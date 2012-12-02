program main

!*****************************************************************************80
!
!! RANDLC_PRB calls sample problems for the RANDLC library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDLC_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the RANDLC library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDLC_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests RANDLC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) randlc
  real ( kind = 8 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  RANDLC computes pseudorandom values '
  write ( *, '(a)' ) '  in the interval [0,1].'

  seed = 123456789.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,f15.0)' ) '  The initial seed is ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I          RANDLC'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.6)' ) i, randlc ( seed )
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests RANDLC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 1000

  integer ( kind = 4 ) i
  real ( kind = 8 ) randlc
  real ( kind = 8 ) seed
  real ( kind = 8 ) seed_in
  real ( kind = 8 ) seed_out
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) u_avg
  real ( kind = 8 ) u_var

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  RANDLC computes a sequence of uniformly'
  write ( *, '(a)' ) '  distributed pseudorandom numbers.'

  seed = 123456789.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,f15.0)' ) '  Initial SEED = ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 10 values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    u(i) = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) &
      i, seed_in, seed_out, u(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) '  Now call RANDLC ', n, ' times.'

  u_avg = 0.0D+00
  do i = 1, n
    u(i) = randlc ( seed )
    u_avg = u_avg + u(i)
  end do

  u_avg = u_avg / real ( n, kind = 8 )

  u_var = 0.0D+00
  do i = 1, n
    u_var = u_var + ( u(i) - u_avg )**2
  end do
  u_var = u_var / real ( n - 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Average value = ', u_avg
  write ( *, '(a,g14.6)' ) '  Expecting       ', 0.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Variance =      ', u_var
  write ( *, '(a,g14.6)' ) '  Expecting       ', 1.0D+00 / 12.0D+00

  return
  end
  subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests RANDLC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) randlc
  real ( kind = 8 ) seed
  real ( kind = 8 ) seed_in
  real ( kind = 8 ) seed_out
  real ( kind = 8 ) seed_save
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  RANDLC computes a sequence of pseudorandom numbers'
  write ( *, '(a)' ) '  but all computations depend on the seed value.'
  write ( *, '(a)' ) '  In this test, we show how a sequence of "random"'
  write ( *, '(a)' ) '  values can be manipulated by accessing the seed.'

  seed = 1066.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,f15.0)' ) '  Set SEED to ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call RANDLC 10 times, and watch SEED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    seed_in = seed

    if ( i == 5 ) then
      seed_save = seed
    end if
    x = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) i, seed_in, seed_out, x
  end do

  seed = seed_save

  write ( *, '(a)' ) ' '
  write ( *, '(a,f15.0)' ) '  Reset SEED to its value at step 5, = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call RANDLC 10 times, and watch how SEED'
  write ( *, '(a)' ) '  and RANDLC restart themselves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    x = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) i, seed_in, seed_out, x
  end do

  seed = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  What happens with an initial zero SEED?'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    x = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) i, seed_in, seed_out, x
  end do

  seed = -123456789.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  What happens with an initial negative SEED?'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    x = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) i, seed_in, seed_out, x
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! RANDLC_TEST04 tests RANDLC_JUMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) klog
  real ( kind = 8 ) randlc
  real ( kind = 8 ) randlc_jump
  real ( kind = 8 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDLC_TEST04'
  write ( *, '(a)' ) '  RANDLC_JUMP jumps directly to the K-th value'
  write ( *, '(a)' ) '  returned by RANDLC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         K X(hard way)     X(jump)'
  write ( *, '(a)' ) ' '

  k = 1

  do klog = 1, 10

    seed = 123456789.0D+00
    do i = 1, k
      x1 = randlc ( seed )
    end do

    seed = 123456789.0D+00
    x2 = randlc_jump ( seed, k )

    write ( *, '(2x,i8,2x,f10.6,2x,f10.6)' ) k, x1, x2

    k = k * 2

  end do

  return
end
