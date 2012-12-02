program main

!*****************************************************************************80
!
!! MAIN is the main program for HALTON_PRB.
!
!  Discussion:
!
!    HALTON_PRB calls a set of problems for HALTON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HALTON_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HALTON library.'

  call test01 ( )
  call test0125 ( )
  call test0126 ( )
  call test02 ( )
  call test025 ( )
  call test03 ( )
  call test04 ( )
  call test045 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )

  call test11 ( )
  call test12 ( )
  call test13 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HALTON_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests HALTON, HALTON_STEP_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_max = 3
  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) base(dim_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_max)
  integer ( kind = 4 ) seed(dim_max)
  integer ( kind = 4 ), dimension ( test_num ) :: step_vec = (/ 0, 5, 1000, 1000000 /)
  integer ( kind = 4 ) step
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  HALTON computes the next element of a Halton sequence.'
  write ( *, '(a)' ) '  HALTON_STEP_SET sets the step.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we try several values of STEP.'
  write ( *, '(a)' ) '  We repeat the test for several dimensions.'
  write ( *, '(a)' ) '  We assume defaults SEED, LEAP and BASE.'

  do dim_num = 1, dim_max

    do test = 1, test_num

      call halton_dim_num_set ( dim_num )
      n = 11
      step = step_vec(test)
      call halton_step_set ( step )
      seed(1:dim_num) = 0
      call halton_seed_set ( seed )
      do i = 1, dim_num
        base(i) = prime(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
      write ( *, '(a,i12)' ) '  N =    ', n
      write ( *, '(a,i12)' ) '  STEP = ', step
      call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
      call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          STEP  Halton'
      write ( *, '(a)' ) ' '
      do j = 1, n
        call halton ( dim_num, r )
        write ( *, '(2x,i12,4f12.8)' ) step+j-1, r(1:dim_num)
      end do

    end do

  end do

  return
end
subroutine test0125 ( )

!*****************************************************************************80
!
!! TEST0125 tests I4_TO_HALTON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_max = 3

  integer ( kind = 4 ) base(dim_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_max)
  integer ( kind = 4 ) seed(dim_max)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0125'
  write ( *, '(a)' ) '  I4_TO_HALTON computes a Halton sequence.'
  write ( *, '(a)' ) '  The user specifies all data explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we call I4_TO_HALTON repeatedly.'
  write ( *, '(a)' ) '  We use distinct primes as bases.'

  do dim_num = 1, dim_max

    n = 11
    step = 0
    seed(1:dim_num) = 0
    leap(1:dim_num) = 1
    do i = 1, dim_num
      base(i) = prime ( i )
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
    write ( *, '(a,i12)' ) '  N =    ', n
    write ( *, '(a,i12)' ) '  STEP = ', step
    call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
    call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
    call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    STEP     Halton'
    write ( *, '(a)' ) ' '
    do j = 1, n
      step = j-1
      call i4_to_halton ( dim_num, step, seed, leap, base, r )
      write ( *, '(2x,i6,2x,4f10.6)' ) step, r(1:dim_num)
    end do

  end do

  return
end
subroutine test0126 ( )

!*****************************************************************************80
!
!! TEST0126 tests I4_TO_HALTON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0126'
  write ( *, '(a)' ) '  I4_TO_HALTON computes a Halton sequence.'
  write ( *, '(a)' ) '  The user specifies all data explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we call I4_TO_HALTON repeatedly.'
  write ( *, '(a)' ) '  We use the same value for all the bases.'

  n = 11
  seed(1:dim_num) = 0
  leap(1:dim_num) = 1
  base(1:dim_num) = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP      Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    step = j-1
    call i4_to_halton ( dim_num, step, seed, leap, base, r )
    write ( *, '(2x,i6,4f10.6)' ) step, r(1:dim_num)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests HALTON_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) j
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  HALTON_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Halton sequence on a single call.'

  call halton_dim_num_set ( dim_num )
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,3g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests I4_TO_HALTON_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 )  j
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  I4_TO_HALTON_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Halton sequence on a single call.'
  write ( *, '(a)' ) '  All arguments are specified explicitly.'

  do test = 1, 4

    if ( test == 1 ) then

      step = 0
      seed(1:dim_num) = 0
      leap(1:dim_num) = 1
      do i = 1, dim_num
        base(i) = prime ( i )
      end do

    else if ( test == 2 ) then

      step = 0
      do i = 1, dim_num
        seed(i) = i
      end do
      leap(1:dim_num) = 1
      do i = 1, dim_num
        base(i) = prime ( i )
      end do

    else if ( test == 3 ) then

      step = 0
      seed(1:dim_num) = 1
      leap(1:dim_num) = 3
      do i = 1, dim_num
        base(i) = prime ( i )
      end do

    else if ( test == 4 ) then

      step = 0
      seed(1:dim_num) = (/ 1, 2, 3 /)
      leap(1:dim_num) = 1
      base(1:dim_num) = (/ 2, 2, 2 /)

    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
    write ( *, '(a,i12)' ) '  N =    ', n
    write ( *, '(a,i12)' ) '  STEP = ', step
    call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
    call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
    call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

    call i4_to_halton_sequence ( dim_num, n, step, seed, leap, base, r )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    STEP   Halton'
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(2x,i6,3g14.6)' ) step+j-1, r(1:dim_num,j)
    end do

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests HALTON_SEQUENCE, HALTON_STEP_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 1
  integer ( kind = 4 ), parameter :: nmax = 11

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(dim_num,nmax)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HALTON_STEP_SET specifies the next element of'
  write ( *, '(a)' ) '    the Halton sequence to compute.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we demonstrate how resetting'
  write ( *, '(a)' ) '  STEP determines the next element computed.'

  call halton_dim_num_set ( dim_num )
  n = 11
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )
  base(1:dim_num) = 2
  call halton_base_set ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP  Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  n = 11
  step = 6
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP  Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  n = 6
  step = 0
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP  Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  n = 5
  step = 100
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP  Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests HALTON, HALTON_BASE_GET, HALTON_BASE_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 1

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  HALTON_BASE_GET gets the current Halton bases.'
  write ( *, '(a)' ) '  HALTON_BASE_SET sets the current Halton bases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compute a Halton sequence'
  write ( *, '(a)' ) '  with the default base, then change the base,'
  write ( *, '(a)' ) '  reset the seed, and recompute the sequence.'

  call halton_dim_num_set ( dim_num )
  n = 10
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )
  call halton_base_get ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N    = ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    call halton ( dim_num, r )
    write ( *, '(2x,i6,g14.6)' ) step+j-1, r(1:dim_num)
  end do

  n = 10
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )
  base(1) = 3
  call halton_base_set ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N    = ', n
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    call halton ( dim_num, r )
    write ( *, '(2x,i6,g14.6)' ) step+j-1, r(1:dim_num)
  end do

  n = 10
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )
  base(1) = 4
  call halton_base_set ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N    = ', n
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    call halton ( dim_num, r )
    write ( *, '(2x,i6,g14.6)' ) step+j-1, r(1:dim_num)
  end do

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests HALTON_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 101
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045'
  write ( *, '(a)' ) '  HALTON_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Halton sequence on a single call.'

  call halton_dim_num_set ( dim_num )
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )
  do i = 1, dim_num
    base(i) = prime ( i )
  end do
  call halton_base_set ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N    = ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2f7.4)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests HALTON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  HALTON computes the elements of a vector '
  write ( *, '(a)' ) '  Halton sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each call produces the next value.  By default,'
  write ( *, '(a)' ) '  the bases are the first DIM_NUM primes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we call HALTON several times,'
  write ( *, '(a)' ) '  with the default bases.'

  call halton_dim_num_set ( dim_num )
  n = 11
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )
  do i = 1, dim_num
    base(i) = prime ( i )
  end do
  call halton_base_set ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    call halton ( dim_num, r )
    write ( *, '(2x,i6,4g14.6)' ) step+j-1, r(1:dim_num)
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests HALTON_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) j
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  HALTON_SEQUENCE computes the next N elements'
  write ( *, '(a)' ) '  of a vector Halton sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each call produces the next value.  By default,'
  write ( *, '(a)' ) '  the bases are the first DIM_NUM primes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we demonstrate how one call can compute'
  write ( *, '(a)' ) '  many successive vector elements of the sequence.'

  call halton_dim_num_set ( dim_num )
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests HALTON_STEP_SET, HALTON_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 4
  integer ( kind = 4 ), parameter :: nmax = 10

  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(dim_num,nmax)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  HALTON_STEP_SET specifies which element of the'
  write ( *, '(a)' ) '    Halton subsequence to compute.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we show how STEP chooses the next element.'

  call halton_dim_num_set ( dim_num )
  n = 10
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  n = 10
  step = 6
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  n = 6
  step = 0
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  n = 5
  step = 100
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  call halton_sequence ( dim_num, n, r )

  do j = 1, n
    write ( *, '(2x,i6,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests HALTON_BASE_GET, HALTON_BASE_SET, HALTON_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  HALTON_BASE_GET gets the current bases.'
  write ( *, '(a)' ) '  HALTON_BASE_SET sets the current bases.'
  write ( *, '(a)' ) '  HALTON_SEQUENCE computes the next N elements'
  write ( *, '(a)' ) '  of a vector Halton sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compute the first 10 elements of the'
  write ( *, '(a)' ) '  default sequence, then change bases, reset the seed'
  write ( *, '(a)' ) '  and recompute the first 10 elements.'

  call halton_dim_num_set ( dim_num )
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )
  do i = 1, dim_num
    base(i) = prime ( i )
  end do
  call halton_base_set ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  step = 0
  call halton_step_set ( step )
  do i = 1, dim_num
    base(i) = prime ( 2 * i )
  end do
  call halton_base_set ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call halton_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests U1_TO_SPHERE_UNIT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 1
  integer ( kind = 4 ), parameter :: dim_num2 = 2

  real ( kind = 8 ) average(dim_num2)
  real ( kind = 8 ) dot_average
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) seed_val
  integer ( kind = 4 ) step
  real ( kind = 8 ) u(dim_num)
  real ( kind = 8 ) v(dim_num2)
  real ( kind = 8 ) x(dim_num2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For the unit sphere in 2 dimensions (the circle):'
  write ( *, '(a)' ) '  HALTON generates "U1" points,'
  write ( *, '(a)' ) '  U1_TO_SPHERE_UNIT_2D samples the circle;'

  call halton_dim_num_set ( dim_num )
  n = 5
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    call halton ( dim_num, u )
    call u1_to_sphere_unit_2d ( u, x )
    write ( *, '(2x,2f8.4)' ) x(1:dim_num2)
  end do

  n = 1000
  step = 0
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  average(1:dim_num2) = 0.0D+00

  do j = 1, n
    call halton ( dim_num, u )
    call u1_to_sphere_unit_2d ( u, x )
    average(1:dim_num2) = average(1:dim_num2) + x(1:dim_num2)
  end do

  average(1:dim_num2) = average(1:dim_num2) / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as N increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  Average:        ', average(1:dim_num2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now choose a random direction, sample the same'
  write ( *, '(a)' ) '  number of points, and compute the dot product with'
  write ( *, '(a)' ) '  the direction.'
  write ( *, '(a)' ) '  Take the absolute value of each dot product '
  write ( *, '(a)' ) '  and sum and average.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We expect a value near 2 / PI = 0.6366...'

  do j2 = 1, 5

    call get_seed ( seed_val )
    step = seed_val + 111 * j2
    call halton_step_set ( step )

    call halton ( dim_num, u )
    call u1_to_sphere_unit_2d ( u, v )

    step = 0
    call halton_step_set ( step )

    dot_average = 0.0D+00

    do j = 1, n
      call halton ( dim_num, u )
      call u1_to_sphere_unit_2d ( u, x )
      dot_average = dot_average + abs ( dot_product ( x(1:dim_num2), v(1:dim_num2) ) )
    end do

    dot_average = dot_average / real ( n, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2f8.4)' ) '  Random V:         ', v(1:dim_num2)
    write ( *, '(a, f8.4)' ) '  Average |(XdotV)| ', dot_average

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests U2_TO_BALL_UNIT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) atan4
  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) average_r
  real ( kind = 8 ) average_theta
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step
  real ( kind = 8 ) theta
  real ( kind = 8 ) u(dim_num)
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For the unit ball in 2 dimensions (the disk):'
  write ( *, '(a)' ) '  U2_TO_BALL_UNIT_2D samples;'

  call halton_dim_num_set ( dim_num )
  n = 5
  step = 0
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N    = ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    call halton ( dim_num, u )
    call u2_to_ball_unit_2d ( u, x )
    write ( *, '(2x,2f8.4)' ) x(1:dim_num)
  end do

  n = 1000
  step = 0
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  average(1:dim_num) = 0.0D+00

  do j = 1, n
    call halton ( dim_num, u )
    call u2_to_ball_unit_2d ( u, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as N increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  Average:        ', average(1:dim_num)

  step = 0
  call halton_step_set ( step )

  average_r = 0.0D+00

  do j = 1, n
    call halton ( dim_num, u )
    call u2_to_ball_unit_2d ( u, x )
    r = sqrt ( sum ( x(1:dim_num)**2 ) )
    average_r = average_r + r
  end do

  average_r = average_r / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the distance of the points from'
  write ( *, '(a,f8.4)' ) '  the center, which should be DIM_NUM/(DIM_NUM+1) = ', &
   real ( dim_num, kind = 8 ) / real ( dim_num + 1, kind = 8 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_r

  step = 0
  call halton_step_set ( step )

  average_theta = 0.0D+00

  do j = 1, n
    call halton ( dim_num, u )
    call u2_to_ball_unit_2d ( u, x )
    theta = atan4 ( x(2), x(1) )
    average_theta = average_theta + theta
  end do

  average_theta = average_theta / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the angle THETA,'
  write ( *, '(a)' ) '  which should approach PI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_theta

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests U2_TO_SPHERE_UNIT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: dim_num2 = 3

  real ( kind = 8 ) average(dim_num2)
  real ( kind = 8 ) dot_average
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) seed_val
  integer ( kind = 4 ) step
  real ( kind = 8 ) u(dim_num)
  real ( kind = 8 ) v(dim_num2)
  real ( kind = 8 ) x(dim_num2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  For the unit sphere in 3 dimensions:'
  write ( *, '(a)' ) '  U2_TO_SPHERE_UNIT_3D samples;'

  call halton_dim_num_set ( dim_num )
  n = 5
  step = 123456789
  call halton_step_set ( step )
  seed(1:dim_num) = 0
  call halton_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    call halton ( dim_num, u )
    call u2_to_sphere_unit_3d ( u, x )
    write ( *, '(2x,3f8.4)' ) x(1:dim_num2)
  end do

  n = 1000
  step = 0
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  average(1:dim_num2) = 0.0D+00

  do j = 1, n
    call halton ( dim_num, u )
    call u2_to_sphere_unit_3d ( u, x )
    average(1:dim_num2) = average(1:dim_num2) + x(1:dim_num2)
  end do

  average(1:dim_num2) = average(1:dim_num2) / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as N increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  Average:        ', average(1:dim_num2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now choose a random direction, sample the same'
  write ( *, '(a)' ) '  number of points, and compute the dot product with'
  write ( *, '(a)' ) '  the direction.'
  write ( *, '(a)' ) '  Take the absolute value of each dot product '
  write ( *, '(a)' ) '  and sum and average.'

  do j2 = 1, 5

    call get_seed ( seed_val )
    step = seed_val + 111 * j2
    call halton_step_set ( step )

    call halton ( dim_num, u )
    call u2_to_sphere_unit_3d ( u, v )

    step = 0
    call halton_step_set ( step )

    dot_average = 0.0D+00

    do j = 1, n
      call halton ( dim_num, u )
      call u2_to_sphere_unit_3d ( u, x )
      dot_average = dot_average + abs ( dot_product ( x(1:dim_num2), v(1:dim_num2) ) )
    end do

    dot_average = dot_average / real ( n, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3f8.4)' ) '  Random V:         ', v(1:dim_num2)
    write ( *, '(a, f8.4)' ) '  Average |(XdotV)| ', dot_average

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests U3_TO_BALL_UNIT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) atan4
  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) average_phi
  real ( kind = 8 ) average_r
  real ( kind = 8 ) average_theta
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) phi
  real ( kind = 8 ) r
  integer ( kind = 4 ) step
  real ( kind = 8 ) theta
  real ( kind = 8 ) u(dim_num)
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  For the unit ball in 3 dimensions:'
  write ( *, '(a)' ) '  U3_TO_BALL_UNIT_3D samples;'

  call halton_dim_num_set ( dim_num )
  n = 5
  step = 0
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    call halton ( dim_num, u )
    call u3_to_ball_unit_3d ( u, x )
    write ( *, '(2x,3f8.4)' ) x(1:dim_num)
  end do

  n = 1000
  step = 0
  call halton_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step

  average(1:dim_num) = 0.0D+00

  do j = 1, n
    call halton ( dim_num, u )
    call u3_to_ball_unit_3d ( u, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as N increases.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  Average:        ', average(1:dim_num)

  step = 0
  call halton_step_set ( step )

  average_r = 0.0D+00

  do j = 1, n
    call halton ( dim_num, u )
    call u3_to_ball_unit_3d ( u, x )
    r = sqrt ( sum ( x(1:dim_num)**2 ) )
    average_r = average_r + r
  end do

  average_r = average_r / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the distance of the points from'
  write ( *, '(a,f8.4)' ) '  the center, which should be DIM_NUM/(DIM_NUM+1) = ', &
   real ( dim_num, kind = 8 ) / real ( dim_num + 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_r

  step = 0
  call halton_step_set ( step )

  average_theta = 0.0D+00

  do j = 1, n
    call halton ( dim_num, u )
    call u3_to_ball_unit_3d ( u, x )
    theta = atan4 ( x(2), x(1) )
    average_theta = average_theta + theta
  end do

  average_theta = average_theta / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the angle THETA,'
  write ( *, '(a)' ) '  which should approach PI.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_theta

  step = 0
  call halton_step_set ( step )

  average_phi = 0.0D+00

  do j = 1, n
    call halton ( dim_num, u )
    call u3_to_ball_unit_3d ( u, x )
    r = sqrt ( sum ( x(1:dim_num)**2 ) )
    if ( r == 0.0D+00 ) then
      phi = 0.0D+00
    else
      phi = acos ( x(3) / r )
    end if
    average_phi = average_phi + phi
  end do

  average_phi = average_phi / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the angle PHI,'
  write ( *, '(a)' ) '  which should approach PI/2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_phi

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests HALHAM_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) base(dim_num)
  character ( len = 80 ) :: file_name = 'halton_03_00010.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  HALHAM_WRITE writes a Halton or Hammersley'
  write ( *, '(a)' ) '  dataset to a file'

  step = 0
  seed(1:dim_num) = 0
  leap(1:dim_num) = 1
  do i = 1, dim_num
    base(i) = prime ( i )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call i4_to_halton_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Halton'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,3g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  call halham_write ( dim_num, n, step, seed, leap, base, r, file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was written to "' // trim ( file_name ) // '".'

  return
end
