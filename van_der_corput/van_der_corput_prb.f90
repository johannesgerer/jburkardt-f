program main

!*****************************************************************************80
!
!! MAIN is the main program for VAN_DER_CORPUT_PRB.
!
!  Discussion:
!
!    VAN_DER_CORPUT_PRB calls a set of problems for VAN_DER_CORPUT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VAN_DER_CORPUT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the VAN_DER_CORPUT library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test045 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VAN_DER_CORPUT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests VAN_DER_CORPUT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  VAN_DER_CORPUT computes the elements of a'
  write ( *, '(a)' ) '  van der Corput sequence.'
  write ( *, '(a)' ) '  Each call produces the next value.  By default,'
  write ( *, '(a)' ) '  the base is 2, and the sequence starts at element 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we call VAN_DER_CORPUT several times.'
  write ( *, '(a)' ) '  We also get and print the seed, although this is'
  write ( *, '(a)' ) '  generally of little interest.'

  base = 2
  call van_der_corput_base_set ( base )
  seed = 0
  call van_der_corput_seed_set ( seed )
  n = 11

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Seed   van der Corput'
  write ( *, '(a)' ) ' '
  do i = 1, n
    call van_der_corput ( r )
    write ( *, '(i6,g16.8)' ) seed+i-1, r
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests VAN_DER_CORPUT_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  VAN_DER_CORPUT_SEQUENCE computes several elements of '
  write ( *, '(a)' ) '  a van der Corput sequence on a single call.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  By default, the base is 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we call VAN_DER_CORPUT_SEQUENCE once.'

  base = 2
  call van_der_corput_base_set ( base )

  seed = 0
  call van_der_corput_seed_set ( seed )

  call van_der_corput_sequence ( n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Element   van der Corput'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,g16.8)' ) seed+i-1, r(i)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests VAN_DER_CORPUT_SEQUENCE, *_SEED_SET, *_GET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 11

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(nmax)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  VAN_DER_CORPUT_SEED_SET specifies the next element of'
  write ( *, '(a)' ) '    the van der Corput sequence to compute.'
  write ( *, '(a)' ) '  VAN_DER_CORPUT_SEED_GET reports the next element of the'
  write ( *, '(a)' ) '    van der Corput sequence that will be computed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  By default, the sequence starts at element 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we demonstrate computing elements'
  write ( *, '(a)' ) '  affects the seed, and how resetting the seed determines'
  write ( *, '(a)' ) '  the next element computed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We compute 11 elements.'
  write ( *, '(a)' ) ' '

  base = 2
  call van_der_corput_base_set ( base )

  seed = 0
  call van_der_corput_seed_set ( seed )

  n = 11
  call van_der_corput_sequence ( n, r )

  do i = 1, n
    write ( *, '(i6,g16.8)' ) seed+i-1, r(i)
  end do

  call van_der_corput_seed_get ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  VAN_DER_CORPUT_SEED_GET: current seed is ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We jump back to element 6 and compute 10 elements.'
  write ( *, '(a)' ) ' '

  seed = 6
  call van_der_corput_seed_set ( seed )

  n = 10
  call van_der_corput_sequence ( n, r )

  do i = 1, n
    write ( *, '(i6,g16.8)' ) seed+i-1, r(i)
  end do

  call van_der_corput_seed_get ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  VAN_DER_CORPUT_SEED_GET: the current seed is ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We restart at element 0 and compute 6 elements.'
  write ( *, '(a)' ) ' '

  seed = 0
  call van_der_corput_seed_set ( seed )

  n = 6
  call van_der_corput_sequence ( n, r )

  do i = 1, n
    write ( *, '(i6,g16.8)' ) seed+i-1, r(i)
  end do

  call van_der_corput_seed_get ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  VAN_DER_CORPUT_SEED_GET: the current seed is ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We jump to element 100 and compute 5 elements.'
  write ( *, '(a)' ) ' '

  seed = 100
  call van_der_corput_seed_set ( seed )

  n = 5
  call van_der_corput_sequence ( n, r )

  do i = 1, n
    write ( *, '(i6,g16.8)' ) seed+i-1, r(i)
  end do

  call van_der_corput_seed_get ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  VAN_DER_CORPUT_SEED_GET: the current seed is ', seed

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests VAN_DER_CORPUT, *_BASE_GET, *_BASE_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  VAN_DER_CORPUT_BASE_GET gets the current base.'
  write ( *, '(a)' ) '  VAN_DER_CORPUT_BASE_SET sets the current base.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The van der Corput base is usually prime, but this is'
  write ( *, '(a)' ) '  not required.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compute a van der Corput sequence'
  write ( *, '(a)' ) '  with the default base, then change the base,'
  write ( *, '(a)' ) '  reset the seed, and recompute the sequence.'

  base = 2
  call van_der_corput_base_set ( base )

  seed = 0
  call van_der_corput_seed_set ( seed )
  n = 10

  call van_der_corput_base_get ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  VAN_DER_CORPUT_BASE_GET: Current base is ', base
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Seed   van der Corput'
  write ( *, '(a)' ) ' '
  do i = 1, n
    call van_der_corput ( r )
    write ( *, '(i6,g16.8)' ) seed+i-1, r
  end do

  base = 3
  call van_der_corput_base_set ( base )

  seed = 0
  call van_der_corput_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Reset base to ', base
  write ( *, '(a,i12)' ) '  Reset seed to ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Seed   van der Corput '
  write ( *, '(a)' ) ' '
  do i = 1, n
    call van_der_corput ( r )
    write ( *, '(i6,g16.8)' ) seed+i-1, r
  end do

  base = 4
  call van_der_corput_base_set ( base )

  seed = 0
  call van_der_corput_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Reset base to ', base
  write ( *, '(a,i12)' ) '  Reset seed to ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Seed   van der Corput'
  write ( *, '(a)' ) ' '
  do i = 1, n
    call van_der_corput ( r )
    write ( *, '(i6,g16.8)' ) seed+i-1, r
  end do

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests VAN_DER_CORPUT, *_SEED_GET, *_SEED_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045'
  write ( *, '(a)' ) '  VAN_DER_CORPUT_SEED_GET gets the current seed.'
  write ( *, '(a)' ) '  VAN_DER_CORPUT_SEED_SET sets the current seed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compute a van der Corput sequence'
  write ( *, '(a)' ) '  with the default seed, then reset the seed,'
  write ( *, '(a)' ) '  and recompute the sequence.'

  base = 2
  call van_der_corput_base_set ( base )
  seed = 0
  call van_der_corput_seed_set ( seed )
  n = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  VAN_DER_CORPUT_BASE_GET: Current base is ', base
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Seed   van der Corput'
  write ( *, '(a)' ) ' '
  do i = 1, n
    call van_der_corput ( r )
    write ( *, '(i6,g16.8)' ) seed+i-1, r
  end do

  seed = 100
  call van_der_corput_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Reset seed to ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Seed   van der Corput '
  write ( *, '(a)' ) ' '
  do i = 1, n
    call van_der_corput ( r )
    write ( *, '(i6,g16.8)' ) seed+i-1, r
  end do

  seed = 3
  call van_der_corput_seed_set ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Reset seed to ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Seed   van der Corput'
  write ( *, '(a)' ) ' '
  do i = 1, n
    call van_der_corput ( r )
    write ( *, '(i6,g16.8)' ) seed+i-1, r
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests I4_TO_VAN_DER_CORPUT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  I4_TO_VAN_DER_CORPUT returns the I-th element'
  write ( *, '(a)' ) '  of a van der Corput sequence to a given base.'

  n = 11

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Base    Seed   R'
  write ( *, '(a)' ) ' '

  do base = 2, 5

    write ( *, '(a)' ) ' '
    write ( *, '(i6)' ) base

    do i = 1, n
      seed = i - 1
      call i4_to_van_der_corput ( seed, base, r )
      write ( *, '(6x,2x,i6,g16.8)' ) seed, r
    end do

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests I4_TO_VAN_DER_CORPUT_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  I4_TO_VAN_DER_CORPUT_SEQUENCE returns N elements'
  write ( *, '(a)' ) '  of a van der Corput sequence to a given base.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Base    Seed   R'
  write ( *, '(a)' ) ' '

  do base = 2, 5

    write ( *, '(a)' ) ' '
    write ( *, '(i6)' ) base

    seed = 0

    call i4_to_van_der_corput_sequence ( seed, base, n, r )

    do i = 1, n
      write ( *, '(6x,2x,i6,g16.8)' ) seed+i-1, r(i)
    end do

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests CIRCLE_UNIT_VAN_DER_CORPUT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) average(ndim)
  integer ( kind = 4 ) base
  real ( kind = 8 ) dot_average
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v(ndim)
  real ( kind = 8 ) x(ndim)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  For the unit circle:'
  write ( *, '(a)' ) '  CIRCLE_UNIT_VAN_DER_CORPUT samples;'
  write ( *, '(a)' ) ' '

  base = 2
  call van_der_corput_base_set ( base )
  seed = 0
  call van_der_corput_seed_set ( seed )
  n = 13

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    call circle_unit_van_der_corput ( x )
    write ( *, '(i6,2f10.6)' ) seed+i-1, x(1:ndim)
  end do

  n = 1000
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of sample points = ', n

  seed = 0
  call van_der_corput_seed_set ( seed )

  average(1:ndim) = 0.0D+00

  do i = 1, n
    call circle_unit_van_der_corput ( x )
    average(1:ndim) = average(1:ndim) + x(1:ndim)
  end do

  average(1:ndim) = average(1:ndim) / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as N increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f10.6)' ) '  Average:        ', average(1:ndim)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now choose a random direction, sample the same'
  write ( *, '(a)' ) '  number of points, and compute the dot product with'
  write ( *, '(a)' ) '  the direction.'
  write ( *, '(a)' ) '  Take the absolute value of each dot product '
  write ( *, '(a)' ) '  and sum and average.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We expect a value near 2 / PI = 0.6366...'

  do j = 1, 5

    call get_seed ( seed )

    call van_der_corput_seed_set ( seed )

    call circle_unit_van_der_corput ( v )

    seed = 0
    call van_der_corput_seed_set ( seed )

    dot_average = 0.0D+00

    do i = 1, n
      call circle_unit_van_der_corput ( x )
      dot_average = dot_average + abs ( dot_product ( x(1:ndim), v(1:ndim) ) )
    end do

    dot_average = dot_average / real ( n, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2f10.6)' ) '  V:                ', v(1:ndim)
    write ( *, '(a, f10.6)' ) '  Average |(XdotV)| ', dot_average

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests I4_TO_VAN_DER_CORPUT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 75

  integer ( kind = 4 ) base
  integer ( kind = 4 ), dimension ( test_num ) :: base_test = (/ &
     2,   2,   2,   2,   2,  2,  2,  2,  2, &
     3,   3,   3,   3,   3,  3,  3,  3,  3, &
     4,   4,   4,   4,   4,  4,  4,  4,  4, &
     2,   3,   4,   5,   7, 11, 13,        &
     2,   3,   4,   5,   7, 11, 13,        &
     2,   3,   4,   5,   7, 11, 13,        &
     2,   3,   4,   5,   7, 11, 13,        &
    29,  29,  29,  29,  29,                &
    71,  71,  71,  71,  71,                &
   173, 173, 173, 173, 173,                &
   409, 409, 409, 409, 409 /)
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), dimension ( test_num ) :: seed_test = (/ &
        0,    1,      2,     3,     4,     5,     6,     7,     8, &
        0,    1,      2,     3,     4,     5,     6,     7,     8, &
        0,    1,      2,     3,     4,     5,     6,     7,     8, &
       10,    10,    10,    10,    10,    10,    10,               &
      100,   100,   100,   100,   100,   100,   100,               &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,               &
    10000, 10000, 10000, 10000, 10000, 10000, 10000,               &
     1000,  1001,  1002,  1003,  1004,                             &
     1000,  1001,  1002,  1003,  1004,                             &
     1000,  1001,  1002,  1003,  1004,                             &
     1000,  1001,  1002,  1003,  1004 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  I4_TO_VAN_DER_CORPUT computes the I-th element'
  write ( *, '(a)' ) '  of a van der Corput sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we try a variety of seeds and bases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* First, we look at the initial sequence of values'
  write ( *, '(a)' ) '  for bases 2, 3, and 4.  Look at how the values '
  write ( *, '(a)' ) '  vary for a fixed base, and then compare the variation'
  write ( *, '(a)' ) '  for corresponding seeds and different bases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* Second, we look at what happens for larger seeds.'
  write ( *, '(a)' ) '  We compare seeds 10, 100, 1000 and 10000 for a variety '
  write ( *, '(a)' ) '  of small bases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* Third, we take relatively large primes (the 10th,'
  write ( *, '(a)' ) '  20th, 40th and 80th) and look at the values of the'
  write ( *, '(a)' ) '  1000, 1001, 1002, 1003, and 1004th values.  As the'
  write ( *, '(a)' ) '  prime base increases, successive values tend to be'
  write ( *, '(a)' ) '  close, which can cause undesirable correlations.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Base      Seed  VDC(Base,Seed)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    base = base_test(test)
    seed = seed_test(test)

    call van_der_corput_base_set ( base )
    call van_der_corput_seed_set ( seed )

    call i4_to_van_der_corput ( seed, base, r )

    write ( *, '(2x,i8,2x,i8,g24.16)' ) base, seed, r

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests VDC_NUMERATOR_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ), allocatable :: r(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  VDC_NUMERATOR_SEQUENCE returns N elements'
  write ( *, '(a)' ) '  of a van der Corput numerator sequence in base 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N:  Sequence'
  write ( *, '(a)' ) ' '

  do n = 1, 20

    allocate ( r(1:n) )

    call vdc_numerator_sequence ( n, r )

    write ( *, '(2x,i2,a1,20(2x,i2))' ) n, ':', r(1:n)

    deallocate ( r )

  end do

  return
end
