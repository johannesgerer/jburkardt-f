program main

!*****************************************************************************80
!
!! MAIN is the main program for IEEE_UNIFORM_PRB.
!
!  Discussion:
!
!    IEEE_UNIFORM_PRB calls sample problems for the IEEE_UNIFORM library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IEEE_UNIFORM_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the IEEE_UNIFORM library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IEEE_UNIFORM_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests I4_SEED_ADVANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i4_seed_advance
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_new
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  I4_SEED_ADVANCE advances the seed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step        SEED input       SEED output'
  write ( *, '(a)' ) ' '

  seed_new = 12345

  do test = 1, 10

    seed = seed_new
    seed_new = i4_seed_advance ( seed )

    write ( *, '(2x,i4,2x,i16,2x,i16)' ) test, seed, seed_new

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R4_IEEE_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) r4
  real ( kind = 4 ) r4_ieee_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R4_IEEE_UNIFORM computes an IEEE uniform'
  write ( *, '(a)' ) '  real ( kind = 4 ) value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step        SEED input       R4 output'
  write ( *, '(a)' ) ' '

  seed = 12345

  do test = 1, 20

    r4 = r4_ieee_uniform ( seed )

    write ( *, '(2x,i4,2x,i16,2x,g14.6)' ) test, seed, r4

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests R4_IEEE_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 1024

  real ( kind = 4 ) b
  integer ( kind = 4 ) count(9)
  integer ( kind = 4 ) j
  real ( kind = 4 ) r4_ieee_uniform
  real ( kind = 4 ) r4vec(test_num)
  integer ( kind = 4 ) seed
  real ( kind = 4 ) t
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  R4_IEEE_UNIFORM computes an IEEE uniform'
  write ( *, '(a)' ) '  real ( kind = 4 ) value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate a lot of values, and count where they fall.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of values to generate = ', test_num
  write ( *, '(a)' ) ' '

  seed = 12345

  do test = 1, test_num

    r4vec(test) = r4_ieee_uniform ( seed )

  end do
!
!  Sign check.
!
  count(1:3) = 0

  do test = 1, test_num
    if ( r4vec(test) < 0.0E+00 ) then
      count(1) = count(1) + 1
    else if ( r4vec(test) == 0.0E+00 ) then
      count(2) = count(2) + 1
    else if ( 0.0E+00 < r4vec(test) ) then
      count(3) = count(3) + 1
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sign check:'
  write ( *, '(a)' ) ' '
  write ( *, '(i8,a)' ) count(1), ' values less than 0'
  write ( *, '(i8,a)' ) count(2), ' values equal to 0'
  write ( *, '(i8,a)' ) count(3), ' values greater than 0.'
!
!  Exponent check.
!
  count(1) = 0

  b = 1.0E+00 / 16.0E+00**10
  t = 1.0E+00 / 16.0E+00**9

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exponent check'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bottom   Top     #( B <= |X| < T )'
  write ( *, '(a)' ) ' '

  do j = 1, 25
    count(1) = 0
    do test = 1, test_num
      if ( b <= abs ( r4vec(test) ) .and. abs ( r4vec(test) ) < t ) then
        count(1) = count(1) + 1
      end if
    end do
    write ( *, * ) b, t, count(1)
    b = t
    t = 16.0E+00 * t
  end do

  return
end
