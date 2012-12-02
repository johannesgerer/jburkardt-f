program main

!*******************************************************************************
!
!! MAIN is the main program for NORMAL_PRB.
!
!  Discussion:
!
!    NORMAL_PRB calls sample problems for the NORMAL library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NORMAL_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version;'
  write ( *, '(a)' ) '  Test the NORMAL library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  if ( .false. ) then
    call test04 ( )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 is being skipped.'
    write ( *, '(a)' ) '  The routine being tested is misbehaving.'
  end if

  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NORMAL_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*******************************************************************************
!
!! TEST01 tests C4_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) c4_normal_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  C4_NORMAL_01 computes pseudorandom complex values '
  write ( *, '(a)' ) '  normally distributed in the unit circle.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8,2x,f14.8)' ) i, c4_normal_01 ( seed )
  end do

  return
end
subroutine test02 ( )

!*******************************************************************************
!
!! TEST02 tests C8_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 8 ) c8_normal_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  C8_NORMAL_01 computes pseudorandom double precision'
  write ( *, '(a)' ) '  complex values normally distributed in the unit'
  write ( *, '(a)' ) '  circle.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8,2x,f14.8)' ) i, c8_normal_01 ( seed )
  end do

  return
end
subroutine test03 ( )

!*******************************************************************************
!
!! TEST03 tests I4_NORMAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_normal
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  I4_NORMAL computes integer pseudonormal values '
  write ( *, '(a)' ) '  with mean A and standard deviation B.'

  a = 70.0E+00
  b = 10.0E+00
  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The mean A = ', a
  write ( *, '(a,g14.6)' ) '  The standard deviation B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,i8)' ) i, i4_normal ( a, b, seed )
  end do

  return
end
subroutine test04 ( )

!*******************************************************************************
!
!! TEST04 tests I8_NORMAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 8 ) i
  integer ( kind = 8 ) i8_normal
  integer ( kind = 8 ) seed
  integer ( kind = 8 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  I8_NORMAL computes integer pseudonormal values '
  write ( *, '(a)' ) '  with mean A and standard deviation B.'

  a = 70.0D+00
  b = 10.0D+00
  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The mean A = ', a
  write ( *, '(a,g14.6)' ) '  The standard deviation B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,i8)' ) i, i8_normal ( a, b, seed )
  end do

  return
end
subroutine test05 ( )

!*******************************************************************************
!
!! TEST05 tests R4_NORMAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_normal
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  R4_NORMAL computes real pseudonormal values '
  write ( *, '(a)' ) '  with mean A and standard deviation B.'

  a = 10.0E+00
  b = 2.0E+00
  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The mean A = ', a
  write ( *, '(a,g14.6)' ) '  The standard deviation B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8)' ) i, r4_normal ( a, b, seed )
  end do

  return
end
subroutine test06 ( )

!*******************************************************************************
!
!! TEST06 tests R4_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_normal_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  R4_NORMAL_01 computes normal pseudorandom values '
  write ( *, '(a)' ) '  in the interval [0,1].'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8)' ) i, r4_normal_01 ( seed )
  end do

  return
end
subroutine test07 ( )

!*******************************************************************************
!
!! TEST07 tests R8_NORMAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_normal
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  R8_NORMAL computes pseudonormal values '
  write ( *, '(a)' ) '  with mean A and standard deviation B.'

  a = 10.0D+00
  b = 2.0D+00
  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The mean A = ', a
  write ( *, '(a,g14.6)' ) '  The standard deviation B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8)' ) i, r8_normal ( a, b, seed )
  end do

  return
end
subroutine test08 ( )

!*******************************************************************************
!
!! TEST08 tests R8_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  R8_NORMAL_01 computes pseudonormal values '
  write ( *, '(a)' ) '  with mean 0 and standard deviation 1.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8)' ) i, r8_normal_01 ( seed )
  end do

  return
end
subroutine test09 ( )

!*******************************************************************************
!
!! TEST09 tests R8_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10000

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_out
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) u_avg
  real ( kind = 8 ) u_var

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  R8_NORMAL_01 computes a sequence of '
  write ( *, '(a)' ) '  normally distributed pseudorandom numbers.'
!
!  Start with a given seed.
!
  seed = 12345

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Initial SEED = ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 10 values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I         Input        Output    R8_NORMAL_01'
  write ( *, '(a)' ) '                  SEED          SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    seed_in = seed
    u(i) = r8_normal_01 ( seed )
    seed_out = seed

    write ( *, '(2x,i8,2x,i12,2x,i12,2x,f14.8)' ) i, seed_in, seed_out, u(i)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12,a)' ) '  Now compute ', n, ' elements.'
  write ( *, '(a)' ) ' '

  do i = 1, n
    u(i) = r8_normal_01 ( seed )
  end do

  u_avg = sum ( u(1:n) ) / real ( n, kind = 8 )

  u_var = sum ( ( u(1:n) - u_avg )**2 ) / real ( n - 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.8)' ) '  Average value = ', u_avg
  write ( *, '(a,f14.8)' ) '  Expecting       ', 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.8)' ) '  Variance =      ', u_var
  write ( *, '(a,f14.8)' ) '  Expecting       ', 1.0D+00

  return
end
subroutine test10 ( )

!*******************************************************************************
!
!! TEST10 tests R8_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789
  integer ( kind = 4 ) seed_input
  integer ( kind = 4 ) seed_output
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  R8_NORMAL_01 computes pseudonormal values '
  write ( *, '(a)' ) '  with mean 0 and standard deviation 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that we can change the seed'
  write ( *, '(a)' ) '  and get the desired results.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         I    Seed(in)   Seed(out)   R8_NORMAL_01'
  write ( *, '(a)' ) ' '

  do i = 1, 5

    seed_input = seed
    value = r8_normal_01 ( seed )
    seed_output = seed

    write ( *, '(2x,i8,2x,i12,2x,i12,2x,f14.8)' ) &
      i, seed_input, seed_output, value

  end do

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Resetting seed to repeat, after an ODD number of steps.'
  write ( *, '(a)' ) ' '

  do i = 6, 10

    seed_input = seed
    value = r8_normal_01 ( seed )
    seed_output = seed

    write ( *, '(2x,i8,2x,i12,2x,i12,2x,f14.8)' ) &
      i, seed_input, seed_output, value

  end do

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Resetting seed to repeat, after an EVEN number of steps.'
  write ( *, '(a)' ) ' '

  do i = 11, 15

    seed_input = seed
    value = r8_normal_01 ( seed )
    seed_output = seed

    write ( *, '(2x,i8,2x,i12,2x,i12,2x,f14.8)' ) &
      i, seed_input, seed_output, value

  end do

  return
end
subroutine test11 ( )

!*******************************************************************************
!
!! TEST11 tests R8_NORMAL_01 and R8MAT_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 100
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  R8_NORMAL_01 computes pseudorandom values '
  write ( *, '(a)' ) '    one at a time.'
  write ( *, '(a)' ) '  R8MAT_NORMAL_01 computes a matrix of values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For the same initial seed, the results should '
  write ( *, '(a)' ) '  be identical, but R8MAT_NORMAL_01 might be faster.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

  seed = seed_init

  do j = 1, n
    do i = 1, m
      a(i,j) = r8_normal_01 ( seed )
    end do
  end do

  seed = seed_init;
  call r8mat_normal_01 ( m, n, seed, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J      A(I,J)        B(I,J)'
  write ( *, '(a)' ) '                       (R8_NORMAL_01)  (R8MAT_NORMAL_01)'
  write ( *, '(a)' ) ' '

  do k = 0, 10
    i = ( k * m + ( 10 - k ) * 1 ) / 10
    j = ( k * n + ( 10 - k ) * 1 ) / 10
    write ( *, '(2x,i8,2x,i8,2x,f14.8,2x,f14.8)' ) i, j, a(i,j), b(i,j)
  end do
  
  return
end
subroutine test12 ( )

!*******************************************************************************
!
!! TEST12 tests R8_NORMAL_01 and R8VEC_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  R8_NORMAL_01 computes pseudorandom values '
  write ( *, '(a)' ) '  one at a time.'
  write ( *, '(a)' ) '  R8VEC_NORMAL_01 computes a vector of values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For the same initial seed, the results should '
  write ( *, '(a)' ) '  be identical, but R8VEC_NORMAL_01 might be faster.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

  seed = seed_init

  do i = 1, n
    a(i) = r8_normal_01 ( seed )
  end do

  seed = seed_init;
  call r8vec_normal_01 ( n, seed, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I      A(I)            B(I)'
  write ( *, '(a)' ) '             (R8_NORMAL_01)  (R8VEC_NORMAL_01)'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8,2x,f14.8)' ) i, a(i), b(i)
  end do
  
  return
end
