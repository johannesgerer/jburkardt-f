program main

!*****************************************************************************80
!
!! MAIN is the main program for SFTPACK_PRB.
!
!  Discussion:
!
!    SFTPACK_PRB tests the SFTPACK routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SFTPACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SFTPACK library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SFTPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests R8VEC_SCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 256

  real ( kind = 8 ), parameter :: ahi = 5.0D+00
  real ( kind = 8 ), parameter :: alo = 0.0D+00
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For slow cosine transforms,'
  write ( *, '(a)' ) '  R8VEC_SCT does a forward or backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  seed = 123456789

  call r8vec_uniform ( n, alo, ahi, seed, c )

  call r8vec_print_part ( n, c, 10, '  The original data:' )
!
!  Compute the coefficients.
!
  call r8vec_sct ( n, c, d )

  call r8vec_print_part ( n, d, 10, '  The cosine coefficients:' )
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.
!
  call r8vec_sct ( n, d, e )

  e(1:n) = e(1:n) / real ( 2 * n, kind = 8 )

  call r8vec_print_part ( n, e, 10, '  The retrieved data:' )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R8VEC_SFTB and R8VEC_SFTF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 36

  real ( kind = 8 ) a(n/2)
  real ( kind = 8 ), parameter :: ahi = 5.0D+00
  real ( kind = 8 ), parameter :: alo = 0.0D+00
  real ( kind = 8 ) azero
  real ( kind = 8 ) b(n/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For real slow Fourier transforms,'
  write ( *, '(a)' ) '  R8VEC_SFTF computes the forward transform.'
  write ( *, '(a)' ) '  R8VEC_SFTB computes the backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data values, N = ', n

  seed = 123456789

  call r8vec_uniform ( n, alo, ahi, seed, x )

  call r8vec_print_part ( n, x, 10, '  The original data:' )
!
!  Compute the slow Fourier transform of the data.
!
  call r8vec_sftf ( n, x, azero, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A (cosine) coefficients:'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,i3,g14.6)' ) 0, azero

  do i = 1, n/2
    write ( *, '(2x,i3,g14.6)' ) i, a(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  B (sine) coefficients:'
  write ( *, '(a)' ) ' '

  do i = 1, n/2
    write ( *, '(2x,i3,g14.6)' ) i, b(i)
  end do
!
!  Now try to retrieve the data from the coefficients.
!
  call r8vec_sftb ( n, azero, a, b, x )

  call r8vec_print_part ( n, x, 10, '  The retrieved data:' )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests R8VEC_SHT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 17

  real ( kind = 8 ), parameter :: ahi = 5.0D+00
  real ( kind = 8 ), parameter :: alo = 0.0D+00
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For real slow Hartley transforms,'
  write ( *, '(a)' ) '  R8VEC_SHT does a forward or backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  seed = 123456789

  call r8vec_uniform ( n, alo, ahi, seed, c )

  call r8vec_print_part ( n, c, 10, '  The original data:' )
!
!  Compute the coefficients.
!
  call r8vec_sht ( n, c, d )

  call r8vec_print_part ( n, d, 10, '  The Hartley coefficients:' )
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.
!
  call r8vec_sht ( n, d, e )

  call r8vec_print_part ( n, e, 10, '  The retrieved data:' )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests R8VEC_SQCTB and R8VEC_SQCTF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 256

  real ( kind = 8 ), parameter :: ahi = 5.0D+00
  real ( kind = 8 ), parameter :: alo = 0.0D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For real slow quarter wave cosine transforms,'
  write ( *, '(a)' ) '  R8VEC_SQCTF does a forward transform;'
  write ( *, '(a)' ) '  R8VEC_SQCTB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  seed = 123456789

  call r8vec_uniform ( n, alo, ahi, seed, x )

  call r8vec_print_part ( n, x, 10, '  The original data:' )
!
!  Compute the coefficients.
!
  call r8vec_sqctf ( n, x, y )

  call r8vec_print_part ( n, y, 10, '  The cosine coefficients:' )
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.
!
  call r8vec_sqctb ( n, y, x )

  call r8vec_print_part ( n, x, 10, '  The retrieved data:' )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests R8VEC_SQSTB and R8VEC_SQSTF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 256

  real ( kind = 8 ), parameter :: ahi = 5.0D+00
  real ( kind = 8 ), parameter :: alo = 0.0D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For real slow quarter wave sine transforms,'
  write ( *, '(a)' ) '  R8VEC_SQSTF does a forward transform;'
  write ( *, '(a)' ) '  R8VEC_SQSTB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  seed = 123456789

  call r8vec_uniform ( n, alo, ahi, seed, x )

  call r8vec_print_part ( n, x, 10, '  The original data:' )
!
!  Compute the coefficients.
!
  call r8vec_sqstf ( n, x, y )

  call r8vec_print_part ( n, y, 10, '  The sine coefficients:' )
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.
!
  call r8vec_sqstb ( n, y, x )

  call r8vec_print_part ( n, x, 10, '  The retrieved data:' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests R8VEC_SST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 256

  real ( kind = 8 ), parameter :: ahi = 5.0D+00
  real ( kind = 8 ), parameter :: alo = 0.0D+00
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For slow sine transforms,'
  write ( *, '(a)' ) '  R8VEC_SST does a forward or backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  seed = 123456789

  call r8vec_uniform ( n, alo, ahi, seed, c )

  call r8vec_print_part ( n, c, 10, '  The original data:' )
!
!  Compute the coefficients.
!
  call r8vec_sst ( n, c, d )

  call r8vec_print_part ( n, d, 10, '  The sine coefficients:' )
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.
!
  call r8vec_sst ( n, d, e )

  e(1:n) = e(1:n) / real ( 2 * ( n + 1 ), kind = 8 )

  call r8vec_print_part ( n, e, 10, '  The retrieved data:' )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests C4VEC_SFTB and C4VEC_SFTF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 36

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(n)
  complex ( kind = 4 ) x2(n)
  complex ( kind = 4 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  For complex slow Fourier transforms,'
  write ( *, '(a)' ) '  C4VEC_SFTF computes the forward transform.'
  write ( *, '(a)' ) '  C4VEC_SFTB computes the backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data values, N = ', n

  seed = 123456789

  call c4vec_uniform_01 ( n, seed, x )

  call c4vec_print_part ( n, x, 10, '  The data X:' )
!
!  Compute the slow Fourier transform of the data.
!
  call c4vec_sftf ( n, x, y )

  call c4vec_print_part ( n, y, 10, '  The Fourier coefficients Y:' )

  call c4vec_sftb ( n, y, x2 )

  call c4vec_print_part ( n, x2, 10, '  The recovered data:' )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests C8VEC_SFTB and C8VEC_SFTF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 36

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) x2(n)
  complex ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  For complex slow Fourier transforms,'
  write ( *, '(a)' ) '  C8VEC_SFTF computes the forward transform.'
  write ( *, '(a)' ) '  C8VEC_SFTB computes the backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data values, N = ', n

  seed = 123456789

  call c8vec_uniform_01 ( n, seed, x )

  call c8vec_print_part ( n, x, 10, '  The data X:' )
!
!  Compute the slow Fourier transform of the data.
!
  call c8vec_sftf ( n, x, y )

  call c8vec_print_part ( n, y, 10, '  The Fourier coefficients Y:' )

  call c8vec_sftb ( n, y, x2 )

  call c8vec_print_part ( n, x2, 10, '  The recovered data:' )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests R4VEC_SFTB and R4VEC_SFTF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 36

  real ( kind = 4 ) a(n/2)
  real ( kind = 4 ), parameter :: ahi = 5.0E+00
  real ( kind = 4 ), parameter :: alo = 0.0E+00
  real ( kind = 4 ) azero
  real ( kind = 4 ) b(n/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For real slow Fourier transforms,'
  write ( *, '(a)' ) '  R4VEC_SFTF computes the forward transform.'
  write ( *, '(a)' ) '  R4VEC_SFTB computes the backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data values, N = ', n

  seed = 123456789

  call r4vec_uniform ( n, alo, ahi, seed, x )

  call r4vec_print_part ( n, x, 10, '  The original data:' )
!
!  Compute the slow Fourier transform of the data.
!
  call r4vec_sftf ( n, x, azero, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A (cosine) coefficients:'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,i3,g14.6)' ) 0, azero

  do i = 1, n/2
    write ( *, '(2x,i3,g14.6)' ) i, a(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  B (sine) coefficients:'
  write ( *, '(a)' ) ' '

  do i = 1, n/2
    write ( *, '(2x,i3,g14.6)' ) i, b(i)
  end do
!
!  Now try to retrieve the data from the coefficients.
!
  call r4vec_sftb ( n, azero, a, b, x )

  call r4vec_print_part ( n, x, 10, '  The retrieved data:' )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests C4MAT_SFTB and C4MAT_SFTF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 10
  integer ( kind = 4 ), parameter :: n2 = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(n1,n2)
  complex ( kind = 4 ) x2(n1,n2)
  complex ( kind = 4 ) y(n1,n2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For complex slow Fourier transforms,'
  write ( *, '(a)' ) '  C4MAT_SFTF computes the forward transform.'
  write ( *, '(a)' ) '  C4MAT_SFTB computes the backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) &
    '  The data has dimension N1 = ', n1, ' by N2 = ', n2

  seed = 123456789

  call c4mat_uniform_01 ( n1, n2, seed, x )

  call c4mat_print_some ( n1, n2, x, 1, 1, 10, 10, '  The data X:' )
!
!  Compute the slow Fourier transform of the data.
!
  call c4mat_sftf ( n1, n2, x, y )

  call c4mat_print_some ( n1, n2, y, 1, 1, 10, 10, '  The Fourier coefficients Y:' )

  call c4mat_sftb ( n1, n2, y, x2 )

  call c4mat_print_some ( n1, n2, x2, 1, 1, 10, 10, '  The recovered data:' )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests C8MAT_SFTB and C8MAT_SFTF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 10
  integer ( kind = 4 ), parameter :: n2 = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n1,n2)
  complex ( kind = 8 ) x2(n1,n2)
  complex ( kind = 8 ) y(n1,n2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  For complex slow Fourier transforms,'
  write ( *, '(a)' ) '  C8MAT_SFTF computes the forward transform.'
  write ( *, '(a)' ) '  C8MAT_SFTB computes the backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) &
    '  The data has dimension N1 = ', n1, ' by N2 = ', n2

  seed = 123456789

  call c8mat_uniform_01 ( n1, n2, seed, x )

  call c8mat_print_some ( n1, n2, x, 1, 1, 10, 10, '  The data X:' )
!
!  Compute the slow Fourier transform of the data.
!
  call c8mat_sftf ( n1, n2, x, y )

  call c8mat_print_some ( n1, n2, y, 1, 1, 10, 10, '  The Fourier coefficients Y:' )

  call c8mat_sftb ( n1, n2, y, x2 )

  call c8mat_print_some ( n1, n2, x2, 1, 1, 10, 10, '  The recovered data:' )

  return
end

