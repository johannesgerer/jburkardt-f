program main

!*****************************************************************************80
!
!! MAIN is the main program for BLAS1_C_PRB.
!
!  Discussion:
!
!    BLAS1_C_PRB tests the BLAS1 single precision complex routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLAS1_C_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BLAS1_C library.'

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
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLAS1_C_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CABS1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) c
  complex ( kind = 4 ) c4_uniform_01
  real ( kind = 4 ) c_norm
  real ( kind = 4 ) cabs1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) &
    '  CABS1 returns the L1 norm of a single precision complex number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Real      Imaginary              '
  write ( *, '(a)' ) '      Part      Part           CABS1(Z)'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    c = 5.0E+00 * c4_uniform_01 ( seed )

    c_norm = cabs1 ( c )

    write ( *, '(2x,2f10.4,5x,f10.4)' ) c, c_norm

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CABS2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) c
  complex ( kind = 4 ) c4_uniform_01
  real ( kind = 4 ) c_norm
  real ( kind = 4 ) cabs2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) &
    '  CABS2 returns the L2 norm of a single precision complex number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Real      Imaginary              '
  write ( *, '(a)' ) '      Part      Part           CABS2(Z)'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    c = 5.0E+00 * c4_uniform_01 ( seed )

    c_norm = cabs2 ( c )

    write ( *, '(2x,2f10.4,5x,f10.4)' ) c, c_norm

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests CAXPY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  complex ( kind = 4 ) s
  complex ( kind = 4 ), dimension ( n ) :: x = (/ &
    (  2.0E+00, -1.0E+00 ), &
    ( -4.0E+00, -2.0E+00 ), &
    (  3.0E+00,  1.0E+00 ), &
    (  2.0E+00,  2.0E+00 ), &
    ( -1.0E+00, -1.0E+00 ) /)
  complex ( kind = 4 ), dimension ( n ) :: y = (/ &
    ( -1.0E+00,  0.0E+00 ), &
    (  0.0E+00, -3.0E+00 ), &
    (  4.0E+00,  0.0E+00 ), &
    ( -3.0E+00,  4.0E+00 ), &
    ( -2.0E+00,  0.0E+00 ) /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  CAXPY adds a multiple of '
  write ( *, '(a)' ) '  one single precision complex vector to another.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, y(i)
  end do

  s = ( 0.50E+00, -1.00E+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  The scalar multiplier is: ', s

  call caxpy ( n, s, x, 1, y, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A * X + Y = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2f10.6)' ) i, y(i)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests CCOPY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) a(5,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 4 ) x(10)
  complex ( kind = 4 ) y(10)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  CCOPY copies a single precision complex vector.'

  do i = 1, 10
    x(i) = cmplx ( 10 * i, i )
  end do

  do i = 1, 10
    y(i) = cmplx ( 20 * i, 2 * i )
  end do

  do i = 1, 5
    do j = 1, 5
      a(i,j) = cmplx ( 10 * i,  j )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,2g14.6)' ) i, x(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y = '
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,2g14.6)' ) i, y(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A = '
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,10f7.1)' ) a(i,1:5)
  end do

  call ccopy ( 5, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  CCOPY ( 5, X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,2g14.6)' ) i, y(i)
  end do

  do i = 1, 10
    y(i) = cmplx ( 20 * i, 2 * i )
  end do

  call ccopy ( 3, x, 2, y, 3 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  CCOPY ( 3, X, 2, Y, 3 )'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,2g14.6)' ) i, y(i)
  end do

  call ccopy ( 5, x, 1, a, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  CCOPY ( 5, X, 1, A, 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A = '
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,10f7.1)' ) a(i,1:5)
  end do

  do i = 1, 5
    do j = 1, 5
      a(i,j) = cmplx ( 10 * i,  j )
    end do
  end do

  call ccopy ( 5, x, 2, a, 5 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  CCOPY ( 5, X, 2, A, 5 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A = '
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,10f7.1)' ) a(i,1:5)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests CDOTC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) i
  complex ( kind = 4 ) x_norm
  complex ( kind = 4 ) xy_dot
  complex ( kind = 4 ), dimension ( n ) :: x = (/ &
    (  2.0E+00, -1.0E+00 ), &
    ( -4.0E+00, -2.0E+00 ), &
    (  3.0E+00,  1.0E+00 ), &
    (  2.0E+00,  2.0E+00 ), &
    ( -1.0E+00, -1.0E+00 ) /)
  complex ( kind = 4 ), dimension ( n ) :: y = (/ &
    ( -1.0E+00,  0.0E+00 ), &
    (  0.0E+00, -3.0E+00 ), &
    (  4.0E+00,  0.0E+00 ), &
    ( -3.0E+00,  4.0E+00 ), &
    ( -2.0E+00,  0.0E+00 ) /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  CDOTC computes the conjugated dot product of '
  write ( *, '(a)' ) '  two single precision complex vectors.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  x_norm = cdotc ( n, x, 1, x, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The square of the norm of X, computed as'
  write ( *, '(a,f10.4,2x,f10.4)' ) '  CDOTC(X,X) = ', x_norm

  xy_dot = cdotc ( n, x, 1, y, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, y(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,2x,f10.4)' ) '  The dot product X.Y* is ', xy_dot

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests CDOTU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  complex ( kind = 4 ) cdotu
  integer ( kind = 4 ) i
  complex ( kind = 4 ) x_norm
  complex ( kind = 4 ) xy_dot
  complex ( kind = 4 ), dimension ( n ) :: x = (/ &
    (  2.0E+00, -1.0E+00 ), &
    ( -4.0E+00, -2.0E+00 ), &
    (  3.0E+00,  1.0E+00 ), &
    (  2.0E+00,  2.0E+00 ), &
    ( -1.0E+00, -1.0E+00 ) /)
  complex ( kind = 4 ), dimension ( n ) :: y = (/ &
    ( -1.0E+00,  0.0E+00 ), &
    (  0.0E+00, -3.0E+00 ), &
    (  4.0E+00,  0.0E+00 ), &
    ( -3.0E+00,  4.0E+00 ), &
    ( -2.0E+00,  0.0E+00 ) /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  CDOTU computes the unconjugated dot product of '
  write ( *, '(a)' ) '  two single precision complex vectors.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  x_norm = cdotu ( n, x, 1, x, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The unconjugated dot product ( X dot X )'
  write ( *, '(a)' ) '  (which is NOT the square of the norm of X!):'
  write ( *, '(a,f10.4,2x,f10.4)' ) '  CDOTU(X,X) = ', x_norm

  xy_dot = cdotu ( n, x, 1, y, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, y(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,2x,f10.4)' ) '  The dot product ( X dot Y ) is ', xy_dot

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests CMACH.
!
!  Discussion:
!
!    The CMACH routine is not part of the official BLAS release.
!    It was used for the testing routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) cmach
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  CMACH computes several machine-dependent'
  write ( *, '(a)' ) '  single precision complex arithmetic parameters.'

  write ( *, '(a)' ) ' '
  write ( *, * ) '  CMACH(1)  = machine epsilon = ', cmach ( 1 )
  write ( *, * ) '  CMACH(2)  = a tiny value    = ', cmach ( 2 )
  write ( *, * ) '  CMACH(3)  = a huge value    = ', cmach ( 3 )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests CROTG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) a
  complex ( kind = 4 ) b
  real ( kind = 4 ) c
  complex ( kind = 4 ) c4_uniform_01
  complex ( kind = 4 ) r
  complex ( kind = 4 ) s
  complex ( kind = 4 ) sa
  complex ( kind = 4 ) sb
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) &
    '  CROTG generates a single precision complex Givens rotation'
  write ( *, '(a)' ) '    (  C  S ) * ( A ) = ( R )'
  write ( *, '(a)' ) '    ( -S  C )   ( B )   ( 0 )'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    a = c4_uniform_01 ( seed )
    b = c4_uniform_01 ( seed )

    sa = a
    sb = b

    call crotg ( sa, sb, c, s )

    r = sa

    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  A =  ', a
    write ( *, '(a,2g14.6)' ) '  B =  ', b
    write ( *, '(a, g14.6)' ) '  C =  ', c
    write ( *, '(a,2g14.6)' ) '  S =  ', s
    write ( *, '(a,2g14.6)' ) '  R =  ', r
    write ( *, '(a,2g14.6)' ) '         C *A+S*B = ',          c   * a + s * b
    write ( *, '(a,2g14.6)' ) '  -conjg(S)*A+C*B = ', -conjg ( s ) * a + c * b

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests CSCAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  complex ( kind = 4 ) da
  integer ( kind = 4 ) i
  complex ( kind = 4 ) x(n)

  do i = 1, n
    x(i) = cmplx ( 10 * i, i )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  CSCAL multiplies a single precision complex scalar'
  write ( *, '(a)' ) '  times a single precision vector.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  da = cmplx ( 5.0E+00, 0.0E+00 )
  call cscal ( n, da, x, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f8.4,a)' ) '  CSCAL ( N, (', da, '), X, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  do i = 1, n
    x(i) = cmplx ( 10 * i, i )
  end do

  da = cmplx ( -2.0E+00, 1.0E+00 )
  call cscal ( 3, da, x, 2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f8.4,a)' ) '  CSCAL ( 3, (', da, '), X, 2 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests CSIGN1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) c1
  complex ( kind = 4 ) c2
  complex ( kind = 4 ) c3
  complex ( kind = 4 ) c4_uniform_01
  complex ( kind = 4 ) csign1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  CSIGN1 ( C1, C2 ) transfers the sign of '
  write ( *, '(a)' ) '  a single precision complex C2'
  write ( *, '(a)' ) '  to the CABS1 magnitude of a single precision complex C1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '           C1                    C2                    C3'
  write ( *, '(a)' ) &
    '  --------------------  --------------------  --------------------'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    c1 = 5.0E+00 * c4_uniform_01 ( seed )
    c2 = 5.0E+00 * c4_uniform_01 ( seed )
    c3 = csign1 ( c1, c2 )

    write ( *, '(2x,2f10.4,2x,2f10.4,2x,2f10.4)' ) c1, c2, c3

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests CSIGN2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) c1
  complex ( kind = 4 ) c2
  complex ( kind = 4 ) c3
  complex ( kind = 4 ) c4_uniform_01
  complex ( kind = 4 ) csign2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  CSIGN2 ( C1, C2 ) transfers the sign of'
  write ( *, '(a)' ) '  a single precision complex C2'
  write ( *, '(a)' ) '  to the CABS2 magnitude of a single precision complex C1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '           C1                    C2                    C3'
  write ( *, '(a)' ) &
    '  --------------------  --------------------  --------------------'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    c1 = 5.0E+00 * c4_uniform_01 ( seed )
    c2 = 5.0E+00 * c4_uniform_01 ( seed )
    c3 = csign2 ( c1, c2 )

    write ( *, '(2x,2f10.4,2x,2f10.4,2x,2f10.4)' ) c1, c2, c3

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests CSROT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 4 ) c
  integer ( kind = 4 ) i
  real ( kind = 4 ) s
  complex ( kind = 4 ) x(n)
  complex ( kind = 4 ) y(n)

  do i = 1, n
    x(i) = cmplx ( 10 * i, i )
  end do

  do i = 1, n
    y(i) = cmplx ( 20 * i, 2 * i )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  CSROT carries out a Givens rotation'
  write ( *, '(a)' ) '  on a single precision complex vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,2f10.1,2x,2f10.1)' ) i, x(i), y(i)
  end do

  c = 0.5E+00
  s = sqrt ( 1.0E+00 - c * c )
  call csrot ( n, x, 1, y, 1, c, s )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a,f8.4,a)' ) '  CSROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,2f10.1,2x,2f10.1)' ) i, x(i), y(i)
  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests CSSCAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 4 ) da
  integer ( kind = 4 ) i
  complex ( kind = 4 ) x(n)

  do i = 1, n
    x(i) = cmplx ( 10 * i, i )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  CSSCAL multiplies a single precision real scalar'
  write ( *, '(a)' ) '  times a single precision complex vector.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  da = 5.0E+00
  call csscal ( n, da, x, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  CSSCAL ( N, ', da, ', X, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  do i = 1, n
    x(i) = cmplx ( 10 * i, i )
  end do

  da = -2.0E+00
  call csscal ( 3, da, x, 2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  CSSCAL ( 3, ', da, ', X, 2 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests CSWAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  complex ( kind = 4 ) x(n)
  complex ( kind = 4 ) y(n)

  do i = 1, n
    x(i) = cmplx ( 10 * i, i )
  end do

  do i = 1, n
    y(i) = cmplx ( 20 * i, 2 * i )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  CSWAP swaps two single precision complex vectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,2f7.1,2x,2f7.1)' ) i, x(i), y(i)
  end do

  call cswap ( n, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CSWAP ( N, X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,2f7.1,2x,2f7.1)' ) i, x(i), y(i)
  end do

  do i = 1, n
    x(i) = cmplx ( 10 * i, i )
  end do

  do i = 1, n
    y(i) = cmplx ( 20 * i, 2 * i )
  end do

  call cswap ( 3, x, 2, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  CSWAP ( 3, X, 2, Y, 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,2f7.1,2x,2f7.1)' ) i, x(i), y(i)
  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests ICAMAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 4 ) cabs1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) icamax
  complex ( kind = 4 ), dimension ( n ) :: x = (/ &
    (  2.0E+00, -1.0E+00 ), &
    ( -4.0E+00, -2.0E+00 ), &
    (  3.0E+00,  1.0E+00 ), &
    (  2.0E+00,  2.0E+00 ), &
    ( -1.0E+00, -1.0E+00 ) /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  ICAMAX returns the index of the entry of '
  write ( *, '(a)' ) '  maximum magnitude in a single precision complex vector.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The entries and CABS1 magnitudes:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '( 2x,i6, 2f8.4,5x,f8.4 )' ) i, x(i), cabs1 ( x(i) )
  end do

  incx = 1

  i = icamax ( n, x, incx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)'    ) '  The index of maximum magnitude = ', i

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests SCASUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ma = 5
  integer ( kind = 4 ), parameter :: na = 4
  integer ( kind = 4 ), parameter :: nx = 8

  complex ( kind = 4 ), dimension ( ma, na ) :: a = reshape ( (/ &
    ( -3.0E+00,  4.0E+00 ), &
    (  2.0E+00,  0.0E+00 ), &
    (  3.0E+00, -4.0E+00 ), &
    (  2.0E+00,  0.0E+00 ), &
    (  2.0E+00, -1.0E+00 ), &
    ( -1.0E+00,  1.0E+00 ), &
    (  0.0E+00,  5.0E+00 ), &
    ( -4.0E+00, -2.0E+00 ), &
    ( -4.0E+00,  1.0E+00 ), &
    ( -4.0E+00, -3.0E+00 ), &
    (  0.0E+00, -2.0E+00 ), &
    (  1.0E+00,  3.0E+00 ), &
    ( -3.0E+00,  3.0E+00 ), &
    ( -3.0E+00,  3.0E+00 ), &
    ( -1.0E+00, -2.0E+00 ), &
    ( -1.0E+00,  2.0E+00 ), &
    (  2.0E+00, -4.0E+00 ), &
    (  0.0E+00, -1.0E+00 ), &
    (  0.0E+00, -1.0E+00 ), &
    ( -2.0E+00,  4.0E+00 ) /), (/ ma, na /) )
  integer ( kind = 4 ) i
  real ( kind = 4 ) scasum
  complex ( kind = 4 ), dimension ( nx ) :: x = (/ &
    (  2.0E+00, -1.0E+00 ), &
    ( -4.0E+00, -2.0E+00 ), &
    (  3.0E+00,  1.0E+00 ), &
    (  2.0E+00,  2.0E+00 ), &
    ( -1.0E+00, -1.0E+00 ), &
    ( -1.0E+00,  0.0E+00 ), &
    (  0.0E+00, -3.0E+00 ), &
    (  4.0E+00,  0.0E+00 ) /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  SCASUM adds the absolute values of elements '
  write ( *, '(a)' ) '  of a single precision complex vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, nx
    write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  SCASUM ( NX,   X, 1    ) = ', &
    scasum ( nx,   x, 1 )
  write ( *, '(a,g14.6)' ) '  SCASUM ( NX/2, X, 2    ) = ', &
    scasum ( nx/2, x, 2 )
  write ( *, '(a,g14.6)' ) '  SCASUM ( 2,    X, NX/2 ) = ', &
    scasum ( 2,    x, nx/2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate with a matrix A:'
  write ( *, '(a)' ) ' '
  do i = 1, ma
    write ( *, '(4(2x,f6.1,2x,f6.1))' ) a(i,1:na)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  SCASUM ( MA, A(1,2), 1 )   = ', &
    scasum ( ma, a(1,2), 1 )
  write ( *, '(a,g14.6)' ) '  SCASUM ( NA, A(2,1), MA ) = ', &
    scasum ( na, a(2,1), ma )

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests SCNRM2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  real ( kind = 4 ) norm
  real ( kind = 4 ) scnrm2
  complex ( kind = 4 ), dimension ( n ) :: x = (/ &
    (  2.0E+00, -1.0E+00 ), &
    ( -4.0E+00, -2.0E+00 ), &
    (  3.0E+00,  1.0E+00 ), &
    (  2.0E+00,  2.0E+00 ), &
    ( -1.0E+00, -1.0E+00 ) /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  SCNRM2 returns the Euclidean norm of a'
  write ( *, '(a)' ) '  single precision complex vector.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vector X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '( 2x, i6, 2x, f6.1, 2x, f6.1 )' ) i, x(i)
  end do

  incx = 1
  norm = scnrm2 ( n, x, incx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)'    ) '  The L2 norm of X is ', norm

  return
end
function c4_uniform_01 ( seed )

!*****************************************************************************80
!
!! C4_UNIFORM_01 returns a unit pseudorandom C4.
!
!  Discussion:
!
!    The angle should be uniformly distributed between 0 and 2 * PI,
!    the square root of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, complex c4_uniform_01, a pseudorandom complex value.
!
  implicit none

  complex ( kind = 4 ) c4_uniform_01
  real ( kind = 4 ) r
  integer ( kind = 4 ) k
  real ( kind = 4 ), parameter :: pi = 3.1415926E+00
  integer ( kind = 4 ) seed
  real ( kind = 4 ) theta

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = sqrt ( real ( real ( seed, kind = 8 ) * 4.656612875D-10 ) )

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  theta = 2.0E+00 * pi * real ( real ( seed, kind = 8 ) * 4.656612875D-10 )

  c4_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ) )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
