program main

!*****************************************************************************80
!
!! MAIN is the main program for BERNSTEIN_PRB.
!
!  Discussion:
!
!    BERNSTEIN_PRB calls the BERNSTEIN test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BERNSTEIN_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BERNSTEIN library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BERNSTEIN_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BERNSTEIN_POLY and BERNSTEIN_POLY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable, dimension ( : ) :: bvec
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  BERNSTEIN_POLY evaluates the Bernstein polynomials'
  write ( *, '(a)' ) '  based on the interval [0,1].'
  write ( *, '(a)' ) '  BERNSTEIN_POLY_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     K     X       Exact         BP01(N,K)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernstein_poly_values ( n_data, n, k, x, b )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( bvec(0:n) )

    call bernstein_poly ( n, x, bvec )

    write ( *, '(2x,i4,2x,i4,2x,f7.4,g14.6,2x,g14.6)' ) n, k, x, b, bvec(k)

    deallocate ( bvec )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BPAB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bern(0:n)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  BPAB evaluates Bernstein polynomials over an'
  write ( *, '(a)' ) '  arbitrary interval [A,B].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we demonstrate that '
  write ( *, '(a)' ) '    BPAB(N,K,A1,B1)(X1) = BPAB(N,K,A2,B2)(X2)'
  write ( *, '(a)' ) '  provided only that'
  write ( *, '(a)' ) '    (X1-A1)/(B1-A1) = (X2-A2)/(B2-A2).'

  x = 0.3D+00
  a = 0.0D+00
  b = 1.0D+00
  call bpab ( n, a, b, x, bern )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     K     A        B        X       BPAB(N,K,A,B)(X)'
  write ( *, '(a)' ) ' ' 
  do k = 0, n
    write ( *, '(2x,i4,2x,i4,2x,f7.4,2x,f7.4,2x,f7.4,2x,g14.6)' ) n, k, a, b, x, bern(k)
  end do
 
  x = 1.3D+00
  a = 1.0D+00
  b = 2.0D+00
  call bpab ( n, a, b, x, bern )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     K     A        B        X       BPAB(N,K,A,B)(X)'
  write ( *, '(a)' ) ' ' 
  do k = 0, n
    write ( *, '(2x,i4,2x,i4,2x,f7.4,2x,f7.4,2x,f7.4,2x,g14.6)' ) n, k, a, b, x, bern(k)
  end do

  x = 2.6D+00
  a = 2.0D+00
  b = 4.0D+00
  call bpab ( n, a, b, x, bern )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     K     A        B        X       BPAB(N,K,A,B)(X)'
  write ( *, '(a)' ) ' '
 
  do k = 0, n
    write ( *, '(2x,i4,2x,i4,2x,f7.4,2x,f7.4,2x,f7.4,2x,g14.6)' ) n, k, a, b, x, bern(k)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests the Partition-of-Unity property.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( : ) :: bvec
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  BERNSTEIN_POLY evaluates the Bernstein polynomials'
  write ( *, '(a)' ) '  based on the interval [0,1].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we test the partition of unity property.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     X          Sum ( 0 <= K <= N ) BP01(N,K)(X)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do n = 0, 10

    allocate ( bvec(0:n) )

    x = r8_uniform_01 ( seed )

    call bernstein_poly ( n, x, bvec )

    write ( *, '(2x,i4,2x,f7.4,2x,g14.6)' ) n, x, sum ( bvec(0:n) )

    deallocate ( bvec )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests BPAB_APPROX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxdata = 20
  integer ( kind = 4 ), parameter :: nval = 501

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) error_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ndata
  integer ( kind = 4 ) nsample
  real ( kind = 8 ) xdata(0:maxdata)
  real ( kind = 8 ) xval(nval)
  real ( kind = 8 ) ydata(0:maxdata)
  real ( kind = 8 ) yval(nval)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  BPAB_APPROX evaluates the Bernstein polynomial'
  write ( *, '(a)' ) '  approximant to a function F(X).'

  a = 1.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N      Max Error'
  write ( *, '(a)' ) ' '

  do ndata = 0, maxdata
!
!  Generate data values.
!
    do i = 0, ndata

      if ( ndata == 0 ) then
        xdata(i) = 0.5D+00 * ( a + b )
      else
        xdata(i) = ( real ( ndata - i, kind = 8 ) * a   &
                   + real (         i, kind = 8 ) * b ) &
                   / real ( ndata,     kind = 8 )
      end if

      ydata(i) = sin ( xdata(i) )

    end do
!
!  Compare the true function and the approximant.
!
    call r8vec_linspace ( nval, a, b, xval )

    error_max = 0.0D+00

    call bpab_approx ( ndata, a, b, ydata, nval, xval, yval )

    error_max = maxval ( abs ( yval(1:nval) - sin ( xval(1:nval) ) ) )

    write ( *, '(2x,i4,2x,g14.6)' ) ndata, error_max

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests BERNSTEIN_MATRIX and BERNSTEIN_MATRIX_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) a_norm_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: b
  real ( kind = 8 ) b_norm_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: c
  real ( kind = 8 ) error_norm_frobenius
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8mat_norm_fro

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  BERNSTEIN_MATRIX returns a matrix A which transforms a'
  write ( *, '(a)' ) '  polynomial coefficient vector from the power basis to'
  write ( *, '(a)' ) '  the Bernstein basis.'
  write ( *, '(a)' ) '  BERNSTEIN_MATRIX_INVERSE computes the inverse B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     ||A||            ||B||      ||I-A*B||'
  write ( *, '(a)' ) ' '

  do n = 5, 15

    allocate ( a(1:n,1:n) )
    allocate ( b(1:n,1:n) )
    allocate ( c(1:n,1:n) )

    call bernstein_matrix ( n, a )
    a_norm_frobenius = r8mat_norm_fro ( n, n, a )

    call bernstein_matrix_inverse ( n, b )
    b_norm_frobenius = r8mat_norm_fro ( n, n, b )

    c = matmul ( a, b )
    call r8mat_is_identity ( n, c, error_norm_frobenius )

    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      n, a_norm_frobenius, b_norm_frobenius, error_norm_frobenius

    deallocate ( a )
    deallocate ( b )
    deallocate ( c )

  end do

  return
end

