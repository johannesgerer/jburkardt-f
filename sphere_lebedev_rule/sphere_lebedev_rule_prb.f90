program main

!*****************************************************************************80
!
!! MAIN is the main program for SPHERE_LEBEDEV_RULE_PRB.
!
!  Discussion:
!
!    SPHERE_LEBEDEV_RULE_PRB tests routines from the SPHERE_LEBEDEV_RULE library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_LEBEDEV_RULE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPHERE_LEBEDEV_RULE library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_LEBEDEV_RULE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests AVAILABLE_TABLE, ORDER_TABLE, PRECISION_TABLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!   12 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) available
  integer ( kind = 4 ) available_table
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_table
  integer ( kind = 4 ) precision
  integer ( kind = 4 ) precision_table
  integer ( kind = 4 ) rule
  integer ( kind = 4 ), parameter :: rule_max = 65

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  List Lebedev rule properties.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rule Avail Order  Prec'
  write ( *, '(a)' ) ' '
  do rule = 1, rule_max
    available = available_table ( rule )
    order = order_table ( rule )
    precision = precision_table ( rule )
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) rule, available, order, precision
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!c TEST02 tests the SPHERE_LEBEDEV_RULE functions.
!
!  Modified:
!
!    13 September 2010
!
!  Author:
!
!    Dmitri Laikov
!
!  Reference:
!
!    Vyacheslav Lebedev, Dmitri Laikov,
!    A quadrature formula for the sphere of the 131st
!    algebraic order of accuracy,
!    Russian Academy of Sciences Doklady Mathematics,
!    Volume 59, Number 3, 1999, pages 477-481.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 65
  integer ( kind = 4 ), parameter :: mmax = &
    ( ( nmax * 2 + 3 ) * ( nmax * 2 + 3 ) / 3 )

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) available
  integer ( kind = 4 ) available_table
  real ( kind = 8 ) beta
  real ( kind = 8 ) err
  real ( kind = 8 ) err_max
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral_exact
  real ( kind = 8 ) integral_approx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_table
  integer ( kind = 4 ) precision_table
  real ( kind = 8 ) s(0:nmax+1)
  real ( kind = 8 ) w(mmax)
  real ( kind = 8 ) x(mmax)
  real ( kind = 8 ) xn(mmax,0:nmax)
  real ( kind = 8 ) y(mmax)
  real ( kind = 8 ) yn(mmax,0:nmax)
  real ( kind = 8 ) z(mmax)
  real ( kind = 8 ) zn(mmax,0:nmax)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Generate each available rule and test for accuracy.'

  do n = 1, nmax

    available = available_table ( n )

    if ( available == 1 ) then

      order = order_table ( n )

      call ld_by_order ( order, x, y, z, w )

      s(0) = 1.0D+00
      do k = 1, n + 1
        s(k) = ( 2 * k - 1 ) * s(k-1)
      end do
!
!  For each abscissa X(M), compute the values 1, X(M)^2, X(M)^4, ..., X(M)^2*N.
!
      do m = 1, order
        xn(m,0) = 1.0D+00
        yn(m,0) = 1.0D+00
        zn(m,0) = 1.0D+00
        do k = 1, n
          xn(m,k) = xn(m,k-1) * x(m) * x(m)
          yn(m,k) = yn(m,k-1) * y(m) * y(m)
          zn(m,k) = zn(m,k-1) * z(m) * z(m)
        end do
      end do

      err_max = 0.0D+00
      do i = 0, n
        do j = 0, n - i
          k = n - i - j
!
!  Apply Lebedev rule to x^2i y^2j z^2k.
!
          integral_approx = 0.0D+00
          do m = 1, order
            integral_approx = integral_approx &
              + w(m) * xn(m,i) * yn(m,j) * zn(m,k)
          end do
!
!  Compute exact value of integral (aside from factor of 4 pi!).
!
          integral_exact = s(i) * s(j) * s(k) / s(1+i+j+k)
!
!  Record the maximum error for this rule.
!
          err = abs ( ( integral_approx - integral_exact ) / integral_exact )
          err_max = max ( err_max, err )

        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4,a,i4,a,g14.6)' ) &
        '  Order = ', order, &
        '  LMAXW = ', precision_table ( n ), &
        '  max error = ', err_max
!
!  Convert (X,Y,Z) to (Theta,Phi) and print the data.
!
      if ( order <= 50 ) then
        do m = 1, order
          call xyz_to_tp ( x(m), y(m), z(m), alpha, beta )
          write ( *, '(g24.15,g24.15,g24.15)' ) alpha, beta, w(m)
        end do
      end if

    end if

  end do

  return
end
