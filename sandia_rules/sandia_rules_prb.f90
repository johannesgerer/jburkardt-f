program main

!*****************************************************************************80
!
!! MAIN is the main program for SANDIA_RULES_PRB.
!
!  Discussion:
!
!    SANDIA_RULES_PRB calls a set of tests for the SANDIA_RULES library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) r

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SANDIA_RULES_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SANDIA_RULES library.'

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
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test21 ( )
  call test22 ( )
!
!  TEST225 takes an input argument R, a rule index between 1 and 10.
!
  r = 1
  call test23 ( r )
  r = 3
  call test23 ( r )
  r = 4
  call test23 ( r )
  r = 11
  call test23 ( r )

  call test24 ( )
  call test25 ( )
  call test26 ( )
  call test27 ( )
  call test28 ( )
  call test285 ( )
  call test29 ( )

  call test30 ( )
  call test31 ( )
  call test32 ( )
  call test33 ( )
  call test34 ( )
  call test35 ( )
  call test36 ( )
  call test37 ( )
  call test38 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SANDIA_RULES_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CHEBYSHEV1_COMPUTE against CHEBYSHEV1_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) / sqrt ( 1 - x^2 ) dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CHEBYSHEV1_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = 2*ORDER-1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Estimate       Exact            Error'

  do order = 1, order_max

    write ( *, '(a)' ) ' '

    allocate ( f(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call chebyshev1_compute ( order, x, w )

    do n = 0, 2 * order + 2

      call chebyshev1_integral ( n, exact )

      if ( n == 0 ) then

        do i = 1, order
          f(i) = 1.0D+00
        end do

      else

        do i = 1, order
          f(i) = x(i)**n
        end do

      end if

      estimate = dot_product ( w(1:order), f(1:order) )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      order, n, estimate, exact, error

    end do

    deallocate ( f )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CHEBYSHEV1_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) / sqrt(1-x^2) dx.'

  do order = 1, order_max

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Order = ', order

    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call chebyshev1_compute ( order, x, w )

    do i = 1, order
      write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests CHEBYSHEV2_COMPUTE against CHEBYSHEV2_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) * sqrt ( 1 - x^2 ) dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CHEBYSHEV2_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = 2*ORDER-1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Estimate       Exact            Error'

  do order = 1, order_max

    write ( *, '(a)' ) ' '

    allocate ( f(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call chebyshev2_compute ( order, x, w )

    do n = 0, 2 * order + 2

      call chebyshev2_integral ( n, exact )

      if ( n == 0 ) then

        do i = 1, order
          f(i) = 1.0D+00
        end do

      else

        do i = 1, order
          f(i) = x(i)**n
        end do

      end if

      estimate = dot_product ( w(1:order), f(1:order) )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      order, n, estimate, exact, error

    end do

    deallocate ( f )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests CHEBYSHEV2_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) * sqrt(1-x^2) dx.'

  do order = 1, order_max

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Order = ', order

    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call chebyshev2_compute ( order, x, w )

    do i = 1, order
      write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests CLENSHAW_CURTIS_COMPUTE against LEGENDRE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_hi
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEGENDRE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '    N = ORDER+1 if ORDER is odd, or'
  write ( *, '(a)' ) '    N = ORDER   if ORDER is even.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Estimate       Exact            Error'

  do order = 1, order_max

    write ( *, '(a)' ) ' '

    allocate ( f(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call clenshaw_curtis_compute ( order, x, w )

    if ( mod ( order, 2 ) == 0 ) then
      n_hi = order + 2
    else
      n_hi = order + 3
    end if

    do n = 0, n_hi

      call legendre_integral ( n, exact )

      if ( n == 0 ) then

        f(1:order) = 1.0D+00

      else

        do i = 1, order
          f(i) = x(i)**n
        end do

      end if

      estimate = dot_product ( w(1:order), f(1:order) )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      order, n, estimate, exact, error

    end do

    deallocate ( f )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test06( )

!*****************************************************************************80
!
!! TEST06 tests CLENSHAW_CURTIS_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   15 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) dx.'

  do order = 1, order_max

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Order = ', order

    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call clenshaw_curtis_compute ( order, x, w )

    do i = 1, order
      write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests FEJER2_COMPUTE against LEGENDRE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_hi
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  FEJER2_COMPUTE computes a Fejer Type 2 rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= 1.0 ) f(x) dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEGENDRE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '    N = ORDER+1 if ORDER is odd, or'
  write ( *, '(a)' ) '    N = ORDER   if ORDER is even.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Estimate       Exact            Error'

  do order = 1, order_max

    write ( *, '(a)' ) ' '

    allocate ( f(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call fejer2_compute ( order, x, w )

    if ( mod ( order, 2 ) == 0 ) then
      n_hi = order + 2
    else
      n_hi = order + 3
    end if

    do n = 0, n_hi

      call legendre_integral ( n, exact )

      if ( n == 0 ) then

        f(1:order) = 1.0D+00

      else

        do i = 1, order
          f(i) = x(i)**n
        end do

      end if

      estimate = dot_product ( w(1:order), f(1:order) )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      order, n, estimate, exact, error

    end do

    deallocate ( f )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests FEJER2_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   15 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  FEJER2_COMPUTE computes a Fejer Type 2 rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) dx.'

  do order = 1, order_max

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Order = ', order

    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call fejer2_compute ( order, x, w )

    do i = 1, order
      write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests GEGENBAUER_COMPUTE against GEGENBAUER_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test =  (/ &
    0.5D+00, 1.0D+00, 2.5D+00 /)
  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  GEGENBAUER_COMPUTE computes a generalized Gauss-Gegenbauer rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( 0 <= x < +oo ) f(x) (1-x^2)^alpha dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GEGENBAUER_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = 2*ORDER-1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Alpha           Estimate       Exact            Error'

  do test = 1, test_num

    alpha = alpha_test(test)

    do order = 1, order_max

      write ( *, '(a)' ) ' '

      allocate ( f(1:order) )
      allocate ( w(1:order) )
      allocate ( x(1:order) )

      call gegenbauer_compute ( order, alpha, x, w )

      do n = 0, 2 * order + 2

        call gegenbauer_integral ( n, alpha, exact )

        if ( n == 0 ) then

          do i = 1, order
            f(i) = 1.0D+00
          end do

        else

          do i = 1, order
            f(i) = x(i)**n
          end do

        end if

        estimate = dot_product ( w(1:order), f(1:order) )

        error = abs ( exact - estimate )

        write ( *, '(2x,i8,2x,i8,2x,g10.4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
          order, n, alpha, estimate, exact, error

      end do

      deallocate ( f )
      deallocate ( w )
      deallocate ( x )

    end do
  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests GEGENBAUER_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test =  (/ &
    0.5D+00, 1.0D+00, 2.5D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  GEGENBAUER_COMPUTE computes a generalized Gauss-Gegenbauer rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) (1-x^2)^alpha dx.'

  do test = 1, test_num

    alpha = alpha_test(test)

    do order = 1, order_max

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Order = ', order
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha

      allocate ( w(1:order) )
      allocate ( x(1:order) )

      call gegenbauer_compute ( order, alpha, x, w )

      do i = 1, order
        write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
      end do

      deallocate ( w )
      deallocate ( x )

    end do

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests GEN_HERMITE_COMPUTE against GEN_HERMITE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test =  (/ &
    0.5D+00, 1.0D+00, 2.5D+00 /)
  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GEN_HERMITE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = 2*ORDER-1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Alpha           Estimate       Exact            Error'

  do test = 1, test_num

    alpha = alpha_test(test)

    do order = 1, order_max

      write ( *, '(a)' ) ' '

      allocate ( f(1:order) )
      allocate ( w(1:order) )
      allocate ( x(1:order) )

      call gen_hermite_compute ( order, alpha, x, w )

      do n = 0, 2 * order + 2

        call gen_hermite_integral ( n, alpha, exact )

        if ( n == 0 ) then

          do i = 1, order
            f(i) = 1.0D+00
          end do

        else

          do i = 1, order
            f(i) = x(i)**n
         end do

        end if

        estimate = dot_product ( w(1:order), f(1:order) )

        error = abs ( exact - estimate )

        write ( *, '(2x,i8,2x,i8,2x,g10.4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
          order, n, alpha, estimate, exact, error

      end do

      deallocate ( f )
      deallocate ( w )
      deallocate ( x )

    end do
  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests GEN_HERMITE_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test =  (/ &
    0.5D+00, 1.0D+00, 2.5D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.'

  do test = 1, test_num

    alpha = alpha_test(test)

    do order = 1, order_max

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Order = ', order
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha

      allocate ( w(1:order) )
      allocate ( x(1:order) )

      call gen_hermite_compute ( order, alpha, x, w )

      do i = 1, order
        write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
      end do

      deallocate ( w )
      deallocate ( x )

    end do
  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests GEN_LAGUERRE_COMPUTE against GEN_LAGUERRE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test =  (/ &
    0.5, 1.0, 2.5 /)
  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GEN_LAGUERRE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = 2*ORDER-1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Alpha           Estimate       Exact            Error'

  do test = 1, test_num

    alpha = alpha_test(test)

    do order = 1, order_max

      write ( *, '(a)' ) ' '

      allocate ( f(1:order) )
      allocate ( w(1:order) )
      allocate ( x(1:order) )

      call gen_laguerre_compute ( order, alpha, x, w )

      do n = 0, 2 * order + 2

        call gen_laguerre_integral ( n, alpha, exact )

        if ( n == 0 ) then

          do i = 1, order
            f(i) = 1.0D+00
          end do

        else

          do i = 1, order
            f(i) = x(i)**n
          end do

        end if

        estimate = dot_product ( w(1:order), f(1:order) )

        error = abs ( exact - estimate )

        write ( *, '(2x,i8,2x,i8,2x,g10.4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
          order, n, alpha, estimate, exact, error

      end do

      deallocate ( f )
      deallocate ( w )
      deallocate ( x )

    end do
  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests GEN_LAGUERRE_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test =  (/ &
    0.5D+00, 1.0D+00, 2.5D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.'

  do test = 1, test_num

    alpha = alpha_test(test)

    do order = 1, order_max

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Order = ', order
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha

      allocate ( w(1:order) )
      allocate ( x(1:order) )

      call gen_laguerre_compute ( order, alpha, x, w )

       do i = 1, order
        write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
      end do

      deallocate ( w )
      deallocate ( x )

    end do
  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests HERMITE_COMPUTE and HERMITE_SS_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) estimate1
  real ( kind = 8 ) estimate2
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  real ( kind = 8 ), allocatable, dimension ( : ) :: w1
  real ( kind = 8 ), allocatable, dimension ( : ) :: w2
  real ( kind = 8 ), allocatable, dimension ( : ) :: x1
  real ( kind = 8 ), allocatable, dimension ( : ) :: x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  HERMITE_COMPUTE and HERMITE_SS_COMPUTE'
  write ( *, '(a)' ) '  compute a Gauss-Hermite rule for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HERMITE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = 2*ORDER-1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Exact           Estimate1      Error1           Estimate2       Error2'

  do order = 1, order_max

    write ( *, '(a)' ) ' '

    allocate ( f(1:order) )
    allocate ( w1(1:order) )
    allocate ( w2(1:order) )
    allocate ( x1(1:order) )
    allocate ( x2(1:order) )

    call hermite_compute ( order, x1, w1 )
    call hermite_ss_compute ( order, x2, w2 )

    do n = 0, 2 * order + 2

      call hermite_integral ( n, exact )

      if ( n == 0 ) then
        f(1:order) = 1.0D+00
      else
        f(1:order) = x1(1:order)**n
      end if

      estimate1 = dot_product ( w1(1:order), f(1:order) )

      error1 = abs ( exact - estimate1 )

      if ( n == 0 ) then
        f(1:order) = 1.0D+00
      else
        f(1:order) = x2(1:order)**n
      end if

      estimate2 = dot_product ( w2(1:order), f(1:order) )

      error2 = abs ( exact - estimate2 )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        order, n, exact, estimate1, error1, estimate2, error2

    end do

    deallocate ( f )
    deallocate ( w1 )
    deallocate ( w2 )
    deallocate ( x1 )
    deallocate ( x2 )

  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests HERMITE_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  HERMITE_COMPUTE computes a Gauss-Hermite rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.'

  do order = 1, order_max

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Order = ', order

    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call hermite_compute ( order, x, w )

    do i = 1, order
      write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests JACOBI_COMPUTE against JACOBI_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test =  (/ &
    0.5D+00, 1.0D+00, 2.5D+00 /)
  real ( kind = 8 ) beta
  real ( kind = 8 ), dimension ( test_num ) :: beta_test =  (/ &
    0.5D+00, 1.0D+00, 2.5D+00 /)
  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test1
  integer ( kind = 4 ) test2
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  JACOBI_COMPUTE computes a Gauss-Jacobi rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JACOBI_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = 2*ORDER-1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Alpha           Beta            Estimate       Exact            Error'

  do test1 = 1, test_num

    alpha = alpha_test(test1)

    do test2 = 1, test_num

      beta = beta_test(test2)

      do order = 1, order_max

        write ( *, '(a)' ) ' '

        allocate ( f(1:order) )
        allocate ( w(1:order) )
        allocate ( x(1:order) )

        call jacobi_compute ( order, alpha, beta, x, w )

        do n = 0, 2 * order + 2

          call jacobi_integral ( n, alpha, beta, exact )

          if ( n == 0 ) then

            do i = 1, order
              f(i) = 1.0D+00
            end do

          else

            do i = 1, order
              f(i) = x(i)**n
            end do

          end if

          estimate = dot_product ( w(1:order), f(1:order) )

          error = abs ( exact - estimate )

         write ( *, '(2x,i8,2x,i8,2x,g10.4,2x,g10.4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
          order, n, alpha, beta, estimate, exact, error

        end do

        deallocate ( f )
        deallocate ( w )
        deallocate ( x )

      end do
    end do
  end do

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests JACOBI_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test =  (/ &
    0.5D+00, 1.0D+00, 2.5D+00 /)
  real ( kind = 8 ) beta
  real ( kind = 8 ), dimension ( test_num ) :: beta_test =  (/ &
    0.5D+00, 1.0D+00, 2.5D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test1
  integer ( kind = 4 ) test2
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  JACOBI_COMPUTE computes a Gauss-Jacobi rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.'

  do test1 = 1, test_num

    alpha = alpha_test(test1)

    do test2 = 1, test_num

      beta = beta_test(test2)

      do order = 1, order_max

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Order = ', order
        write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
        write ( *, '(a,g14.6)' ) '  BETA =  ', beta

        allocate ( w(1:order) )
        allocate ( x(1:order) )

        call jacobi_compute ( order, alpha, beta, x, w )

        do i = 1, order
          write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
        end do

        deallocate ( w )
        deallocate ( x )

      end do
    end do
  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests LAGUERRE_COMPUTE against LAGUERRE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  LAGUERRE_COMPUTE computes a Gauss-Laguerre rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LAGUERRE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = 2*ORDER-1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Estimate       Exact            Error'

  do order = 1, order_max

    write ( *, '(a)' ) ' '

    allocate ( f(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call laguerre_compute ( order, x, w )

    do n = 0, 2 * order + 2

      call laguerre_integral ( n, exact )

      if ( n == 0 ) then

        do i = 1, order
          f(i) = 1.0D+00
        end do

      else

        do i = 1, order
          f(i) = x(i)**n
        end do

      end if

      estimate = dot_product ( w(1:order), f(1:order) )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        order, n, estimate, exact, error

    end do

    deallocate ( f )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests LAGUERRE_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.'

  do order = 1, order_max

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Order = ', order

    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call laguerre_compute ( order, x, w )

    do i = 1, order
      write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
    end do

    deallocate ( w )
    deallocate ( x )
  end do

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests LEGENDRE_COMPUTE against LEGENDRE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  LEGENDRE_COMPUTE computes a Gauss-Legendre rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x)  dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEGENDRE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = 2*ORDER-1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Estimate       Exact            Error'

  do order = 1, order_max

    write ( *, '(a)' ) ' '

    allocate ( f(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call legendre_compute ( order, x, w )

    do n = 0, 2 * order + 2

      call legendre_integral ( n, exact )

      if ( n == 0 ) then

        do i = 1, order
          f(i) = 1.0D+00
        end do

      else

        do i = 1, order
          f(i) = x(i)**n
        end do

      end if

      estimate = dot_product ( w(1:order), f(1:order) )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        order, n, estimate, exact, error

    end do

    deallocate ( f )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests LEGENDRE_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: order_max = 10
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  LEGENDRE_COMPUTE computes a Gauss-Legendre rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x) dx.'

  do order = 1, order_max

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Order = ', order

    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call legendre_compute ( order, x, w )

    do i = 1, order
      write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) i, x(i), w(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test23 ( r )

!*****************************************************************************80
!
!! TEST23 tests LEVEL_GROWTH_TO_ORDER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 11

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) g
  integer ( kind = 4 ) growth(dim_num)
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rule(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  LEVEL_GROWTH_TO_ORDER uses rule, level and growth'
  write ( *, '(a)' ) '  to determine the order of a 1D rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we examine rule ', r
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     LEVEL   0     1     2     3     4     5     6     7     8     9    10'
  write ( *, '(a)' ) 'GROWTH'

  rule(1:dim_num) = r

  do g = 0, 6

    if ( r == 3 .or. r == 10 ) then
      if ( g == 1 .or. g == 2 .or. g == 3 ) then
        cycle
      end if
    end if

    growth(1:dim_num) = g

    do dim = 1, dim_num
      level(dim) = dim - 1
    end do

    call level_growth_to_order ( dim_num, level, rule, growth, order )

    write ( *, '(2x,i4,2x,11(2x,i4))' ) g, order(1:dim_num)

  end do

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests LEVEL_TO_ORDER_DEFAULT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 11

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rule(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  LEVEL_TO_ORDER_DEFAULT uses a default rule to'
  write ( *, '(a)' ) '  determine the order of a rule from its level.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10'
  write ( *, '(a)' ) ' '

  do r = 1, 16

    rule(1:dim_num) = r

    do dim = 1, dim_num
      level(dim) = dim - 1
    end do

    call level_to_order_default ( dim_num, level, rule, order )

    write ( *, '(2x,i4,2x,11(2x,i4))' ) r, order(1:dim_num)

  end do

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests LEVEL_TO_ORDER_EXPONENTIAL
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 11

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rule(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  LEVEL_TO_ORDER_EXPONENTIAL uses an exponential rule to'
  write ( *, '(a)' ) '  determine the order of a rule from its level.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10'
  write ( *, '(a)' ) ' '

  do r = 1, 10

    rule(1:dim_num) = r

    do dim = 1, dim_num
      level(dim) = dim - 1
    end do

    call level_to_order_exponential ( dim_num, level, rule, order )

    write ( *, '(2x,i4,2x,11i6)' ) r, order(1:dim_num)

  end do

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests LEVEL_TO_ORDER_EXPONENTIAL_SLOW
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 11

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rule(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  LEVEL_TO_ORDER_EXPONENTIAL_SLOW uses a slow exponential'
  write ( *, '(a)' ) '  rule to determine the order of a rule from its level.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Since it is really only useful for fully nested rules,'
  write ( *, '(a)' ) '  we only consider rules 11, 12 and 13.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10'
  write ( *, '(a)' ) ' '

  do r = 11, 13

    rule(1:dim_num) = r

    do dim = 1, dim_num
      level(dim) = dim - 1
    end do

    call level_to_order_exponential_slow ( dim_num, level, rule, order )

    write ( *, '(2x,i4,2x,11i6)' ) r, order(1:dim_num)

  end do

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests LEVEL_TO_ORDER_LINEAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 11

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rule(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  LEVEL_TO_ORDER_LINEAR uses a linear rule to'
  write ( *, '(a)' ) '  determine the order of a rule from its level.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10'
  write ( *, '(a)' ) ' '

  do r = 1, 10

    rule(1:dim_num) = r

    do dim = 1, dim_num
      level(dim) = dim - 1
    end do

    call level_to_order_linear ( dim_num, level, rule, order )

    write ( *, '(2x,i4,2x,11i6)' ) r, order(1:dim_num)

  end do

  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests PATTERSON_LOOKUP against LEGENDRE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) level
  integer ( kind = 4 ) :: level_max = 5
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) p
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  PATTERSON_LOOKUP looks up a Gauss-Patterson rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    Integral ( -1 <= x <= +1 ) f(x)  dx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEGENDRE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^n.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A rule of order ORDER should be exact for monomials X^N'
  write ( *, '(a)' ) '  up to N = (3*ORDER+1)/2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, for each order, the LAST THREE estimates'
  write ( *, '(a)' ) '  are made on monomials that exceed the exactness limit for the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order         N       Estimate       Exact            Error'

  do level = 0, level_max

    order = 2**( level + 1 ) - 1

    write ( *, '(a)' ) ' '

    allocate ( f(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:order) )

    call patterson_lookup ( order, x, w )

    if ( level == 0 ) then
      p = 1
    else
      p = ( 3 * order + 1 ) / 2
    end if

    do n = 0, p + 3

      call legendre_integral ( n, exact )

      if ( n == 0 ) then

        do i = 1, order
          f(i) = 1.0D+00
        end do

      else

        do i = 1, order
          f(i) = x(i)**n
        end do

      end if

      estimate = dot_product ( w(1:order), f(1:order) )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        order, n, estimate, exact, error

    end do

    deallocate ( f )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test285 ( )

!*****************************************************************************80
!
!  Purpose:
!
!    TEST285 tests POINT_RADIAL_TOL_UNIQUE_INDEX_INC1, INC2 and INC3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    0
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n1 = 11
  integer ( kind = 4 ), parameter :: n2 = 8

  real ( kind = 8 ) :: a1(m,n1) = reshape ( (/ &
    0.0, 0.0, &
    0.5, 0.5, &
    1.0, 0.0, &
    0.0, 1.0, &
    1.0, 1.0, &
    0.500000001, 0.5, &
    0.0, 0.0, &
    0.0, 0.5, &
    0.5, 0.0, &
    1.0, 0.5, &
    0.5, 1.0 /), (/ m, n1 /) )
  real ( kind = 8 ) :: a2(m,n2) = reshape ( (/ &
    0.4999999999, 0.5, &
    0.75,        0.25, &
    0.500000001, 0.9999999999, &
    0.500000001, 0.0000000001, &
    0.25,        0.75, &
    0.75,        0.25, &
    0.250000001, 0.7499999999, &
    0.75,        0.75 /), (/ m, n2 /) )
  real ( kind = 8 ), allocatable :: a3(:,:)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ), allocatable :: indx1(:)
  integer ( kind = 4 ), allocatable :: indx2(:)
  integer ( kind = 4 ), allocatable :: indx3(:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n3
  real ( kind = 8 ), allocatable :: r1(:)
  real ( kind = 8 ), allocatable :: r2(:)
  real ( kind = 8 ), allocatable :: r3(:)
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx_value
  integer ( kind = 4 ), allocatable :: undx1(:)
  integer ( kind = 4 ), allocatable :: undx2(:)
  integer ( kind = 4 ), allocatable :: undx3(:)
  integer ( kind = 4 ) unique_num1
  integer ( kind = 4 ) unique_num2
  integer ( kind = 4 ) unique_num3
  logical, allocatable :: unique1(:)
  logical, allocatable :: unique2(:)
  logical, allocatable :: unique3(:)
  integer ( kind = 4 ), allocatable :: xdnu1(:)
  integer ( kind = 4 ), allocatable :: xdnu2(:)
  integer ( kind = 4 ), allocatable :: xdnu3(:)
  real ( kind = 8 ), allocatable :: z(:)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST285'
  write ( *, '(a)' ) '  POINT_RADIAL_TOL_UNIQUE_INDEX_INC1 can index unique'
  write ( *, '(a)' ) '  points in a "permanent" point set'
  write ( *, '(a)' ) '  POINT_RADIAL_TOL_UNIQUE_INDEX_INC2 can incremented by'
  write ( *, '(a)' ) '  "temporary" points.'
  write ( *, '(a)' ) '  POINT_RADIAL_TOL_UNIQUE_INDEX_INC3 can merge permanent'
  write ( *, '(a)' ) '  and temporary points.'

  tol = sqrt ( r8_epsilon ( ) )
  write ( *, '(a,g14.6)' ) '  Using tolerance TOL = ', tol
!
!  Step 1
!
  allocate ( indx1(n1) )
  allocate ( r1(n1) )
  allocate ( undx1(n1) )
  allocate ( unique1(n1) )
  allocate ( xdnu1(n1) )
  allocate ( z(m) )

  call point_radial_tol_unique_index_inc1 ( m, n1, a1, tol, seed, z, r1, &
    indx1, unique1, unique_num1, undx1, xdnu1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  UNIQUE_NUM1 = ', unique_num1
  write ( *, '(a,i8)' ) '  Expected   =  ', 9

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Item I1, unique index XDNU1(I1), representative location UNDX1(XDNU1(I1)):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           I1  XDNU1  UNDX1'
  write ( *, '(a)' ) ' '
  do i1 = 1, n1
    write ( *, '(2x,4x,2x,i4,x,i4,2x,i4)' ), i1, xdnu1(i1), undx1(xdnu1(i1))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unique item I1, location UNDX1(I1), value A1(:,UNDX1(I1)):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          I1 UNDX1  --A1(1,*)---  --A1(2,*)---'
  write ( *, '(a)' ) ' '
  do i1 = 1, unique_num1
    write ( *, '(2x,4x,2x,i4,2x,i4,2x,f12.4,2x,f12.4)' ) &
      i1, undx1(i1), a1(1,undx1(i1)), a1(2,undx1(i1))
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unique item I1, location UNDX1(I1), value A1(:,UNDX1(I1)):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          I1 UNIQUE1       R1     --A1(1,I1)--  --A1(2,I1)--'
  write ( *, '(a)' ) ' '
  do i1 = 1, n1
    write ( *, '(2x,4x,2x,i4,2x,l4,2x,f12.4,2x,f12.4,2x,f12.4)' ) &
      i1, unique1(i1), r1(i1), a1(1,i1), a1(2,i1)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          I1   INDX1   R1(INDX1)'
  write ( *, '(a)' ) ' '
  do i1 = 1, n1
    write ( *, '(2x,4x,2x,i4,2x,i4,2x,f12.4)' ) i1, indx1(i1), r1(indx1(i1))
  end do
!
!  Step 2
!
  allocate ( indx2(n2) )
  allocate ( r2(n2) )
  allocate ( undx2(n2) )
  allocate ( unique2(n2) )
  allocate ( xdnu2(n2) )

  call point_radial_tol_unique_index_inc2 ( m, n1, a1, n2, a2, tol, z, &
    r1, indx1, unique1, unique_num1, undx1, xdnu1, &
    r2, indx2, unique2, unique_num2, undx2, xdnu2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  UNIQUE_NUM2 = ', unique_num2
  write ( *, '(a,i8)' ) '  Expected   =  ', 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Item I2, unique index XDNU2(I2), representative location UNDX2(XDNU2(I2)):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I2 XDNU2 UNDX2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Temporary data)'
  write ( *, '(a)' ) ' '
  do i2 = 1, n2
    if ( xdnu2(i2) < unique_num1 ) then
      undx_value = undx1(xdnu2(i2))
    else
      undx_value = undx2( xdnu2(i2) - unique_num1 )
    end if
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) &
      i2, i2 + n1, xdnu2(i2), undx_value
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unique item I2, location UNDX2(I2), value A2(:,UNDX2(I2)):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I2 UNDX2  --A2(1,*)---  --A2(2,*)---'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Temporary data)'
  write ( *, '(a)' ) ' '
  do i2 = 1, unique_num2
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,f12.4,2x,f12.4)' ) &
      i2, i2 + unique_num1, undx2(i2), a2(1,undx2(i2)-n1), a2(2,undx2(i2)-n1)
  end do
!
!  Step 3.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Merge the temporary data with the permanent data.'

  n3 = n1 + n2

  allocate ( a3(m,n3) )
  allocate ( indx3(n3) )
  allocate ( r3(n3) )
  allocate ( undx3(n3) )
  allocate ( unique3(n3) )
  allocate ( xdnu3(n3) )

  call point_radial_tol_unique_index_inc3 ( m, &
     n1, a1, r1, indx1, unique1, unique_num1, undx1, xdnu1, &
     n2, a2, r2, indx2, unique2, unique_num2, undx2, xdnu2, &
     n3, a3, r3, indx3, unique3, unique_num3, undx3, xdnu3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  UNIQUE_NUM3 = ', unique_num3
  write ( *, '(a)' ) '  Expected   =  12'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Item I3, unique index XDNU3(I3), representative location UNDX3(XDNU3(I3)):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          I3 XDNU3 UNDX3'
  write ( *, '(a)' ) ' '
  do i3 = 1, n3
    write ( *, '(2x,4x,2x,i4,2x,i4,2x,i4)' ) &
      i3, xdnu3(i3), undx3(xdnu3(i3))
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unique item I3, location UNDX3(I3), value A3(:,UNDX3(I3)):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          I3 UNDX3  --A3(1,*)---  --A3(2,*)---'
  write ( *, '(a)' ) ' '
  do i3 = 1, unique_num3
    write ( *, '(2x,4x,2x,i4,2x,i4,2x,f12.4,2x,f12.4)' ) &
      i3, undx3(i3), a3(1,undx3(i3)), a3(2,undx3(i3))
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unique item I3, location UNDX3(I3), value A3(:,UNDX3(I3)):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          I3 UNIQUE3       R3     --A3(1,I3)--  --A3(2,I3)--'
  write ( *, '(a)' ) ' '
  do i3 = 1, n3
    write ( *, '(2x,4x,2x,i4,2x,l4,2x,f12.4,2x,f12.4,2x,f12.4)' ) &
      i3, unique3(i3), r3(i3), a3(1,i3), a3(2,i3)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          I3   INDX3   R3(INDX3)'
  write ( *, '(a)' ) ' '
  do i3 = 1, n3
    write ( *, '(2x,4x,2x,i4,2x,i4,2x,f12.4)' ) &
      i3, indx3(i3), r3(indx3(i3))
  end do

  deallocate ( a3 )
  deallocate ( indx1 )
  deallocate ( indx2 )
  deallocate ( indx3 )
  deallocate ( r1 )
  deallocate ( r2 )
  deallocate ( r3 )
  deallocate ( undx1 )
  deallocate ( undx2 )
  deallocate ( undx3 )
  deallocate ( unique1 )
  deallocate ( unique2 )
  deallocate ( unique3 )
  deallocate ( xdnu1 )
  deallocate ( xdnu2 )
  deallocate ( xdnu3 )
  deallocate ( z )

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests R8COL_TOL_UNDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 22

  real ( kind = 8 ) :: a(m,n) = reshape ( (/ &
    1.9,  0.0, 10.0, &
    2.0,  6.0, 10.0, &
    4.0,  8.0, 12.0, &
    1.0,  5.0,  9.0, &
    3.0,  7.0, 11.0, &
    2.0,  6.0,  0.0, &
    2.0,  0.0, 10.1, &
    2.0,  0.1, 10.0, &
    3.0,  4.0, 18.0, &
    1.9,  8.0, 10.0, &
    0.0,  0.0,  0.0, &
    0.0,  6.0, 10.0, &
    2.1,  0.0, 10.0, &
    2.0,  6.0, 10.0, &
    3.0,  7.0, 11.0, &
    2.0,  0.0, 10.0, &
    2.0,  0.0, 10.0, &
    2.0,  6.0, 10.0, &
    1.0,  5.0,  9.0, &
    2.0,  0.0, 10.1, &
    1.0,  5.0,  9.1, &
    1.0,  5.1,  9.0 /), (/ m, n /) )
  real ( kind = 8 ), allocatable :: au(:,:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n_unique
  real ( kind = 8 ) tol
  integer ( kind = 4 ), allocatable :: undx(:)
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ), allocatable :: xdnu(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  R8COL_TOL_UNDEX produces index vectors which create a sorted'
  write ( *, '(a)' ) '  list of the tolerably unique columns of an R8COL,'
  write ( *, '(a)' ) '  and a map from the original R8COL to the (implicit)'
  write ( *, '(a)' ) '  R8COL of sorted tolerably unique elements.'

  call r8mat_transpose_print ( m, n, a, '  The unsorted R8COL (transposed):' )

  tol = 0.25

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Using tolerance = ', tol

  call r8col_tol_unique_count ( m, n, a, tol, n_unique )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of tolerably unique columns is ', n_unique

  allocate ( au(m,n_unique) )
  allocate ( undx(n_unique) )
  allocate ( xdnu(n) )

  call r8col_tol_undex ( m, n, a, n_unique, tol, undx, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XDNU points to the representative for each item.'
  write ( *, '(a)' ) '  UNDX selects the representatives.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  XDNU  UNDX'
  write ( *, '(a)' ) ' '
  do i = 1, n_unique
    write ( *, '(2x,i4,2x,i4,2x,i4)' ) i, xdnu(i), undx(i)
  end do
  do i = n_unique + 1, n
    write ( *, '(2x,i4,2x,i4)' ) i, xdnu(i)
  end do

  do j = 1, n_unique
    au(1:m,j) = a(1:m,undx(j))
  end do

  call r8mat_transpose_print ( m, n_unique, au, &
    '  The tolerably unique R8COL (transposed):' )

  deallocate ( au )
  deallocate ( undx )
  deallocate ( xdnu )

  return
end
subroutine test30 ( )

!****************************************************************************80
!
!  Purpose:
!
!    TEST30 tests R8VEC_SORT_HEAP_INDEX_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  R8VEC_SORT_HEAP_INDEX_A creates an ascending'
  write ( *, '(a)' ) '  sort index for a R8VEC.'

  seed = 123456789

  call r8vec_uniform_01 ( n, seed, a )

  call r8vec_print ( n, a, '  The unsorted array:' )

  call r8vec_sort_heap_index_a ( n, a, indx )

  call i4vec_print ( n, indx, '  The index vector:' )

  b(1:n) = a(indx(1:n))

  call r8vec_print ( n, b, '  The sorted array A(INDX(:)):' )

  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests HCE_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 11
  integer ( kind = 4 ), parameter :: n = 2 * m

  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31:'
  write ( *, '(a)' ) '  HCE_COMPUTE returns a quadrature rule'
  write ( *, '(a)' ) '  for piecewise Hermite cubic splines which are based'
  write ( *, '(a)' ) '  on equally spaced function and derivative data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we compute a rule of order N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I        X(I)        W(I)'
  write ( *, '(a)' ) ' '
  call hce_compute ( n, x, w )

  do i = 1, n
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, x(i), w(i)
  end do

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests HCC_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 11
  integer ( kind = 4 ), parameter :: n = 2 * m

  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32:'
  write ( *, '(a)' ) '  HCC_COMPUTE returns a quadrature rule'
  write ( *, '(a)' ) '  for piecewise Hermite cubic splines which are based'
  write ( *, '(a)' ) '  on Chebyshev-spaced function and derivative data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we compute a rule of order N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I        X(I)        W(I)'
  write ( *, '(a)' ) ' '
  call hcc_compute ( n, x, w )

  do i = 1, n
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, x(i), w(i)
  end do

  return
end
subroutine test33 ( )

!*****************************************************************************80
!
!! TEST33 tests HC_COMPUTE_WEIGHTS_FROM_POINTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ) dn(n)
  real ( kind = 8 ) fn(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) q
  real ( kind = 8 ) q_exact
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  integer ( kind = 4 ) test
  real ( kind = 8 ) w(2,n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33:'
  write ( *, '(a)' ) '  HC_COMPUTE_WEIGHTS_FROM_POINTS returns quadrature'
  write ( *, '(a)' ) '  weights for a Hermite cubic spline given the points.'

  seed = 123456789

  do test = 1, 3

    call r8vec_uniform_01 ( n, seed, r )

    x(1) = r(1)
    do j = 2, n
      x(j) = x(j-1) + r(j)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i1)' ) '  Test #', test
    write ( *, '(a)' ) '  Random spacing'
    write ( *, '(a,i8)' ) '  Number of points N = ', n
    write ( *, '(a,g14.6,a,g14.6,a)' ) '  Interval = [', x(1), ',', x(n), ']'

    call hc_compute_weights_from_points ( n, x, w )

    do j = 1, n
      call cubic_value ( x(j), fn(j), dn(j), s, t )
    end do

    q = dot_product ( w(1,1:n), fn(1:n) ) + dot_product ( w(2,1:n), dn(1:n) )

    call cubic_integrate ( x(1), x(n), q_exact )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Q         = ', q
    write ( *, '(a,g14.6)' ) '  Q (exact) = ', q_exact

  end do

  return
end
subroutine test34 ( )

!*****************************************************************************80
!
!! TEST34 uses HERMITE_INTERPOLANT on the Runge function at equally spaced points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) max_dif
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ndp
  integer ( kind = 4 ) ns
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xdp(:)
  real ( kind = 8 ), allocatable :: xs(:)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xt
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: ydp(:)
  real ( kind = 8 ), allocatable :: yp(:)
  real ( kind = 8 ), allocatable :: ys(:)
  real ( kind = 8 ), allocatable :: ysp(:)
  real ( kind = 8 ) yt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  HERMITE_INTERPOLANT computes the Hermite interpolant to data.'
  write ( *, '(a)' ) '  Here, f(x) is the Runge function'
  write ( *, '(a)' ) '  and the data is evaluated at equally spaced points.'
  write ( *, '(a)' ) '  As N increases, the maximum error grows.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     Max | F(X) - H(F(X)) |'
  write ( *, '(a)' ) ' '

  do n = 3, 15, 2

    allocate ( x(1:n) )
    allocate ( y(1:n) )
    allocate ( yp(1:n) )

    nd = 2 * n

    allocate ( xd(1:nd) )
    allocate ( yd(1:nd) )

    ndp = 2 * n - 1

    allocate ( xdp(1:ndp) )
    allocate ( ydp(1:ndp) )

    ns = 10 * ( n - 1 ) + 1
    allocate ( xs(1:ns) )
    allocate ( ys(1:ns) )
    allocate ( ysp(1:ns) )

    xlo = -5.0D+00
    xhi = +5.0D+00
    call r8vec_linspace ( n, xlo, xhi, x )

    y(1:n) = 1.0D+00 / ( 1.0D+00 + x(1:n)**2 )
    yp(1:n) = - 2.0D+00 * x(1:n) / ( 1.0D+00 + x(1:n)**2 )**2

    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
!
!  Compare exact and interpolant at sample points.
!
    call r8vec_linspace ( ns, xlo, xhi, xs )

    call hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, ys, ysp )

    max_dif = 0.0D+00
    do i = 1, ns
      xt = xs(i)
      yt = 1.0D+00 / ( 1.0D+00 + xt * xt )
      max_dif = max ( max_dif, abs ( ys(i) - yt ) )
    end do

    write ( *, '(2x,i4,2x,g14.6)' ) n, max_dif

    deallocate ( x )
    deallocate ( xd )
    deallocate ( xdp )
    deallocate ( xs )
    deallocate ( y )
    deallocate ( yd )
    deallocate ( ydp )
    deallocate ( yp )
    deallocate ( ys )
    deallocate ( ysp )

  end do

  return
end
subroutine test35 ( )

!*****************************************************************************80
!
!! TEST35 uses HERMITE_INTERPOLANT on the Runge function at Chebyshev points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) max_dif
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ndp
  integer ( kind = 4 ) ns
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xdp(:)
  real ( kind = 8 ), allocatable :: xs(:)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xt
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: ydp(:)
  real ( kind = 8 ), allocatable :: yp(:)
  real ( kind = 8 ), allocatable :: ys(:)
  real ( kind = 8 ), allocatable :: ysp(:)
  real ( kind = 8 ) yt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  HERMITE_INTERPOLANT computes the Hermite interpolant to data.'
  write ( *, '(a)' ) '  Here, f(x) is the Runge function'
  write ( *, '(a)' ) '  and the data is evaluated at Chebyshev points.'
  write ( *, '(a)' ) '  As N increases, the maximum error goes down.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     Max | F(X) - H(F(X)) |'
  write ( *, '(a)' ) ' '

  do n = 3, 15, 2

    allocate ( x(1:n) )
    allocate ( y(1:n) )
    allocate ( yp(1:n) )

    nd = 2 * n
    allocate ( xd(1:nd) )
    allocate ( yd(1:nd) )

    ndp = 2 * n - 1
    allocate ( xdp(1:ndp) )
    allocate ( ydp(1:ndp) )

    ns = 10 * ( n - 1 ) + 1
    allocate ( xs(1:ns) )
    allocate ( ys(1:ns) )
    allocate ( ysp(1:ns) )

    xlo = -5.0D+00
    xhi = +5.0D+00
    call r8vec_chebyshev ( n, xlo, xhi, x )

    y(1:n) = 1.0D+00 / ( 1.0D+00 + x(1:n)**2 )
    yp(1:n) = - 2.0D+00 * x(1:n) / ( 1.0D+00 + x(1:n)**2 )**2

    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
!
!  Compare exact and interpolant at sample points.
!
    call r8vec_linspace ( ns, xlo, xhi, xs )

    call hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, ys, ysp )

    max_dif = 0.0D+00
    do i = 1, ns
      xt = xs(i)
      yt = 1.0D+00 / ( 1.0D+00 + xt * xt )
      max_dif = max ( max_dif, abs ( ys(i) - yt ) )
    end do

    write ( *, '(2x,i4,2x,g14.6)' ) n, max_dif

    deallocate ( x )
    deallocate ( xd )
    deallocate ( xdp )
    deallocate ( xs )
    deallocate ( y )
    deallocate ( yd )
    deallocate ( ydp )
    deallocate ( yp )
    deallocate ( ys )
    deallocate ( ysp )

  end do

  return
end
subroutine test36 ( )

!*****************************************************************************80
!
!! TEST36 tests HERMITE_INTERPOLANT_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) e
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) q
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36:'
  write ( *, '(a)' ) '  HERMITE_INTERPOLANT_RULE'
  write ( *, '(a)' ) '  is given a set of N abscissas for a Hermite interpolant'
  write ( *, '(a)' ) '  and returns N pairs of quadrature weights'
  write ( *, '(a)' ) '  for function and derivative values at the abscissas.'
!
!  1: Behavior with increasing N.
!
  a = 0.0D+00
  b = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Observe behavior of quadrature weights for increasing N'
  write ( *, '(a,g14.6,a,g14.6)' ) '  We are working in ', a, ' <= X <= ', b

  do n = 3, 11, 2

    allocate ( x(1:n) )
    allocate ( w(1:2*n) )

    call r8vec_linspace ( n, a, b, x )
    call hermite_interpolant_rule ( n, a, b, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I       X               W(F(X))        W(F''(X))'
    write ( *, '(a)' ) ' '
    k = 1
    do i = 1, n
      write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w(k), w(k+1)
      k = k + 2
    end do

    deallocate ( x )
    deallocate ( w )

  end do
!
!  2: Integral estimates with equally spaced points.
!
  a = -5.0D+00
  b = 5.0D+00
  n = 11

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a)' ) '  Use the rule with N = ', n, ' to estimate integrals.'
  write ( *, '(a)' ) '  Points are equally spaced.'
  write ( *, '(a,g14.6,a,g14.6)' ) '  We are working in ', a, ' <= X <= ', b

  allocate ( x(1:n) )
  allocate ( w(1:2*n) )
  call r8vec_linspace ( n, a, b, x )
  call hermite_interpolant_rule ( n, a, b, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X               W(F(X))        W(F''(X))'
  write ( *, '(a)' ) ' '
  k = 1
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w(k), w(k+1)
    k = k + 2
  end do

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * 1 + w(k+1) * 0.0D+00
    k = k + 2
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i) + w(k+1) * 1.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X^2 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) / ( 1.0D+00 + x(i)**2 ) - w(k+1) * 2.0D+00 * x(i) / ( 1.0D+00 + x(i)**2 )**2
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1/(1+x^2) = ', q

  deallocate ( w )
  deallocate ( x )
!
!  3: Integral estimates with Chebyshev spaced points.
!
  a = -5.0D+00
  b = 5.0D+00
  n = 11

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a)' ) '  Use the rule with N = ', n, ' to estimate integrals.'
  write ( *, '(a)' ) '  Points are Chebyshev spaced.'
  write ( *, '(a,g14.6,a,g14.6)' ) '  We are working in ', a, ' <= X <= ', b

  allocate ( x(1:n) )
  allocate ( w(1:2*n) )
  call r8vec_chebyshev ( n, a, b, x )
  call hermite_interpolant_rule ( n, a, b, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X               W(F(X))        W(F''(X))'
  write ( *, '(a)' ) ' '
  k = 1
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w(k), w(k+1)
    k = k + 2
  end do

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * 1 + w(k+1) * 0.0D+00
    k = k + 2
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i) + w(k+1) * 1.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X^2 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) / ( 1.0D+00 + x(i)**2 ) - w(k+1) * 2.0D+00 * x(i) / ( 1.0D+00 + x(i)**2 )**2
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1/(1+x^2) = ', q

  deallocate ( w )
  deallocate ( x )
!
!  4: Integral estimates with Legendre spaced points.
!
  a = -5.0D+00
  b = 5.0D+00
  n = 11

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a)' ) '  Use the rule with N = ', n, ' to estimate integrals.'
  write ( *, '(a)' ) '  Points are Legendre spaced.'
  write ( *, '(a,g14.6,a,g14.6)' ) '  We are working in ', a, ' <= X <= ', b

  allocate ( x(1:n) )
  allocate ( w(1:2*n) )
  call r8vec_legendre ( n, a, b, x )
  call hermite_interpolant_rule ( n, a, b, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X               W(F(X))        W(F''(X))'
  write ( *, '(a)' ) ' '
  k = 1
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w(k), w(k+1)
    k = k + 2
  end do

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * 1 + w(k+1) * 0.0D+00
    k = k + 2
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i) + w(k+1) * 1.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X^2 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) / ( 1.0D+00 + x(i)**2 ) - w(k+1) * 2.0D+00 * x(i) / ( 1.0D+00 + x(i)**2 )**2
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1/(1+x^2) = ', q

  deallocate ( w )
  deallocate ( x )
!
!  5: Integral estimates with Chebyshev spaced points on 1/(1+x^2), increasing N.
!
  a = -5.0D+00
  b = 5.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a)' ) '  Approximate integral of 1/(1+x^2) with increasing N.'
  write ( *, '(a)' ) '  Points are Chebyshev spaced.'
  write ( *, '(a,g14.6,a,g14.6)' ) '  We are working in ', a, ' <= X <= ', b
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     Estimate         Error'
  write ( *, '(a)' ) ' '

  do n = 2, 11

    allocate ( x(1:n) )
    allocate ( w(1:2*n) )
    call r8vec_chebyshev ( n, a, b, x )
    call hermite_interpolant_rule ( n, a, b, x, w )
    q = 0.0D+00
    k = 1
    do i = 1, n
      q = q + w(k) / ( 1.0D+00 + x(i)**2 ) &
            - w(k+1) * 2.0D+00 * x(i) / ( 1.0D+00 + x(i)**2 )**2
      k = k + 2
    end do

    exact = atan ( b ) - atan ( a )
    error = abs ( q - exact )
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) n, q, error

    deallocate ( w )
    deallocate ( x )

  end do
!
!  6: Integral estimates, with Chebyshev spaced points, for monomials, using N = 11.
!
  a = -1.0D+00
  b = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a)' ) '  Compute integral estimates for monomials X^0 through X^15.'
  write ( *, '(a)' ) '  Use N = 5, 9, 13, 17, 21 point rules.'
  write ( *, '(a)' ) '  Points are Chebyshev spaced.'
  write ( *, '(a,g14.6,a,g14.6)' ) '  We are working in ', a, ' <= X <= ', b

  do n = 5, 21, 4

    allocate ( x(1:n) )
    allocate ( w(1:2*n) )
    call r8vec_chebyshev ( n, a, b, x )
    call hermite_interpolant_rule ( n, a, b, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i2)' ) '  Estimates are made using N = ', n
    write ( *, '(a)' ) '  F(X)         Integral        Estimate           Error'
    write ( *, '(a)' ) ' '

    do e = 0, 15
      q = 0.0D+00
      k = 1
      do i = 1, n
        if ( e == 0 ) then
          q = q + w(k)
        else
          q = q + w(k) * x(i)**e + w(k+1) * e * x(i)**(e-1)
        end if
        k = k + 2
      end do
      exact = ( b**(e+1) - a**(e+1) ) / real ( e + 1, kind = 8 )
      write ( *, '(a,i7,2x,g14.6,2x,g14.6,2x,g14.6)' ) '  X^', e, exact, q, abs ( exact - q )
    end do

    q = 0.0D+00
    k = 1
    do i = 1, n
      q = q + w(k) / ( 1.0D+00 + x(i)**2 ) - w(k+1) * 2.0D+00 * x(i) / ( 1.0D+00 + x(i)**2 )**2
      k = k + 2
    end do
    exact = atan ( b ) - atan ( a )
    write ( *, '(a,2x,g14.6,2x,g14.6,2x,g14.6)' ) '  1/(1+x^2)', exact, q, abs ( exact - q )

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test37 ( )

!*****************************************************************************80
!
!! TEST37 checks that the HGK weights are correctly scaled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) o
  integer ( kind = 4 ), dimension ( 8 ) :: order = (/ 1, 3, 9, 19, 35, 37, 41, 43 /)
  real ( kind = 8 ), parameter :: sqrtpi = 1.7724538509055159D+00
  integer ( kind = 4 ) rule;
  real ( kind = 8 ) s
  real ( kind = 8 ), allocatable :: w(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS looks up weights'
  write ( *, '(a)' ) '  for Genz-Keister quadrature rules for the Hermite weight function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This test simply checks that, for each rule, the quadrature'
  write ( *, '(a)' ) '  weights correctly sum to sqrt(pi).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Index     Order         Sum of Weights'
  write ( *, '(a)' ) ' '

  do rule = 1, 8

    o = order(rule)

    allocate ( w(1:o) )

    call hermite_genz_keister_lookup_weights ( o, w )

    s = sum ( w(1:o) )

    write ( *, '(2x,i4,2x,i8,2x,g14.6)' ) rule, o, s

    deallocate ( w )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) ' Correct sum:     ', sqrtpi

  return
end
subroutine test38 ( )

!*****************************************************************************80
!
!! TEST38 tabulates the Hermite interpolant and its derivative. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ndp
  integer ( kind = 4 ) ns
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xdp(:)
  real ( kind = 8 ), allocatable :: xs(:)
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: ydp(:)
  real ( kind = 8 ), allocatable :: yp(:)
  real ( kind = 8 ), allocatable :: ys(:)
  real ( kind = 8 ), allocatable :: ysp(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST38'
  write ( *, '(a)' ) '  HERMITE_INTERPOLANT sets up the Hermite interpolant.'
  write ( *, '(a)' ) '  HERMITE_INTERPOLANT_VALUE evaluates it.'
  write ( *, '(a)' ) '  Consider data for y=sin(x) at x=0,1,2,3,4.'

  n = 5
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  allocate ( yp(1:n) )

  nd = 2 * n
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )

  ndp = 2 * n - 1
  allocate ( xdp(1:ndp) )
  allocate ( ydp(1:ndp) )

  call r8vec_linspace ( n, 0.0D+00, 4.0D+00, x )
  y(1:n) = sin ( x(1:n) )
  yp(1:n) = cos ( x(1:n) ) 

  call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
!
!  Now sample the interpolant at NS points, which include data values.
!
  ns = 4 * ( n - 1 ) + 1
  allocate ( xs(1:ns) )
  allocate ( ys(1:ns) )
  allocate ( ysp(1:ns) )

  call r8vec_linspace ( ns, 0.0D+00, 4.0D+00, xs )

  call hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, ys, ysp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, there should be perfect'
  write ( *, '(a)' ) '  agreement between F and H, and F'' and H'''
  write ( *, '(a)' ) '  at the data points X = 0, 1, 2, 3, and 4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In between, H and H'' approximate F and F''.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X(I)          F(X(I))         H(X(I))' // &
    '        F''(X(I))        H''(X(I))'
  write ( *, '(a)' ) ' '
  do i = 1, ns
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, xs(i), sin ( xs(i) ), ys(i), &
      cos ( xs(i) ), ysp(i)
  end do

  deallocate ( x )
  deallocate ( xd )
  deallocate ( xdp )
  deallocate ( xs )
  deallocate ( y )
  deallocate ( yd )
  deallocate ( ydp )
  deallocate ( yp )
  deallocate ( ys )
  deallocate ( ysp )

  return
end
function cubic_antiderivative ( x )

!*****************************************************************************80
!
!! CUBIC_ANTIDERIVATIVE evaluates the antiderivative function of a cubic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) CUBIC_ANTIDERIVATIVE, the value.
!
  implicit none

  real ( kind = 8 ) cubic_antiderivative
  real ( kind = 8 ) x

  cubic_antiderivative = x * x * ( 5.0D+00 + x * ( - 7.0D+00 / 3.0D+00 &
    + x * 1.0D+00 / 4.0D+00 ) )

  return
end
subroutine cubic_integrate ( a, b, q )

!*****************************************************************************80
!
!! CUBIC_INTEGRATE integrates the cubic from A to B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the integration interval.
!
!    Output, real ( kind = 8 ) Q, the integral from A to B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cubic_antiderivative
  real ( kind = 8 ) q

  q = cubic_antiderivative ( b ) - cubic_antiderivative ( a )

  return
end
subroutine cubic_value ( x, f, d, s, t )

!*****************************************************************************80
!
!! CUBIC_VALUE evaluates a cubic function.
!
!  Discussion:
!
!    f(x) =   x^3 -  7 x^2 + 10 x
!    d(x) = 3 x^2 - 14 x   + 10
!    s(x) = 6 x   - 14
!    t(x) = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) F, D, S, T, the value and first three
!    derivatives of the cubic function.
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  f = 0.0D+00 + x * ( 10.0D+00 + x * (  - 7.0D+00 + x * 1.0D+00 ) )
  d =                 10.0D+00 + x * ( - 14.0D+00 + x * 3.0D+00 )
  s =                                  - 14.0D+00 + x * 6.0D+00
  t =                                                   6.0D+00

  return
end
