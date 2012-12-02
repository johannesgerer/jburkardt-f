program main

!*****************************************************************************80
!
!! MAIN is the main program for QUADRULE_PRB.
!
!  Discussion:
!
!    QUADRULE_PRB calls a set of tests for the QUADRULE library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) alpha
  integer ( kind = 4 ) n

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUADRULE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the QUADRULE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test0725 ( )
  call test075 ( )
  call test076 ( )
  call test078 ( )

  n = 5
  alpha = +0.250D+00
  call test079 ( n, alpha )

  n = 5
  alpha = -0.5D+00
  call test079 ( n, alpha )

  call test08 ( )
  call test085 ( )
  
  n = 31
  call test087 ( n )
  
  n = 63
  call test087 ( n )
   
  call test09 ( )
  call test095 ( )
  call test096 ( )

  call test10 ( )
  call test105 ( )
  call test108 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )

  n = 11
  alpha = 0.0D+00
  call test165 ( n, alpha )

  n = 11
  alpha = 0.5D+00
  call test165 ( n, alpha )

  n = 11
  alpha = 2.0D+00
  call test165 ( n, alpha )

  call test17 ( )
!
!  Compare computed and lookup versions of Gauss-Legendre rules.
!
  n = 31
  call test18 ( n )
  n = 64
  call test18 ( n )
  n = 129
  call test18 ( n )
  n = 255
  call test18 ( n )

  call test185 ( )
  call test19 ( )

  call test20 ( )

  n = 127
  call test205 ( n )

  n = 255
  call test205 ( n )

  call test21 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
  call test25 ( )
  call test26 ( )
  call test27 ( )
  call test275 ( )
  call test28 ( )
  call test29 ( )

  call test30 ( )
  call test31 ( )
  call test32 ( )
  call test33 ( )
  call test34 ( )
  call test345 ( )
  call test35 ( )
  call test36 ( )
  call test37 ( )
  call test38 ( )
  call test39 ( )

  call test40 ( )
  call test401 ( )
  call test402 ( )
  call test403 ( )
  call test404 ( )
  call test41 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUADRULE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )
 
  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BASHFORTH_SET and SUMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: n_max = 10
  integer ( kind = 4 ) order
  real ( kind = 8 ), allocatable, dimension ( : ) :: result
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  BASHFORTH_SET sets up an Adams-Bashforth rule;'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is [0,1].'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do n = 1, n_max

      allocate ( w(1:n) ) 
      allocate ( x(1:n) )

      call bashforth_set ( n, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call summer ( function_value, n, x, w, result(i) )
 
      end do

      write ( *, '(2x,i2,2x,5f14.8)' ) order, result(ilo:ihi)

      deallocate ( w )
      deallocate ( x )

    end do

  end do

  deallocate ( result )
 
  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BDFC_SET and BDF_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: n_max = 10

  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) n
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(n_max)
  real      ( kind = 8 ) x(n_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) &
    '  BDFC_SET sets up a Backward Difference Corrector rule;'
  write ( *, '(a)' ) '  BDF_SUM carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is [0,1].'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do n = 1, n_max

      call bdfc_set ( n, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call bdf_sum ( function_value, n, x, w, result(i) )

      end do

      write ( *, '(2x,i2,2x,5f14.8)' ) n, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests BDFP_SET and BDF_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: n_max = 10

  real      ( kind = 8 ) diftab(n_max)
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) n
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ), allocatable, dimension ( : ) :: w
  real      ( kind = 8 ), allocatable, dimension ( : ) :: x

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) &
    '  BDFP_SET sets up a Backward Difference Predictor rule;'
  write ( *, '(a)' ) '  BDF_SUM carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is [0,1].'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do n = 1, n_max

      allocate ( w(1:n) )
      allocate ( x(1:n) )

      call bdfp_set ( n, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call bdf_sum ( function_value, n, x, w, result(i) )

      end do

      write ( *, '(2x,i2,2x,5f14.8)' ) n, result(ilo:ihi)

      deallocate ( w )
      deallocate ( x )

    end do

  end do

  deallocate ( result )
 
  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests CHEB_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: n_max = 9

  real      ( kind = 8 ) a
  real      ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) nsub
  integer   ( kind = 4 ) n
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(n_max)
  real      ( kind = 8 ) x(n_max)
  real      ( kind = 8 ) xhi
  real      ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a =  0.0D+00
  b =  1.0D+00

  nsub = 1

  xlo = -1.0D+00
  xhi =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  CHEB_SET sets up a Chebyshev rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do n = 1, n_max

      if ( n /= 8 ) then

        call cheb_set ( n, x, w )

        do i = ilo, ihi

          call function_set ( 'SET', i )

          call sum_sub ( function_value, a, b, nsub, n, xlo, xhi, &
            x, w, result(i) )
   
        end do

        write ( *, '(2x,i2,2x,5f14.8)' ) n, result(ilo:ihi)

      end if

    end do
 
  end do

  deallocate ( result )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests CHEBYSHEV1_COMPUTE and SUMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: n_max = 6

  real      ( kind = 8 ) a
  real      ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) nsub
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(n_max)
  real      ( kind = 8 ) x(n_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = -1.0D+00
  b =  1.0D+00
  nsub = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  CHEBYSHEV1_COMPUTE sets a Gauss-Chebyshev type 1 rule,'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) '  The weight function is 1 / sqrt ( 1 - X**2 )'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '
    
    do n = 1, n_max
 
      call chebyshev1_compute ( n, x, w )

      do i = ilo, ihi
 
        call function_set ( 'SET', i )

        call summer ( function_value, n, x, w, result(i) )

      end do

      write ( *, '(2x,i2,2x,5f14.8)' ) n, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests CHEBYSHEV2_COMPUTE and SUMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: n_max = 4

  real      ( kind = 8 ) a
  real      ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) nsub
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(n_max)
  real      ( kind = 8 ) x(n_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = - 1.0D+00
  b =   1.0D+00
  nsub = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  CHEBYSHEV2_COMPUTE sets a Gauss-Chebyshev type 2 rule;'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) '  The weight function is sqrt ( 1 - X**2 )'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '
 
    do n = 1, n_max
 
      call chebyshev2_compute ( n, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call summer ( function_value, n, x, w, result(i) )

      end do

      write ( *, '(2x,i2,2x,5f14.8)' ) n, result(ilo:ihi)

    end do
 
  end do

  deallocate ( result )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests CHEBYSHEV3_COMPUTE and SUMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: order_max = 6

  real      ( kind = 8 ) a
  real      ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) nsub
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = -1.0D+00
  b =  1.0D+00
  nsub = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  CHEBYSHEV3_COMPUTE sets a Gauss-Chebyshev type 3 rule,'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) '  The weight function is 1 / sqrt ( 1 - X**2 )'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max
 
      call chebyshev3_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call summer ( function_value, order, x, w, result(i) )

      end do

      write ( *, '(2x,i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test0725 ( )

!*****************************************************************************80
!
!! TEST0725 tests CLENSHAW_CURTIS_COMPUTE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  real    ( kind = 8 ) w(n_max)
  real    ( kind = 8 ) x(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0725'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_COMPUTE computes'
  write ( *, '(a)' ) '  a Clenshaw-Curtis quadrature rule over [-1,1]'
  write ( *, '(a)' ) '  of given order.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order       W                       X'
  write ( *, '(a)' ) ' '

  do n = 1, 10

    call clenshaw_curtis_compute ( n, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8)' ) n

    do i = 1, n
      write ( *, '(10x,2x,g25.16,2x,f24.16)' ) w(i), x(i)
    end do

  end do

  return
end
subroutine test075 ( )

!*****************************************************************************80
!
!! TEST075 tests CLENSHAW_CURTIS_SET and SUMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: order_max = 16

  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST075'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_SET sets up a Clenshaw-Curtis rule;'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is [-1,1].'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '
  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max
 
      call clenshaw_curtis_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call summer ( function_value, order, x, w, result(i) )
 
      end do

      write ( *, '(2x,i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )
 
  return
end
subroutine test076 ( )

!*****************************************************************************80
!
!! TEST076 compares FEJER1_COMPUTE and FEJER1_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real    ( kind = 8 ) w1(order_max)
  real    ( kind = 8 ) w2(order_max)
  real    ( kind = 8 ) x1(order_max)
  real    ( kind = 8 ) x2(order_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST076'
  write ( *, '(a)' ) '  FEJER1_COMPUTE computes a Fejer type 1 rule.'
  write ( *, '(a)' ) '  FEJER1_SET looks up a Fejer type 1 rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare:'
  write ( *, '(a)' ) '    (W1,X1) from FEJER1_SET,'
  write ( *, '(a)' ) '    (W2,X2) from FEJER1_COMPUTE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    Order        W1              W2              X1             X2'
  write ( *, '(a)' ) ' '

  do order = 1, order_max

    call fejer1_set ( order, x1, w1 )
    call fejer1_compute ( order, x2, w2 )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8)' ) order

    do i = 1, order
      write ( *, '(10x,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        w1(i), w2(i), x1(i), x2(i)
    end do

  end do

  return
end
subroutine test078 ( )

!*****************************************************************************80
!
!! TEST078 compares FEJER2_COMPUTE and FEJER2_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real    ( kind = 8 ) w1(order_max)
  real    ( kind = 8 ) w2(order_max)
  real    ( kind = 8 ) x1(order_max)
  real    ( kind = 8 ) x2(order_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST078'
  write ( *, '(a)' ) '  FEJER2_COMPUTE computes a Fejer type 2 rule.'
  write ( *, '(a)' ) '  FEJER2_SET looks up a Fejer type 2 rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare:'
  write ( *, '(a)' ) '    (W1,X1) from FEJER2_SET,'
  write ( *, '(a)' ) '    (W2,X2) from FEJER2_COMPUTE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    Order        W1              W2              X1             X2'
  write ( *, '(a)' ) ' '

  do order = 1, order_max

    call fejer2_set ( order, x1, w1 )
    call fejer2_compute ( order, x2, w2 )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8)' ) order

    do i = 1, order
      write ( *, '(10x,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        w1(i), w2(i), x1(i), x2(i)
    end do

  end do

  return
end
subroutine test079 ( n, alpha )

!*****************************************************************************80
!
!! TEST079 tests GEGENBAUER_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!
  implicit none

  real    ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real    ( kind = 8 ), allocatable, dimension ( : ) :: w
  real    ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST079'
  write ( *, '(a)' ) '  GEGENBAUER_COMPUTE computes a Gauss-Gegenbauer rule;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The printed output of this test can be inserted into'
  write ( *, '(a)' ) '  a FORTRAN program.'

  allocate ( w(n) )
  allocate ( x(n) )

  call gegenbauer_compute ( n, alpha, x, w )
 
  write ( *, '(a)'       ) '!'
  write ( *, '(a)'       ) '!  Abscissas X and weights W for a Gauss Gegenbauer rule'
  write ( *, '(a,i8)'    ) '!  of order = ', n
  write ( *, '(a,g14.6)' ) '!  with ALPHA = ', alpha
  write ( *, '(a)'       ) '!'

  do i = 1, n
    write ( *, '(a,i3,a,g24.16)' ) '    x(', i, ') = ', x(i)
  end do
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(a,i3,a,g24.16)' ) '    w(', i, ') = ', w(i)
  end do

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests HERMITE_EK_COMPUTE and SUMMER.
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

  integer   ( kind = 4 ), parameter :: order_max = 20

  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  HERMITE_EK_COMPUTE computes a Gauss-Hermite rule;'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is ( -oo, +oo ).'
  write ( *, '(a)' ) '  The weight function is exp ( - x * x )'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '
 
    do order = 1, order_max

      call hermite_ek_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call summer ( function_value, order, x, w, result(i) )
 
      end do
 
      write ( *, '(2x,i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do
 
  deallocate ( result )

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests HERMITE_EK_COMPUTE against HERMITE_INTEGRAL.
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

  integer ( kind = 4 ), parameter :: n_max = 10

  real    ( kind = 8 ) error
  real    ( kind = 8 ) estimate
  real    ( kind = 8 ) exact
  real    ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real    ( kind = 8 ), allocatable, dimension ( : ) :: w
  real    ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  HERMITE_EK_COMPUTE computes a Gauss-Hermite rule'
  write ( *, '(a)' ) '  which is appropriate for integrands of the form'
  write ( *, '(a)' ) '    f(x) * exp(-x*x) from -oo to +oo.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HERMITE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  this integal when f(x) = x^m.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         M         N       Estimate       Exact            Error'

  do m = 0, 10, 2
 
    call hermite_integral ( m, exact )

    write ( *, '(a)' ) ' '

    do n = 1, n_max

      allocate ( f(n) )
      allocate ( w(n) )
      allocate ( x(n) )

      call hermite_ek_compute ( n, x, w )
 
      if ( m == 0 ) then
        f(1:n) = 1.0D+00
      else
        f(1:n) = x(1:n)**m
      end if

      estimate = dot_product ( w(1:n), f(1:n) )

      error = abs ( exact - estimate )
  
      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        m, n, estimate, exact, error

      deallocate ( f )
      deallocate ( w )
      deallocate ( x)

    end do

  end do
 
  return
end
subroutine test087 ( order )

!*****************************************************************************80
!
!! TEST087 tests HERMITE_EK_COMPUTE.
!
!  Discussion:
!
!    I used this test to generate tabular values of weights and
!    abscissas for Gauss-Hermite quadrature.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) order

  integer ( kind = 4 ) i
  real    ( kind = 8 ) w(order)
  real    ( kind = 8 ) x(order)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST087'
  write ( *, '(a)' ) '  HERMITE_EK_COMPUTE computes a Gauss-Hermite rule;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Compute the data for ORDER = ', order

  call hermite_ek_compute ( order, x, w )
 
  write ( *, '(a)' ) ' '
  do i = 1, order
    write ( *, '(a,i3,a,g40.32)' ) '    x(', i, ') = ', x(i)
  end do
  write ( *, '(a)' ) ' '
  do i = 1, order
    write ( *, '(a,i3,a,g40.32)' ) '    w(', i, ') = ', w(i)
  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests HERMITE_SET and SUMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: order_max = 20

  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  HERMITE_SET sets up a Gauss-Hermite rule;'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is ( -oo, +oo ).'
  write ( *, '(a)' ) '  The weight function is exp ( - x * x )'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '
 
    do order = 1, order_max
 
      call hermite_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call summer ( function_value, order, x, w, result(i) )
 
      end do

      write ( *, '(2x,i2,2x,5f14.8)' )  order, result(ilo:ihi)
 
    end do
 
  end do

  deallocate ( result )

  return
end
subroutine test095 ( )

!*****************************************************************************80
!
!! TEST095 tests HERMITE_GK16_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l_max = 8

  real    ( kind = 8 ) error
  real    ( kind = 8 ) estimate
  real    ( kind = 8 ) exact
  real    ( kind = 8 ), allocatable :: f(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: n_list(0:l_max) = (/ &
    1, 3, 7, 9, 17, 19, 31, 33, 35 /)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) :: p_list(0:l_max) = (/ &
    1, 5, 7, 15, 17, 29, 31, 33, 51 /)
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST095'
  write ( *, '(a)' ) '  HERMITE_GK16_SET sets up a nested rule'
  write ( *, '(a)' ) '  for the Hermite integration problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is ( -oo, +oo ).'
  write ( *, '(a)' ) '  The weight function is exp ( - x * x )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HERMITE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  the integal when f(x) = x^m.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         M         N       Estimate       Exact            Error'

  do l = 0, l_max

    write ( *, '(a)' ) ' '
    n = n_list(l)
    allocate ( f(1:n) )
    allocate ( x(1:n) )
    allocate ( w(1:n) )

    p = p_list(l)

    call hermite_gk16_set ( n, x, w )

    do m = 0, min ( p + 2, 20 ), 2

      call hermite_integral ( m, exact )

      if ( m == 0 ) then
        f(1:n) = 1.0D+00
      else
        f(1:n) = x(1:n)**m
      end if

      estimate = dot_product ( w(1:n), f(1:n) )

      error = abs ( exact - estimate )
  
      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        m, n, estimate, exact, error

    end do

    deallocate ( f )
    deallocate ( w )
    deallocate ( x )
 
  end do

  return
end
subroutine test096 ( )

!*****************************************************************************80
!
!! TEST096 tests HERMITE_GK**_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) error
  real    ( kind = 8 ) estimate
  real    ( kind = 8 ) exact
  real    ( kind = 8 ), allocatable :: f(:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST096'
  write ( *, '(a)' ) '  HERMITE_GK**_SET sets up a nested rule'
  write ( *, '(a)' ) '  for the Hermite integration problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is ( -oo, +oo ).'
  write ( *, '(a)' ) '  The weight function is exp ( - x * x )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HERMITE_INTEGRAL determines the exact value of'
  write ( *, '(a)' ) '  the integal when f(x) = x^m.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we just test the highest order rule.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HERMITE_GK16_SET:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         M         N       Estimate       Exact            Error'
  write ( *, '(a)' ) ' '

  n = 35
  p = 51

  allocate ( f(1:n) )
  allocate ( w(1:n) )
  allocate ( x(1:n) )

  call hermite_gk16_set ( n, x, w )

  do m = 0, min ( p + 2, 20 ), 2

    call hermite_integral ( m, exact )

    if ( m == 0 ) then
      f(1:n) = 1.0D+00
    else
      f(1:n) = x(1:n)**m
    end if

    estimate = dot_product ( w(1:n), f(1:n) )

    error = abs ( exact - estimate )
  
    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      m, n, estimate, exact, error

  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HERMITE_GK18_SET:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         M         N       Estimate       Exact            Error'
  write ( *, '(a)' ) ' '
  
  n = 37
  p = 55

  allocate ( f(1:n) )
  allocate ( w(1:n) )
  allocate ( x(1:n) )

  call hermite_gk18_set ( n, x, w )

  do m = 0, min ( p + 2, 20 ), 2

    call hermite_integral ( m, exact )

    if ( m == 0 ) then
      f(1:n) = 1.0D+00
    else
      f(1:n) = x(1:n)**m
    end if

    estimate = dot_product ( w(1:n), f(1:n) )

    error = abs ( exact - estimate )
  
    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      m, n, estimate, exact, error

  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HERMITE_GK22_SET:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         M         N       Estimate       Exact            Error'
  write ( *, '(a)' ) ' '
  n = 41
  p = 63

  allocate ( f(1:n) )
  allocate ( w(1:n) )
  allocate ( x(1:n) )

  call hermite_gk22_set ( n, x, w )

  do m = 0, min ( p + 2, 20 ), 2

    call hermite_integral ( m, exact )

    if ( m == 0 ) then
      f(1:n) = 1.0D+00
    else
      f(1:n) = x(1:n)**m
    end if

    estimate = dot_product ( w(1:n), f(1:n) )

    error = abs ( exact - estimate )
  
    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      m, n, estimate, exact, error

  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HERMITE_GK24_SET:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         M         N       Estimate       Exact            Error'
  write ( *, '(a)' ) ' '

  n = 43
  p = 67

  allocate ( f(1:n) )
  allocate ( w(1:n) )
  allocate ( x(1:n) )

  call hermite_gk24_set ( n, x, w )

  do m = 0, min ( p + 2, 20 ), 2

    call hermite_integral ( m, exact )

    if ( m == 0 ) then
      f(1:n) = 1.0D+00
    else
      f(1:n) = x(1:n)**m
    end if

    estimate = dot_product ( w(1:n), f(1:n) )

    error = abs ( exact - estimate )
  
    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      m, n, estimate, exact, error

  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests JACOBI_EK_COMPUTE and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 10

  real    ( kind = 8 ) a
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  integer ( kind = 4 ) test
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = -1.0D+00
  b =  1.0D+00

  nsub = 1
  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  JACOBI_EK_COMPUTE computes a Gauss-Jacobi rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'

  do test = 1, 2

    if ( test == 1 ) then
      alpha = 0.0D+00
      beta = 0.0D+00
    else if ( test == 2 ) then
      alpha = 1.0D+00
      beta = 0.0D+00
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
    write ( *, '(a,g14.6)' ) '  BETA =  ', beta

    do ilo = 1, function_num, 5

      ihi = min ( ilo + 4, function_num )

      write ( *, '(a)' ) ' '
      write ( *, '(2x,a7, 5(a10,4x) )' ) &
        'Order  ', ( function_name ( i ), i = ilo, ihi )
      write ( *, '(a)' ) ' '

      do order = 1, order_max

        call jacobi_ek_compute ( order, alpha, beta, x, w )

        do i = ilo, ihi

          call function_set ( 'SET', i )

          call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
            x, w, result(i) )

        end do

        write ( *, '(2x,i2,2x,5f14.8)' ) order, result(ilo:ihi)

      end do

    end do

  end do 

  deallocate ( result )

  return
end
subroutine test105 ( )

!*****************************************************************************80
!
!! TEST105 tests JACOBI_EK_COMPUTE and JACOBI_SS_COMPUTE.
!
!  Discussion:
!
!    Compare with tabular values on page 178 of Stroud and Secrest.
!
!     In particular,
!
!             X              W
!
!     1  -0.9833999115   0.4615276287E-03
!     2  -0.9447138932   0.2732249104E-02
!     3  -0.8849310847   0.8045830455E-02
!    ..  .............   ................
!    19   0.9656375637   0.7613987785E-01
!    20   0.9934477866   0.3348337670E-01
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 20

  real    ( kind = 8 ) a
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  real    ( kind = 8 ) w(order)
  real    ( kind = 8 ) x(order)

  a = -1.0D+00
  b =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105'
  write ( *, '(a)' ) '  JACOBI_EK_COMPUTE computes a Gauss-Jacobi rule;'
  write ( *, '(a)' ) '  JACOBI_SS_COMPUTE computes a Gauss-Jacobi rule;'
  write ( *, '(a)' ) '  Here, we simply compute a single rule and'
  write ( *, '(a)' ) '  print it, to check for accuracy.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'

  alpha = 1.5D+00
  beta = 0.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', order
  write ( *, '(a,f14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,f14.6)' ) '  BETA =  ', beta

  call jacobi_ek_compute ( order, alpha, beta, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JACOBI_EK_COMPUTE:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I        X(I)            W(I)'
  write ( *, '(a)' ) ' '

  do i = 1, order
    write ( *, '(2x,i4,2x,f14.8,2x,f14.8)' ) i, x(i), w(i)
  end do

  call jacobi_ss_compute ( order, alpha, beta, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JACOBI_SS_COMPUTE:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I        X(I)            W(I)'
  write ( *, '(a)' ) ' '

  do i = 1, order
    write ( *, '(2x,i4,2x,f14.8,2x,f14.8)' ) i, x(i), w(i)
  end do

  return
end
subroutine test108 ( )

!*****************************************************************************80
!
!! TEST108 tests JACOBI_EK_COMPUTE.
!
!  Discussion:
!
!    I used this test to generate tabular values of weights and
!    abscissas for Gauss-Jacobi quadrature.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 10 

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  real    ( kind = 8 ) w(order)
  real    ( kind = 8 ) x(order)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST108'
  write ( *, '(a)' ) '  JACOBI_EK_COMPUTE computes a Gauss-Jacobi rule;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The printed output of this test can be inserted into'
  write ( *, '(a)' ) '  a FORTRAN program.'

  alpha = 0.5D+00
  beta  = 2.0D+00

  call jacobi_ek_compute ( order, alpha, beta, x, w )
 
  write ( *, '(a)'       ) '!'
  write ( *, '(a)'       ) '!  Abscissas X and weights W for a Gauss Jacobi rule'
  write ( *, '(a,i8)'    ) '!  of order = ', order
  write ( *, '(a,g14.6)' ) '!  with ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '!  and BETA = ', beta
  write ( *, '(a)'       ) '!'

  do i = 1, order
    write ( *, '(a,i2,a,g24.16)' ) '    x(', i, ') = ', x(i)
  end do

  write ( *, '(a)' ) ' '

  do i = 1, order
    write ( *, '(a,i2,a,g24.16)' ) '    w(', i, ') = ', w(i)
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests KRONROD_SET, LEGENDRE_SET and SUMMER_GK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: orderg = 10
  integer ( kind = 4 ), parameter :: orderk = 2 * orderg + 1

  real    ( kind = 8 ), external :: fx2sd1
  real    ( kind = 8 ) resultg
  real    ( kind = 8 ) resultk
  real    ( kind = 8 ) wg(orderg)
  real    ( kind = 8 ) wk(orderk)
  real    ( kind = 8 ) xg(orderg)
  real    ( kind = 8 ) xk(orderk)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  KRONROD_SET sets up a Kronrod rule;'
  write ( *, '(a)' ) '  LEGENDRE_SET sets up a Gauss-Legendre rule;'
  write ( *, '(a)' ) '  SUMMER_GK carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is [-1, 1].'
  write ( *, '(a)' ) '  Integrand is x * x / SQRT ( 1.1 - x * x ).'
  write ( *, '(a)' ) ' '

  call legendre_set ( orderg, xg, wg )

  call kronrod_set ( orderk, xk, wk )

  call summer_gk ( fx2sd1, orderg, wg, resultg, orderk, xk, &
    wk, resultk )

  write ( *, '(2x,i2,2x,g16.10)' ) orderg, resultg
  write ( *, '(2x,i2,2x,g16.10)' ) orderk, resultk
  write ( *, '(2x,2x,2x,g16.10)' )          resultg-resultk

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests KRONROD_SET, LEGENDRE_SET and SUM_SUB_GK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: orderg = 7
  integer ( kind = 4 ), parameter :: orderk = 2 * orderg + 1

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) error
  real    ( kind = 8 ), external :: fx2sd1
  integer ( kind = 4 ) nsub
  real    ( kind = 8 ) resultg
  real    ( kind = 8 ) resultk
  real    ( kind = 8 ) wg(orderg)
  real    ( kind = 8 ) wk(orderk)
  real    ( kind = 8 ) xg(orderg)
  real    ( kind = 8 ) xk(orderk)

  a = - 1.0D+00
  b =   1.0D+00
  nsub = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  KRONROD_SET sets up a Kronrod rule;'
  write ( *, '(a)' ) '  LEGENDRE_SET sets up a Gauss-Legendre rule;'
  write ( *, '(a)' ) '  SUM_SUB_GK carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Integrand is X**2 / SQRT ( 1.1 - X**2 ).'
  write ( *, '(a)' ) ' '

  call legendre_set ( orderg, xg, wg )

  call kronrod_set ( orderk, xk, wk )

  call sum_sub_gk ( fx2sd1, a, b, nsub, orderg, wg, resultg, &
    orderk, xk, wk, resultk, error )

  write ( *, '(2x,i2,2x,g16.10)' ) orderg, resultg
  write ( *, '(2x,i2,2x,g16.10)' ) orderk, resultk
  write ( *, '(2x,2x,2x,g16.10)' )          error

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests LAGUERRE_EK_COMPUTE and LAGUERRE_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20

  real      ( kind = 8 ) a
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  LAGUERRE_EK_COMPUTE computes a Gauss-Laguerre rule;'
  write ( *, '(a)' ) '  LAGUERRE_SUM carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) '  The weight function is exp ( - x )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a)' ) '  The integration interval is [ ', &
    a, ', +oo ).'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call laguerre_ek_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call laguerre_sum ( function_value, a, order, x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)
 
    end do

  end do

  deallocate ( result )

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests LAGUERRE_EK_COMPUTE and LAGUERRE_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20

  real      ( kind = 8 ) a
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  LAGUERRE_EK_COMPUTE computes a Gauss-Laguerre rule;'
  write ( *, '(a)' ) '  LAGUERRE_SUM carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) '  The weight function is exp ( - x )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a)' ) '  The integration interval is [ ',  &
    a, ', +oo ).'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call laguerre_ek_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call laguerre_sum ( function_value, a, order, x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)
 
    end do

  end do

  deallocate ( result )

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests LAGUERRE_EK_COMPUTE and LAGUERRE_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20

  real      ( kind = 8 ) a
  real      ( kind = 8 ) alpha
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  LAGUERRE_EK_COMPUTE computes a Gauss-Laguerre rule;'
  write ( *, '(a)' ) '  LAGUERRE_SUM carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a,g14.6,a)' ) '  The integration interval is [ ', &
    a, ', +oo ).'
  write ( *, '(a)' ) '  The weight function is exp ( - x )'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5
  
    ihi = min ( ilo + 4, function_num )
    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call laguerre_ek_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call laguerre_sum ( function_value, a, order, x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5g14.6)' ) order, result(ilo:ihi)
 
    end do

  end do

  deallocate ( result )

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests GEN_LAGUERRE_EK_COMPUTE and LAGUERRE_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20

  real      ( kind = 8 ) a
  real      ( kind = 8 ) alpha
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  alpha = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  GEN_LAGUERRE_EK_COMPUTE computes a Gauss-Laguerre rule;'
  write ( *, '(a)' ) '  LAGUERRE_SUM carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a,g14.6)' ) &
    '  The weight function is exp ( - x ) * x ^ ', alpha
  write ( *, '(a,g14.6,a)' ) '  The integration interval is [ ', &
    a, ', +oo ).'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call gen_laguerre_ek_compute ( order, alpha, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call laguerre_sum ( function_value, a, order, x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5g14.6)' ) order, result(ilo:ihi)
 
    end do

  end do

  deallocate ( result )

  return
end
subroutine test165 ( order, alpha )

!*****************************************************************************80
!
!! TEST165 tests GEN_LAGUERRE_EK_COMPUTE.
!
!  Discussion:
!
!    I used this test to generate tabular values of weights and
!    abscissas for generalized Gauss-Laguerre quadrature.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!
  implicit none

  real    ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: w
  real    ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST165'
  write ( *, '(a)' ) '  GEN_LAGUERRE_EK_COMPUTE computes a generalized Gauss-Laguerre rule;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The printed output of this test can be inserted into'
  write ( *, '(a)' ) '  a FORTRAN program.'

  allocate ( w(order) )
  allocate ( x(order) )

  call gen_laguerre_ek_compute ( order, alpha, x, w )
 
  write ( *, '(a)'       ) '!'
  write ( *, '(a)'       ) '!  Abscissas X and weights W for a Gauss Laguerre rule'
  write ( *, '(a,i8)'    ) '!  of order = ', order
  write ( *, '(a,g14.6)' ) '!  with ALPHA = ', alpha
  write ( *, '(a)'       ) '!'

  do i = 1, order
    write ( *, '(a,i3,a,g24.16)' ) '    x(', i, ') = ', x(i)
  end do
  write ( *, '(a)' ) ' '
  do i = 1, order
    write ( *, '(a,i3,a,g24.16)' ) '    w(', i, ') = ', w(i)
  end do

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests LAGUERRE_SET and LAGUERRE_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: order_max = 20

  real      ( kind = 8 ) a
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  LAGUERRE_SET sets up a Gauss-Laguerre rule;'
  write ( *, '(a)' ) '  LAGUERRE_SUM carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a)' ) '  The integration interval is [ ', &
    a, ', +oo ).'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) '  The weight function is EXP ( - X ).'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call laguerre_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call laguerre_sum ( function_value, a, order, x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)
 
    end do

  end do
 
  deallocate ( result )

  return
end
subroutine test18 ( n )

!*****************************************************************************80
!
!! TEST18 compares LEGENDRE_EK_COMPUTE and LEGENDRE_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwdifmax
  integer ( kind = 4 ) ixdifmax
  real    ( kind = 8 ) p1(1)
  real    ( kind = 8 ) p2(1)
  real    ( kind = 8 ) wdifmax
  real    ( kind = 8 ) w1(n)
  real    ( kind = 8 ) w2(n)
  real    ( kind = 8 ) xdifmax
  real    ( kind = 8 ) x1(n)
  real    ( kind = 8 ) x2(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  LEGENDRE_EK_COMPUTE computes a Gauss-Legendre rule;'
  write ( *, '(a)' ) '  LEGENDRE_SET sets up a Gauss-Legendre rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Compare the data for N = ', n

  call legendre_ek_compute ( n, x1, w1 )
  call legendre_set ( n, x2, w2 )
 
  xdifmax = 0.0D+00
  ixdifmax = -1

  wdifmax = 0.0D+00
  iwdifmax = -1

  do i = 1, n

    if ( xdifmax < abs ( x1(i) - x2(i) ) ) then
      xdifmax = abs ( x1(i) - x2(i) )
      ixdifmax = i
    end if

    if ( wdifmax < abs ( w1(i) - w2(i) ) ) then
      wdifmax = abs ( w1(i) - w2(i) )
      iwdifmax = i
    end if

  end do

  if ( 0 < ixdifmax ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,g25.18)' ) '  Maximum abscissa difference is ', xdifmax
    write ( *, '(a,i8)'     ) '  for index I = ', ixdifmax
    write ( *, '(a,g25.18)' ) '  Computed X1:', x1(ixdifmax)
    write ( *, '(a,g25.18)' ) '  Stored   X2:', x2(ixdifmax)
    call legendre_polynomial_value ( 1, n, x1(ixdifmax), p1 )
    call legendre_polynomial_value ( 1, n, x2(ixdifmax), p2 )
    write ( *, '(a,g25.18)' ) '  Computed P(X1):', p1(1)
    write ( *, '(a,g25.18)' ) '  Stored   P(X2):', p2(1)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The computed and stored abscissas are identical.'
  end if

  if ( 0 < iwdifmax ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,g25.18)' ) '  Maximum weight difference is   ', wdifmax
    write ( *, '(a,i8)'     ) '  for index I = ', iwdifmax
    write ( *, '(a,g25.18)' ) '  Computed:', w1(iwdifmax)
    write ( *, '(a,g25.18)' ) '  Stored:  ', w2(iwdifmax)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The computed and stored weights are identical.'
  end if
 
  return
end
subroutine test185 ( )

!*****************************************************************************80
!
!! TEST185 tests LEGENDRE_EK_COMPUTE.
!
!  Discussion:
!
!    I used this test to generate tabular values of weights and
!    abscissas for Gauss-Legendre quadrature.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 31

  integer ( kind = 4 ) i
  real    ( kind = 8 ) w(n)
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST185'
  write ( *, '(a)' ) '  LEGENDRE_EK_COMPUTE computes a Gauss-Legendre rule;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Compute the data for N = ', n

  call legendre_ek_compute ( n, x, w )
 
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(a,i3,a,g24.16)' ) '    x(', i, ') = ', x(i)
  end do
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(a,i3,a,g24.16)' ) '    w(', i, ') = ', w(i)
  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests LEGENDRE_EK_COMPUTE and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ), external :: fx2sd1
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) nsub
  real    ( kind = 8 ) result
  real    ( kind = 8 ) w(n)
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  a = 0.0D+00
  b = 1.0D+00

  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  LEGENDRE_EK_COMPUTE computes a Gauss-Legendre rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out over subintervals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Here, we use a fixed order N = ', n
  write ( *, '(a)' ) '  and use more and more subintervals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NSUB     Integral'
  write ( *, '(a)' ) ' '
 
  call legendre_ek_compute ( n, x, w )
 
  do iexp = 0, 9

    nsub = 2**iexp

    call sum_sub ( fx2sd1, a, b, nsub, n, xlo, xhi, x, w, result )

    write ( *, '(i4,g16.8)' ) nsub, result

  end do
 
  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests LEGENDRE_EK_COMPUTE and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 10

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1

  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  LEGENDRE_EK_COMPUTE computes a Gauss-Legendre rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call legendre_ek_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )
 
  return
end
subroutine test205 ( n )

!*****************************************************************************80
!
!! TEST205: LEGENDRE_**_COMPUTE and LEGENDRE_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) p1(n)
  real    ( kind = 8 ) p2(n)
  real    ( kind = 8 ) p3(n)
  real    ( kind = 8 ) p4(n)
  real    ( kind = 8 ) p5(n)
  real    ( kind = 8 ) t0
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) t3
  real    ( kind = 8 ) t4
  real    ( kind = 8 ) t5
  real    ( kind = 8 ) w1(n)
  real    ( kind = 8 ) w2(n)
  real    ( kind = 8 ) w3(n)
  real    ( kind = 8 ) w4(n)
  real    ( kind = 8 ) w5(n)
  real    ( kind = 8 ) x1(n)
  real    ( kind = 8 ) x2(n)
  real    ( kind = 8 ) x3(n)
  real    ( kind = 8 ) x4(n)
  real    ( kind = 8 ) x5(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST205'
  write ( *, '(a)' ) '  Compare timings and accuracy of'
  write ( *, '(a)' ) '  LEGENDRE_DR_COMPUTE,'
  write ( *, '(a)' ) '  LEGENDRE_EK_COMPUTE,'
  write ( *, '(a)' ) '  LEGENDRE_GW_COMPUTE,'
  write ( *, '(a)' ) '  LEGENDRE_SS_COMPUTE,'
  write ( *, '(a)' ) '  LEGENDRE_SET.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Compare the data for N = ', n

  call cpu_time ( t0 )
  call legendre_dr_compute ( n, x1, w1 )
  call cpu_time ( t1 )
  call legendre_ek_compute ( n, x2, w2 )
  call cpu_time ( t2 )
  call legendre_gw_compute ( n, x3, w3 )
  call cpu_time ( t3 )
  call legendre_ss_compute ( n, x4, w4 )
  call cpu_time ( t4 )
  call legendre_set ( n, x5, w5 )
  call cpu_time ( t5 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Timings in seconds:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LEGENDRE_DR_COMPUTE:  ', t1 - t0
  write ( *, '(a,g14.6)' ) '  LEGENDRE_EK_COMPUTE:  ', t2 - t1
  write ( *, '(a,g14.6)' ) '  LEGENDRE_GW_COMPUTE:  ', t3 - t2
  write ( *, '(a,g14.6)' ) '  LEGENDRE_SS_COMPUTE:  ', t4 - t3
  write ( *, '(a,g14.6)' ) '  LEGENDRE_SET:         ', t5 - t4

  call legendre_polynomial_value ( n, n, x1, p1 )
  call legendre_polynomial_value ( n, n, x2, p2 )
  call legendre_polynomial_value ( n, n, x3, p3 )
  call legendre_polynomial_value ( n, n, x4, p4 )
  call legendre_polynomial_value ( n, n, x5, p5 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Accuracy: Max ( Abs ( L(N,X) ):'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LEGENDRE_DR_COMPUTE:  ', maxval ( abs ( p1(1:n) ) )
  write ( *, '(a,g14.6)' ) '  LEGENDRE_EK_COMPUTE:  ', maxval ( abs ( p2(1:n) ) )
  write ( *, '(a,g14.6)' ) '  LEGENDRE_GW_COMPUTE:  ', maxval ( abs ( p3(1:n) ) )
  write ( *, '(a,g14.6)' ) '  LEGENDRE_SS_COMPUTE:  ', maxval ( abs ( p4(1:n) ) )
  write ( *, '(a,g14.6)' ) '  LEGENDRE_SET:         ', maxval ( abs ( p5(1:n) ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Accuracy: Abs ( 2 - sum ( W ):'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LEGENDRE_DR_COMPUTE:  ', abs ( 2.0D+00 - sum ( w1(1:n) ) )
  write ( *, '(a,g14.6)' ) '  LEGENDRE_EK_COMPUTE:  ', abs ( 2.0D+00 - sum ( w2(1:n) ) )
  write ( *, '(a,g14.6)' ) '  LEGENDRE_GW_COMPUTE:  ', abs ( 2.0D+00 - sum ( w3(1:n) ) )
  write ( *, '(a,g14.6)' ) '  LEGENDRE_SS_COMPUTE:  ', abs ( 2.0D+00 - sum ( w4(1:n) ) )
  write ( *, '(a,g14.6)' ) '  LEGENDRE_SET:         ', abs ( 2.0D+00 - sum ( w5(1:n) ) )
 
  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests LEGENDRE_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1

  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  LEGENDRE_SET sets up a Gauss-Legendre rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call legendre_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, x, &
          w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do
 
  end do

  deallocate ( result )

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests LEGENDRE_SET, LEGENDRE_X0_01_SET and RULE_ADJUST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 5

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  integer ( kind = 4 ) i
  real    ( kind = 8 ) w1(order)
  real    ( kind = 8 ) w2(order)
  real    ( kind = 8 ) w3(order)
  real    ( kind = 8 ) x1(order)
  real    ( kind = 8 ) x2(order)
  real    ( kind = 8 ) x3(order)

  a = -1.0D+00
  b = +1.0D+00
  c =  0.0D+00
  d =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  LEGENDRE_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating F(X) over [-1,1];'
  write ( *, '(a)' ) '  RULE_ADJUST adjusts a rule for a new interval.'
  write ( *, '(a)' ) '  LEGENDRE_X0_01_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating F(X) over [0,1];'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will use LEGENDRE_SET to get a rule for [-1,1],'
  write ( *, '(a)' ) '  adjust it to [0,1] using RULE_ADJUST,'
  write ( *, '(a)' ) '  and compare the results of LEGENDRE_X0_01_SET.'
  write ( *, '(a)' ) ' '

  call legendre_set ( order, x1, w1 )

  x2(1:order) = x1(1:order)
  w2(1:order) = w1(1:order)

  call rule_adjust ( a, b, c, d, order, x2, w2 )

  call legendre_x0_01_set ( order, x3, w3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Abscissas:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Original        Adjusted        Stored'
  write ( *, '(a)' ) ' '

  do i = 1, order
    write ( *, '(i2,3f16.12)' ) i, x1(i), x2(i), x3(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Original        Adjusted        Stored'
  write ( *, '(a)' ) ' '

  do i = 1, order
    write ( *, '(i2,3f16.12)' ) i, w1(i), w2(i), w3(i)
  end do

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests LEGENDRE_COS_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ) :: pi = 3.141592653589793D+00
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  a = -0.5D+00 * pi
  b = +0.5D+00 * pi

  nsub = 1

  xlo = -0.5D+00 * pi
  xhi = +0.5D+00 * pi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  LEGENDRE_COS_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    over [-PI/2,PI/2] with weight function COS(X);'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do iexp = 0, 4

      order = 2**iexp

      call legendre_cos_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, x, &
          w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests LEGENDRE_SQRTX_01_SET and SUMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: order_max = 20

  real      ( kind = 8 ) a
  real      ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer   ( kind = 4 ) function_num
  real      ( kind = 8 ), external :: function_value
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) iexp
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) order
  real      ( kind = 8 ), allocatable, dimension ( : ) :: result
  real      ( kind = 8 ) w(order_max)
  real      ( kind = 8 ) x(order_max)

  a = 0.0D+00
  b = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  LEGENDRE_SQRTX_01_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    over [0,1] with weight function SQRT(X);'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do iexp = 0, 3

      order = 2**iexp

      call legendre_sqrtx_01_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call summer ( function_value, order, x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests LEGENDRE_SQRTX2_01_SET and SUMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)

  a = 0.0D+00
  b = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  LEGENDRE_SQRTX2_01_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    over [0,1] with weight function 1/SQRT(X);'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do iexp = 0, 3

      order = 2**iexp

      call legendre_sqrtx2_01_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call summer ( function_value, order, x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests LEGENDRE_COS2_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ) :: pi = 3.141592653589793D+00
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  a = -0.5D+00 * pi
  b = +0.5D+00 * pi

  nsub = 1
  xlo = -0.5D+00 * pi
  xhi = +0.5D+00 * pi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  LEGENDRE_COS2_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    over [0,PI/2] with weight function COS(X);'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do iexp = 1, 4

      order = 2**iexp

      call legendre_cos2_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, x, &
          w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests LEGENDRE_LOG_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20
  integer ( kind = 4 ), parameter :: test_num = 9

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) order
  integer ( kind = 4 ), dimension ( test_num ) :: order_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 16 /)
  integer ( kind = 4 ) nsub
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  integer ( kind = 4 ) test
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1

  xlo = 0.0D+00
  xhi = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  LEGENDRE_LOG_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    to integrate -LOG(X)*F(X) over [0,1];'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do test = 1, test_num

      order = order_test(test)
 
      call legendre_log_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order,  xlo, xhi, &
          x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do
 
  end do

  deallocate ( result )

  return
end
subroutine test275 ( )

!*****************************************************************************80
!
!! TEST275 tests LEGENDRE_LOG_COMPUTE and SUM_SUB.
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

  integer ( kind = 4 ), parameter :: order_max = 20
  integer ( kind = 4 ), parameter :: test_num = 13

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) order
  integer ( kind = 4 ), dimension ( test_num ) :: order_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 /)
  integer ( kind = 4 ) nsub
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  integer ( kind = 4 ) test
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1

  xlo = 0.0D+00
  xhi = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST275'
  write ( *, '(a)' ) '  LEGENDRE_LOG_COMPUTE computes a Gauss-Legendre rule'
  write ( *, '(a)' ) '    to integrate -LOG(X)*F(X) over [0,1];'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do test = 1, test_num

      order = order_test(test)
 
      call legendre_log_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order,  xlo, xhi, &
          x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do
 
  end do

  deallocate ( result )
!
!  DEBUG
!
  do test = 1, test_num

    order = order_test(test)
 
    call legendre_log_compute ( order, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Rule for order ', order
    write ( *, '(a)' ) ' '

    do i = 1, order
      write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, x(i), w(i)
    end do
 
  end do

  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests LEGENDRE_X0_01_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 8

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1

  xlo =  0.0D+00
  xhi =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  LEGENDRE_X0_01_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating F(X) over [0,1];'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call legendre_x0_01_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, & 
          x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests LEGENDRE_X1_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 9

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a =  0.0D+00
  b =  1.0D+00

  nsub = 1

  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  LEGENDRE_X1_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating ( 1 + X ) * F(X) over [-1,1];'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call legendre_x1_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test30 ( )

!*****************************************************************************80
!
!! TEST30 tests LEGENDRE_X1_SET, LEGENDRE_X1_01_SET and RULE_ADJUST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 5

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  integer ( kind = 4 ) i
  real    ( kind = 8 ) w1(order)
  real    ( kind = 8 ) w2(order)
  real    ( kind = 8 ) w3(order)
  real    ( kind = 8 ) x1(order)
  real    ( kind = 8 ) x2(order)
  real    ( kind = 8 ) x3(order)

  a = -1.0D+00
  b = +1.0D+00
  c =  0.0D+00
  d =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  LEGENDRE_X1_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating ( 1 + X ) * F(X) over [-1,1];'
  write ( *, '(a)' ) '  RULE_ADJUST adjusts a rule for a new interval.'
  write ( *, '(a)' ) '  LEGENDRE_X1_01_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating X * F(X) over [0,1];'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will use LEGENDRE_X1_SET to get a rule for [-1,1],'
  write ( *, '(a)' ) '  adjust it to [0,1] using RULE_ADJUST,'
  write ( *, '(a)' ) '  make further adjustments because the weight function'
  write ( *, '(a)' ) '  is not 1, '
  write ( *, '(a)' ) '  and compare the results of LEGENDRE_X1_01_SET.'
  write ( *, '(a)' ) ' '

  call legendre_x1_set ( order, x1, w1 )

  x2(1:order) = x1(1:order)
  w2(1:order) = w1(1:order)

  call rule_adjust ( a, b, c, d, order, x2, w2 )
!
!  Because the weight function W(X) is not 1, we need to do more
!  adjustments to the weight vector.
!
  w2(1:order) = w2(1:order) / 2.0D+00

  call legendre_x1_01_set ( order, x3, w3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Abscissas:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original        Adjusted         Stored'
  write ( *, '(a)' ) ' '

  do i = 1, order
    write ( *, '(i2,3f16.12)' ) i, x1(i), x2(i), x3(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)') '  Weights:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original  Adjusted Stored'
  write ( *, '(a)' ) ' '

  do i = 1, order
    write ( *, '(i2,3f16.12)' ) i, w1(i), w2(i), w3(i)
  end do

  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests LEGENDRE_X1_01_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 8

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1

  xlo = 0.0D+00
  xhi = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  LEGENDRE_X1_01_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating X * F(X) over [0,1];'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call legendre_x1_01_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests LEGENDRE_X2_SET and SUM_SUB
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 9

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1
  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  LEGENDRE_X2_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating (1+X)**2 * F(X) over [-1,1];'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call legendre_x2_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test33 ( )

!*****************************************************************************80
!
!! TEST33 tests LEGENDRE_X2_SET, LEGENDRE_X2_01_SET and RULE_ADJUST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 5

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  integer ( kind = 4 ) i
  real    ( kind = 8 ) w1(order)
  real    ( kind = 8 ) w2(order)
  real    ( kind = 8 ) w3(order)
  real    ( kind = 8 ) x1(order)
  real    ( kind = 8 ) x2(order)
  real    ( kind = 8 ) x3(order)

  a = -1.0D+00
  b = +1.0D+00
  c =  0.0D+00
  d =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  LEGENDRE_X2_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating ( 1 + X )^2 * F(X) over [-1,1];'
  write ( *, '(a)' ) '  RULE_ADJUST adjusts a rule for a new interval.'
  write ( *, '(a)' ) '  LEGENDRE_X2_01_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating X^2 * F(X) over [0,1];'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will use LEGENDRE_X2_SET to get a rule for [-1,1],'
  write ( *, '(a)' ) '  adjust it to [0,1] using RULE_ADJUST,'
  write ( *, '(a)' ) '  make further adjustments because the weight function'
  write ( *, '(a)' ) '  is not 1, '
  write ( *, '(a)' ) '  and compare the results of LEGENDRE_X2_01_SET.'
  write ( *, '(a)' ) ' '

  call legendre_x2_set ( order, x1, w1 )

  x2(1:order) = x1(1:order)
  w2(1:order) = w1(1:order)

  call rule_adjust ( a, b, c, d, order, x2, w2 )
!
!  Because the weight function W(X) is not 1, we need to do more
!  adjustments to the weight vector.
!
  w2(1:order) = w2(1:order) / 4.0D+00

  call legendre_x2_01_set ( order, x3, w3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Abscissas:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original  Adjusted Stored'
  write ( *, '(a)' ) ' '

  do i = 1, order
    write ( *, '(i2,3f16.12)' ) i, x1(i), x2(i), x3(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original  Adjusted Stored'
  write ( *, '(a)' ) ' '

  do i = 1, order
    write ( *, '(i2,3f16.12)' ) i, w1(i), w2(i), w3(i)
  end do

  return
end
subroutine test34 ( )

!*****************************************************************************80
!
!! TEST34 tests LEGENDRE_X2_01_SET and SUM_SUB
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 8

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1
  xlo = 0.0D+00
  xhi = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  LEGENDRE_X2_01_SET sets up a Gauss-Legendre rule'
  write ( *, '(a)' ) '    for integrating X*X * F(X) over [0,1];'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '
  
    do order = 1, order_max

      call legendre_x2_01_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test345 ( )

!*****************************************************************************80
!
!! TEST345 tests LOBATTO_COMPUTE and LOBATTO_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real    ( kind = 8 ), allocatable, dimension ( : ) :: w1
  real    ( kind = 8 ), allocatable, dimension ( : ) :: w2
  real    ( kind = 8 ), allocatable, dimension ( : ) :: x1
  real    ( kind = 8 ), allocatable, dimension ( : ) :: x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST345'
  write ( *, '(a)' ) '  LOBATTO_COMPUTE computes a Lobatto rule;'
  write ( *, '(a)' ) '  LOBATTO_SET sets a rule from a table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         I      X1            X2            W1            W2'
 
  do n = 4, 12, 3

    allocate ( w1(1:n) )
    allocate ( w2(1:n) )
    allocate ( x1(1:n) )
    allocate ( x2(1:n) )

    call lobatto_compute ( n, x1, w1 )
    call lobatto_set ( n, x2, w2 )

    write ( * , '(a)' ) ' '
    do i = 1, n
      write (  *, '(2x,i8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8)' ) &
        i, x1(i), x2(i), w1(i), w2(i)
    end do

    deallocate ( w1 )
    deallocate ( w2 )
    deallocate ( x1 )
    deallocate ( x2 )

  end do
 
  return
end
subroutine test35 ( )

!*****************************************************************************80
!
!! TEST35 tests LOBATTO_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 20

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = -1.0D+00
  b =  1.0D+00

  nsub = 1
  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  LOBATTO_SET sets up a Lobatto rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '
  
    do order = 1, order_max
 
      if ( order == 1 ) then
        cycle
      end if

      call lobatto_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do
 
  deallocate ( result )

  return
end
subroutine test36 ( )

!*****************************************************************************80
!
!! TEST36 tests MOULTON_SET and SUMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 10

  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36'
  write ( *, '(a)' ) '  MOULTON_SET sets up an Adams-Moulton rule;'
  write ( *, '(a)' ) '  SUMMER carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is [0,1].'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max
 
      call moulton_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call summer ( function_value, order, x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do

  deallocate ( result )

  return
end
subroutine test37 ( )

!*****************************************************************************80
!
!! TEST37 tests NCC_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 21

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1

  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  NCC_SET sets up a closed Newton-Cotes rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call ncc_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do
 
  end do

  deallocate ( result )

  return
end
subroutine test38 ( )

!*****************************************************************************80
!
!! TEST38 tests NCC_COMPUTE and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 21

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1
  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST38'
  write ( *, '(a)' ) '  NCC_COMPUTE computes a closed Newton-Cotes rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '
  
    do order = 1, order_max
 
      call ncc_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do
 
  deallocate ( result )

  return
end
subroutine test39 ( )

!*****************************************************************************80
!
!! TEST39 tests NCO_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 9

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1
  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST39'
  write ( *, '(a)' ) '  NCO_SET sets up an open Newton-Cotes rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)') '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max
 
      if ( order <= 7 .or. order == 9 ) then

        call nco_set ( order, x, w )

        do i = ilo, ihi

          call function_set ( 'SET', i )
 
          call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
            x, w, result(i) )
 
        end do

        write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

      end if

    end do

  end do
 
  deallocate ( result )

  return
end
subroutine test40 ( )

!*****************************************************************************80
!
!! TEST40 tests NCO_COMPUTE and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 9

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1
  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST40'
  write ( *, '(a)' ) '  NCO_COMPUTE computes an open Newton-Cotes rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max
 
      call nco_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do
 
  deallocate ( result )

  return
end
subroutine test401 ( )

!*****************************************************************************80
!
!! TEST401 tests NCOH_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 10

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1
  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST401'
  write ( *, '(a)' ) '  NCOH_SET sets up an open half Newton-Cotes rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)') '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max
 
      call ncoh_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do
 
  deallocate ( result )

  return
end
subroutine test402 ( )

!*****************************************************************************80
!
!! TEST402 tests NCOH_COMPUTE and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 10

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1
  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST402'
  write ( *, '(a)' ) '  NCOH_COMPUTE computes an open half Newton-Cotes rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max
 
      call ncoh_compute ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )
 
        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )
 
      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do

  end do
 
  deallocate ( result )

  return
end
subroutine test403 ( )

!*****************************************************************************80
!
!! TEST403 tests PATTERSON_SET and SUM_SUB.
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

  integer ( kind = 4 ), parameter :: order_max = 255

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) level
  integer ( kind = 4 ), parameter :: level_max = 8
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1

  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST403'
  write ( *, '(a)' ) '  PATTERSON_SET sets up a Gauss-Patterson rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(1x,a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do level = 1, level_max

      order = ( 2**level ) - 1

      call patterson_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, x, &
          w, result(i) )

      end do

      write ( *, '(i3,2x,5f14.8)' ) order, result(ilo:ihi)

    end do
 
  end do

  deallocate ( result )

  return
end
subroutine test404 ( )

!*****************************************************************************80
!
!! TEST404 tests RADAU_COMPUTE and RADAU_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real    ( kind = 8 ), allocatable, dimension ( : ) :: w1
  real    ( kind = 8 ), allocatable, dimension ( : ) :: w2
  real    ( kind = 8 ), allocatable, dimension ( : ) :: x1
  real    ( kind = 8 ), allocatable, dimension ( : ) :: x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST404'
  write ( *, '(a)' ) '  RADAU_COMPUTE computes a Radau rule;'
  write ( *, '(a)' ) '  RADAU_SET sets a rule from a table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         I      X1            X2            W1            W2'
 
  do n = 4, 12, 3

    allocate ( w1(1:n) )
    allocate ( w2(1:n) )
    allocate ( x1(1:n) )
    allocate ( x2(1:n) )

    call radau_compute ( n, x1, w1 )
    call radau_set ( n, x2, w2 )

    write ( * , '(a)' ) ' '
    do i = 1, n
      write (  *, '(2x,i8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8)' ) &
        i, x1(i), x2(i), w1(i), w2(i)
    end do

    deallocate ( w1 )
    deallocate ( w2 )
    deallocate ( x1 )
    deallocate ( x2 )

  end do
 
  return
end
subroutine test41 ( )

!*****************************************************************************80
!
!! TEST41 tests RADAU_SET and SUM_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 15

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  character ( len = 10 ) function_name
  integer ( kind = 4 ) function_num
  real    ( kind = 8 ), external :: function_value
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable, dimension ( : ) :: result
  real    ( kind = 8 ) w(order_max)
  real    ( kind = 8 ) x(order_max)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo

  call function_set ( 'COUNT', function_num )

  allocate ( result(1:function_num) )

  a = 0.0D+00
  b = 1.0D+00

  nsub = 1
  xlo = -1.0D+00
  xhi = +1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST41'
  write ( *, '(a)' ) '  RADAU_SET sets up a Radau rule;'
  write ( *, '(a)' ) '  SUM_SUB carries it out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.4,a,f10.4,a)' ) &
    '  The integration interval is [', a, ',', b, ']'
  write ( *, '(a,i8)' ) '  Number of subintervals is ', nsub
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '

  do ilo = 1, function_num, 5

    ihi = min ( ilo + 4, function_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a7, 5(a10,4x) )' ) &
      'Order  ', ( function_name ( i ), i = ilo, ihi )
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call radau_set ( order, x, w )

      do i = ilo, ihi

        call function_set ( 'SET', i )

        call sum_sub ( function_value, a, b, nsub, order, xlo, xhi, &
          x, w, result(i) )

      end do

      write ( *, '(i2,2x,5f14.8)' ) order, result(ilo:ihi)

    end do
 
  end do

  deallocate ( result )

  return
end
function f1sd1 ( x )

!*****************************************************************************80
!
!! F1SD1 evaluates the function 1.0D+00/ sqrt ( 1.1 - x**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) F1SD1, the value of the function.
!
  implicit none

  real    ( kind = 8 ) f1sd1
  real    ( kind = 8 ) x

  f1sd1 = 1.0D+00 / sqrt ( 1.1D+00 - x**2 )
 
  return
end
function fxsd1 ( x )

!*****************************************************************************80
!
!! FXSD1 evaluates the function x / sqrt ( 1.1 - x * x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FXSD1, the value of the function.
!
  implicit none

  real    ( kind = 8 ) fxsd1
  real    ( kind = 8 ) x

  fxsd1 = x / sqrt ( 1.1D+00 - x * x )
 
  return
end
function fx2sd1 ( x )

!*****************************************************************************80
!
!! FX2SD1 evaluates the function x**2 / sqrt ( 1.1 - x**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX2SD1, the value of the function.
!
  implicit none

  real    ( kind = 8 ) fx2sd1
  real    ( kind = 8 ) x

  fx2sd1 = x**2 / sqrt ( 1.1D+00 - x**2 )
 
  return
end
function function_value ( x )

!*****************************************************************************80
!
!! FUNCTION_VALUE evaluates a function of X, as chosen by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FUNCTION_VALUE, the value of the function.
!
  implicit none

  integer ( kind = 4 ) function_index
  real    ( kind = 8 ) function_value
  real    ( kind = 8 ) x

  call function_set ( 'GET', function_index )

  if ( function_index == 1 ) then
    function_value = 1.0D+00
  else if ( function_index == 2 ) then
    function_value = x
  else if ( function_index == 3 ) then
    function_value = x**2
  else if ( function_index == 4 ) then
    function_value = x**3
  else if ( function_index == 5 ) then
    function_value = x**4
  else if ( function_index == 6 ) then
    function_value = x**5
  else if ( function_index == 7 ) then
    function_value = x**6
  else if ( function_index == 8 ) then
    function_value = x**7
  else if ( function_index == 9 ) then
    function_value = sin ( x )
  else if ( function_index == 10 ) then
    function_value = exp ( x )
  else if ( function_index == 11 ) then
    function_value = sqrt ( abs ( x ) )
  else
    function_value = 0.0D+00
  end if

  return
end
subroutine function_set ( action, i )

!*****************************************************************************80
!
!! FUNCTION_SET sets the function to be returned by FUNCTION_VALUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action to be carried out.
!    'COUNT' means the call is made to count the number of functions available.
!    'GET' means the call is made to find out the current function index.
!    'SET' means the call is made to set the current function index.
!
!    Input/output, integer I.
!    For 'COUNT', I is output as the number of functions available;
!    For 'GET', I is output as the currently chosen function;
!    For 'SET', I is input as the user's new choice for the function.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: function_index = -1

  if ( action == 'COUNT' ) then
    i = 11
  else if ( action == 'GET' ) then
    i = function_index
  else if ( action == 'SET' ) then
    function_index = i
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNCTION_SET - Warning!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
  end if

  return
end
function function_name ( function_index )

!*****************************************************************************80
!
!! FUNCTION_NAME returns the name of the function evaluated in FUNCTION_VALUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FUNCTION_INDEX, the index of the function.
!
!    Output, character ( len = 10 ) FUNCTION_NAME, the name of the function.
!
  implicit none

  integer ( kind = 4 ) function_index
  character ( len = 10 ) function_name

  if ( function_index == 1 ) then
    function_name = '         1'
  else if ( function_index == 2 ) then
    function_name = '         X'
  else if ( function_index == 3 ) then
    function_name = '       X^2'
  else if ( function_index == 4 ) then
    function_name = '       X^3'
  else if ( function_index == 5 ) then
    function_name = '       X^4'
  else if ( function_index == 6 ) then
    function_name = '       X^5'
  else if ( function_index == 7 ) then
    function_name = '       X^6'
  else if ( function_index == 8 ) then
    function_name = '       X^7'
  else if ( function_index == 9 ) then
    function_name = '    SIN(X)'
  else if ( function_index == 10 ) then
    function_name = '    EXP(X)'
  else if ( function_index == 11 ) then
    function_name = ' SQRT(|X|)'
  else
    function_name = '??????????'
  end if

  return
end
