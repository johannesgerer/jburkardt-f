program main

!*****************************************************************************80
!
!! MAIN is the main program for NINTLIB_PRB.
!
!  Discussion:
!
!    NINTLIB_PRB runs the NINTLIB tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), dimension ( test_num ) :: dim_num_test = (/ 2, 3, 4 /)
  real ( kind = 8 ), external :: f1dn
  real ( kind = 8 ), external :: fbdn
  real ( kind = 8 ), external :: fedn
  real ( kind = 8 ), external :: fxdn
  real ( kind = 8 ), external :: fx2dn
  real ( kind = 8 ), external :: fx3dn
  integer ( kind = 4 ) test

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NINTLIB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NINTLIB library.'

  a = 0.0D+00
  b = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TESTND'
  write ( *, '(a)' ) '  Test routines for estimating the integral of'
  write ( *, '(a)' ) '  of F(X) in the hypercube [A,B]**DIM_NUM.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    dim_num = dim_num_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a,g20.12)' ) '  A(1:DIM_NUM) = ', a
    write ( *, '(a,g20.12)' ) '  B(1:DIM_NUM) = ', b

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  F(X(1:DIM_NUM)) = 1'
    write ( *, '(a)' ) ' '

    call testnd ( dim_num, f1dn )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  F(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM) )'
    write ( *, '(a)' ) ' '

    call testnd ( dim_num, fxdn )

    write ( *, '(a)') ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  F(X(1:DIM_NUM)) = sum( X(1:DIM_NUM)^2 )'
    write ( *, '(a)' ) ' '

    call testnd ( dim_num, fx2dn )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  F(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM)^3 )'
    write ( *, '(a)' ) ' '

    call testnd ( dim_num, fx3dn )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  F(X(1:DIM_NUM)) = exp(sum(X(1:DIM_NUM)))'
    write ( *, '(a)' ) ' '

    call testnd ( dim_num, fedn )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  F(X(1:DIM_NUM)) = 1/(1+sum(X(1:DIM_NUM)^2))'
    write ( *, '(a)' ) ' '

    call testnd ( dim_num, fbdn )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NINTLIB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine testnd ( dim_num, func )

!*****************************************************************************80
!
!! TESTND tests the integrators on a particular function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function
!    to be integrated.
!
  implicit none

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) dim_num

  call test01 ( dim_num, func )
  call test02 ( dim_num, func )
  call test03 ( dim_num, func )
  call test04 ( dim_num, func )
!
!  TEST05 is only set up for DIM_NUM = 2.
!
  if ( dim_num == 2 ) then
    call test05 ( dim_num, func )
  end if

  call test06 ( dim_num, func )

  return
end
subroutine test01 ( dim_num, func )

!*****************************************************************************80
!
!! TEST01 tests BOX_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function
!    to be integrated.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), parameter :: order = 5

  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) result
  real ( kind = 8 ), dimension ( order ) :: wtab = (/ &
    0.236926885056189087514264040720D+00, &
    0.478628670499366468041291514836D+00, &
    0.568888888888888888888888888889D+00, &
    0.478628670499366468041291514836D+00, &
    0.236926885056189087514264040720D+00 /)
  real ( kind = 8 ), dimension ( order ) :: wtab2
  real ( kind = 8 ), dimension ( order ) :: xtab = (/ &
    -0.906179845938663992797626878299D+00, &
    -0.538469310105683091036314420700D+00, &
     0.0D+00, &
     0.538469310105683091036314420700D+00, &
     0.906179845938663992797626878299D+00 /)
  real ( kind = 8 ), dimension ( order ) :: xtab2
!
!  Adjust the quadrature rule from [-1,1] to [0,1]:
!
  xtab2(1:order) = ( xtab(1:order) + 1.0D+00 ) / 2.0D+00
  wtab2(1:order) = 0.5D+00 * wtab(1:order)

  call box_nd ( func, dim_num, order, xtab2, wtab2, result, eval_num )

  write ( *, '(a,g20.12,2x,i8)' ) '  BOX_ND:         ', result, eval_num

  return
end
subroutine test02 ( dim_num, func )

!*****************************************************************************80
!
!! TEST02 tests P5_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function
!    to be integrated.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) result
!
!  Set the integration limits.
!
  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  call p5_nd ( func, dim_num, a, b, result, eval_num )

  write ( *, '(a,g20.12,2x,i8)' ) '  P5_ND:          ', result, eval_num

  return
end
subroutine test03 ( dim_num, func )

!*****************************************************************************80
!
!! TEST03 tests ROMBERG_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function
!    to be integrated.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) ind
  integer ( kind = 4 ), parameter :: it_max = 3
  real ( kind = 8 ) result
  integer ( kind = 4 ) sub_num(dim_num)
  real ( kind = 8 ) tol
!
!  Set the integration limits.
!
  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  tol = 0.001D+00

  sub_num(1:dim_num) = 10

  call romberg_nd ( func, a, b, dim_num, sub_num, it_max, tol, &
    result, ind, eval_num )

  write ( *, '(a,g20.12,2x,i8)' ) '  ROMBERG_ND:     ', result, eval_num

  return
end
subroutine test04 ( dim_num, func )

!*****************************************************************************80
!
!! TEST04 tests SAMPLE_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function
!    to be integrated.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ), parameter :: k2 = 4

  real ( kind = 8 ) dev1(k2)
  real ( kind = 8 ) dev2(k2)
  real ( kind = 8 ) err1(k2)
  real ( kind = 8 ) est1(k2)
  real ( kind = 8 ) est2(k2)
  real ( kind = 8 ) err2(k2)
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) k1

  k1 = 1

  call sample_nd ( func, k1, k2, dim_num, est1, err1, dev1, est2, err2, &
    dev2, eval_num )

  write ( *, '(a,g20.12,2x,i8)' ) '  SAMPLE_ND:      ', est2(k2), eval_num

  return
end
subroutine test05 ( dim_num, func )

!*****************************************************************************80
!
!! TEST05 demonstrates how to refine N-dimensional integration results.
!
!  Discussion:
!
!    This routine is only set up for DIM_NUM = 2 for now.
!
!    We are given a routine, NDP5, which will integrate over a
!    DIM_NUM dimensional hypercube using a fixed method.  In order to
!    improve the approximation to an integral, we can subdivide
!    the hypercube and call NDP5 to integrate again over each of
!    these regions.
!
!    The information that we gather can be used to tell us when
!    to expect that we have achieved a certain degree of accuracy.
!
!    With a little more work, we could make this code adaptive.
!    That is, it would only refine SOME of the subregions, where
!    the approximation to the integral was still not good enough.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function
!    to be integrated.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) eval_num
  integer ( kind = 4 ) eval_total
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ngrid
  real ( kind = 8 ) result
  real ( kind = 8 ) result_total
  real ( kind = 8 ) xlo(dim_num)
  real ( kind = 8 ) xhi(dim_num)

  xlo(1:dim_num) = 0.0D+00
  xhi(1:dim_num) = 1.0D+00

  do igrid = 1, 6

    ngrid = 2**( igrid - 1 )

    result_total = 0.0D+00
    eval_total = 0

    do i = 1, ngrid

      a(1) = ( real ( ngrid - i + 1, kind = 8 ) * xlo(1)   &
             + real (         i - 1, kind = 8 ) * xhi(1) ) &
             / real ( ngrid,         kind = 8 )

      b(1) = ( real ( ngrid - i, kind = 8 ) * xlo(1)   &
             + real (         i, kind = 8 ) * xhi(1) ) &
             / real ( ngrid,     kind = 8 )

      do j = 1, ngrid

        a(2) = ( real ( ngrid - j + 1, kind = 8 ) * xlo(2)   &
               + real (         j - 1, kind = 8 ) * xhi(2) ) &
               / real ( ngrid,         kind = 8 )

        b(2) = ( real ( ngrid - j, kind = 8 ) * xlo(2)   &
               + real (         j, kind = 8 ) * xhi(2) ) &
               / real ( ngrid,     kind = 8 )

        call p5_nd ( func, dim_num, a, b, result, eval_num )

        result_total = result_total + result
        eval_total = eval_total + eval_num

      end do

    end do

    write ( *, '(a,g20.12,2x,i8)' ) &
      '  P5_ND+:         ', result_total, eval_total

  end do

  return
end
subroutine test06 ( dim_num, func )

!*****************************************************************************80
!
!! TEST06 tests MONTE_CARLO_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function
!    to be integrated.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test

  seed = 123456789
!
!  Set the integration limits.
!
  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  do test = 1, test_num

    eval_num = 8**test * 10000

    call monte_carlo_nd ( func, dim_num, a, b, eval_num, seed, result )

    write ( *, '(a,g20.12,2x,i8)' ) '  MONTE_CARLO_ND: ', result, eval_num

  end do

  return
end
function f1dn ( dim_num, x )

!*****************************************************************************80
!
!! F1DN(X(1:DIM_NUM)) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the argument.
!
!    Output, real ( kind = 8 ) F1DN, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) f1dn
  real ( kind = 8 ) x(dim_num)

  f1dn = 1.0D+00

  return
end
function fbdn ( dim_num, x )

!*****************************************************************************80
!
!! FBDN(X(1:DIM_NUM)) = 1 / ( 1 + sum ( X(1:DIM_NUM)**2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the argument.
!
!    Output, real ( kind = 8 ) FBDN, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) fbdn
  real ( kind = 8 ) x(dim_num)

  fbdn = 1.0D+00 / ( 1.0D+00 + sum ( x(1:dim_num)**2 ) )

  return
end
function fedn ( dim_num, x )

!*****************************************************************************80
!
!! FEDN(X(1:DIM_NUM)) = EXP ( sum ( X(1:DIM_NUM) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the argument.
!
!    Output, real ( kind = 8 ) FEDN, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) fedn
  real ( kind = 8 ) x(dim_num)

  fedn = exp ( sum ( x(1:dim_num) ) )

  return
end
function fxdn ( dim_num, x )

!*****************************************************************************80
!
!! FXDN(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the argument.
!
!    Output, real ( kind = 8 ) FXDN, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) fxdn
  real ( kind = 8 ) x(dim_num)

  fxdn = sum ( x(1:dim_num) )

  return
end
function fx2dn ( dim_num, x )

!*****************************************************************************80
!
!! FX2DN(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM)**2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the argument.
!
!    Output, real ( kind = 8 ) FX2DN, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) fx2dn
  real ( kind = 8 ) x(dim_num)

  fx2dn = sum ( x(1:dim_num)**2 )

  return
end
function fx3dn ( dim_num, x )

!*****************************************************************************80
!
!! FX3DN(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM)**3 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the argument.
!
!    Output, real ( kind = 8 ) FX3DN, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) fx3dn
  real ( kind = 8 ) x(dim_num)

  fx3dn = sum ( x(1:dim_num)**3 )

  return
end
