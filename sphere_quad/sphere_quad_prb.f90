program main

!*****************************************************************************80
!
!! MAIN is the main program for SPHERE_QUAD_PRB.
!
!  Discussion:
!
!    SPHERE_QUAD_PRB tests SPHERE_QUAD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_QUAD_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPHERE_QUAD library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_QUAD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SPHERE01_QUAD_LL*.
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

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: polyterm_value_3d
  real ( kind = 8 ) h
  integer ( kind = 4 ) h_test
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n_llc
  integer ( kind = 4 ) n_llm
  integer ( kind = 4 ) n_llv
  integer ( kind = 4 ) n_mc
  real ( kind = 8 ) result_llc
  real ( kind = 8 ) result_llm
  real ( kind = 8 ) result_llv
  real ( kind = 8 ) result_mc
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Approximate the integral of a function on the unit sphere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SPHERE01_QUAD_MC uses a Monte Carlo method.'
  write ( *, '(a)' ) '  SPHERE01_QUAD_LLC uses centroids of spherical triangles.'
  write ( *, '(a)' ) '  SPHERE01_QUAD_LLM uses midsides of spherical triangles.'
  write ( *, '(a)' ) '  SPHERE01_QUAD_LLV uses vertices of spherical triangles.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  H              QUAD_MC       QUAD_LLC      QUAD_LLM      QUAD_LLV         EXACT'

  do i = 0, 17

    if ( i == 0 ) then
      e(1:3) = (/ 0, 0, 0 /)
    else if ( i == 1 ) then
      e(1:3) = (/ 0, 0, 0 /)
    else if ( i == 2 ) then
      e(1:3) = (/ 1, 0, 0 /)
    else if ( i == 3 ) then
      e(1:3) = (/ 0, 1, 0 /)
    else if ( i == 4 ) then
      e(1:3) = (/ 0, 0, 1 /)
    else if ( i == 5 ) then
      e(1:3) = (/ 2, 0, 0 /)
    else if ( i == 6 ) then
      e(1:3) = (/ 0, 2, 2 /)
    else if ( i == 7 ) then
      e(1:3) = (/ 2, 2, 2 /)
    else if ( i == 8 ) then
      e(1:3) = (/ 0, 2, 4 /)
    else if ( i == 9 ) then
      e(1:3) = (/ 0, 0, 6 /)
    else if ( i == 10 ) then
      e(1:3) = (/ 1, 2, 4 /)
    else if ( i == 11 ) then
      e(1:3) = (/ 2, 4, 2 /)
    else if ( i == 12 ) then
      e(1:3) = (/ 6, 2, 0 /)
    else if ( i == 13 ) then
      e(1:3) = (/ 0, 0, 8 /)
    else if ( i == 14 ) then
      e(1:3) = (/ 6, 0, 4 /)
    else if ( i == 15 ) then
      e(1:3) = (/ 4, 6, 2 /)
    else if ( i == 16 ) then
      e(1:3) = (/ 2, 4, 8 /)
    else if ( i == 17 ) then
      e(1:3) = (/ 16, 0, 0 /)
    end if

    call polyterm_exponent ( 'SET', e )

    if ( i == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Point counts per method:'
    else
      call polyterm_exponent ( 'PRINT', e )
    end if

    do h_test = 1, 3

      if ( h_test == 1 ) then
        h = 1.0D+00
      else if (  h_test == 2 ) then
        h = 0.1D+00
      else if (  h_test == 3 ) then
        h = 0.01D+00
      end if

      call sphere01_quad_mc_size ( h, n_mc )

      call sphere01_quad_mc ( polyterm_value_3d, h, seed, n_mc, result_mc )

      call sphere01_quad_llc ( polyterm_value_3d, h, n_llc, result_llc )

      call sphere01_quad_llm ( polyterm_value_3d, h, n_llm, result_llm )

      call sphere01_quad_llv ( polyterm_value_3d, h, n_llv, result_llv )

      call sphere01_monomial_integral ( e, exact )

      if ( i == 0 ) then
        write ( *, '(g14.6,5i14)' ) h, n_mc, n_llc, n_llm, n_llv
      else
        write ( *, '(6g14.6)' ) &
          h, result_mc, result_llc, result_llm, result_llv, exact

      end if

    end do

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests SPHERE01_QUAD_ICOS1C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: polyterm_value_3d
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) factor_log
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Approximate the integral of a function on the unit sphere.'
  write ( *, '(a)' ) '  SPHERE01_QUAD_ICOS1C uses centroids of spherical triangles.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'FACTOR         N        QUAD          EXACT         ERROR'

  do i = 1, 17

    if ( i == 1 ) then
      e(1:3) = (/ 0, 0, 0 /)
    else if ( i == 2 ) then
      e(1:3) = (/ 1, 0, 0 /)
    else if ( i == 3 ) then
      e(1:3) = (/ 0, 1, 0 /)
    else if ( i == 4 ) then
      e(1:3) = (/ 0, 0, 1 /)
    else if ( i == 5 ) then
      e(1:3) = (/ 2, 0, 0 /)
    else if ( i == 6 ) then
      e(1:3) = (/ 0, 2, 2 /)
    else if ( i == 7 ) then
      e(1:3) = (/ 2, 2, 2 /)
    else if ( i == 8 ) then
      e(1:3) = (/ 0, 2, 4 /)
    else if ( i == 9 ) then
      e(1:3) = (/ 0, 0, 6 /)
    else if ( i == 10 ) then
      e(1:3) = (/ 1, 2, 4 /)
    else if ( i == 11 ) then
      e(1:3) = (/ 2, 4, 2 /)
    else if ( i == 12 ) then
      e(1:3) = (/ 6, 2, 0 /)
    else if ( i == 13 ) then
      e(1:3) = (/ 0, 0, 8 /)
    else if ( i == 14 ) then
      e(1:3) = (/ 6, 0, 4 /)
    else if ( i == 15 ) then
      e(1:3) = (/ 4, 6, 2 /)
    else if ( i == 16 ) then
      e(1:3) = (/ 2, 4, 8 /)
    else if ( i == 17 ) then
      e(1:3) = (/ 16, 0, 0 /)
    end if

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    factor = 1
    do factor_log = 0, 5

      call sphere01_quad_icos1c ( factor, polyterm_value_3d, n, result )

      call sphere01_monomial_integral ( e, exact )

      error = abs ( exact - result )

      write ( *, '(2x,i4,2x,i8,2x,3g14.6)' ) &
          factor, n, result, exact, error

      factor = factor * 2

    end do

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests SPHERE01_QUAD_ICOS1M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: polyterm_value_3d
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) factor_log
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Approximate the integral of a function on the unit sphere.'
  write ( *, '(a)' ) '  SPHERE01_QUAD_ICOS1M uses midsides of spherical triangles.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'FACTOR         N        QUAD          EXACT         ERROR'

  do i = 1, 17
 
    if ( i == 1 ) then
      e(1:3) = (/ 0, 0, 0 /)
    else if ( i == 2 ) then
      e(1:3) = (/ 1, 0, 0 /)
    else if ( i == 3 ) then
      e(1:3) = (/ 0, 1, 0 /)
    else if ( i == 4 ) then
      e(1:3) = (/ 0, 0, 1 /)
    else if ( i == 5 ) then
      e(1:3) = (/ 2, 0, 0 /)
    else if ( i == 6 ) then
      e(1:3) = (/ 0, 2, 2 /)
    else if ( i == 7 ) then
      e(1:3) = (/ 2, 2, 2 /)
    else if ( i == 8 ) then
      e(1:3) = (/ 0, 2, 4 /)
    else if ( i == 9 ) then
      e(1:3) = (/ 0, 0, 6 /)
    else if ( i == 10 ) then
      e(1:3) = (/ 1, 2, 4 /)
    else if ( i == 11 ) then
      e(1:3) = (/ 2, 4, 2 /)
    else if ( i == 12 ) then
      e(1:3) = (/ 6, 2, 0 /)
    else if ( i == 13 ) then
      e(1:3) = (/ 0, 0, 8 /)
    else if ( i == 14 ) then
      e(1:3) = (/ 6, 0, 4 /)
    else if ( i == 15 ) then
      e(1:3) = (/ 4, 6, 2 /)
    else if ( i == 16 ) then
      e(1:3) = (/ 2, 4, 8 /)
    else if ( i == 17 ) then
      e(1:3) = (/ 16, 0, 0 /)
    end if

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    factor = 1
    do factor_log = 0, 5

      call sphere01_quad_icos1m ( factor, polyterm_value_3d, n, result )

      call sphere01_monomial_integral ( e, exact )

      error = abs ( exact - result )

      write ( *, '(2x,i4,2x,i8,2x,3g14.6)' ) &
          factor, n, result, exact, error

      factor = factor * 2

    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests SPHERE01_QUAD_ICOS1V.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: polyterm_value_3d
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) factor_log
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Approximate the integral of a function on the unit sphere.'
  write ( *, '(a)' ) '  SPHERE01_QUAD_ICOS1V uses vertices of spherical triangles.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'FACTOR         N        QUAD          EXACT         ERROR'

  do i = 1, 17
 
    if ( i == 1 ) then
      e(1:3) = (/ 0, 0, 0 /)
    else if ( i == 2 ) then
      e(1:3) = (/ 1, 0, 0 /)
    else if ( i == 3 ) then
      e(1:3) = (/ 0, 1, 0 /)
    else if ( i == 4 ) then
      e(1:3) = (/ 0, 0, 1 /)
    else if ( i == 5 ) then
      e(1:3) = (/ 2, 0, 0 /)
    else if ( i == 6 ) then
      e(1:3) = (/ 0, 2, 2 /)
    else if ( i == 7 ) then
      e(1:3) = (/ 2, 2, 2 /)
    else if ( i == 8 ) then
      e(1:3) = (/ 0, 2, 4 /)
    else if ( i == 9 ) then
      e(1:3) = (/ 0, 0, 6 /)
    else if ( i == 10 ) then
      e(1:3) = (/ 1, 2, 4 /)
    else if ( i == 11 ) then
      e(1:3) = (/ 2, 4, 2 /)
    else if ( i == 12 ) then
      e(1:3) = (/ 6, 2, 0 /)
    else if ( i == 13 ) then
      e(1:3) = (/ 0, 0, 8 /)
    else if ( i == 14 ) then
      e(1:3) = (/ 6, 0, 4 /)
    else if ( i == 15 ) then
      e(1:3) = (/ 4, 6, 2 /)
    else if ( i == 16 ) then
      e(1:3) = (/ 2, 4, 8 /)
    else if ( i == 17 ) then
      e(1:3) = (/ 16, 0, 0 /)
    end if

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    factor = 1
    do factor_log = 0, 5

      call sphere01_quad_icos1v ( factor, polyterm_value_3d, n, result )

      call sphere01_monomial_integral ( e, exact )

      error = abs ( exact - result )

      write ( *, '(2x,i4,2x,i8,2x,3g14.6)' ) &
          factor, n, result, exact, error

      factor = factor * 2

    end do

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests SPHERE01_QUAD_ICOS2V.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: polyterm_value_3d
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) factor_log
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Approximate the integral of a function on the unit sphere.'
  write ( *, '(a)' ) '  SPHERE01_QUAD_ICOS2V uses vertices of spherical triangles.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'FACTOR         N        QUAD          EXACT         ERROR'

  do i = 1, 17

    if ( i == 1 ) then
      e(1:3) = (/ 0, 0, 0 /)
    else if ( i == 2 ) then
      e(1:3) = (/ 1, 0, 0 /)
    else if ( i == 3 ) then
      e(1:3) = (/ 0, 1, 0 /)
    else if ( i == 4 ) then
      e(1:3) = (/ 0, 0, 1 /)
    else if ( i == 5 ) then
      e(1:3) = (/ 2, 0, 0 /)
    else if ( i == 6 ) then
      e(1:3) = (/ 0, 2, 2 /)
    else if ( i == 7 ) then
      e(1:3) = (/ 2, 2, 2 /)
    else if ( i == 8 ) then
      e(1:3) = (/ 0, 2, 4 /)
    else if ( i == 9 ) then
      e(1:3) = (/ 0, 0, 6 /)
    else if ( i == 10 ) then
      e(1:3) = (/ 1, 2, 4 /)
    else if ( i == 11 ) then
      e(1:3) = (/ 2, 4, 2 /)
    else if ( i == 12 ) then
      e(1:3) = (/ 6, 2, 0 /)
    else if ( i == 13 ) then
      e(1:3) = (/ 0, 0, 8 /)
    else if ( i == 14 ) then
      e(1:3) = (/ 6, 0, 4 /)
    else if ( i == 15 ) then
      e(1:3) = (/ 4, 6, 2 /)
    else if ( i == 16 ) then
      e(1:3) = (/ 2, 4, 8 /)
    else if ( i == 17 ) then
      e(1:3) = (/ 16, 0, 0 /)
    end if

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    factor = 1
    do factor_log = 0, 5

      call sphere01_quad_icos2v ( factor, polyterm_value_3d, n, result )

      call sphere01_monomial_integral ( e, exact )

      error = abs ( exact - result )

      write ( *, '(2x,i4,2x,i8,2x,3g14.6)' ) &
          factor, n, result, exact, error

      factor = factor * 2

    end do

  end do

  return
end
subroutine polyterm_exponent ( action, e )

!*****************************************************************************80
!
!! POLYTERM_EXPONENT gets or sets the exponents for the polynomial term.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 3 ) ACTION.
!    'GET' asks the routine to return the current values in E.
!    'SET' asks the routine to set the current values to E.
!
!    Input/output, integer ( kind = 4 ) E(3), storage used to set or get values.
!
  implicit none

  character ( len = * )  action
  integer   ( kind = 4 ) e(3)
  integer   ( kind = 4 ), save, dimension ( 3 ) :: e_save = (/ 0, 0, 0 /)
  character ( len = 80 ) text
  character ( len = 80 ) text2

  if ( action(1:1) == 'G' ) then

    e(1:3) = e_save(1:3)

  else if ( action(1:1) == 'P' ) then

    write ( *, '(a)' ) ' '

    if ( all ( e_save(1:3) == 0 ) ) then

      text = 'P(X,Y,Z) = 1'

    else

      text = 'P(X,Y,Z) = '

      if ( e_save(1) == 0 ) then

      else if ( e_save(1) == 1 ) then

        call s_cat ( text, ' X', text )

      else

        call s_cat ( text, ' X^', text )

        write ( text2, '(i2)' ) e_save(1)
        text2 = adjustl ( text2 )
        call s_cat ( text, text2, text )

      end if

      if ( e_save(2) == 0 ) then

      else if ( e_save(2) == 1 ) then

        call s_cat ( text, ' Y', text )

      else

        call s_cat ( text, ' Y^', text )

        write ( text2, '(i2)' ) e_save(2)
        text2 = adjustl ( text2 )
        call s_cat ( text, text2, text )

      end if
       
      if ( e_save(3) == 0 ) then

      else if ( e_save(3) == 1 ) then

        call s_cat ( text, ' Z', text )

      else

        call s_cat ( text, ' Z^', text )

        write ( text2, '(i2)' ) e_save(3)
        text2 = adjustl ( text2 )
        call s_cat ( text, text2, text )

      end if
 
    end if

    write ( *, '(a)' ) trim ( text )
    
  else if ( action(1:1) == 'S' ) then

    e_save(1:3) = e(1:3)

  end if

  return
end
subroutine polyterm_value_3d ( n, x, f )

!*****************************************************************************80
!
!! POLYTERM_VALUE_3D evaluates a single polynomial term in 3D.
!
!  Discussion:
!
!    The polynomial term has the form:
!
!      F(X) = X(1)^E(1) * X(2)^E(2) * X(3)^E(3)
!
!    The exponents E(1:3) are set by calling POLYTERM_EXPONENT.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(3,N), the points where the polynomial term 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the polynomial term.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(3,n)

  call polyterm_exponent ( 'GET', e )

  f(1:n) = 1.0D+00

  do i = 1, 3

    if ( e(i) /= 0 ) then
      f(1:n) = f(1:n) * x(i,1:n)**e(i)
    end if

  end do
  
  return
end
