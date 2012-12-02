program main

!*****************************************************************************80
!
!! MAIN is the main program for SPHERE_DESIGN_RULE_PRB.
!
!  Discussion:
!
!    DESIGN_PRB tests the DESIGN routines for integrals on a sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_DESIGN_RULE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPHERE_DESIGN_RULE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_DESIGN_RULE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DESIGN_QUAD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ), external :: polyterm_value_3d
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_max
  real ( kind = 8 ) quad
  real ( kind = 8 ), parameter :: r = 1.0D+00
  real ( kind = 8 ) sphere_area
  real ( kind = 8 ), dimension ( 3 ) :: xc = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  DESIGN_QUAD can return the average value of a'
  write ( *, '(a)' ) '  function F(X,Y,Z) at the points of a spherical'
  write ( *, '(a)' ) '  design.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this test, we will use single polynomial terms.'

  call sphere_area_3d ( r, sphere_area )

  call design_order_max ( order_max )

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

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' Order   Quad          Integral'
    write ( *, '(a)' ) ' '

    do order = 0, order_max

      if ( order == 0 ) then

        call sphere_monomial_int_3d ( r, e, integral )

        write ( *, '(a6,14x,g14.6)' ) 'Exact:', integral

      else

        call design_quad ( xc, r, polyterm_value_3d, order, quad )
 
        integral = quad * sphere_area

        write ( *, '(i6,2g14.6)' ) order, quad, integral

      end if

    end do
 
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DESIGN_QUAD on a shifted scaled sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ), external :: polyterm_value_3d
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_max
  real ( kind = 8 ) quad
  real ( kind = 8 ), parameter :: r = 2.0D+00
  real ( kind = 8 ) sphere_area
  real ( kind = 8 ), dimension ( 3 ) :: xc = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  DESIGN_QUAD can return the average value of a'
  write ( *, '(a)' ) '  function F(X,Y,Z) at the points of a spherical'
  write ( *, '(a)' ) '  design.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this test, we will use single polynomial terms.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The sphere will have a non-unit radius of R = ', r
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sphere will have a nonzero center of:'
  write ( *, '(3g14.6)' ) xc(1:3)

  call sphere_area_3d ( r, sphere_area )

  call design_order_max ( order_max )

  do i = 1, 15

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
      e(1:3) = (/ 1, 2, 4 /)
    else if ( i == 10 ) then
      e(1:3) = (/ 2, 4, 2 /)
    else if ( i == 11 ) then
      e(1:3) = (/ 6, 2, 0 /)
    else if ( i == 12 ) then
      e(1:3) = (/ 0, 0, 8 /)
    else if ( i == 13 ) then
      e(1:3) = (/ 6, 0, 4 /)
    else if ( i == 14 ) then
      e(1:3) = (/ 4, 6, 2 /)
    else if ( i == 15 ) then
      e(1:3) = (/ 2, 4, 8 /)
    end if

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' Order   Quad          Integral'
    write ( *, '(a)' ) ' '

    do order = 1, order_max

      call design_quad ( xc, r, polyterm_value_3d, order, quad )

      integral = quad * sphere_area

      write ( *, '(i6,2g14.6)' ) order, quad, integral

    end do
 
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 writes a sphere design rule to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) design_filename
  real      ( kind = 8 ) x(3,180)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  R8MAT_WRITE can write a sphere design rule to a file.'

  call design_18_180_3d ( x )

  design_filename = 'sphere_design_rule_18.txt'

  call r8mat_write ( design_filename, 3, 180, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sphere design rule of 180 points written to "' &
    // trim ( design_filename ) // '".'

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

    write ( *, * ) ' '

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
function polyterm_value_3d ( x )

!*****************************************************************************80
!
!! POLYTERM_VALUE_3D evaluates a single polynomial term in 3D.
!
!  Discussion:
!
!    The polynomial term has the form:
!
!      X(1)**E(1) * X(2)**E(2) * X(3)**E(3)
!
!    The exponents E(1:3) are set by calling POLYTERM_EXPONENT_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(3), the point where the polynomial 
!    term is to be evaluated.
!
!    Output, real ( kind = 8 ) POLYTERM_VALUE_3D, the value of the 
!    polynomial term.
!
  implicit none

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ) polyterm_value_3d
  real ( kind = 8 ) x(3)

  call polyterm_exponent ( 'GET', e )

  if ( all ( e(1:3) == 0 )  ) then
    polyterm_value_3d = 1
  else if ( any ( x(1:3) == 0.0D+00 ) ) then
    polyterm_value_3d = 0.0D+00
  else
    polyterm_value_3d = x(1)**e(1) * x(2)**e(2) * x(3)**e(3)
  end if
  
  return
end

