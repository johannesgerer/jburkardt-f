function c8_abs ( z )

!*****************************************************************************80
!
!! C8_ABS evaluates the absolute value of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the absolute value of a C8 with the ABS function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, real ( kind = 8 ) C8_ABS, the function value.
!
  implicit none

  real ( kind = 8 ) c8_abs
  complex ( kind = 8 ) z

  c8_abs = sqrt ( ( real ( z, kind = 8 ) )**2 &
                + ( aimag ( z ) )**2 )

  return
end
function c8_acos ( z )

!*****************************************************************************80
!
!! C8_ACOS evaluates the inverse cosine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the inverse cosine of a C8.
!
!    Here we use the relationship:
!
!       C8_ACOS ( Z ) = pi/2 - C8_ASIN ( Z ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_ACOS, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_acos
  complex ( kind = 8 ) c8_asin
  real ( kind = 8 ), parameter :: r8_pi_half = 1.57079632679489661923D+00
  complex ( kind = 8 ) z

  c8_acos = r8_pi_half - c8_asin ( z )

  return
end
function c8_acosh ( z )

!*****************************************************************************80
!
!! C8_ACOSH evaluates the inverse hyperbolic cosine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the inverse hyperbolic cosine of a C8.
!
!    Here we use the relationship:
!
!      C8_ACOSH ( Z ) = i * C8_ACOS ( Z ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_ACOSH, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_acos
  complex ( kind = 8 ) c8_acosh
  complex ( kind = 8 ) c8_i
  complex ( kind = 8 ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )

  c8_acosh = c8_i * c8_acos ( z )

  return
end
function c8_add ( z1, z2 )

!*****************************************************************************80
!
!! C8_ADD adds two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports addition of C8's with the "+" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z1, Z2, the values to add.
!
!    Output, complex ( kind = 8 ) C8_ADD, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_add
  complex ( kind = 8 ) z1
  complex ( kind = 8 ) z2

  c8_add = z1 + z2

  return
end
function c8_arg ( x )

!*****************************************************************************80
!
!! C8_ARG returns the argument of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the argument of a C8.
!
!    By convention, the argument of a C8 is expected to lie between
!    -PI and PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the complex number.
!
!    Output, real ( kind = 8 ) C8_ARG, the function value.
!
  implicit none

  real ( kind = 8 ) c8_arg
  complex ( kind = 8 ) x

  if ( aimag ( x )           == 0.0D+00 .and. &
       real  ( x, kind = 8 ) == 0.0D+00 ) then

    c8_arg = 0.0D+00

  else

    c8_arg = atan2 ( aimag ( x ), real ( x ) )

  end if

  return
end
function c8_asin ( z )

!*****************************************************************************80
!
!! C8_ASIN evaluates the inverse sine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the inverse sine of a C8.
!
!    Here we use the relationship:
!
!      C8_ASIN ( Z ) = - i * log ( i * z + sqrt ( 1 - z * z ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_ASIN, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_asin
  complex ( kind = 8 ) c8_i
  complex ( kind = 8 ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )

  c8_asin = - c8_i * log ( c8_i * z + sqrt ( 1.0D+00 - z * z ) )

  return
end
function c8_asinh ( z )

!*****************************************************************************80
!
!! C8_ASINH evaluates the inverse hyperbolic sine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the inverse hyperbolic sine of a C8.
!
!    Here we use the relationship:
!
!      C8_ASINH ( Z ) = - i * C8_ASIN ( i * Z ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_ASINH, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_asin
  complex ( kind = 8 ) c8_asinh
  complex ( kind = 8 ) c8_i
  complex ( kind = 8 ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )

  c8_asinh = - c8_i * c8_asin ( c8_i * z )

  return
end
function c8_atan ( z )

!*****************************************************************************80
!
!! C8_ATAN evaluates the inverse tangent of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the inverse tangent of a C8.
!
!    Here we use the relationship:
!
!      C8_ATAN ( Z ) = ( i / 2 ) * log ( ( 1 - i * z ) / ( 1 + i * z ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_ATAN, the function value.
!
  implicit none

  complex ( kind = 8 ) arg
  complex ( kind = 8 ) c8_atan
  complex ( kind = 8 ) c8_i
  complex ( kind = 8 ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )

  arg = ( 1.0D+00 - c8_i * z ) / ( 1.0D+00 + c8_i * z )

  c8_atan = 0.5D+00 * c8_i * log ( arg )

  return
end
function c8_atanh ( z )

!*****************************************************************************80
!
!! C8_ATANH evaluates the inverse hyperbolic tangent of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the inverse hyperbolic tangent of a C8.
!
!    Here we use the relationship:
!
!      C8_ATANH ( Z ) = - i * C8_ATAN ( i * Z ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_ATANH, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_atan
  complex ( kind = 8 ) c8_atanh
  complex ( kind = 8 ) c8_i
  complex ( kind = 8 ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )

  c8_atanh = - c8_i * c8_atan ( c8_i * z )

  return
end
function c8_conj ( z )

!*****************************************************************************80
!
!! C8_CONJ evaluates the conjugate of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the conjugate of a C8 with the CONJG function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_CONJ, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_conj
  complex ( kind = 8 ) z

  c8_conj = cmplx ( real ( z, kind = 8 ), - aimag ( z ) )

  return
end
function c8_copy ( z )

!*****************************************************************************80
!
!! C8_COPY copies a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the copy of a C8 with the "=" operator.
!
!    The order of the arguments may seem unnatural, but it is arranged so
!    that the call
!
!      c8_copy ( c1, c2 )
!
!    mimics the assignment
!
!      c1 = c2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_COPY, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_copy
  complex ( kind = 8 ) z

  c8_copy = z

  return
end
function c8_cos ( z )

!*****************************************************************************80
!
!! C8_COS evaluates the cosine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the cosine of a C8 with the COS function.
!
!    We use the relationship:
!
!      C8_COS ( C ) = ( C8_EXP ( i * C ) + C8_EXP ( - i * C ) ) / 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_COS, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_cos
  complex ( kind = 8 ) c8_exp
  complex ( kind = 8 ) c8_i
  complex ( kind = 8 ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00 )

  c8_cos = ( c8_exp ( c8_i * z ) + c8_exp ( - c8_i * z ) ) / 2.0D+00

  return
end
function c8_cosh ( z )

!*****************************************************************************80
!
!! C8_COSH evaluates the hyperbolic cosine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the hyperbolic cosine of a C8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_COSH, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_cosh
  complex ( kind = 8 ) c8_exp
  complex ( kind = 8 ) z

  c8_cosh = ( c8_exp ( z ) + c8_exp ( - z ) ) / 2.0D+00

  return
end
function c8_cube_root ( x )

!*****************************************************************************80
!
!! C8_CUBE_ROOT returns the principal cube root of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the cube root of a C8 through the "**(1.0/3.0)" operator.
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
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the argument.
!
!    Output, complex ( kind = 8 ) C8_CUBE_ROOT, the function value.
!
  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) c8_arg
  complex ( kind = 8 ) c8_cube_root
  real ( kind = 8 ) c8_mag
  real ( kind = 8 ) mag
  complex ( kind = 8 ) x

  arg = c8_arg ( x )
  mag = c8_mag ( x )

  if ( mag == 0.0D+00 ) then

    c8_cube_root = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  else

    c8_cube_root = mag**( 1.0D+00 / 3.0D+00 ) &
      * cmplx ( cos ( arg / 3.0D+00 ), &
                sin ( arg / 3.0D+00 ), kind = 8 )

  end if

  return
end
function c8_div ( z1, z2 )

!*****************************************************************************80
!
!! C8_DIV divides two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports division of C8's with the "/" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z1, Z2, the arguments.
!
!    Output, complex ( kind = 8 ) C8_DIV, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_div
  complex ( kind = 8 ) z1
  complex ( kind = 8 ) z2

  c8_div = z1 / z2

  return
end
function c8_div_r8 ( z1, r )

!*****************************************************************************80
!
!! C8_DIV_R8 divides a C8 by an R8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!    An R8 is a real ( kind = 8 ) value.
!
!    FORTRAN90 supports division of a C8 by an R8 with the "/" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z1, the value to be divided.
!
!    Input, real ( kind = 8 ) R, the divisor.
!
!    Output, complex ( kind = 8 ) C8_DIV_R8, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_div_r8
  real ( kind = 8 ) r
  complex ( kind = 8 ) z1

  c8_div_r8 = z1 / r

  return
end
function c8_exp ( z )

!*****************************************************************************80
!
!! C8_EXP evaluates the exponential of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the exponential of a C8 with the EXP function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double complex Z, the argument.
!
!    Output, double complex C8_EXP, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_exp
  complex ( kind = 8 ) z
  real ( kind = 8 ) zi
  real ( kind = 8 ) zr

  zr = real ( z, kind = 8 )
  zi = aimag ( z )

  c8_exp = exp ( zr ) * cmplx ( cos ( zi ), sin ( zi ), kind = 8 )

  return
end
function c8_i ( )

!*****************************************************************************80
!
!! C8_I returns the imaginary unit, i as a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, complex ( kind = 8 ) C8_I, the value of complex i.
!
  implicit none

  complex ( kind = 8 ) c8_i

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )

  return
end
function c8_imag ( z )

!*****************************************************************************80
!
!! C8_IMAG evaluates the imaginary part of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the imaginary part of a C8 with the AIMAG function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, real ( kind = 8 ) C8_IMAG, the function value.
!
  implicit none

  real ( kind = 8 ) c8_imag
  complex ( kind = 8 ) z

  c8_imag = aimag ( z )

  return
end
function c8_inv ( z )

!*****************************************************************************80
!
!! C8_INV evaluates the inverse of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the inverse of a C8 with the "1/" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_INV, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_inv
  complex ( kind = 8 ) z
  real ( kind = 8 ) z_imag
  real ( kind = 8 ) z_norm
  real ( kind = 8 ) z_real

  z_real = real ( z, kind = 8 )
  z_imag = aimag ( z )

  z_norm = sqrt ( z_real * z_real + z_imag * z_imag )

  c8_inv = cmplx ( z_real, - z_imag ) / z_norm / z_norm

  return
end
function c8_le_l1 ( x, y )

!*****************************************************************************80
!
!! C8_LE_L1 := X <= Y for C8 values, and the L1 norm.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    The L1 norm can be defined here as:
!
!      C8_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, Y, the values to be compared.
!
!    Output, logical C8_LE_L1, is TRUE if X <= Y.
!
  implicit none

  logical c8_le_l1
  complex ( kind = 8 ) x
  complex ( kind = 8 ) y

  if ( abs ( real ( x, kind = 8 ) ) + abs ( aimag ( x ) ) <= &
       abs ( real ( y, kind = 8 ) ) + abs ( aimag ( y ) ) ) then
    c8_le_l1 = .true.
  else
    c8_le_l1 = .false.
  end if

  return
end
function c8_le_l2 ( x, y )

!*****************************************************************************80
!
!! C8_LE_L2 := X <= Y for C8 values, and the L2 norm.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    The L2 norm can be defined here as:
!
!      C8_NORM_L2(X) = sqrt ( ( real (X) )**2 + ( imag (X) )**2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, Y, the values to be compared.
!
!    Output, logical C8_LE_L2, is TRUE if X <= Y.
!
  implicit none

  logical c8_le_l2
  complex ( kind = 8 ) x
  complex ( kind = 8 ) y

  if ( ( real ( x, kind = 8 ) )**2 + ( aimag ( x ) )**2 <= &
       ( real ( y, kind = 8 ) )**2 + ( aimag ( y ) )**2 ) then
    c8_le_l2 = .true.
  else
    c8_le_l2 = .false.
  end if

  return
end
function c8_le_li ( x, y )

!*****************************************************************************80
!
!! C8_LE_LI := X <= Y for C8 values, and the L Infinity norm.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    The L Infinity norm can be defined here as:
!
!      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, Y, the values to be compared.
!
!    Output, logical C8_LE_LI, is TRUE if X <= Y.
!
  implicit none

  logical c8_le_li
  complex ( kind = 8 ) x
  complex ( kind = 8 ) y

  if ( max ( abs ( real ( x, kind = 8 ) ), abs ( aimag ( x ) ) ) <= &
       max ( abs ( real ( y, kind = 8 ) ), abs ( aimag ( y ) ) ) ) then
    c8_le_li = .true.
  else
    c8_le_li = .false.
  end if

  return
end
function c8_log ( z )

!*****************************************************************************80
!
!! C8_LOG evaluates the logarithm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the logarithm of a C8 with the LOG function.
!
!    Here we use the relationship:
!
!      C8_LOG ( Z ) = LOG ( MAG ( Z ) ) + i * ARG ( Z )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_LOG, the function value.
!
  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) c8_arg
  complex ( kind = 8 ) c8_i
  complex ( kind = 8 ) c8_log
  real ( kind = 8 ) c8_mag
  real ( kind = 8 ) mag
  complex ( kind = 8 ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00 )

  arg = c8_arg ( z )
  mag = c8_mag ( z )

  c8_log = log ( mag ) + c8_i * arg
 
  return
end
function c8_mag ( x )

!*****************************************************************************80
!
!! C8_MAG returns the magnitude of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the magnitude of a C8 with the ABS function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) C8_MAG, the function value.
!
  implicit none

  real ( kind = 8 ) c8_mag
  complex ( kind = 8 ) x

  c8_mag = sqrt ( ( real ( x, kind = 8 ) )**2 + ( aimag ( x ) )**2 )

  return
end
function c8_mul ( z1, z2 )

!*****************************************************************************80
!
!! C8_MUL multiplies two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports multiplication of C8's with the "*" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z1, Z2, the values to multiply.
!
!    Output, complex ( kind = 8 ) C8_MUL, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_mul
  complex ( kind = 8 ) z1
  complex ( kind = 8 ) z2

  c8_mul = z1 * z2

  return
end
function c8_neg ( c1 )

!*****************************************************************************80
!
!! C8_NEG returns the negative of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports negation of a C8 with the "-" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) C1, the value to be negated.
!
!    Output, complex ( kind = 8 ) C8_NEG, the function value.
!
  implicit none

  complex ( kind = 8 ) c1
  complex ( kind = 8 ) c8_neg

  c8_neg = - c1

  return
end
function c8_norm_l1 ( x )

!*****************************************************************************80
!
!! C8_NORM_L1 evaluates the L1 norm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    Numbers of equal norm lie along diamonds centered at (0,0).
!
!    The L1 norm can be defined here as:
!
!      C8_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the value whose norm is desired.
!
!    Output, real ( kind = 8 ) C8_NORM_L1, the norm of X.
!
  implicit none

  real ( kind = 8 ) c8_norm_l1
  complex ( kind = 8 ) x

  c8_norm_l1 = abs ( real ( x, kind = 8 ) ) + abs ( aimag ( x ) )

  return
end
function c8_norm_l2 ( x )

!*****************************************************************************80
!
!! C8_NORM_L2 evaluates the L2 norm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    Numbers of equal norm lie on circles centered at (0,0).
!
!    The L2 norm can be defined here as:
!
!      C8_NORM_L2(X) = sqrt ( ( real (X) )**2 + ( imag ( X ) )**2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the value whose norm is desired.
!
!    Output, real ( kind = 8 ) C8_NORM_L2, the 2-norm of X.
!
  implicit none

  real ( kind = 8 ) c8_norm_l2
  complex ( kind = 8 ) x

  c8_norm_l2 = sqrt ( ( real ( x, kind = 8 ) )**2 &
                   + ( aimag ( x ) )**2 )

  return
end
function c8_norm_li ( x )

!*****************************************************************************80
!
!! C8_NORM_LI evaluates the L-infinity norm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    Numbers of equal norm lie along squares whose centers are at (0,0).
!
!    The L-infinity norm can be defined here as:
!
!      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the value whose norm is desired.
!
!    Output, real ( kind = 8 ) C8_NORM_LI, the infinity norm of X.
!
  implicit none

  real ( kind = 8 ) c8_norm_li
  complex ( kind = 8 ) x

  c8_norm_li = max ( abs ( real ( x, kind = 8 ) ), abs ( aimag ( x ) ) )

  return
end
function c8_normal_01 ( seed )

!*****************************************************************************80
!
!! C8_NORMAL_01 returns a unit pseudonormal C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, complex ( kind = 8 ) C8_NORMAL_01, a unit pseudornormal value.
!
  implicit none

  complex ( kind = 8 ) c8_normal_01
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) x_c
  real ( kind = 8 ) x_r

  v1 = r8_uniform_01 ( seed )
  v2 = r8_uniform_01 ( seed )

  x_r = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * r8_pi * v2 )
  x_c = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * r8_pi * v2 )

  c8_normal_01 = cmplx ( x_r, x_c, kind = 8 )

  return
end
function c8_one ( )

!*****************************************************************************80
!
!! C8_ONE returns the value of 1 as a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, complex ( kind = 8 ) C8_ONE, the value of complex 1.
!
  implicit none

  complex ( kind = 8 ) c8_one

  c8_one = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )

  return
end
subroutine c8_print ( a, title )

!*****************************************************************************80
!
!! C8_PRINT prints a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) A, the value to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  complex ( kind = 8 ) a
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a,2x,a,g14.6,a,g14.6,a)' ) &
      trim ( title ), '(', real ( a ), ',', imag ( a ), ')'
  else
    write ( *, '(a,g14.6,a,g14.6,a)' ) &
      '(', real ( a ), ',', imag ( a ), ')'
  end if

  return
end
function c8_real ( z )

!*****************************************************************************80
!
!! C8_REAL evaluates the real part of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the real part of a C8 with the REAL function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, real ( kind = 8 ) C8_REAL, the function value.
!
  implicit none

  real ( kind = 8 ) c8_real
  complex ( kind = 8 ) z

  c8_real = real ( z )

  return
end
function c8_sin ( z )

!*****************************************************************************80
!
!! C8_SIN evaluates the sine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the sine of a C8 with the SIN function.
!
!    We use the relationship:
!
!      C8_SIN ( C ) = - i * ( C8_EXP ( i * C ) - C8_EXP ( - i * C ) ) / 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_SIN, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_exp
  complex ( kind = 8 ) c8_i
  complex ( kind = 8 ) c8_sin
  complex ( kind = 8 ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00 )

  c8_sin = - c8_i * ( c8_exp ( c8_i * z ) - c8_exp ( - c8_i * z ) ) / 2.0D+00

  return
end
function c8_sinh ( z )

!*****************************************************************************80
!
!! C8_SINH evaluates the hyperbolic sine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the hyperbolic sine of a C8.
!
!    We use the relationship:
!
!      C8_SINH ( C ) = ( C8_EXP ( C ) - C8_EXP ( - C ) ) / 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_SINH, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_exp
  complex ( kind = 8 ) c8_sinh
  complex ( kind = 8 ) z

  c8_sinh = ( c8_exp ( z ) - c8_exp ( - z ) ) / 2.0D+00

  return
end
function c8_sqrt ( x )

!*****************************************************************************80
!
!! C8_SQRT returns the principal square root of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 supports the square root of a C8 with the SQRT function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the argument.
!
!    Output, complex ( kind = 8 ) C8_SQRT, the function value.
!
  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) c8_arg
  real ( kind = 8 ) c8_mag
  complex ( kind = 8 ) c8_sqrt
  real ( kind = 8 ) mag
  complex ( kind = 8 ) x

  arg = c8_arg ( x )
  mag = c8_mag ( x )

  if ( mag == 0.0D+00 ) then

    c8_sqrt = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  else

    c8_sqrt = sqrt ( mag ) &
      * cmplx ( cos ( arg / 2.0D+00 ), &
                sin ( arg / 2.0D+00 ), kind = 8 )

  end if

  return
end
function c8_sub ( z1, z2 )

!*****************************************************************************80
!
!! C8_SUB subtracts two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 directly supports C8 subtraction with the "-" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z1, Z2, the values to subtract.
!
!    Output, complex ( kind = 8 ) C8_SUB, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_sub
  complex ( kind = 8 ) z1
  complex ( kind = 8 ) z2

  c8_sub = z1 - z2

  return
end
subroutine c8_swap ( x, y )

!*****************************************************************************80
!
!! C8_SWAP swaps two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, complex ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  complex ( kind = 8 ) x
  complex ( kind = 8 ) y
  complex ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
function c8_tan ( z )

!*****************************************************************************80
!
!! C8_TAN evaluates the tangent of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the tangent of a C8.
!
!    We use the relationship:
!
!      C8_TAN ( C ) = - i * ( C8_EXP ( i * C ) - C8_EXP ( - i * C ) ) 
!                         / ( C8_EXP ( I * C ) + C8_EXP ( - i * C ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_TAN, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_exp
  complex ( kind = 8 ) c8_i
  complex ( kind = 8 ) c8_tan
  complex ( kind = 8 ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00 )

  c8_tan =  - c8_i * ( c8_exp ( c8_i * z ) - c8_exp ( - c8_i * z ) ) &
         /           ( c8_exp ( c8_i * z ) + c8_exp ( - c8_i * z ) )

  return
end
function c8_tanh ( z )

!*****************************************************************************80
!
!! C8_TANH evaluates the hyperbolic tangent of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    FORTRAN90 does not support the hyperbolic tangent of a C8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) C8_TANH, the function value.
!
  implicit none

  complex ( kind = 8 ) c8_exp
  complex ( kind = 8 ) c8_tanh
  complex ( kind = 8 ) z

  c8_tanh = ( c8_exp ( z ) - c8_exp ( - z ) ) &
          / ( c8_exp ( z ) + c8_exp ( - z ) )

  return
end
subroutine c8_to_cartesian ( z, x, y )

!*****************************************************************************80
!
!! C8_TO_CARTESIAN converts a C8 to Cartesian form.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, real ( kind = 8 ) X, Y, the Cartesian form.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  complex ( kind = 8 ) z

  x = real ( z )
  y = aimag ( z )

  return
end
subroutine c8_to_polar ( z, r, theta )

!*****************************************************************************80
!
!! C8_TO_POLAR converts a C8 to polar form.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, real ( kind = 8 ) R, THETA, the polar form.
!
  implicit none

  real ( kind = 8 ) c8_arg
  real ( kind = 8 ) c8_mag
  real ( kind = 8 ) r
  real ( kind = 8 ) theta
  complex ( kind = 8 ) z

  r = c8_mag ( z )
  theta = c8_arg ( z )

  return
end
function c8_uniform_01 ( seed )

!*****************************************************************************80
!
!! C8_UNIFORM_01 returns a unit pseudorandom C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C8_UNIFORM_01, a pseudorandom complex value.
!
  implicit none

  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  theta = 2.0D+00 * r8_pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

  c8_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

  return
end
function c8_zero ( )

!*****************************************************************************80
!
!! C8_ZERO returns the value of 0 as a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, complex ( kind = 8 ) C8_ZERO, the value of complex 0.
!
  implicit none

  complex ( kind = 8 ) c8_zero

  c8_zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  return
end
subroutine c8mat_identity ( n, a )

!*****************************************************************************80
!
!! C8MAT_IDENTITY sets a C8MAT to the identity.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, complex ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i

  a(1:n,1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  do i = 1, n
    a(i,i) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
  end do

  return
end
subroutine c8mat_indicator ( m, n, a )

!*****************************************************************************80
!
!! C8MAT_INDICATOR returns the C8MAT indicator matrix.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2008
!
!  Author
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Output, complex ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n
    do i = 1, m
      a(i,j) = cmplx ( i, j, kind = 8 )
    end do
  end do

  return
end
subroutine c8mat_nint ( m, n, a )

!*****************************************************************************80
!
!! C8MAT_NINT rounds the entries of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input/output, complex ( kind = 8 ) A(M,N), the matrix to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n
    do i = 1, m
      a(i,j) = cmplx ( nint ( real ( a(i,j) ) ), &
                       nint ( imag ( a(i,j) ) ), kind = 8 )
    end do
  end do

  return
end
subroutine c8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! C8MAT_PRINT prints a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call c8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8MAT_PRINT_SOME prints some of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 4
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)' 
    return
  end if
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == zero ) then
          ctemp(j2) = '       0.0          '
        else if ( imag ( a(i,j) ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      theta = 2.0D+00 * r8_pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

    end do

  end do

  return
end
subroutine c8vec_indicator ( n, a )

!*****************************************************************************80
!
!! C8VEC_INDICATOR sets a C8VEC to the indicator vector.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = cmplx ( i, -i, kind = 8 )
  end do

  return
end
subroutine c8vec_nint ( n, a )

!*****************************************************************************80
!
!! C8VEC_NINT rounds the entries of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, complex ( kind = 8 ) A(N), the vector to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = cmplx ( nint ( real ( a(i) ) ), &
                   nint ( imag ( a(i) ) ), kind = 8 )
  end do

  return
end
function c8vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! C8VEC_NORM_L2 returns the L2 norm of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The vector L2 norm is defined as:
!
!      C8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, complex ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) C8VEC_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  real ( kind = 8 ) c8vec_norm_l2

  c8vec_norm_l2 = sqrt ( sum ( conjg ( a(1:n) ) * a(1:n) ) )

  return
end
function c8vec_norm_squared ( n, a )

!*****************************************************************************80
!
!! C8VEC_NORM_SQUARED returns the square of the L2 norm of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The square of the vector L2 norm is defined as:
!
!      C8VEC_NORM_SQUARED = sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, complex ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) C8VEC_NORM_SQUARED, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  real ( kind = 8 ) c8vec_norm_squared

  c8vec_norm_squared = sum ( conjg ( a(1:n) ) * a(1:n) )

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! C8VEC_PRINT prints a C8VEC, with an optional title.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, complex ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, a(i)
  end do

  return
end
subroutine c8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! C8VEC_PRINT_PART prints "part" of a C8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) max_print
  character ( len = * )  title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), &
      '...more entries...'

  end if

  return
end
subroutine c8vec_print_some ( n, x, i_lo, i_hi, title )

!*****************************************************************************80
!
!! C8VEC_PRINT_SOME prints some of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last entries
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i_hi
  integer   ( kind = 4 ) i_lo
  character ( len = * )  title
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = max ( 1, i_lo ), min ( n, i_hi )
    write ( *, '(2x,i8,2x,2g14.6)' ) i, x(i)
  end do

  return
end
subroutine c8vec_sort_a_l1 ( n, x )

!*****************************************************************************80
!
!! C8VEC_SORT_A_L1 ascending sorts a C8VEC by L1 norm.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The L1 norm of A+Bi is abs(A) + abs(B).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, complex ( kind = 8 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical              c8_le_l1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c8_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c8_le_l1 ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_sort_a_l2 ( n, x )

!*****************************************************************************80
!
!! C8VEC_SORT_A_L2 ascending sorts a C8VEC by L2 norm.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The L2 norm of A+Bi is sqrt ( A**2 + B**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, complex ( kind = 8 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical              c8_le_l2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c8_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c8_le_l2 ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_sort_a_li ( n, x )

!*****************************************************************************80
!
!! C8VEC_SORT_A_LI ascending sorts a C8VEC by L-infinity norm.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The L infinity norm of A+Bi is max ( abs ( A ), abs ( B ) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, complex ( kind = 8 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical              c8_le_li
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c8_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c8_le_li ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_spiral ( n, m, c1, c2, c )

!*****************************************************************************80
!
!! C8VEC_SPIRAL returns N points on a spiral between C1 and C2.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    Let the polar form of C1 be ( R1, T1 ) and the polar form of C2 
!    be ( R2, T2 ) where, if necessary, we increase T2 by 2*PI so that T1 <= T2.
!    
!    Then the polar form of the I-th point C(I) is:
!
!      R(I) = ( ( N - I     ) * R1 
!             + (     I - 1 ) * R2 ) 
!              / ( N    - 1 )
!
!      T(I) = ( ( N - I     ) * T1 
!             + (     I - 1 ) * ( T2 + M * 2 * PI ) ) 
!             / ( N     - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points on the spiral.
!
!    Input, integer ( kind = 4 ) M, the number of full circuits the spiral makes.
!
!    Input, complex ( kind = 8 ) C1, C2, the first and last points on the spiral.
!
!    Output, complex ( kind = 8 ) C(N), the points.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  complex ( kind = 8 ) c1
  complex ( kind = 8 ) c2
  real ( kind = 8 ) c8_arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) ri
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) ti

  r1 = abs ( c1 )
  r2 = abs ( c2 )

  t1 = c8_arg ( c1 )
  t2 = c8_arg ( c2 )

  if ( m == 0 ) then

    if ( t2 < t1 ) then
      t2 = t2 + 2.0D+00 * r8_pi
    end if

  else if ( 0 < m ) then

    if ( t2 < t1 ) then
      t2 = t2 + 2.0D+00 * r8_pi
    end if

    t2 = t2 + real ( m, kind = 8 ) * 2.0D+00 * r8_pi

  else if ( m < 0 ) then

    if ( t1 < t2 ) then
      t2 = t2 - 2.0D+00 * r8_pi
    end if

    t2 = t2 - real ( m, kind = 8 ) * 2.0D+00 * r8_pi

  end if

  do i = 1, n

    ri = ( real ( n - i,     kind = 8 ) * r1 &
         + real (     i - 1, kind = 8 ) * r2 ) &
         / real ( n     - 1, kind = 8 )

    ti = ( real ( n - i,     kind = 8 ) * t1 &
         + real (     i - 1, kind = 8 ) * t2 ) &
         / real ( n     - 1, kind = 8 )

    call polar_to_c8 ( ri, ti, c(i) )

  end do

  return
end
subroutine c8vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
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
!    Input, integer ( kind = 4 ) N, the number of values to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    theta = 2.0D+00 * r8_pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

  end do

  return
end
subroutine c8vec_unity ( n, a )

!*****************************************************************************80
!
!! C8VEC_UNITY returns the N roots of unity.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
!
!    X(1:N)**N = ( (1,0), (1,0), ..., (1,0) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 8 ) A(N), the N roots of unity.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) theta

  do i = 1, n
    theta = r8_pi * real ( 2 * ( i - 1 ), kind = 8 ) / real ( n, kind = 8 )
    a(i) = cmplx ( cos ( theta ), sin ( theta ), kind = 8 )
  end do

  return
end
subroutine cartesian_to_c8 ( x, y, z )

!*****************************************************************************80
!
!! CARTESIAN_TO_C8 converts a Cartesian form to a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the Cartesian form.
!
!    Output, complex ( kind = 8 ) Z, the complex number.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  complex ( kind = 8 ) z

  z = cmplx ( x, y )

  return
end
subroutine polar_to_c8 ( r, theta, z )

!*****************************************************************************80
!
!! POLAR_TO_C8 converts a polar form to a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, THETA, the polar form.
!
!    Output, complex ( kind = 8 ) Z, the complex number.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) theta
  complex ( kind = 8 ) z

  z = r * cmplx ( cos ( theta ), sin ( theta ) )

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
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
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements 
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
