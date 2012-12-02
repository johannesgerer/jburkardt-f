function c4_cos ( z )

!*****************************************************************************80
!
!! C4_COS evaluates the cosine of a C4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_COS, the cosine of Z.
!
  implicit none

  complex ( kind = 4 ) c4_cos
  real ( kind = 4 ) cs
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  complex ( kind = 4 ) z

  x = real ( z, kind = 4 )
  y = imag ( z )

  cs = cos ( x )

  c4_cos = cmplx ( cs * cosh ( y ), - sin ( x ) * sinh ( y ) )

  return
end
function c4_sin ( z )

!*****************************************************************************80
!
!! C4_SIN evaluates the sine of a C4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_SIN, the sine of Z.
!
  implicit none

  complex ( kind = 4 ) c4_sin
  real ( kind = 4 ) sn
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  complex ( kind = 4 ) z

  x = real ( z, kind = 4 )
  y = imag ( z )

  sn = sin ( x )

  c4_sin = cmplx ( sn * cosh ( y ), cos ( x ) * sinh ( y ) )

  return
end
function i4_mach ( i )

!*****************************************************************************80
!
!! I4_MACH returns integer machine dependent constants.
!
!  Discussion:
!
!    Input/output unit numbers.
!
!      I4_MACH(1) = the standard input unit.
!      I4_MACH(2) = the standard output unit.
!      I4_MACH(3) = the standard punch unit.
!      I4_MACH(4) = the standard error message unit.
!
!    Words.
!
!      I4_MACH(5) = the number of bits per integer storage unit.
!      I4_MACH(6) = the number of characters per integer storage unit.
!
!    Integers.
!
!    Assume integers are represented in the S digit base A form:
!
!      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
!
!    where 0 <= X(1:S-1) < A.
!
!      I4_MACH(7) = A, the base.
!      I4_MACH(8) = S, the number of base A digits.
!      I4_MACH(9) = A^S-1, the largest integer.
!
!    Floating point numbers
!
!    Assume floating point numbers are represented in the T digit 
!    base B form:
!
!      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
!
!    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
!
!      I4_MACH(10) = B, the base.
!
!    Single precision
!
!      I4_MACH(11) = T, the number of base B digits.
!      I4_MACH(12) = EMIN, the smallest exponent E.
!      I4_MACH(13) = EMAX, the largest exponent E.
!
!    Double precision
!
!      I4_MACH(14) = T, the number of base B digits.
!      I4_MACH(15) = EMIN, the smallest exponent E.
!      I4_MACH(16) = EMAX, the largest exponent E.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528,
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 16.
!
!    Output, integer ( kind = 4 ) I4_MACH, the value of the chosen parameter.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_mach

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
    write ( *, '(a,i12)' ) '  I = ', i
    i4_mach = 0
    stop
  else if ( i == 1 ) then
    i4_mach = 5
  else if ( i == 2 ) then
    i4_mach = 6
  else if ( i == 3 ) then
    i4_mach = 7
  else if ( i == 4 ) then
    i4_mach = 6
  else if ( i == 5 ) then
    i4_mach = 32
  else if ( i == 6 ) then
    i4_mach = 4
  else if ( i == 7 ) then
    i4_mach = 2
  else if ( i == 8 ) then
    i4_mach = 31
  else if ( i == 9 ) then
    i4_mach = 2147483647
  else if ( i == 10 ) then
    i4_mach = 2
  else if ( i == 11 ) then
    i4_mach = 24
  else if ( i == 12 ) then
    i4_mach = -125
  else if ( i == 13 ) then
    i4_mach = 128
  else if ( i == 14 ) then
    i4_mach = 53
  else if ( i == 15 ) then
    i4_mach = -1021
  else if ( i == 16 ) then
    i4_mach = 1024
  else if ( 16 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
    write ( *, '(a,i12)' ) '  I = ', i
    i4_mach = 0
    stop
  end if

  return
end
function r4_acos ( x )

!*****************************************************************************80
!
!! R4_ACOS evaluates the arc-cosine of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_ACOS, the arc-cosine of X.
!
  implicit none

  real ( kind = 4 ), parameter :: pi2 = 1.57079632679489661923E+00
  real ( kind = 4 ) r4_acos
  real ( kind = 4 ) r4_asin
  real ( kind = 4 ) x

  r4_acos = pi2 - r4_asin ( x )

  return
end
function r4_acosh ( x )

!*****************************************************************************80
!
!! R4_ACOSH evaluates the arc-hyperbolic cosine of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_ACOSH, the arc-hyperbolic cosine of X.
!
  implicit none

  real ( kind = 4 ), save :: aln2 = 0.69314718055994530942E+00
  real ( kind = 4 ) r4_acosh
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ), save :: xmax = 0.0E+00

  if ( xmax == 0.0E+00 ) then
    xmax = 1.0E+00 / sqrt ( r4_mach ( 3 ) )
  end if

  if ( x < 1.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_ACOSH - Fatal error!'
    write ( *, '(a)' ) '  X < 1.' 
    stop
  end if

  if ( x < xmax ) then
    value = log ( x + sqrt ( x * x - 1.0E+00 ) )
  else
    value = aln2 + log ( x )
  end if

  r4_acosh = value

  return
end
subroutine r4_admp ( x, ampl, phi )

!*****************************************************************************80
!
!! R4_ADMP: modulus and phase of the derivative of the Airy function.
!
!  Description:
!
!    This function must only be called when X <= -1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) AMPL, PHI, the modulus and phase of the 
!    derivative of the Airy function.
!
  implicit none

  real ( kind = 4 ) ampl
  real ( kind = 4 ) an20cs(16)
  real ( kind = 4 ) an21cs(24)
  real ( kind = 4 ) an22cs(33)
  real ( kind = 4 ) aph0cs(15)
  real ( kind = 4 ) aph1cs(22)
  real ( kind = 4 ) aph2cs(32)
  real ( kind = 4 ) eta
  integer ( kind = 4 ) nan20
  integer ( kind = 4 ) nan21
  integer ( kind = 4 ) nan22
  integer ( kind = 4 ) naph0
  integer ( kind = 4 ) naph1
  integer ( kind = 4 ) naph2
  real ( kind = 4 ) phi
  real ( kind = 4 ) pi34
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sqrtx
  real ( kind = 4 ) x
  real ( kind = 4 ) xsml
  real ( kind = 4 ) z

  save an20cs
  save an21cs
  save an22cs
  save aph0cs
  save aph1cs
  save aph2cs
  save nan20
  save nan21
  save nan22
  save naph0
  save naph1
  save naph2
  save pi34
  save xsml

  data an22cs(  1) /     0.0537418629629794329E+00/
  data an22cs(  2) /    -0.0126661435859883193E+00/
  data an22cs(  3) /    -0.0011924334106593007E+00/
  data an22cs(  4) /    -0.0002032327627275655E+00/
  data an22cs(  5) /    -0.0000446468963075164E+00/
  data an22cs(  6) /    -0.0000113359036053123E+00/
  data an22cs(  7) /    -0.0000031641352378546E+00/
  data an22cs(  8) /    -0.0000009446708886149E+00/
  data an22cs(  9) /    -0.0000002966562236472E+00/
  data an22cs( 10) /    -0.0000000969118892024E+00/
  data an22cs( 11) /    -0.0000000326822538653E+00/
  data an22cs( 12) /    -0.0000000113144618964E+00/
  data an22cs( 13) /    -0.0000000040042691002E+00/
  data an22cs( 14) /    -0.0000000014440333684E+00/
  data an22cs( 15) /    -0.0000000005292853746E+00/
  data an22cs( 16) /    -0.0000000001967763374E+00/
  data an22cs( 17) /    -0.0000000000740800096E+00/
  data an22cs( 18) /    -0.0000000000282016314E+00/
  data an22cs( 19) /    -0.0000000000108440066E+00/
  data an22cs( 20) /    -0.0000000000042074801E+00/
  data an22cs( 21) /    -0.0000000000016459150E+00/
  data an22cs( 22) /    -0.0000000000006486827E+00/
  data an22cs( 23) /    -0.0000000000002574095E+00/
  data an22cs( 24) /    -0.0000000000001027889E+00/
  data an22cs( 25) /    -0.0000000000000412846E+00/
  data an22cs( 26) /    -0.0000000000000166711E+00/
  data an22cs( 27) /    -0.0000000000000067657E+00/
  data an22cs( 28) /    -0.0000000000000027585E+00/
  data an22cs( 29) /    -0.0000000000000011296E+00/
  data an22cs( 30) /    -0.0000000000000004645E+00/
  data an22cs( 31) /    -0.0000000000000001917E+00/
  data an22cs( 32) /    -0.0000000000000000794E+00/
  data an22cs( 33) /    -0.0000000000000000330E+00/

  data an21cs(  1) /     0.0198313155263169394E+00/
  data an21cs(  2) /    -0.0029376249067087533E+00/
  data an21cs(  3) /    -0.0001136260695958196E+00/
  data an21cs(  4) /    -0.0000100554451087156E+00/
  data an21cs(  5) /    -0.0000013048787116563E+00/
  data an21cs(  6) /    -0.0000002123881993151E+00/
  data an21cs(  7) /    -0.0000000402270833384E+00/
  data an21cs(  8) /    -0.0000000084996745953E+00/
  data an21cs(  9) /    -0.0000000019514839426E+00/
  data an21cs( 10) /    -0.0000000004783865344E+00/
  data an21cs( 11) /    -0.0000000001236733992E+00/
  data an21cs( 12) /    -0.0000000000334137486E+00/
  data an21cs( 13) /    -0.0000000000093702824E+00/
  data an21cs( 14) /    -0.0000000000027130128E+00/
  data an21cs( 15) /    -0.0000000000008075954E+00/
  data an21cs( 16) /    -0.0000000000002463214E+00/
  data an21cs( 17) /    -0.0000000000000767656E+00/
  data an21cs( 18) /    -0.0000000000000243883E+00/
  data an21cs( 19) /    -0.0000000000000078831E+00/
  data an21cs( 20) /    -0.0000000000000025882E+00/
  data an21cs( 21) /    -0.0000000000000008619E+00/
  data an21cs( 22) /    -0.0000000000000002908E+00/
  data an21cs( 23) /    -0.0000000000000000993E+00/
  data an21cs( 24) /    -0.0000000000000000343E+00/

  data an20cs(  1) /     0.0126732217145738027E+00/
  data an20cs(  2) /    -0.0005212847072615621E+00/
  data an20cs(  3) /    -0.0000052672111140370E+00/
  data an20cs(  4) /    -0.0000001628202185026E+00/
  data an20cs(  5) /    -0.0000000090991442687E+00/
  data an20cs(  6) /    -0.0000000007438647126E+00/
  data an20cs(  7) /    -0.0000000000795494752E+00/
  data an20cs(  8) /    -0.0000000000104050944E+00/
  data an20cs(  9) /    -0.0000000000015932426E+00/
  data an20cs( 10) /    -0.0000000000002770648E+00/
  data an20cs( 11) /    -0.0000000000000535343E+00/
  data an20cs( 12) /    -0.0000000000000113062E+00/
  data an20cs( 13) /    -0.0000000000000025772E+00/
  data an20cs( 14) /    -0.0000000000000006278E+00/
  data an20cs( 15) /    -0.0000000000000001621E+00/
  data an20cs( 16) /    -0.0000000000000000441E+00/

  data aph2cs(  1) /    -0.2057088719781465107E+00/
  data aph2cs(  2) /     0.0422196961357771922E+00/
  data aph2cs(  3) /     0.0020482560511207275E+00/
  data aph2cs(  4) /     0.0002607800735165006E+00/
  data aph2cs(  5) /     0.0000474824268004729E+00/
  data aph2cs(  6) /     0.0000105102756431612E+00/
  data aph2cs(  7) /     0.0000026353534014668E+00/
  data aph2cs(  8) /     0.0000007208824863499E+00/
  data aph2cs(  9) /     0.0000002103236664473E+00/
  data aph2cs( 10) /     0.0000000644975634555E+00/
  data aph2cs( 11) /     0.0000000205802377264E+00/
  data aph2cs( 12) /     0.0000000067836273921E+00/
  data aph2cs( 13) /     0.0000000022974015284E+00/
  data aph2cs( 14) /     0.0000000007961306765E+00/
  data aph2cs( 15) /     0.0000000002813860610E+00/
  data aph2cs( 16) /     0.0000000001011749057E+00/
  data aph2cs( 17) /     0.0000000000369306738E+00/
  data aph2cs( 18) /     0.0000000000136615066E+00/
  data aph2cs( 19) /     0.0000000000051142751E+00/
  data aph2cs( 20) /     0.0000000000019351689E+00/
  data aph2cs( 21) /     0.0000000000007393607E+00/
  data aph2cs( 22) /     0.0000000000002849792E+00/
  data aph2cs( 23) /     0.0000000000001107281E+00/
  data aph2cs( 24) /     0.0000000000000433412E+00/
  data aph2cs( 25) /     0.0000000000000170801E+00/
  data aph2cs( 26) /     0.0000000000000067733E+00/
  data aph2cs( 27) /     0.0000000000000027017E+00/
  data aph2cs( 28) /     0.0000000000000010835E+00/
  data aph2cs( 29) /     0.0000000000000004367E+00/
  data aph2cs( 30) /     0.0000000000000001769E+00/
  data aph2cs( 31) /     0.0000000000000000719E+00/
  data aph2cs( 32) /     0.0000000000000000294E+00/

  data aph1cs(  1) /    -0.1024172908077571694E+00/
  data aph1cs(  2) /     0.0071697275146591248E+00/
  data aph1cs(  3) /     0.0001209959363122329E+00/
  data aph1cs(  4) /     0.0000073361512841220E+00/
  data aph1cs(  5) /     0.0000007535382954272E+00/
  data aph1cs(  6) /     0.0000001041478171741E+00/
  data aph1cs(  7) /     0.0000000174358728519E+00/
  data aph1cs(  8) /     0.0000000033399795033E+00/
  data aph1cs(  9) /     0.0000000007073075174E+00/
  data aph1cs( 10) /     0.0000000001619187515E+00/
  data aph1cs( 11) /     0.0000000000394539982E+00/
  data aph1cs( 12) /     0.0000000000101192282E+00/
  data aph1cs( 13) /     0.0000000000027092778E+00/
  data aph1cs( 14) /     0.0000000000007523806E+00/
  data aph1cs( 15) /     0.0000000000002156369E+00/
  data aph1cs( 16) /     0.0000000000000635283E+00/
  data aph1cs( 17) /     0.0000000000000191757E+00/
  data aph1cs( 18) /     0.0000000000000059143E+00/
  data aph1cs( 19) /     0.0000000000000018597E+00/
  data aph1cs( 20) /     0.0000000000000005950E+00/
  data aph1cs( 21) /     0.0000000000000001934E+00/
  data aph1cs( 22) /     0.0000000000000000638E+00/

  data aph0cs(  1) /    -0.0855849241130933257E+00/
  data aph0cs(  2) /     0.0011214378867065261E+00/
  data aph0cs(  3) /     0.0000042721029353664E+00/
  data aph0cs(  4) /     0.0000000817607381483E+00/
  data aph0cs(  5) /     0.0000000033907645000E+00/
  data aph0cs(  6) /     0.0000000002253264423E+00/
  data aph0cs(  7) /     0.0000000000206284209E+00/
  data aph0cs(  8) /     0.0000000000023858763E+00/
  data aph0cs(  9) /     0.0000000000003301618E+00/
  data aph0cs( 10) /     0.0000000000000527010E+00/
  data aph0cs( 11) /     0.0000000000000094555E+00/
  data aph0cs( 12) /     0.0000000000000018709E+00/
  data aph0cs( 13) /     0.0000000000000004024E+00/
  data aph0cs( 14) /     0.0000000000000000930E+00/
  data aph0cs( 15) /     0.0000000000000000229E+00/

  data nan20 / 0 /
  data nan21 / 0 /
  data nan22 / 0 /
  data naph0 / 0 /
  data naph1 / 0 /
  data naph2 / 0 /
  data pi34 / 2.3561944901923449E+00 /
  data xsml / 0.0E+00 /

  if ( nan20 == 0 ) then
    eta = 0.1E+00 * r4_mach ( 3 )
    nan20 = r4_inits ( an20cs, 16, eta )
    nan21 = r4_inits ( an21cs, 24, eta )
    nan22 = r4_inits ( an22cs, 33, eta )
    naph0 = r4_inits ( aph0cs, 15, eta )
    naph1 = r4_inits ( aph1cs, 22, eta )
    naph2 = r4_inits ( aph2cs, 32, eta )
    xsml = - ( 128.0E+00 / r4_mach ( 3 ) ) ** 0.3333E+00
  end if

  if ( x <= xsml ) then
    z = 1.0E+00
    ampl = 0.3125E+00 + r4_csevl ( z, an20cs, nan20 )
    phi = - 0.625E+00 + r4_csevl ( z, aph0cs, naph0 )
  else if ( x < - 4.0E+00 ) then
    z = 128.0E+00 / x / x / x + 1.0E+00
    ampl = 0.3125E+00 + r4_csevl ( z, an20cs, nan20 )
    phi = - 0.625E+00 + r4_csevl ( z, aph0cs, naph0 )
  else if ( x < - 2.0E+00 ) then
    z = ( 128.0E+00 / x / x / x + 9.0E+00 ) / 7.0E+00
    ampl = 0.3125E+00 + r4_csevl ( z, an21cs, nan21 )
    phi = - 0.625E+00 + r4_csevl ( z, aph1cs, naph1 )
  else if ( x <= - 1.0E+00 ) then
    z = ( 16.0E+00 / x / x / x + 9.0E+00 ) / 7.0E+00
    ampl = 0.3125E+00 + r4_csevl ( z, an22cs, nan22 )
    phi = - 0.625E+00 + r4_csevl ( z, aph2cs, naph2 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_ADMP - Fatal error!'
    write ( *, '(a)' ) '  - 1.0 < x.'
    stop
  end if

  sqrtx = sqrt ( - x )
  ampl = sqrt ( ampl * sqrtx )
  phi = pi34 - x * sqrtx * phi

  return
end
function r4_ai ( x )

!*****************************************************************************80
!
!! R4_AI evaluates the Airy function Ai of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_AI, the Airy function Ai of X.
!
  implicit none

  real ( kind = 4 ) aifcs(9)
  real ( kind = 4 ) aigcs(8)
  integer ( kind = 4 ) naif
  integer ( kind = 4 ) naig
  real ( kind = 4 ) r4_ai
  real ( kind = 4 ) r4_aie
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) theta
  real ( kind = 4 ) x
  real ( kind = 4 ) x3sml
  real ( kind = 4 ) xm
  real ( kind = 4 ) xmax
  real ( kind = 4 ) z

  save aifcs
  save aigcs
  save naif
  save naig
  save x3sml
  save xmax

  data aifcs( 1) /   -0.03797135849666999750E+00 /
  data aifcs( 2) /    0.05919188853726363857E+00 /
  data aifcs( 3) /    0.00098629280577279975E+00 /
  data aifcs( 4) /    0.00000684884381907656E+00 /
  data aifcs( 5) /    0.00000002594202596219E+00 /
  data aifcs( 6) /    0.00000000006176612774E+00 /
  data aifcs( 7) /    0.00000000000010092454E+00 /
  data aifcs( 8) /    0.00000000000000012014E+00 /
  data aifcs( 9) /    0.00000000000000000010E+00 /

  data aigcs( 1) /    0.01815236558116127E+00 /
  data aigcs( 2) /    0.02157256316601076E+00 /
  data aigcs( 3) /    0.00025678356987483E+00 /
  data aigcs( 4) /    0.00000142652141197E+00 /
  data aigcs( 5) /    0.00000000457211492E+00 /
  data aigcs( 6) /    0.00000000000952517E+00 /
  data aigcs( 7) /    0.00000000000001392E+00 /
  data aigcs( 8) /    0.00000000000000001E+00 /

  data naif / 0 /
  data naig / 0 /
  data x3sml / 0.0E+00 / 
  data xmax / 0.0E+00 /

  if ( naif == 0 ) then
    naif = r4_inits ( aifcs, 9, 0.1E+00 * r4_mach ( 3 ) )
    naig = r4_inits ( aigcs, 8, 0.1E+00 * r4_mach ( 3 ) )
    x3sml = r4_mach ( 3 )**0.3334E+00
    xmax = ( - 1.5E+00 * log ( r4_mach ( 1 ) ) )**0.6667E+00
    xmax = xmax - xmax * log ( xmax ) &
      / ( 4.0E+00 * xmax * sqrt ( xmax ) + 1.0E+00 ) - 0.01E+00
  end if

  if ( x < - 1.0E+00 ) then
    call r4_aimp ( x, xm, theta )
    r4_ai = xm * cos ( theta )
  else if ( abs ( x ) <= x3sml ) then
    z = 0.0E+00
    r4_ai = 0.375E+00 + ( r4_csevl ( z, aifcs, naif ) &
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) )
  else if ( x <= 1.0E+00 ) then
    z = x * x * x
    r4_ai = 0.375E+00 + ( r4_csevl ( z, aifcs, naif ) &
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) )
  else if ( x <= xmax ) then
    r4_ai = r4_aie ( x ) * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  else
    r4_ai = 0.0E+00
  end if

  return
end
function r4_aid ( x )

!*****************************************************************************80
!
!! R4_AID evaluates the derivative of the Airy function Ai of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_AID, the derivative of the Airy function Ai.
!
  implicit none

  real ( kind = 4 ) aifcs(8)
  real ( kind = 4 ) aigcs(9)
  integer ( kind = 4 ) naif
  integer ( kind = 4 ) naig
  real ( kind = 4 ) phi
  real ( kind = 4 ) r4_aid
  real ( kind = 4 ) r4_aide
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) x2
  real ( kind = 4 ) x2sml
  real ( kind = 4 ) x3
  real ( kind = 4 ) x3sml
  real ( kind = 4 ) xn

  save aifcs
  save aigcs
  save naif
  save naig
  save x2sml
  save x3sml

  data aifcs(  1) /     0.10527461226531408809E+00 /
  data aifcs(  2) /     0.01183613628152997844E+00 /
  data aifcs(  3) /     0.00012328104173225664E+00 /
  data aifcs(  4) /     0.00000062261225638140E+00 /
  data aifcs(  5) /     0.00000000185298887844E+00 /
  data aifcs(  6) /     0.00000000000363328873E+00 /
  data aifcs(  7) /     0.00000000000000504622E+00 /
  data aifcs(  8) /     0.00000000000000000522E+00 /

  data aigcs(  1) /     0.021233878150918666852E+00 /
  data aigcs(  2) /     0.086315930335214406752E+00 /
  data aigcs(  3) /     0.001797594720383231358E+00 /
  data aigcs(  4) /     0.000014265499875550693E+00 /
  data aigcs(  5) /     0.000000059437995283683E+00 /
  data aigcs(  6) /     0.000000000152403366479E+00 /
  data aigcs(  7) /     0.000000000000264587660E+00 /
  data aigcs(  8) /     0.000000000000000331562E+00 /
  data aigcs(  9) /     0.000000000000000000314E+00 /

  data naif / 0 /
  data naig / 0 /
  data x2sml / 0.0E+00 /
  data x3sml / 0.0E+00 /

  if ( naif == 0 ) then
    naif = r4_inits ( aifcs, 8, 0.1E+00 * r4_mach ( 3 ) )
    naig = r4_inits ( aigcs, 9, 0.1E+00 * r4_mach ( 3 ) )
    x3sml = r4_mach ( 3 ) ** 0.3334E+00
    x2sml = sqrt ( r4_mach ( 3 ) )
  end if

  if ( x < - 1.0E+00 ) then
    call r4_admp ( x, xn, phi )
    r4_aid = xn * cos ( phi )
  else if ( abs ( x ) <= x2sml ) then
    x2 = 0.0E+00
    x3 = 0.0E+00
    r4_aid = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) &
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00
  else if ( abs ( x ) <= x3sml ) then
    x2 = x * x
    x3 = 0.0E+00
    r4_aid = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) &
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00
  else if ( x <= 1.0E+00 ) then
    x2 = x * x
    x3 = x * x * x
    r4_aid = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) &
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00
  else
    r4_aid = r4_aide ( x ) * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  end if

  return
end
function r4_aide ( x )

!*****************************************************************************80
!
!! R4_AIDE: exponentially scaled derivative, Airy function Ai of an R4 argument.
!
!  Discussion:
!
!    if X <= 0,
!      R4_AIDE ( X ) = R4_AID ( X )
!    else
!      R4_AIDE ( X ) = R4_AID ( X ) * exp ( 2/3 * X**(3/2) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_AIDE, the exponentially scaled derivative of 
!    the Airy function Ai of X.
!
  implicit none

  real ( kind = 4 ) aifcs(8)
  real ( kind = 4 ) aigcs(9)
  real ( kind = 4 ) aip1cs(25)
  real ( kind = 4 ) aip2cs(15)
  real ( kind = 4 ) eta
  integer ( kind = 4 ) naif
  integer ( kind = 4 ) naig
  integer ( kind = 4 ) naip1
  integer ( kind = 4 ) naip2
  real ( kind = 4 ) phi
  real ( kind = 4 ) r4_aide
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sqrtx
  real ( kind = 4 ) x
  real ( kind = 4 ) x2
  real ( kind = 4 ) x2sml
  real ( kind = 4 ) x3
  real ( kind = 4 ) x32sml
  real ( kind = 4 ) x3sml
  real ( kind = 4 ) xn
  real ( kind = 4 ) xbig
  real ( kind = 4 ) z

  save aifcs
  save aigcs
  save aip1cs
  save aip2cs
  save naif
  save naig
  save naip1
  save naip2
  save x2sml
  save x32sml
  save x3sml
  save xbig

  data aifcs(  1) /     0.10527461226531408809E+00 /
  data aifcs(  2) /     0.01183613628152997844E+00 /
  data aifcs(  3) /     0.00012328104173225664E+00 /
  data aifcs(  4) /     0.00000062261225638140E+00 /
  data aifcs(  5) /     0.00000000185298887844E+00 /
  data aifcs(  6) /     0.00000000000363328873E+00 /
  data aifcs(  7) /     0.00000000000000504622E+00 /
  data aifcs(  8) /     0.00000000000000000522E+00 /

  data aigcs(  1) /     0.021233878150918666852E+00 /
  data aigcs(  2) /     0.086315930335214406752E+00 /
  data aigcs(  3) /     0.001797594720383231358E+00 /
  data aigcs(  4) /     0.000014265499875550693E+00 /
  data aigcs(  5) /     0.000000059437995283683E+00 /
  data aigcs(  6) /     0.000000000152403366479E+00 /
  data aigcs(  7) /     0.000000000000264587660E+00 /
  data aigcs(  8) /     0.000000000000000331562E+00 /
  data aigcs(  9) /     0.000000000000000000314E+00 /

  data aip2cs(  1) /     0.0065457691989713757E+00 /
  data aip2cs(  2) /     0.0023833724120774592E+00 /
  data aip2cs(  3) /    -0.0000430700770220586E+00 /
  data aip2cs(  4) /     0.0000015629125858629E+00 /
  data aip2cs(  5) /    -0.0000000815417186163E+00 /
  data aip2cs(  6) /     0.0000000054103738057E+00 /
  data aip2cs(  7) /    -0.0000000004284130883E+00 /
  data aip2cs(  8) /     0.0000000000389497963E+00 /
  data aip2cs(  9) /    -0.0000000000039623161E+00 /
  data aip2cs( 10) /     0.0000000000004428184E+00 /
  data aip2cs( 11) /    -0.0000000000000536297E+00 /
  data aip2cs( 12) /     0.0000000000000069650E+00 /
  data aip2cs( 13) /    -0.0000000000000009620E+00 /
  data aip2cs( 14) /     0.0000000000000001403E+00 /
  data aip2cs( 15) /    -0.0000000000000000215E+00 /

  data aip1cs(  1) /     0.0358865097808301538E+00 /
  data aip1cs(  2) /     0.0114668575627764899E+00 /
  data aip1cs(  3) /    -0.0007592073583861400E+00 /
  data aip1cs(  4) /     0.0000869517610893841E+00 /
  data aip1cs(  5) /    -0.0000128237294298592E+00 /
  data aip1cs(  6) /     0.0000022062695681038E+00 /
  data aip1cs(  7) /    -0.0000004222295185921E+00 /
  data aip1cs(  8) /     0.0000000874686415726E+00 /
  data aip1cs(  9) /    -0.0000000192773588418E+00 /
  data aip1cs( 10) /     0.0000000044668460054E+00 /
  data aip1cs( 11) /    -0.0000000010790108052E+00 /
  data aip1cs( 12) /     0.0000000002700029447E+00 /
  data aip1cs( 13) /    -0.0000000000696480108E+00 /
  data aip1cs( 14) /     0.0000000000184489907E+00 /
  data aip1cs( 15) /    -0.0000000000050027817E+00 /
  data aip1cs( 16) /     0.0000000000013852243E+00 /
  data aip1cs( 17) /    -0.0000000000003908218E+00 /
  data aip1cs( 18) /     0.0000000000001121536E+00 /
  data aip1cs( 19) /    -0.0000000000000326862E+00 /
  data aip1cs( 20) /     0.0000000000000096619E+00 /
  data aip1cs( 21) /    -0.0000000000000028935E+00 /
  data aip1cs( 22) /     0.0000000000000008770E+00 /
  data aip1cs( 23) /    -0.0000000000000002688E+00 /
  data aip1cs( 24) /     0.0000000000000000832E+00 /
  data aip1cs( 25) /    -0.0000000000000000260E+00 /

  data naif / 0 /
  data naig / 0 /
  data naip1 / 0 /
  data naip2 / 0 /
  data x2sml / 0.0E+00 /
  data x32sml / 0.0E+00 /
  data x3sml / 0.0E+00 /
  data xbig / 0.0E+00 /

  if ( naif == 0 ) then
    eta = 0.1E+00 * r4_mach ( 3 )
    naif = r4_inits ( aifcs, 8, eta )
    naig = r4_inits ( aigcs, 9, eta )
    naip1 = r4_inits ( aip1cs, 25, eta )
    naip2 = r4_inits ( aip2cs, 15, eta )
    x2sml = sqrt ( eta )
    x3sml = eta ** 0.3333E+00
    x32sml = 1.3104E+00 * x3sml * x3sml
    xbig = r4_mach ( 2 ) ** 0.6666E+00
  end if

  if ( x < - 1.0E+00 ) then
    call r4_admp ( x, xn, phi )
    r4_aide = xn * cos ( phi )
  else if ( abs ( x ) <= x2sml ) then
    x2 = 0.0E+00
    x3 = 0.0E+00
    r4_aide = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) &
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00
  else if ( abs ( x ) <= x3sml ) then
    x2 = x * x
    x3 = 0.0E+00
    r4_aide = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) &
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00
  else if ( abs ( x ) <= x32sml ) then
    x2 = x * x
    x3 = x * x * x
    r4_aide = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) &
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00
  else if ( x <= 1.0E+00 ) then
    x2 = x * x
    x3 = x * x * x
    r4_aide = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) &
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00
    r4_aide = r4_aide * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  else if ( x <= 4.0E+00 ) then
    sqrtx = sqrt ( x )
    z = ( 16.0E+00 / ( x * sqrtx ) - 9.0E+00 ) / 7.0E+00
    r4_aide = ( - 0.28125E+00 &
      - r4_csevl ( z, aip1cs, naip1 ) ) * sqrt ( sqrtx )
  else if ( x < xbig ) then
    sqrtx = sqrt ( x )
    z = 16.0E+00 / ( x * sqrtx ) - 1.0E+00
    r4_aide = ( - 0.28125E+00 &
      - r4_csevl ( z, aip2cs, naip2 ) ) * sqrt ( sqrtx )
  else
    sqrtx = sqrt ( x )
    z = - 1.0E+00
    r4_aide = ( - 0.28125E+00 &
      - r4_csevl ( z, aip2cs, naip2 ) ) * sqrt ( sqrtx )
  end if

  return
end
function r4_aie ( x )

!*****************************************************************************80
!
!! R4_AIE evaluates the exponential scaled Airy function Ai of an R4 argument.
!
!  Discussion:
!
!    If X <= 0
!      R4_AIE ( X ) = R4_AI ( X )
!    else
!      R4_AIE ( X ) = R4_AI ( X ) * exp ( 2/3 X^(3/2) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_AIE, the Airy function Ai of X.
!
  implicit none

  real ( kind = 4 ) aifcs(9)
  real ( kind = 4 ) aigcs(8)
  real ( kind = 4 ) aipcs(34)
  real ( kind = 4 ) eta
  integer ( kind = 4 ) naif
  integer ( kind = 4 ) naig
  integer ( kind = 4 ) naip
  real ( kind = 4 ) r4_aie
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sqrtx
  real ( kind = 4 ) theta
  real ( kind = 4 ) x
  real ( kind = 4 ) x32sml
  real ( kind = 4 ) x3sml
  real ( kind = 4 ) xbig
  real ( kind = 4 ) xm
  real ( kind = 4 ) z

  save aifcs
  save aigcs
  save aipcs
  save naif
  save naig
  save naip
  save x32sml
  save x3sml
  save xbig

  data aifcs( 1) /   -0.03797135849666999750E+00 /
  data aifcs( 2) /    0.05919188853726363857E+00 /
  data aifcs( 3) /    0.00098629280577279975E+00 /
  data aifcs( 4) /    0.00000684884381907656E+00 /
  data aifcs( 5) /    0.00000002594202596219E+00 /
  data aifcs( 6) /    0.00000000006176612774E+00 /
  data aifcs( 7) /    0.00000000000010092454E+00 /
  data aifcs( 8) /    0.00000000000000012014E+00 /
  data aifcs( 9) /    0.00000000000000000010E+00 /

  data aigcs( 1) /    0.01815236558116127E+00 /
  data aigcs( 2) /    0.02157256316601076E+00 /
  data aigcs( 3) /    0.00025678356987483E+00 /
  data aigcs( 4) /    0.00000142652141197E+00 /
  data aigcs( 5) /    0.00000000457211492E+00 /
  data aigcs( 6) /    0.00000000000952517E+00 /
  data aigcs( 7) /    0.00000000000001392E+00 /
  data aigcs( 8) /    0.00000000000000001E+00 /

  data aipcs( 1) /   -0.0187519297793868E+00 /
  data aipcs( 2) /   -0.0091443848250055E+00 /
  data aipcs( 3) /    0.0009010457337825E+00 /
  data aipcs( 4) /   -0.0001394184127221E+00 /
  data aipcs( 5) /    0.0000273815815785E+00 /
  data aipcs( 6) /   -0.0000062750421119E+00 /
  data aipcs( 7) /    0.0000016064844184E+00 /
  data aipcs( 8) /   -0.0000004476392158E+00 /
  data aipcs( 9) /    0.0000001334635874E+00 /
  data aipcs(10) /   -0.0000000420735334E+00 /
  data aipcs(11) /    0.0000000139021990E+00 /
  data aipcs(12) /   -0.0000000047831848E+00 /
  data aipcs(13) /    0.0000000017047897E+00 /
  data aipcs(14) /   -0.0000000006268389E+00 /
  data aipcs(15) /    0.0000000002369824E+00 /
  data aipcs(16) /   -0.0000000000918641E+00 /
  data aipcs(17) /    0.0000000000364278E+00 /
  data aipcs(18) /   -0.0000000000147475E+00 /
  data aipcs(19) /    0.0000000000060851E+00 /
  data aipcs(20) /   -0.0000000000025552E+00 /
  data aipcs(21) /    0.0000000000010906E+00 /
  data aipcs(22) /   -0.0000000000004725E+00 /
  data aipcs(23) /    0.0000000000002076E+00 /
  data aipcs(24) /   -0.0000000000000924E+00 /
  data aipcs(25) /    0.0000000000000417E+00 /
  data aipcs(26) /   -0.0000000000000190E+00 /
  data aipcs(27) /    0.0000000000000087E+00 /
  data aipcs(28) /   -0.0000000000000040E+00 /
  data aipcs(29) /    0.0000000000000019E+00 /
  data aipcs(30) /   -0.0000000000000009E+00 /
  data aipcs(31) /    0.0000000000000004E+00 /
  data aipcs(32) /   -0.0000000000000002E+00 /
  data aipcs(33) /    0.0000000000000001E+00 /
  data aipcs(34) /   -0.0000000000000000E+00 /

  data naif / 0 /
  data naig / 0 /
  data naip / 0 /
  data x3sml / 0.0E+00 /
  data x32sml / 0.0E+00 /
  data xbig / 0.0E+00 /

  if ( naif == 0 ) then
    eta = 0.1E+00 * r4_mach ( 3 )
    naif = r4_inits ( aifcs, 9, eta )
    naig = r4_inits ( aigcs, 8, eta )
    naip = r4_inits ( aipcs, 34, eta )
    x3sml = eta ** 0.3333E+00
    x32sml = 1.3104E+00 * x3sml * x3sml
    xbig = r4_mach ( 2 ) ** 0.6666E+00
  end if

  if ( x < - 1.0E+00 ) then
    call r4_aimp ( x, xm, theta )
    r4_aie = xm * cos ( theta )
  else if ( abs ( x ) <= x32sml ) then
    z = 0.0E+00
    r4_aie = 0.375E+00 + ( r4_csevl ( z, aifcs, naif ) &
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) )
  else if ( abs ( x ) <= x3sml ) then
    z = 0.0E+00
    r4_aie = 0.375E+00 + ( r4_csevl ( z, aifcs, naif ) &
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) )
    r4_aie = r4_aie * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  else if ( x <= 1.0E+00 ) then
    z = x * x * x
    r4_aie = 0.375E+00 + ( r4_csevl ( z, aifcs, naif ) &
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) )
    r4_aie = r4_aie * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  else if ( x < xbig ) then
    sqrtx = sqrt ( x )
    z = 2.0E+00 / ( x * sqrtx ) - 1.0E+00
    r4_aie = ( 0.28125E+00 + r4_csevl ( z, aipcs, naip ) ) / sqrt ( sqrtx )
  else
    sqrtx = sqrt ( x )
    z = - 1.0E+00
    r4_aie = ( 0.28125E+00 + r4_csevl ( z, aipcs, naip ) ) / sqrt ( sqrtx )
  end if

  return
end
subroutine r4_aimp ( x, ampl, theta )

!*****************************************************************************80
!
!! R4_AIMP evaluates the modulus and phase of the Airy function.
!
!  Description:
!
!    This function must only be called when X <= -1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) AMPL, PHI, the modulus and phase of the 
!    Airy function.
!
  implicit none

  real ( kind = 4 ) am21cs(40)
  real ( kind = 4 ) am22cs(33)
  real ( kind = 4 ) ampl
  real ( kind = 4 ) ath1cs(36)
  real ( kind = 4 ) ath2cs(32)
  real ( kind = 4 ) eta
  integer ( kind = 4 ) nam21
  integer ( kind = 4 ) nam22
  integer ( kind = 4 ) nath1
  integer ( kind = 4 ) nath2
  real ( kind = 4 ) phi
  real ( kind = 4 ) pi4
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sqrtx
  real ( kind = 4 ) theta
  real ( kind = 4 ) x
  real ( kind = 4 ) xsml
  real ( kind = 4 ) z

  save am21cs
  save am22cs
  save ath1cs
  save ath2cs
  save nam21
  save nam22
  save nath1
  save nath2
  save pi4
  save xsml

  data am21cs( 1) /   +0.0065809191761485E+00 /
  data am21cs( 2) /   +0.0023675984685722E+00 /
  data am21cs( 3) /   +0.0001324741670371E+00 /
  data am21cs( 4) /   +0.0000157600904043E+00 /
  data am21cs( 5) /   +0.0000027529702663E+00 /
  data am21cs( 6) /   +0.0000006102679017E+00 /
  data am21cs( 7) /   +0.0000001595088468E+00 /
  data am21cs( 8) /   +0.0000000471033947E+00 /
  data am21cs( 9) /   +0.0000000152933871E+00 /
  data am21cs(10) /   +0.0000000053590722E+00 /
  data am21cs(11) /   +0.0000000020000910E+00 /
  data am21cs(12) /   +0.0000000007872292E+00 /
  data am21cs(13) /   +0.0000000003243103E+00 /
  data am21cs(14) /   +0.0000000001390106E+00 /
  data am21cs(15) /   +0.0000000000617011E+00 /
  data am21cs(16) /   +0.0000000000282491E+00 /
  data am21cs(17) /   +0.0000000000132979E+00 /
  data am21cs(18) /   +0.0000000000064188E+00 /
  data am21cs(19) /   +0.0000000000031697E+00 /
  data am21cs(20) /   +0.0000000000015981E+00 /
  data am21cs(21) /   +0.0000000000008213E+00 /
  data am21cs(22) /   +0.0000000000004296E+00 /
  data am21cs(23) /   +0.0000000000002284E+00 /
  data am21cs(24) /   +0.0000000000001232E+00 /
  data am21cs(25) /   +0.0000000000000675E+00 /
  data am21cs(26) /   +0.0000000000000374E+00 /
  data am21cs(27) /   +0.0000000000000210E+00 /
  data am21cs(28) /   +0.0000000000000119E+00 /
  data am21cs(29) /   +0.0000000000000068E+00 /
  data am21cs(30) /   +0.0000000000000039E+00 /
  data am21cs(31) /   +0.0000000000000023E+00 /
  data am21cs(32) /   +0.0000000000000013E+00 /
  data am21cs(33) /   +0.0000000000000008E+00 /
  data am21cs(34) /   +0.0000000000000005E+00 /
  data am21cs(35) /   +0.0000000000000003E+00 /
  data am21cs(36) /   +0.0000000000000001E+00 /
  data am21cs(37) /   +0.0000000000000001E+00 /
  data am21cs(38) /   +0.0000000000000000E+00 /
  data am21cs(39) /   +0.0000000000000000E+00 /
  data am21cs(40) /   +0.0000000000000000E+00 /

  data ath1cs( 1) /   -0.07125837815669365E+00 /
  data ath1cs( 2) /   -0.00590471979831451E+00 /
  data ath1cs( 3) /   -0.00012114544069499E+00 /
  data ath1cs( 4) /   -0.00000988608542270E+00 /
  data ath1cs( 5) /   -0.00000138084097352E+00 /
  data ath1cs( 6) /   -0.00000026142640172E+00 /
  data ath1cs( 7) /   -0.00000006050432589E+00 /
  data ath1cs( 8) /   -0.00000001618436223E+00 /
  data ath1cs( 9) /   -0.00000000483464911E+00 /
  data ath1cs(10) /   -0.00000000157655272E+00 /
  data ath1cs(11) /   -0.00000000055231518E+00 /
  data ath1cs(12) /   -0.00000000020545441E+00 /
  data ath1cs(13) /   -0.00000000008043412E+00 /
  data ath1cs(14) /   -0.00000000003291252E+00 /
  data ath1cs(15) /   -0.00000000001399875E+00 /
  data ath1cs(16) /   -0.00000000000616151E+00 /
  data ath1cs(17) /   -0.00000000000279614E+00 /
  data ath1cs(18) /   -0.00000000000130428E+00 /
  data ath1cs(19) /   -0.00000000000062373E+00 /
  data ath1cs(20) /   -0.00000000000030512E+00 /
  data ath1cs(21) /   -0.00000000000015239E+00 /
  data ath1cs(22) /   -0.00000000000007758E+00 /
  data ath1cs(23) /   -0.00000000000004020E+00 /
  data ath1cs(24) /   -0.00000000000002117E+00 /
  data ath1cs(25) /   -0.00000000000001132E+00 /
  data ath1cs(26) /   -0.00000000000000614E+00 /
  data ath1cs(27) /   -0.00000000000000337E+00 /
  data ath1cs(28) /   -0.00000000000000188E+00 /
  data ath1cs(29) /   -0.00000000000000105E+00 /
  data ath1cs(30) /   -0.00000000000000060E+00 /
  data ath1cs(31) /   -0.00000000000000034E+00 /
  data ath1cs(32) /   -0.00000000000000020E+00 /
  data ath1cs(33) /   -0.00000000000000011E+00 /
  data ath1cs(34) /   -0.00000000000000007E+00 /
  data ath1cs(35) /   -0.00000000000000004E+00 /
  data ath1cs(36) /   -0.00000000000000002E+00 /

  data am22cs( 1) /   -0.01562844480625341E+00 /
  data am22cs( 2) /   +0.00778336445239681E+00 /
  data am22cs( 3) /   +0.00086705777047718E+00 /
  data am22cs( 4) /   +0.00015696627315611E+00 /
  data am22cs( 5) /   +0.00003563962571432E+00 /
  data am22cs( 6) /   +0.00000924598335425E+00 /
  data am22cs( 7) /   +0.00000262110161850E+00 /
  data am22cs( 8) /   +0.00000079188221651E+00 /
  data am22cs( 9) /   +0.00000025104152792E+00 /
  data am22cs(10) /   +0.00000008265223206E+00 /
  data am22cs(11) /   +0.00000002805711662E+00 /
  data am22cs(12) /   +0.00000000976821090E+00 /
  data am22cs(13) /   +0.00000000347407923E+00 /
  data am22cs(14) /   +0.00000000125828132E+00 /
  data am22cs(15) /   +0.00000000046298826E+00 /
  data am22cs(16) /   +0.00000000017272825E+00 /
  data am22cs(17) /   +0.00000000006523192E+00 /
  data am22cs(18) /   +0.00000000002490471E+00 /
  data am22cs(19) /   +0.00000000000960156E+00 /
  data am22cs(20) /   +0.00000000000373448E+00 /
  data am22cs(21) /   +0.00000000000146417E+00 /
  data am22cs(22) /   +0.00000000000057826E+00 /
  data am22cs(23) /   +0.00000000000022991E+00 /
  data am22cs(24) /   +0.00000000000009197E+00 /
  data am22cs(25) /   +0.00000000000003700E+00 /
  data am22cs(26) /   +0.00000000000001496E+00 /
  data am22cs(27) /   +0.00000000000000608E+00 /
  data am22cs(28) /   +0.00000000000000248E+00 /
  data am22cs(29) /   +0.00000000000000101E+00 /
  data am22cs(30) /   +0.00000000000000041E+00 /
  data am22cs(31) /   +0.00000000000000017E+00 /
  data am22cs(32) /   +0.00000000000000007E+00 /
  data am22cs(33) /   +0.00000000000000002E+00 /

  data ath2cs( 1) /   +0.00440527345871877E+00 /
  data ath2cs( 2) /   -0.03042919452318455E+00 /
  data ath2cs( 3) /   -0.00138565328377179E+00 /
  data ath2cs( 4) /   -0.00018044439089549E+00 /
  data ath2cs( 5) /   -0.00003380847108327E+00 /
  data ath2cs( 6) /   -0.00000767818353522E+00 /
  data ath2cs( 7) /   -0.00000196783944371E+00 /
  data ath2cs( 8) /   -0.00000054837271158E+00 /
  data ath2cs( 9) /   -0.00000016254615505E+00 /
  data ath2cs(10) /   -0.00000005053049981E+00 /
  data ath2cs(11) /   -0.00000001631580701E+00 /
  data ath2cs(12) /   -0.00000000543420411E+00 /
  data ath2cs(13) /   -0.00000000185739855E+00 /
  data ath2cs(14) /   -0.00000000064895120E+00 /
  data ath2cs(15) /   -0.00000000023105948E+00 /
  data ath2cs(16) /   -0.00000000008363282E+00 /
  data ath2cs(17) /   -0.00000000003071196E+00 /
  data ath2cs(18) /   -0.00000000001142367E+00 /
  data ath2cs(19) /   -0.00000000000429811E+00 /
  data ath2cs(20) /   -0.00000000000163389E+00 /
  data ath2cs(21) /   -0.00000000000062693E+00 /
  data ath2cs(22) /   -0.00000000000024260E+00 /
  data ath2cs(23) /   -0.00000000000009461E+00 /
  data ath2cs(24) /   -0.00000000000003716E+00 /
  data ath2cs(25) /   -0.00000000000001469E+00 /
  data ath2cs(26) /   -0.00000000000000584E+00 /
  data ath2cs(27) /   -0.00000000000000233E+00 /
  data ath2cs(28) /   -0.00000000000000093E+00 /
  data ath2cs(29) /   -0.00000000000000037E+00 /
  data ath2cs(30) /   -0.00000000000000015E+00 /
  data ath2cs(31) /   -0.00000000000000006E+00 /
  data ath2cs(32) /   -0.00000000000000002E+00 /

  data nam21 / 0 /
  data nam22 / 0 /
  data nath1 / 0 /
  data nath2 / 0 /
  data pi4 / 0.78539816339744831E+00 /
  data xsml / 0.0E+00 /

  if ( nam21 == 0 ) then
    eta = 0.1E+00 * r4_mach ( 3 )
    nam21 = r4_inits ( am21cs, 40, eta )
    nath1 = r4_inits ( ath1cs, 36, eta )
    nam22 = r4_inits ( am22cs, 33, eta )
    nath2 = r4_inits ( ath2cs, 32, eta )
    xsml = - ( 16.0E+00 / r4_mach ( 3 ) ) ** 0.3333E+00
  end if

  if ( x <= xsml ) then
    z = 1.0E+00
    ampl = 0.3125E+00 + r4_csevl ( z, am21cs, nam21 )
    theta = - 0.625E+00 + r4_csevl ( z, ath1cs, nath1 )
  else if ( x < - 2.0E+00 ) then
    z = 16.0 / x / x / x + 1.0
    ampl = 0.3125E+00 + r4_csevl ( z, am21cs, nam21 )
    theta = - 0.625E+00 + r4_csevl ( z, ath1cs, nath1 )
  else if ( x <= - 1.0E+00 ) then
    z = ( 16.0E+00 / x / x / x + 9.0E+00 ) / 7.0E+00
    ampl = 0.3125E+00 + r4_csevl ( z, am22cs, nam22 )
    theta = - 0.625E+00 + r4_csevl ( z, ath2cs, nath2 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_AIMP - Fatal error!'
    write ( *, '(a)' ) '  - 1.0 < X.'
    stop
  end if

  sqrtx = sqrt ( - x )
  ampl = sqrt ( ampl / sqrtx )
  theta = pi4 - x * sqrtx * theta

  return
end
function r4_aint ( x )

!*****************************************************************************80
!
!! R4_AINT truncates an R4 argument to an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_AINT, the truncated version of X.
!
  implicit none

  real ( kind = 4 ) r4_aint
  real ( kind = 4 ) value
  real ( kind = 4 ) x

  if ( x < 0.0E+00 ) then
    value = - real ( int ( abs ( x ) ), kind = 4 )
  else
    value =   real ( int ( abs ( x ) ), kind = 4 )
  end if

  r4_aint = value

  return
end
function r4_asin ( x )

!*****************************************************************************80
!
!! R4_ASIN evaluates the arc-sine of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_ASIN, the arc-sine of X.
!
  implicit none

  real ( kind = 4 ) asincs(20)
  integer ( kind = 4 ) nterms
  real ( kind = 4 ), parameter :: pi2 = 1.57079632679489661923E+00
  real ( kind = 4 ) r4_asin
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  save asincs
  save nterms

  data asincs( 1) / 0.10246391753227159E+00 /
  data asincs( 2) / 0.054946487221245833E+00 /
  data asincs( 3) / 0.004080630392544969E+00 /
  data asincs( 4) / 0.000407890068546044E+00 /
  data asincs( 5) / 0.000046985367432203E+00 /
  data asincs( 6) / 0.000005880975813970E+00 /
  data asincs( 7) / 0.000000777323124627E+00 /
  data asincs( 8) / 0.000000106774233400E+00 /
  data asincs( 9) / 0.000000015092399536E+00 /
  data asincs(10) / 0.000000002180972408E+00 /
  data asincs(11) / 0.000000000320759842E+00 /
  data asincs(12) / 0.000000000047855369E+00 /
  data asincs(13) / 0.000000000007225128E+00 /
  data asincs(14) / 0.000000000001101833E+00 /
  data asincs(15) / 0.000000000000169476E+00 /
  data asincs(16) / 0.000000000000026261E+00 /
  data asincs(17) / 0.000000000000004095E+00 /
  data asincs(18) / 0.000000000000000642E+00 /
  data asincs(19) / 0.000000000000000101E+00 /
  data asincs(20) / 0.000000000000000016E+00 /

  data nterms / 0 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( asincs, 20, 0.1E+00 * r4_mach ( 3 ) )
    sqeps = sqrt ( 6.0E+00 * r4_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( x < - 1.0E+00 - sqeps ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_ASIN - Fatal error!'
    write ( *, '(a)' ) '  X < - 1.0'
    stop

  else if ( x < - 1.0E+00 ) then

    value = - pi2

  else if ( x < 1.0E+00 ) then

    z = 0.0E+00
    if ( sqeps < y ) then
      z = y * y
    end if

    if ( z <= 0.5E+00 ) then
      value = x * ( 1.0E+00 + r4_csevl ( 4.0E+00 * z - 1.0E+00, &
        asincs, nterms ) )
    else
      value = pi2 - sqrt ( 1.0E+00 - z ) * ( 1.0E+00 + &
        r4_csevl ( 3.0E+00 - 4.0E+00 * z, asincs, nterms ) )
    end if

    if ( x < 0.0E+00 ) then
      value = - abs ( value )
    else if ( 0.0E+00 < x ) then
      value = + abs ( value )
    end if

  else if ( x < 1.0E+00 + sqeps ) then

    value = pi2

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_ASIN - Fatal error!'
    write ( *, '(a)' ) '  1.0 < X'
    stop

  end if

  r4_asin = value

  return
end
function r4_asinh ( x )

!*****************************************************************************80
!
!! R4_ASINH evaluates the arc-sine of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_ASINH, the arc-hyperbolic sine of X.
!
  implicit none

  real ( kind = 4 ) aln2
  real ( kind = 4 ) asnhcs(20)
  integer ( kind = 4 ) nterms
  real ( kind = 4 ) r4_asinh
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) y

  save aln2
  save asnhcs
  save nterms
  save xmax

  data aln2 / 0.69314718055994530942E+00 /

  data asnhcs( 1) /   -0.12820039911738186E+00 /
  data asnhcs( 2) /   -0.058811761189951768E+00 /
  data asnhcs( 3) /    0.004727465432212481E+00 /
  data asnhcs( 4) /   -0.000493836316265361E+00 /
  data asnhcs( 5) /    0.000058506207058557E+00 /
  data asnhcs( 6) /   -0.000007466998328931E+00 /
  data asnhcs( 7) /    0.000001001169358355E+00 /
  data asnhcs( 8) /   -0.000000139035438587E+00 /
  data asnhcs( 9) /    0.000000019823169483E+00 /
  data asnhcs(10) /   -0.000000002884746841E+00 /
  data asnhcs(11) /    0.000000000426729654E+00 /
  data asnhcs(12) /   -0.000000000063976084E+00 /
  data asnhcs(13) /    0.000000000009699168E+00 /
  data asnhcs(14) /   -0.000000000001484427E+00 /
  data asnhcs(15) /    0.000000000000229037E+00 /
  data asnhcs(16) /   -0.000000000000035588E+00 /
  data asnhcs(17) /    0.000000000000005563E+00 /
  data asnhcs(18) /   -0.000000000000000874E+00 /
  data asnhcs(19) /    0.000000000000000138E+00 /
  data asnhcs(20) /   -0.000000000000000021E+00 /

  data nterms / 0 /
  data xmax / 0.0E+00 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( asnhcs, 20, 0.1E+00 * r4_mach ( 3 ) )
    sqeps = sqrt ( r4_mach ( 3 ) )
    xmax = 1.0E+00 / sqeps
  end if

  y = abs ( x )

  if ( y <= 1.0E+00 ) then

    value = x
    if ( sqeps < y ) then
      value = x * ( 1.0E+00 &
        + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, asnhcs, nterms ) )
    end if

  else

    if ( y < xmax ) then
      value = log ( y + sqrt ( y * y + 1.0E+00 ) )
    else
      value = aln2 + log ( y )
    end if

    if ( x < 0.0E+00 ) then 
      value = - value
    end if

  end if

  r4_asinh = value

  return
end
function r4_atan ( x )

!*****************************************************************************80
!
!! R4_ATAN evaluates the arc-tangent of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_ATAN, the arc-tangent of X.
!
  implicit none

  real ( kind = 4 ) atancs(9)
  real ( kind = 4 ) conpi8(4)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nterms
  real ( kind = 4 ) pi8(4)
  real ( kind = 4 ) r4_atan
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) t
  real ( kind = 4 ) tanp8(3)
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) xbig
  real ( kind = 4 ) xbnd1
  real ( kind = 4 ) xbnd2
  real ( kind = 4 ) xbnd3
  real ( kind = 4 ) xbnd4
  real ( kind = 4 ) y

  save atancs
  save conpi8
  save nterms
  save pi8
  save tanp8
  save xbig
  save xbnd1
  save xbnd2
  save xbnd3
  save xbnd4

  data atancs( 1) / 0.48690110349241406E+00 /
  data atancs( 2) / -0.006510831636717464E+00 /
  data atancs( 3) / 0.000038345828265245E+00 /
  data atancs( 4) / -0.000000268722128762E+00 /
  data atancs( 5) / 0.000000002050093098E+00 /
  data atancs( 6) / -0.000000000016450717E+00 /
  data atancs( 7) / 0.000000000000136509E+00 /
  data atancs( 8) / -0.000000000000001160E+00 /
  data atancs( 9) / 0.000000000000000010E+00 /

  data xbnd1 / +0.198912367379658006E+00 /
  data xbnd2 / +0.668178637919298919E+00 /
  data xbnd3 / +1.49660576266548901E+00 /
  data xbnd4 / +5.02733949212584810E+00 /

  data tanp8(1) / 0.414213562373095048E+00 /
  data tanp8(2) / 1.0E+00 /
  data tanp8(3) / 2.41421356237309504E+00 /

  data conpi8(1) / 0.375E+00 /
  data conpi8(2) / 0.75E+00 /
  data conpi8(3) / 1.125E+00 /
  data conpi8(4) / 1.5E+00 /

  data pi8(1) / +0.176990816987241548E-01 /
  data pi8(2) / +0.353981633974483096E-01 /
  data pi8(3) / +0.530972450961724644E-01 /
  data pi8(4) / 0.0707963267948966192E+00 /

  data nterms / 0 /
  data xbig / 0.0E+00 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( atancs, 9, 0.1E+00 * r4_mach ( 3 ) )
    sqeps = sqrt ( 6.0E+00 * r4_mach ( 3 ) )
    xbig = 1.0E+00 / r4_mach ( 3 )
  end if

  y = abs ( x )

  if ( y <= xbnd1 ) then

    value = x
    if ( sqeps < y ) then
      value = x * ( 0.75E+00 + r4_csevl ( &
        50.0E+00 * y * y - 1.0E+00, atancs, nterms ) )
    end if

  else if ( y <= xbnd4 ) then

    if ( xbnd3 < y ) then
      n = 3
    else if ( xbnd2 < y ) then
      n = 2
    else
      n = 1
    end if

    t = ( y - tanp8(n) ) / ( 1.0E+00 + y * tanp8(n) )

    value = conpi8(n) + ( pi8(n) + t * ( 0.75E+00 + &
      r4_csevl ( 50.0E+00 * t * t - 1.0E+00, atancs, nterms ) ) )

  else

    value = conpi8(4) + pi8(4)

    if ( y < xbig ) then
      value = conpi8(4) + ( pi8(4) - ( 0.75E+00 + &
        r4_csevl ( 50.0E+00 / y / y - 1.0E+00, atancs, &
        nterms ) ) / y )
    end if

  end if

  if ( x < 0.0E+00 ) then
    value = - abs ( value )
  else
    value = + abs ( value )
  end if

  r4_atan = value

  return
end
function r4_atan2 ( sn, cs )

!*****************************************************************************80
!
!! R4_ATAN2 evaluates the arc-tangent of two R4 arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) SN, CS, the Y and X coordinates of a point
!    on the angle.
!
!    Output, real ( kind = 4 ) R4_ATAN2, the arc-tangent of the angle.
!
  implicit none

  real ( kind = 4 ) abscs
  real ( kind = 4 ) abssn
  real ( kind = 4 ) big
  real ( kind = 4 ) cs
  real ( kind = 4 ), save :: pi = 3.14159265358979323846E+00
  real ( kind = 4 ) r4_atan2
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sml
  real ( kind = 4 ) sn

  save big
  save sml

  data sml / 0.0E+00 /
  data big / 0.0E+00 /

  if ( sml == 0.0E+00 ) then
    sml = r4_mach ( 1 )
    big = r4_mach ( 2 )
  end if
!
!  We now make sure SN can be divided by CS.  It is painful.
!
  abssn = abs ( sn )
  abscs = abs ( cs )

  if ( abscs <= abssn ) then

    if ( abscs < 1.0E+00 .and. abscs * big <= abssn ) then

      if ( sn < 0.0E+00 ) then
        r4_atan2 = - 0.5E+00 * pi
      else if ( sn == 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_ATAN2 - Fatal error!'
        write ( *, '(a)' ) '  Both arguments are 0.'
        stop
      else
        r4_atan2 = 0.5E+00 * pi
      end if

      return

    end if

  else

    if ( 1.0E+00 < abscs .and. abssn <= abscs * sml ) then

      if ( 0.0E+00 <= cs ) then
        r4_atan2 = 0.0E+00
      else
        r4_atan2 = pi
      end if

      return

    end if

  end if

  r4_atan2 = atan ( sn / cs )

  if ( cs < 0.0E+00 ) then
    r4_atan2 = r4_atan2 + pi
  end if

  if ( pi < r4_atan2 ) then
    r4_atan2 = r4_atan2 - 2.0E+00 * pi
  end if

  return
end
function r4_atanh ( x )

!*****************************************************************************80
!
!! R4_ATANH evaluates the arc-hyperbolic tangent of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_ATANH, the arc-hyperbolic tangent of X.
!
  implicit none

  real ( kind = 4 ) atnhcs(15)
  real ( kind = 4 ) dxrel
  integer ( kind = 4 ) nterms
  real ( kind = 4 ) r4_atanh
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) y

  save atnhcs
  save dxrel
  save nterms

  data atnhcs( 1) /    0.094395102393195492E+00 /
  data atnhcs( 2) /    0.049198437055786159E+00 /
  data atnhcs( 3) /    0.002102593522455432E+00 /
  data atnhcs( 4) /    0.000107355444977611E+00 /
  data atnhcs( 5) /    0.000005978267249293E+00 /
  data atnhcs( 6) /    0.000000350506203088E+00 /
  data atnhcs( 7) /    0.000000021263743437E+00 /
  data atnhcs( 8) /    0.000000001321694535E+00 /
  data atnhcs( 9) /    0.000000000083658755E+00 /
  data atnhcs(10) /    0.000000000005370503E+00 /
  data atnhcs(11) /    0.000000000000348665E+00 /
  data atnhcs(12) /    0.000000000000022845E+00 /
  data atnhcs(13) /    0.000000000000001508E+00 /
  data atnhcs(14) /    0.000000000000000100E+00 /
  data atnhcs(15) /    0.000000000000000006E+00 /

  data nterms / 0 /
  data dxrel / 0.0E+00 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( atnhcs, 15, 0.1E+00 * r4_mach ( 3 ) )
    dxrel = sqrt ( r4_mach ( 4 ) )
    sqeps = sqrt ( 3.0E+00 * r4_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    value = x
  else if ( y <= 0.5E+00 ) then
    value = x * ( 1.0E+00 &
      + r4_csevl ( 8.0E+00 * x * x - 1.0E+00, atnhcs, nterms ) )
  else if ( y < 1.0E+00 ) then
    value = 0.5E+00 * log ( ( 1.0E+00 + x ) / ( 1.0E+00 - x ) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_ATANH - Fatal error!'
    write ( *, '(a)' ) '  1 <= |X|.'
    stop
  end if

  r4_atanh = value

  return
end
function r4_besi0 ( x )

!*****************************************************************************80
!
!! R4_BESI0 evaluates the Bessel function I of order 0 of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESI0, the Bessel function I of order 0 of X.
!
  implicit none

  real ( kind = 4 ) bi0cs(12)
  integer ( kind = 4 ) nti0
  real ( kind = 4 ) r4_besi0
  real ( kind = 4 ) r4_besi0e
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save bi0cs
  save nti0
  save xmax
  save xsml

  data bi0cs( 1) /   -0.07660547252839144951E+00/
  data bi0cs( 2) /    1.927337953993808270E+00/
  data bi0cs( 3) /    0.2282644586920301339E+00/
  data bi0cs( 4) /    0.01304891466707290428E+00/
  data bi0cs( 5) /    0.00043442709008164874E+00/
  data bi0cs( 6) /    0.00000942265768600193E+00/
  data bi0cs( 7) /    0.00000014340062895106E+00/
  data bi0cs( 8) /    0.00000000161384906966E+00/
  data bi0cs( 9) /    0.00000000001396650044E+00/
  data bi0cs(10) /    0.00000000000009579451E+00/
  data bi0cs(11) /    0.00000000000000053339E+00/
  data bi0cs(12) /    0.00000000000000000245E+00/

  data nti0 / 0 /
  data xsml / 0.0E+00 / 
  data xmax / 0.0E+00 /

  if ( nti0 == 0 ) then
    nti0 = r4_inits ( bi0cs, 12, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) )
    xmax = log ( r4_mach ( 2 ) )
  end if

  y = abs ( x )

  if ( y <= xsml ) then
    r4_besi0 = 1.0E+00
  else if ( y <= 3.0E+00 ) then
    r4_besi0 = 2.75E+00 + r4_csevl ( y * y / 4.5E+00 - 1.0E+00, &
      bi0cs, nti0 )
  else if ( y <= xmax ) then
    r4_besi0 = exp ( y ) * r4_besi0e ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESI0 - Fatal error!'
    write ( *, '(a)' ) '  Result overflows.'
    write ( *, '(a)' ) '  |X| too large.'
    stop
  end if

  return
end
function r4_besi0e ( x )

!*****************************************************************************80
!
!! R4_BESI0E evaluates the exponentially scaled Bessel function I0(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESI0E, the exponentially scaled Bessel 
!    function I0(X).
!
  implicit none

  real ( kind = 4 ) ai02cs(22)
  real ( kind = 4 ) ai0cs(21)
  real ( kind = 4 ) bi0cs(12)
  integer ( kind = 4 ) ntai0
  integer ( kind = 4 ) ntai02
  integer ( kind = 4 ) nti0
  real ( kind = 4 ) r4_besi0e
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save ai02cs
  save ai0cs
  save bi0cs
  save ntai0
  save ntai02
  save nti0
  save xsml

  data bi0cs( 1) /   -0.07660547252839144951E+00 /
  data bi0cs( 2) /    1.927337953993808270E+00 /
  data bi0cs( 3) /    0.2282644586920301339E+00 /
  data bi0cs( 4) /    0.01304891466707290428E+00 /
  data bi0cs( 5) /    0.00043442709008164874E+00 /
  data bi0cs( 6) /    0.00000942265768600193E+00 /
  data bi0cs( 7) /    0.00000014340062895106E+00 /
  data bi0cs( 8) /    0.00000000161384906966E+00 /
  data bi0cs( 9) /    0.00000000001396650044E+00 /
  data bi0cs(10) /    0.00000000000009579451E+00 /
  data bi0cs(11) /    0.00000000000000053339E+00 /
  data bi0cs(12) /    0.00000000000000000245E+00 /

  data ai0cs( 1) /    0.07575994494023796E+00 /
  data ai0cs( 2) /    0.00759138081082334E+00 /
  data ai0cs( 3) /    0.00041531313389237E+00 /
  data ai0cs( 4) /    0.00001070076463439E+00 /
  data ai0cs( 5) /   -0.00000790117997921E+00 /
  data ai0cs( 6) /   -0.00000078261435014E+00 /
  data ai0cs( 7) /    0.00000027838499429E+00 /
  data ai0cs( 8) /    0.00000000825247260E+00 /
  data ai0cs( 9) /   -0.00000001204463945E+00 /
  data ai0cs(10) /    0.00000000155964859E+00 /
  data ai0cs(11) /    0.00000000022925563E+00 /
  data ai0cs(12) /   -0.00000000011916228E+00 /
  data ai0cs(13) /    0.00000000001757854E+00 /
  data ai0cs(14) /    0.00000000000112822E+00 /
  data ai0cs(15) /   -0.00000000000114684E+00 /
  data ai0cs(16) /    0.00000000000027155E+00 /
  data ai0cs(17) /   -0.00000000000002415E+00 /
  data ai0cs(18) /   -0.00000000000000608E+00 /
  data ai0cs(19) /    0.00000000000000314E+00 /
  data ai0cs(20) /   -0.00000000000000071E+00 /
  data ai0cs(21) /    0.00000000000000007E+00 /

  data ai02cs( 1) /    0.05449041101410882E+00 /
  data ai02cs( 2) /    0.00336911647825569E+00 /
  data ai02cs( 3) /    0.00006889758346918E+00 /
  data ai02cs( 4) /    0.00000289137052082E+00 /
  data ai02cs( 5) /    0.00000020489185893E+00 /
  data ai02cs( 6) /    0.00000002266668991E+00 /
  data ai02cs( 7) /    0.00000000339623203E+00 /
  data ai02cs( 8) /    0.00000000049406022E+00 /
  data ai02cs( 9) /    0.00000000001188914E+00 /
  data ai02cs(10) /   -0.00000000003149915E+00 /
  data ai02cs(11) /   -0.00000000001321580E+00 /
  data ai02cs(12) /   -0.00000000000179419E+00 /
  data ai02cs(13) /    0.00000000000071801E+00 /
  data ai02cs(14) /    0.00000000000038529E+00 /
  data ai02cs(15) /    0.00000000000001539E+00 /
  data ai02cs(16) /   -0.00000000000004151E+00 /
  data ai02cs(17) /   -0.00000000000000954E+00 /
  data ai02cs(18) /    0.00000000000000382E+00 /
  data ai02cs(19) /    0.00000000000000176E+00 /
  data ai02cs(20) /   -0.00000000000000034E+00 /
  data ai02cs(21) /   -0.00000000000000027E+00 /
  data ai02cs(22) /    0.00000000000000003E+00 /

  data ntai0 / 0 /
  data ntai02 / 0 /
  data nti0 / 0 /
  data xsml / 0.0E+00 /

  if ( nti0 == 0 ) then
    nti0 = r4_inits ( bi0cs, 12, 0.1E+00 * r4_mach ( 3 ) )
    ntai0 = r4_inits ( ai0cs, 21, 0.1E+00 * r4_mach ( 3 ) )
    ntai02 = r4_inits ( ai02cs, 22, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= xsml ) then
    r4_besi0e = 1.0E+00
  else if ( y <= 3.0E+00 ) then
    r4_besi0e = exp ( - y ) * ( 2.75E+00 + &
      r4_csevl ( y * y / 4.5E+00 - 1.0E+00, bi0cs, nti0 ) )
  else if ( y <= 8.0E+00 ) then
    r4_besi0e = ( 0.375E+00 + r4_csevl &
      ( ( 48.0E+00 / y - 11.0E+00 ) / 5.0E+00, ai0cs, ntai0 ) ) &
      / sqrt ( y )
  else
    r4_besi0e = ( 0.375E+00 + r4_csevl &
      ( 16.0E+00 / y - 1.0E+00, ai02cs, ntai02 ) ) / sqrt ( y )
  end if

  return
end
function r4_besi1 ( x )

!*****************************************************************************80
!
!! R4_BESI1 evaluates the Bessel function I of order 1 of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESI1, the Bessel function I of order 1 of X.
!
  implicit none

  real ( kind = 4 ) bi1cs(11)
  integer ( kind = 4 ) nti1
  real ( kind = 4 ) r4_besi1
  real ( kind = 4 ) r4_besi1e
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save bi1cs
  save nti1
  save xmax
  save xmin
  save xsml

  data bi1cs( 1) /   -0.001971713261099859E+00 /
  data bi1cs( 2) /    0.40734887667546481E+00 /
  data bi1cs( 3) /    0.034838994299959456E+00 /
  data bi1cs( 4) /    0.001545394556300123E+00 /
  data bi1cs( 5) /    0.000041888521098377E+00 /
  data bi1cs( 6) /    0.000000764902676483E+00 /
  data bi1cs( 7) /    0.000000010042493924E+00 /
  data bi1cs( 8) /    0.000000000099322077E+00 /
  data bi1cs( 9) /    0.000000000000766380E+00 /
  data bi1cs(10) /    0.000000000000004741E+00 /
  data bi1cs(11) /    0.000000000000000024E+00 /

  data nti1 / 0 /
  data xmax / 0.0E+00 /
  data xmin / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( nti1 == 0 ) then
    nti1 = r4_inits ( bi1cs, 11, 0.1E+00 * r4_mach ( 3 ) )
    xmin = 2.0E+00 * r4_mach ( 1 )
    xsml = sqrt ( 8.0E+00 * r4_mach ( 3 ) )
    xmax = log ( r4_mach ( 2 ) )
  end if

  y = abs ( x )

  if ( y <= xmin ) then
    r4_besi1 = 0.0E+00
  else if ( y <= xsml ) then
    r4_besi1 = 0.5E+00 * x
  else if ( y <= 3.0E+00 ) then
    r4_besi1 = x * ( 0.875E+00 + r4_csevl &
      ( y * y / 4.5E+00 - 1.0E+00, bi1cs, nti1 ) )
  else if ( y <= xmax ) then
    r4_besi1 = exp ( y ) * r4_besi1e ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESI1 - Fatal error!'
    write ( *, '(a)' ) '  Result overflows.'
    stop
  end if

  return
end
function r4_besi1e ( x )

!*****************************************************************************80
!
!! R4_BESI1E: exponentially scaled Bessel function I of order 1, R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESI1E, the exponentially scaled Bessel 
!    function I of order 1 of X.
!
  implicit none

  real ( kind = 4 ) ai1cs(21)
  real ( kind = 4 ) ai12cs(22)
  real ( kind = 4 ) bi1cs(11)
  integer ( kind = 4 ) ntai1
  integer ( kind = 4 ) ntai12
  integer ( kind = 4 ) nti1
  real ( kind = 4 ) r4_besi1e
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save ai1cs
  save ai12cs
  save bi1cs
  save ntai1
  save ntai12
  save nti1
  save xmin
  save xsml

  data bi1cs( 1) /   -0.001971713261099859E+00 /
  data bi1cs( 2) /    0.40734887667546481E+00 /
  data bi1cs( 3) /    0.034838994299959456E+00 /
  data bi1cs( 4) /    0.001545394556300123E+00 /
  data bi1cs( 5) /    0.000041888521098377E+00 /
  data bi1cs( 6) /    0.000000764902676483E+00 /
  data bi1cs( 7) /    0.000000010042493924E+00 /
  data bi1cs( 8) /    0.000000000099322077E+00 /
  data bi1cs( 9) /    0.000000000000766380E+00 /
  data bi1cs(10) /    0.000000000000004741E+00 /
  data bi1cs(11) /    0.000000000000000024E+00 /

  data ai1cs( 1) /   -0.02846744181881479E+00 /
  data ai1cs( 2) /   -0.01922953231443221E+00 /
  data ai1cs( 3) /   -0.00061151858579437E+00 /
  data ai1cs( 4) /   -0.00002069971253350E+00 /
  data ai1cs( 5) /    0.00000858561914581E+00 /
  data ai1cs( 6) /    0.00000104949824671E+00 /
  data ai1cs( 7) /   -0.00000029183389184E+00 /
  data ai1cs( 8) /   -0.00000001559378146E+00 /
  data ai1cs( 9) /    0.00000001318012367E+00 /
  data ai1cs(10) /   -0.00000000144842341E+00 /
  data ai1cs(11) /   -0.00000000029085122E+00 /
  data ai1cs(12) /    0.00000000012663889E+00 /
  data ai1cs(13) /   -0.00000000001664947E+00 /
  data ai1cs(14) /   -0.00000000000166665E+00 /
  data ai1cs(15) /    0.00000000000124260E+00 /
  data ai1cs(16) /   -0.00000000000027315E+00 /
  data ai1cs(17) /    0.00000000000002023E+00 /
  data ai1cs(18) /    0.00000000000000730E+00 /
  data ai1cs(19) /   -0.00000000000000333E+00 /
  data ai1cs(20) /    0.00000000000000071E+00 /
  data ai1cs(21) /   -0.00000000000000006E+00 /

  data ai12cs( 1) /    0.02857623501828014E+00 /
  data ai12cs( 2) /   -0.00976109749136147E+00 /
  data ai12cs( 3) /   -0.00011058893876263E+00 /
  data ai12cs( 4) /   -0.00000388256480887E+00 /
  data ai12cs( 5) /   -0.00000025122362377E+00 /
  data ai12cs( 6) /   -0.00000002631468847E+00 /
  data ai12cs( 7) /   -0.00000000383538039E+00 /
  data ai12cs( 8) /   -0.00000000055897433E+00 /
  data ai12cs( 9) /   -0.00000000001897495E+00 /
  data ai12cs(10) /    0.00000000003252602E+00 /
  data ai12cs(11) /    0.00000000001412580E+00 /
  data ai12cs(12) /    0.00000000000203564E+00 /
  data ai12cs(13) /   -0.00000000000071985E+00 /
  data ai12cs(14) /   -0.00000000000040836E+00 /
  data ai12cs(15) /   -0.00000000000002101E+00 /
  data ai12cs(16) /    0.00000000000004273E+00 /
  data ai12cs(17) /    0.00000000000001041E+00 /
  data ai12cs(18) /   -0.00000000000000382E+00 /
  data ai12cs(19) /   -0.00000000000000186E+00 /
  data ai12cs(20) /    0.00000000000000033E+00 /
  data ai12cs(21) /    0.00000000000000028E+00 /
  data ai12cs(22) /   -0.00000000000000003E+00 /

  data ntai1 / 0 /
  data ntai12 / 0 /
  data nti1 / 0 /
  data xmin / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( nti1 == 0 ) then
    nti1 = r4_inits ( bi1cs, 11, 0.1E+00 * r4_mach ( 3 ) )
    ntai1 = r4_inits ( ai1cs, 21, 0.1E+00 * r4_mach ( 3 ) )
    ntai12 = r4_inits ( ai12cs, 22, 0.1E+00 * r4_mach ( 3 ) )
    xmin = 2.0E+00 * r4_mach ( 1 )
    xsml = sqrt ( 8.0E+00 * r4_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( x == 0.0E+00 ) then
    r4_besi1e = 0.0E+00
  else if ( y <= xmin ) then
    r4_besi1e = 0.0E+00
  else if ( y <= xsml ) then
    r4_besi1e = 0.5E+00 * x
    r4_besi1e = exp ( - y ) * r4_besi1e
  else if ( y <= 3.0E+00 ) then
    r4_besi1e = x * ( 0.875E+00 &
      + r4_csevl ( y * y / 4.5E+00 - 1.0E+00, bi1cs, nti1 ) )
    r4_besi1e = exp ( - y ) * r4_besi1e
  else if ( y <= 8.0E+00 ) then
    r4_besi1e = ( 0.375E+00 &
      + r4_csevl ( ( 48.0E+00 / y - 11.0E+00 ) / 5.0E+00, &
      ai1cs, ntai1) ) / sqrt ( y )
    if ( x < 0.0E+00 ) then
      r4_besi1e = - r4_besi1e
    end if
  else
    r4_besi1e = ( 0.375E+00 &
      + r4_csevl ( 16.0E+00 / y - 1.0E+00, ai12cs, ntai12 ) ) &
      / sqrt ( y )
    if ( x < 0.0E+00 ) then
      r4_besi1e = - r4_besi1e
    end if
  end if

  return
end
function r4_besj0 ( x )

!*****************************************************************************80
!
!! R4_BESJ0 evaluates the Bessel function J of order 0 of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESJ0, the Bessel function J of order 0 of X.
!
  implicit none

  real ( kind = 4 ) ampl
  real ( kind = 4 ) bj0cs(13)
  real ( kind = 4 ) bm0cs(21)
  real ( kind = 4 ) bth0cs(24)
  integer ( kind = 4 ) ntj0
  integer ( kind = 4 ) ntm0
  integer ( kind = 4 ) ntth0
  real ( kind = 4 ) pi4
  real ( kind = 4 ) r4_besj0
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) theta
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  save bj0cs
  save bm0cs
  save bth0cs
  save ntj0
  save ntm0
  save ntth0
  save pi4
  save xmax
  save xsml

  data bj0cs( 1) /    0.100254161968939137E+00 /
  data bj0cs( 2) /   -0.665223007764405132E+00 /
  data bj0cs( 3) /    0.248983703498281314E+00 /
  data bj0cs( 4) /   -0.0332527231700357697E+00 /
  data bj0cs( 5) /    0.0023114179304694015E+00 /
  data bj0cs( 6) /   -0.0000991127741995080E+00 /
  data bj0cs( 7) /    0.0000028916708643998E+00 /
  data bj0cs( 8) /   -0.0000000612108586630E+00 /
  data bj0cs( 9) /    0.0000000009838650793E+00 /
  data bj0cs(10) /   -0.0000000000124235515E+00 /
  data bj0cs(11) /    0.0000000000001265433E+00 /
  data bj0cs(12) /   -0.0000000000000010619E+00 /
  data bj0cs(13) /    0.0000000000000000074E+00 /

  data bm0cs( 1) /    0.09284961637381644E+00 /
  data bm0cs( 2) /   -0.00142987707403484E+00 /
  data bm0cs( 3) /    0.00002830579271257E+00 /
  data bm0cs( 4) /   -0.00000143300611424E+00 /
  data bm0cs( 5) /    0.00000012028628046E+00 /
  data bm0cs( 6) /   -0.00000001397113013E+00 /
  data bm0cs( 7) /    0.00000000204076188E+00 /
  data bm0cs( 8) /   -0.00000000035399669E+00 /
  data bm0cs( 9) /    0.00000000007024759E+00 /
  data bm0cs(10) /   -0.00000000001554107E+00 /
  data bm0cs(11) /    0.00000000000376226E+00 /
  data bm0cs(12) /   -0.00000000000098282E+00 /
  data bm0cs(13) /    0.00000000000027408E+00 /
  data bm0cs(14) /   -0.00000000000008091E+00 /
  data bm0cs(15) /    0.00000000000002511E+00 /
  data bm0cs(16) /   -0.00000000000000814E+00 /
  data bm0cs(17) /    0.00000000000000275E+00 /
  data bm0cs(18) /   -0.00000000000000096E+00 /
  data bm0cs(19) /    0.00000000000000034E+00 /
  data bm0cs(20) /   -0.00000000000000012E+00 /
  data bm0cs(21) /    0.00000000000000004E+00 /

  data bth0cs( 1) /   -0.24639163774300119E+00 /
  data bth0cs( 2) /    0.001737098307508963E+00 /
  data bth0cs( 3) /   -0.000062183633402968E+00 /
  data bth0cs( 4) /    0.000004368050165742E+00 /
  data bth0cs( 5) /   -0.000000456093019869E+00 /
  data bth0cs( 6) /    0.000000062197400101E+00 /
  data bth0cs( 7) /   -0.000000010300442889E+00 /
  data bth0cs( 8) /    0.000000001979526776E+00 /
  data bth0cs( 9) /   -0.000000000428198396E+00 /
  data bth0cs(10) /    0.000000000102035840E+00 /
  data bth0cs(11) /   -0.000000000026363898E+00 /
  data bth0cs(12) /    0.000000000007297935E+00 /
  data bth0cs(13) /   -0.000000000002144188E+00 /
  data bth0cs(14) /    0.000000000000663693E+00 /
  data bth0cs(15) /   -0.000000000000215126E+00 /
  data bth0cs(16) /    0.000000000000072659E+00 /
  data bth0cs(17) /   -0.000000000000025465E+00 /
  data bth0cs(18) /    0.000000000000009229E+00 /
  data bth0cs(19) /   -0.000000000000003448E+00 /
  data bth0cs(20) /    0.000000000000001325E+00 /
  data bth0cs(21) /   -0.000000000000000522E+00 /
  data bth0cs(22) /    0.000000000000000210E+00 /
  data bth0cs(23) /   -0.000000000000000087E+00 /
  data bth0cs(24) /    0.000000000000000036E+00 /

  data pi4 / 0.78539816339744831E+00 /
  data ntj0 / 0 /
  data ntm0 / 0 /
  data ntth0 / 0 / 
  data xsml / 0.0E+00 / 
  data xmax / 0.0E+00 /

  if ( ntj0 == 0 ) then
    ntj0 = r4_inits ( bj0cs, 13, 0.1E+00 * r4_mach ( 3 ) )
    ntm0 = r4_inits ( bm0cs, 21, 0.1E+00 * r4_mach ( 3 ) )
    ntth0 = r4_inits ( bth0cs, 24, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) )
    xmax = 1.0E+00 / r4_mach ( 4 )
  end if

  y = abs ( x )
 
  if ( y <= xsml ) then
    r4_besj0 = 1.0E+00
  else if ( y <= 4.0E+00 ) then
    r4_besj0 = r4_csevl ( 0.125E+00 * y * y - 1.0E+00, bj0cs, ntj0 )
  else
    z = 32.0E+00 / y / y - 1.0E+00
    ampl = ( 0.75E+00 + r4_csevl ( z, bm0cs, ntm0 ) ) / sqrt ( y )
    theta = y - pi4 + r4_csevl ( z, bth0cs, ntth0 ) / y
    r4_besj0 = ampl * cos ( theta )
  end if

  return
end
function r4_besj1 ( x )

!*****************************************************************************80
!
!! R4_BESJ1 evaluates the Bessel function J of order 1 of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESJ1, the Bessel function J of order 1 of X.
!
  implicit none

  real ( kind = 4 ) ampl
  real ( kind = 4 ) bj1cs(12)
  real ( kind = 4 ) bm1cs(21)
  real ( kind = 4 ) bth1cs(24)
  integer ( kind = 4 ) ntj1
  integer ( kind = 4 ) ntm1
  integer ( kind = 4 ) ntth1
  real ( kind = 4 ) pi4
  real ( kind = 4 ) r4_besj1
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) theta
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  save bj1cs
  save bm1cs
  save bth1cs
  save ntj1
  save ntm1
  save ntth1
  save pi4
  save xmax
  save xmin
  save xsml

  data bj1cs( 1) /   -0.11726141513332787E+00 /
  data bj1cs( 2) /   -0.25361521830790640E+00 /
  data bj1cs( 3) /   +0.050127080984469569E+00 /
  data bj1cs( 4) /   -0.004631514809625081E+00 /
  data bj1cs( 5) /   +0.000247996229415914E+00 /
  data bj1cs( 6) /   -0.000008678948686278E+00 /
  data bj1cs( 7) /   +0.000000214293917143E+00 /
  data bj1cs( 8) /   -0.000000003936093079E+00 /
  data bj1cs( 9) /   +0.000000000055911823E+00 /
  data bj1cs(10) /   -0.000000000000632761E+00 /
  data bj1cs(11) /   +0.000000000000005840E+00 /
  data bj1cs(12) /   -0.000000000000000044E+00 /

  data bm1cs( 1) /   +0.1047362510931285E+00 /
  data bm1cs( 2) /   +0.00442443893702345E+00 /
  data bm1cs( 3) /   -0.00005661639504035E+00 /
  data bm1cs( 4) /   +0.00000231349417339E+00 /
  data bm1cs( 5) /   -0.00000017377182007E+00 /
  data bm1cs( 6) /   +0.00000001893209930E+00 /
  data bm1cs( 7) /   -0.00000000265416023E+00 /
  data bm1cs( 8) /   +0.00000000044740209E+00 /
  data bm1cs( 9) /   -0.00000000008691795E+00 /
  data bm1cs(10) /   +0.00000000001891492E+00 /
  data bm1cs(11) /   -0.00000000000451884E+00 /
  data bm1cs(12) /   +0.00000000000116765E+00 /
  data bm1cs(13) /   -0.00000000000032265E+00 /
  data bm1cs(14) /   +0.00000000000009450E+00 /
  data bm1cs(15) /   -0.00000000000002913E+00 /
  data bm1cs(16) /   +0.00000000000000939E+00 /
  data bm1cs(17) /   -0.00000000000000315E+00 /
  data bm1cs(18) /   +0.00000000000000109E+00 /
  data bm1cs(19) /   -0.00000000000000039E+00 /
  data bm1cs(20) /   +0.00000000000000014E+00 /
  data bm1cs(21) /   -0.00000000000000005E+00 /

  data bth1cs( 1) /   +0.74060141026313850E+00 /
  data bth1cs( 2) /   -0.004571755659637690E+00 /
  data bth1cs( 3) /   +0.000119818510964326E+00 /
  data bth1cs( 4) /   -0.000006964561891648E+00 /
  data bth1cs( 5) /   +0.000000655495621447E+00 /
  data bth1cs( 6) /   -0.000000084066228945E+00 /
  data bth1cs( 7) /   +0.000000013376886564E+00 /
  data bth1cs( 8) /   -0.000000002499565654E+00 /
  data bth1cs( 9) /   +0.000000000529495100E+00 /
  data bth1cs(10) /   -0.000000000124135944E+00 /
  data bth1cs(11) /   +0.000000000031656485E+00 /
  data bth1cs(12) /   -0.000000000008668640E+00 /
  data bth1cs(13) /   +0.000000000002523758E+00 /
  data bth1cs(14) /   -0.000000000000775085E+00 /
  data bth1cs(15) /   +0.000000000000249527E+00 /
  data bth1cs(16) /   -0.000000000000083773E+00 /
  data bth1cs(17) /   +0.000000000000029205E+00 /
  data bth1cs(18) /   -0.000000000000010534E+00 /
  data bth1cs(19) /   +0.000000000000003919E+00 /
  data bth1cs(20) /   -0.000000000000001500E+00 /
  data bth1cs(21) /   +0.000000000000000589E+00 /
  data bth1cs(22) /   -0.000000000000000237E+00 /
  data bth1cs(23) /   +0.000000000000000097E+00 /
  data bth1cs(24) /   -0.000000000000000040E+00 /

  data pi4 / 0.78539816339744831E+00 /
  data ntj1 / 0 /
  data ntm1 / 0 /
  data ntth1 / 0 /
  data xsml / 0.0E+00 /
  data xmin / 0.0E+00 /
  data xmax / 0.0E+00 /

  if ( ntj1 == 0 ) then
    ntj1 = r4_inits ( bj1cs, 12, 0.1E+00 * r4_mach ( 3 ) )
    ntm1 = r4_inits ( bm1cs, 21, 0.1E+00 * r4_mach ( 3 ) )
    ntth1 = r4_inits ( bth1cs, 24, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( 8.0E+00 * r4_mach ( 3 ) )
    xmin = 2.0E+00 * r4_mach ( 1 )
    xmax = 1.0E+00 / r4_mach ( 4 )
  end if

  y = abs ( x )

  if ( y <= xmin ) then
    r4_besj1 = 0.0E+00
  else if ( y <= xsml ) then
    r4_besj1 = 0.5E+00 * x
  else if ( y <= 4.0E+00 ) then
    r4_besj1 = x * ( 0.25E+00 &
      + r4_csevl ( 0.125E+00 * y * y - 1.0E+00, bj1cs, ntj1 ) )
  else
    z = 32.0E+00 / y / y - 1.0E+00
    ampl = ( 0.75E+00 + r4_csevl ( z, bm1cs, ntm1 ) ) / sqrt ( y )
    theta = y - 3.0E+00 * pi4 &
      + r4_csevl ( z, bth1cs, ntth1 ) / y
    if ( x < 0.0E+00 ) then
      r4_besj1 = - ampl * cos ( theta )
    else
      r4_besj1 = + ampl * cos ( theta )
    end if
  end if

  return
end
function r4_besk ( nu, x )

!*****************************************************************************80
!
!! R4_BESK evaluates the Bessel function K of order NU of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) NU, the order.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESK, the Bessel function K of order NU at X.
!
  implicit none

  real ( kind = 4 ), allocatable :: bke(:)
  integer ( kind = 4 ) nin
  real ( kind = 4 ) nu
  real ( kind = 4 ) r4_besk
  real ( kind = 4 ) x
  real ( kind = 4 ) xnu

  xnu = nu - int ( nu )
  nin = int ( nu ) + 1
  allocate ( bke(1:nin) )

  call r4_besks ( xnu, x, nin, bke )

  r4_besk = bke(nin)

  deallocate ( bke )

  return
end
function r4_besk0 ( x )

!*****************************************************************************80
!
!! R4_BESK0 evaluates the Bessel function K of order 0 of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESK0, the Bessel function K of order 0 of X.
!
  implicit none

  real ( kind = 4 ) bk0cs(11)
  integer ( kind = 4 ) ntk0
  real ( kind = 4 ) r4_besi0
  real ( kind = 4 ) r4_besk0
  real ( kind = 4 ) r4_besk0e
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save bk0cs
  save ntk0
  save xmax
  save xsml

  data bk0cs( 1) /   -0.03532739323390276872E+00 /
  data bk0cs( 2) /    0.3442898999246284869E+00 /
  data bk0cs( 3) /    0.03597993651536150163E+00 /
  data bk0cs( 4) /    0.00126461541144692592E+00 /
  data bk0cs( 5) /    0.00002286212103119451E+00 /
  data bk0cs( 6) /    0.00000025347910790261E+00 /
  data bk0cs( 7) /    0.00000000190451637722E+00 /
  data bk0cs( 8) /    0.00000000001034969525E+00 /
  data bk0cs( 9) /    0.00000000000004259816E+00 /
  data bk0cs(10) /    0.00000000000000013744E+00 /
  data bk0cs(11) /    0.00000000000000000035E+00 /

  data ntk0 / 0 /
  data xsml / 0.0E+00 /
  data xmax / 0.0E+00 /

  if ( ntk0 == 0 ) then
    ntk0 = r4_inits ( bk0cs, 11, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) )
    xmax = - log ( r4_mach ( 1 ) )
    xmax = xmax - 0.5E+00 * xmax * log ( xmax ) &
      / ( xmax + 0.5E+00 ) - 0.01E+00
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESK0 = Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0E+00
    r4_besk0 = - log ( 0.5E+00 * x ) * r4_besi0 ( x ) &
      - 0.25E+00 + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk0cs, ntk0 )
  else if ( x <= 2.0E+00 ) then
    y = x * x
    r4_besk0 = - log ( 0.5E+00 * x ) * r4_besi0 ( x ) &
      - 0.25E+00 + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk0cs, ntk0 )
  else if ( x <= xmax ) then
    r4_besk0 = exp ( - x ) * r4_besk0e ( x )
  else
    r4_besk0 = 0.0E+00
  end if

  return
end
function r4_besk0e ( x )

!*****************************************************************************80
!
!! R4_BESK0E evaluates the exponentially scaled Bessel function K0(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESK0E, the exponentially scaled Bessel 
!    function K0(X).
!
  implicit none

  real ( kind = 4 ) ak02cs(14)
  real ( kind = 4 ) ak0cs(17)
  real ( kind = 4 ) bk0cs(11)
  integer ( kind = 4 ) ntak0
  integer ( kind = 4 ) ntak02
  integer ( kind = 4 ) ntk0
  real ( kind = 4 ) r4_besi0
  real ( kind = 4 ) r4_besk0e
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save ak02cs
  save ak0cs
  save bk0cs
  save ntak0
  save ntak02
  save ntk0
  save xsml

  data bk0cs( 1) /   -0.03532739323390276872E+00 /
  data bk0cs( 2) /   +0.3442898999246284869E+00 /
  data bk0cs( 3) /   +0.03597993651536150163E+00 /
  data bk0cs( 4) /   +0.00126461541144692592E+00 /
  data bk0cs( 5) /   +0.00002286212103119451E+00 /
  data bk0cs( 6) /   +0.00000025347910790261E+00 /
  data bk0cs( 7) /   +0.00000000190451637722E+00 /
  data bk0cs( 8) /   +0.00000000001034969525E+00 /
  data bk0cs( 9) /   +0.00000000000004259816E+00 /
  data bk0cs(10) /   +0.00000000000000013744E+00 /
  data bk0cs(11) /   +0.00000000000000000035E+00 /

  data ak0cs( 1) /   -0.07643947903327941E+00 /
  data ak0cs( 2) /   -0.02235652605699819E+00 /
  data ak0cs( 3) /   +0.00077341811546938E+00 /
  data ak0cs( 4) /   -0.00004281006688886E+00 /
  data ak0cs( 5) /   +0.00000308170017386E+00 /
  data ak0cs( 6) /   -0.00000026393672220E+00 /
  data ak0cs( 7) /   +0.00000002563713036E+00 /
  data ak0cs( 8) /   -0.00000000274270554E+00 /
  data ak0cs( 9) /   +0.00000000031694296E+00 /
  data ak0cs(10) /   -0.00000000003902353E+00 /
  data ak0cs(11) /   +0.00000000000506804E+00 /
  data ak0cs(12) /   -0.00000000000068895E+00 /
  data ak0cs(13) /   +0.00000000000009744E+00 /
  data ak0cs(14) /   -0.00000000000001427E+00 /
  data ak0cs(15) /   +0.00000000000000215E+00 /
  data ak0cs(16) /   -0.00000000000000033E+00 /
  data ak0cs(17) /   +0.00000000000000005E+00 /

  data ak02cs( 1) /   -0.01201869826307592E+00 /
  data ak02cs( 2) /   -0.00917485269102569E+00 /
  data ak02cs( 3) /   +0.00014445509317750E+00 /
  data ak02cs( 4) /   -0.00000401361417543E+00 /
  data ak02cs( 5) /   +0.00000015678318108E+00 /
  data ak02cs( 6) /   -0.00000000777011043E+00 /
  data ak02cs( 7) /   +0.00000000046111825E+00 /
  data ak02cs( 8) /   -0.00000000003158592E+00 /
  data ak02cs( 9) /   +0.00000000000243501E+00 /
  data ak02cs(10) /   -0.00000000000020743E+00 /
  data ak02cs(11) /   +0.00000000000001925E+00 /
  data ak02cs(12) /   -0.00000000000000192E+00 /
  data ak02cs(13) /   +0.00000000000000020E+00 /
  data ak02cs(14) /   -0.00000000000000002E+00 /

  data ntak0 / 0 /
  data ntak02 / 0 /
  data ntk0 / 0 /
  data xsml / 0.0E+00 /

  if ( ntk0 == 0 ) then
    ntk0 = r4_inits ( bk0cs, 11, 0.1E+00 * r4_mach ( 3 ) )
    ntak0 = r4_inits ( ak0cs, 17, 0.1E+00 * r4_mach ( 3 ) )
    ntak02 = r4_inits ( ak02cs, 14, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) )
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESK0E - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0E+00
    r4_besk0e = exp ( x ) * ( - log ( 0.5E+00 * x ) &
      * r4_besi0 ( x ) - 0.25E+00 &
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk0cs, ntk0 ) )
  else if ( x <= 2.0E+00 ) then
    y = x * x
    r4_besk0e = exp ( x ) * ( - log ( 0.5E+00 * x ) &
      * r4_besi0 ( x ) - 0.25E+00 &
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk0cs, ntk0 ) )
  else if ( x <= 8.0E+00 ) then
    r4_besk0e = ( 1.25E+00 &
      + r4_csevl ( ( 16.0E+00 / x - 5.0E+00 ) / 3.0E+00, &
      ak0cs, ntak0 ) ) / sqrt ( x )
  else
    r4_besk0e = ( 1.25E+00 &
      + r4_csevl ( 16.0E+00 / x - 1.0E+00, ak02cs, ntak02 ) ) &
      / sqrt ( x )
  end if

  return
end
function r4_besk1 ( x )

!*****************************************************************************80
!
!! R4_BESK1 evaluates the Bessel function K of order 1 of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESK1, the Bessel function K of order 1 of X.
!
  implicit none

  real ( kind = 4 ) bk1cs(11)
  integer ( kind = 4 ) ntk1
  real ( kind = 4 ) r4_besi1
  real ( kind = 4 ) r4_besk1
  real ( kind = 4 ) r4_besk1e
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save bk1cs
  save ntk1
  save xmax
  save xmin
  save xsml

  data bk1cs( 1) /    0.0253002273389477705E+00 /
  data bk1cs( 2) /   -0.353155960776544876E+00 /
  data bk1cs( 3) /   -0.122611180822657148E+00 /
  data bk1cs( 4) /   -0.0069757238596398643E+00 /
  data bk1cs( 5) /   -0.0001730288957513052E+00 /
  data bk1cs( 6) /   -0.0000024334061415659E+00 /
  data bk1cs( 7) /   -0.0000000221338763073E+00 /
  data bk1cs( 8) /   -0.0000000001411488392E+00 /
  data bk1cs( 9) /   -0.0000000000006666901E+00 /
  data bk1cs(10) /   -0.0000000000000024274E+00 /
  data bk1cs(11) /   -0.0000000000000000070E+00 /

  data ntk1 / 0 /
  data xmax / 0.0E+00 /
  data xmin / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( ntk1 == 0 ) then
    ntk1 = r4_inits ( bk1cs, 11, 0.1E+00 * r4_mach ( 3 ) )
    xmin = exp ( max ( log ( r4_mach ( 1 ) ), &
      - log ( r4_mach ( 2 ) ) ) + 0.01E+00 )
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) )
    xmax = - log ( r4_mach ( 1 ) )
    xmax = xmax - 0.5E+00 * xmax * log ( xmax ) &
      / ( xmax + 0.5E+00 ) - 0.01E+00
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESK1 = Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0E+00
    r4_besk1 = log ( 0.5E+00 * x ) * r4_besi1 ( x ) + ( 0.75E+00 &
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk1cs, ntk1 ) ) / x
  else if ( x <= 2.0E+00 ) then
    y = x * x
    r4_besk1 = log ( 0.5E+00 * x ) * r4_besi1 ( x ) + ( 0.75E+00 &
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk1cs, ntk1 ) ) / x
  else if ( x <= xmax ) then
    r4_besk1 = exp ( - x ) * r4_besk1e ( x )
  else
    r4_besk1 = 0.0E+00
  end if

  return
end
function r4_besk1e ( x )

!*****************************************************************************80
!
!! R4_BESK1E evaluates the exponentially scaled Bessel function K1(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESK1E, the exponentially scaled Bessel 
!    function K1(X).
!
  implicit none

  real ( kind = 4 ) ak12cs(14)
  real ( kind = 4 ) ak1cs(17)
  real ( kind = 4 ) bk1cs(11)
  integer ( kind = 4 ) ntak1
  integer ( kind = 4 ) ntak12
  integer ( kind = 4 ) ntk1
  real ( kind = 4 ) r4_besi1
  real ( kind = 4 ) r4_besk1e
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save ak12cs
  save ak1cs
  save bk1cs
  save ntak1
  save ntak12
  save ntk1
  save xmin
  save xsml

  data bk1cs( 1) /   +0.0253002273389477705E+00 /
  data bk1cs( 2) /   -0.353155960776544876E+00 /
  data bk1cs( 3) /   -0.122611180822657148E+00 /
  data bk1cs( 4) /   -0.0069757238596398643E+00 /
  data bk1cs( 5) /   -0.0001730288957513052E+00 /
  data bk1cs( 6) /   -0.0000024334061415659E+00 /
  data bk1cs( 7) /   -0.0000000221338763073E+00 /
  data bk1cs( 8) /   -0.0000000001411488392E+00 /
  data bk1cs( 9) /   -0.0000000000006666901E+00 /
  data bk1cs(10) /   -0.0000000000000024274E+00 /
  data bk1cs(11) /   -0.0000000000000000070E+00 /

  data ak1cs( 1) /   +0.2744313406973883E+00 /
  data ak1cs( 2) /   +0.07571989953199368E+00 /
  data ak1cs( 3) /   -0.00144105155647540E+00 /
  data ak1cs( 4) /   +0.00006650116955125E+00 /
  data ak1cs( 5) /   -0.00000436998470952E+00 /
  data ak1cs( 6) /   +0.00000035402774997E+00 /
  data ak1cs( 7) /   -0.00000003311163779E+00 /
  data ak1cs( 8) /   +0.00000000344597758E+00 /
  data ak1cs( 9) /   -0.00000000038989323E+00 /
  data ak1cs(10) /   +0.00000000004720819E+00 /
  data ak1cs(11) /   -0.00000000000604783E+00 /
  data ak1cs(12) /   +0.00000000000081284E+00 /
  data ak1cs(13) /   -0.00000000000011386E+00 /
  data ak1cs(14) /   +0.00000000000001654E+00 /
  data ak1cs(15) /   -0.00000000000000248E+00 /
  data ak1cs(16) /   +0.00000000000000038E+00 /
  data ak1cs(17) /   -0.00000000000000006E+00 /

  data ak12cs( 1) /   +0.06379308343739001E+00 /
  data ak12cs( 2) /   +0.02832887813049721E+00 /
  data ak12cs( 3) /   -0.00024753706739052E+00 /
  data ak12cs( 4) /   +0.00000577197245160E+00 /
  data ak12cs( 5) /   -0.00000020689392195E+00 /
  data ak12cs( 6) /   +0.00000000973998344E+00 /
  data ak12cs( 7) /   -0.00000000055853361E+00 /
  data ak12cs( 8) /   +0.00000000003732996E+00 /
  data ak12cs( 9) /   -0.00000000000282505E+00 /
  data ak12cs(10) /   +0.00000000000023720E+00 /
  data ak12cs(11) /   -0.00000000000002176E+00 /
  data ak12cs(12) /   +0.00000000000000215E+00 /
  data ak12cs(13) /   -0.00000000000000022E+00 /
  data ak12cs(14) /   +0.00000000000000002E+00 /

  data ntak1 / 0 /
  data ntak12 / 0 /
  data ntk1 / 0 /
  data xmin / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( ntk1 == 0 ) then
    ntk1 = r4_inits ( bk1cs, 11, 0.1E+00 * r4_mach ( 3 ) )
    ntak1 = r4_inits ( ak1cs, 17, 0.1E+00 * r4_mach ( 3 ) )
    ntak12 = r4_inits ( ak12cs, 14, 0.1E+00 * r4_mach ( 3 ) )
    xmin = exp ( max ( log ( r4_mach ( 1 ) ), &
      - log ( r4_mach ( 2 ) ) ) + 0.01E+00 )
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) )
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESK1E = Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0E+00
    r4_besk1e = exp ( x ) * ( log ( 0.5E+00 * x ) * r4_besi1 ( x ) &
      + ( 0.75E+00 &
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk1cs, ntk1 ) ) / x )
  else if ( x <= 2.0E+00 ) then
    y = x * x
    r4_besk1e = exp ( x ) * ( log ( 0.5E+00 * x ) * r4_besi1 ( x ) &
      + ( 0.75E+00 &
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk1cs, ntk1 ) ) / x )
  else if ( x <= 8.0E+00 ) then
    r4_besk1e = ( 1.25E+00 &
      + r4_csevl ( ( 16.0E+00 / x - 5.0E+00 ) / 3.0E+00, &
      ak1cs, ntak1 ) ) / sqrt ( x )
  else
    r4_besk1e = ( 1.25E+00 &
      + r4_csevl ( 16.E+00 / x - 1.0E+00, ak12cs, ntak12 ) ) &
      / sqrt ( x )
  end if

  return
end
subroutine r4_beskes ( xnu, x, nin, bke )

!*****************************************************************************80
!
!! R4_BESKES: a sequence of exponentially scaled K Bessel functions at X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) XNU, the order of the first function.
!    |XNU| < 1.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Input, integer ( kind = 4 ) NIN, the absolute value of NIN indicates the 
!    number of terms to compute.
!    If NIN < 0, successive values of NU count DOWN from XNU.
!    If NIN > 0, successive values of NU count UP from XNU.
!
!    Output, real ( kind = 4 ) BKE(abs(NIN)), the exponentially scaled 
!    K Bessel functions.
!
  implicit none

  integer ( kind = 4 ) nin

  real ( kind = 4 ) bke(abs(nin))
  real ( kind = 4 ) bknu1
  real ( kind = 4 ) direct
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iswtch
  integer ( kind = 4 ) n
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) v
  real ( kind = 4 ) vend
  real ( kind = 4 ) vincr
  real ( kind = 4 ) x
  real ( kind = 4 ) xnu

  v = abs ( xnu )
  n = abs ( nin )

  if ( 1.0E+00 <= v ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESKES - Fatal error!'
    write ( *, '(a)' ) '  |XNU| must be less than 1.'
    stop
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESKES - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  if ( n == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESKES - Fatal error!'
    write ( *, '(a)' ) '  N = 0.'
    stop
  end if

  call r4_knus ( v, x, bke(1), bknu1, iswtch )

  if ( n == 1 ) then
    return
  end if

  if ( nin < 0 ) then
    vincr = - 1.0E+00
  else
    vincr = + 1.0E+00
  end if

  if ( xnu < 0.0E+00 ) then
    direct = - vincr
  else
    direct = vincr
  end if

  bke(2) = bknu1

  if ( direct < 0.0E+00 ) then
    call r4_knus ( abs ( xnu + vincr ), x, bke(2), bknu1, iswtch )
  end if

  if ( n == 2 ) then
    return
  end if

  vend = abs ( xnu + real ( nin, kind = 4 ) ) - 1.0E+00

  v = xnu
  do i = 3, n
    v = v + vincr
    bke(i) = 2.0E+00 * v * bke(i-1) / x + bke(i-2)
  end do

  return
end
subroutine r4_besks ( xnu, x, nin, bk )

!*****************************************************************************80
!
!! R4_BESKS evaluates a sequence of K Bessel functions at X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) XNU, the order of the first function.
!    |XNU| < 1.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Input, integer ( kind = 4 ) NIN, the absolute value of NIN indicates 
!    the number of terms to compute.
!    If NIN < 0, successive values of NU count DOWN from XNU.
!    If NIN > 0, successive values of NU count UP from XNU.
!
!    Output, real ( kind = 4 ) BK(abs(NIN)), the K Bessel functions.
!
  implicit none

  integer ( kind = 4 ) nin

  real ( kind = 4 ) bk(abs(nin))
  real ( kind = 4 ) expxi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xnu

  save xmax

  data xmax / 0.0E+00 /

  if ( xmax == 0.0E+00 ) then
    xmax = - log ( r4_mach ( 1 ) )
    xmax = xmax + 0.5E+00 * log ( 3.14E+00 * 0.5E+00 / xmax )
  end if

  call r4_beskes ( xnu, x, nin, bk )

  expxi = exp ( - x )
  n = abs ( nin )

  do i = 1, n
    bk(i) = expxi * bk(i)
  end do

  return
end
function r4_besy0 ( x )

!*****************************************************************************80
!
!! R4_BESY0 evaluates the Bessel function Y of order 0 of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESY0, the Bessel function Y of order 0 of X.
!
  implicit none

  real ( kind = 4 ) alnhaf
  real ( kind = 4 ) ampl
  real ( kind = 4 ) bm0cs(21)
  real ( kind = 4 ) bth0cs(24)
  real ( kind = 4 ) by0cs(13)
  integer ( kind = 4 ) ntm0
  integer ( kind = 4 ) ntth0
  integer ( kind = 4 ) nty0
  real ( kind = 4 ) pi4
  real ( kind = 4 ) r4_besj0
  real ( kind = 4 ) r4_besy0
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) theta
  real ( kind = 4 ) twodpi
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  save alnhaf
  save bm0cs
  save bth0cs
  save by0cs
  save ntm0
  save ntth0
  save nty0
  save pi4
  save twodpi
  save xmax
  save xsml

  data by0cs( 1) /   -0.011277839392865573E+00 /
  data by0cs( 2) /   -0.12834523756042035E+00 /
  data by0cs( 3) /   -0.10437884799794249E+00 /
  data by0cs( 4) /   +0.023662749183969695E+00 /
  data by0cs( 5) /   -0.002090391647700486E+00 /
  data by0cs( 6) /   +0.000103975453939057E+00 /
  data by0cs( 7) /   -0.000003369747162423E+00 /
  data by0cs( 8) /   +0.000000077293842676E+00 /
  data by0cs( 9) /   -0.000000001324976772E+00 /
  data by0cs(10) /   +0.000000000017648232E+00 /
  data by0cs(11) /   -0.000000000000188105E+00 /
  data by0cs(12) /   +0.000000000000001641E+00 /
  data by0cs(13) /   -0.000000000000000011E+00 /

  data bm0cs( 1) /   +0.09284961637381644E+00 /
  data bm0cs( 2) /   -0.00142987707403484E+00 /
  data bm0cs( 3) /   +0.00002830579271257E+00 /
  data bm0cs( 4) /   -0.00000143300611424E+00 /
  data bm0cs( 5) /   +0.00000012028628046E+00 /
  data bm0cs( 6) /   -0.00000001397113013E+00 /
  data bm0cs( 7) /   +0.00000000204076188E+00 /
  data bm0cs( 8) /   -0.00000000035399669E+00 /
  data bm0cs( 9) /   +0.00000000007024759E+00 /
  data bm0cs(10) /   -0.00000000001554107E+00 /
  data bm0cs(11) /   +0.00000000000376226E+00 /
  data bm0cs(12) /   -0.00000000000098282E+00 /
  data bm0cs(13) /   +0.00000000000027408E+00 /
  data bm0cs(14) /   -0.00000000000008091E+00 /
  data bm0cs(15) /   +0.00000000000002511E+00 /
  data bm0cs(16) /   -0.00000000000000814E+00 /
  data bm0cs(17) /   +0.00000000000000275E+00 /
  data bm0cs(18) /   -0.00000000000000096E+00 /
  data bm0cs(19) /   +0.00000000000000034E+00 /
  data bm0cs(20) /   -0.00000000000000012E+00 /
  data bm0cs(21) /   +0.00000000000000004E+00 /

  data bth0cs( 1) /   -0.24639163774300119E+00 /
  data bth0cs( 2) /   +0.001737098307508963E+00 /
  data bth0cs( 3) /   -0.000062183633402968E+00 /
  data bth0cs( 4) /   +0.000004368050165742E+00 /
  data bth0cs( 5) /   -0.000000456093019869E+00 /
  data bth0cs( 6) /   +0.000000062197400101E+00 /
  data bth0cs( 7) /   -0.000000010300442889E+00 /
  data bth0cs( 8) /   +0.000000001979526776E+00 /
  data bth0cs( 9) /   -0.000000000428198396E+00 /
  data bth0cs(10) /   +0.000000000102035840E+00 /
  data bth0cs(11) /   -0.000000000026363898E+00 /
  data bth0cs(12) /   +0.000000000007297935E+00 /
  data bth0cs(13) /   -0.000000000002144188E+00 /
  data bth0cs(14) /   +0.000000000000663693E+00 /
  data bth0cs(15) /   -0.000000000000215126E+00 /
  data bth0cs(16) /   +0.000000000000072659E+00 /
  data bth0cs(17) /   -0.000000000000025465E+00 /
  data bth0cs(18) /   +0.000000000000009229E+00 /
  data bth0cs(19) /   -0.000000000000003448E+00 /
  data bth0cs(20) /   +0.000000000000001325E+00 /
  data bth0cs(21) /   -0.000000000000000522E+00 /
  data bth0cs(22) /   +0.000000000000000210E+00 /
  data bth0cs(23) /   -0.000000000000000087E+00 /
  data bth0cs(24) /   +0.000000000000000036E+00 /

  data twodpi / 0.63661977236758134E+00 /
  data pi4 / 0.78539816339744831E+00 /
  data alnhaf / -0.693147180559945309E+00 /
  data nty0 / 0 /
  data ntm0 / 0 /
  data ntth0 / 0 /
  data xsml / 0.0E+00 /
  data xmax / 0.0E+00 /

  if ( nty0 == 0 ) then
    nty0 = r4_inits ( by0cs, 13, 0.1E+00 * r4_mach ( 3 ) )
    ntm0 = r4_inits ( bm0cs, 21, 0.1E+00 * r4_mach ( 3 ) )
    ntth0 = r4_inits ( bth0cs, 24, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) )
    xmax = 1.0E+00 / r4_mach ( 4 )
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESY0 - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0E+00
    r4_besy0 = twodpi * ( alnhaf + log ( x ) ) &
      * r4_besj0 ( x ) + 0.375E+00 &
      + r4_csevl ( 0.125E+00 * y - 1.0E+00, by0cs, nty0 )
  else if ( x <= 4.0E+00 ) then
    y = x * x
    r4_besy0 = twodpi * ( alnhaf + log ( x ) ) &
      * r4_besj0 ( x ) + 0.375E+00 &
      + r4_csevl ( 0.125E+00 * y - 1.0E+00, by0cs, nty0 )
  else 
    z = 32.0E+00 / x / x - 1.0E+00
    ampl = ( 0.75E+00 + r4_csevl ( z, bm0cs, ntm0 ) ) / sqrt ( x )
    theta = x - pi4 + r4_csevl ( z, bth0cs, ntth0 ) / x
    r4_besy0 = ampl * sin ( theta )
  end if

  return
end
function r4_besy1 ( x )

!*****************************************************************************80
!
!! R4_BESY1 evaluates the Bessel function Y of order 1 of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BESY1, the Bessel function Y of order 1 of X.
!
  implicit none

  real ( kind = 4 ) ampl
  real ( kind = 4 ) bm1cs(21)
  real ( kind = 4 ) bth1cs(24)
  real ( kind = 4 ) by1cs(14)
  integer ( kind = 4 ) ntm1
  integer ( kind = 4 ) ntth1
  integer ( kind = 4 ) nty1
  real ( kind = 4 ) pi4
  real ( kind = 4 ) r4_besj1
  real ( kind = 4 ) r4_besy1
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) theta
  real ( kind = 4 ) twodpi
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  save bm1cs
  save bth1cs
  save by1cs
  save ntm1
  save ntth1
  save nty1
  save pi4
  save twodpi
  save xmax
  save xmin
  save xsml

  data by1cs( 1) /   +0.03208047100611908629E+00 /
  data by1cs( 2) /   +1.262707897433500450E+00 /
  data by1cs( 3) /   +0.00649996189992317500E+00 /
  data by1cs( 4) /   -0.08936164528860504117E+00 /
  data by1cs( 5) /   +0.01325088122175709545E+00 /
  data by1cs( 6) /   -0.00089790591196483523E+00 /
  data by1cs( 7) /   +0.00003647361487958306E+00 /
  data by1cs( 8) /   -0.00000100137438166600E+00 /
  data by1cs( 9) /   +0.00000001994539657390E+00 /
  data by1cs(10) /   -0.00000000030230656018E+00 /
  data by1cs(11) /   +0.00000000000360987815E+00 /
  data by1cs(12) /   -0.00000000000003487488E+00 /
  data by1cs(13) /   +0.00000000000000027838E+00 /
  data by1cs(14) /   -0.00000000000000000186E+00 /

  data bm1cs( 1) /   +0.1047362510931285E+00 /
  data bm1cs( 2) /   +0.00442443893702345E+00 /
  data bm1cs( 3) /   -0.00005661639504035E+00 /
  data bm1cs( 4) /   +0.00000231349417339E+00 /
  data bm1cs( 5) /   -0.00000017377182007E+00 /
  data bm1cs( 6) /   +0.00000001893209930E+00 /
  data bm1cs( 7) /   -0.00000000265416023E+00 /
  data bm1cs( 8) /   +0.00000000044740209E+00 /
  data bm1cs( 9) /   -0.00000000008691795E+00 /
  data bm1cs(10) /   +0.00000000001891492E+00 /
  data bm1cs(11) /   -0.00000000000451884E+00 /
  data bm1cs(12) /   +0.00000000000116765E+00 /
  data bm1cs(13) /   -0.00000000000032265E+00 /
  data bm1cs(14) /   +0.00000000000009450E+00 /
  data bm1cs(15) /   -0.00000000000002913E+00 /
  data bm1cs(16) /   +0.00000000000000939E+00 /
  data bm1cs(17) /   -0.00000000000000315E+00 /
  data bm1cs(18) /   +0.00000000000000109E+00 /
  data bm1cs(19) /   -0.00000000000000039E+00 /
  data bm1cs(20) /   +0.00000000000000014E+00 /
  data bm1cs(21) /   -0.00000000000000005E+00 /

  data bth1cs( 1) /   +0.74060141026313850E+00 /
  data bth1cs( 2) /   -0.004571755659637690E+00 /
  data bth1cs( 3) /   +0.000119818510964326E+00 /
  data bth1cs( 4) /   -0.000006964561891648E+00 /
  data bth1cs( 5) /   +0.000000655495621447E+00 /
  data bth1cs( 6) /   -0.000000084066228945E+00 /
  data bth1cs( 7) /   +0.000000013376886564E+00 /
  data bth1cs( 8) /   -0.000000002499565654E+00 /
  data bth1cs( 9) /   +0.000000000529495100E+00 /
  data bth1cs(10) /   -0.000000000124135944E+00 /
  data bth1cs(11) /   +0.000000000031656485E+00 /
  data bth1cs(12) /   -0.000000000008668640E+00 /
  data bth1cs(13) /   +0.000000000002523758E+00 /
  data bth1cs(14) /   -0.000000000000775085E+00 /
  data bth1cs(15) /   +0.000000000000249527E+00 /
  data bth1cs(16) /   -0.000000000000083773E+00 /
  data bth1cs(17) /   +0.000000000000029205E+00 /
  data bth1cs(18) /   -0.000000000000010534E+00 /
  data bth1cs(19) /   +0.000000000000003919E+00 /
  data bth1cs(20) /   -0.000000000000001500E+00 /
  data bth1cs(21) /   +0.000000000000000589E+00 /
  data bth1cs(22) /   -0.000000000000000237E+00 /
  data bth1cs(23) /   +0.000000000000000097E+00 /
  data bth1cs(24) /   -0.000000000000000040E+00 /

  data ntm1 / 0 /
  data ntth1 / 0 /
  data nty1 / 0 /
  data pi4 / 0.78539816339744831E+00 /
  data twodpi / 0.63661977236758134E+00 /
  data xmax / 0.0E+00 /
  data xmin / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( nty1 == 0 ) then
    nty1 = r4_inits ( by1cs, 14, 0.1E+00 * r4_mach ( 3 ) )
    ntm1 = r4_inits ( bm1cs, 21, 0.1E+00 * r4_mach ( 3 ) )
    ntth1 = r4_inits ( bth1cs, 24, 0.1E+00 * r4_mach ( 3 ) )
    xmin = 1.571E+00 * exp ( max ( log ( r4_mach ( 1 ) ), &
      - log ( r4_mach ( 2 ) ) ) + 0.01E+00 )
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) )
    xmax = 1.0E+00 / r4_mach ( 4 )
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BESY1 - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0E+00
    r4_besy1 = twodpi * log ( 0.5E+00 * x ) * r4_besj1 ( x ) &
      + ( 0.5E+00 + r4_csevl ( 0.125E+00 * y - 1.0E+00, by1cs, &
      nty1 ) ) / x
  else if ( x <= 4.0E+00 ) then
    y = x * x
    r4_besy1 = twodpi * log ( 0.5E+00 * x ) * r4_besj1 ( x ) &
      + ( 0.5E+00 + r4_csevl ( 0.125E+00 * y - 1.0E+00, by1cs, &
      nty1 ) ) / x
  else
    z = 32.0E+00 / x / x - 1.0E+00
    ampl = ( 0.75E+00 + r4_csevl ( z, bm1cs, ntm1 ) ) / sqrt ( x )
    theta = x - 3.0E+00 * pi4 + r4_csevl ( z, bth1cs, ntth1 ) / x
    r4_besy1 = ampl * sin ( theta )
  end if

  return
end
function r4_beta ( a, b )

!*****************************************************************************80
!
!! R4_BETA evaluates the beta function of R4 arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the arguments.
!
!    Output, real ( kind = 4 ) R4_BETA, the beta function of A and B.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) alnsml
  real ( kind = 4 ) b
  real ( kind = 4 ) r4_beta
  real ( kind = 4 ) r4_gamma
  real ( kind = 4 ) r4_lbeta
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin

  save alnsml
  save xmax

  data alnsml / 0.0E+00 /
  data xmax / 0.0E+00 /

  if ( xmax == 0.0E+00 ) then
    call r4_gaml ( xmin, xmax )
    alnsml = log ( r4_mach ( 1 ) )
  end if

  if ( a <= 0.0E+00 .or. b <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BETA - Fatal error!'
    write ( *, '(a)' ) '  A and B must be greater than 0.'
    write ( *, '(a,g14.6)' ) '  A = ', a
    write ( *, '(a,g14.6)' ) '  B = ', b
    stop
  end if

  if ( a + b < xmax ) then
    r4_beta = r4_gamma ( a ) * r4_gamma ( b ) / r4_gamma ( a + b )
    return
  end if

  r4_beta = r4_lbeta ( a, b )

  r4_beta = exp ( r4_beta )

  return
end
function r4_betai ( x, pin, qin )

!*****************************************************************************80
!
!! R4_BETAI evaluates the incomplete beta ratio of R4 arguments.
!
!  Discussion:
!
!    The incomplete Beta function ratio is the probability that a
!    random variable from a beta distribution having parameters
!    P and Q will be less than or equal to X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Nancy Bosten, EL Battiste,
!    Remark on Algorithm 179: 
!    Incomplete Beta Ratio,
!    Communications of the ACM,
!    Volume 17, Number 3, March 1974, pages 156-157.
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the upper limit of integration.
!    0.0 <= X <= 1.0.
!
!    Input, real ( kind = 4 ) PIN, the first distribution parameter.
!    0.0 < PIN.
!
!    Input, real ( kind = 4 ) QIN, the second distribution parameter.
!    0.0 < QIN.
!
!    Output, real ( kind = 4 ) R4_BETAI, the incomplete beta function ratio.
!
  implicit none

  real ( kind = 4 ) alneps
  real ( kind = 4 ) alnsml
  real ( kind = 4 ) c
  real ( kind = 4 ) eps
  real ( kind = 4 ) finsum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) n
  real ( kind = 4 ) p
  real ( kind = 4 ) p1
  real ( kind = 4 ) pin
  real ( kind = 4 ) ps
  real ( kind = 4 ) q
  real ( kind = 4 ) qin
  real ( kind = 4 ) r4_betai
  real ( kind = 4 ) r4_lbeta
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sml
  real ( kind = 4 ) term
  real ( kind = 4 ) x
  real ( kind = 4 ) xb
  real ( kind = 4 ) y

  save alneps
  save alnsml
  save eps
  save sml

  data alneps / 0.0E+00 /
  data alnsml / 0.0E+00 /
  data eps / 0.0E+00 /
  data sml / 0.0E+00 /

  if ( eps == 0.0E+00 ) then
    eps = r4_mach ( 3 )
    alneps = log ( eps )
    sml = r4_mach ( 1 )
    alnsml = log ( sml )
  end if

  if ( x < 0.0E+00 .or. 1.0E+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BETAI - Fatal error!'
    write ( *, '(a)' ) '  0 <= X <= 1 is required.'
    stop
  end if

  if ( pin <= 0.0E+00 .or. qin <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BETAI - Fatal error!'
    write ( *, '(a)' ) '  P or Q <= 0.0.'
    stop
  end if

  y = x
  p = pin
  q = qin

  if ( p < q .or. 0.8E+00 <= x ) then

    if ( 0.2E+00 <= x ) then
      y = 1.0E+00 - y
      p = qin
      q = pin
    end if

  end if

  if ( ( p + q ) * y / ( p + 1.0E+00 ) < eps ) then

    r4_betai = 0.0E+00

    xb = p * log ( max ( y, sml ) ) - log ( p ) - r4_lbeta ( p, q )

    if ( alnsml < xb .and. y /= 0.0E+00 ) then
      r4_betai = exp ( xb )
    end if

    if ( y /= x .or. p /= pin ) then
      r4_betai = 1.0E+00 - r4_betai
    end if

    return

  end if
!
!  Evaluate the infinite sum first.
!  TERM will equal y^p/beta(ps,p) * (1.-ps)i * y^i / fac(i)
!
  ps = q - aint ( q )
  if ( ps == 0.0E+00 ) then
    ps = 1.0E+00
  end if

  xb = p * log ( y ) - r4_lbeta ( ps, p ) - log ( p )

  if ( xb < alnsml ) then

    r4_betai = 0.0E+00

  else

    r4_betai = exp ( xb )
    term = r4_betai * p

    if ( ps /= 1.0E+00 ) then

      n = int ( max ( alneps / log ( y ), 4.0E+00 ) )
      do i = 1, n
        term = term * ( real ( i, kind = 4 ) - ps ) * y / real ( i, kind = 4 )
        r4_betai = r4_betai + term / ( p + real ( i, kind = 4 ) )
      end do

    end if

  end if
!
!  Now evaluate the finite sum.
!
  if  ( 1.0E+00 < q ) then

    xb = p * log ( y ) + q * log ( 1.0E+00 - y ) - r4_lbeta ( p, q ) - log ( q )
    ib = int ( max ( xb / alnsml, 0.0E+00 ) )
    term = exp ( xb - real ( ib, kind = 4 ) * alnsml )
    c = 1.0E+00 / ( 1.0E+00 - y )
    p1 = q * c / ( p + q - 1.0E+00 )

    finsum = 0.0E+00
    n = int ( q )
    if ( q == real ( n, kind = 4 ) ) then
      n = n - 1
    end if

    do i = 1, n

      if ( p1 <= 1.0E+00 .and. term / eps <= finsum ) then
        exit
      end if

      term = ( q - real ( i - 1, kind = 4 ) ) * c * term &
        / ( p + q - real ( i, kind = 4 ) )

      if ( 1.0E+00 < term ) then
        ib = ib - 1
        term = term * sml
      end if

      if ( ib == 0 ) then
        finsum = finsum + term
      end if

    end do

    r4_betai = r4_betai + finsum

  end if

  if ( y /= x .or. p /= pin ) then
    r4_betai = 1.0E+00 - r4_betai
  end if

  if ( r4_betai < 0.0E+00 ) then
    r4_betai =  0.0E+00
  end if

  if ( 1.0E+00 < r4_betai ) then
    r4_betai = 1.0E+00
  end if

  return
end
function r4_bi ( x )

!*****************************************************************************80
!
!! R4_BI evaluates the Airy function Bi of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BI, the Airy function Bi of X.
!
  implicit none

  real ( kind = 4 ) bif2cs(10)
  real ( kind = 4 ) bifcs(9)
  real ( kind = 4 ) big2cs(10)
  real ( kind = 4 ) bigcs(8)
  real ( kind = 4 ) eta
  integer ( kind = 4 ) nbif
  integer ( kind = 4 ) nbif2
  integer ( kind = 4 ) nbig
  integer ( kind = 4 ) nbig2
  real ( kind = 4 ) r4_bi
  real ( kind = 4 ) r4_bie
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) theta
  real ( kind = 4 ) x
  real ( kind = 4 ) x3sml
  real ( kind = 4 ) xm
  real ( kind = 4 ) xmax
  real ( kind = 4 ) z

  save bif2cs
  save bifcs
  save big2cs
  save bigcs
  save nbif
  save nbif2
  save nbig
  save nbig2
  save x3sml
  save xmax

  data bifcs( 1) /   -0.01673021647198664948E+00 /
  data bifcs( 2) /    0.1025233583424944561E+00 /
  data bifcs( 3) /    0.00170830925073815165E+00 /
  data bifcs( 4) /    0.00001186254546774468E+00 /
  data bifcs( 5) /    0.00000004493290701779E+00 /
  data bifcs( 6) /    0.00000000010698207143E+00 /
  data bifcs( 7) /    0.00000000000017480643E+00 /
  data bifcs( 8) /    0.00000000000000020810E+00 /
  data bifcs( 9) /    0.00000000000000000018E+00 /

  data bigcs( 1) /    0.02246622324857452E+00 /
  data bigcs( 2) /    0.03736477545301955E+00 /
  data bigcs( 3) /    0.00044476218957212E+00 /
  data bigcs( 4) /    0.00000247080756363E+00 /
  data bigcs( 5) /    0.00000000791913533E+00 /
  data bigcs( 6) /    0.00000000001649807E+00 /
  data bigcs( 7) /    0.00000000000002411E+00 /
  data bigcs( 8) /    0.00000000000000002E+00 /

  data bif2cs( 1) /    0.09984572693816041E+00 /
  data bif2cs( 2) /    0.478624977863005538E+00 /
  data bif2cs( 3) /    0.0251552119604330118E+00 /
  data bif2cs( 4) /    0.0005820693885232645E+00 /
  data bif2cs( 5) /    0.0000074997659644377E+00 /
  data bif2cs( 6) /    0.0000000613460287034E+00 /
  data bif2cs( 7) /    0.0000000003462753885E+00 /
  data bif2cs( 8) /    0.0000000000014288910E+00 /
  data bif2cs( 9) /    0.0000000000000044962E+00 /
  data bif2cs(10) /    0.0000000000000000111E+00 /

  data big2cs( 1) /    0.033305662145514340E+00 /
  data big2cs( 2) /    0.161309215123197068E+00 /
  data big2cs( 3) /    0.0063190073096134286E+00 /
  data big2cs( 4) /    0.0001187904568162517E+00 /
  data big2cs( 5) /    0.0000013045345886200E+00 /
  data big2cs( 6) /    0.0000000093741259955E+00 /
  data big2cs( 7) /    0.0000000000474580188E+00 /
  data big2cs( 8) /    0.0000000000001783107E+00 /
  data big2cs( 9) /    0.0000000000000005167E+00 /
  data big2cs(10) /    0.0000000000000000011E+00 /

  data nbif / 0 /
  data nbif2 / 0 /
  data nbig / 0 /
  data nbig2 / 0 /
  data x3sml / 0.0E+00 /
  data xmax / 0.0E+00 /

  if ( nbif == 0 ) then
    eta = 0.1E+00 * r4_mach ( 3 )
    nbif = r4_inits ( bifcs, 9, eta )
    nbig = r4_inits ( bigcs, 8, eta )
    nbif2 = r4_inits ( bif2cs, 10, eta )
    nbig2 = r4_inits ( big2cs, 10, eta )
    x3sml = eta ** 0.3333E+00
    xmax = ( 1.5E+00 * log ( r4_mach ( 2 ) ) ) ** 0.6666E+00
  end if

  if ( x <= - 1.0E+00 ) then

    call r4_aimp ( x, xm, theta )
    r4_bi = xm * sin ( theta )

  else if ( abs ( x ) <= x3sml ) then

    z = 0.0E+00

    r4_bi = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) &
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) )

  else if ( x <= 1.0E+00 ) then

    z = x * x * x
    r4_bi = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) &
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) )

  else if ( x <= 2.0E+00 ) then

    z = ( 2.0E+00 * x * x * x - 9.0E+00 ) / 7.0E+00
    r4_bi = 1.125E+00 + r4_csevl ( z, bif2cs, nbif2 ) &
      + x * ( 0.625E+00 + r4_csevl ( z, big2cs, nbig2 ) )

  else

    r4_bi = r4_bie ( x ) &
      * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )

  end if

  return
end
function r4_bid ( x )

!*****************************************************************************80
!
!! R4_BID evaluates the derivative of the Airy function Bi of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BID, the derivative of the Airy function Bi.
!
  implicit none

  real ( kind = 4 ) bif2cs(10)
  real ( kind = 4 ) bifcs(8)
  real ( kind = 4 ) big2cs(10)
  real ( kind = 4 ) bigcs(9)
  real ( kind = 4 ) eta
  integer ( kind = 4 ) nbif
  integer ( kind = 4 ) nbif2
  integer ( kind = 4 ) nbig
  integer ( kind = 4 ) nbig2
  real ( kind = 4 ) phi
  real ( kind = 4 ) r4_bid
  real ( kind = 4 ) r4_bide
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) x2
  real ( kind = 4 ) x2sml
  real ( kind = 4 ) x3
  real ( kind = 4 ) x3sml
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xn
  real ( kind = 4 ) z

  save bif2cs
  save bifcs
  save big2cs
  save bigcs
  save nbif
  save nbif2
  save nbig
  save nbig2
  save x2sml
  save x3sml
  save xmax

  data bifcs(  1) /     0.1153536790828570243E+00 /
  data bifcs(  2) /     0.0205007894049192875E+00 /
  data bifcs(  3) /     0.0002135290278902876E+00 /
  data bifcs(  4) /     0.0000010783960614677E+00 /
  data bifcs(  5) /     0.0000000032094708833E+00 /
  data bifcs(  6) /     0.0000000000062930407E+00 /
  data bifcs(  7) /     0.0000000000000087403E+00 /
  data bifcs(  8) /     0.0000000000000000090E+00 /

  data bigcs(  1) /    -0.097196440416443537390E+00 /
  data bigcs(  2) /     0.149503576843167066571E+00 /
  data bigcs(  3) /     0.003113525387121326042E+00 /
  data bigcs(  4) /     0.000024708570579821297E+00 /
  data bigcs(  5) /     0.000000102949627731379E+00 /
  data bigcs(  6) /     0.000000000263970373987E+00 /
  data bigcs(  7) /     0.000000000000458279271E+00 /
  data bigcs(  8) /     0.000000000000000574283E+00 /
  data bigcs(  9) /     0.000000000000000000544E+00 /

  data bif2cs(  1) /     0.323493987603522033521E+00 /
  data bif2cs(  2) /     0.086297871535563559139E+00 /
  data bif2cs(  3) /     0.002994025552655397426E+00 /
  data bif2cs(  4) /     0.000051430528364661637E+00 /
  data bif2cs(  5) /     0.000000525840250036811E+00 /
  data bif2cs(  6) /     0.000000003561751373958E+00 /
  data bif2cs(  7) /     0.000000000017146864007E+00 /
  data bif2cs(  8) /     0.000000000000061663520E+00 /
  data bif2cs(  9) /     0.000000000000000171911E+00 /
  data bif2cs( 10) /     0.000000000000000000382E+00 /

  data big2cs(  1) /     1.6062999463621294578E+00 /
  data big2cs(  2) /     0.7449088819876088652E+00 /
  data big2cs(  3) /     0.0470138738610277380E+00 /
  data big2cs(  4) /     0.0012284422062548239E+00 /
  data big2cs(  5) /     0.0000173222412256624E+00 /
  data big2cs(  6) /     0.0000001521901652368E+00 /
  data big2cs(  7) /     0.0000000009113560249E+00 /
  data big2cs(  8) /     0.0000000000039547918E+00 /
  data big2cs(  9) /     0.0000000000000130017E+00 /
  data big2cs( 10) /     0.0000000000000000335E+00 /

  data nbif / 0/
  data nbif2 / 0 /
  data nbig / 0 /
  data nbig2 / 0 /
  data x2sml / 0.0E+00 /
  data x3sml / 0.0E+00 /
  data xmax / 0.0E+00 /

  if ( nbif == 0 ) then
    eta = 0.1E+00 * r4_mach ( 3 )
    nbif = r4_inits ( bifcs, 8, eta )
    nbig = r4_inits ( bigcs, 9, eta )
    nbif2 = r4_inits ( bif2cs, 10, eta )
    nbig2 = r4_inits ( big2cs, 10, eta )
    x2sml = sqrt ( eta )
    x3sml = eta ** 0.3333E+00
    xmax = ( 1.5E+00 * log ( r4_mach ( 2 ) ) ) ** 0.6666E+00
  end if

  if ( x < - 1.0E+00 ) then
    call r4_admp ( x, xn, phi )
    r4_bid = xn * sin ( phi )
  else if ( abs ( x ) <= x2sml ) then
    x2 = 0.0E+00
    x3 = 0.0E+00
    r4_bid = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) &
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00
  else if ( abs ( x ) <= x3sml ) then
    x2 = x * x
    x3 = 0.0E+00
    r4_bid = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) &
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00
  else if ( x <= 1.0E+00 ) then
    x2 = x * x
    x3 = x * x * x
    r4_bid = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) &
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00
  else if ( x <= 2.0E+00 ) then
    z = ( 2.0E+00 * x * x * x - 9.0E+00 ) / 7.0E+00
    r4_bid = x * x * ( r4_csevl ( z, bif2cs, nbif2 ) + 0.25E+00 ) &
      + r4_csevl ( z, big2cs, nbig2 ) + 0.5E+00
  else
    r4_bid = r4_bide ( x ) &
      * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  end if

  return
end
function r4_bide ( x )

!*****************************************************************************80
!
!! R4_BIDE: exponentially scaled derivative, Airy function Bi of an R4 argument.
!
!  Discussion:
!
!    if X < 0,
!      R4_BIDE ( X ) = R4_BID ( X )
!    else
!      R4_BIDE ( X ) = R4_BID ( X ) * exp ( - 2/3 * X**(3/2) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BIDE, the exponentially scaled derivative of 
!    the Airy function Bi of X.
!
  implicit none

  real ( kind = 4 ) atr
  real ( kind = 4 ) bif2cs(10)
  real ( kind = 4 ) bifcs(8)
  real ( kind = 4 ) big2cs(10)
  real ( kind = 4 ) bigcs(9)
  real ( kind = 4 ) bip1cs(24)
  real ( kind = 4 ) bip2cs(29)
  real ( kind = 4 ) btr
  real ( kind = 4 ) eta
  integer ( kind = 4 ) nbif
  integer ( kind = 4 ) nbif2
  integer ( kind = 4 ) nbig
  integer ( kind = 4 ) nbig2
  integer ( kind = 4 ) nbip1
  integer ( kind = 4 ) nbip2
  real ( kind = 4 ) phi
  real ( kind = 4 ) r4_bide
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sqrtx
  real ( kind = 4 ) x
  real ( kind = 4 ) x2
  real ( kind = 4 ) x2sml
  real ( kind = 4 ) x3
  real ( kind = 4 ) x3sml
  real ( kind = 4 ) x32sml
  real ( kind = 4 ) xbig
  real ( kind = 4 ) xn
  real ( kind = 4 ) z

  save atr
  save bif2cs
  save bifcs
  save big2cs
  save bigcs
  save bip1cs
  save bip2cs
  save btr
  save nbif
  save nbif2
  save nbig
  save nbig2
  save nbip1
  save nbip2
  save x2sml
  save x3sml
  save x32sml
  save xbig

  data bifcs(  1) /     0.1153536790828570243E+00 /
  data bifcs(  2) /     0.0205007894049192875E+00 /
  data bifcs(  3) /     0.0002135290278902876E+00 /
  data bifcs(  4) /     0.0000010783960614677E+00 /
  data bifcs(  5) /     0.0000000032094708833E+00 /
  data bifcs(  6) /     0.0000000000062930407E+00 /
  data bifcs(  7) /     0.0000000000000087403E+00 /
  data bifcs(  8) /     0.0000000000000000090E+00 /

  data bigcs(  1) /    -0.097196440416443537390E+00 /
  data bigcs(  2) /     0.149503576843167066571E+00 /
  data bigcs(  3) /     0.003113525387121326042E+00 /
  data bigcs(  4) /     0.000024708570579821297E+00 /
  data bigcs(  5) /     0.000000102949627731379E+00 /
  data bigcs(  6) /     0.000000000263970373987E+00 /
  data bigcs(  7) /     0.000000000000458279271E+00 /
  data bigcs(  8) /     0.000000000000000574283E+00 /
  data bigcs(  9) /     0.000000000000000000544E+00 /

  data bif2cs(  1) /     0.323493987603522033521E+00 /
  data bif2cs(  2) /     0.086297871535563559139E+00 /
  data bif2cs(  3) /     0.002994025552655397426E+00 /
  data bif2cs(  4) /     0.000051430528364661637E+00 /
  data bif2cs(  5) /     0.000000525840250036811E+00 /
  data bif2cs(  6) /     0.000000003561751373958E+00 /
  data bif2cs(  7) /     0.000000000017146864007E+00 /
  data bif2cs(  8) /     0.000000000000061663520E+00 /
  data bif2cs(  9) /     0.000000000000000171911E+00 /
  data bif2cs( 10) /     0.000000000000000000382E+00 /

  data big2cs(  1) /     1.6062999463621294578E+00 /
  data big2cs(  2) /     0.7449088819876088652E+00 /
  data big2cs(  3) /     0.0470138738610277380E+00 /
  data big2cs(  4) /     0.0012284422062548239E+00 /
  data big2cs(  5) /     0.0000173222412256624E+00 /
  data big2cs(  6) /     0.0000001521901652368E+00 /
  data big2cs(  7) /     0.0000000009113560249E+00 /
  data big2cs(  8) /     0.0000000000039547918E+00 /
  data big2cs(  9) /     0.0000000000000130017E+00 /
  data big2cs( 10) /     0.0000000000000000335E+00 /

  data bip2cs(  1) /    -0.13269705443526630495E+00 /
  data bip2cs(  2) /    -0.00568443626045977481E+00 /
  data bip2cs(  3) /    -0.00015643601119611610E+00 /
  data bip2cs(  4) /    -0.00001136737203679562E+00 /
  data bip2cs(  5) /    -0.00000143464350991284E+00 /
  data bip2cs(  6) /    -0.00000018098531185164E+00 /
  data bip2cs(  7) /     0.00000000926177343611E+00 /
  data bip2cs(  8) /     0.00000001710005490721E+00 /
  data bip2cs(  9) /     0.00000000476698163504E+00 /
  data bip2cs( 10) /    -0.00000000035195022023E+00 /
  data bip2cs( 11) /    -0.00000000058890614316E+00 /
  data bip2cs( 12) /    -0.00000000006678499608E+00 /
  data bip2cs( 13) /     0.00000000006395565102E+00 /
  data bip2cs( 14) /     0.00000000001554529427E+00 /
  data bip2cs( 15) /    -0.00000000000792397000E+00 /
  data bip2cs( 16) /    -0.00000000000258326243E+00 /
  data bip2cs( 17) /     0.00000000000121655048E+00 /
  data bip2cs( 18) /     0.00000000000038707207E+00 /
  data bip2cs( 19) /    -0.00000000000022487045E+00 /
  data bip2cs( 20) /    -0.00000000000004953477E+00 /
  data bip2cs( 21) /     0.00000000000004563782E+00 /
  data bip2cs( 22) /     0.00000000000000332998E+00 /
  data bip2cs( 23) /    -0.00000000000000921750E+00 /
  data bip2cs( 24) /     0.00000000000000094157E+00 /
  data bip2cs( 25) /     0.00000000000000167154E+00 /
  data bip2cs( 26) /    -0.00000000000000055134E+00 /
  data bip2cs( 27) /    -0.00000000000000022369E+00 /
  data bip2cs( 28) /     0.00000000000000017487E+00 /
  data bip2cs( 29) /     0.00000000000000000207E+00 /

  data bip1cs(  1) /    -0.1729187351079553719E+00 /
  data bip1cs(  2) /    -0.0149358492984694364E+00 /
  data bip1cs(  3) /    -0.0005471104951678566E+00 /
  data bip1cs(  4) /     0.0001537966292958408E+00 /
  data bip1cs(  5) /     0.0000154353476192179E+00 /
  data bip1cs(  6) /    -0.0000065434113851906E+00 /
  data bip1cs(  7) /     0.0000003728082407879E+00 /
  data bip1cs(  8) /     0.0000002072078388189E+00 /
  data bip1cs(  9) /    -0.0000000658173336470E+00 /
  data bip1cs( 10) /     0.0000000074926746354E+00 /
  data bip1cs( 11) /     0.0000000011101336884E+00 /
  data bip1cs( 12) /    -0.0000000007265140553E+00 /
  data bip1cs( 13) /     0.0000000001782723560E+00 /
  data bip1cs( 14) /    -0.0000000000217346352E+00 /
  data bip1cs( 15) /    -0.0000000000020302035E+00 /
  data bip1cs( 16) /     0.0000000000019311827E+00 /
  data bip1cs( 17) /    -0.0000000000006044953E+00 /
  data bip1cs( 18) /     0.0000000000001209450E+00 /
  data bip1cs( 19) /    -0.0000000000000125109E+00 /
  data bip1cs( 20) /    -0.0000000000000019917E+00 /
  data bip1cs( 21) /     0.0000000000000015154E+00 /
  data bip1cs( 22) /    -0.0000000000000004977E+00 /
  data bip1cs( 23) /     0.0000000000000001155E+00 /
  data bip1cs( 24) /    -0.0000000000000000186E+00 /

  data atr / 8.7506905708484345E+00 /
  data btr /-2.0938363213560543E+00 /
  data nbif / 0 /
  data nbif2 / 0 /
  data nbig / 0 /
  data nbig2 / 0 /
  data nbip1 / 0 /
  data nbip2 / 0 /
  data x2sml / 0.0E+00 /
  data x3sml / 0.0E+00 /
  data x32sml / 0.0E+00 /
  data xbig / 0.0E+00 /

  if ( nbif == 0 ) then
    eta = 0.1E+00 * r4_mach ( 3 )
    nbif = r4_inits ( bifcs, 8, eta )
    nbig = r4_inits ( bigcs, 9, eta )
    nbif2 = r4_inits ( bif2cs, 10, eta )
    nbig2 = r4_inits ( big2cs, 10, eta )
    nbip1 = r4_inits ( bip1cs, 24, eta )
    nbip2 = r4_inits ( bip2cs, 29, eta )
    x2sml = sqrt ( eta )
    x3sml = eta ** 0.3333E+00
    x32sml = 1.3104 * x3sml * x3sml
    xbig = r4_mach ( 2 ) ** 0.6666E+00
  end if

  if ( x <= - 1.0E+00 ) then
    call r4_admp ( x, xn, phi )
    r4_bide = xn * sin ( phi )
  else if ( 0.0E+00 <= x .and. x <= x32sml ) then
    x2 = 0.0E+00
    x3 = 0.0E+00
    r4_bide = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) &
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00
  else if ( abs ( x ) <= x2sml ) then
    x2 = 0.0E+00
    x3 = 0.0E+00
    r4_bide = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) &
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00
    r4_bide = r4_bide * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  else if ( x <= x3sml ) then
    x2 = x * x
    x3 = 0.0E+00
    r4_bide = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) &
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00
    r4_bide = r4_bide * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  else if ( x <= 1.0E+00 ) then
    x2 = x * x
    x3 = x * x * x
    r4_bide = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) &
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00
    r4_bide = r4_bide * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  else if ( x <= 2.0E+00 ) then
    z = ( 2.0E+00 * x * x * x - 9.0E+00 ) / 7.0E+00
    r4_bide = exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 ) &
      * ( x * x * ( 0.25E+00 + r4_csevl ( z, bif2cs, nbif2 ) ) &
      + 0.5E+00 + r4_csevl ( z, big2cs, nbig2 ) )
  else if ( x <= 4.0E+00 ) then
    sqrtx = sqrt ( x )
    z = atr / ( x * sqrtx ) + btr
    r4_bide = ( 0.625E+00 &
      + r4_csevl ( z, bip1cs, nbip1 ) ) * sqrt ( sqrtx )
  else if ( x < xbig ) then
    sqrtx = sqrt ( x )
    z = 16.0E+00 / ( x * sqrtx ) - 1.0E+00
    r4_bide = ( 0.625E+00 + r4_csevl ( z, bip2cs, nbip2 ) ) &
      * sqrt ( sqrtx )
  else
    sqrtx = sqrt ( x )
    z = - 1.0E+00
    r4_bide = ( 0.625E+00 + r4_csevl ( z, bip2cs, nbip2 ) ) * sqrt ( sqrtx )
  end if

  return
end
function r4_bie ( x )

!*****************************************************************************80
!
!! R4_BIE evaluates the exponentially scaled Airy function Bi of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_BIE, the exponentially scaled Airy 
!    function Bi of X.
!
  implicit none

  real ( kind = 4 ) atr
  real ( kind = 4 ) bif2cs(10)
  real ( kind = 4 ) bifcs(9)
  real ( kind = 4 ) big2cs(10)
  real ( kind = 4 ) bigcs(8)
  real ( kind = 4 ) bip2cs(29)
  real ( kind = 4 ) bipcs(24)
  real ( kind = 4 ) btr
  real ( kind = 4 ) eta
  integer ( kind = 4 ) nbif
  integer ( kind = 4 ) nbif2
  integer ( kind = 4 ) nbig
  integer ( kind = 4 ) nbig2
  integer ( kind = 4 ) nbip
  integer ( kind = 4 ) nbip2
  real ( kind = 4 ) r4_bie
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sqrtx
  real ( kind = 4 ) theta
  real ( kind = 4 ) x
  real ( kind = 4 ) x32sml
  real ( kind = 4 ) x3sml
  real ( kind = 4 ) xbig
  real ( kind = 4 ) xm
  real ( kind = 4 ) z

  save atr
  save bif2cs
  save bifcs
  save big2cs
  save bigcs
  save bip2cs
  save bipcs
  save btr
  save nbif
  save nbif2
  save nbig
  save nbig2
  save nbip
  save nbip2
  save x32sml
  save x3sml
  save xbig

  data bifcs( 1) /   -0.01673021647198664948E+00 /
  data bifcs( 2) /    0.1025233583424944561E+00 /
  data bifcs( 3) /    0.00170830925073815165E+00 /
  data bifcs( 4) /    0.00001186254546774468E+00 /
  data bifcs( 5) /    0.00000004493290701779E+00 /
  data bifcs( 6) /    0.00000000010698207143E+00 /
  data bifcs( 7) /    0.00000000000017480643E+00 /
  data bifcs( 8) /    0.00000000000000020810E+00 /
  data bifcs( 9) /    0.00000000000000000018E+00 /

  data bigcs( 1) /    0.02246622324857452E+00 /
  data bigcs( 2) /    0.03736477545301955E+00 /
  data bigcs( 3) /    0.00044476218957212E+00 /
  data bigcs( 4) /    0.00000247080756363E+00 /
  data bigcs( 5) /    0.00000000791913533E+00 /
  data bigcs( 6) /    0.00000000001649807E+00 /
  data bigcs( 7) /    0.00000000000002411E+00 /
  data bigcs( 8) /    0.00000000000000002E+00 /

  data bif2cs( 1) /    0.09984572693816041E+00 /
  data bif2cs( 2) /    0.478624977863005538E+00 /
  data bif2cs( 3) /    0.0251552119604330118E+00 /
  data bif2cs( 4) /    0.0005820693885232645E+00 /
  data bif2cs( 5) /    0.0000074997659644377E+00 /
  data bif2cs( 6) /    0.0000000613460287034E+00 /
  data bif2cs( 7) /    0.0000000003462753885E+00 /
  data bif2cs( 8) /    0.0000000000014288910E+00 /
  data bif2cs( 9) /    0.0000000000000044962E+00 /
  data bif2cs(10) /    0.0000000000000000111E+00 /

  data big2cs( 1) /    0.033305662145514340E+00 /
  data big2cs( 2) /    0.161309215123197068E+00 /
  data big2cs( 3) /    0.0063190073096134286E+00 /
  data big2cs( 4) /    0.0001187904568162517E+00 /
  data big2cs( 5) /    0.0000013045345886200E+00 /
  data big2cs( 6) /    0.0000000093741259955E+00 /
  data big2cs( 7) /    0.0000000000474580188E+00 /
  data big2cs( 8) /    0.0000000000001783107E+00 /
  data big2cs( 9) /    0.0000000000000005167E+00 /
  data big2cs(10) /    0.0000000000000000011E+00 /

  data bipcs( 1) /   -0.08322047477943447E+00 /
  data bipcs( 2) /    0.01146118927371174E+00 /
  data bipcs( 3) /    0.00042896440718911E+00 /
  data bipcs( 4) /   -0.00014906639379950E+00 /
  data bipcs( 5) /   -0.00001307659726787E+00 /
  data bipcs( 6) /    0.00000632759839610E+00 /
  data bipcs( 7) /   -0.00000042226696982E+00 /
  data bipcs( 8) /   -0.00000019147186298E+00 /
  data bipcs( 9) /    0.00000006453106284E+00 /
  data bipcs(10) /   -0.00000000784485467E+00 /
  data bipcs(11) /   -0.00000000096077216E+00 /
  data bipcs(12) /    0.00000000070004713E+00 /
  data bipcs(13) /   -0.00000000017731789E+00 /
  data bipcs(14) /    0.00000000002272089E+00 /
  data bipcs(15) /    0.00000000000165404E+00 /
  data bipcs(16) /   -0.00000000000185171E+00 /
  data bipcs(17) /    0.00000000000059576E+00 /
  data bipcs(18) /   -0.00000000000012194E+00 /
  data bipcs(19) /    0.00000000000001334E+00 /
  data bipcs(20) /    0.00000000000000172E+00 /
  data bipcs(21) /   -0.00000000000000145E+00 /
  data bipcs(22) /    0.00000000000000049E+00 /
  data bipcs(23) /   -0.00000000000000011E+00 /
  data bipcs(24) /    0.00000000000000001E+00 /

  data bip2cs( 1) /   -0.113596737585988679E+00 /
  data bip2cs( 2) /    0.0041381473947881595E+00 /
  data bip2cs( 3) /    0.0001353470622119332E+00 /
  data bip2cs( 4) /    0.0000104273166530153E+00 /
  data bip2cs( 5) /    0.0000013474954767849E+00 /
  data bip2cs( 6) /    0.0000001696537405438E+00 /
  data bip2cs( 7) /   -0.0000000100965008656E+00 /
  data bip2cs( 8) /   -0.0000000167291194937E+00 /
  data bip2cs( 9) /   -0.0000000045815364485E+00 /
  data bip2cs(10) /    0.0000000003736681366E+00 /
  data bip2cs(11) /    0.0000000005766930320E+00 /
  data bip2cs(12) /    0.0000000000621812650E+00 /
  data bip2cs(13) /   -0.0000000000632941202E+00 /
  data bip2cs(14) /   -0.0000000000149150479E+00 /
  data bip2cs(15) /    0.0000000000078896213E+00 /
  data bip2cs(16) /    0.0000000000024960513E+00 /
  data bip2cs(17) /   -0.0000000000012130075E+00 /
  data bip2cs(18) /   -0.0000000000003740493E+00 /
  data bip2cs(19) /    0.0000000000002237727E+00 /
  data bip2cs(20) /    0.0000000000000474902E+00 /
  data bip2cs(21) /   -0.0000000000000452616E+00 /
  data bip2cs(22) /   -0.0000000000000030172E+00 /
  data bip2cs(23) /    0.0000000000000091058E+00 /
  data bip2cs(24) /   -0.0000000000000009814E+00 /
  data bip2cs(25) /   -0.0000000000000016429E+00 /
  data bip2cs(26) /    0.0000000000000005533E+00 /
  data bip2cs(27) /    0.0000000000000002175E+00 /
  data bip2cs(28) /   -0.0000000000000001737E+00 /
  data bip2cs(29) /   -0.0000000000000000010E+00 /

  data atr / 8.7506905708484345E+00 /
  data btr / -2.093836321356054E+00 /
  data nbif / 0 /
  data nbif2 / 0 /
  data nbig / 0 /
  data nbig2 / 0 /
  data nbip / 0 /
  data nbip2 / 0 /
  data x3sml / 0.0E+00 /
  data x32sml / 0.0E+00 /
  data xbig / 0.0E+00 /

  if ( nbif == 0 ) then
    eta = 0.1E+00 * r4_mach ( 3 )
    nbif = r4_inits ( bifcs, 9, eta )
    nbig = r4_inits ( bigcs, 8, eta )
    nbif2 = r4_inits ( bif2cs, 10, eta )
    nbig2 = r4_inits ( big2cs, 10, eta )
    nbip  = r4_inits ( bipcs, 24, eta )
    nbip2 = r4_inits ( bip2cs, 29, eta )
    x3sml = eta ** 0.3333E+00
    x32sml = 1.3104E+00 * x3sml * x3sml
    xbig = r4_mach ( 2 ) ** 0.6666E+00
  end if

  if ( x < -1.0E+00 ) then
    call r4_aimp ( x, xm, theta )
    r4_bie = xm * sin ( theta )
  else if ( abs ( x ) <= x32sml ) then
    z = 0.0E+00
    r4_bie = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) &
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) )
  else if ( abs ( x ) <= x3sml ) then
    z = 0.0E+00
    r4_bie = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) &
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) )
    r4_bie = r4_bie * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  else if ( x <= 1.0E+00 ) then
    z = x * x * x
    r4_bie = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) &
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) )
    r4_bie = r4_bie * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 )
  else if ( x <= 2.0E+00 ) then
    z = ( 2.0E+00 * x * x * x - 9.0E+00 ) / 7.0E+00
    r4_bie = exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 ) &
      * ( 1.125E+00 + r4_csevl ( z, bif2cs, nbif2 )&
      + x * ( 0.625E+00 + r4_csevl ( z, big2cs, nbig2 ) ) )
  else if ( x <= 4.0E+00 ) then
    sqrtx = sqrt ( x )
    z = atr / ( x * sqrtx ) + btr
    r4_bie = ( 0.625E+00 + r4_csevl ( z, bipcs, nbip ) ) / sqrt ( sqrtx )
  else if ( x <= xbig ) then
    sqrtx = sqrt ( x )
    z = 16.0E+00 / ( x * sqrtx ) - 1.0E+00
    r4_bie = ( 0.625E+00 + r4_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx )
  else
    sqrtx = sqrt ( x )
    z = - 1.0E+00
    r4_bie = ( 0.625E+00 + r4_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx )
  end if

  return
end
function r4_binom ( n, m )

!*****************************************************************************80
!
!! R4_BINOM evaluates the binomial coefficient using R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, M, the arguments.
!
!    Output, real ( kind = 4 ) R4_BINOM, the binomial coefficient.
!
  implicit none

  real ( kind = 4 ) bilnmx
  real ( kind = 4 ) corr
  real ( kind = 4 ) fintmx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ) r4_binom
  real ( kind = 4 ) r4_lgmc
  real ( kind = 4 ) r4_lnrel
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sq2pil
  real ( kind = 4 ) xk
  real ( kind = 4 ) xn
  real ( kind = 4 ) xnk

  save sq2pil
  save bilnmx
  save fintmx

  data sq2pil / 0.91893853320467274E+00 /
  data bilnmx / 0.0E+00 /
  data fintmx / 0.0E+00 /

  if ( bilnmx == 0.0E+00 ) then
    bilnmx = log ( r4_mach ( 2 ) )
    fintmx = 0.9E+00 / r4_mach ( 3 )
  end if

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BINOM - Fatal error!'
    write ( *, '(a)' ) '  N < 0.'
    stop
  end if

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BINOM - Fatal error!'
    write ( *, '(a)' ) '  M < 0.'
    stop
  end if

  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_BINOM - Fatal error!'
    write ( *, '(a)' ) '  N < M.'
    stop
  end if

  k = min ( m, n - m )

  if ( k <= 20 .and. real ( k, kind = 4 ) &
    * log ( real ( max ( n, 1 ), kind = 4 ) ) <= bilnmx ) then

    r4_binom = 1.0E+00

    do i = 1, k
      r4_binom = r4_binom * real ( n - i + 1, kind = 4 ) / real ( i, kind = 4 )
    end do

  else

    if ( k < 9 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_BINOM - Fatal error!'
      write ( *, '(a)' ) '  Result overflows.'
      write ( *, '(a)' ) '  N or M is too big.'
      stop
    end if

    xn = real ( n + 1, kind = 4 )
    xk = real ( k + 1, kind = 4 )
    xnk = real ( n - k + 1, kind = 4 )

    corr = r4_lgmc ( xn ) - r4_lgmc ( xk ) - r4_lgmc ( xnk )

    r4_binom = xk * log ( xnk / xk ) &
      - xn * r4_lnrel ( - ( xk - 1.0E+00 ) / xn ) &
      - 0.5E+00 * log ( xn * xnk / xk ) + 1.0E+00 - sq2pil + corr

    if ( bilnmx < r4_binom ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_BINOM - Fatal error!'
      write ( *, '(a)' ) '  Result overflows.'
      write ( *, '(a)' ) '  N or M is too big.'
      stop
    end if

    r4_binom = exp ( r4_binom )

  end if

  if ( r4_binom < fintmx ) then
    r4_binom = aint ( r4_binom + 0.5E+00 )
  end if

  return
end
function r4_cbrt ( x )

!*****************************************************************************80
!
!! R4_CBRT computes the cube root of an R4.
!
!  Discussion:
!
!    The approximation is a generalized Chebyshev series converted
!    to polynomial form.  The approximation is nearly best in the 
!    sense of relative error with 4.085 digits accuracy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the number whose square root is desired.
!
!    Output, real ( kind = 4 ) R4_CBRT, the cube root of X.
!
  implicit none

  real ( kind = 4 ) cbrt2(5)
  integer ( kind = 4 ) irem
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) ixpnt
  integer ( kind = 4 ) n
  integer ( kind = 4 ) niter
  real ( kind = 4 ) r4_cbrt
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_pak
  real ( kind = 4 ) value
  real ( kind = 4 ) vsq
  real ( kind = 4 ) x
  real ( kind = 4 ) y

  save cbrt2
  save niter

  data cbrt2(1) / 0.62996052494743658E+00 /
  data cbrt2(2) / 0.79370052598409974E+00 /
  data cbrt2(3) / 1.0E+00 /
  data cbrt2(4) / 1.25992104989487316E+00 /
  data cbrt2(5) / 1.58740105196819947E+00 /

  data niter / 0 /

  if ( niter == 0 ) then
    niter = int ( 1.443E+00 * log ( -0.106E+00 &
      * log ( 0.1E+00 * r4_mach ( 3 ) ) ) + 1.0E+00 )
  end if

  value = 0.0E+00

  if ( x /= 0.0E+00 ) then

    call r4_upak ( abs ( x ), y, n )
    ixpnt = n / 3
    irem = n - 3 * ixpnt + 3

    value = 0.439581E+00 + y * ( &
            0.928549E+00 + y * ( &
          - 0.512653E+00 + y * &
            0.144586E+00 ) )

    do iter = 1, niter
      vsq = value * value
      value = value + ( y - value * vsq ) / ( 3.0E+00 * vsq )
    end do

    if ( x < 0.0E+00 ) then
      value = - abs ( value )
    else
      value = + abs ( value )
    end if

    value = r4_pak ( cbrt2(irem) * value, ixpnt )

  end if

  r4_cbrt = value

  return
end
function r4_chi ( x )

!*****************************************************************************80
!
!! R4_CHI evaluates the hyperbolic cosine integral of an R4 argument.
!
!  Discussion:
!
!    The hyperbolic cosine integral is defined by
!
!      CHI(X) = gamma + log ( x ) 
!        + integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
!
!    where gamma is Euler's constant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_CHI, the hyperbolic cosine integral
!    evaluated at X.
!
  implicit none

  real ( kind = 4 ) r4_chi
  real ( kind = 4 ) r4_e1
  real ( kind = 4 ) r4_ei
  real ( kind = 4 ) x

  r4_chi = 0.5E+00 * ( r4_ei ( x ) - r4_e1 ( x ) )

  return
end
function r4_chu ( a, b, x )

!*****************************************************************************80
!
!! R4_CHU evaluates the confluent hypergeometric function of R4 arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the parameters.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_CHU, the function value.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) a0
  real ( kind = 4 ) aintb
  real ( kind = 4 ) alnx
  real ( kind = 4 ) b
  real ( kind = 4 ) b0
  real ( kind = 4 ) beps
  real ( kind = 4 ) c0
  real ( kind = 4 ), save :: eps = 0.0E+00
  real ( kind = 4 ) factor
  real ( kind = 4 ) gamri1
  real ( kind = 4 ) gamrni
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istrt
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), save :: pi = 3.14159265358979324E+00
  real ( kind = 4 ) pch1ai
  real ( kind = 4 ) pch1i
  real ( kind = 4 ) pochai
  real ( kind = 4 ) r4_chu
  real ( kind = 4 ) r4_chu_scaled
  real ( kind = 4 ) r4_exprel
  real ( kind = 4 ) r4_gamma
  real ( kind = 4 ) r4_gamr
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_mop
  real ( kind = 4 ) r4_poch
  real ( kind = 4 ) r4_poch1
  real ( kind = 4 ) sum
  real ( kind = 4 ) t
  real ( kind = 4 ) x
  real ( kind = 4 ) xeps1
  real ( kind = 4 ) xi
  real ( kind = 4 ) xi1
  real ( kind = 4 ) xn
  real ( kind = 4 ) xtoeps

  real ( kind = 4 ) temp

  if ( eps == 0.0E+00 ) then
    eps = r4_mach ( 3 )
  end if

  if ( x < 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_CHU - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop
  end if

  if ( x == 0.0E+00 ) then
    if ( 1.0E+00 <= b ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_CHU - Fatal error!'
      write ( *, '(a)' ) '  X = 0 and 1 <= B.'
      stop
    end if
    r4_chu = r4_gamma ( 1.0E+00 - b ) / r4_gamma ( 1.0E+00 + a - b )
    return
  end if

  if ( max ( abs ( a ), 1.0E+00 ) * max ( abs ( 1.0E+00 + a - b ), 1.0E+00 ) &
    < 0.99E+00 * abs ( x ) ) then
    r4_chu = x ** ( - a ) * r4_chu_scaled ( a, b, x )
    return
  end if
!
!  The ascending series will be used, because the descending rational
!  approximation (which is based on the asymptotic series) is unstable.
!
  if ( b < 0.0E+00 ) then
    aintb = aint ( b - 0.5E+00 )
  else
    aintb = aint ( b + 0.5E+00 )
  end if
  beps = b - aintb
  n = aintb

  alnx = log ( x )
  xtoeps = exp ( - beps * alnx )
!
!  Evaluate the finite sum.
!
!  Consider the case b < 1.0 first.
!
  if ( n < 1 ) then

    sum = 1.0E+00
    t = 1.0E+00
    m = - n
    do i = 1, m
      xi1 = real ( i - 1, kind = 4 )
      t = t * ( a + xi1 ) * x / ( ( b + xi1 ) * ( xi1 + 1.0E+00 ) )
      sum = sum + t
    end do

    sum = r4_poch ( 1.0E+00 + a - b, - a ) * sum
!
!  Now consider the case b .ge. 1.0.
!
  else

    sum = 0.0E+00
    m = n - 2

    if ( 0 <= m ) then

      t = 1.0E+00
      sum = 1.0E+00

      do i = 1, m
        xi = real ( i, kind = 4 )
        t = t * ( a - b + xi ) * x / ( ( 1.0E+00 - b + xi ) * xi )
        sum = sum + t
      end do

      sum = r4_gamma ( b - 1.0E+00 ) * r4_gamr ( a ) &
        * x ** ( 1 - n ) * xtoeps * sum

    end if

  end if
!
!  Now evaluate the infinite sum.
!
  if ( n < 1 ) then
    istrt = 1 - n
  else
    istrt = 0
  end if

  xi = real ( istrt, kind = 4 )

  factor = r4_mop ( n ) * r4_gamr ( 1.0E+00 + a - b ) * x ** istrt

  if ( beps /= 0.0E+00 ) then
    factor = factor * beps * pi / sin ( beps * pi )
  end if

  pochai = r4_poch ( a, xi )
  gamri1 = r4_gamr ( xi + 1.0E+00 )
  gamrni = r4_gamr ( aintb + xi )
  b0 = factor * r4_poch ( a, xi - beps ) * gamrni &
    * r4_gamr ( xi + 1.0E+00 - beps )
!
!  x^(-beps) is close to 1.0, so we must be careful in evaluating
!  the differences.
!
  if ( abs ( xtoeps - 1.0E+00 ) <= 0.5E+00 ) then

    pch1ai = r4_poch1 ( a + xi, - beps )
    pch1i = r4_poch1 ( xi + 1.0E+00 - beps, beps )
    c0 = factor * pochai * gamrni * gamri1 * ( &
      - r4_poch1 ( b + xi, -beps ) + pch1ai &
      - pch1i + beps * pch1ai * pch1i )
!
!  xeps1 = (1.0 - x^(-beps)) / beps
!
    xeps1 = alnx * r4_exprel ( - beps * alnx )

    r4_chu = sum + c0 + xeps1 * b0
    xn = real ( n, kind = 4 )

    do i = 1, 1000
      xi = real ( istrt + i, kind = 4 )
      xi1 = real ( istrt + i - 1, kind = 4 )
      b0 = ( a + xi1 - beps ) * b0 * x &
        / ( ( xn + xi1 ) * ( xi - beps ) )
      c0 = ( a + xi1 ) * c0 * x / ( ( b + xi1 ) * xi ) &
        - ( ( a - 1.0E+00 ) * ( xn + 2.0E+00 * xi - 1.0E+00 )&
        + xi * ( xi - beps ) ) * b0 &
        / ( xi * ( b + xi1 ) * ( a + xi1 - beps ) )
      t = c0 + xeps1 * b0
      r4_chu = r4_chu + t
      if ( abs ( t ) < eps * abs ( r4_chu ) ) then
        return
      end if
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_CHU - Fatal error!'
    write ( *, '(a)' ) '  No convergence in 1000 terms.'
    stop

  end if
!
!  x^(-beps) is very different from 1.0, so the straightforward
!  formulation is stable.
!
  a0 = factor * pochai * r4_gamr ( b + xi ) * gamri1 / beps
  b0 = xtoeps * b0 / beps

  r4_chu = sum + a0 - b0

  do i = 1, 1000
    xi = real ( istrt + i, kind = 4 )
    xi1 = real ( istrt + i - 1, kind = 4 )
    a0 = ( a + xi1 ) * a0 * x / ( ( b + xi1 ) * xi )
    b0 = ( a + xi1 - beps ) * b0 * x &
      / ( ( aintb + xi1 ) * ( xi - beps ) )
    t = a0 - b0
    r4_chu = r4_chu + t
    if ( abs ( t ) < eps * abs ( r4_chu ) ) then
      return
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R4_CHU - Fatal error!'
  write ( *, '(a)' ) '  No convergence in 1000 terms.'

  stop
end
function r4_chu_scaled ( a, b, z )

!*****************************************************************************80
!
!! R4_CHU_SCALED: scaled confluent hypergeometric function of R4 arguments.
!
!  Discussion:
!
!    Evaluate, for large z, z**a * u(a,b,z)  where U is the logarithmic
!    confluent hypergeometric function.  A rational approximation due to
!    Y L Luke is used.  When U is not in the asymptotic region, that is, when A
!    or B is large compared with Z, considerable significance loss occurs.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the parameters.
!
!    Input, real ( kind = 4 ) Z, the argument.
!
!    Output, real ( kind = 4 ) R4_CHU_SCALED, the function value.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) aa(4)
  real ( kind = 4 ) ab
  real ( kind = 4 ) anbn
  real ( kind = 4 ) b
  real ( kind = 4 ) bb(4)
  real ( kind = 4 ) bp
  real ( kind = 4 ) c2
  real ( kind = 4 ) ct1
  real ( kind = 4 ) ct2
  real ( kind = 4 ) ct3
  real ( kind = 4 ) d1z
  real ( kind = 4 ) eps
  real ( kind = 4 ) g1
  real ( kind = 4 ) g2
  real ( kind = 4 ) g3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_chu_scaled
  real ( kind = 4 ) sab
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) x2i1
  real ( kind = 4 ) z

  save eps

  data eps / 0.0E+00 /

  if ( eps == 0.0E+00 ) then
    eps = 4.0E+00 * r4_mach ( 4 )
    sqeps = sqrt ( r4_mach ( 4 ) )
  end if

  bp = 1.0E+00 + a - b
  ab = a * bp
  ct2 = 2.0E+00 * ( z - ab )
  sab = a + bp

  bb(1) = 1.0E+00
  aa(1) = 1.0E+00

  ct3 = sab + 1.0E+00 + ab
  bb(2) = 1.0E+00 + 2.0E+00 * z / ct3
  aa(2) = 1.0E+00 + ct2 / ct3

  anbn = ct3 + sab + 3.0E+00
  ct1 = 1.0E+00 + 2.0E+00 * z / anbn
  bb(3) = 1.0E+00 + 6.0E+00 * ct1 * z / ct3
  aa(3) = 1.0E+00 + 6.0E+00 * ab / anbn + 3.0E+00 * ct1 * ct2 / ct3

  do i = 4, 300

    x2i1 = real ( 2 * i - 3, kind = 4 )
    ct1 = x2i1 / ( x2i1 - 2.0 )
    anbn = anbn + x2i1 + sab
    ct2 = ( x2i1 - 1.0E+00 ) / anbn
    c2 = x2i1 * ct2 - 1.0E+00
    d1z = x2i1 * 2.0E+00 * z / anbn

    ct3 = sab * ct2
    g1 = d1z + ct1 * ( c2 + ct3 )
    g2 = d1z - c2
    g3 = ct1 * ( 1.0E+00 - ct3 - 2.0E+00 * ct2 )

    bb(4) = g1 * bb(3) + g2 * bb(2) + g3 * bb(1)
    aa(4) = g1 * aa(3) + g2 * aa(2) + g3 * aa(1)

    r4_chu_scaled = aa(4) / bb(4)

    if ( abs ( r4_chu_scaled - aa(1) / bb(1) ) &
      < eps * abs ( r4_chu_scaled ) ) then
      return
    end if
!
!  If overflows or underflows prove to be a problem, the statements
!  below could be altered to incorporate a dynamically adjusted scale
!  factor.
!
    do j = 1, 3
      bb(j) = bb(j+1)
      aa(j) = aa(j+1)
    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R4_CHU_SCALED - Fatal error!'
  write ( *, '(a)' ) '  No convergence in 300 terms.'

  stop
end
function r4_ci ( x )

!*****************************************************************************80
!
!! R4_CI evaluates the cosine integral Ci of an R4 argument.
!
!  Discussion:
!
!    The cosine integral is defined by
!
!      CI(X) = - integral ( X <= T < Infinity ) ( cos ( T ) ) / T  dT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_CI, the cosine integral Ci evaluated at X.
!
  implicit none

  real ( kind = 4 ) cics(13)
  real ( kind = 4 ) f
  real ( kind = 4 ) g
  integer ( kind = 4 ) nci
  real ( kind = 4 ) r4_ci
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sinx
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) xsml

  save cics
  save nci
  save xsml

  data cics(  1) /    -0.34004281856055363156E+00 /
  data cics(  2) /    -1.03302166401177456807E+00 /
  data cics(  3) /     0.19388222659917082877E+00 /
  data cics(  4) /    -0.01918260436019865894E+00 /
  data cics(  5) /     0.00110789252584784967E+00 /
  data cics(  6) /    -0.00004157234558247209E+00 /
  data cics(  7) /     0.00000109278524300229E+00 /
  data cics(  8) /    -0.00000002123285954183E+00 /
  data cics(  9) /     0.00000000031733482164E+00 /
  data cics( 10) /    -0.00000000000376141548E+00 /
  data cics( 11) /     0.00000000000003622653E+00 /
  data cics( 12) /    -0.00000000000000028912E+00 /
  data cics( 13) /     0.00000000000000000194E+00 /

  data nci / 0 /
  data xsml / 0.0E+00 /

  if ( nci == 0 ) then
    nci = r4_inits ( cics, 13, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( r4_mach ( 3 ) )
  end if

  if ( x <= 0.0E+00 ) then
    
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_CI - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.0.'
    stop

  else if ( x <= xsml ) then

    y = - 1.0E+00
    r4_ci = log ( x ) - 0.5E+00 + r4_csevl ( y, cics, nci )

  else if ( x <= 4.0E+00 ) then

    y = ( x * x - 8.0E+00 ) * 0.125E+00
    r4_ci = log ( x ) - 0.5E+00 + r4_csevl ( y, cics, nci )

  else

    call r4_sifg ( x, f, g )
    sinx = sin ( x )
    r4_ci = f * sinx - g * cos ( x )

  end if

  return
end
function r4_cin ( x )

!*****************************************************************************80
!
!! R4_CIN evaluates the alternate cosine integral Cin of an R4 argument.
!
!  Discussion:
!
!    CIN(X) = gamma + log(X) 
!      + integral ( 0 <= T <= X ) ( cos ( T ) - 1 ) / T  dT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_CIN, the cosine integral Cin evaluated at X.
!
  implicit none

  real ( kind = 4 ) absx
  real ( kind = 4 ) cincs(11)
  real ( kind = 4 ) eul
  real ( kind = 4 ) f
  real ( kind = 4 ) g
  integer ( kind = 4 ) ncin
  real ( kind = 4 ) r4_cin
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sinx
  real ( kind = 4 ) x
  real ( kind = 4 ) xmin

  save cincs
  save eul
  save ncin
  save xmin

  data cincs(  1) /     0.3707450175090968874E+00 /
  data cincs(  2) /    -0.0589357489636444683E+00 /
  data cincs(  3) /     0.0053818964211356912E+00 /
  data cincs(  4) /    -0.0002986005284196214E+00 /
  data cincs(  5) /     0.0000109557257532162E+00 /
  data cincs(  6) /    -0.0000002840545487735E+00 /
  data cincs(  7) /     0.0000000054697399488E+00 /
  data cincs(  8) /    -0.0000000000812418746E+00 /
  data cincs(  9) /     0.0000000000009586859E+00 /
  data cincs( 10) /    -0.0000000000000092027E+00 /
  data cincs( 11) /     0.0000000000000000733E+00 /

  data eul / 0.57721566490153286E+00 /
  data ncin / 0 /
  data xmin / 0.0E+00 /

  if ( ncin == 0 ) then
    ncin = r4_inits ( cincs, 11, 0.1E+00 * r4_mach ( 3 ) )
    xmin = sqrt ( r4_mach ( 1 ) )
  end if

  absx = abs ( x )

  if ( absx <= xmin ) then

    r4_cin = 0.0E+00

  else if ( x <= 4.0E+00 ) then

    r4_cin = x * x * r4_csevl ( ( x * x - 8.0E+00 ) * 0.125E+00, &
      cincs, ncin )

  else

    call r4_sifg ( x, f, g )
    sinx = sin ( absx )
    r4_cin = - f * sinx + g * cos ( absx ) + log ( absx ) + eul

  end if

  return
end
function r4_cinh ( x )

!*****************************************************************************80
!
!! R4_CINH: alternate hyperbolic cosine integral Cinh of an R4 argument.
!
!  Discussion:
!
!    Cinh ( x ) = Integral ( 0 <= t <= x ) ( cosh ( t ) - 1 ) dt / t
!
!    The original text of this program had a mistake:
!      y = x * x / 9.0E+00 - 1.0E+00
!    has been corrected to
!      y = x * x / 4.5E+00 - 1.0E+00
!    JVB, 27 March 2010
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_CINH, the hyperbolic cosine integral Cinh
!    evaluated at X.
!
  implicit none

  real ( kind = 4 ) absx
  real ( kind = 4 ) cinhcs(10)
  real ( kind = 4 ) eul
  integer ( kind = 4 ) ncinh
  real ( kind = 4 ) r4_chi
  real ( kind = 4 ) r4_cinh
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save cinhcs
  save eul
  save ncinh
  save xmin
  save xsml

  data cinhcs(  1) /     0.1093291636520734431E+00 /
  data cinhcs(  2) /     0.0573928847550379676E+00 /
  data cinhcs(  3) /     0.0028095756978830353E+00 /
  data cinhcs(  4) /     0.0000828780840721357E+00 /
  data cinhcs(  5) /     0.0000016278596173914E+00 /
  data cinhcs(  6) /     0.0000000227809519256E+00 /
  data cinhcs(  7) /     0.0000000002384484842E+00 /
  data cinhcs(  8) /     0.0000000000019360830E+00 /
  data cinhcs(  9) /     0.0000000000000125454E+00 /
  data cinhcs( 10) /     0.0000000000000000664E+00 /

  data eul / 0.57721566490153286E+00 /
  data ncinh / 0 /
  data xmin / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( ncinh == 0 ) then
    ncinh = r4_inits ( cinhcs, 10, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( r4_mach ( 3 ) )
    xmin = 2.0E+00 * sqrt ( r4_mach ( 1 ) )
  end if

  absx = abs ( x )

  if ( x == 0.0E+00 ) then
    r4_cinh = 0.0E+00
  else if ( absx <= xmin ) then
    r4_cinh = 0.0E+00
  else if ( x <= xsml ) then
    y = - 1.0E+00
    r4_cinh = x * x * ( 0.25E+00 + r4_csevl ( y, cinhcs, ncinh ) )
  else if ( x <= 3.0E+00 ) then
    y = x * x / 4.5E+00 - 1.0E+00
    r4_cinh = x * x * ( 0.25E+00 + r4_csevl ( y, cinhcs, ncinh ) )
  else
    r4_cinh = r4_chi ( absx ) - eul - log ( absx )
  end if

  return
end
function r4_cos ( x )

!*****************************************************************************80
!
!! R4_COS evaluates the cosine of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_COS, the cosine of X.
!
  implicit none

  real ( kind = 4 ) absx
  real ( kind = 4 ) f
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) ntsn
  real ( kind = 4 ), parameter :: pi2 = 1.57079632679489661923E+00
  real ( kind = 4 ) pi2rec
  real ( kind = 4 ) pihi
  real ( kind = 4 ) pilo
  real ( kind = 4 ) pirec
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_cos
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) sincs(10)
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xn
  real ( kind = 4 ) xsml
  real ( kind = 4 ) xwarn
  real ( kind = 4 ) y

  save ntsn
  save pi2rec
  save pihi
  save pilo
  save pirec
  save sincs
  save xmax
  save xsml
  save xwarn

  data sincs(  1) / -0.374991154955873175840E+00 /
  data sincs(  2) / -0.181603155237250201864E+00 /
  data sincs(  3) / +0.005804709274598633559E+00 /
  data sincs(  4) / -0.000086954311779340757E+00 /
  data sincs(  5) / +0.000000754370148088851E+00 /
  data sincs(  6) / -0.000000004267129665056E+00 /
  data sincs(  7) / +0.000000000016980422945E+00 /
  data sincs(  8) / -0.000000000000050120579E+00 /
  data sincs(  9) / +0.000000000000000114101E+00 /
  data sincs( 10) / -0.000000000000000000206E+00 /
!
!  pihi + pilo = pi.  pihi is exactly representable on all machines
!  with at least 8 bits of precision.  whether it is exactly
!  represented depends on the compiler.  this routine is more
!  accurate if it is exactly represented.
!
  data ntsn / 0 /
  data pi2rec / 0.636619772367581343E+00 /
  data pihi / 3.140625E+00 /
  data pilo / 9.6765358979323846E-04 /
  data pirec / 0.31830988618379067E+00 /
  data xmax / 0.0E+00 /
  data xsml / 0.0E+00 /
  data xwarn / 0.0E+00 /

  if ( ntsn == 0 ) then
    ntsn = r4_inits ( sincs, 10, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( 2.0E+00 * r4_mach ( 3 ) )
    xmax = 1.0E+00 / r4_mach ( 4 )
    xwarn = sqrt ( xmax )
  end if
  
  absx = abs ( x )
  y = absx + pi2

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_COS - Warning!'
    write ( *, '(a)' ) '  No precision because |X| is big.'
    r4_cos = 0.0E+00
    return
  end if

  if ( xwarn < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_COS - Warning!'
    write ( *, '(a)' ) '  Answer < half precision because |X| is big.'
  end if

  r4_cos = 1.0E+00

  if ( absx < xsml ) then
    return
  end if

  xn = aint ( y * pirec + 0.5E+00 )
  n2 = int ( mod ( xn, 2.0E+00 ) + 0.5E+00 )
  xn = xn - 0.5E+00
  f = ( absx - xn * pihi ) - xn * pilo

  xn = 2.0E+00 * ( f * pi2rec ) * ( f * pi2rec ) - 1.0E+00
  r4_cos = f + f * r4_csevl ( xn, sincs, ntsn )

  if ( n2 /= 0 ) then
    r4_cos = - r4_cos
  end if

  if ( r4_cos < - 1.0E+00 ) then
    r4_cos = - 1.0E+00
  else if ( 1.0E+00 < r4_cos ) then
    r4_cos = + 1.0E+00
  end if

  return
end
function r4_cos_deg ( x )

!*****************************************************************************80
!
!! R4_COS_DEG evaluates the cosine of an R4 argument in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument in degrees.
!
!    Output, real ( kind = 4 ) R4_COS_DEG, the cosine of X.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 4 ) r4_cos_deg
  real ( kind = 4 ), parameter :: raddeg = 0.017453292519943296E+00
  real ( kind = 4 ) value
  real ( kind = 4 ) x

  value = cos ( raddeg * x )

  if ( mod ( x, 90.0E+00 ) == 0.0E+00 ) then

    n = int ( abs ( x ) / 90.0E+00 + 0.5E+00 )
    n = mod ( n, 2 )

    if ( n == 1 ) then
      value = 0.0E+00
    else if ( value < 0.0E+00 ) then
      value = - 1.0E+00
    else
      value = + 1.0E+00
    end if

  end if

  r4_cos_deg = value

  return
end
function r4_cosh ( x )

!*****************************************************************************80
!
!! R4_COSH evaluates the hyperbolic cosine of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_COSH, the hyperbolic cosine of X.
!
  implicit none

  real ( kind = 4 ) r4_cosh
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) ymax

  save ymax

  data ymax / 0.0E+00 /

  if ( ymax == 0.0E+00 ) then 
    ymax = 1.0E+00 / sqrt ( r4_mach ( 3 ) )
  end if

  y = exp ( abs ( x ) )

  value = 0.5E+00 * y

  if ( y < ymax ) then
    value = 0.5E+00 * ( y + 1.0E+00 / y )
  end if

  r4_cosh = value

  return
end
function r4_cot ( x )

!*****************************************************************************80
!
!! R4_COT evaluates the cotangent of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_COT, the cotangent of X.
!
  implicit none

  real ( kind = 4 ) ainty
  real ( kind = 4 ) ainty2
  real ( kind = 4 ) cotcs(8)
  integer ( kind = 4 ) ifn
  integer ( kind = 4 ) nterms
  real ( kind = 4 ) pi2rec
  real ( kind = 4 ) prodbg
  real ( kind = 4 ) r4_cot
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y
  real ( kind = 4 ) yrem

  save cotcs
  save nterms
  save pi2rec
  save xmax
  save xmin
  save xsml

  data cotcs( 1) /    0.24025916098295630E+00 /
  data cotcs( 2) /   -0.016533031601500228E+00 /
  data cotcs( 3) /   -0.000042998391931724E+00 /
  data cotcs( 4) /   -0.000000159283223327E+00 /
  data cotcs( 5) /   -0.000000000619109313E+00 /
  data cotcs( 6) /   -0.000000000002430197E+00 /
  data cotcs( 7) /   -0.000000000000009560E+00 /
  data cotcs( 8) /   -0.000000000000000037E+00 /

  data nterms / 0 /
  data pi2rec / 0.0116197723675813430E+00 /
  data xmax / 0.0E+00 /
  data xmin / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( cotcs, 8, 0.1E+00 * r4_mach ( 3 ) )
    xmax = 1.0E+00 / r4_mach ( 4 )
    xsml = sqrt ( 3.0E+00 * r4_mach ( 3 ) )
    xmin = exp ( max ( log ( r4_mach ( 1 ) ), &
      - log ( r4_mach ( 2 ) ) )  + 0.01E+00 )
    sqeps = sqrt ( r4_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( y < xmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_COT - Fatal error!'
    write ( *, '(a)' ) '  |X| is too small.'
    stop
  end if

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_COT - Fatal error!'
    write ( *, '(a)' ) '  |X| is too big.'
    stop
  end if
!
!  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
!  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
!  = aint(.625*y) + aint(z) + rem(z)
!
  ainty = aint ( y )
  yrem = y - ainty
  prodbg = 0.625E+00 * ainty
  ainty = aint ( prodbg )
  y = ( prodbg - ainty ) + 0.625E+00 * yrem + y * pi2rec
  ainty2 = aint ( y )
  ainty = ainty + ainty2
  y = y - ainty2

  ifn = int ( mod ( ainty, 2.0E+00 ) )

  if ( ifn == 1 ) then
    y = 1.0E+00 - y
  end if

  if ( 0.5E+00 < abs ( x ) .and. y < abs ( x ) * sqeps ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_COT - Warning!'
    write ( *, '(a)' ) '  Answer less than half precision.'
    write ( *, '(a)' ) &
      '  |X| too big, or X nearly a nonzero multiple of pi.'
  end if

  if ( y == 0.0E+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_COT - Fatal error!'
    write ( *, '(a)' ) '  X is a multiple of pi.'
    stop

  else if ( y <= xsml ) then

    r4_cot = 1.0E+00 / y 

  else if ( y <= 0.25E+00 ) then

    r4_cot = ( 0.5E+00 &
      + r4_csevl ( 32.0E+00 * y * y - 1.0E+00, cotcs, nterms ) ) / y

  else if ( y <= 0.5E+00 ) then

    r4_cot = ( 0.5E+00 + r4_csevl ( 8.0E+00 * y * y - 1.0E+00, &
      cotcs, nterms ) ) / ( 0.5E+00 * y )

    r4_cot = ( r4_cot * r4_cot - 1.0E+00 ) * 0.5E+00 / r4_cot

  else

    r4_cot = ( 0.5E+00 + r4_csevl ( 2.0E+00 * y * y - 1.0E+00, &
      cotcs, nterms ) ) / ( 0.25E+00 * y )
    r4_cot = ( r4_cot * r4_cot - 1.0E+00 ) * 0.5E+00 / r4_cot
    r4_cot = ( r4_cot * r4_cot - 1.0E+00 ) * 0.5E+00 / r4_cot

  end if

  if ( x < 0.0E+00 ) then
    r4_cot = - abs ( r4_cot )
  else
    r4_cot = + abs ( r4_cot )
  end if

  if ( ifn == 1 ) then
    r4_cot = - r4_cot
  end if

  return
end
function r4_csevl ( x, cs, n )

!*****************************************************************************80
!
!! R4_CSEVL evaluates a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the evaluation point.
!
!    Input, real ( kind = 4 ) CS(N), the Chebyshev coefficients.
!
!    Input, integer ( kind = 4 ) N, the number of Chebyshev coefficients.
!
!    Output, real ( kind = 4 ) R4_CSEVL, the Chebyshev series evaluated at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) b0
  real ( kind = 4 ) b1
  real ( kind = 4 ) b2
  real ( kind = 4 ) cs(n)
  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_csevl
  real ( kind = 4 ) twox
  real ( kind = 4 ) x

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms <= 0.'
    stop
  end if

  if ( 1000 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms > 1000.'
    stop
  end if

  if ( x < -1.1E+00 .or. 1.0E+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_CSEVL - Fatal error!'
    write ( *, '(a)' ) '  X outside (-1,+1).'
    write ( *, '(a,g14.6)' ) '  X = ', x
    stop
  end if

  b1 = 0.0E+00
  b0 = 0.0E+00
  twox = 2.0E+00 * x

  do i = n, 1, -1
    b2 = b1
    b1 = b0
    b0 = twox * b1 - b2 + cs(i)
  end do

  r4_csevl = 0.5E+00 * ( b0 - b2 )

  return
end
function r4_dawson ( x )

!*****************************************************************************80
!
!! R4_DAWSON evaluates Dawson's integral of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_DAWSON, the value of Dawson's integral at X.
!
  implicit none

  real ( kind = 4 ) daw2cs(29)
  real ( kind = 4 ) dawacs(26)
  real ( kind = 4 ) dawcs(13)
  real ( kind = 4 ) eps
  integer ( kind = 4 ) ntdaw
  integer ( kind = 4 ) ntdaw2
  integer ( kind = 4 ) ntdawa
  real ( kind = 4 ) r4_csevl
  real ( kind = 4 ) r4_dawson
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xbig
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save daw2cs
  save dawacs
  save dawcs
  save ntdaw
  save ntdaw2
  save ntdawa
  save xbig
  save xmax
  save xsml

  data dawcs( 1) /   -0.006351734375145949E+00 /
  data dawcs( 2) /   -0.22940714796773869E+00 /
  data dawcs( 3) /    0.022130500939084764E+00 /
  data dawcs( 4) /   -0.001549265453892985E+00 /
  data dawcs( 5) /    0.000084973277156849E+00 /
  data dawcs( 6) /   -0.000003828266270972E+00 /
  data dawcs( 7) /    0.000000146285480625E+00 /
  data dawcs( 8) /   -0.000000004851982381E+00 /
  data dawcs( 9) /    0.000000000142146357E+00 /
  data dawcs(10) /   -0.000000000003728836E+00 /
  data dawcs(11) /    0.000000000000088549E+00 /
  data dawcs(12) /   -0.000000000000001920E+00 /
  data dawcs(13) /    0.000000000000000038E+00 /

  data daw2cs( 1) /   -0.056886544105215527E+00 /
  data daw2cs( 2) /   -0.31811346996168131E+00 /
  data daw2cs( 3) /    0.20873845413642237E+00 /
  data daw2cs( 4) /   -0.12475409913779131E+00 /
  data daw2cs( 5) /    0.067869305186676777E+00 /
  data daw2cs( 6) /   -0.033659144895270940E+00 /
  data daw2cs( 7) /    0.015260781271987972E+00 /
  data daw2cs( 8) /   -0.006348370962596214E+00 /
  data daw2cs( 9) /    0.002432674092074852E+00 /
  data daw2cs(10) /   -0.000862195414910650E+00 /
  data daw2cs(11) /    0.000283765733363216E+00 /
  data daw2cs(12) /   -0.000087057549874170E+00 /
  data daw2cs(13) /    0.000024986849985481E+00 /
  data daw2cs(14) /   -0.000006731928676416E+00 /
  data daw2cs(15) /    0.000001707857878557E+00 /
  data daw2cs(16) /   -0.000000409175512264E+00 /
  data daw2cs(17) /    0.000000092828292216E+00 /
  data daw2cs(18) /   -0.000000019991403610E+00 /
  data daw2cs(19) /    0.000000004096349064E+00 /
  data daw2cs(20) /   -0.000000000800324095E+00 /
  data daw2cs(21) /    0.000000000149385031E+00 /
  data daw2cs(22) /   -0.000000000026687999E+00 /
  data daw2cs(23) /    0.000000000004571221E+00 /
  data daw2cs(24) /   -0.000000000000751873E+00 /
  data daw2cs(25) /    0.000000000000118931E+00 /
  data daw2cs(26) /   -0.000000000000018116E+00 /
  data daw2cs(27) /    0.000000000000002661E+00 /
  data daw2cs(28) /   -0.000000000000000377E+00 /
  data daw2cs(29) /    0.000000000000000051E+00 /

  data dawacs( 1) /    0.01690485637765704E+00 /
  data dawacs( 2) /    0.00868325227840695E+00 /
  data dawacs( 3) /    0.00024248640424177E+00 /
  data dawacs( 4) /    0.00001261182399572E+00 /
  data dawacs( 5) /    0.00000106645331463E+00 /
  data dawacs( 6) /    0.00000013581597947E+00 /
  data dawacs( 7) /    0.00000002171042356E+00 /
  data dawacs( 8) /    0.00000000286701050E+00 /
  data dawacs( 9) /   -0.00000000019013363E+00 /
  data dawacs(10) /   -0.00000000030977804E+00 /
  data dawacs(11) /   -0.00000000010294148E+00 /
  data dawacs(12) /   -0.00000000000626035E+00 /
  data dawacs(13) /    0.00000000000856313E+00 /
  data dawacs(14) /    0.00000000000303304E+00 /
  data dawacs(15) /   -0.00000000000025236E+00 /
  data dawacs(16) /   -0.00000000000042106E+00 /
  data dawacs(17) /   -0.00000000000004431E+00 /
  data dawacs(18) /    0.00000000000004911E+00 /
  data dawacs(19) /    0.00000000000001235E+00 /
  data dawacs(20) /   -0.00000000000000578E+00 /
  data dawacs(21) /   -0.00000000000000228E+00 /
  data dawacs(22) /    0.00000000000000076E+00 /
  data dawacs(23) /    0.00000000000000038E+00 /
  data dawacs(24) /   -0.00000000000000011E+00 /
  data dawacs(25) /   -0.00000000000000006E+00 /
  data dawacs(26) /    0.00000000000000002E+00 /

  data ntdaw / 0 /
  data ntdaw2 / 0 /
  data ntdawa / 0 /
  data xbig / 0.0E+00 /
  data xmax / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( ntdaw == 0 ) then

    eps = r4_mach ( 3 )
    ntdaw  = r4_inits ( dawcs,  13, 0.1E+00 * eps )
    ntdaw2 = r4_inits ( daw2cs, 29, 0.1E+00 * eps )
    ntdawa = r4_inits ( dawacs, 26, 0.1E+00 * eps )

    xsml = sqrt ( 1.5E+00 * eps )
    xbig = sqrt ( 0.5E+00 / eps )
    xmax = exp ( min ( - log ( 2.0E+00 * r4_mach ( 1 ) ), &
      log ( r4_mach ( 2 ) ) ) - 0.01E+00 )

  end if

  y = abs ( x )

  if ( y <= xsml ) then
    r4_dawson = x
  else if ( y <= 1.0E+00 ) then
    r4_dawson = x * ( 0.75E+00 &
      + r4_csevl ( 2.0E+00 * y * y - 1.0E+00, dawcs, ntdaw ) )
  else if ( y <= 4.0E+00 ) then
    r4_dawson = x * ( 0.25E+00 &
      + r4_csevl ( 0.125E+00 * y * y - 1.0E+00, daw2cs, ntdaw2 ) )
  else if ( y < xbig ) then
    r4_dawson = ( 0.5E+00 &
      + r4_csevl ( 32.0E+00 / y / y - 1.0E+00, dawacs, ntdawa ) ) / x
  else if ( y <= xmax ) then
    r4_dawson = 0.5E+00 / x
  else
    r4_dawson = 0.0E+00
  end if

  return
end
function r4_e1 ( x )

!*****************************************************************************80
!
!! R4_E1 evaluates the exponential integral E1 for an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_E1, the exponential integral E1.
!
  implicit none

  real ( kind = 4 ) ae11cs(39)
  real ( kind = 4 ) ae12cs(25)
  real ( kind = 4 ) ae13cs(25)
  real ( kind = 4 ) ae14cs(26)
  real ( kind = 4 ) e11cs(19)
  real ( kind = 4 ) e12cs(16)
  real ( kind = 4 ) eta
  integer ( kind = 4 ) ntae11
  integer ( kind = 4 ) ntae12
  integer ( kind = 4 ) ntae13
  integer ( kind = 4 ) ntae14
  integer ( kind = 4 ) nte11
  integer ( kind = 4 ) nte12
  real ( kind = 4 ) r4_csevl
  real ( kind = 4 ) r4_e1
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax

  save ae11cs
  save ae12cs
  save ae13cs
  save ae14cs
  save e11cs
  save e12cs
  save ntae11
  save ntae12
  save ntae13
  save ntae14
  save nte11
  save nte12
  save xmax

  data ae11cs( 1) /    0.12150323971606579E+00 /
  data ae11cs( 2) /   -0.065088778513550150E+00 /
  data ae11cs( 3) /    0.004897651357459670E+00 /
  data ae11cs( 4) /   -0.000649237843027216E+00 /
  data ae11cs( 5) /    0.000093840434587471E+00 /
  data ae11cs( 6) /    0.000000420236380882E+00 /
  data ae11cs( 7) /   -0.000008113374735904E+00 /
  data ae11cs( 8) /    0.000002804247688663E+00 /
  data ae11cs( 9) /    0.000000056487164441E+00 /
  data ae11cs(10) /   -0.000000344809174450E+00 /
  data ae11cs(11) /    0.000000058209273578E+00 /
  data ae11cs(12) /    0.000000038711426349E+00 /
  data ae11cs(13) /   -0.000000012453235014E+00 /
  data ae11cs(14) /   -0.000000005118504888E+00 /
  data ae11cs(15) /    0.000000002148771527E+00 /
  data ae11cs(16) /    0.000000000868459898E+00 /
  data ae11cs(17) /   -0.000000000343650105E+00 /
  data ae11cs(18) /   -0.000000000179796603E+00 /
  data ae11cs(19) /    0.000000000047442060E+00 /
  data ae11cs(20) /    0.000000000040423282E+00 /
  data ae11cs(21) /   -0.000000000003543928E+00 /
  data ae11cs(22) /   -0.000000000008853444E+00 /
  data ae11cs(23) /   -0.000000000000960151E+00 /
  data ae11cs(24) /    0.000000000001692921E+00 /
  data ae11cs(25) /    0.000000000000607990E+00 /
  data ae11cs(26) /   -0.000000000000224338E+00 /
  data ae11cs(27) /   -0.000000000000200327E+00 /
  data ae11cs(28) /   -0.000000000000006246E+00 /
  data ae11cs(29) /    0.000000000000045571E+00 /
  data ae11cs(30) /    0.000000000000016383E+00 /
  data ae11cs(31) /   -0.000000000000005561E+00 /
  data ae11cs(32) /   -0.000000000000006074E+00 /
  data ae11cs(33) /   -0.000000000000000862E+00 /
  data ae11cs(34) /    0.000000000000001223E+00 /
  data ae11cs(35) /    0.000000000000000716E+00 /
  data ae11cs(36) /   -0.000000000000000024E+00 /
  data ae11cs(37) /   -0.000000000000000201E+00 /
  data ae11cs(38) /   -0.000000000000000082E+00 /
  data ae11cs(39) /    0.000000000000000017E+00 /

  data ae12cs( 1) /    0.58241749513472674E+00 /
  data ae12cs( 2) /   -0.15834885090578275E+00 /
  data ae12cs( 3) /   -0.006764275590323141E+00 /
  data ae12cs( 4) /    0.005125843950185725E+00 /
  data ae12cs( 5) /    0.000435232492169391E+00 /
  data ae12cs( 6) /   -0.000143613366305483E+00 /
  data ae12cs( 7) /   -0.000041801320556301E+00 /
  data ae12cs( 8) /   -0.000002713395758640E+00 /
  data ae12cs( 9) /    0.000001151381913647E+00 /
  data ae12cs(10) /    0.000000420650022012E+00 /
  data ae12cs(11) /    0.000000066581901391E+00 /
  data ae12cs(12) /    0.000000000662143777E+00 /
  data ae12cs(13) /   -0.000000002844104870E+00 /
  data ae12cs(14) /   -0.000000000940724197E+00 /
  data ae12cs(15) /   -0.000000000177476602E+00 /
  data ae12cs(16) /   -0.000000000015830222E+00 /
  data ae12cs(17) /    0.000000000002905732E+00 /
  data ae12cs(18) /    0.000000000001769356E+00 /
  data ae12cs(19) /    0.000000000000492735E+00 /
  data ae12cs(20) /    0.000000000000093709E+00 /
  data ae12cs(21) /    0.000000000000010707E+00 /
  data ae12cs(22) /   -0.000000000000000537E+00 /
  data ae12cs(23) /   -0.000000000000000716E+00 /
  data ae12cs(24) /   -0.000000000000000244E+00 /
  data ae12cs(25) /   -0.000000000000000058E+00 /

  data e11cs( 1) /  -16.113461655571494026E+00 /
  data e11cs( 2) /    7.7940727787426802769E+00 /
  data e11cs( 3) /   -1.9554058188631419507E+00 /
  data e11cs( 4) /    0.37337293866277945612E+00 /
  data e11cs( 5) /   -0.05692503191092901938E+00 /
  data e11cs( 6) /    0.00721107776966009185E+00 /
  data e11cs( 7) /   -0.00078104901449841593E+00 /
  data e11cs( 8) /    0.00007388093356262168E+00 /
  data e11cs( 9) /   -0.00000620286187580820E+00 /
  data e11cs(10) /    0.00000046816002303176E+00 /
  data e11cs(11) /   -0.00000003209288853329E+00 /
  data e11cs(12) /    0.00000000201519974874E+00 /
  data e11cs(13) /   -0.00000000011673686816E+00 /
  data e11cs(14) /    0.00000000000627627066E+00 /
  data e11cs(15) /   -0.00000000000031481541E+00 /
  data e11cs(16) /    0.00000000000001479904E+00 /
  data e11cs(17) /   -0.00000000000000065457E+00 /
  data e11cs(18) /    0.00000000000000002733E+00 /
  data e11cs(19) /   -0.00000000000000000108E+00 /

  data e12cs( 1) /   -0.037390214792202795E+00 /
  data e12cs( 2) /    0.042723986062209577E+00 /
  data e12cs( 3) /   -0.1303182079849700544E+00 /
  data e12cs( 4) /    0.01441912402469889073E+00 /
  data e12cs( 5) /   -0.00134617078051068022E+00 /
  data e12cs( 6) /    0.00010731029253063780E+00 /
  data e12cs( 7) /   -0.00000742999951611943E+00 /
  data e12cs( 8) /    0.00000045377325690753E+00 /
  data e12cs( 9) /   -0.00000002476417211390E+00 /
  data e12cs(10) /    0.00000000122076581374E+00 /
  data e12cs(11) /   -0.00000000005485141480E+00 /
  data e12cs(12) /    0.00000000000226362142E+00 /
  data e12cs(13) /   -0.00000000000008635897E+00 /
  data e12cs(14) /    0.00000000000000306291E+00 /
  data e12cs(15) /   -0.00000000000000010148E+00 /
  data e12cs(16) /    0.00000000000000000315E+00 /

  data ae13cs( 1) /   -0.60577324664060346E+00 /
  data ae13cs( 2) /   -0.11253524348366090E+00 /
  data ae13cs( 3) /    0.013432266247902779E+00 /
  data ae13cs( 4) /   -0.001926845187381145E+00 /
  data ae13cs( 5) /    0.000309118337720603E+00 /
  data ae13cs( 6) /   -0.000053564132129618E+00 /
  data ae13cs( 7) /    0.000009827812880247E+00 /
  data ae13cs( 8) /   -0.000001885368984916E+00 /
  data ae13cs( 9) /    0.000000374943193568E+00 /
  data ae13cs(10) /   -0.000000076823455870E+00 /
  data ae13cs(11) /    0.000000016143270567E+00 /
  data ae13cs(12) /   -0.000000003466802211E+00 /
  data ae13cs(13) /    0.000000000758754209E+00 /
  data ae13cs(14) /   -0.000000000168864333E+00 /
  data ae13cs(15) /    0.000000000038145706E+00 /
  data ae13cs(16) /   -0.000000000008733026E+00 /
  data ae13cs(17) /    0.000000000002023672E+00 /
  data ae13cs(18) /   -0.000000000000474132E+00 /
  data ae13cs(19) /    0.000000000000112211E+00 /
  data ae13cs(20) /   -0.000000000000026804E+00 /
  data ae13cs(21) /    0.000000000000006457E+00 /
  data ae13cs(22) /   -0.000000000000001568E+00 /
  data ae13cs(23) /    0.000000000000000383E+00 /
  data ae13cs(24) /   -0.000000000000000094E+00 /
  data ae13cs(25) /    0.000000000000000023E+00 /

  data ae14cs( 1) /   -0.1892918000753017E+00 /
  data ae14cs( 2) /   -0.08648117855259871E+00 /
  data ae14cs( 3) /    0.00722410154374659E+00 /
  data ae14cs( 4) /   -0.00080975594575573E+00 /
  data ae14cs( 5) /    0.00010999134432661E+00 /
  data ae14cs( 6) /   -0.00001717332998937E+00 /
  data ae14cs( 7) /    0.00000298562751447E+00 /
  data ae14cs( 8) /   -0.00000056596491457E+00 /
  data ae14cs( 9) /    0.00000011526808397E+00 /
  data ae14cs(10) /   -0.00000002495030440E+00 /
  data ae14cs(11) /    0.00000000569232420E+00 /
  data ae14cs(12) /   -0.00000000135995766E+00 /
  data ae14cs(13) /    0.00000000033846628E+00 /
  data ae14cs(14) /   -0.00000000008737853E+00 /
  data ae14cs(15) /    0.00000000002331588E+00 /
  data ae14cs(16) /   -0.00000000000641148E+00 /
  data ae14cs(17) /    0.00000000000181224E+00 /
  data ae14cs(18) /   -0.00000000000052538E+00 /
  data ae14cs(19) /    0.00000000000015592E+00 /
  data ae14cs(20) /   -0.00000000000004729E+00 /
  data ae14cs(21) /    0.00000000000001463E+00 /
  data ae14cs(22) /   -0.00000000000000461E+00 /
  data ae14cs(23) /    0.00000000000000148E+00 /
  data ae14cs(24) /   -0.00000000000000048E+00 /
  data ae14cs(25) /    0.00000000000000016E+00 /
  data ae14cs(26) /   -0.00000000000000005E+00 /

  data ntae11 / 0 /
  data ntae12 / 0 /
  data nte11 / 0 /
  data nte12 / 0 /
  data ntae13 / 0 /
  data ntae14 / 0 /
  data xmax / 0.0E+00 /

  if ( ntae11 == 0 ) then
    eta = 0.1E+00 * r4_mach ( 3 )
    ntae11 = r4_inits ( ae11cs, 39, eta )
    ntae12 = r4_inits ( ae12cs, 25, eta )
    nte11 = r4_inits ( e11cs, 19, eta )
    nte12 = r4_inits ( e12cs, 16, eta )
    ntae13 = r4_inits ( ae13cs, 25, eta )
    ntae14 = r4_inits ( ae14cs, 26, eta )
    xmax = - log ( r4_mach ( 1 ) )
    xmax = xmax - log ( xmax )
  end if

  if ( x <= - 10.0E+00 ) then

    r4_e1 = exp ( - x ) / x * ( 1.0E+00 &
      + r4_csevl ( 20.0E+00 / x + 1.0E+00, ae11cs, ntae11 ) )

  else if ( x <= - 4.0E+00 ) then

    r4_e1 = exp ( - x ) / x * ( 1.0E+00 + r4_csevl ( &
      ( 40.0E+00 / x + 7.0E+00 ) / 3.0E+00, ae12cs, ntae12 ) )

  else if ( x <= - 1.0E+00 ) then

    r4_e1 = - log ( abs ( x ) ) + r4_csevl ( &
      ( 2.0E+00 * x + 5.0E+00 ) / 3.0E+00, e11cs, nte11 )

  else if ( x == 0.0E+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_E1 - Fatal error!'
    write ( *, '(a)' ) '  X is zero.'
    stop

  else if ( x <= 1.0E+00 ) then

    r4_e1 = ( - log ( abs ( x ) ) - 0.6875E+00 + x ) &
      + r4_csevl ( x, e12cs, nte12 )

  else if ( x <= 4.0E+00 ) then

    r4_e1 = exp ( - x ) / x * ( 1.0E+00 + r4_csevl ( &
      ( 8.0E+00 / x - 5.0E+00 ) / 3.0E+00, ae13cs, ntae13 ) )

  else if ( x <= xmax ) then

    r4_e1 = exp ( - x ) / x * ( 1.0E+00 + r4_csevl ( &
      8.0E+00 / x - 1.0E+00, ae14cs, ntae14 ) )

  else

    r4_e1 = 0.0E+00

  end if

  return
end
function r4_ei ( x )

!*****************************************************************************80
!
!! R4_EI evaluates the exponential integral Ei for an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_EI, the exponential integral Ei..
!
  implicit none

  real ( kind = 4 ) r4_e1
  real ( kind = 4 ) r4_ei
  real ( kind = 4 ) x

  r4_ei = - r4_e1 ( - x )

  return
end
function r4_erf ( x )

!*****************************************************************************80
!
!! R4_ERF evaluates the error function of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_ERF, the error function of X.
!
  implicit none

  real ( kind = 4 ) erfcs(13)
  integer ( kind = 4 ) nterf
  real ( kind = 4 ) r4_csevl
  real ( kind = 4 ) r4_erf
  real ( kind = 4 ) r4_erfc
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) sqrtpi
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) xbig
  real ( kind = 4 ) y

  save erfcs
  save nterf
  save sqrtpi
  save xbig

  data erfcs( 1) /   -0.049046121234691808E+00 /
  data erfcs( 2) /   -0.14226120510371364E+00 /
  data erfcs( 3) /    0.010035582187599796E+00 /
  data erfcs( 4) /   -0.000576876469976748E+00 /
  data erfcs( 5) /    0.000027419931252196E+00 /
  data erfcs( 6) /   -0.000001104317550734E+00 /
  data erfcs( 7) /    0.000000038488755420E+00 /
  data erfcs( 8) /   -0.000000001180858253E+00 /
  data erfcs( 9) /    0.000000000032334215E+00 /
  data erfcs(10) /   -0.000000000000799101E+00 /
  data erfcs(11) /    0.000000000000017990E+00 /
  data erfcs(12) /   -0.000000000000000371E+00 /
  data erfcs(13) /    0.000000000000000007E+00 /

  data sqrtpi / 1.7724538509055160E+00 /
  data nterf / 0 /
  data xbig / 0.0E+00 /

  if ( nterf == 0 ) then
    nterf = r4_inits ( erfcs, 13, 0.1E+00 * r4_mach ( 3 ) )
    xbig = sqrt ( - log ( sqrtpi * r4_mach ( 3 ) ) )
    sqeps = sqrt ( 2.0E+00 * r4_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    value = 2.0E+00 * x / sqrtpi
  else if ( y <= 1.0E+00 ) then
    value = x * ( 1.0E+00 &
      + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, erfcs, nterf ) )
  else if ( y <= xbig ) then
    value = 1.0E+00 - r4_erfc ( y )
    if ( x < 0.0E+00 ) then
      value = - value
    end if
  else
    value = 1.0E+00
    if ( x < 0.0E+00 ) then
      value = - value
    end if
  end if

  r4_erf = value

  return
end
function r4_erfc ( x )

!*****************************************************************************80
!
!! R4_ERFC evaluates the co-error function of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_ERFC, the co-error function of X.
!
  implicit none

  real ( kind = 4 ) erc2cs(23)
  real ( kind = 4 ) erfccs(24)
  real ( kind = 4 ) erfcs(13)
  real ( kind = 4 ) eta
  integer ( kind = 4 ) nterc2
  integer ( kind = 4 ) nterf
  integer ( kind = 4 ) nterfc
  real ( kind = 4 ) r4_csevl
  real ( kind = 4 ) r4_erfc
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) sqrtpi
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save erfccs
  save erfcs
  save erc2cs
  save nterc2
  save nterf
  save nterfc
  save sqrtpi
  save xmax
  save xsml

  data erfcs( 1) /   -0.049046121234691808E+00 /
  data erfcs( 2) /   -0.14226120510371364E+00 /
  data erfcs( 3) /    0.010035582187599796E+00 /
  data erfcs( 4) /   -0.000576876469976748E+00 /
  data erfcs( 5) /    0.000027419931252196E+00 /
  data erfcs( 6) /   -0.000001104317550734E+00 /
  data erfcs( 7) /    0.000000038488755420E+00 /
  data erfcs( 8) /   -0.000000001180858253E+00 /
  data erfcs( 9) /    0.000000000032334215E+00 /
  data erfcs(10) /   -0.000000000000799101E+00 /
  data erfcs(11) /    0.000000000000017990E+00 /
  data erfcs(12) /   -0.000000000000000371E+00 /
  data erfcs(13) /    0.000000000000000007E+00 /

  data erc2cs( 1) /   -0.069601346602309501E+00 /
  data erc2cs( 2) /   -0.041101339362620893E+00 /
  data erc2cs( 3) /    0.003914495866689626E+00 /
  data erc2cs( 4) /   -0.000490639565054897E+00 /
  data erc2cs( 5) /    0.000071574790013770E+00 /
  data erc2cs( 6) /   -0.000011530716341312E+00 /
  data erc2cs( 7) /    0.000001994670590201E+00 /
  data erc2cs( 8) /   -0.000000364266647159E+00 /
  data erc2cs( 9) /    0.000000069443726100E+00 /
  data erc2cs(10) /   -0.000000013712209021E+00 /
  data erc2cs(11) /    0.000000002788389661E+00 /
  data erc2cs(12) /   -0.000000000581416472E+00 /
  data erc2cs(13) /    0.000000000123892049E+00 /
  data erc2cs(14) /   -0.000000000026906391E+00 /
  data erc2cs(15) /    0.000000000005942614E+00 /
  data erc2cs(16) /   -0.000000000001332386E+00 /
  data erc2cs(17) /    0.000000000000302804E+00 /
  data erc2cs(18) /   -0.000000000000069666E+00 /
  data erc2cs(19) /    0.000000000000016208E+00 /
  data erc2cs(20) /   -0.000000000000003809E+00 /
  data erc2cs(21) /    0.000000000000000904E+00 /
  data erc2cs(22) /   -0.000000000000000216E+00 /
  data erc2cs(23) /    0.000000000000000052E+00 /

  data erfccs( 1) /    0.0715179310202925E+00 /
  data erfccs( 2) /   -0.026532434337606719E+00 /
  data erfccs( 3) /    0.001711153977920853E+00 /
  data erfccs( 4) /   -0.000163751663458512E+00 /
  data erfccs( 5) /    0.000019871293500549E+00 /
  data erfccs( 6) /   -0.000002843712412769E+00 /
  data erfccs( 7) /    0.000000460616130901E+00 /
  data erfccs( 8) /   -0.000000082277530261E+00 /
  data erfccs( 9) /    0.000000015921418724E+00 /
  data erfccs(10) /   -0.000000003295071356E+00 /
  data erfccs(11) /    0.000000000722343973E+00 /
  data erfccs(12) /   -0.000000000166485584E+00 /
  data erfccs(13) /    0.000000000040103931E+00 /
  data erfccs(14) /   -0.000000000010048164E+00 /
  data erfccs(15) /    0.000000000002608272E+00 /
  data erfccs(16) /   -0.000000000000699105E+00 /
  data erfccs(17) /    0.000000000000192946E+00 /
  data erfccs(18) /   -0.000000000000054704E+00 /
  data erfccs(19) /    0.000000000000015901E+00 /
  data erfccs(20) /   -0.000000000000004729E+00 /
  data erfccs(21) /    0.000000000000001432E+00 /
  data erfccs(22) /   -0.000000000000000439E+00 /
  data erfccs(23) /    0.000000000000000138E+00 /
  data erfccs(24) /   -0.000000000000000048E+00 /

  data sqrtpi / 1.7724538509055160E+00 /
  data nterf / 0 /
  data nterfc / 0 /
  data nterc2 / 0 /
  data xsml / 0.0E+00 /
  data xmax / 0.0E+00 /

  if ( nterf == 0 ) then

    eta = 0.1E+00 * r4_mach ( 3 )
    nterf = r4_inits ( erfcs, 13, eta )
    nterfc = r4_inits ( erfccs, 24, eta )
    nterc2 = r4_inits ( erc2cs, 23, eta )

    xsml = - sqrt ( - log ( sqrtpi * r4_mach ( 3 ) ) )
    xmax = sqrt ( - log ( sqrtpi * r4_mach ( 1 ) ) )
    xmax = xmax - 0.5E+00 * log ( xmax ) / xmax - 0.01E+00
    sqeps = sqrt ( 2.0E+00 * r4_mach ( 3 ) )

  end if

  if ( x <= xsml ) then

    r4_erfc = 2.0E+00
    return

  end if

  if ( xmax < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_ERFC - Warning!'
    write ( *, '(a)' ) '  X so big that ERFC underflows.'
    r4_erfc = 0.0E+00
    return
  end if

  y = abs ( x )

  if ( y < sqeps ) then
    r4_erfc = 1.0E+00 - 2.0E+00 * x / sqrtpi
    return
  else if ( y <= 1.0E+00 ) then
    r4_erfc = 1.0E+00 - x * ( 1.0E+00 &
      + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, erfcs, nterf ) )
    return
  end if

  y = y * y

  if ( y <= 4.0E+00 ) then
    r4_erfc = exp ( - y ) / abs ( x ) * ( 0.5E+00 &
      + r4_csevl ( ( 8.0E+00 / y - 5.0E+00 ) / 3.0E+00, erc2cs, nterc2 ) )
  else 
    r4_erfc = exp ( - y ) / abs ( x ) * ( 0.5E+00 &
      + r4_csevl ( 8.0E+00 / y - 1.0E+00, erfccs, nterfc ) )
  end if

  if ( x < 0.0E+00 ) then
    r4_erfc = 2.0E+00 - r4_erfc
  end if

  return
end
function r4_exp ( x )

!*****************************************************************************80
!
!! R4_EXP evaluates the exponential of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_EXP, the exponential of X.
!
  implicit none

  real ( kind = 4 ) aln216
  real ( kind = 4 ) expcs(8)
  real ( kind = 4 ) f
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n16
  integer ( kind = 4 ) ndx
  integer ( kind = 4 ) nterms
  real ( kind = 4 ) r4_csevl
  real ( kind = 4 ) r4_exp
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_pak
  real ( kind = 4 ) twon16(17)
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) xint
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) y

  save aln216
  save expcs
  save nterms
  save twon16
  save xmax
  save xmin

  data expcs( 1) / 0.086656949331498571E+00 /
  data expcs( 2) / 0.000938494869299839E+00 /
  data expcs( 3) / 0.000006776039709981E+00 /
  data expcs( 4) / 0.000000036693120039E+00 /
  data expcs( 5) / 0.000000000158959053E+00 /
  data expcs( 6) / 0.000000000000573859E+00 /
  data expcs( 7) / 0.000000000000001775E+00 /
  data expcs( 8) / 0.000000000000000004E+00 /

  data twon16( 1) / 0.0E+00 /
  data twon16( 2) / 0.44273782427413840E-01 /
  data twon16( 3) / 0.90507732665257659E-01 /
  data twon16( 4) / 0.13878863475669165E+00 /
  data twon16( 5) / 0.18920711500272107E+00 /
  data twon16( 6) / 0.24185781207348405E+00 /
  data twon16( 7) / 0.29683955465100967E+00 /
  data twon16( 8) / 0.35425554693689273E+00 /
  data twon16( 9) / 0.41421356237309505E+00 /
  data twon16(10) / 0.47682614593949931E+00 /
  data twon16(11) / 0.54221082540794082E+00 /
  data twon16(12) / 0.61049033194925431E+00 /
  data twon16(13) / 0.68179283050742909E+00 /
  data twon16(14) / 0.75625216037329948E+00 /
  data twon16(15) / 0.83400808640934246E+00 /
  data twon16(16) / 0.91520656139714729E+00 /
  data twon16(17) / 1.0E+00 /

  data aln216 / 0.083120654223414518E+00 /
  data nterms / 0 /
  data xmax / 0.0E+00 /
  data xmin / 0.0E+00 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( expcs, 8, 0.1E+00 * r4_mach ( 3 ) )
    xmin = log ( r4_mach ( 1 ) ) + 0.01E+00
    xmax = log ( r4_mach ( 2 ) ) - 0.001E+00
  end if

  if ( x < xmin ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_EXP - Warning!'
    write ( *, '(a)' ) '  X so small that exp(X) underflows.'
    value = 0.0E+00

  else if ( x <= xmax ) then

    xint = aint ( x )
    y = x - xint

    y = 23.0E+00 * y + x * aln216
    n = int ( y )
    f = y - real ( n, kind = 4 )
    n = 23.0E+00 * xint + real ( n, kind = 4 )
    n16 = n / 16
    if ( n < 0 ) then
      n16 = n16 - 1
    end if
    ndx = n - 16 * n16 + 1

    value = 1.0E+00 + ( twon16(ndx) &
      + f * ( 1.0E+00 + twon16(ndx) ) &
      * r4_csevl ( f, expcs, nterms ) )

    value = r4_pak ( value, n16 )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_EXP - Fatal error!'
    write ( *, '(a)' ) '  X so large that exp(X) overflows.'
    stop

  end if

  r4_exp = value

  return
end
function r4_exprel ( x )

!*****************************************************************************80
!
!! R4_EXPREL evaluates the exponential relative error term of an R4 argument.
!
!  Discussion:
!
!    The relative error term is ( exp ( x ) - 1 ) / x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_EXPREL, the exponential relative error term
!    at X.
!
  implicit none

  real ( kind = 4 ) absx
  real ( kind = 4 ) alneps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nterms
  real ( kind = 4 ) r4_exprel
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xbnd
  real ( kind = 4 ) xln
  real ( kind = 4 ) xn

  save nterms
  save xbnd

  data nterms / 0 /
  data xbnd / 0.0E+00 /

  if ( nterms == 0 ) then
    alneps = log ( r4_mach ( 3 ) )
    xn = 3.72E+00 - 0.3E+00 * alneps
    xln = log ( ( xn + 1.0E+00 ) / 1.36E+00 )
    nterms = int ( xn - ( xn * xln + alneps ) / ( xln + 1.36E+00 ) + 1.5E+00 )
    xbnd = r4_mach ( 3 )
  end if

  absx = abs ( x )

  if ( absx < xbnd ) then
    r4_exprel = 1.0E+00
  else if ( absx <= 0.5E+00 ) then
    r4_exprel = 0.0E+00
    do i = 1, nterms
      r4_exprel = 1.0E+00 + r4_exprel * x / real ( nterms + 2 - i, kind = 4 )
    end do
  else
    r4_exprel = ( exp ( x ) - 1.0E+00 ) / x
  end if

  return
end
function r4_fac ( n )

!*****************************************************************************80
!
!! R4_FAC evaluates the factorial of an I4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument.
!
!    Output, real ( kind = 4 ) R4_FAC, the factorial of N.
!
  implicit none

  real ( kind = 4 ) facn(26)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nmax
  real ( kind = 4 ) r4_fac
  real ( kind = 4 ) r4_lgmc
  real ( kind = 4 ) sq2pil
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin

  save facn
  save nmax
  save sq2pil

  data facn( 1) / 1.0E+00 /
  data facn( 2) / 1.0E+00 /
  data facn( 3) / 2.0E+00 /
  data facn( 4) / 6.0E+00 /
  data facn( 5) / 24.0E+00 /
  data facn( 6) / 120.0E+00 /
  data facn( 7) / 720.0E+00 /
  data facn( 8) / 5040.0E+00 /
  data facn( 9) / 40320.0E+00 /
  data facn(10) / 362880.0E+00 /
  data facn(11) / 3628800.0E+00 /
  data facn(12) / 39916800.0E+00 /
  data facn(13) / 479001600.0E+00 /
  data facn(14) / 6227020800.0E+00 /
  data facn(15) / 87178291200.0E+00 /
  data facn(16) / 1307674368000.0E+00 /
  data facn(17) / 20922789888000.0E+00 /
  data facn(18) / 355687428096000.0E+00 /
  data facn(19) / 6402373705728000.0E+00 /
  data facn(20) /  0.12164510040883200E+18 /
  data facn(21) /  0.24329020081766400E+19 /
  data facn(22) /  0.51090942171709440E+20 /
  data facn(23) /  0.11240007277776077E+22 /
  data facn(24) /  0.25852016738884977E+23 /
  data facn(25) /  0.62044840173323944E+24 /
  data facn(26) /  0.15511210043330986E+26 /

  data nmax / 0 /
  data sq2pil / 0.91893853320467274E+00 /

  if ( nmax == 0 ) then
    call r4_gaml ( xmin, xmax )
    nmax = int ( xmax - 1.0E+00 )
  end if

  if ( n < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_FAC - Fatal error!'
    write ( *, '(a)' ) '  Input argument is negative.'
    stop

  else if ( n <= 25 ) then

    r4_fac = facn(n+1)
 
  else if ( n <= nmax ) then

    x = real ( n + 1, kind = 4 )
    r4_fac = exp ( ( x - 0.5E+00 ) * log ( x ) - x + sq2pil + r4_lgmc ( x ) )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_FAC - Fatal error!'
    write ( *, '(a)' ) '  Factorial overflows.'
    stop

  end if

  return
end
function r4_gami ( a, x )

!*****************************************************************************80
!
!! R4_GAMI evaluates the incomplete gamma function for an R4 argument.
!
!  Discussion:
!
!    GAMI = Integral ( 0 <= T <= X ) exp ( - t ) * t^( a - 1 )  dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the parameter.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_GAMI, the value of the incomplete 
!    gamma function.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) factor
  real ( kind = 4 ) r4_gami
  real ( kind = 4 ) r4_gamit
  real ( kind = 4 ) r4_lngam
  real ( kind = 4 ) x

  if ( a <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_GAMI - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    stop
  end if

  if ( x < 0.0E+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_GAMI - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop

  else if ( x == 0.0E+00 ) then

    r4_gami = 0.0E+00

  else

    factor = exp ( r4_lngam ( a ) + a * log ( x ) )

    r4_gami = factor * r4_gamit ( a, x )

  end if

  return
end
function r4_gamic ( a, x )

!*****************************************************************************80
!
!! R4_GAMIC evaluates the complementary incomplete gamma function.
!
!  Discussion:
!
!    GAMIC = integral ( x <= t < oo ) exp(-t) * t^(a-1) dt
!
!    GAMIC is evaluated for arbitrary real values of A and non-negative
!    values X (even though GAMIC is defined for X < 0.0), except that
!    for X = 0 and A <= 0.0, GAMIC is undefined.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!    Walter Gautschi,
!    A Computational Procedure for Incomplete Gamma Functions,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 4, December 1979, pages 466-481.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the parameter.
!
!    Input, real ( kind = 4 ) X, the evaluation point.
!
!    Output, real ( kind = 4 ) R4_GAMIC, the value of the incomplete 
!    gamma function.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) aeps
  real ( kind = 4 ) algap1
  real ( kind = 4 ) alneps
  real ( kind = 4 ) alngs
  real ( kind = 4 ) alx
  real ( kind = 4 ) bot
  real ( kind = 4 ) e
  real ( kind = 4 ) eps
  real ( kind = 4 ) fm
  real ( kind = 4 ) gstar
  real ( kind = 4 ) h
  integer ( kind = 4 ) izero
  integer ( kind = 4 ) ma
  real ( kind = 4 ) r4_gamic
  real ( kind = 4 ) r4_gmic
  real ( kind = 4 ) r4_gmit
  real ( kind = 4 ) r4_lgic
  real ( kind = 4 ) r4_lgit
  real ( kind = 4 ) r4_lngam
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_sign
  real ( kind = 4 ) sga
  real ( kind = 4 ) sgng
  real ( kind = 4 ) sgngam
  real ( kind = 4 ) sgngs
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) t
  real ( kind = 4 ) x

  save alneps
  save bot
  save eps

  data alneps / 0.0E+00 /
  data bot / 0.0E+00 /
  data eps / 0.0E+00 /

  if ( eps == 0.0E+00 ) then
    eps = 0.5E+00 * r4_mach ( 3 )
    sqeps = sqrt ( r4_mach ( 4 ) )
    alneps = - log ( r4_mach ( 3 ) )
    bot = log ( r4_mach ( 1 ) )
  end if

  if ( x < 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_GAMIC - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop
  end if

  if ( x == 0.0E+00 ) then

    if ( a <= 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_GAMIC - Fatal error!'
      write ( *, '(a)' ) '  X = 0 and A <= 0.'
      stop
    end if

    r4_gamic = exp ( r4_lngam ( a + 1.0E+00 ) - log ( a ) )

    return
  end if

  alx = log ( x )
  if ( a < 0.0E+00 ) then
    sga = - 1.0E+00
  else
    sga = + 1.0E+00
  end if

  ma = int ( a + 0.5E+00 * sga )
  aeps = a - real ( ma, kind = 4 )

  izero = 0

  if ( x < 1.0E+00 ) then

    if ( a <= 0.5E+00 .and. abs ( aeps ) <= 0.001E+00 ) then

      fm = - real ( ma, kind = 4 )

      if ( fm <= 1.0E+00 ) then
        e = 2.0E+00
      else
        e = 2.0E+00 * ( fm + 2.0E+00 ) / ( fm * fm - 1.0E+00 )
      end if

      e = e - alx * x ** ( - 0.001E+00 )

      if ( e * abs ( aeps ) <= eps ) then
        r4_gamic = r4_gmic ( a, x, alx )
        return
      end if

    end if

    call r4_lgams ( a + 1.0E+00, algap1, sgngam )
    gstar = r4_gmit ( a, x, algap1, sgngam, alx )

    if ( gstar == 0.0E+00 ) then
      izero = 1
    else
      alngs = log ( abs ( gstar ) )
      sgngs = r4_sign ( gstar )
    end if

  else

    if ( a < x ) then
      r4_gamic = exp ( r4_lgic ( a, x, alx ) )
      return
    end if

    sgngam = 1.0E+00
    algap1 = r4_lngam ( a + 1.0E+00 )
    sgngs = 1.0E+00
    alngs = r4_lgit ( a, x, algap1 )

  end if

  h = 1.0E+00

  if ( izero /= 1 ) then

    t = a * alx + alngs

    if ( alneps < t ) then
      sgng = - sgngs * sga * sgngam
      t = t + algap1 - log ( abs ( a ) )
      r4_gamic = sgng * exp ( t )
      return
    end if

    if ( - alneps < t ) then
      h = 1.0E+00 - sgngs * exp ( t )
    end if

  end if

  sgng = r4_sign ( h ) * sga * sgngam
  t = log ( abs ( h ) ) + algap1 - log ( abs ( a ) )
  r4_gamic = sgng * exp ( t )

  return
end
function r4_gamit ( a, x )

!*****************************************************************************80
!
!! R4_GAMIT evaluates Tricomi's incomplete gamma function for an R4 argument.
!
!  Discussion:
!
!      GAMIT = x^(-a) / gamma(a) 
!        * Integral ( 0 <= t <= x ) exp(-t) * t^(a-1) dt
!
!    with analytic continuation for a <= 0.0.  Gamma(x) is the complete
!    gamma function of X.  GAMIT is evaluated for arbitrary real values of
!    A and for non-negative values of X (even though GAMIT is defined for
!    X < 0.0).
!
!    A slight deterioration of 2 or 3 digits accuracy will occur when
!    gamit is very large or very small in absolute value, because log-
!    arithmic variables are used.  Also, if the parameter A is very close
!    to a negative integer (but not a negative integer), there is a loss
!    of accuracy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!    Walter Gautschi,
!    A Computational Procedure for Incomplete Gamma Functions,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 4, December 1979, pages 466-481.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the parameter.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_GAMIT, the function value.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) aeps
  real ( kind = 4 ) ainta
  real ( kind = 4 ) algap1
  real ( kind = 4 ) alneps
  real ( kind = 4 ) alng
  real ( kind = 4 ) alx
  real ( kind = 4 ) bot
  real ( kind = 4 ) h
  real ( kind = 4 ) r4_gamit
  real ( kind = 4 ) r4_gamr
  real ( kind = 4 ) r4_gmit
  real ( kind = 4 ) r4_lgic
  real ( kind = 4 ) r4_lgit
  real ( kind = 4 ) r4_lngam
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sga
  real ( kind = 4 ) sgngam
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) t
  real ( kind = 4 ) x

  save alneps
  save bot

  data alneps / 0.0E+00 /
  data bot / 0.0E+00 /

  if ( alneps == 0.0E+00 ) then
    alneps = - log ( r4_mach ( 3 ) )
    sqeps = sqrt ( r4_mach ( 4 ) )
    bot = log ( r4_mach ( 1 ) )
  end if

  if ( x < 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_GAMIT - Fatal error!'
    write ( *, '(a)' ) '  X is negative.'
    stop
  else if ( x == 0.0E+00 ) then
    alx = 0.0E+00
  else
    alx = log ( x )
  end if

  if ( a < 0.0E+00 ) then
    sga = - 1.0E+00
  else
    sga = + 1.0E+00
  end if

  ainta = aint ( a + 0.5E+00 * sga )
  aeps = a - ainta

  if ( x == 0.0E+00 ) then
    if ( 0.0E+00 < ainta .or. aeps /= 0.0E+00 ) then
      r4_gamit = r4_gamr ( a + 1.0E+00 )
    else
      r4_gamit = 0.0E+00
    end if
    return
  end if

  if ( x <= 1.0E+00 ) then
    if ( - 0.5E+00 <= a .or. aeps /= 0.0E+00 ) then
      call r4_lgams ( a + 1.0E+00, algap1, sgngam )
    end if
    r4_gamit = r4_gmit ( a, x, algap1, sgngam, alx )
    return
  end if

  if ( x <= a ) then
    t = r4_lgit ( a, x, r4_lngam ( a + 1.0E+00 ) )
    r4_gamit = exp ( t )
    return
  end if

  alng = r4_lgic ( a, x, alx )
!
!  Evaluate in terms of alog(r4_gamic(a,x))
!
  h = 1.0E+00

  if ( aeps /= 0.0E+00 .or. 0.0E+00 < ainta ) then

    call r4_lgams ( a + 1.0E+00, algap1, sgngam )
    t = log ( abs ( a ) ) + alng - algap1

    if ( alneps < t ) then
      t = t - a * alx
      r4_gamit = - sga * sgngam * exp ( t )
      return
    end if

    if ( - alneps < t ) then
      h = 1.0E+00 - sga * sgngam * exp ( t )
    end if

  end if

  t = - a * alx + log ( abs ( h ) )

  if ( h < 0.0E+00 ) then
    r4_gamit = - exp ( t )
  else
    r4_gamit = + exp ( t )
  end if

  return
end
subroutine r4_gaml ( xmin, xmax )

!*****************************************************************************80
!
!! R4_GAML evaluates bounds for an R4 argument of the gamma function.
!
!  Discussion:
!
!    This function calculates the minimum and maximum legal bounds 
!    for X in the evaluation of GAMMA ( X ).
!
!    XMIN and XMAX are not the only bounds, but they are the only 
!    non-trivial ones to calculate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Output, real ( kind = 4 ) XMIN, XMAX, the bounds.
!
  implicit none

  real ( kind = 4 ) alnbig
  real ( kind = 4 ) alnsml
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) xln
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xold

  alnsml = log ( r4_mach ( 1 ) )
  xmin = - alnsml

  do i = 1, 10

    xold = xmin
    xln = log ( xmin )
    xmin = xmin - xmin * ( ( xmin + 0.5E+00 ) * xln - xmin &
      - 0.2258E+00 + alnsml ) / ( xmin * xln + 0.5E+00 )

    if ( abs ( xmin - xold ) < 0.005E+00 ) then

      xmin = - xmin + 0.01E+00

      alnbig = log ( r4_mach ( 2 ) )
      xmax = alnbig

      do j = 1, 10

        xold = xmax
        xln = log ( xmax )
        xmax = xmax - xmax * ( ( xmax - 0.5E+00 ) * xln - xmax &
          + 0.9189E+00 - alnbig ) / ( xmax * xln - 0.5E+00 )

        if ( abs ( xmax - xold ) < 0.005E+00 ) then
          xmax = xmax - 0.01E+00
          xmin = max ( xmin, - xmax + 1.0E+00 )
          return
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_GAML - Fatal error!'
      write ( *, '(a)' ) '  Unable to find XMAX.'
      stop

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R4_GAML - Fatal error!'
  write ( *, '(a)' ) '  Unable to find XMIN.'

  stop
end
function r4_gamma ( x )

!*****************************************************************************80
!
!! R4_GAMMA evaluates the gamma function of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_GAMMA, the gamma function of X.
!
  implicit none

  real ( kind = 4 ) dxrel
  real ( kind = 4 ) gcs(23)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ngcs
  real ( kind = 4 ), save :: pi = 3.14159265358979323846E+00
  real ( kind = 4 ) r4_csevl
  real ( kind = 4 ) r4_gamma
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_lgmc
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sinpiy
  real ( kind = 4 ) sq2pil
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y

  save dxrel
  save gcs
  save ngcs
  save sq2pil
  save xmax
  save xmin
  save xsml

  data gcs( 1) / 0.008571195590989331E+00 /
  data gcs( 2) / 0.004415381324841007E+00 /
  data gcs( 3) / 0.05685043681599363E+00 /
  data gcs( 4) /-0.004219835396418561E+00 /
  data gcs( 5) / 0.001326808181212460E+00 /
  data gcs( 6) /-0.0001893024529798880E+00 /
  data gcs( 7) / 0.0000360692532744124E+00 /
  data gcs( 8) /-0.0000060567619044608E+00 /
  data gcs( 9) / 0.0000010558295463022E+00 /
  data gcs(10) /-0.0000001811967365542E+00 /
  data gcs(11) / 0.0000000311772496471E+00 /
  data gcs(12) /-0.0000000053542196390E+00 /
  data gcs(13) / 0.0000000009193275519E+00 /
  data gcs(14) /-0.0000000001577941280E+00 /
  data gcs(15) / 0.0000000000270798062E+00 /
  data gcs(16) /-0.0000000000046468186E+00 /
  data gcs(17) / 0.0000000000007973350E+00 /
  data gcs(18) /-0.0000000000001368078E+00 /
  data gcs(19) / 0.0000000000000234731E+00 /
  data gcs(20) /-0.0000000000000040274E+00 /
  data gcs(21) / 0.0000000000000006910E+00 /
  data gcs(22) /-0.0000000000000001185E+00 /
  data gcs(23) / 0.0000000000000000203E+00 /

  data dxrel / 0.0E+00 /
  data ngcs / 0 /
  data sq2pil / 0.91893853320467274E+00 /
  data xmax / 0.0E+00 /
  data xmin / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( ngcs == 0 ) then
    ngcs = r4_inits ( gcs, 23, 0.1E+00 * r4_mach ( 3 ) )
    call r4_gaml ( xmin, xmax )
    xsml = exp ( max ( log ( r4_mach ( 1 ) ), &
      - log ( r4_mach ( 2 ) ) ) + 0.01E+00 )
    dxrel = sqrt ( r4_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( y <= 10.0E+00 ) then

    n = int ( x )
    if ( x < 0.0E+00 ) then
      n = n - 1
    end if
    y = x - real ( n, kind = 4 )
    n = n - 1
    r4_gamma = 0.9375E+00 + r4_csevl ( 2.0E+00 * y - 1.0E+00, gcs, ngcs )

    if ( n == 0 ) then

      return

    else if ( n < 0 ) then

      n = - n

      if ( x == 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_GAMMA - Fatal error!'
        write ( *, '(a)' ) '  X is 0.'
        stop
      end if

      if ( x < 0.0E+00 .and. &
        x + real ( n - 2, kind = 4 ) == 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_GAMMA - Fatal error!'
        write ( *, '(a)' ) '  X is a negative integer.'
        stop
      end if

      if ( x < - 0.5E+00 .and. &
        abs ( ( x - aint ( x - 0.5E+00 ) ) / x ) < dxrel ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_GAMMA - Warning!'
        write ( *, '(a)' ) '  X too near a negative integer,'
        write ( *, '(a)' ) '  answer is half precision.'
      end if

      if ( y < xsml ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_GAMMA - Fatal error!'
        write ( *, '(a)' ) '  X is so close to zero that Gamma overflows.'
        write ( *, * ) x
        stop
      end if

      do i = 1, n
        r4_gamma = r4_gamma / ( x + real ( i - 1, kind = 4 ) )
      end do

    else if ( n == 0 ) then

    else

      do i = 1, n
        r4_gamma = ( y + real ( i, kind = 4 ) ) * r4_gamma
      end do

    end if

  else

    if ( xmax < x ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_GAMMA - Fatal error!'
      write ( *, '(a)' ) '  X so big that Gamma overflows.'
      stop
    end if
!
!  Underflow.
!
    if ( x < xmin ) then
      r4_gamma = 0.0E+00
      return
    end if

    r4_gamma = exp ( ( y - 0.5E+00 ) * log ( y ) - y + sq2pil &
      + r4_lgmc ( y ) )

    if ( 0.0E+00 < x ) then
      return
    end if

    if ( abs ( ( x - aint ( x - 0.5E+00 ) ) / x ) < dxrel ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_GAMMA - Warning!'
      write ( *, '(a)' ) '  X too near a negative integer,'
      write ( *, '(a)' ) '  answer is half precision.'
    end if

    sinpiy = sin ( pi * y )

    if ( sinpiy == 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_GAMMA - Fatal error!'
      write ( *, '(a)' ) '  X is a negative integer.'
      stop
    end if

    r4_gamma = - pi / ( y * sinpiy * r4_gamma )

  end if

  return
end
function r4_gamr ( x )

!*****************************************************************************80
!
!! R4_GAMR evaluates the reciprocal gamma function of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_GAMR, the value of the reciprocal gamma
!    function at X.
!
  implicit none

  real ( kind = 4 ) alngx
  real ( kind = 4 ) r4_gamma
  real ( kind = 4 ) r4_gamr
  real ( kind = 4 ) sgngx
  real ( kind = 4 ) x

  if ( x <= 0.0E+00 .and. aint ( x ) == x ) then

    r4_gamr = 0.0E+00
    
  else if ( abs ( x ) <= 10.0E+00 ) then

    r4_gamr = 1.0E+00 / r4_gamma ( x )

  else

    call r4_lgams ( x, alngx, sgngx )
    r4_gamr = sgngx * exp ( - alngx )

  end if

  return
end
function r4_gmic ( a, x, alx )

!*****************************************************************************80
!
!! R4_GMIC: complementary incomplete gamma, small X, A near negative integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the parameter.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Input, real ( kind = 4 ) ALX, the logarithm of X.
!
!    Output, real ( kind = 4 ) R4_GMIC, the complementary incomplete 
!    gamma function.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) alng
  real ( kind = 4 ) alx
  real ( kind = 4 ) bot
  logical converged
  real ( kind = 4 ) eps
  real ( kind = 4 ) euler
  real ( kind = 4 ) fk
  real ( kind = 4 ) fkp1
  real ( kind = 4 ) fm
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ma
  integer ( kind = 4 ) mm1
  real ( kind = 4 ) r4_gmic
  real ( kind = 4 ) r4_lngam
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) s
  real ( kind = 4 ) sgng
  real ( kind = 4 ) t
  real ( kind = 4 ) te
  real ( kind = 4 ) x

  save bot
  save eps
  save euler

  data bot / 0.0E+00 /
  data eps / 0.0E+00 /
  data euler / 0.5772156649015329E+00 /

  if ( eps == 0.0E+00 ) then
    eps = 0.5E+00 * r4_mach ( 3 )
    bot = log ( r4_mach ( 1 ) )
  end if

  if ( 0.0E+00 < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_GMIC - Fatal error!'
    write ( *, '(a)' ) '  A must be near a negative integer.'
    stop
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_GMIC - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  ma = int ( a - 0.5E+00 )
  fm = - real ( ma, kind = 4 )
  m = - ma

  te = 1.0E+00
  t = 1.0E+00
  s = t
  converged = .false.
  do k = 1, 200
    fkp1 = real ( k + 1, kind = 4 )
    te = - x * te / ( fm + fkp1 )
    t = te / fkp1
    s = s + t
    if ( abs ( t ) < eps * s ) then
      converged = .true.
      exit
    end if
  end do

  if ( .not. converged ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_GMIC - Fatal error!'
    write ( *, '(a)' ) '  No convergence after 200 iterations.'
    stop
  end if

  r4_gmic = - alx - euler + x * s / ( fm + 1.0E+00 )

  if ( m == 0 ) then
    return
  else if ( m == 1 ) then
    r4_gmic = - r4_gmic - 1.0E+00 + 1.0E+00 / x
    return
  end if

  te = fm
  t = 1.0E+00
  s = t
  mm1 = m - 1
  do k = 1, mm1
    fk = real ( k, kind = 4 )
    te = - x * te / fk
    t = te / ( fm - fk )
    s = s + t
    if ( abs ( t ) < eps * abs ( s ) ) then
      exit
    end if
  end do

  do k = 1, m
    r4_gmic = r4_gmic + 1.0E+00 / real ( k, kind = 4 )
  end do

  if ( mod ( m, 2 ) == 1 ) then
    sgng = - 1.0E+00
  else
    sgng = + 1.0E+00
  end if

  alng = log ( r4_gmic ) - r4_lngam ( fm + 1.0E+00 )

  if ( bot < alng ) then
    r4_gmic = sgng * exp ( alng )
  else
    r4_gmic = 0.0E+00
  end if

  if ( s /= 0.0E+00 ) then
    r4_gmic = r4_gmic &
      + sign ( exp ( - fm * alx + log ( abs ( s ) / fm ) ), s )
  end if

  return
end
function r4_gmit ( a, x, algap1, sgngam, alx )

!*****************************************************************************80
!
!! R4_GMIT: Tricomi's incomplete gamma function for small X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the parameter.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Input, real ( kind = 4 ) ALGAP1, the logarithm of Gamma ( A + 1 ).
!
!    Input, real ( kind = 4 ) SGNGAM, the sign of Gamma ( A + 1 ).
!
!    Input, real ( kind = 4 ) ALX, the logarithm of X.
!
!    Output, real ( kind = 4 ) R4_GMIT, the Tricomi incomplete gamma function.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) ae
  real ( kind = 4 ) aeps
  real ( kind = 4 ) alg2
  real ( kind = 4 ) algap1
  real ( kind = 4 ) algs
  real ( kind = 4 ) alx
  real ( kind = 4 ) bot
  logical converged
  real ( kind = 4 ) eps
  real ( kind = 4 ) fk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ma
  real ( kind = 4 ) r4_gmit
  real ( kind = 4 ) r4_lngam
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_sign
  real ( kind = 4 ) s
  real ( kind = 4 ) sgng2
  real ( kind = 4 ) sgngam
  real ( kind = 4 ) t
  real ( kind = 4 ) te
  real ( kind = 4 ) x

  save bot
  save eps

  data bot / 0.0E+00 /
  data eps / 0.0E+00 /

  if ( eps == 0.0E+00 ) then
    eps = 0.5E+00 * r4_mach ( 3 )
    bot = log ( r4_mach ( 1 ) )
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_GMIT - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  if ( a < 0.0E+00 ) then
    ma = int ( a - 0.5E+00 )
  else
    ma = int ( a + 0.5E+00 )
  end if

  aeps = a - real ( ma, kind = 4 )

  if ( a < - 0.5E+00 ) then
    ae = aeps
  else
    ae = a
  end if

  t = 1.0E+00
  te = ae
  s = t
  converged = .false.
  do k = 1, 200
    fk = real ( k, kind = 4 )
    te = - x * te / fk
    t = te / ( ae + fk )
    s = s + t
    if ( abs ( t ) < eps * abs ( s ) ) then
      converged = .true.
      exit
    end if
  end do

  if ( .not. converged ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_GMIT - Fatal error!'
    write ( *, '(a)' ) '  No convergence in 200 iterations.'
    stop
  end if

  if ( - 0.5E+00 <= a ) then
    algs = - algap1 + log ( s )
    r4_gmit = exp ( algs )
    return
  end if

  algs = - r4_lngam ( 1.0E+00 + aeps ) + log ( s )
  s = 1.0E+00
  m = - ma - 1
  t = 1.0E+00

  do k = 1, m
    t = x * t / ( aeps - real ( m + 1 - k, kind = 4 ) )
    s = s + t
    if ( abs ( t ) < eps * abs ( s ) ) then
      exit
    end if
  end do

  r4_gmit = 0.0E+00
  algs = - real ( ma, kind = 4 ) * log ( x ) + algs

  if ( s == 0.0E+00 .or. aeps == 0.0E+00 ) then
    r4_gmit = exp ( algs )
    return
  end if

  sgng2 = sgngam * r4_sign ( s )
  alg2 = - x - algap1 + log ( abs ( s ) )

  if ( bot < alg2 ) then
    r4_gmit = sgng2 * exp ( alg2 )
  end if

  if ( bot < algs ) then
    r4_gmit = r4_gmit + exp ( algs )
  end if

  return
end
function r4_inits ( os, nos, eta )

!*****************************************************************************80
!
!! R4_INITS initializes a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DOS(NOS), the Chebyshev coefficients.
!
!    Input, integer ( kind = 4 ) NOS, the number of coefficients.
!
!    Input, real ( kind = 8 ) ETA, the desired accuracy.
!
!    Output, integer ( kind = 4 ) R8_INITS, the number of terms of the series 
!    needed to ensure the requested accuracy.
!
  implicit none

  integer ( kind = 4 ) nos

  real ( kind = 4 ) err
  real ( kind = 4 ) eta
  integer ( kind = 4 ) i
  real ( kind = 4 ) os(nos)
  integer ( kind = 4 ) r4_inits

  if ( nos < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_INITS - Fatal error!'
    write ( *, '(a)' ) '  Number of coefficients < 1.'
    stop
  end if

  err = 0.0E+00
  do i = nos, 1, -1
    err = err + abs ( os(i) )
    if ( eta < err ) then
      r4_inits = i
      return
    end if
  end do

  r4_inits = nos
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R4_INITS - Warning!'
  write ( *, '(a)' ) '  ETA may be too small.'

  return
end
function r4_int ( x )

!*****************************************************************************80
!
!! R4_INT returns the integer part of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_INT, the integer part of X.
!
  implicit none

  integer ( kind = 4 ) expo
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_mach
  integer ( kind = 4 ) ibase
  integer ( kind = 4 ) ipart
  integer ( kind = 4 ) npart
  real ( kind = 4 ) r4_int
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) part
  real ( kind = 4 ) scale
  real ( kind = 4 ) x
  real ( kind = 4 ) xbig
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xscl

  save npart
  save scale
  save xbig
  save xmax

  data npart / 0 /
  data scale / 0.0E+00 /
  data xbig / 0.0E+00 /
  data xmax / 0.0E+00 /

  if ( npart == 0 ) then
    ibase = i4_mach ( 10 )
    xmax = 1.0E+00 / r4_mach ( 4 )
    xbig = min ( real ( i4_mach ( 9 ), kind = 4 ), xmax )
    expo = int ( log ( xbig ) / log ( real ( ibase, kind = 4 ) ) - 0.5E+00 )
    scale = ibase ** expo
    npart = int ( log ( xmax ) / log ( scale ) + 1.0E+00 )
  end if

  if ( x < - xmax ) then

    r4_int = x

  else if ( x < - xbig ) then

    xscl = - x

    do i = 1, npart
      xscl = xscl / scale
    end do

    r4_int = 0.0E+00
    do i = 1, npart
      xscl = xscl * scale
      ipart = int ( xscl )
      part = real ( ipart, kind = 4 )
      xscl = xscl - part
      r4_int = r4_int * scale + part
    end do

    r4_int = - r4_int

  else if ( x < + xbig ) then

    r4_int = int ( x )

  else if ( x < + xmax ) then

    xscl = x

    do i = 1, npart
      xscl = xscl / scale
    end do

    r4_int = 0.0E+00
    do i = 1, npart
      xscl = xscl * scale
      ipart = int ( xscl )
      part = real ( ipart, kind = 4 )
      xscl = xscl - part
      r4_int = r4_int * scale + part
    end do

  else

    r4_int = x

  end if

  return
end
subroutine r4_knus ( xnu, x, bknu, bknu1, iswtch )

!*****************************************************************************80
!
!! R4_KNUS computes a sequence of K Bessel functions.
!
!  Discussion:
!
!    This routine computes Bessel functions 
!      exp(x) * k-sub-xnu (x)  
!    and
!      exp(x) * k-sub-xnu+1 (x) 
!    for 0.0 <= xnu < 1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) XNU, the order parameter.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) BKNU, BKNU1, the two K Bessel functions.
!
!    Output, integer ( kind = 4 ) ISWTCH, ?
!
  implicit none

  real ( kind = 4 ) a(15)
  real ( kind = 4 ) a0
  real ( kind = 4 ) aln2
  real ( kind = 4 ) alnbig
  real ( kind = 4 ) alneps
  real ( kind = 4 ) alnsml
  real ( kind = 4 ) alnz
  real ( kind = 4 ) alpha(15)
  real ( kind = 4 ) an
  real ( kind = 4 ) b0
  real ( kind = 4 ) beta(15)
  real ( kind = 4 ) bknu
  real ( kind = 4 ) bknu0
  real ( kind = 4 ) bknu1
  real ( kind = 4 ) bknud
  real ( kind = 4 ) bn
  real ( kind = 4 ) c0
  real ( kind = 4 ) c0kcs(16)
  real ( kind = 4 ) euler
  real ( kind = 4 ) expx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inu
  integer ( kind = 4 ) iswtch
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntc0k
  integer ( kind = 4 ) nterms
  integer ( kind = 4 ) ntznu1
  real ( kind = 4 ) p1
  real ( kind = 4 ) p2
  real ( kind = 4 ) p3
  real ( kind = 4 ) qq
  real ( kind = 4 ) r4_csevl
  real ( kind = 4 ) r4_gamma
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) result
  real ( kind = 4 ) sqpi2
  real ( kind = 4 ) sqrtx
  real ( kind = 4 ) v
  real ( kind = 4 ) vlnz
  real ( kind = 4 ) x
  real ( kind = 4 ) x2n
  real ( kind = 4 ) x2tov
  real ( kind = 4 ) xi
  real ( kind = 4 ) xmu
  real ( kind = 4 ) xnu
  real ( kind = 4 ) xnusml
  real ( kind = 4 ) xsml
  real ( kind = 4 ) z
  real ( kind = 4 ) znu1cs(12)
  real ( kind = 4 ) ztov

  save aln2
  save alnbig
  save alneps
  save alnsml
  save c0kcs
  save euler
  save ntc0k
  save ntznu1
  save sqpi2
  save xnusml
  save xsml
  save znu1cs

  data c0kcs( 1) /    0.060183057242626108E+00 /
  data c0kcs( 2) /   -0.15364871433017286E+00 /
  data c0kcs( 3) /   -0.011751176008210492E+00 /
  data c0kcs( 4) /   -0.000852487888919795E+00 /
  data c0kcs( 5) /   -0.000061329838767496E+00 /
  data c0kcs( 6) /   -0.000004405228124551E+00 /
  data c0kcs( 7) /   -0.000000316312467283E+00 /
  data c0kcs( 8) /   -0.000000022710719382E+00 /
  data c0kcs( 9) /   -0.000000001630564460E+00 /
  data c0kcs(10) /   -0.000000000117069392E+00 /
  data c0kcs(11) /   -0.000000000008405206E+00 /
  data c0kcs(12) /   -0.000000000000603466E+00 /
  data c0kcs(13) /   -0.000000000000043326E+00 /
  data c0kcs(14) /   -0.000000000000003110E+00 /
  data c0kcs(15) /   -0.000000000000000223E+00 /
  data c0kcs(16) /   -0.000000000000000016E+00 /

  data znu1cs( 1) /    0.20330675699419173E+00 /
  data znu1cs( 2) /    0.14007793341321977E+00 /
  data znu1cs( 3) /    0.007916796961001613E+00 /
  data znu1cs( 4) /    0.000339801182532104E+00 /
  data znu1cs( 5) /    0.000011741975688989E+00 /
  data znu1cs( 6) /    0.000000339357570612E+00 /
  data znu1cs( 7) /    0.000000008425941769E+00 /
  data znu1cs( 8) /    0.000000000183336677E+00 /
  data znu1cs( 9) /    0.000000000003549698E+00 /
  data znu1cs(10) /    0.000000000000061903E+00 /
  data znu1cs(11) /    0.000000000000000981E+00 /
  data znu1cs(12) /    0.000000000000000014E+00 /

  data aln2 / 0.69314718055994531E+00 /
  data alnbig / 0.0E+00 /
  data alneps / 0.0E+00 /
  data alnsml / 0.0E+00 /
  data euler / 0.57721566490153286E+00 /
  data ntc0k / 0 /
  data ntznu1 / 0 /
  data sqpi2 / 1.2533141373155003E+00 /
  data xnusml / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( ntc0k == 0 ) then
    ntc0k = r4_inits ( c0kcs, 16, 0.1E+00 * r4_mach ( 3 ) )
    ntznu1 = r4_inits ( znu1cs, 12, 0.1E+00 * r4_mach ( 3 ) )
    xnusml = sqrt ( r4_mach ( 3 ) / 8.0E+00 )
    xsml = 0.1E+00 * r4_mach ( 3 )
    alnsml = log ( r4_mach ( 1 ) )
    alnbig = log ( r4_mach ( 2 ) )
    alneps = log ( 0.1E+00 * r4_mach ( 3 ) )
  end if

  if ( xnu < 0.0E+00 .or. 1.0E+00 <= xnu ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_KNUS - Fatal error!'
    write ( *, '(a)' ) '  XNU < 0 or. 1 <= XNU.'
    stop
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_KNUS - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  iswtch = 0
!
!  X is small.  Compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
!  then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-.5,+.5)
!  then to (0., .5), because k of negative order (-nu) = k of positive
!  order (+nu).
!
  if ( x <= 2.0E+00 ) then

    if ( 0.5E+00 < xnu ) then
      v = 1.0E+00 - xnu
    else
      v = xnu
    end if
!
!  Carefully find (x/2)^xnu and z^xnu where z = x*x/4.
!
    alnz = 2.0E+00 * ( log ( x ) - aln2 )

    if ( x <= xnu ) then

      if ( alnbig < - 0.5E+00 * xnu * alnz - aln2 - log ( xnu ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_KNUS - Fatal error!'
        write ( *, '(a)' ) '  Small X causing overflow.'
        stop
      end if

    end if

    vlnz = v * alnz
    x2tov = exp ( 0.5E+00 * vlnz )

    if ( alnsml < vlnz ) then
      ztov = x2tov * x2tov
    else
      ztov = 0.0E+00
    end if

    a0 = 0.5E+00 * r4_gamma ( 1.0E+00 + v )
    b0 = 0.5E+00 * r4_gamma ( 1.0E+00 - v )
    c0 = - euler

    if ( 0.5E+00 < ztov .and. xnusml < v ) then
      c0 = - 0.75E+00 + r4_csevl ( ( 8.0E+00 * v ) * v - 1.0E+00, c0kcs, ntc0k )
    end if

    if ( ztov <= 0.5E+00 ) then
      alpha(1) = ( a0 - ztov * b0 ) / v
    else
      alpha(1) = c0 - alnz * ( 0.75E+00 + &
        r4_csevl ( vlnz / 0.35E+00 + 1.0E+00, znu1cs, ntznu1 ) ) * b0
    end if

    beta(1) = - 0.5E+00 * ( a0 + ztov * b0 )

    if ( xsml < x ) then
      z = 0.25E+00 * x * x
    else
      z = 0.0E+00
    end if

    nterms = max ( 2, 11 + int ( ( 8.0E+00 * alnz - 25.19E+00 - alneps ) &
      / ( 4.28E+00 - alnz ) ) )

    do i = 2, nterms
      xi = real ( i - 1, kind = 4 )
      a0 = a0 / ( xi * ( xi - v) )
      b0 = b0 / ( xi * ( xi + v) )
      alpha(i) = ( alpha(i-1) + 2.0E+00 * xi * a0 ) / ( xi * ( xi + v ) )
      beta(i) = ( xi - 0.5E+00 * v ) * alpha(i) - ztov * b0
    end do

    bknu = alpha(nterms)
    bknud = beta(nterms)
    do ii = 2, nterms
      i = nterms + 1 - ii
      bknu = alpha(i) + bknu * z
      bknud = beta(i) + bknud * z
    end do

    expx = exp ( x )
    bknu = expx * bknu / x2tov

    if ( alnbig < - 0.5E+00 * ( xnu + 1.0E+00 ) * alnz &
      - 2.0E+00 * aln2 ) then
      iswtch = 1
      return
    end if

    bknud = expx * bknud * 2.0E+00 / ( x2tov * x )

    if ( xnu <= 0.5E+00 ) then
      bknu1 = v * bknu / x - bknud
      return
    end if

    bknu0 = bknu
    bknu = - v * bknu / x - bknud
    bknu1 = 2.0E+00 * xnu * bknu / x + bknu0
!
!  X is large.  Find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke's
!  rational expansion.
!
  else

    sqrtx = sqrt ( x )

    if ( 1.0E+00 / xsml < x ) then
      bknu = sqpi2 / sqrtx
      bknu1 = bknu
      return
    end if

    an = - 1.56E+00 + 4.0E+00 / x
    bn = - 0.29E+00 - 0.22E+00 / x
    nterms = min ( 15, max ( 3, int ( an + bn * alneps ) ) )

    do inu = 1, 2

      if ( inu == 1 ) then
        if ( xnusml < xnu ) then
          xmu = ( 4.0E+00 * xnu ) * xnu
        else
          xmu = 0.0E+00
        end if
      else
        xmu = 4.0E+00 * ( abs ( xnu ) + 1.0E+00 ) ** 2
      end if

      a(1) = 1.0E+00 - xmu
      a(2) = 9.0E+00 - xmu
      a(3) = 25.0E+00 - xmu

      if ( a(2) == 0.0E+00 ) then

        result = sqpi2 * ( 16.0E+00 * x + xmu + 7.0E+00 ) &
          / ( 16.0E+00 * x * sqrtx )

      else

        alpha(1) = 1.0E+00
        alpha(2) = ( 16.0E+00 * x + a(2) ) / a(2)
        alpha(3) = ( ( 768.0E+00 * x + 48.0E+00 * a(3) ) * x &
          + a(2) * a(3) ) / ( a(2) * a(3) )

        beta(1) = 1.0E+00
        beta(2) = ( 16.0E+00 * x + ( xmu + 7.0E+00 ) ) / a(2)
        beta(3) = ( ( 768.0E+00 * x &
          + 48.0E+00 * ( xmu + 23.0E+00 ) ) * x &
          + ( ( xmu + 62.0E+00 ) * xmu + 129.0E+00 ) ) &
          / ( a(2) * a(3) )

        do i = 4, nterms

          n = i - 1
          x2n = real ( 2 * n - 1, kind = 4 )

          a(i) = ( x2n + 2.0E+00 ) ** 2 - xmu
          qq = 16.0E+00 * x2n / a(i)
          p1 = - x2n * ( real ( 12 * n * n - 20 * n, kind = 4 ) &
            - a(1) ) / ( ( x2n - 2.0E+00 ) * a(i) ) - qq * x
          p2 = ( real ( 12 * n * n - 28 * n + 8, kind = 4 ) &
            - a(1) ) / a(i) - qq * x
          p3 = - x2n * a(i-3) / ( ( x2n - 2.0E+00 ) * a(i) )

          alpha(i) = - p1 * alpha(i-1) &
                     - p2 * alpha(i-2) &
                     - p3 * alpha(i-3)

          beta(i) =  - p1 * beta(i-1) &
                     - p2 * beta(i-2) &
                     - p3 * beta(i-3)

        end do

        result = sqpi2 * beta(nterms) / ( sqrtx * alpha(nterms) )

      end if

      if ( inu == 1 ) then
        bknu = result
      else
        bknu1 = result
      end if

    end do

  end if

  return
end
function r4_lbeta ( a, b )

!*****************************************************************************80
!
!! R4_LBETA evaluates the logarithm of the beta function of R4 arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the arguments.
!
!    Output, real ( kind = 4 ) R4_LBETA, the logarithm of the beta function
!    of A and B.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) corr
  real ( kind = 4 ) p
  real ( kind = 4 ) q
  real ( kind = 4 ) r4_gamma
  real ( kind = 4 ) r4_lbeta
  real ( kind = 4 ) r4_lgmc
  real ( kind = 4 ) r4_lngam
  real ( kind = 4 ) r4_lnrel
  real ( kind = 4 ) sq2pil

  save sq2pil

  data sq2pil / 0.91893853320467274E+00 /

  p = min ( a, b )
  q = max ( a, b )

  if ( p <= 0.0E+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LBETA - Fatal error!'
    write ( *, '(a)' ) '  Both arguments must be greater than 0.'
    stop

  else if ( p < 10.0E+00 .and. q <= 10.0E+00 ) then

    r4_lbeta = log ( r4_gamma ( p ) &
      * ( r4_gamma ( q ) / r4_gamma ( p + q ) ) )

  else if ( p < 10.0E+00 ) then

    corr = r4_lgmc ( q ) - r4_lgmc ( p + q )

    r4_lbeta = r4_lngam ( p ) + corr + p - p * log ( p + q ) + &
      ( q - 0.5E+00 ) * r4_lnrel ( - p / ( p + q ) )

  else

    corr = r4_lgmc ( p ) + r4_lgmc ( q ) - r4_lgmc ( p + q )

    r4_lbeta = - 0.5E+00 * log ( q ) + sq2pil + corr &
      + ( p - 0.5E+00 ) * log ( p / ( p + q ) ) &
      + q * r4_lnrel ( - p / ( p + q ) )

  end if

  return
end
subroutine r4_lgams ( x, algam, sgngam )

!*****************************************************************************80
!
!! R4_LGAMS evaluates the log of |gamma(x)| and sign, for an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) ALGAM, the logarithm of the absolute value of
!    gamma ( X ).
!
!    Output, real ( kind = 4 ) SGNGAM, the sign (+1 or -1 ) of gamma ( X ).
!
  implicit none

  real ( kind = 4 ) algam
  integer ( kind = 4 ) k
  real ( kind = 4 ) r4_lngam
  real ( kind = 4 ) sgngam
  real ( kind = 4 ) x

  algam = r4_lngam ( x )
  sgngam = 1.0E+00

  if ( x <= 0.0E+00 ) then

    k = int ( mod ( - aint ( x ), 2.0E+00 ) + 0.1E+00 )

    if ( k == 0 ) then
      sgngam = - 1.0E+00
    end if

  end if

  return
end
function r4_lgic ( a, x, alx )

!*****************************************************************************80
!
!! R4_LGIC: log complementary incomplete gamma function for large X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the parameter.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Input, real ( kind = 4 ) ALX, the logarithm of X.
!
!    Output, real ( kind = 4 ) R4_LGIC, the log complementary incomplete 
!    gamma function.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) alx
  real ( kind = 4 ) eps
  real ( kind = 4 ) fk
  integer ( kind = 4 ) k
  real ( kind = 4 ) p
  real ( kind = 4 ) r
  real ( kind = 4 ) r4_lgic
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) s
  real ( kind = 4 ) t
  real ( kind = 4 ) x
  real ( kind = 4 ) xma
  real ( kind = 4 ) xpa

  save eps

  data eps / 0.0E+00 /

  if ( eps == 0.0E+00 ) then
    eps = 0.5E+00 * r4_mach ( 3 )
  end if

  xpa = x + 1.0E+00 - a
  xma = x - 1.0E+00 - a

  r = 0.0E+00
  p = 1.0E+00
  s = p
  do k = 1, 200
    fk = real ( k, kind = 4 )
    t = fk * ( a - fk ) * ( 1.0E+00 + r )
    r = - t / ( ( xma + 2.0E+00 * fk ) * ( xpa + 2.0E+00 * fk ) + t )
    p = r * p
    s = s + p
    if ( abs ( p ) < eps * s ) then
      r4_lgic = a * alx - x + log ( s / xpa )
      return
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R4_LGIC - Fatal error!'
  write ( *, '(a)' ) '  No convergence in 200 iterations.'

  stop
end
function r4_lgit ( a, x, algap1 )

!*****************************************************************************80
!
!! R4_LGIT evaluates the log of Tricomi's incomplete gamma function.
!
!  Discussion:
!
!    Perron's continued fraction is used for large X and X <= A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the parameter.
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Input, real ( kind = 4 ) ALGAP1, the logarithm of A+1.
!
!    Output, real ( kind = 4 ) R4_LGIT, the log of Tricomi's incomplete
!    gamma function.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) a1x
  real ( kind = 4 ) algap1
  real ( kind = 4 ) ax
  real ( kind = 4 ) eps
  real ( kind = 4 ) fk
  real ( kind = 4 ) hstar
  integer ( kind = 4 ) k
  real ( kind = 4 ) p
  real ( kind = 4 ) r
  real ( kind = 4 ) r4_lgit
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) s
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) t
  real ( kind = 4 ) x

  save eps

  data eps / 0.0E+00 /

  if ( eps == 0.0E+00 ) then
    eps = 0.5E+00 * r4_mach ( 3 )
    sqeps = sqrt ( r4_mach ( 4 ) )
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LGIT - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  if ( a < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LGIT - Fatal error!'
    write ( *, '(a)' ) '  A < X.'
    stop
  end if

  ax = a + x
  a1x = ax + 1.0E+00
  r = 0.0E+00
  p = 1.0E+00
  s = p
  do k = 1, 200
    fk = real ( k, kind = 4 )
    t = ( a + fk ) * x * ( 1.0E+00 + r )
    r = t / ( ( ax + fk ) * ( a1x + fk ) - t )
    p = r * p
    s = s + p
    if ( abs ( p ) < eps * s ) then
      hstar = 1.0E+00 - x * s / a1x
      r4_lgit = - x - algap1 - log ( hstar )
      return
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R4_LGIT - Fatal error!'
  write ( *, '(a)' ) '  No convergence after 200 iterations.'
  stop
end
function r4_lgmc ( x )

!*****************************************************************************80
!
!! R4_LGMC evaluates the log gamma correction factor for an R4 argument.
!
!  Discussion:
!
!    For 10 <= X, compute the log gamma correction factor so that
!
!      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
!                          + ( x - 0.5 ) * log ( x ) - x 
!                          + r4_lgmc ( x )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_LGMC, the correction factor.
!
  implicit none

  real ( kind = 4 ) algmcs(6)
  integer ( kind = 4 ) nalgm
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_lgmc
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xbig
  real ( kind = 4 ) xmax

  save algmcs
  save nalgm
  save xbig
  save xmax

  data algmcs( 1) /    0.166638948045186E+00 /
  data algmcs( 2) /   -0.0000138494817606E+00 /
  data algmcs( 3) /    0.0000000098108256E+00 /
  data algmcs( 4) /   -0.0000000000180912E+00 /
  data algmcs( 5) /    0.0000000000000622E+00 /
  data algmcs( 6) /   -0.0000000000000003E+00 /

  data nalgm / 0 /
  data xbig / 0.0E+00 /
  data xmax / 0.0E+00 /

  if ( nalgm == 0 ) then
    nalgm = r4_inits ( algmcs, 6, r4_mach ( 3 ) )
    xbig = 1.0E+00 / sqrt ( r4_mach ( 3 ) )
    xmax = exp ( min ( log ( r4_mach ( 2 ) / 12.0E+00 ), &
      - log ( 12.0E+00 * r4_mach ( 1 ) ) ) )
  end if

  if ( x < 10.0E+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LGMC - Fatal error!'
    write ( *, '(a)' ) '  X must be at least 10.'
    stop

  else if ( x < xbig ) then

    r4_lgmc = r4_csevl ( 2.0E+00 * ( 10.0E+00 / x ) &
      * ( 10.0E+00 / x ) - 1.0E+00, algmcs, nalgm ) / x

  else if ( x < xmax ) then

    r4_lgmc = 1.0E+00 / ( 12.0E+00 * x )

  else

    r4_lgmc = 0.0E+00

  end if

  return
end
function r4_li ( x )

!*****************************************************************************80
!
!! R4_LI evaluates the logarithmic integral for an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_LI, the logarithmic integral evaluated at X.
!
  implicit none

  real ( kind = 4 ) r4_ei
  real ( kind = 4 ) r4_li
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) x

  if ( sqeps == 0.0E+00 ) then
    sqeps = sqrt ( r4_mach ( 3 ) )
  end if

  if ( x < 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LI - Fatal error!'
    write ( *, '(a)' ) '  Function undefined for X <= 0.'
    stop
  end if

  if ( x == 0.0E+00 ) then
    r4_li = 0.0E+00
    return
  end if

  if ( x == 1.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LI - Fatal error!'
    write ( *, '(a)' ) '  Function undefined for X = 1.'
    stop
  end if

  if ( abs ( 1.0E+00 - x ) < sqeps ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LI - Warning!'
    write ( *, '(a)' ) '  Answer less than half precision.'
    write ( *, '(a)' ) '  X is too close to 1.'
  end if

  r4_li = r4_ei ( log ( x ) )

  return
end
function r4_lngam ( x )

!*****************************************************************************80
!
!! R4_LNGAM evaluates the log of the absolute value of gamma of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_LNGAM, the logarithm of the absolute value of
!    the gamma function of X.
!
  implicit none

  real ( kind = 4 ) dxrel
  real ( kind = 4 ), save :: pi = 3.14159265358979323846E+00
  real ( kind = 4 ) r4_gamma
  real ( kind = 4 ) r4_lgmc
  real ( kind = 4 ) r4_lngam
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) sinpiy
  real ( kind = 4 ) sq2pil
  real ( kind = 4 ) sqpi2l
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) y

  save dxrel
  save sq2pil
  save sqpi2l
  save xmax
  
  data dxrel / 0.0E+00 /
  data sq2pil / 0.91893853320467274E+00 /
  data sqpi2l / 0.22579135264472743E+00 /
  data xmax / 0.0E+00 /

  if ( xmax == 0.0E+00 ) then
    xmax = r4_mach ( 2 ) / log ( r4_mach ( 2 ) )
    dxrel = sqrt ( r4_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( y <= 10.0E+00 ) then
    r4_lngam = log ( abs ( r4_gamma ( x ) ) )
    return
  end if

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LNGAM - Fatal error!'
    write ( *, '(a)' ) '  Result overflows, |X| too big.'
    stop
  end if

  if ( 0.0E+00 < x ) then
    r4_lngam = sq2pil + ( x - 0.5E+00 ) * log ( x ) - x + r4_lgmc ( y )
    return
  end if

  sinpiy = abs ( sin ( pi * y ) )

  if ( sinpiy == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LNGAM - Fatal error!'
    write ( *, '(a)' ) '  X is a negative integer.'
    stop
  end if

  r4_lngam = sqpi2l + ( x - 0.5E+00 ) * log ( y ) - x &
    - log ( sinpiy ) - r4_lgmc ( y )

  if ( abs ( ( x - aint ( x - 0.5E+00 ) ) * r4_lngam / x ) < dxrel ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LNGAM - Warning!'
    write ( *, '(a)' ) '  Result is half precision because'
    write ( *, '(a)' ) '  X is too near a negative integer.'
  end if

  return
end
function r4_lnrel ( x )

!*****************************************************************************80
!
!! R4_LNREL evaluates log ( 1 + X ) for an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_LNREL, the value of LOG ( 1 + X ).
!
  implicit none

  real ( kind = 4 ) alnrcs(23)
  integer ( kind = 4 ) nlnrel
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_lnrel
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) x
  real ( kind = 4 ) xmin

  save alnrcs
  save nlnrel
  save xmin

  data alnrcs( 1) /    1.0378693562743770E+00 /
  data alnrcs( 2) /   -0.13364301504908918E+00 /
  data alnrcs( 3) /    0.019408249135520563E+00 /
  data alnrcs( 4) /   -0.003010755112753577E+00 /
  data alnrcs( 5) /    0.000486946147971548E+00 /
  data alnrcs( 6) /   -0.000081054881893175E+00 /
  data alnrcs( 7) /    0.000013778847799559E+00 /
  data alnrcs( 8) /   -0.000002380221089435E+00 /
  data alnrcs( 9) /    0.000000416404162138E+00 /
  data alnrcs(10) /   -0.000000073595828378E+00 /
  data alnrcs(11) /    0.000000013117611876E+00 /
  data alnrcs(12) /   -0.000000002354670931E+00 /
  data alnrcs(13) /    0.000000000425227732E+00 /
  data alnrcs(14) /   -0.000000000077190894E+00 /
  data alnrcs(15) /    0.000000000014075746E+00 /
  data alnrcs(16) /   -0.000000000002576907E+00 /
  data alnrcs(17) /    0.000000000000473424E+00 /
  data alnrcs(18) /   -0.000000000000087249E+00 /
  data alnrcs(19) /    0.000000000000016124E+00 /
  data alnrcs(20) /   -0.000000000000002987E+00 /
  data alnrcs(21) /    0.000000000000000554E+00 /
  data alnrcs(22) /   -0.000000000000000103E+00 /
  data alnrcs(23) /    0.000000000000000019E+00 /

  data nlnrel / 0 /
  data xmin / 0.0E+00 /

  if ( nlnrel == 0 ) then
    nlnrel = r4_inits ( alnrcs, 23, 0.1E+00 * r4_mach ( 3 ) )
    xmin = - 1.0E+00 + sqrt ( r4_mach ( 4 ) )
  end if

  if ( x <= - 1.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LNREL - Fatal error!'
    write ( *, '(a)' ) '  X <= -1.'
    stop
  else if ( x < xmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LNREL - Warning!'
    write ( *, '(a)' ) '  Result is less than half precision.'
    write ( *, '(a)' ) '  X is too close to - 1.'
  end if

  if ( abs ( x ) <= 0.375E+00 ) then
    r4_lnrel = x * ( 1.0E+00 &
      - x * r4_csevl ( x / 0.375E+00, alnrcs, nlnrel ) )
  else
    r4_lnrel = log ( 1.0E+00 + x )
  end if
  return
end
function r4_log ( x )

!*****************************************************************************80
!
!! R4_LOG evaluates the logarithm of an R4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the evaluation point.
!
!    Output, real ( kind = 4 ) R4_LOG, the logarithm of X.
!
  implicit none

  real ( kind = 4 ) aln2
  real ( kind = 4 ) alncen(5)
  real ( kind = 4 ) alncs(6)
  real ( kind = 4 ) center(4)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nterms
  integer ( kind = 4 ) ntrval
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_log
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) t
  real ( kind = 4 ) t2
  real ( kind = 4 ) x
  real ( kind = 4 ) xn
  real ( kind = 4 ) y

  save aln2
  save alncen
  save alncs
  save center
  save nterms

  data alncs(1) / 1.3347199877973882E+00 /
  data alncs(2) / 0.000693756283284112E+00 /
  data alncs(3) / 0.000000429340390204E+00 /
  data alncs(4) / 0.000000000289338477E+00 /
  data alncs(5) / 0.000000000000205125E+00 /
  data alncs(6) / 0.000000000000000150E+00 /

  data center(1) / 1.0E+00 /
  data center(2) / 1.25E+00 /
  data center(3) / 1.50E+00 /
  data center(4) / 1.75E+00 /

  data alncen(1) / 0.0E+00 /
  data alncen(2) / +0.223143551314209755E+00 /
  data alncen(3) / +0.405465108108164381E+00 /
  data alncen(4) / +0.559615787935422686E+00 /
  data alncen(5) / +0.693147180559945309E+00 /

  data aln2 / 0.068147180559945309E+00 /
  data nterms / 0 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( alncs, 6, 28.9E+00 * r4_mach ( 3 ) )
  end if

  if ( x <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_LOG - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.0'
    stop
  end if

  call r4_upak ( x, y, n )

  xn = real ( n - 1, kind = 4 )
  y = 2.0E+00 * y
  ntrval = int ( 4.0E+00 * y - 2.5E+00 )

  if ( ntrval == 5 ) then
    t = ( ( y - 1.0E+00 ) - 1.0E+00 ) / ( y + 2.0E+00 )
  else if ( ntrval < 5 ) then
    t = ( y - center(ntrval) ) / ( y + center(ntrval) )
  end if

  t2 = t * t

  r4_log = 0.625E+00 * xn + ( aln2 * xn + alncen(ntrval) + 2.0E+00 * t &
    + t * t2 * r4_csevl ( 578.0E+00 * t2 - 1.0, alncs, nterms ) )

  return
end
function r4_log10 ( x )

!*****************************************************************************80
!
!! R4_LOG10 evaluates the logarithm, base 10, of an R4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the evaluation point.
!
!    Output, real ( kind = 4 ) R4_LOG10, the logarithm, base 10, of X.
!
  implicit none

  real ( kind = 4 ) aloge
  real ( kind = 4 ) r4_log10
  real ( kind = 4 ) x

  save aloge

  data aloge / 0.43429448190325182765E+00 /

  r4_log10 = aloge * log ( x )

  return
end
function r4_mach ( i )

!*****************************************************************************80
!
!! R4_MACH returns single precision real machine constants.
!
!  Discussion:
!
!    Assume that single precision real numbers are stored with a mantissa 
!    of T digits in base B, with an exponent whose value must lie 
!    between EMIN and EMAX.  Then for values of I between 1 and 5, 
!    R4_MACH will return the following values:
!
!      R4_MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!      R4_MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
!      R4_MACH(3) = B**(-T), the smallest relative spacing.
!      R4_MACH(4) = B**(1-T), the largest relative spacing.
!      R4_MACH(5) = log10(B)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528,
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 5.
!
!    Output, real ( kind = 4 ) R4_MACH, the value of the chosen parameter.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_mach

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r4_mach = 0.0E+00
    stop
  else if ( i == 1 ) then
    r4_mach = 1.1754944E-38
  else if ( i == 2 ) then
    r4_mach = 3.4028235E+38
  else if ( i == 3 ) then
    r4_mach = 5.9604645E-08
  else if ( i == 4 ) then
    r4_mach = 1.1920929E-07
  else if ( i == 5 ) then
    r4_mach = 0.3010300E+00
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r4_mach = 0.0E+00
    stop
  end if

  return
end
subroutine r4_machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, &
  minexp, maxexp, eps, epsneg, xmin, xmax )

!*****************************************************************************80
!
!! R4_MACHAR determines single precision machine constants.
!
!  Discussion:
!
!    This routine determines the parameters of the floating-point 
!    arithmetic system specified below.  The determination of the first 
!    three uses an extension of an algorithm due to Malcolm, 
!    incorporating some of the improvements suggested by Gentleman and 
!    Marovich.  
!
!    This routine appeared as ACM algorithm 665.
!
!    An earlier version of this program was published in Cody and Waite.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
!    machine parameters,
!    ACM Transactions on Mathematical Software,
!    Volume 14, Number 4, pages 303-311, 1988.
!
!    William Cody, William Waite,
!    Software Manual for the Elementary Functions,
!    Prentice Hall, 1980.
!
!    Morven Gentleman, Scott Marovich,
!    Communications of the ACM,
!    Volume 17, pages 276-277, 1974.
!
!    Michael Malcolm,
!    Communications of the ACM,
!    Volume 15, pages 949-951, 1972.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IBETA, the radix for the floating-point 
!    representation.
!
!    Output, integer ( kind = 4 ) IT, the number of base IBETA digits 
!    in the floating-point significand.
!
!    Output, integer ( kind = 4 ) IRND:
!    0, if floating-point addition chops.
!    1, if floating-point addition rounds, but not in the IEEE style.
!    2, if floating-point addition rounds in the IEEE style.
!    3, if floating-point addition chops, and there is partial underflow.
!    4, if floating-point addition rounds, but not in the IEEE style, and 
!      there is partial underflow.
!    5, if floating-point addition rounds in the IEEE style, and there is 
!      partial underflow.
!
!    Output, integer ( kind = 4 ) NGRD, the number of guard digits for 
!    multiplication with truncating arithmetic.  It is
!    0, if floating-point arithmetic rounds, or if it truncates and only 
!      IT base IBETA digits participate in the post-normalization shift of the
!      floating-point significand in multiplication;
!    1, if floating-point arithmetic truncates and more than IT base IBETA
!      digits participate in the post-normalization shift of the floating-point
!      significand in multiplication.
!
!    Output, integer ( kind = 4 ) MACHEP, the largest negative integer such that
!      1.0 + real ( IBETA ) ** MACHEP /= 1.0, 
!    except that MACHEP is bounded below by - ( IT + 3 ).
!
!    Output, integer ( kind = 4 ) NEGEPS, the largest negative integer such that
!      1.0 - real ( IBETA ) ** NEGEPS /= 1.0, 
!    except that NEGEPS is bounded below by - ( IT + 3 ).
!
!    Output, integer ( kind = 4 ) IEXP, the number of bits (decimal places 
!    if IBETA = 10) reserved for the representation of the exponent (including 
!    the bias or sign) of a floating-point number.
!
!    Output, integer ( kind = 4 ) MINEXP, the largest in magnitude negative 
!    integer such that
!      real ( IBETA ) ** MINEXP 
!    is positive and normalized.
!
!    Output, integer ( kind = 4 ) MAXEXP, the smallest positive power of 
!    BETA that overflows.
!
!    Output, real ( kind = 4 ) EPS, the smallest positive floating-point 
!    number such that  
!      1.0 + EPS /= 1.0. 
!    in particular, if either IBETA = 2  or IRND = 0, 
!      EPS = real ( IBETA ) ** MACHEP.
!    Otherwise,  
!      EPS = ( real ( IBETA ) ** MACHEP ) / 2.
!
!    Output, real ( kind = 4 ) EPSNEG, a small positive floating-point number 
!    such that
!      1.0 - EPSNEG /= 1.0. 
!    In particular, if IBETA = 2 or IRND = 0, 
!      EPSNEG = real ( IBETA ) ** NEGEPS.
!    Otherwise,  
!      EPSNEG = ( real ( IBETA ) ** NEGEPS ) / 2.  
!    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
!    smallest number that can alter 1.0 by subtraction.
!
!    Output, real ( kind = 4 ) XMIN, the smallest non-vanishing normalized 
!    floating-point power of the radix:
!      XMIN = real ( IBETA ) ** MINEXP
!
!    Output, real ( kind = 4 ) XMAX, the largest finite floating-point number.
!    In particular,
!      XMAX = ( 1.0 - EPSNEG ) * real ( IBETA ) ** MAXEXP
!    On some machines, the computed value of XMAX will be only the second, 
!    or perhaps third, largest number, being too small by 1 or 2 units in 
!    the last digit of the significand.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) beta
  real ( kind = 4 ) betah
  real ( kind = 4 ) betain
  real ( kind = 4 ) eps
  real ( kind = 4 ) epsneg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  integer ( kind = 4 ) nxres
  real ( kind = 4 ) one
  real ( kind = 4 ) t
  real ( kind = 4 ) temp
  real ( kind = 4 ) temp1
  real ( kind = 4 ) tempa
  real ( kind = 4 ) two
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) y
  real ( kind = 4 ) z
  real ( kind = 4 ) zero

  one = real ( 1, kind = 4 )
  two = one + one
  zero = one - one
!
!  Determine IBETA and BETA ala Malcolm.
!
  a = one

  do
  
    a = a + a
    temp = a + one
    temp1 = temp - a

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  b = one

  do

    b = b + b
    temp = a + b
    itemp = int ( temp - a )

    if ( itemp /= 0 ) then
      exit
    end if

  end do

  ibeta = itemp
  beta = real ( ibeta, kind = 4 )
!
!  Determine IT and IRND.
!
  it = 0
  b = one

  do

    it = it + 1
    b = b * beta
    temp = b + one
    temp1 = temp - b

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  irnd = 0
  betah = beta / two
  temp = a + betah

  if ( temp - a /= zero ) then
    irnd = 1
  end if

  tempa = a + beta
  temp = tempa + betah

  if ( irnd == 0 .and. temp - tempa /= zero ) then
    irnd = 2
  end if
!
!  Determine NEGEP and EPSNEG.
!
  negep = it + 3
  betain = one / beta
  a = one
  do i = 1, negep
    a = a * betain
  end do

  b = a

  do

    temp = one - a

    if ( temp - one /= zero ) then
      exit
    end if

    a = a * beta
    negep = negep - 1

  end do

  negep = -negep
  epsneg = a

  if ( ibeta /= 2 .and. irnd /= 0 ) then

    a = ( a * ( one + a ) ) / two
    temp = one - a

    if ( temp - one /= zero ) then
      epsneg = a
    end if

  end if
!
!  Determine MACHEP and EPS.
!
  machep = -it - 3
  a = b

  do

    temp = one + a

    if ( temp - one /= zero ) then
      exit
    end if

    a = a * beta
    machep = machep + 1

  end do

  eps = a
  temp = tempa + beta * ( one + eps )

  if ( ibeta /= 2 .and. irnd /= 0 ) then

    a = ( a * ( one + a ) ) / two
    temp = one + a

    if ( temp - one /= zero ) then
      eps = a
    end if

  end if
!
!  Determine NGRD.
!
  ngrd = 0
  temp = one + eps

  if ( irnd == 0 .and. temp * one - one /= zero ) then
    ngrd = 1
  end if
!
!  Determine IEXP, MINEXP and XMIN.
!
!  Loop to determine largest I and K = 2**I such that (1/BETA) ** (2**(I))
!  does not underflow.  Exit from loop is signaled by an underflow.
!
  i = 0
  k = 1
  z = betain
  t = one + eps
  nxres = 0

  do

    y = z
    z = y * y

    a = z * one
    temp = z * t

    if ( a + a == zero .or. y <= abs ( z ) ) then
      exit
    end if

    temp1 = temp * betain

    if ( temp1 * beta == z ) then
      exit
    end if

    i = i + 1
    k = k + k

  end do
!
!  This segment is for nondecimal machines.
!
  if ( ibeta /= 10 ) then

    iexp = i + 1
    mx = k + k
!
!  This segment is for decimal machines only.
!
  else

    iexp = 2
    iz = ibeta

    do

      if ( k < iz ) then
        exit
      end if

      iz = iz * ibeta
      iexp = iexp + 1

    end do

    mx = iz + iz - 1

  end if
!
!  Loop to determine MINEXP, XMIN.
!  Exit from loop is signaled by an underflow.
!
  do

    xmin = y
    y = y * betain

    a = y * one
    temp = y * t

    if ( a + a == zero .or. xmin <= abs ( y ) ) then
      exit
    end if

    k = k + 1
    temp1 = temp * betain

    if ( temp1 * beta == y ) then
      nxres = 3
      xmin = y
      exit
    end if

  end do

  minexp = -k
!
!  Determine MAXEXP and XMAX.
!
  if ( mx <= k + k - 3 .and. ibeta /= 10 ) then
    mx = mx + mx
    iexp = iexp + 1
  end if

  maxexp = mx + minexp
!
!  Adjust IRND to reflect partial underflow.
!
  irnd = irnd + nxres
!
!  Adjust for IEEE-style machines.
!
  if ( irnd == 2 .or. irnd == 5 ) then
    maxexp = maxexp - 2
  end if
!
!  Adjust for non-IEEE machines with partial underflow.
!
  if ( irnd == 3 .or. irnd == 4 ) then
    maxexp = maxexp - it
  end if
!
!  Adjust for machines with implicit leading bit in binary significand, 
!  and machines with radix point at extreme right of significand.
!
  i = maxexp + minexp

  if ( ibeta == 2 .and. i == 0 ) then
    maxexp = maxexp - 1
  end if

  if ( 20 < i ) then
    maxexp = maxexp - 1
  end if

  if ( a /= y ) then
    maxexp = maxexp - 2
  end if

  xmax = one - epsneg

  if ( xmax * one /= xmax ) then
    xmax = one - beta * epsneg
  end if

  xmax = xmax / ( beta * beta * beta * xmin )

  i = maxexp + minexp + 3

  do j = 1, i

    if ( ibeta == 2 ) then
      xmax = xmax + xmax
    else
      xmax = xmax * beta
    end if

  end do

  return
end
function r4_mop ( i )

!*****************************************************************************80
!
!! R4_MOP returns the I-th power of -1 as an R4.
!
!  Discussion:
!
!    An R4 is a real value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, real ( kind = 4 ) R4_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_mop

  if ( mod ( i, 2 ) == 0 ) then
    r4_mop = + 1.0E+00
  else
    r4_mop = - 1.0E+00
  end if

  return
end
function r4_pak ( y, n )

!*****************************************************************************80
!
!! R4_PAK packs a base 2 exponent into an R4.
!
!  Discussion:
!
!    This routine is almost the inverse of R4_UPAK.  It is not exactly 
!    the inverse, because abs ( x ) need not be between 0.5 and 1.0.  
!    If both R4_PAK and 2.0^n were known to be in range, we could compute
!    R4_PAK = x * 2.0^n .
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) Y, the mantissa.
!
!    Input, integer ( kind = 4 ) N, the exponent.
!
!    Output, real ( kind = 4 ) R4_PAK, the packed value.
!
  implicit none

  real ( kind = 4 ) aln210
  real ( kind = 4 ) aln2b
  integer ( kind = 4 ) i4_mach
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) nsum
  integer ( kind = 4 ) ny
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_pak
  real ( kind = 4 ) value
  real ( kind = 4 ) y

  save aln210
  save nmax
  save nmin

  data aln210 / 3.321928094887362E+00 /
  data nmax / 0 /
  data nmin / 0 /

  if ( nmin == 0 ) then
    aln2b = 1.0E+00
    if ( i4_mach ( 10 ) /= 2 ) then
      aln2b = r4_mach ( 5 ) * aln210
    end if
    nmin = aln2b * real ( i4_mach ( 12 ), kind = 4 )
    nmax = aln2b * real ( i4_mach ( 13 ), kind = 4 )
  end if

  call r4_upak ( y, value, ny )

  nsum = n + ny

  if ( nsum < nmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_PAK - Warning!'
    write ( *, '(a)' ) '  Packed number underflows.'
    r4_pak = 0.0E+00
    return
  end if

  if ( nmax < nsum ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_PAK - Fatal error!'
    write ( *, '(a)' ) '  Packed number overflows.'
    stop
  end if

  do while ( nsum < 0 )
    value = 0.5E+00 * value
    nsum = nsum + 1
  end do

  do while ( 0 < nsum )
    value = 2.0E+00 * value
    nsum = nsum - 1
  end do

  r4_pak = value

  return
end
function r4_poch ( a, x )

!*****************************************************************************80
!
!! R4_POCH evaluates Pochhammer's function of R4 arguments.
!
!  Discussion:
!
!    POCH ( A, X ) = Gamma ( A + X ) / Gamma ( A ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, X, the arguments.
!
!    Output, real ( kind = 4 ) R4_POCH, the Pochhammer function of A and X.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) absa
  real ( kind = 4 ) absax
  real ( kind = 4 ) alnga
  real ( kind = 4 ) alngax
  real ( kind = 4 ) ax
  real ( kind = 4 ) b
  real ( kind = 4 ) cospia
  real ( kind = 4 ) cospix
  real ( kind = 4 ) den
  real ( kind = 4 ) eps
  real ( kind = 4 ) err
  real ( kind = 4 ) errpch
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 4 ), save :: pi = 3.14159265358979323846E+00
  real ( kind = 4 ) r4_fac
  real ( kind = 4 ) r4_gamma
  real ( kind = 4 ) r4_gamr
  real ( kind = 4 ) r4_lgmc
  real ( kind = 4 ) r4_lnrel
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_mop
  real ( kind = 4 ) r4_poch
  real ( kind = 4 ) sgnga
  real ( kind = 4 ) sgngax
  real ( kind = 4 ) sinpia
  real ( kind = 4 ) sinpix
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) x

  save eps

  data eps / 0.0E+00 /

  if ( eps == 0.0E+00 ) then
    eps = r4_mach ( 4 )
    sqeps = sqrt ( eps )
  end if

  ax = a + x

  if ( ax <= 0.0E+00 .and. aint ( ax ) == ax ) then

    if ( 0.0E+00 < a .or. aint ( a ) /= a ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_POCH - Fatal error!'
      write ( *, '(a)' ) '  A + X is nonpositive integer,'
      write ( *, '(a)' ) '  but A is not.'
      stop
    end if
!
!  We know here that both A+X and A are non-positive integers.
!
    if ( x == 0.0E+00 ) then
      r4_poch = 1.0E+00
    else if ( - 20.0E+00 <= min ( a + x, a ) ) then
      n = int ( x )
      r4_poch = r4_mop ( n ) * r4_fac ( - int ( a ) ) &
        / r4_fac ( - int ( a ) - n )
    else
      n = int ( x )
      r4_poch = r4_mop ( n ) * exp ( ( a - 0.5E+00 ) &
        * r4_lnrel ( x / ( a - 1.0E+00 ) ) &
        + x * log ( - a + 1.0E+00 - x ) - x &
        + r4_lgmc ( - a + 1.0E+00 ) &
        - r4_lgmc ( - a - x + 1.0E+00 ) )
    end if

    return

  end if
!
!  Here we know A+X is not zero or a negative integer.
!
  if ( a <= 0.0E+00 .and. aint ( a ) == a ) then
    r4_poch = 0.0E+00
    return
  end if

  n = abs ( x )
!
!  X is a small non-positive integer, presummably a common case.
!
  if ( real ( n, kind = 4 ) == x .and. n <= 20 ) then
    r4_poch = 1.0E+00
    do i = 1, n
      r4_poch = r4_poch * ( a + real ( i - 1, kind = 4 ) )
    end do
    return
  end if

  absax = abs ( a + x )
  absa = abs ( a )

  if ( max ( absax, absa ) <= 20.0E+00 ) then
    r4_poch = r4_gamma ( a + x ) * r4_gamr ( a )
    return
  end if

  if ( 0.5E+00 * absa < abs ( x ) ) then
    call r4_lgams ( a + x, alngax, sgngax )
    call r4_lgams ( a, alnga, sgnga )
    r4_poch = sgngax * sgnga * exp ( alngax - alnga )
    return
  end if
!
!  Here abs ( x ) is small and both abs(a+x) and abs(a) are large.  Thus,
!  a+x and a must have the same sign.  For negative a, we use
!  gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) *
!  sin(pi*a)/sin(pi*(a+x))
!
  if ( a < 0.0E+00 ) then
    b = - a - x + 1.0E+00
  else
    b = a
  end if

  r4_poch = exp ( ( b - 0.5E+00 ) * r4_lnrel ( x / b ) &
    + x * log ( b + x ) - x + r4_lgmc ( b + x ) - r4_lgmc ( b ) )

  if ( 0.0E+00 <= a .or. r4_poch == 0.0E+00 ) then
    return
  end if

  cospix = cos ( pi * x )
  sinpix = sin ( pi * x )
  cospia = cos ( pi * a )
  sinpia = sin ( pi * a )

  errpch = abs ( x ) * ( 1.0E+00 + log ( b ) )
  den = cospix + cospia * sinpix / sinpia
  err = ( abs ( x ) * ( abs ( sinpix ) &
    + abs ( cospia * cospix / sinpia ) ) &
    + abs ( a * sinpix ) / sinpia ** 2 ) * pi
  err = errpch + err / abs ( den )

  r4_poch = r4_poch / den

  return
end
function r4_poch1 ( a, x )

!*****************************************************************************80
!
!! R4_POCH1 evaluates a quantity related to Pochhammer's symbol.
!
!  Discussion:
!
!    Evaluate a generalization of Pochhammer's symbol for special
!    situations that require especially accurate values when x is small in
!      poch1(a,x) = (poch(a,x)-1)/x
!                 = (gamma(a+x)/gamma(a) - 1.0)/x .
!    This specification is particularly suited for stably computing
!    expressions such as
!      (gamma(a+x)/gamma(a) - gamma(b+x)/gamma(b))/x
!           = poch1(a,x) - poch1(b,x)
!    Note that poch1(a,0.0) = psi(a)
!
!    When abs ( x ) is so small that substantial cancellation will occur if
!    the straightforward formula is used, we  use an expansion due
!    to fields and discussed by y. l. luke, the special functions and their
!    approximations, vol. 1, academic press, 1969, page 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the parameter.
!
!    Input, real ( kind = 4 ) X, the evaluation point.
!
!    Output, real ( kind = 4 ) R4_POCH1, the value of the function.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) absa
  real ( kind = 4 ) absx
  real ( kind = 4 ), save :: alneps = 0.0E+00
  real ( kind = 4 ) alnvar
  real ( kind = 4 ) b
  real ( kind = 4 ) bern(9)
  real ( kind = 4 ) binv
  real ( kind = 4 ) bp
  real ( kind = 4 ) gbern(10)
  real ( kind = 4 ) gbk
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ndx
  integer ( kind = 4 ) nterms
  real ( kind = 4 ), save :: pi = 3.14159265358979324E+00
  real ( kind = 4 ) poly1
  real ( kind = 4 ) q
  real ( kind = 4 ) r4_cot
  real ( kind = 4 ) r4_exprel
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_poch
  real ( kind = 4 ) r4_poch1
  real ( kind = 4 ) r4_psi
  real ( kind = 4 ) rho
  real ( kind = 4 ) sinpx2
  real ( kind = 4 ) sinpxx
  real ( kind = 4 ), save :: sqtbig = 0.0
  real ( kind = 4 ) term
  real ( kind = 4 ) trig
  real ( kind = 4 ) var
  real ( kind = 4 ) var2
  real ( kind = 4 ) x

  save bern

  data bern( 1) /   0.83333333333333333E-01 /
  data bern( 2) /  -0.13888888888888889E-02 /
  data bern( 3) /   0.33068783068783069E-04 /
  data bern( 4) /  -0.82671957671957672E-06 /
  data bern( 5) /   0.20876756987868099E-07 /
  data bern( 6) /  -0.52841901386874932E-09 /
  data bern( 7) /   0.13382536530684679E-10 /
  data bern( 8) /  -0.33896802963225829E-12 /
  data bern( 9) /   0.85860620562778446E-14 /

  if ( sqtbig == 0.0E+00 ) then
    sqtbig = 1.0E+00 / sqrt ( 24.0E+00 * r4_mach ( 1 ) )
    alneps = log ( r4_mach ( 3 ) )
  end if

  if ( x == 0.0E+00 ) then
    r4_poch1 = r4_psi ( a )
    return
  end if

  absx = abs ( x )
  absa = abs ( a )

  if ( 0.1E+00 * absa < absx .or. &
    0.1E+00 < absx * log ( max ( absa, 2.0E+00 ) ) ) then
    r4_poch1 = r4_poch ( a, x )
    r4_poch1 = ( r4_poch1 - 1.0E+00 ) / x
    return
  end if

  if ( a < - 0.5E+00 ) then
    bp = 1.0E+00 - a - x
  else
    bp = a
  end if

  if ( bp < 10.0E+00 ) then
    incr = 11.0E+00 - bp
  else
    incr = 0
  end if

  b = bp + real ( incr, kind = 4 )

  var = b + 0.5E+00 * ( x - 1.0E+00 )
  alnvar = log ( var )
  q = x * alnvar

  poly1 = 0.0E+00

  if ( var < sqtbig ) then

    var2 = 1.0E+00 / var / var

    rho = 0.5E+00 * ( x + 1.0E+00 )
    gbern(1) = 1.0E+00
    gbern(2) = - rho / 12.0E+00
    term = var2
    poly1 = gbern(2) * term

    nterms = int ( - 0.5E+00 * alneps / alnvar + 1.0E+00 )

    if ( 9 < nterms ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_POCH1 - Fatal error!'
      write ( *, '(a)' ) '  9 < NTERMS.'
      stop
    end if 

    do k = 2, nterms
      gbk = 0.0E+00
      do j = 1, k
        ndx = k - j + 1
        gbk = gbk + bern(ndx) * gbern(j)
      end do
      gbern(k+1) = - rho * gbk / real ( k, kind = 4 )
      term = term * ( real ( 2 * k - 2, kind = 4 ) - x ) &
        * ( real ( 2 * k - 1, kind = 4 ) - x ) * var2
      poly1 = poly1 + gbern(k+1) * term
    end do

  end if

  poly1 = ( x - 1.0E+00 ) * poly1
  r4_poch1 = r4_exprel ( q ) * ( alnvar + q * poly1 ) + poly1
!
!  We have poch1(b,x).  but bp is small, so we use backwards recursion
!  to obtain poch1(bp,x).
!
  do ii = 1, incr
    i = incr - ii
    binv = 1.0E+00 / ( bp + real ( i, kind = 4 ) )
    r4_poch1 = ( r4_poch1 - binv ) / ( 1.0E+00 + x * binv )
  end do

  if ( bp == a ) then
    return
  end if
!
!  We have poch1(bp,x), but a is lt -0.5.  we therefore use a reflection
!  formula to obtain poch1(a,x).
!
  sinpxx = sin ( pi * x ) / x
  sinpx2 = sin ( 0.5E+00 * pi * x )
  trig = sinpxx * r4_cot ( pi * b ) - 2.0E+00 * sinpx2 * ( sinpx2 / x )

  r4_poch1 = trig + ( 1.0E+00 + x * trig ) * r4_poch1

  return
end
function r4_pow ( x, y )

!*****************************************************************************80
!
!! R4_POW computes a power of an R4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    1 September 2011
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the base.
!
!    Input, real ( kind = 4 ) Y, the power.
!
!    Output, real ( kind = 4 ) R4_POW, the value of X^Y.
!
  implicit none

  real ( kind = 4 ) r4_pow
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) y

  value = x ** y

  r4_pow = value

  return
end
function r4_psi ( x )

!*****************************************************************************80
!
!! R4_PSI evaluates the psi function of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_PSI, the psi function of X.
!
  implicit none

  real ( kind = 4 ) apsics(16)
  real ( kind = 4 ) aux
  real ( kind = 4 ) dxrel
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntapsi
  integer ( kind = 4 ) ntpsi
  real ( kind = 4 ), save :: pi = 3.14159265358979323846E+00
  real ( kind = 4 ) psics(23)
  real ( kind = 4 ) r4_cot
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_psi
  real ( kind = 4 ) x
  real ( kind = 4 ) xbig
  real ( kind = 4 ) y

  save apsics
  save dxrel
  save ntapsi
  save ntpsi
  save psics
  save xbig

  data psics( 1) /   -0.038057080835217922E+00 /
  data psics( 2) /    0.49141539302938713E+00 /
  data psics( 3) /   -0.056815747821244730E+00 /
  data psics( 4) /    0.008357821225914313E+00 /
  data psics( 5) /   -0.001333232857994342E+00 /
  data psics( 6) /    0.000220313287069308E+00 /
  data psics( 7) /   -0.000037040238178456E+00 /
  data psics( 8) /    0.000006283793654854E+00 /
  data psics( 9) /   -0.000001071263908506E+00 /
  data psics(10) /    0.000000183128394654E+00 /
  data psics(11) /   -0.000000031353509361E+00 /
  data psics(12) /    0.000000005372808776E+00 /
  data psics(13) /   -0.000000000921168141E+00 /
  data psics(14) /    0.000000000157981265E+00 /
  data psics(15) /   -0.000000000027098646E+00 /
  data psics(16) /    0.000000000004648722E+00 /
  data psics(17) /   -0.000000000000797527E+00 /
  data psics(18) /    0.000000000000136827E+00 /
  data psics(19) /   -0.000000000000023475E+00 /
  data psics(20) /    0.000000000000004027E+00 /
  data psics(21) /   -0.000000000000000691E+00 /
  data psics(22) /    0.000000000000000118E+00 /
  data psics(23) /   -0.000000000000000020E+00 /

  data apsics( 1) /   -0.0204749044678185E+00 /
  data apsics( 2) /   -0.0101801271534859E+00 /
  data apsics( 3) /    0.0000559718725387E+00 /
  data apsics( 4) /   -0.0000012917176570E+00 /
  data apsics( 5) /    0.0000000572858606E+00 /
  data apsics( 6) /   -0.0000000038213539E+00 /
  data apsics( 7) /    0.0000000003397434E+00 /
  data apsics( 8) /   -0.0000000000374838E+00 /
  data apsics( 9) /    0.0000000000048990E+00 /
  data apsics(10) /   -0.0000000000007344E+00 /
  data apsics(11) /    0.0000000000001233E+00 /
  data apsics(12) /   -0.0000000000000228E+00 /
  data apsics(13) /    0.0000000000000045E+00 /
  data apsics(14) /   -0.0000000000000009E+00 /
  data apsics(15) /    0.0000000000000002E+00 /
  data apsics(16) /   -0.0000000000000000E+00 /

  data dxrel / 0.0E+00 /
  data ntpsi / 0 /
  data ntapsi / 0 /
  data xbig / 0.0E+00 /

  if ( ntpsi == 0 ) then
    ntpsi = r4_inits ( psics, 23, 0.1E+00 * r4_mach ( 3 ) )
    ntapsi = r4_inits ( apsics, 16, 0.1E+00 * r4_mach ( 3 ) )
    xbig = 1.0E+00 / sqrt ( r4_mach ( 3 ) )
    dxrel = sqrt ( r4_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( y < 2.0E+00 ) then

    n = int ( x )
    if ( x < 0.0E+00 ) then
      n = n - 1
    end if
    y = x - real ( n, kind = 4 )
    n = n - 1
    r4_psi = r4_csevl ( 2.0E+00 * y - 1.0E+00, psics, ntpsi )

    if ( n == 0 ) then
      return
    end if

    n = - n

    if ( x == 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_PSI - Fatal error!'
      write ( *, '(a)' ) '  X is zero.'
      stop
    end if

    if ( x < 0.0E+00 .and. x + real ( n - 2, kind = 4 ) == 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_PSI - Fatal error!'
      write ( *, '(a)' ) '  X is a negative integer.'
      stop
    end if

    if ( x < - 0.5E+00 .and. &
      abs ( ( x - aint ( x - 0.5E+00 ) ) / x ) < dxrel ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_PSI - Warning!'
      write ( *, '(a)' ) '  Answer is less than half precision'
      write ( *, '(a)' ) '  because X is near a negative integer.'
    end if

    do i = 1, n
      r4_psi = r4_psi - 1.0E+00 / ( x + real ( i - 1, kind = 4 ) )
    end do

  else

    if ( y < xbig ) then
      aux = r4_csevl ( 8.0E+00 / y / y - 1.0E+00, apsics, ntapsi )
    else
      aux = 0.0E+00
    end if

    if ( x < 0.0E+00 ) then

      r4_psi = log ( abs ( x ) ) - 0.5E+00 / x + aux &
        - pi * r4_cot ( pi * x )

    else if ( 0.0E+00 < x ) then

      r4_psi = log ( x ) - 0.5E+00 / x + aux

    end if

  end if

  return
end
function r4_rand ( r )

!*****************************************************************************80
!
!! R4_RAND is a portable pseudorandom number generator.
!
!  Discussion:
!
!    This pseudo-random number generator is portable amoung a wide
!    variety of computers.  It is undoubtedly not as good as many
!    readily available installation dependent versions, and so this
!    routine is not recommended for widespread usage.  Its redeeming
!    feature is that the exact same random numbers (to within final round-
!    off error) can be generated from machine to machine.  Thus, programs
!    that make use of random numbers can be easily transported to and
!    checked in a new environment.
!
!    The random numbers are generated by the linear congruential
!    method described by Knuth in seminumerical methods (p.9),
!    addison-wesley, 1969.  Given the i-th number of a pseudo-random
!    sequence, the i+1 -st number is generated from
!      x(i+1) = (a*x(i) + c) mod m,
!    where here m = 2^22 = 4194304, c = 1731 and several suitable values
!    of the multiplier a are discussed below.  Both the multiplier a and
!    random number x are represented in real ( kind = 8 ) as two 11-bit
!    words.  The constants are chosen so that the period is the maximum
!    possible, 4194304.
!
!    In order that the same numbers be generated from machine to
!    machine, it is necessary that 23-bit integers be reducible modulo
!    2^11 exactly, that 23-bit integers be added exactly, and that 11-bit
!    integers be multiplied exactly.  Furthermore, if the restart option
!    is used (where r is between 0 and 1), then the product r*2^22 =
!    r*4194304 must be correct to the nearest integer.
!
!    The first four random numbers should be 
!
!      0.0004127026,
!      0.6750836372, 
!      0.1614754200, 
!      0.9086198807.
!
!    The tenth random number is 
!
!      0.5527787209.
!
!    The hundredth random number is 
!
!      0.3600893021.  
!
!    The thousandth number should be 
!
!      0.2176990509.
!
!    In order to generate several effectively independent sequences
!    with the same generator, it is necessary to know the random number
!    for several widely spaced calls.  The I-th random number times 2^22,
!    where I=K*P/8 and P is the period of the sequence (P = 2^22), is
!    still of the form L*P/8.  In particular we find the I-th random
!    number multiplied by 2^22 is given by
!      I   =  0  1*p/8  2*p/8  3*p/8  4*p/8  5*p/8  6*p/8  7*p/8  8*p/8
!      RAND=  0  5*p/8  2*p/8  7*p/8  4*p/8  1*p/8  6*p/8  3*p/8  0
!    thus the 4*P/8 = 2097152 random number is 2097152/2^22.
!
!    Several multipliers have been subjected to the spectral test
!    (see Knuth, p. 82).  Four suitable multipliers roughly in order of
!    goodness according to the spectral test are
!      3146757 = 1536*2048 + 1029 = 2^21 + 2^20 + 2^10 + 5
!      2098181 = 1024*2048 + 1029 = 2^21 + 2^10 + 5
!      3146245 = 1536*2048 +  517 = 2^21 + 2^20 + 2^9 + 5
!      2776669 = 1355*2048 + 1629 = 5^9 + 7^7 + 1
!
!    In the table below log10(NU(I)) gives roughly the number of
!    random decimal digits in the random numbers considered I at a time.
!    C is the primary measure of goodness.  In both cases bigger is better.
!
!                     log10 nu(i)              c(i)
!         a       i=2  i=3  i=4  i=5    i=2  i=3  i=4  i=5
!
!      3146757    3.3  2.0  1.6  1.3    3.1  1.3  4.6  2.6
!      2098181    3.3  2.0  1.6  1.2    3.2  1.3  4.6  1.7
!      3146245    3.3  2.2  1.5  1.1    3.2  4.2  1.1  0.4
!      2776669    3.3  2.1  1.6  1.3    2.5  2.0  1.9  2.6
!     best
!      possible   3.3  2.3  1.7  1.4    3.6  5.9  9.7  14.9
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) R, determines the action.
!    * R = 0.0, the next random number of the sequence is generated.
!    * R < 0.0, the last generated number will be returned for
!    possible use in a restart procedure.
!    * R > 0.0, the sequence of random numbers will start with the
!    seed ( R mod 1 ).  This seed is also returned as the value of
!    R4_RAND provided the arithmetic is done exactly.
!
!    Output, real ( kind = 4 ) R4_RAND, a pseudo-random number between 
!    0.0 and 1.0.
!
  implicit none

  integer ( kind = 4 ), save :: ia0 = 1029
  integer ( kind = 4 ), save :: ia1 = 1536
  integer ( kind = 4 ), save :: ia1ma0 = 507
  integer ( kind = 4 ), save :: ic = 1731
  integer ( kind = 4 ), save :: ix0 = 0
  integer ( kind = 4 ), save :: ix1 = 0
  integer ( kind = 4 ) iy0
  integer ( kind = 4 ) iy1
  real ( kind = 4 ) r
  real ( kind = 4 ) r4_rand

  if ( r == 0.0E+00 ) then
    iy0 = ia0 * ix0
    iy1 = ia1 * ix1 + ia1ma0 * ( ix0 - ix1 ) + iy0
    iy0 = iy0 + ic
    ix0 = mod ( iy0, 2048 )
    iy1 = iy1 + ( iy0 - ix0 ) / 2048
    ix1 = mod ( iy1, 2048 )
  end if

  if ( 0.0 < r ) then
    ix1 = int ( mod ( r, 1.0E+00 ) * 4194304.0E+00 + 0.5E+00 )
    ix0 = mod ( ix1, 2048 )
    ix1 = ( ix1 - ix0 ) / 2048
  end if

  r4_rand = real ( ix1 * 2048 + ix0, kind = 4 )
  r4_rand = r4_rand / 4194304.0E+00
 
  return
end
function r4_randgs ( xmean, sd )

!*****************************************************************************80
!
!! R4_RANDGS generates a normally distributed random number.
!
!  Discussion:
!
!    This function generate a normally distributed random number, that is, 
!    it generates random numbers with a Gaussian distribution.  These 
!    random numbers are not exceptionally good, especially in the tails 
!    of the distribution, but this implementation is simple and suitable 
!    for most applications.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Hamming,
!    Numerical Methods for Scientists and Engineers,
!    Dover, 1986,
!    ISBN: 0486652416,
!    LC: QA297.H28.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) XMEAN, the mean of the Gaussian distribution.
!
!    Input, real ( kind = 4 ) SD, the standard deviation of the 
!    Gaussian function.
!
!    Output, real ( kind = 4 ) R4_RANDGS, a normally distributed random number.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_rand
  real ( kind = 4 ) r4_randgs
  real ( kind = 4 ) sd
  real ( kind = 4 ) xmean

  r4_randgs = - 6.0E+00
  do i = 1, 12
    r4_randgs = r4_randgs + r4_rand ( 0.0E+00 )
  end do

  r4_randgs = xmean + sd * r4_randgs

  return
end
function r4_random ( t, n )

!*****************************************************************************80
!
!! R4_RANDOM is a portable pseudorandom number generator.
!
!  Discussion:
!
!    This random number generator is portable amoung a wide variety of
!    computers.  It generates a random number between 0.0 and 1.0 
!    according to the algorithm presented by Bays and Durham.
!
!    The motivation for using this scheme, which resembles the
!    Maclaren-Marsaglia method, is to greatly increase the period of the
!    random sequence.  If the period of the basic generator is P,
!    then the expected mean period of the sequence generated by this
!    generator is given by
!
!      new mean P = sqrt ( pi * factorial ( N ) / ( 8 * P ) ),
!
!    where factorial ( N ) must be much greater than P in this 
!    asymptotic formula.  Generally, N should be 16 to maybe 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Carter Bays, Stephen Durham,
!    Improving a Poor Random Number Generator,
!    ACM Transactions on Mathematical Software,
!    Volume 2, Number 1, March 1976, pages 59-64.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N.  The absolute value of N is the number 
!    of random numbers in an auxiliary table.  Note though that abs(N)+1 is 
!    the number of items in array T.  If N is positive and differs from its 
!    value in the previous invocation, then the table is initialized for 
!    the new value of N.  If N is negative, abs(N) is the number of items 
!    in an auxiliary table, but the tables are now assumed already to
!    be initialized.  This option enables the user to save the table T at 
!    the end of a long computer run and to restart with the same sequence.  
!    Normally, this function would be called at most once with negative N.  
!    Subsequent invocations would have N positive and of the correct magnitude.
!
!    Input/output, real ( kind = 4 ) T(abs(N)+1), an array of random numbers 
!    from a previous invocation of this function.  Whenever N is positive 
!    and differs from the old N, the table is initialized.  The first 
!    abs(N) numbers are the table discussed in the reference, and the 
!    last value is Y.  This array may be saved in order to restart a sequence.
!
!    Output, real ( kind = 4 ) R4_RANDOM, a random number between 0.0 and 1.0.
!
  implicit none

  real ( kind = 4 ) dummy
  real ( kind = 4 ), save :: floatn = - 1.0E+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: nold = - 1
  real ( kind = 4 ) r4_rand
  real ( kind = 4 ) r4_random
  real ( kind = 4 ) t(*)

  if ( n /= nold ) then

    nold = abs ( n )
    floatn = real ( nold, kind = 4 )
 
    if ( n < 0 ) then

      dummy = r4_rand ( t(nold+1) )

    else

      do i = 1, nold
        t(i) = r4_rand ( 0.0E+00 )
      end do
      t(nold+1) = r4_rand ( 0.0E+00 )

    end if

  end if

  j = int ( t(nold+1) * floatn + 1.0E+00 )
  t(nold+1) = t(j)
  r4_random = t(j)
  t(j) = r4_rand ( 0.0E+00 )

  return
end
function r4_ranf ( sw )

!*****************************************************************************80
!
!! R4_RANF is a driver for R4_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Carter Bays, Stephen Durham,
!    Improving a Poor Random Number Generator,
!    ACM Transactions on Mathematical Software,
!    Volume 2, Number 1, March 1976, pages 59-64.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) SW, chooses the action.
!    0.0 <= SW, compute and return the next random number.
!    0.0 > SW, print the internal table, and return the current (old)
!    random number.
!
!    Output, real ( kind = 4 ) R4_RANF, the random value.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_random
  real ( kind = 4 ) r4_ranf
  real ( kind = 4 ) ranold
  real ( kind = 4 ) sw
  real ( kind = 4 ) t(33)

  save ranold
  save t

  data ranold / 0.0E+00 /

  if ( 0.0E+00 <= sw .or. ranold == 0.0E+00 ) then

    r4_ranf = r4_random ( t, 32 )
    ranold = r4_ranf

  end if

  if ( sw < 0.0E+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Current random number table:'
    write ( *, '(a)' ) ' '
    do i = 1, 33
      write ( *, '(2x,i2,2x,f15.10)' ) i, t(i)
    end do

    r4_ranf = ranold

  end if

  return
end
function r4_ren ( )

!*****************************************************************************80
!
!! R4_REN is a simple random number generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Malcolm Pike, David Hill,
!    Algorithm 266:
!    Pseudo-Random Numbers,
!    Communications of the ACM,
!    Volume 8, Number 10, October 1965, page 605.
!
!  Parameters:
!
!    Output, real ( kind = 4 ) R4_REN, the random value.
!
  implicit none

  integer ( kind = 4 ) iy
  real ( kind = 4 ) r4_ren

  save iy

  data iy / 100001 /

  iy = iy * 125
  iy = iy - ( iy / 2796203 ) * 2796203
  r4_ren = real ( iy, kind = 4 ) / 2796203.0E+00

  return
end
function r4_shi ( x )

!*****************************************************************************80
!
!! R4_SHI evaluates the hyperbolic sine integral Shi of an R4 argument.
!
!  Discussion:
!
!    Shi ( x ) = Integral ( 0 <= t <= x ) sinh ( t ) dt / t
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_SHI, the hyperbolic sine integral Shi 
!    evaluated at X.
!
  implicit none

  real ( kind = 4 ) absx
  integer ( kind = 4 ) nshi
  real ( kind = 4 ) r4_csevl
  real ( kind = 4 ) r4_e1
  real ( kind = 4 ) r4_ei
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_shi
  real ( kind = 4 ) shics(7)
  real ( kind = 4 ) x
  real ( kind = 4 ) xsml

  save nshi
  save shics
  save xsml

  data shics(  1) /     0.0078372685688900950695E+00 /
  data shics(  2) /     0.0039227664934234563973E+00 /
  data shics(  3) /     0.0000041346787887617267E+00 /
  data shics(  4) /     0.0000000024707480372883E+00 /
  data shics(  5) /     0.0000000000009379295591E+00 /
  data shics(  6) /     0.0000000000000002451817E+00 /
  data shics(  7) /     0.0000000000000000000467E+00 /

  data nshi / 0 /
  data xsml / 0.0E+00 /

  if ( nshi == 0 ) then
    nshi = r4_inits ( shics, 7, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( r4_mach ( 3 ) )
  end if

  absx = abs ( x )

  if ( absx <= xsml ) then
    r4_shi = x
  else if ( absx <= 0.375E+00 ) then
    r4_shi = x * ( 1.0E+00 &
      + r4_csevl ( 128.0E+00 * x * x / 9.0E+00 - 1.0E+00, shics, nshi ) )
  else
    r4_shi = 0.5E+00 * ( r4_ei ( x ) + r4_e1 ( x ) )
  end if

  return
end
function r4_si ( x )

!*****************************************************************************80
!
!! R4_SI evaluates the sine integral Si of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_SI, the sine integral Si evaluated at X.
!
  implicit none

  real ( kind = 4 ) absx
  real ( kind = 4 ) cosx
  real ( kind = 4 ) eps
  real ( kind = 4 ) f
  real ( kind = 4 ) g
  integer ( kind = 4 ) nsi
  real ( kind = 4 ), parameter :: pi2 = 1.57079632679489661923E+00
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_si
  real ( kind = 4 ) sics(12)
  real ( kind = 4 ) x
  real ( kind = 4 ) xsml

  save nsi
  save sics
  save xsml

  data sics(  1) /    -0.1315646598184841929E+00 /
  data sics(  2) /    -0.2776578526973601892E+00 /
  data sics(  3) /     0.0354414054866659180E+00 /
  data sics(  4) /    -0.0025631631447933978E+00 /
  data sics(  5) /     0.0001162365390497009E+00 /
  data sics(  6) /    -0.0000035904327241606E+00 /
  data sics(  7) /     0.0000000802342123706E+00 /
  data sics(  8) /    -0.0000000013562997693E+00 /
  data sics(  9) /     0.0000000000179440722E+00 /
  data sics( 10) /    -0.0000000000001908387E+00 /
  data sics( 11) /     0.0000000000000016670E+00 /
  data sics( 12) /    -0.0000000000000000122E+00 /

  data nsi / 0 /
  data xsml / 0.0E+00 /

  if ( nsi == 0 ) then
    nsi = r4_inits ( sics, 12, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( eps )
  end if

  absx = abs ( x )

  if ( absx < xsml ) then

    r4_si = x

  else if ( absx <= 4.0E+00 ) then

    r4_si = x * ( 0.75E+00 &
      + r4_csevl ( ( x * x - 8.0E+00 ) * 0.125E+00, sics, nsi ) )

  else

    call r4_sifg ( absx, f, g )
    cosx = cos ( absx )

    if ( x < 0.0E+00 ) then
      r4_si = - pi2 + f * cosx + g * sin ( x )
    else
      r4_si = pi2 - f * cosx - g * sin ( x )
    end if

  end if

  return
end
subroutine r4_sifg ( x, f, g )

!*****************************************************************************80
!
!! R4_SIFG is a utility routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) F, G.
!
  implicit none

  real ( kind = 4 ) f
  real ( kind = 4 ) f1cs(20)
  real ( kind = 4 ) f2cs(29)
  real ( kind = 4 ) g
  real ( kind = 4 ) g1cs(21)
  real ( kind = 4 ) g2cs(34)
  integer ( kind = 4 ) nf1
  integer ( kind = 4 ) nf2
  integer ( kind = 4 ) ng1
  integer ( kind = 4 ) ng2
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) tol
  real ( kind = 4 ) x
  real ( kind = 4 ) xbig
  real ( kind = 4 ) xbnd
  real ( kind = 4 ) xmaxf
  real ( kind = 4 ) xmaxg

  save f1cs
  save f2cs
  save g1cs
  save g2cs
  save nf1
  save nf2
  save ng1
  save ng2
  save xbig
  save xbnd
  save xmaxf
  save xmaxg

  data f1cs(  1) /    -0.1191081969051363610E+00 /
  data f1cs(  2) /    -0.0247823144996236248E+00 /
  data f1cs(  3) /     0.0011910281453357821E+00 /
  data f1cs(  4) /    -0.0000927027714388562E+00 /
  data f1cs(  5) /     0.0000093373141568271E+00 /
  data f1cs(  6) /    -0.0000011058287820557E+00 /
  data f1cs(  7) /     0.0000001464772071460E+00 /
  data f1cs(  8) /    -0.0000000210694496288E+00 /
  data f1cs(  9) /     0.0000000032293492367E+00 /
  data f1cs( 10) /    -0.0000000005206529618E+00 /
  data f1cs( 11) /     0.0000000000874878885E+00 /
  data f1cs( 12) /    -0.0000000000152176187E+00 /
  data f1cs( 13) /     0.0000000000027257192E+00 /
  data f1cs( 14) /    -0.0000000000005007053E+00 /
  data f1cs( 15) /     0.0000000000000940241E+00 /
  data f1cs( 16) /    -0.0000000000000180014E+00 /
  data f1cs( 17) /     0.0000000000000035063E+00 /
  data f1cs( 18) /    -0.0000000000000006935E+00 /
  data f1cs( 19) /     0.0000000000000001391E+00 /
  data f1cs( 20) /    -0.0000000000000000282E+00 /

  data f2cs(  1) /    -0.0348409253897013234E+00 /
  data f2cs(  2) /    -0.0166842205677959686E+00 /
  data f2cs(  3) /     0.0006752901241237738E+00 /
  data f2cs(  4) /    -0.0000535066622544701E+00 /
  data f2cs(  5) /     0.0000062693421779007E+00 /
  data f2cs(  6) /    -0.0000009526638801991E+00 /
  data f2cs(  7) /     0.0000001745629224251E+00 /
  data f2cs(  8) /    -0.0000000368795403065E+00 /
  data f2cs(  9) /     0.0000000087202677705E+00 /
  data f2cs( 10) /    -0.0000000022601970392E+00 /
  data f2cs( 11) /     0.0000000006324624977E+00 /
  data f2cs( 12) /    -0.0000000001888911889E+00 /
  data f2cs( 13) /     0.0000000000596774674E+00 /
  data f2cs( 14) /    -0.0000000000198044313E+00 /
  data f2cs( 15) /     0.0000000000068641396E+00 /
  data f2cs( 16) /    -0.0000000000024731020E+00 /
  data f2cs( 17) /     0.0000000000009226360E+00 /
  data f2cs( 18) /    -0.0000000000003552364E+00 /
  data f2cs( 19) /     0.0000000000001407606E+00 /
  data f2cs( 20) /    -0.0000000000000572623E+00 /
  data f2cs( 21) /     0.0000000000000238654E+00 /
  data f2cs( 22) /    -0.0000000000000101714E+00 /
  data f2cs( 23) /     0.0000000000000044259E+00 /
  data f2cs( 24) /    -0.0000000000000019634E+00 /
  data f2cs( 25) /     0.0000000000000008868E+00 /
  data f2cs( 26) /    -0.0000000000000004074E+00 /
  data f2cs( 27) /     0.0000000000000001901E+00 /
  data f2cs( 28) /    -0.0000000000000000900E+00 /
  data f2cs( 29) /     0.0000000000000000432E+00 /

  data g1cs(  1) /    -0.3040578798253495954E+00 /
  data g1cs(  2) /    -0.0566890984597120588E+00 /
  data g1cs(  3) /     0.0039046158173275644E+00 /
  data g1cs(  4) /    -0.0003746075959202261E+00 /
  data g1cs(  5) /     0.0000435431556559844E+00 /
  data g1cs(  6) /    -0.0000057417294453025E+00 /
  data g1cs(  7) /     0.0000008282552104503E+00 /
  data g1cs(  8) /    -0.0000001278245892595E+00 /
  data g1cs(  9) /     0.0000000207978352949E+00 /
  data g1cs( 10) /    -0.0000000035313205922E+00 /
  data g1cs( 11) /     0.0000000006210824236E+00 /
  data g1cs( 12) /    -0.0000000001125215474E+00 /
  data g1cs( 13) /     0.0000000000209088918E+00 /
  data g1cs( 14) /    -0.0000000000039715832E+00 /
  data g1cs( 15) /     0.0000000000007690431E+00 /
  data g1cs( 16) /    -0.0000000000001514697E+00 /
  data g1cs( 17) /     0.0000000000000302892E+00 /
  data g1cs( 18) /    -0.0000000000000061400E+00 /
  data g1cs( 19) /     0.0000000000000012601E+00 /
  data g1cs( 20) /    -0.0000000000000002615E+00 /
  data g1cs( 21) /     0.0000000000000000548E+00 /

  data g2cs(  1) /    -0.0967329367532432218E+00 /
  data g2cs(  2) /    -0.0452077907957459871E+00 /
  data g2cs(  3) /     0.0028190005352706523E+00 /
  data g2cs(  4) /    -0.0002899167740759160E+00 /
  data g2cs(  5) /     0.0000407444664601121E+00 /
  data g2cs(  6) /    -0.0000071056382192354E+00 /
  data g2cs(  7) /     0.0000014534723163019E+00 /
  data g2cs(  8) /    -0.0000003364116512503E+00 /
  data g2cs(  9) /     0.0000000859774367886E+00 /
  data g2cs( 10) /    -0.0000000238437656302E+00 /
  data g2cs( 11) /     0.0000000070831906340E+00 /
  data g2cs( 12) /    -0.0000000022318068154E+00 /
  data g2cs( 13) /     0.0000000007401087359E+00 /
  data g2cs( 14) /    -0.0000000002567171162E+00 /
  data g2cs( 15) /     0.0000000000926707021E+00 /
  data g2cs( 16) /    -0.0000000000346693311E+00 /
  data g2cs( 17) /     0.0000000000133950573E+00 /
  data g2cs( 18) /    -0.0000000000053290754E+00 /
  data g2cs( 19) /     0.0000000000021775312E+00 /
  data g2cs( 20) /    -0.0000000000009118621E+00 /
  data g2cs( 21) /     0.0000000000003905864E+00 /
  data g2cs( 22) /    -0.0000000000001708459E+00 /
  data g2cs( 23) /     0.0000000000000762015E+00 /
  data g2cs( 24) /    -0.0000000000000346151E+00 /
  data g2cs( 25) /     0.0000000000000159996E+00 /
  data g2cs( 26) /    -0.0000000000000075213E+00 /
  data g2cs( 27) /     0.0000000000000035970E+00 /
  data g2cs( 28) /    -0.0000000000000017530E+00 /
  data g2cs( 29) /     0.0000000000000008738E+00 /
  data g2cs( 30) /    -0.0000000000000004487E+00 /
  data g2cs( 31) /     0.0000000000000002397E+00 /
  data g2cs( 32) /    -0.0000000000000001347E+00 /
  data g2cs( 33) /     0.0000000000000000801E+00 /
  data g2cs( 34) /    -0.0000000000000000501E+00 /

  data nf1 / 0 /
  data nf2 / 0 /
  data ng1 / 0 /
  data ng2 / 0 /
  data xbnd / 0.0E+00 /
  data xbig / 0.0E+00 /
  data xmaxf / 0.0E+00 /
  data xmaxg / 0.0E+00 /

  if ( nf1 == 0 ) then
    tol = 0.1E+00 * r4_mach ( 3 )
    nf1 = r4_inits ( f1cs, 20, tol )
    nf2 = r4_inits ( f2cs, 29, tol )
    ng1 = r4_inits ( g1cs, 21, tol )
    ng2 = r4_inits ( g2cs, 34, tol )
    xbig = sqrt ( 1.0E+00 / r4_mach ( 3 ) )
    xmaxf = exp ( min ( - log ( r4_mach ( 1 ) ), &
      log ( r4_mach ( 2 ) ) ) - 0.01E+00 )
    xmaxg = 1.0E+00 / sqrt ( r4_mach ( 1 ) )
    xbnd = sqrt ( 50.0E+00 )
  end if

  if ( x < 4.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_SIFG - Fatal error!'
    write ( *, '(a)' ) '  Approximation invalid for X < 4.'
    stop
  end if

  if ( x <= xbnd ) then

    f = ( 1.0E+00 + r4_csevl ( ( 1.0E+00 / x / x - 0.04125E+00 ) &
      / 0.02125E+00, f1cs, nf1 ) ) / x
    g = ( 1.0E+00 + r4_csevl ( ( 1.0E+00 / x / x - 0.04125E+00 ) &
      / 0.02125E+00, g1cs, ng1 ) ) / x / x

  else if ( x <= xbig ) then

    f = ( 1.0E+00 + r4_csevl ( 100.0E+00 / x / x - 1.0E+00, &
      f2cs, nf2 ) ) / x
    g = ( 1.0E+00 + r4_csevl ( 100.0E+00 / x / x - 1.0E+00, &
      g2cs, ng2) ) / x / x

  else

    if ( x < xmaxf ) then
      f = 1.0E+00 / x
    else
      f = 0.0E+00
    end if

    if ( x < xmaxg ) then
      g = 1.0E+00 / x / x
    else
      g = 0.0E+00
    end if

  end if

  return
end
function r4_sign ( x )

!*****************************************************************************80
!
!! R4_SIGN returns the sign of an R4.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value =  0 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 4 ) R4_SIGN, the sign of X:
!
  implicit none

  real ( kind = 4 ) r4_sign
  real ( kind = 4 ) x

  if ( x < 0.0E+00 ) then
    r4_sign = -1.0E+00
  else
    r4_sign = +1.0E+00
  end if

  return
end
function r4_sin ( x )

!*****************************************************************************80
!
!! R4_SIN evaluates the sine of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_SIN, the sine of X.
!
  implicit none

  real ( kind = 4 ) f
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) ntsn
  real ( kind = 4 ) pi2rec
  real ( kind = 4 ) pihi
  real ( kind = 4 ) pilo
  real ( kind = 4 ) pirec
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_sin
  real ( kind = 4 ) sincs(10)
  real ( kind = 4 ) sgn
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xn
  real ( kind = 4 ) xsml
  real ( kind = 4 ) xwarn
  real ( kind = 4 ) y

  save ntsn
  save pi2rec
  save pihi
  save pilo
  save pirec
  save sincs
  save xmax
  save xsml
  save xwarn

  data sincs(  1) / -0.374991154955873175840E+00 /
  data sincs(  2) / -0.181603155237250201864E+00 /
  data sincs(  3) / +0.005804709274598633559E+00 /
  data sincs(  4) / -0.000086954311779340757E+00 /
  data sincs(  5) / +0.000000754370148088851E+00 /
  data sincs(  6) / -0.000000004267129665056E+00 /
  data sincs(  7) / +0.000000000016980422945E+00 /
  data sincs(  8) / -0.000000000000050120579E+00 /
  data sincs(  9) / +0.000000000000000114101E+00 /
  data sincs( 10) / -0.000000000000000000206E+00 /
!
!  pihi + pilo = pi.  pihi is exactly representable on all machines
!  with at least 8 bits of precision.  whether it is exactly
!  represented depends on the compiler.  this routine is more
!  accurate if it is exactly represented.
!
  data ntsn / 0 /
  data pihi / 3.140625E+00 /
  data pilo / 9.6765358979323846E-04 /
  data pirec / 0.31830988618379067E+00 /
  data pi2rec / 0.636619772367581343E+00 /
  data xmax / 0.0E+00 /
  data xsml / 0.0E+00 /
  data xwarn / 0.0E+00 /

  if ( ntsn == 0 ) then
    ntsn = r4_inits ( sincs, 10, 0.1E+00 * r4_mach ( 3 ) )
    xsml = sqrt ( 6.0E+00 * r4_mach ( 3 ) )
    xmax = 1.0E+00 / r4_mach ( 4 )
    xwarn = sqrt ( xmax )
  end if

  y = abs ( x )

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_SIN - Warning!'
    write ( *, '(a)' ) '  No precision because |X| is big.'
    r4_sin = 0.0E+00
    return
  end if

  if ( xwarn < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_SIN - Warning!'
    write ( *, '(a)' ) '  Answer < half precision because |X| is big.'
  end if

  r4_sin = x

  if ( y < xsml ) then
    return
  end if

  xn = aint ( y * pirec + 0.5E+00 )
  n2 = int ( mod ( xn, 2.0E+00 ) + 0.5E+00 )
  sgn = x
  if ( n2 /= 0 ) then
    sgn = - sgn
  end if

  f = ( y - xn * pihi ) - xn * pilo
  xn = 2.0E+00 * ( f * pi2rec ) * ( f * pi2rec ) - 1.0E+00

  r4_sin = f + f * r4_csevl ( xn, sincs, ntsn )

  if ( sgn < 0.0E+00 ) then
    r4_sin = - r4_sin
  end if

  if ( r4_sin < - 1.0E+00 ) then
    r4_sin = - 1.0E+00
  else if ( 1.0E+00 < r4_sin ) then
    r4_sin = 1.0E+00
  end if

  return
end
function r4_sin_deg ( x )

!*****************************************************************************80
!
!! R4_SIN_DEG evaluates the sine of an R4 argument in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument in degrees.
!
!    Output, real ( kind = 4 ) R4_SIN_DEG, the sine of X.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 4 ) r4_sin_deg
  real ( kind = 4 ), parameter :: raddeg = 0.017453292519943296E+00
  real ( kind = 4 ) value
  real ( kind = 4 ) x

  value = sin ( raddeg * x )

  if ( mod ( x, 90.0E+00 ) == 0.0E+00 ) then

    n = int ( abs ( x ) / 90.0E+00 + 0.5E+00 )
    n = mod ( n, 2 )

    if ( n == 0 ) then
      value = 0.0E+00
    else if ( value < 0.0E+00 ) then
      value = - 1.0E+00
    else
      value = + 1.0E+00
    end if

  end if

  r4_sin_deg = value

  return
end
function r4_sinh ( x )

!*****************************************************************************80
!
!! R4_SINH evaluates the hyperbolic sine of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_SINH, the hyperbolic sine of X.
!
  implicit none

  integer ( kind = 4 ) nterms
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_sinh
  real ( kind = 4 ) sinhcs(8)
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) ymax

  save nterms
  save sinhcs
  save ymax

  data sinhcs( 1) /    0.1730421940471796E+00 /
  data sinhcs( 2) /    0.08759422192276048E+00 /
  data sinhcs( 3) /    0.00107947777456713E+00 /
  data sinhcs( 4) /    0.00000637484926075E+00 /
  data sinhcs( 5) /    0.00000002202366404E+00 /
  data sinhcs( 6) /    0.00000000004987940E+00 /
  data sinhcs( 7) /    0.00000000000007973E+00 /
  data sinhcs( 8) /    0.00000000000000009E+00 /

  data nterms / 0 /
  data ymax / 0.0E+00 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( sinhcs, 8, 0.1E+00 * r4_mach ( 3 ) )
    sqeps = sqrt ( 6.0E+00 * r4_mach ( 3 ) )
    ymax = 1.0E+00 / sqrt ( r4_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then

    value = x

  else if ( y <= 1.0E+00 ) then

    value = x * ( 1.0E+00 &
      + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, sinhcs, nterms ) )

  else

    y = exp ( y )

    if ( ymax <= y ) then 
      value = 0.5E+00 * y
    else
      value = 0.5E+00 * ( y - 1.0E+00 / y )
    end if

    if ( x < 0.0E+00 ) then
      value = - value
    end if

  end if

  r4_sinh = value

  return
end
function r4_spence ( x )

!*****************************************************************************80
!
!! R4_SPENCE evaluates a form of Spence's function for an R4 argument.
!
!  Discussion:
!
!    This function evaluates a form of Spence's function defined by
!
!      f(x) = Integral ( 0 <= y <= x ) - log ( abs ( 1 - y ) ) / y dy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions, page 1004,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!    K Mitchell,
!    Tables of the function Integral ( 0 < y < x ) - log | 1 - y | dy / y
!    with an account of some properties of this and related functions,
!    Philosophical Magazine,
!    Volume 40, pages 351-368, 1949.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_SPENCE, Spence's function evaluated at X.
!
  implicit none

  real ( kind = 4 ) aln
  integer ( kind = 4 ) nspenc
  real ( kind = 4 ) pi26
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_spence
  real ( kind = 4 ) spencs(19)
  real ( kind = 4 ) x
  real ( kind = 4 ) xbig

  save nspenc
  save pi26
  save spencs
  save xbig

  data spencs( 1) /    0.1527365598892406E+00 /
  data spencs( 2) /    0.08169658058051014E+00 /
  data spencs( 3) /    0.00581415714077873E+00 /
  data spencs( 4) /    0.00053716198145415E+00 /
  data spencs( 5) /    0.00005724704675185E+00 /
  data spencs( 6) /    0.00000667454612164E+00 /
  data spencs( 7) /    0.00000082764673397E+00 /
  data spencs( 8) /    0.00000010733156730E+00 /
  data spencs( 9) /    0.00000001440077294E+00 /
  data spencs(10) /    0.00000000198444202E+00 /
  data spencs(11) /    0.00000000027940058E+00 /
  data spencs(12) /    0.00000000004003991E+00 /
  data spencs(13) /    0.00000000000582346E+00 /
  data spencs(14) /    0.00000000000085767E+00 /
  data spencs(15) /    0.00000000000012768E+00 /
  data spencs(16) /    0.00000000000001918E+00 /
  data spencs(17) /    0.00000000000000290E+00 /
  data spencs(18) /    0.00000000000000044E+00 /
  data spencs(19) /    0.00000000000000006E+00 /

  data pi26 / 1.644934066848226E+00 /
  data nspenc / 0 /
  data xbig / 0.0E+00 /

  if ( nspenc == 0 ) then
    nspenc = r4_inits ( spencs, 19, 0.1E+00 * r4_mach ( 3 ) )
    xbig = 1.0E+00 / r4_mach ( 3 )
  end if

  if ( x <= - xbig ) then

    aln = log ( 1.0E+00 - x )
    r4_spence = - pi26 &
      - 0.5E+00 * aln * ( 2.0E+00 * log ( - x ) - aln )

  else if ( x <= - 1.0E+00 ) then

    aln = log ( 1.0E+00 - x )

    r4_spence = - pi26 - 0.5E+00 * aln * ( 2.0E+00 &
      * log ( - x ) - aln ) + ( 1.0E+00 + r4_csevl ( &
      4.0E+00 / ( 1.0E+00 - x ) - 1.0E+00, spencs, nspenc ) ) &
   / ( 1.0E+00 - x )

  else if ( x <= 0.0E+00 ) then

    r4_spence = - 0.5E+00 * log ( 1.0E+00 - x ) &
      * log ( 1.0E+00 - x ) - x * ( 1.0E+00 + r4_csevl ( &
      4.0E+00 * x / ( x - 1.0E+00 ) - 1.0E+00, spencs, nspenc ) ) &
   / ( x - 1.0E+00 )

  else if ( x <= 0.5E+00 ) then

    r4_spence = x * ( 1.0E+00 + r4_csevl ( 4.0E+00 * x - 1.0E+00, &
      spencs, nspenc ) )

  else if ( x < 1.0E+00 ) then

    r4_spence = pi26 - log ( x ) * log ( 1.0E+00 - x ) &
      - ( 1.0E+00 - x ) * ( 1.0E+00 + r4_csevl ( 4.0E+00 &
      * ( 1.0E+00 - x ) - 1.0E+00, spencs, nspenc ) )

  else if ( x == 1.0E+00 ) then

    r4_spence = pi26

  else if ( x <= 2.0E+00 ) then

    r4_spence = pi26 - 0.5E+00 * log ( x ) &
      * log ( ( x - 1.0E+00 ) * ( x - 1.0E+00 ) / x ) &
      + ( x - 1.0E+00 ) * ( 1.0E+00 + r4_csevl ( 4.0E+00 &
      * ( x - 1.0E+00 ) / x - 1.0E+00, spencs, nspenc ) ) / x

  else if ( x < xbig ) then

    r4_spence = 2.0E+00 * pi26 - 0.5E+00 * log ( x ) * log ( x ) &
      - ( 1.0E+00 + r4_csevl ( 4.0E+00 / x - 1.0E+00, spencs, &
      nspenc ) ) / x

  else

    r4_spence = 2.0E+00 * pi26 - 0.5E+00 * log ( x ) * log ( x )

  end if

  return
end
function r4_sqrt ( x )

!*****************************************************************************80
!
!! R4_SQRT computes the square root of an R4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the number whose square root is desired.
!
!    Output, real ( kind = 4 ) R4_SQRT, the square root of X.
!
  implicit none

  integer ( kind = 4 ) irem
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) ixpnt
  integer ( kind = 4 ) n
  integer ( kind = 4 ) niter
  real ( kind = 4 ) r4_log
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_pak
  real ( kind = 4 ) r4_sqrt
  real ( kind = 4 ) sqrt2(3)
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) y

  save niter
  save sqrt2

  data sqrt2(1) / 0.70710678118654752E+00 /
  data sqrt2(2) / 1.0E+00 /
  data sqrt2(3) / 1.41421356237309505E+00 /

  data niter / 0 /

  if ( niter == 0 ) then
    niter = 1.443E+00 * r4_log ( - 0.104E+00 &
      * r4_log ( 0.1E+00 * r4_mach ( 3 ) ) ) + 1.0E+00
  end if

  if ( x < 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_SQRT - Fatal error!'
    write ( *, '(a)' ) '  X < 0.0'
    stop
  else if ( x == 0.0E+00 ) then
    r4_sqrt = 0.0E+00
  else

    call r4_upak ( x, y, n )
    ixpnt = n / 2
    irem = n - 2 * ixpnt + 2
    value = 0.261599E+00 + y * ( 1.114292E+00 &
      + y * ( -0.516888E+00 + y * 0.141067E+00 ) )

    do iter = 1, niter
      value = value + 0.5E+00 * ( y - value * value ) / value
    end do

    value = r4_pak ( sqrt2(irem) * value, ixpnt )

    r4_sqrt = value

  end if

  return
end
function r4_tan ( x )

!*****************************************************************************80
!
!! R4_TAN evaluates the tangent of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_TAN, the tangent of X.
!
  implicit none

  real ( kind = 4 ) ainty
  real ( kind = 4 ) ainty2
  integer ( kind = 4 ) ifn
  integer ( kind = 4 ) nterms
  real ( kind = 4 ) pi2rec
  real ( kind = 4 ) prodbg
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_tan
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) tancs(11)
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xsml
  real ( kind = 4 ) y
  real ( kind = 4 ) yrem

  save nterms
  save pi2rec
  save tancs
  save xmax
  save xsml

  data tancs( 1) / 0.226279327631293578E+00 /
  data tancs( 2) / 0.0430179131465489618E+00 /
  data tancs( 3) / 0.0006854461068256508E+00 /
  data tancs( 4) / 0.0000110453269475970E+00 /
  data tancs( 5) / 0.0000001781747790392E+00 /
  data tancs( 6) / 0.0000000028744968582E+00 /
  data tancs( 7) / 0.0000000000463748541E+00 /
  data tancs( 8) / 0.0000000000007481760E+00 /
  data tancs( 9) / 0.0000000000000120704E+00 /
  data tancs(10) / 0.0000000000000001947E+00 /
  data tancs(11) / 0.0000000000000000031E+00 /

  data pi2rec / 0.0116197723675813430E+00 /
  data nterms / 0 /
  data xmax / 0.0E+00 /
  data xsml / 0.0E+00 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( tancs, 11, 0.1E+00 * r4_mach ( 3 ) )
    xmax = 1.0E+00 / r4_mach ( 4 )
    xsml = sqrt ( 3.0E+00 * r4_mach ( 3 ) )
    sqeps = sqrt ( r4_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_TAN - Warning'
    write ( *, '(a)' ) '  No precision because |X| is big.'
    r4_tan = 0.0E+00
    return
  end if
!
!  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
!  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
!  = aint(.625*y) + aint(z) + rem(z)
!
  ainty = aint ( y )
  yrem = y - ainty
  prodbg = 0.625E+00 * ainty
  ainty = aint ( prodbg )
  y = ( prodbg - ainty ) + 0.625E+00 * yrem + y * pi2rec
  ainty2 = aint ( y )
  ainty = ainty + ainty2
  y = y - ainty2

  ifn = int ( mod ( ainty, 2.0E+00 ) )

  if ( ifn == 1 ) then
    y = 1.0E+00 - y
  end if

  if ( 1.0E+00 - y < abs ( x ) * sqeps ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_TAN - Warning!'
    write ( *, '(a)' ) '  Answer < half precision.'
    write ( *, '(a)' ) '  |X| big or X near pi/2 or 3*pi/2.'
  end if

  if ( y == 1.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_TAN - Fatal error!'
    write ( *, '(a)' ) '  X is pi/2 or 3*pi/2.'
    value = 0.0E+00
    stop
  end if

  if ( y <= 0.25E+00 ) then

    value = y
    if ( xsml < y ) then
      value = y * ( 1.5E+00 &
        + r4_csevl ( 32.0E+00 * y * y - 1.0E+00, tancs, nterms ) )
    end if

  else if ( y <= 0.5E+00 ) then

    value = 0.5E+00 * y * ( 1.5E+00 &
      + r4_csevl ( 8.0E+00 * y * y - 1.0E+00, tancs, nterms ) )

    value = 2.0E+00 * value / ( 1.0E+00 - value * value )

  else

    value = 0.25E+00 * y * ( 1.5E+00 &
      + r4_csevl ( 2.0E+00 * y * y - 1.0E+00, tancs, nterms ) )
    value = 2.0E+00 * value / ( 1.0E+00 - value * value )
    value = 2.0E+00 * value / ( 1.0E+00 - value * value )

  end if

  if ( x < 0.0E+00 ) then
    value = - abs ( value )
  else if ( 0.0E+00 < x ) then
    value = + abs ( value )
  end if

  if ( ifn == 1 ) then
    value = - value
  end if

  r4_tan = value

  return
end
function r4_tanh ( x )

!*****************************************************************************80
!
!! R4_TANH evaluates the hyperbolic tangent of an R4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) R4_TANH, the hyperbolic tangent of X.
!
  implicit none

  integer ( kind = 4 ) nterms
  real ( kind = 4 ) r4_csevl
  integer ( kind = 4 ) r4_inits
  real ( kind = 4 ) r4_mach
  real ( kind = 4 ) r4_tanh
  real ( kind = 4 ), save :: sqeps = 0.0E+00
  real ( kind = 4 ) tanhcs(17)
  real ( kind = 4 ) value
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) y

  save nterms
  save tanhcs
  save xmax

  data tanhcs( 1) /   -0.25828756643634710E+00 /
  data tanhcs( 2) /   -0.11836106330053497E+00 /
  data tanhcs( 3) /   +0.009869442648006398E+00 /
  data tanhcs( 4) /   -0.000835798662344582E+00 /
  data tanhcs( 5) /   +0.000070904321198943E+00 /
  data tanhcs( 6) /   -0.000006016424318120E+00 /
  data tanhcs( 7) /   +0.000000510524190800E+00 /
  data tanhcs( 8) /   -0.000000043320729077E+00 /
  data tanhcs( 9) /   +0.000000003675999055E+00 /
  data tanhcs(10) /   -0.000000000311928496E+00 /
  data tanhcs(11) /   +0.000000000026468828E+00 /
  data tanhcs(12) /   -0.000000000002246023E+00 /
  data tanhcs(13) /   +0.000000000000190587E+00 /
  data tanhcs(14) /   -0.000000000000016172E+00 /
  data tanhcs(15) /   +0.000000000000001372E+00 /
  data tanhcs(16) /   -0.000000000000000116E+00 /
  data tanhcs(17) /   +0.000000000000000009E+00 /

  data nterms / 0 /
  data xmax / 0.0E+00 /

  if ( nterms == 0 ) then
    nterms = r4_inits ( tanhcs, 17, 0.1E+00 * r4_mach ( 3 ) )
    sqeps = sqrt ( 3.0E+00 * r4_mach ( 3 ) )
    xmax = - 0.5E+00 * log ( r4_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then

    value = x

  else if ( y <= 1.0E+00 ) then

    value = x * ( 1.0E+00 &
      + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, tanhcs, nterms ) )

  else if ( y <= xmax ) then

    y = exp ( y )
    value = ( y - 1.0E+00 / y ) / ( y + 1.0E+00 / y )
    if ( x < 0.0E+00 ) then
      value = - value
    end if

  else

    if ( x < 0.0E+00 ) then
      value = - 1.0E+00
    else
      value = + 1.0E+00
    end if

  end if

  r4_tanh = value

  return
end
subroutine r4_upak ( x, y, n )

!*****************************************************************************80
!
!! R4_UPAK unpacks an R4 into a mantissa and exponent.
!
!  Discussion:
!
!    This function unpacks a floating point number x so that
!
!      x = y * 2.0^n
!
!    where
!
!      0.5 <= abs ( y ) < 1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the number to be unpacked.
!
!    Output, real ( kind = 4 ) Y, the mantissa.
!
!    Output, integer ( kind = 4 ) N, the exponent.
!
  implicit none

  real ( kind = 4 ) absx
  integer ( kind = 4 ) n
  real ( kind = 4 ) x
  real ( kind = 4 ) y

  absx = abs ( x )
  n = 0
  y = 0.0E+00

  if ( x == 0.0E+00 ) then
    return
  end if

  do while ( absx < 0.5E+00 )
    n = n - 1
    absx = absx * 2.0E+00
  end do

  do while ( 1.0E+00 <= absx )
    n = n + 1
    absx = absx * 0.5E+00
  end do

  if ( x < 0.0E+00 ) then
    y = - absx
  else
    y = + absx
  end if

  return
end
function r8_acos ( x )

!*****************************************************************************80
!
!! R8_ACOS evaluates the arc-cosine of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_ACOS, the arc-cosine of X.
!
  implicit none

  real ( kind = 8 ) pi2
  real ( kind = 8 ) r8_acos
  real ( kind = 8 ) r8_asin
  real ( kind = 8 ) x

  save pi2

  data pi2 / 1.57079632679489661923132169163975D+00 /

  r8_acos = pi2 - r8_asin ( x )

  return
end
function r8_acosh ( x )

!*****************************************************************************80
!
!! R8_ACOSH evaluates the arc-hyperbolic cosine of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_ACOSH, the arc-hyperbolic cosine of X.
!
  implicit none

  real ( kind = 8 ) dln2
  real ( kind = 8 ) r8_acosh
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax

  save dln2
  save xmax

  data dln2 / 0.69314718055994530941723212145818D+00 /
  data xmax / 0.0D+00 /

  if ( xmax == 0.0D+00 ) then
    xmax = 1.0D+00 / sqrt ( r8_mach ( 3 ) )
  end if

  if ( x < 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ACOSH - Fatal error!'
    write ( *, '(a)' ) '  X < 1.0'
    stop
  else if ( x < xmax ) then
    value = log ( x + sqrt ( x * x - 1.0D+00 ) )
  else
    value = dln2 + log ( x )
  end if

  r8_acosh = value

  return
end
subroutine r8_admp ( x, ampl, phi )

!*****************************************************************************80
!
!! R8_ADMP: modulus and phase of the derivative of the Airy function.
!
!  Description:
!
!    This function must only be called when X <= -1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) AMPL, PHI, the modulus and phase of the 
!    derivative of the Airy function.
!
  implicit none

  real ( kind = 8 ) ampl
  real ( kind = 8 ) an20cs(57)
  real ( kind = 8 ) an21cs(60)
  real ( kind = 8 ) an22cs(74)
  real ( kind = 8 ) aph0cs(53)
  real ( kind = 8 ) aph1cs(58)
  real ( kind = 8 ) aph2cs(72)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nan20
  integer ( kind = 4 ) nan21
  integer ( kind = 4 ) nan22
  integer ( kind = 4 ) naph0
  integer ( kind = 4 ) naph1
  integer ( kind = 4 ) naph2
  real ( kind = 8 ) phi
  real ( kind = 8 ) pi34
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sqrtx
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml
  real ( kind = 8 ) z

  save an20cs
  save an21cs
  save an22cs
  save aph0cs
  save aph1cs
  save aph2cs
  save nan20
  save nan21
  save nan22
  save naph0
  save naph1
  save naph2
  save pi34
  save xsml

  data an22cs(  1) /  0.0537418629629794329091103360917783D+00 /
  data an22cs(  2) / -0.0126661435859883193466312085036450D+00 /
  data an22cs(  3) / -0.0011924334106593006840848916913681D+00 /
  data an22cs(  4) / -0.0002032327627275654552687155176363D+00 /
  data an22cs(  5) / -0.0000446468963075163979516164905945D+00 /
  data an22cs(  6) / -0.0000113359036053123490416997893086D+00 /
  data an22cs(  7) / -0.0000031641352378546107356671355827D+00 /
  data an22cs(  8) / -0.0000009446708886148939120888532442D+00 /
  data an22cs(  9) / -0.0000002966562236471765527900905456D+00 /
  data an22cs( 10) / -0.0000000969118892024367799908661433D+00 /
  data an22cs( 11) / -0.0000000326822538653274091533072559D+00 /
  data an22cs( 12) / -0.0000000113144618963583865900447294D+00 /
  data an22cs( 13) / -0.0000000040042691001741501738278050D+00 /
  data an22cs( 14) / -0.0000000014440333683907423778522199D+00 /
  data an22cs( 15) / -0.0000000005292853746152611585663541D+00 /
  data an22cs( 16) / -0.0000000001967763373707889528245726D+00 /
  data an22cs( 17) / -0.0000000000740800095755849858816731D+00 /
  data an22cs( 18) / -0.0000000000282016314294661982842740D+00 /
  data an22cs( 19) / -0.0000000000108440066463128331337590D+00 /
  data an22cs( 20) / -0.0000000000042074800682644236920617D+00 /
  data an22cs( 21) / -0.0000000000016459149670634819724739D+00 /
  data an22cs( 22) / -0.0000000000006486826705121018896077D+00 /
  data an22cs( 23) / -0.0000000000002574095003354105832300D+00 /
  data an22cs( 24) / -0.0000000000001027889029407822132143D+00 /
  data an22cs( 25) / -0.0000000000000412845827195222720128D+00 /
  data an22cs( 26) / -0.0000000000000166711029332862509726D+00 /
  data an22cs( 27) / -0.0000000000000067656696165608023403D+00 /
  data an22cs( 28) / -0.0000000000000027585448232693576823D+00 /
  data an22cs( 29) / -0.0000000000000011296397915297168938D+00 /
  data an22cs( 30) / -0.0000000000000004644848225457314333D+00 /
  data an22cs( 31) / -0.0000000000000001917198035033912928D+00 /
  data an22cs( 32) / -0.0000000000000000794197570111893530D+00 /
  data an22cs( 33) / -0.0000000000000000330116492300368930D+00 /
  data an22cs( 34) / -0.0000000000000000137658057726549714D+00 /
  data an22cs( 35) / -0.0000000000000000057578093720012791D+00 /
  data an22cs( 36) / -0.0000000000000000024152700858632017D+00 /
  data an22cs( 37) / -0.0000000000000000010159301700933666D+00 /
  data an22cs( 38) / -0.0000000000000000004284434955330055D+00 /
  data an22cs( 39) / -0.0000000000000000001811344052168016D+00 /
  data an22cs( 40) / -0.0000000000000000000767602045619422D+00 /
  data an22cs( 41) / -0.0000000000000000000326026346758614D+00 /
  data an22cs( 42) / -0.0000000000000000000138773806682627D+00 /
  data an22cs( 43) / -0.0000000000000000000059191627103729D+00 /
  data an22cs( 44) / -0.0000000000000000000025297256431944D+00 /
  data an22cs( 45) / -0.0000000000000000000010832077293819D+00 /
  data an22cs( 46) / -0.0000000000000000000004646674880404D+00 /
  data an22cs( 47) / -0.0000000000000000000001996797783865D+00 /
  data an22cs( 48) / -0.0000000000000000000000859524108705D+00 /
  data an22cs( 49) / -0.0000000000000000000000370584152073D+00 /
  data an22cs( 50) / -0.0000000000000000000000160027503479D+00 /
  data an22cs( 51) / -0.0000000000000000000000069208124999D+00 /
  data an22cs( 52) / -0.0000000000000000000000029974448994D+00 /
  data an22cs( 53) / -0.0000000000000000000000013000356362D+00 /
  data an22cs( 54) / -0.0000000000000000000000005646100942D+00 /
  data an22cs( 55) / -0.0000000000000000000000002455341103D+00 /
  data an22cs( 56) / -0.0000000000000000000000001069119686D+00 /
  data an22cs( 57) / -0.0000000000000000000000000466095090D+00 /
  data an22cs( 58) / -0.0000000000000000000000000203441579D+00 /
  data an22cs( 59) / -0.0000000000000000000000000088900866D+00 /
  data an22cs( 60) / -0.0000000000000000000000000038891813D+00 /
  data an22cs( 61) / -0.0000000000000000000000000017032637D+00 /
  data an22cs( 62) / -0.0000000000000000000000000007467295D+00 /
  data an22cs( 63) / -0.0000000000000000000000000003277097D+00 /
  data an22cs( 64) / -0.0000000000000000000000000001439618D+00 /
  data an22cs( 65) / -0.0000000000000000000000000000633031D+00 /
  data an22cs( 66) / -0.0000000000000000000000000000278620D+00 /
  data an22cs( 67) / -0.0000000000000000000000000000122743D+00 /
  data an22cs( 68) / -0.0000000000000000000000000000054121D+00 /
  data an22cs( 69) / -0.0000000000000000000000000000023884D+00 /
  data an22cs( 70) / -0.0000000000000000000000000000010549D+00 /
  data an22cs( 71) / -0.0000000000000000000000000000004663D+00 /
  data an22cs( 72) / -0.0000000000000000000000000000002063D+00 /
  data an22cs( 73) / -0.0000000000000000000000000000000913D+00 /
  data an22cs( 74) / -0.0000000000000000000000000000000405D+00 /

  data an21cs(  1) /  0.0198313155263169394420342483165643D+00 /
  data an21cs(  2) / -0.0029376249067087533460593745594484D+00 /
  data an21cs(  3) / -0.0001136260695958195549872611137182D+00 /
  data an21cs(  4) / -0.0000100554451087156009750981645918D+00 /
  data an21cs(  5) / -0.0000013048787116563250421785598252D+00 /
  data an21cs(  6) / -0.0000002123881993150664830666079609D+00 /
  data an21cs(  7) / -0.0000000402270833384269040347850109D+00 /
  data an21cs(  8) / -0.0000000084996745953161799142201792D+00 /
  data an21cs(  9) / -0.0000000019514839426178614099532934D+00 /
  data an21cs( 10) / -0.0000000004783865343840384282992480D+00 /
  data an21cs( 11) / -0.0000000001236733992099450501137105D+00 /
  data an21cs( 12) / -0.0000000000334137486398754232219789D+00 /
  data an21cs( 13) / -0.0000000000093702823540766329897780D+00 /
  data an21cs( 14) / -0.0000000000027130128156139564687240D+00 /
  data an21cs( 15) / -0.0000000000008075953800583479535949D+00 /
  data an21cs( 16) / -0.0000000000002463214304700125252160D+00 /
  data an21cs( 17) / -0.0000000000000767655689109321564410D+00 /
  data an21cs( 18) / -0.0000000000000243882598807354919791D+00 /
  data an21cs( 19) / -0.0000000000000078831466358760308462D+00 /
  data an21cs( 20) / -0.0000000000000025882400995585864077D+00 /
  data an21cs( 21) / -0.0000000000000008619457862945690828D+00 /
  data an21cs( 22) / -0.0000000000000002907994739663128534D+00 /
  data an21cs( 23) / -0.0000000000000000992846796122890484D+00 /
  data an21cs( 24) / -0.0000000000000000342720229187774480D+00 /
  data an21cs( 25) / -0.0000000000000000119511048205515026D+00 /
  data an21cs( 26) / -0.0000000000000000042069729043678359D+00 /
  data an21cs( 27) / -0.0000000000000000014939697762818400D+00 /
  data an21cs( 28) / -0.0000000000000000005348981161589517D+00 /
  data an21cs( 29) / -0.0000000000000000001929877577826238D+00 /
  data an21cs( 30) / -0.0000000000000000000701313701018203D+00 /
  data an21cs( 31) / -0.0000000000000000000256585738509682D+00 /
  data an21cs( 32) / -0.0000000000000000000094475894562734D+00 /
  data an21cs( 33) / -0.0000000000000000000034996401941465D+00 /
  data an21cs( 34) / -0.0000000000000000000013037622466397D+00 /
  data an21cs( 35) / -0.0000000000000000000004883334163346D+00 /
  data an21cs( 36) / -0.0000000000000000000001838477586152D+00 /
  data an21cs( 37) / -0.0000000000000000000000695527324058D+00 /
  data an21cs( 38) / -0.0000000000000000000000264351910209D+00 /
  data an21cs( 39) / -0.0000000000000000000000100918094655D+00 /
  data an21cs( 40) / -0.0000000000000000000000038688924289D+00 /
  data an21cs( 41) / -0.0000000000000000000000014892036525D+00 /
  data an21cs( 42) / -0.0000000000000000000000005754342426D+00 /
  data an21cs( 43) / -0.0000000000000000000000002231725971D+00 /
  data an21cs( 44) / -0.0000000000000000000000000868607480D+00 /
  data an21cs( 45) / -0.0000000000000000000000000339220403D+00 /
  data an21cs( 46) / -0.0000000000000000000000000132910128D+00 /
  data an21cs( 47) / -0.0000000000000000000000000052239309D+00 /
  data an21cs( 48) / -0.0000000000000000000000000020594383D+00 /
  data an21cs( 49) / -0.0000000000000000000000000008142614D+00 /
  data an21cs( 50) / -0.0000000000000000000000000003228473D+00 /
  data an21cs( 51) / -0.0000000000000000000000000001283529D+00 /
  data an21cs( 52) / -0.0000000000000000000000000000511622D+00 /
  data an21cs( 53) / -0.0000000000000000000000000000204451D+00 /
  data an21cs( 54) / -0.0000000000000000000000000000081901D+00 /
  data an21cs( 55) / -0.0000000000000000000000000000032886D+00 /
  data an21cs( 56) / -0.0000000000000000000000000000013235D+00 /
  data an21cs( 57) / -0.0000000000000000000000000000005338D+00 /
  data an21cs( 58) / -0.0000000000000000000000000000002158D+00 /
  data an21cs( 59) / -0.0000000000000000000000000000000874D+00 /
  data an21cs( 60) / -0.0000000000000000000000000000000355D+00 /

  data an20cs(  1) /  0.0126732217145738027154610751034240D+00 /
  data an20cs(  2) / -0.0005212847072615621184780942309478D+00 /
  data an20cs(  3) / -0.0000052672111140370429809074052969D+00 /
  data an20cs(  4) / -0.0000001628202185026483752632460680D+00 /
  data an20cs(  5) / -0.0000000090991442687371386325973075D+00 /
  data an20cs(  6) / -0.0000000007438647126242192890685403D+00 /
  data an20cs(  7) / -0.0000000000795494751591469486122822D+00 /
  data an20cs(  8) / -0.0000000000104050944288303742803960D+00 /
  data an20cs(  9) / -0.0000000000015932425598414551523990D+00 /
  data an20cs( 10) / -0.0000000000002770648272341913946674D+00 /
  data an20cs( 11) / -0.0000000000000535342629237606295104D+00 /
  data an20cs( 12) / -0.0000000000000113061541781728314051D+00 /
  data an20cs( 13) / -0.0000000000000025772190078943167788D+00 /
  data an20cs( 14) / -0.0000000000000006278033116032485076D+00 /
  data an20cs( 15) / -0.0000000000000001621295400189939757D+00 /
  data an20cs( 16) / -0.0000000000000000440992985240675353D+00 /
  data an20cs( 17) / -0.0000000000000000125655516553258972D+00 /
  data an20cs( 18) / -0.0000000000000000037336906988015204D+00 /
  data an20cs( 19) / -0.0000000000000000011524626926724671D+00 /
  data an20cs( 20) / -0.0000000000000000003683081499099144D+00 /
  data an20cs( 21) / -0.0000000000000000001215206965331797D+00 /
  data an20cs( 22) / -0.0000000000000000000412916177724016D+00 /
  data an20cs( 23) / -0.0000000000000000000144177364239347D+00 /
  data an20cs( 24) / -0.0000000000000000000051631842875864D+00 /
  data an20cs( 25) / -0.0000000000000000000018931242668250D+00 /
  data an20cs( 26) / -0.0000000000000000000007096054668569D+00 /
  data an20cs( 27) / -0.0000000000000000000002715406646904D+00 /
  data an20cs( 28) / -0.0000000000000000000001059486979400D+00 /
  data an20cs( 29) / -0.0000000000000000000000421030035685D+00 /
  data an20cs( 30) / -0.0000000000000000000000170233781664D+00 /
  data an20cs( 31) / -0.0000000000000000000000069966677028D+00 /
  data an20cs( 32) / -0.0000000000000000000000029206643813D+00 /
  data an20cs( 33) / -0.0000000000000000000000012373128203D+00 /
  data an20cs( 34) / -0.0000000000000000000000005315871095D+00 /
  data an20cs( 35) / -0.0000000000000000000000002314622618D+00 /
  data an20cs( 36) / -0.0000000000000000000000001020779922D+00 /
  data an20cs( 37) / -0.0000000000000000000000000455706227D+00 /
  data an20cs( 38) / -0.0000000000000000000000000205831071D+00 /
  data an20cs( 39) / -0.0000000000000000000000000094015189D+00 /
  data an20cs( 40) / -0.0000000000000000000000000043405874D+00 /
  data an20cs( 41) / -0.0000000000000000000000000020247792D+00 /
  data an20cs( 42) / -0.0000000000000000000000000009539214D+00 /
  data an20cs( 43) / -0.0000000000000000000000000004537234D+00 /
  data an20cs( 44) / -0.0000000000000000000000000002178016D+00 /
  data an20cs( 45) / -0.0000000000000000000000000001054823D+00 /
  data an20cs( 46) / -0.0000000000000000000000000000515242D+00 /
  data an20cs( 47) / -0.0000000000000000000000000000253763D+00 /
  data an20cs( 48) / -0.0000000000000000000000000000125983D+00 /
  data an20cs( 49) / -0.0000000000000000000000000000063030D+00 /
  data an20cs( 50) / -0.0000000000000000000000000000031771D+00 /
  data an20cs( 51) / -0.0000000000000000000000000000016131D+00 /
  data an20cs( 52) / -0.0000000000000000000000000000008248D+00 /
  data an20cs( 53) / -0.0000000000000000000000000000004246D+00 /
  data an20cs( 54) / -0.0000000000000000000000000000002200D+00 /
  data an20cs( 55) / -0.0000000000000000000000000000001147D+00 /
  data an20cs( 56) / -0.0000000000000000000000000000000602D+00 /
  data an20cs( 57) / -0.0000000000000000000000000000000318D+00 /

  data aph2cs(  1) / -0.2057088719781465106973648665602125D+00 /
  data aph2cs(  2) /  0.0422196961357771921673114980369460D+00 /
  data aph2cs(  3) /  0.0020482560511207275042660577813334D+00 /
  data aph2cs(  4) /  0.0002607800735165005631187879922652D+00 /
  data aph2cs(  5) /  0.0000474824268004728875381750519293D+00 /
  data aph2cs(  6) /  0.0000105102756431611743473630026955D+00 /
  data aph2cs(  7) /  0.0000026353534014667945109314041983D+00 /
  data aph2cs(  8) /  0.0000007208824863499147299790783731D+00 /
  data aph2cs(  9) /  0.0000002103236664473352859749477082D+00 /
  data aph2cs( 10) /  0.0000000644975634555295598437362273D+00 /
  data aph2cs( 11) /  0.0000000205802377264368507978116888D+00 /
  data aph2cs( 12) /  0.0000000067836273920906428963513918D+00 /
  data aph2cs( 13) /  0.0000000022974015284009400168343792D+00 /
  data aph2cs( 14) /  0.0000000007961306765491187534883226D+00 /
  data aph2cs( 15) /  0.0000000002813860609741591719003632D+00 /
  data aph2cs( 16) /  0.0000000001011749056931973922841793D+00 /
  data aph2cs( 17) /  0.0000000000369306737952476559097060D+00 /
  data aph2cs( 18) /  0.0000000000136615066127098031778842D+00 /
  data aph2cs( 19) /  0.0000000000051142751416045045119388D+00 /
  data aph2cs( 20) /  0.0000000000019351688931706516247975D+00 /
  data aph2cs( 21) /  0.0000000000007393606916493224217271D+00 /
  data aph2cs( 22) /  0.0000000000002849792219222743597555D+00 /
  data aph2cs( 23) /  0.0000000000001107280782459648335733D+00 /
  data aph2cs( 24) /  0.0000000000000433412199370134633169D+00 /
  data aph2cs( 25) /  0.0000000000000170800825265670367471D+00 /
  data aph2cs( 26) /  0.0000000000000067733080195631114673D+00 /
  data aph2cs( 27) /  0.0000000000000027016904789262414108D+00 /
  data aph2cs( 28) /  0.0000000000000010834720751810782141D+00 /
  data aph2cs( 29) /  0.0000000000000004367060312970286167D+00 /
  data aph2cs( 30) /  0.0000000000000001768511738053366608D+00 /
  data aph2cs( 31) /  0.0000000000000000719359213093645717D+00 /
  data aph2cs( 32) /  0.0000000000000000293823610002933154D+00 /
  data aph2cs( 33) /  0.0000000000000000120482811525848357D+00 /
  data aph2cs( 34) /  0.0000000000000000049586659491091389D+00 /
  data aph2cs( 35) /  0.0000000000000000020479438315847217D+00 /
  data aph2cs( 36) /  0.0000000000000000008486019944410629D+00 /
  data aph2cs( 37) /  0.0000000000000000003527351765384506D+00 /
  data aph2cs( 38) /  0.0000000000000000001470563996804903D+00 /
  data aph2cs( 39) /  0.0000000000000000000614817826902188D+00 /
  data aph2cs( 40) /  0.0000000000000000000257737706565077D+00 /
  data aph2cs( 41) /  0.0000000000000000000108323903590042D+00 /
  data aph2cs( 42) /  0.0000000000000000000045638898024998D+00 /
  data aph2cs( 43) /  0.0000000000000000000019273635403662D+00 /
  data aph2cs( 44) /  0.0000000000000000000008157668569775D+00 /
  data aph2cs( 45) /  0.0000000000000000000003460202828346D+00 /
  data aph2cs( 46) /  0.0000000000000000000001470726482427D+00 /
  data aph2cs( 47) /  0.0000000000000000000000626356074088D+00 /
  data aph2cs( 48) /  0.0000000000000000000000267261292780D+00 /
  data aph2cs( 49) /  0.0000000000000000000000114246948763D+00 /
  data aph2cs( 50) /  0.0000000000000000000000048923460516D+00 /
  data aph2cs( 51) /  0.0000000000000000000000020985807810D+00 /
  data aph2cs( 52) /  0.0000000000000000000000009016618807D+00 /
  data aph2cs( 53) /  0.0000000000000000000000003880129464D+00 /
  data aph2cs( 54) /  0.0000000000000000000000001672282170D+00 /
  data aph2cs( 55) /  0.0000000000000000000000000721790800D+00 /
  data aph2cs( 56) /  0.0000000000000000000000000311982573D+00 /
  data aph2cs( 57) /  0.0000000000000000000000000135035015D+00 /
  data aph2cs( 58) /  0.0000000000000000000000000058524861D+00 /
  data aph2cs( 59) /  0.0000000000000000000000000025397686D+00 /
  data aph2cs( 60) /  0.0000000000000000000000000011035457D+00 /
  data aph2cs( 61) /  0.0000000000000000000000000004800788D+00 /
  data aph2cs( 62) /  0.0000000000000000000000000002090956D+00 /
  data aph2cs( 63) /  0.0000000000000000000000000000911743D+00 /
  data aph2cs( 64) /  0.0000000000000000000000000000397998D+00 /
  data aph2cs( 65) /  0.0000000000000000000000000000173923D+00 /
  data aph2cs( 66) /  0.0000000000000000000000000000076083D+00 /
  data aph2cs( 67) /  0.0000000000000000000000000000033316D+00 /
  data aph2cs( 68) /  0.0000000000000000000000000000014604D+00 /
  data aph2cs( 69) /  0.0000000000000000000000000000006407D+00 /
  data aph2cs( 70) /  0.0000000000000000000000000000002814D+00 /
  data aph2cs( 71) /  0.0000000000000000000000000000001237D+00 /
  data aph2cs( 72) /  0.0000000000000000000000000000000544D+00 /

  data aph1cs(  1) / -0.1024172908077571694021123321813917D+00 /
  data aph1cs(  2) /  0.0071697275146591248047211649144704D+00 /
  data aph1cs(  3) /  0.0001209959363122328589813856491397D+00 /
  data aph1cs(  4) /  0.0000073361512841219912080297845684D+00 /
  data aph1cs(  5) /  0.0000007535382954271607069982903869D+00 /
  data aph1cs(  6) /  0.0000001041478171741301926885109155D+00 /
  data aph1cs(  7) /  0.0000000174358728518545691858907606D+00 /
  data aph1cs(  8) /  0.0000000033399795033346451660184961D+00 /
  data aph1cs(  9) /  0.0000000007073075174363527083399508D+00 /
  data aph1cs( 10) /  0.0000000001619187515189773266792272D+00 /
  data aph1cs( 11) /  0.0000000000394539981881954889879668D+00 /
  data aph1cs( 12) /  0.0000000000101192281734227133292631D+00 /
  data aph1cs( 13) /  0.0000000000027092778259520332198030D+00 /
  data aph1cs( 14) /  0.0000000000007523806418422548885854D+00 /
  data aph1cs( 15) /  0.0000000000002156368733008966357328D+00 /
  data aph1cs( 16) /  0.0000000000000635282777126068410174D+00 /
  data aph1cs( 17) /  0.0000000000000191756972641501729345D+00 /
  data aph1cs( 18) /  0.0000000000000059143072446464891558D+00 /
  data aph1cs( 19) /  0.0000000000000018597128517275028357D+00 /
  data aph1cs( 20) /  0.0000000000000005950444923946103668D+00 /
  data aph1cs( 21) /  0.0000000000000001934229956430180252D+00 /
  data aph1cs( 22) /  0.0000000000000000637843021489504324D+00 /
  data aph1cs( 23) /  0.0000000000000000213127290087312393D+00 /
  data aph1cs( 24) /  0.0000000000000000072081380656728500D+00 /
  data aph1cs( 25) /  0.0000000000000000024652494144769247D+00 /
  data aph1cs( 26) /  0.0000000000000000008519110570266154D+00 /
  data aph1cs( 27) /  0.0000000000000000002972384468491170D+00 /
  data aph1cs( 28) /  0.0000000000000000001046426648811446D+00 /
  data aph1cs( 29) /  0.0000000000000000000371493036347327D+00 /
  data aph1cs( 30) /  0.0000000000000000000132923247793472D+00 /
  data aph1cs( 31) /  0.0000000000000000000047912837925909D+00 /
  data aph1cs( 32) /  0.0000000000000000000017390619859336D+00 /
  data aph1cs( 33) /  0.0000000000000000000006353585173501D+00 /
  data aph1cs( 34) /  0.0000000000000000000002335643614263D+00 /
  data aph1cs( 35) /  0.0000000000000000000000863643881606D+00 /
  data aph1cs( 36) /  0.0000000000000000000000321123006944D+00 /
  data aph1cs( 37) /  0.0000000000000000000000120031540983D+00 /
  data aph1cs( 38) /  0.0000000000000000000000045091488699D+00 /
  data aph1cs( 39) /  0.0000000000000000000000017020228580D+00 /
  data aph1cs( 40) /  0.0000000000000000000000006453744630D+00 /
  data aph1cs( 41) /  0.0000000000000000000000002457788564D+00 /
  data aph1cs( 42) /  0.0000000000000000000000000939897684D+00 /
  data aph1cs( 43) /  0.0000000000000000000000000360863150D+00 /
  data aph1cs( 44) /  0.0000000000000000000000000139077884D+00 /
  data aph1cs( 45) /  0.0000000000000000000000000053797184D+00 /
  data aph1cs( 46) /  0.0000000000000000000000000020882551D+00 /
  data aph1cs( 47) /  0.0000000000000000000000000008133371D+00 /
  data aph1cs( 48) /  0.0000000000000000000000000003178080D+00 /
  data aph1cs( 49) /  0.0000000000000000000000000001245700D+00 /
  data aph1cs( 50) /  0.0000000000000000000000000000489742D+00 /
  data aph1cs( 51) /  0.0000000000000000000000000000193099D+00 /
  data aph1cs( 52) /  0.0000000000000000000000000000076349D+00 /
  data aph1cs( 53) /  0.0000000000000000000000000000030269D+00 /
  data aph1cs( 54) /  0.0000000000000000000000000000012032D+00 /
  data aph1cs( 55) /  0.0000000000000000000000000000004795D+00 /
  data aph1cs( 56) /  0.0000000000000000000000000000001915D+00 /
  data aph1cs( 57) /  0.0000000000000000000000000000000767D+00 /
  data aph1cs( 58) /  0.0000000000000000000000000000000308D+00 /

  data aph0cs(  1) / -0.0855849241130933256920124260179491D+00 /
  data aph0cs(  2) /  0.0011214378867065260735786722471124D+00 /
  data aph0cs(  3) /  0.0000042721029353664113951573742015D+00 /
  data aph0cs(  4) /  0.0000000817607381483243644018062323D+00 /
  data aph0cs(  5) /  0.0000000033907645000492724207816418D+00 /
  data aph0cs(  6) /  0.0000000002253264422619113939845276D+00 /
  data aph0cs(  7) /  0.0000000000206284209229015251256990D+00 /
  data aph0cs(  8) /  0.0000000000023858762828130887627258D+00 /
  data aph0cs(  9) /  0.0000000000003301618105886705480628D+00 /
  data aph0cs( 10) /  0.0000000000000527009648508328581123D+00 /
  data aph0cs( 11) /  0.0000000000000094555482203813492868D+00 /
  data aph0cs( 12) /  0.0000000000000018709426951344836908D+00 /
  data aph0cs( 13) /  0.0000000000000004023980041825392741D+00 /
  data aph0cs( 14) /  0.0000000000000000930192879258983167D+00 /
  data aph0cs( 15) /  0.0000000000000000229038635402379945D+00 /
  data aph0cs( 16) /  0.0000000000000000059634359822083386D+00 /
  data aph0cs( 17) /  0.0000000000000000016320279659403399D+00 /
  data aph0cs( 18) /  0.0000000000000000004671145658861339D+00 /
  data aph0cs( 19) /  0.0000000000000000001392334415363502D+00 /
  data aph0cs( 20) /  0.0000000000000000000430642670285155D+00 /
  data aph0cs( 21) /  0.0000000000000000000137781416318755D+00 /
  data aph0cs( 22) /  0.0000000000000000000045476710480396D+00 /
  data aph0cs( 23) /  0.0000000000000000000015448420203026D+00 /
  data aph0cs( 24) /  0.0000000000000000000005389770551212D+00 /
  data aph0cs( 25) /  0.0000000000000000000001927726737155D+00 /
  data aph0cs( 26) /  0.0000000000000000000000705659320166D+00 /
  data aph0cs( 27) /  0.0000000000000000000000263985084827D+00 /
  data aph0cs( 28) /  0.0000000000000000000000100791301805D+00 /
  data aph0cs( 29) /  0.0000000000000000000000039228928481D+00 /
  data aph0cs( 30) /  0.0000000000000000000000015547422955D+00 /
  data aph0cs( 31) /  0.0000000000000000000000006268306372D+00 /
  data aph0cs( 32) /  0.0000000000000000000000002568563962D+00 /
  data aph0cs( 33) /  0.0000000000000000000000001068858883D+00 /
  data aph0cs( 34) /  0.0000000000000000000000000451347253D+00 /
  data aph0cs( 35) /  0.0000000000000000000000000193267262D+00 /
  data aph0cs( 36) /  0.0000000000000000000000000083865369D+00 /
  data aph0cs( 37) /  0.0000000000000000000000000036857386D+00 /
  data aph0cs( 38) /  0.0000000000000000000000000016396202D+00 /
  data aph0cs( 39) /  0.0000000000000000000000000007379298D+00 /
  data aph0cs( 40) /  0.0000000000000000000000000003358392D+00 /
  data aph0cs( 41) /  0.0000000000000000000000000001544891D+00 /
  data aph0cs( 42) /  0.0000000000000000000000000000718013D+00 /
  data aph0cs( 43) /  0.0000000000000000000000000000337026D+00 /
  data aph0cs( 44) /  0.0000000000000000000000000000159710D+00 /
  data aph0cs( 45) /  0.0000000000000000000000000000076382D+00 /
  data aph0cs( 46) /  0.0000000000000000000000000000036855D+00 /
  data aph0cs( 47) /  0.0000000000000000000000000000017935D+00 /
  data aph0cs( 48) /  0.0000000000000000000000000000008800D+00 /
  data aph0cs( 49) /  0.0000000000000000000000000000004353D+00 /
  data aph0cs( 50) /  0.0000000000000000000000000000002170D+00 /
  data aph0cs( 51) /  0.0000000000000000000000000000001090D+00 /
  data aph0cs( 52) /  0.0000000000000000000000000000000551D+00 /
  data aph0cs( 53) /  0.0000000000000000000000000000000281D+00 /

  data nan20 / 0 /
  data nan21 / 0 /
  data nan22 / 0 /
  data naph0 / 0 /
  data naph1 / 0 /
  data naph2 / 0 /
  data pi34 / 2.35619449019234492884698253745962716313D+00  /
  data xsml / 0.0D+00 /

  if ( nan20 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nan20 = r8_inits ( an20cs, 57, eta )
    nan21 = r8_inits ( an21cs, 60, eta )
    nan22 = r8_inits ( an22cs, 74, eta )
    naph0 = r8_inits ( aph0cs, 53, eta )
    naph1 = r8_inits ( aph1cs, 58, eta )
    naph2 = r8_inits ( aph2cs, 72, eta )
    xsml = - ( 128.0D+00 / r8_mach ( 3 ) ) ** 0.3333D+00
  end if

  if ( x < xsml ) then
    z = 1.0D+00
    ampl = 0.3125D+00 + r8_csevl ( z, an20cs, nan20 )
    phi = - 0.625D+00 + r8_csevl ( z, aph0cs, naph0 )
  else if ( x < - 4.0D+00 ) then
    z = 128.0D+00 / x / x / x + 1.0D+00
    ampl = 0.3125D+00 + r8_csevl ( z, an20cs, nan20 )
    phi = - 0.625D+00 + r8_csevl ( z, aph0cs, naph0 )
  else if ( x < - 2.0D+00 ) then
    z = ( 128.0D+00 / x / x / x + 9.0D+00 ) / 7.0D+00
    ampl = 0.3125D+00 + r8_csevl ( z, an21cs, nan21 )
    phi = - 0.625D+00 + r8_csevl ( z, aph1cs, naph1 )
  else if ( x <= - 1.0D+00 ) then
    z = ( 16.0D+00 / x / x / x + 9.0D+00 ) / 7.0D+00
    ampl = 0.3125D+00 + r8_csevl ( z, an22cs, nan22 )
    phi = - 0.625D+00 + r8_csevl ( z, aph2cs, naph2 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ADMP - Fatal error!'
    write ( *, '(a)' ) '  - 1.0 < X.'
    stop
  end if

  sqrtx = sqrt ( - x )
  ampl = sqrt ( ampl * sqrtx )
  phi = pi34 - x * sqrtx * phi

  return
end
function r8_ai ( x )

!*****************************************************************************80
!
!! R8_AI evaluates the Airy function Ai of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_AI, the Airy function Ai of X.
!
  implicit none

  real ( kind = 8 ) aifcs(13)
  real ( kind = 8 ) aigcs(13)
  integer ( kind = 4 ) naif
  integer ( kind = 4 ) naig
  real ( kind = 8 ) r8_ai
  real ( kind = 8 ) r8_aie
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) x3sml
  real ( kind = 8 ) xm
  real ( kind = 8 ) xmax
  real ( kind = 8 ) z

  save aifcs
  save aigcs
  save naif
  save naig
  save x3sml
  save xmax

  data aifcs(  1) / -0.37971358496669997496197089469414D-01 /
  data aifcs(  2) / +0.59191888537263638574319728013777D-01 /
  data aifcs(  3) / +0.98629280577279975365603891044060D-03 /
  data aifcs(  4) / +0.68488438190765667554854830182412D-05 /
  data aifcs(  5) / +0.25942025962194713019489279081403D-07 /
  data aifcs(  6) / +0.61766127740813750329445749697236D-10 /
  data aifcs(  7) / +0.10092454172466117901429556224601D-12 /
  data aifcs(  8) / +0.12014792511179938141288033225333D-15 /
  data aifcs(  9) / +0.10882945588716991878525295466666D-18 /
  data aifcs( 10) / +0.77513772196684887039238400000000D-22 /
  data aifcs( 11) / +0.44548112037175638391466666666666D-25 /
  data aifcs( 12) / +0.21092845231692343466666666666666D-28 /
  data aifcs( 13) / +0.83701735910741333333333333333333D-32 /

  data aigcs(  1) / +0.18152365581161273011556209957864D-01 /
  data aigcs(  2) / +0.21572563166010755534030638819968D-01 /
  data aigcs(  3) / +0.25678356987483249659052428090133D-03 /
  data aigcs(  4) / +0.14265214119792403898829496921721D-05 /
  data aigcs(  5) / +0.45721149200180426070434097558191D-08 /
  data aigcs(  6) / +0.95251708435647098607392278840592D-11 /
  data aigcs(  7) / +0.13925634605771399051150420686190D-13 /
  data aigcs(  8) / +0.15070999142762379592306991138666D-16 /
  data aigcs(  9) / +0.12559148312567778822703205333333D-19 /
  data aigcs( 10) / +0.83063073770821340343829333333333D-23 /
  data aigcs( 11) / +0.44657538493718567445333333333333D-26 /
  data aigcs( 12) / +0.19900855034518869333333333333333D-29 /
  data aigcs( 13) / +0.74702885256533333333333333333333D-33 /

  data naif / 0 /
  data naig / 0 /
  data x3sml / 0.0D+00 /
  data xmax / 0.0D+00 /

  if ( naif == 0 ) then
    naif = r8_inits ( aifcs, 13, 0.1D+00 * r8_mach ( 3 ) )
    naig = r8_inits ( aigcs, 13, 0.1D+00 * r8_mach ( 3 ) )
    x3sml = r8_mach ( 3 ) ** 0.3334D+00
    xmax = ( - 1.5D+00 * log ( r8_mach ( 1 ) ) ) ** 0.6667D+00
    xmax = xmax - xmax * log ( xmax ) / &
      ( 4.0D+00 * xmax * sqrt ( xmax ) + 1.0D+00 ) - 0.01D+00
  end if

  if ( x < - 1.0D+00 ) then
    call r8_aimp ( x, xm, theta )
    r8_ai = xm * cos ( theta )
  else if ( abs ( x ) <= x3sml ) then
    z = 0.0D+00
    r8_ai = 0.375D+00 + ( r8_csevl ( z, aifcs, naif ) &
      - x * ( 0.25D+00 + r8_csevl ( z, aigcs, naig ) ) )
  else if ( x <= 1.0D+00 ) then
    z = x * x * x
    r8_ai = 0.375D+00 + ( r8_csevl ( z, aifcs, naif ) &
      - x * ( 0.25D+00 + r8_csevl ( z, aigcs, naig ) ) )
  else if ( x <= xmax ) then
    r8_ai = r8_aie ( x ) &
      * exp ( - 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
  else
    r8_ai = 0.0D+00
  end if

  return
end
function r8_aid ( x )

!*****************************************************************************80
!
!! R8_AID evaluates the derivative of the Airy function Ai of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_AID, the derivative of the Airy function 
!    Ai of X.
!
  implicit none

  real ( kind = 8 ) aifcs(13)
  real ( kind = 8 ) aigcs(13)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) naif
  integer ( kind = 4 ) naig
  real ( kind = 8 ) phi
  real ( kind = 8 ) r8_aid
  real ( kind = 8 ) r8_aide
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) x2sml
  real ( kind = 8 ) x3
  real ( kind = 8 ) x3sml
  real ( kind = 8 ) xn

  save aifcs
  save aigcs
  save naif
  save naig
  save x2sml
  save x3sml

  data aifcs(  1) /  0.105274612265314088088970057325134114D+00 /
  data aifcs(  2) /  0.011836136281529978442889292583980840D+00 /
  data aifcs(  3) /  0.000123281041732256643051689242469164D+00 /
  data aifcs(  4) /  0.000000622612256381399016825658693579D+00 /
  data aifcs(  5) /  0.000000001852988878441452950548140821D+00 /
  data aifcs(  6) /  0.000000000003633288725904357915995625D+00 /
  data aifcs(  7) /  0.000000000000005046217040440664768330D+00 /
  data aifcs(  8) /  0.000000000000000005223816555471480985D+00 /
  data aifcs(  9) /  0.000000000000000000004185745090748989D+00 /
  data aifcs( 10) /  0.000000000000000000000002672887324883D+00 /
  data aifcs( 11) /  0.000000000000000000000000001392128006D+00 /
  data aifcs( 12) /  0.000000000000000000000000000000602653D+00 /
  data aifcs( 13) /  0.000000000000000000000000000000000220D+00 /

  data aigcs(  1) /  0.0212338781509186668523122276848937D+00 /
  data aigcs(  2) /  0.0863159303352144067524942809461604D+00 /
  data aigcs(  3) /  0.0017975947203832313578033963225230D+00 /
  data aigcs(  4) /  0.0000142654998755506932526620687495D+00 /
  data aigcs(  5) /  0.0000000594379952836832010488787064D+00 /
  data aigcs(  6) /  0.0000000001524033664794478945214786D+00 /
  data aigcs(  7) /  0.0000000000002645876603490435305100D+00 /
  data aigcs(  8) /  0.0000000000000003315624296815020591D+00 /
  data aigcs(  9) /  0.0000000000000000003139789757594792D+00 /
  data aigcs( 10) /  0.0000000000000000000002325767379040D+00 /
  data aigcs( 11) /  0.0000000000000000000000001384384231D+00 /
  data aigcs( 12) /  0.0000000000000000000000000000676629D+00 /
  data aigcs( 13) /  0.0000000000000000000000000000000276D+00 /

  data naif / 0 /
  data naig / 0 /
  data x2sml / 0.0D+00 /
  data x3sml / 0.0D+00 /

  if ( naif == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    naif = r8_inits ( aifcs, 13, eta )
    naig = r8_inits ( aigcs, 13, eta )
    x3sml = r8_mach ( 3 ) ** 0.3334D+00
    x2sml = sqrt ( r8_mach ( 3 ) )
  end if

  if ( x < - 1.0D+00 ) then
    call r8_admp ( x, xn, phi )
    r8_aid = xn * cos ( phi )
  else if ( abs ( x ) <= x2sml ) then
    x2 = 0.0D+00
    x3 = 0.0D+00
    r8_aid = ( x2 * ( 0.125D+00 + r8_csevl ( x3, aifcs, naif ) ) &
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25D+00
  else if ( abs ( x ) <= x3sml ) then
    x2 = x * x
    x3 = 0.0D+00
    r8_aid = ( x2 * ( 0.125D+00 + r8_csevl ( x3, aifcs, naif ) ) &
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25D+00
  else if ( x <= 1.0D+00 ) then
    x2 = x * x
    x3 = x * x * x
    r8_aid = ( x2 * ( 0.125D+00 + r8_csevl ( x3, aifcs, naif ) ) &
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25D+00
  else
    r8_aid = r8_aide ( x ) &
      * exp ( - 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
  end if

  return
end
function r8_aide ( x )

!*****************************************************************************80
!
!! R8_AIDE: exponentially scaled derivative, Airy function Ai of an R8 argument.
!
!  Discussion:
!
!    if X <= 0,
!      R8_AIDE ( X ) = R8_AID ( X )
!    else
!      R8_AIDE ( X ) = R8_AID ( X ) * exp ( 2/3 * X^(3/2) )
!
!    Thanks to Aleksandra Piper for pointing out a correction involving 
!    the computation of Z in the last two cases, 02 February 2012.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 2012
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_AIDE, the exponentially scaled derivative of 
!    the Airy function Ai of X.
!
  implicit none

  real ( kind = 8 ) aifcs(13)
  real ( kind = 8 ) aigcs(13)
  real ( kind = 8 ) aip1cs(57)
  real ( kind = 8 ) aip2cs(37)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) naif
  integer ( kind = 4 ) naig
  integer ( kind = 4 ) naip1
  integer ( kind = 4 ) naip2
  real ( kind = 8 ) phi
  real ( kind = 8 ) r8_aide
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sqrtx
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) x2sml
  real ( kind = 8 ) x3
  real ( kind = 8 ) x32sml
  real ( kind = 8 ) x3sml
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xn
  real ( kind = 8 ) z

  save aifcs
  save aigcs
  save aip1cs
  save aip2cs
  save naif
  save naig
  save naip1
  save naip2
  save x2sml
  save x32sml
  save x3sml
  save xbig

  data aifcs(  1) /  0.105274612265314088088970057325134114D+00 /
  data aifcs(  2) /  0.011836136281529978442889292583980840D+00 /
  data aifcs(  3) /  0.000123281041732256643051689242469164D+00 /
  data aifcs(  4) /  0.000000622612256381399016825658693579D+00 /
  data aifcs(  5) /  0.000000001852988878441452950548140821D+00 /
  data aifcs(  6) /  0.000000000003633288725904357915995625D+00 /
  data aifcs(  7) /  0.000000000000005046217040440664768330D+00 /
  data aifcs(  8) /  0.000000000000000005223816555471480985D+00 /
  data aifcs(  9) /  0.000000000000000000004185745090748989D+00 /
  data aifcs( 10) /  0.000000000000000000000002672887324883D+00 /
  data aifcs( 11) /  0.000000000000000000000000001392128006D+00 /
  data aifcs( 12) /  0.000000000000000000000000000000602653D+00 /
  data aifcs( 13) /  0.000000000000000000000000000000000220D+00 /

  data aigcs(  1) /  0.0212338781509186668523122276848937D+00 /
  data aigcs(  2) /  0.0863159303352144067524942809461604D+00 /
  data aigcs(  3) /  0.0017975947203832313578033963225230D+00 /
  data aigcs(  4) /  0.0000142654998755506932526620687495D+00 /
  data aigcs(  5) /  0.0000000594379952836832010488787064D+00 /
  data aigcs(  6) /  0.0000000001524033664794478945214786D+00 /
  data aigcs(  7) /  0.0000000000002645876603490435305100D+00 /
  data aigcs(  8) /  0.0000000000000003315624296815020591D+00 /
  data aigcs(  9) /  0.0000000000000000003139789757594792D+00 /
  data aigcs( 10) /  0.0000000000000000000002325767379040D+00 /
  data aigcs( 11) /  0.0000000000000000000000001384384231D+00 /
  data aigcs( 12) /  0.0000000000000000000000000000676629D+00 /
  data aigcs( 13) /  0.0000000000000000000000000000000276D+00 /

  data aip2cs(  1) /  0.0065457691989713756794276979067064D+00 /
  data aip2cs(  2) /  0.0023833724120774591992772552886923D+00 /
  data aip2cs(  3) / -0.0000430700770220585862775012110584D+00 /
  data aip2cs(  4) /  0.0000015629125858629202330785369063D+00 /
  data aip2cs(  5) / -0.0000000815417186162706965112501015D+00 /
  data aip2cs(  6) /  0.0000000054103738056935918208008783D+00 /
  data aip2cs(  7) / -0.0000000004284130882614696528766222D+00 /
  data aip2cs(  8) /  0.0000000000389497962832286424862198D+00 /
  data aip2cs(  9) / -0.0000000000039623161264979257658071D+00 /
  data aip2cs( 10) /  0.0000000000004428184214405989602353D+00 /
  data aip2cs( 11) / -0.0000000000000536296527150689675318D+00 /
  data aip2cs( 12) /  0.0000000000000069649872139936028200D+00 /
  data aip2cs( 13) / -0.0000000000000009619636286095319210D+00 /
  data aip2cs( 14) /  0.0000000000000001403454967784808032D+00 /
  data aip2cs( 15) / -0.0000000000000000215097136525875715D+00 /
  data aip2cs( 16) /  0.0000000000000000034471230632678283D+00 /
  data aip2cs( 17) / -0.0000000000000000005753907621819442D+00 /
  data aip2cs( 18) /  0.0000000000000000000997001165824168D+00 /
  data aip2cs( 19) / -0.0000000000000000000178811436021458D+00 /
  data aip2cs( 20) /  0.0000000000000000000033110307923551D+00 /
  data aip2cs( 21) / -0.0000000000000000000006315885529506D+00 /
  data aip2cs( 22) /  0.0000000000000000000001238666952364D+00 /
  data aip2cs( 23) / -0.0000000000000000000000249324053394D+00 /
  data aip2cs( 24) /  0.0000000000000000000000051426030999D+00 /
  data aip2cs( 25) / -0.0000000000000000000000010854236402D+00 /
  data aip2cs( 26) /  0.0000000000000000000000002341316852D+00 /
  data aip2cs( 27) / -0.0000000000000000000000000515542099D+00 /
  data aip2cs( 28) /  0.0000000000000000000000000115758841D+00 /
  data aip2cs( 29) / -0.0000000000000000000000000026479669D+00 /
  data aip2cs( 30) /  0.0000000000000000000000000006165328D+00 /
  data aip2cs( 31) / -0.0000000000000000000000000001459931D+00 /
  data aip2cs( 32) /  0.0000000000000000000000000000351331D+00 /
  data aip2cs( 33) / -0.0000000000000000000000000000085863D+00 /
  data aip2cs( 34) /  0.0000000000000000000000000000021297D+00 /
  data aip2cs( 35) / -0.0000000000000000000000000000005358D+00 /
  data aip2cs( 36) /  0.0000000000000000000000000000001367D+00 /
  data aip2cs( 37) / -0.0000000000000000000000000000000353D+00 /

  data aip1cs(  1) /  0.0358865097808301537956710489261688D+00 /
  data aip1cs(  2) /  0.0114668575627764898572700883121766D+00 /
  data aip1cs(  3) / -0.0007592073583861400301335647601603D+00 /
  data aip1cs(  4) /  0.0000869517610893841271948619434021D+00 /
  data aip1cs(  5) / -0.0000128237294298591691789607600486D+00 /
  data aip1cs(  6) /  0.0000022062695681038336934376250420D+00 /
  data aip1cs(  7) / -0.0000004222295185920749486945988432D+00 /
  data aip1cs(  8) /  0.0000000874686415726348479356130376D+00 /
  data aip1cs(  9) / -0.0000000192773588418365388625693417D+00 /
  data aip1cs( 10) /  0.0000000044668460054492719699777137D+00 /
  data aip1cs( 11) / -0.0000000010790108051948168015747466D+00 /
  data aip1cs( 12) /  0.0000000002700029446696248083071434D+00 /
  data aip1cs( 13) / -0.0000000000696480108007915257318929D+00 /
  data aip1cs( 14) /  0.0000000000184489907003246687076806D+00 /
  data aip1cs( 15) / -0.0000000000050027817358071698301149D+00 /
  data aip1cs( 16) /  0.0000000000013852243366012168297298D+00 /
  data aip1cs( 17) / -0.0000000000003908218466657048253473D+00 /
  data aip1cs( 18) /  0.0000000000001121536072524563451273D+00 /
  data aip1cs( 19) / -0.0000000000000326861522579502522443D+00 /
  data aip1cs( 20) /  0.0000000000000096619179010090805752D+00 /
  data aip1cs( 21) / -0.0000000000000028934767442698434271D+00 /
  data aip1cs( 22) /  0.0000000000000008770086661150897069D+00 /
  data aip1cs( 23) / -0.0000000000000002688046261195853754D+00 /
  data aip1cs( 24) /  0.0000000000000000832498823872342992D+00 /
  data aip1cs( 25) / -0.0000000000000000260343254786947057D+00 /
  data aip1cs( 26) /  0.0000000000000000082159528142686287D+00 /
  data aip1cs( 27) / -0.0000000000000000026150406704984940D+00 /
  data aip1cs( 28) /  0.0000000000000000008390563463261051D+00 /
  data aip1cs( 29) / -0.0000000000000000002712685618629660D+00 /
  data aip1cs( 30) /  0.0000000000000000000883333375271942D+00 /
  data aip1cs( 31) / -0.0000000000000000000289603206822333D+00 /
  data aip1cs( 32) /  0.0000000000000000000095562185928676D+00 /
  data aip1cs( 33) / -0.0000000000000000000031727463569051D+00 /
  data aip1cs( 34) /  0.0000000000000000000010595576960768D+00 /
  data aip1cs( 35) / -0.0000000000000000000003558253765402D+00 /
  data aip1cs( 36) /  0.0000000000000000000001201334680517D+00 /
  data aip1cs( 37) / -0.0000000000000000000000407666883800D+00 /
  data aip1cs( 38) /  0.0000000000000000000000139016944446D+00 /
  data aip1cs( 39) / -0.0000000000000000000000047628165730D+00 /
  data aip1cs( 40) /  0.0000000000000000000000016391265551D+00 /
  data aip1cs( 41) / -0.0000000000000000000000005665491354D+00 /
  data aip1cs( 42) /  0.0000000000000000000000001966381969D+00 /
  data aip1cs( 43) / -0.0000000000000000000000000685230229D+00 /
  data aip1cs( 44) /  0.0000000000000000000000000239706939D+00 /
  data aip1cs( 45) / -0.0000000000000000000000000084166831D+00 /
  data aip1cs( 46) /  0.0000000000000000000000000029659364D+00 /
  data aip1cs( 47) / -0.0000000000000000000000000010487947D+00 /
  data aip1cs( 48) /  0.0000000000000000000000000003721150D+00 /
  data aip1cs( 49) / -0.0000000000000000000000000001324570D+00 /
  data aip1cs( 50) /  0.0000000000000000000000000000472976D+00 /
  data aip1cs( 51) / -0.0000000000000000000000000000169405D+00 /
  data aip1cs( 52) /  0.0000000000000000000000000000060855D+00 /
  data aip1cs( 53) / -0.0000000000000000000000000000021924D+00 /
  data aip1cs( 54) /  0.0000000000000000000000000000007920D+00 /
  data aip1cs( 55) / -0.0000000000000000000000000000002869D+00 /
  data aip1cs( 56) /  0.0000000000000000000000000000001042D+00 /
  data aip1cs( 57) / -0.0000000000000000000000000000000379D+00 /

  data naif / 0 /
  data naig / 0 /
  data naip1 / 0 /
  data naip2 / 0 /
  data x2sml / 0.0D+00 /
  data x32sml / 0.0D+00 /
  data x3sml / 0.0D+00 /
  data xbig / 0.0D+00 /

  if ( naif == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    naif = r8_inits ( aifcs, 13, eta )
    naig = r8_inits ( aigcs, 13, eta )
    naip1 = r8_inits ( aip1cs, 57, eta )
    naip2 = r8_inits ( aip2cs, 37, eta )
    x2sml = sqrt ( eta )
    x3sml = eta ** 0.3333D+00
    x32sml = 1.3104D+00 * x3sml * x3sml
    xbig = r8_mach ( 2 ) ** 0.6666D+00
  end if

  if ( x < - 1.0D+00 ) then
    call r8_admp ( x, xn, phi )
    r8_aide = xn * cos ( phi )
  else if ( abs ( x ) < x2sml ) then
    x2 = 0.0D+00
    x3 = 0.0D+00
    r8_aide = ( x2 * ( 0.125D+00 + r8_csevl ( x3, aifcs, naif ) ) &
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25D+00
    if ( x32sml < x ) then
      r8_aide = r8_aide * exp ( 2.0D+00 * x * sqrt ( x ) &
        / 3.0D+00 )
    end if
  else if ( abs ( x ) < x3sml ) then
    x2 = x * x
    x3 = 0.0D+00
    r8_aide = ( x2 * ( 0.125D+00 + r8_csevl ( x3, aifcs, naif ) ) &
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25D+00
    if ( x32sml < x ) then
      r8_aide = r8_aide * exp ( 2.0D+00 * x * sqrt ( x ) &
        / 3.0D+00 )
    end if
  else if ( x <= 1.0D+00 ) then
    x2 = x * x
    x3 = x * x
    r8_aide = ( x2 * ( 0.125D+00 + r8_csevl ( x3, aifcs, naif ) ) &
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25D+00
    if ( x32sml < x ) then
      r8_aide = r8_aide * exp ( 2.0D+00 * x * sqrt ( x ) &
        / 3.0D+00 )
    end if
  else if ( x <= 4.0D+00 ) then
    sqrtx = sqrt ( x )
    z = ( 16.0D+00  / ( x * sqrtx ) - 9.0D+00 ) / 7.0D+00
    r8_aide = ( - 0.28125D+00 - r8_csevl ( z, aip1cs, naip1 ) ) &
      * sqrt ( sqrtx )
  else if ( x < xbig ) then
    sqrtx = sqrt ( x )
    z = 16.0D+00  / ( x * sqrtx ) - 1.0D+00
    r8_aide = ( - 0.28125D+00 - r8_csevl ( z, aip2cs, naip2 ) ) &
      * sqrt ( sqrtx )
  else
    sqrtx = sqrt ( x )
    z = - 1.0D+00
    r8_aide = ( - 0.28125D+00 - r8_csevl ( z, aip2cs, naip2 ) ) &
      * sqrt ( sqrtx )
  end if

  return
end
function r8_aie ( x )

!*****************************************************************************80
!
!! R8_AIE evaluates the exponentially scaled Airy function Ai of an R8 argument.
!
!  Discussion:
!
!    if X <= 0,
!      R8_AIE ( X ) = R8_AI ( X )
!    else
!      R8_AIE ( X ) = R8_AI ( X ) * exp ( 2/3 * X^(3/2) )
!
!    Thanks to Aleksandra Piper for pointing out a correction involving a
!    missing assignment to SQRTX, 27 January 2012.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2012
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_AIE, the exponentially scaled Airy 
!    function Ai of X.
!
  implicit none

  real ( kind = 8 ) aifcs(13)
  real ( kind = 8 ) aigcs(13)
  real ( kind = 8 ) aip1cs(57)
  real ( kind = 8 ) aip2cs(37)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) naif
  integer ( kind = 4 ) naig
  integer ( kind = 4 ) naip1
  integer ( kind = 4 ) naip2
  real ( kind = 8 ) r8_aie
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sqrtx
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) x32sml
  real ( kind = 8 ) x3sml
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xm
  real ( kind = 8 ) z

  save aifcs
  save aigcs
  save aip1cs
  save aip2cs
  save naif
  save naig
  save naip1
  save naip2
  save x32sml
  save x3sml
  save xbig

  data aifcs(  1) / -0.37971358496669997496197089469414D-01 /
  data aifcs(  2) / +0.59191888537263638574319728013777D-01 /
  data aifcs(  3) / +0.98629280577279975365603891044060D-03 /
  data aifcs(  4) / +0.68488438190765667554854830182412D-05 /
  data aifcs(  5) / +0.25942025962194713019489279081403D-07 /
  data aifcs(  6) / +0.61766127740813750329445749697236D-10 /
  data aifcs(  7) / +0.10092454172466117901429556224601D-12 /
  data aifcs(  8) / +0.12014792511179938141288033225333D-15 /
  data aifcs(  9) / +0.10882945588716991878525295466666D-18 /
  data aifcs( 10) / +0.77513772196684887039238400000000D-22 /
  data aifcs( 11) / +0.44548112037175638391466666666666D-25 /
  data aifcs( 12) / +0.21092845231692343466666666666666D-28 /
  data aifcs( 13) / +0.83701735910741333333333333333333D-32 /

  data aigcs(  1) / +0.18152365581161273011556209957864D-01 /
  data aigcs(  2) / +0.21572563166010755534030638819968D-01 /
  data aigcs(  3) / +0.25678356987483249659052428090133D-03 /
  data aigcs(  4) / +0.14265214119792403898829496921721D-05 /
  data aigcs(  5) / +0.45721149200180426070434097558191D-08 /
  data aigcs(  6) / +0.95251708435647098607392278840592D-11 /
  data aigcs(  7) / +0.13925634605771399051150420686190D-13 /
  data aigcs(  8) / +0.15070999142762379592306991138666D-16 /
  data aigcs(  9) / +0.12559148312567778822703205333333D-19 /
  data aigcs( 10) / +0.83063073770821340343829333333333D-23 /
  data aigcs( 11) / +0.44657538493718567445333333333333D-26 /
  data aigcs( 12) / +0.19900855034518869333333333333333D-29 /
  data aigcs( 13) / +0.74702885256533333333333333333333D-33 /

  data aip1cs(  1) / -0.2146951858910538455460863467778D-01 /
  data aip1cs(  2) / -0.7535382535043301166219720865565D-02 /
  data aip1cs(  3) / +0.5971527949026380852035388881994D-03 /
  data aip1cs(  4) / -0.7283251254207610648502368291548D-04 /
  data aip1cs(  5) / +0.1110297130739299666517381821140D-04 /
  data aip1cs(  6) / -0.1950386152284405710346930314033D-05 /
  data aip1cs(  7) / +0.3786973885159515193885319670057D-06 /
  data aip1cs(  8) / -0.7929675297350978279039072879154D-07 /
  data aip1cs(  9) / +0.1762247638674256075568420122202D-07 /
  data aip1cs( 10) / -0.4110767539667195045029896593893D-08 /
  data aip1cs( 11) / +0.9984770057857892247183414107544D-09 /
  data aip1cs( 12) / -0.2510093251387122211349867730034D-09 /
  data aip1cs( 13) / +0.6500501929860695409272038601725D-10 /
  data aip1cs( 14) / -0.1727818405393616515478877107366D-10 /
  data aip1cs( 15) / +0.4699378842824512578362292872307D-11 /
  data aip1cs( 16) / -0.1304675656297743914491241246272D-11 /
  data aip1cs( 17) / +0.3689698478462678810473948382282D-12 /
  data aip1cs( 18) / -0.1061087206646806173650359679035D-12 /
  data aip1cs( 19) / +0.3098414384878187438660210070110D-13 /
  data aip1cs( 20) / -0.9174908079824139307833423547851D-14 /
  data aip1cs( 21) / +0.2752049140347210895693579062271D-14 /
  data aip1cs( 22) / -0.8353750115922046558091393301880D-15 /
  data aip1cs( 23) / +0.2563931129357934947568636168612D-15 /
  data aip1cs( 24) / -0.7950633762598854983273747289822D-16 /
  data aip1cs( 25) / +0.2489283634603069977437281175644D-16 /
  data aip1cs( 26) / -0.7864326933928735569664626221296D-17 /
  data aip1cs( 27) / +0.2505687311439975672324470645019D-17 /
  data aip1cs( 28) / -0.8047420364163909524537958682241D-18 /
  data aip1cs( 29) / +0.2604097118952053964443401104392D-18 /
  data aip1cs( 30) / -0.8486954164056412259482488834184D-19 /
  data aip1cs( 31) / +0.2784706882142337843359429186027D-19 /
  data aip1cs( 32) / -0.9195858953498612913687224151354D-20 /
  data aip1cs( 33) / +0.3055304318374238742247668225583D-20 /
  data aip1cs( 34) / -0.1021035455479477875902177048439D-20 /
  data aip1cs( 35) / +0.3431118190743757844000555680836D-21 /
  data aip1cs( 36) / -0.1159129341797749513376922463109D-21 /
  data aip1cs( 37) / +0.3935772844200255610836268229154D-22 /
  data aip1cs( 38) / -0.1342880980296717611956718989038D-22 /
  data aip1cs( 39) / +0.4603287883520002741659190305314D-23 /
  data aip1cs( 40) / -0.1585043927004064227810772499387D-23 /
  data aip1cs( 41) / +0.5481275667729675908925523755008D-24 /
  data aip1cs( 42) / -0.1903349371855047259064017948945D-24 /
  data aip1cs( 43) / +0.6635682302374008716777612115968D-25 /
  data aip1cs( 44) / -0.2322311650026314307975200986453D-25 /
  data aip1cs( 45) / +0.8157640113429179313142743695359D-26 /
  data aip1cs( 46) / -0.2875824240632900490057489929557D-26 /
  data aip1cs( 47) / +0.1017329450942901435079714319018D-26 /
  data aip1cs( 48) / -0.3610879108742216446575703490559D-27 /
  data aip1cs( 49) / +0.1285788540363993421256640342698D-27 /
  data aip1cs( 50) / -0.4592901037378547425160693022719D-28 /
  data aip1cs( 51) / +0.1645597033820713725812102485333D-28 /
  data aip1cs( 52) / -0.5913421299843501842087920271360D-29 /
  data aip1cs( 53) / +0.2131057006604993303479369509546D-29 /
  data aip1cs( 54) / -0.7701158157787598216982761745066D-30 /
  data aip1cs( 55) / +0.2790533307968930417581783777280D-30 /
  data aip1cs( 56) / -0.1013807715111284006452241367039D-30 /
  data aip1cs( 57) / +0.3692580158719624093658286216533D-31 /

  data aip2cs(  1) / -0.174314496929375513390355844011D-02 /
  data aip2cs(  2) / -0.167893854325541671632190613480D-02 /
  data aip2cs(  3) / +0.359653403352166035885983858114D-04 /
  data aip2cs(  4) / -0.138081860273922835457399383100D-05 /
  data aip2cs(  5) / +0.741122807731505298848699095233D-07 /
  data aip2cs(  6) / -0.500238203900133013130422866325D-08 /
  data aip2cs(  7) / +0.400693917417184240675446866355D-09 /
  data aip2cs(  8) / -0.367331242795905044199318496207D-10 /
  data aip2cs(  9) / +0.376034439592373852439592002918D-11 /
  data aip2cs( 10) / -0.422321332718747538026564938968D-12 /
  data aip2cs( 11) / +0.513509454033657070919618754120D-13 /
  data aip2cs( 12) / -0.669095850390477595651681356676D-14 /
  data aip2cs( 13) / +0.926667545641290648239550724382D-15 /
  data aip2cs( 14) / -0.135514382416070576333397356591D-15 /
  data aip2cs( 15) / +0.208115496312830995299006549335D-16 /
  data aip2cs( 16) / -0.334116499159176856871277570256D-17 /
  data aip2cs( 17) / +0.558578584585924316868032946585D-18 /
  data aip2cs( 18) / -0.969219040152365247518658209109D-19 /
  data aip2cs( 19) / +0.174045700128893206465696557738D-19 /
  data aip2cs( 20) / -0.322640979731130400247846333098D-20 /
  data aip2cs( 21) / +0.616074471106625258533259618986D-21 /
  data aip2cs( 22) / -0.120936347982490059076420676266D-21 /
  data aip2cs( 23) / +0.243632763310138108261570095786D-22 /
  data aip2cs( 24) / -0.502914221497457468943403144533D-23 /
  data aip2cs( 25) / +0.106224175543635689495470626133D-23 /
  data aip2cs( 26) / -0.229284284895989241509856324266D-24 /
  data aip2cs( 27) / +0.505181733929503744986884778666D-25 /
  data aip2cs( 28) / -0.113498123714412404979793920000D-25 /
  data aip2cs( 29) / +0.259765565985606980698374144000D-26 /
  data aip2cs( 30) / -0.605124621542939506172231679999D-27 /
  data aip2cs( 31) / +0.143359777966772800720295253333D-27 /
  data aip2cs( 32) / -0.345147757060899986280721066666D-28 /
  data aip2cs( 33) / +0.843875190213646740427025066666D-29 /
  data aip2cs( 34) / -0.209396142298188169434453333333D-29 /
  data aip2cs( 35) / +0.527008873478945503182848000000D-30 /
  data aip2cs( 36) / -0.134457433014553385789030399999D-30 /
  data aip2cs( 37) / +0.347570964526601147340117333333D-31 /

  data naif / 0 /
  data naig / 0 /
  data naip1 / 0 /
  data naip2 / 0 /
  data x32sml / 0.0D+00 /
  data x3sml / 0.0D+00 /
  data xbig / 0.0D+00 /

  if ( naif == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    naif = r8_inits ( aifcs, 13, eta )
    naig = r8_inits ( aigcs, 13, eta )
    naip1 = r8_inits ( aip1cs, 57, eta )
    naip2 = r8_inits ( aip2cs, 37, eta )
    x3sml = eta ** 0.3333D+00
    x32sml = 1.3104D+00 * x3sml * x3sml
    xbig = r8_mach ( 2 ) ** 0.6666D+00
  end if

  if ( x < - 1.0D+00 ) then
    call r8_aimp ( x, xm, theta )
    r8_aie = xm * cos ( theta )
  else if ( 0.0D+00 <= x .and. x <= x32sml ) then
    z = 0.0D+00
    r8_aie = 0.375D+000 + ( r8_csevl ( z, aifcs, naif ) &
      - x * ( 0.25D+00 + r8_csevl ( z, aigcs, naig ) ) )
  else if ( abs ( x ) <= x3sml ) then
    z = 0.0D+00
    r8_aie = 0.375D+000 + ( r8_csevl ( z, aifcs, naif ) &
      - x * ( 0.25D+00 + r8_csevl ( z, aigcs, naig ) ) )
    r8_aie = r8_aie * exp ( 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
  else if ( x <= 1.0D+00 ) then
    z = x * x * x
    r8_aie = 0.375D+000 + ( r8_csevl ( z, aifcs, naif ) &
      - x * ( 0.25D+00 + r8_csevl ( z, aigcs, naig ) ) )
    r8_aie = r8_aie * exp ( 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
  else if ( x <= 4.0D+00 ) then
    sqrtx = sqrt ( x )
    z = ( 16.0D+00 / ( x * sqrtx ) - 9.0D+00 ) / 7.0D+00
    r8_aie = ( 0.28125D+00 + r8_csevl ( z, aip1cs, naip1 ) ) &
      / sqrt ( sqrtx )
  else if ( x < xbig ) then
    sqrtx = sqrt ( x )
    z = 16.0D+00 / ( x * sqrtx ) - 1.0D+00
    r8_aie = ( 0.28125D+00 + r8_csevl ( z, aip2cs, naip2 ) ) &
      / sqrt ( sqrtx )
  else
    sqrtx = sqrt ( x )
    z = - 1.0D+00
    r8_aie = ( 0.28125D+00 + r8_csevl ( z, aip2cs, naip2 ) ) &
      / sqrt ( sqrtx )
  end if

  return
end
subroutine r8_aimp ( x, ampl, theta )

!*****************************************************************************80
!
!! R8_AIMP evaluates the modulus and phase of the Airy function.
!
!  Description:
!
!    This function must only be called when X <= -1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) AMPL, PHI, the modulus and phase of the 
!    Airy function.
!
  implicit none

  real ( kind = 8 ) am20cs(57)
  real ( kind = 8 ) am21cs(60)
  real ( kind = 8 ) am22cs(74)
  real ( kind = 8 ) ampl
  real ( kind = 8 ) ath0cs(53)
  real ( kind = 8 ) ath1cs(58)
  real ( kind = 8 ) ath2cs(72)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nam20
  integer ( kind = 4 ) nam21
  integer ( kind = 4 ) nam22
  integer ( kind = 4 ) nath0
  integer ( kind = 4 ) nath1
  integer ( kind = 4 ) nath2
  real ( kind = 8 ) pi4
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sqrtx
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml
  real ( kind = 8 ) z

  save am20cs
  save am21cs
  save am22cs
  save ath0cs
  save ath1cs
  save ath2cs
  save nam20
  save nam21
  save nam22
  save nath0
  save nath1
  save nath2
  save pi4
  save xsml

  data am20cs(  1) / +0.108716749086561856615730588125D-01 /
  data am20cs(  2) / +0.369489228982663555091728665146D-03 /
  data am20cs(  3) / +0.440680100484689563667507001327D-05 /
  data am20cs(  4) / +0.143686762361911153929183952833D-06 /
  data am20cs(  5) / +0.824275552390078308670628855353D-08 /
  data am20cs(  6) / +0.684426758893661606173927278180D-09 /
  data am20cs(  7) / +0.739566697282739287731004740213D-10 /
  data am20cs(  8) / +0.974595633696825017638702600847D-11 /
  data am20cs(  9) / +0.150076885829405775650973119497D-11 /
  data am20cs( 10) / +0.262147910221527634206252854802D-12 /
  data am20cs( 11) / +0.508354111376487180357278966914D-13 /
  data am20cs( 12) / +0.107684753358811440492985997070D-13 /
  data am20cs( 13) / +0.246091286618433429335914062617D-14 /
  data am20cs( 14) / +0.600786380358656418436110373550D-15 /
  data am20cs( 15) / +0.155449156102388071150651388384D-15 /
  data am20cs( 16) / +0.423535125035576604426382780182D-16 /
  data am20cs( 17) / +0.120862166289299840154401109189D-16 /
  data am20cs( 18) / +0.359609651214658240861499706423D-17 /
  data am20cs( 19) / +0.111134218386395638261774604677D-17 /
  data am20cs( 20) / +0.355559532432366609893680289225D-18 /
  data am20cs( 21) / +0.117433021600139309998766947387D-18 /
  data am20cs( 22) / +0.399397454661077561389162200966D-19 /
  data am20cs( 23) / +0.139576671528916310425606325640D-19 /
  data am20cs( 24) / +0.500240055309236041393459280716D-20 /
  data am20cs( 25) / +0.183552760958132679184834866457D-20 /
  data am20cs( 26) / +0.688490998179202743197790112404D-21 /
  data am20cs( 27) / +0.263631035611417012359996885105D-21 /
  data am20cs( 28) / +0.102924890237338360287153563785D-21 /
  data am20cs( 29) / +0.409246966671594885489762960571D-22 /
  data am20cs( 30) / +0.165558573406734651039727903828D-22 /
  data am20cs( 31) / +0.680797467063033356116599685727D-23 /
  data am20cs( 32) / +0.284326559934079832419751134476D-23 /
  data am20cs( 33) / +0.120507398348965255097287818819D-23 /
  data am20cs( 34) / +0.517961243287505217976613610424D-24 /
  data am20cs( 35) / +0.225622613427562816303268640887D-24 /
  data am20cs( 36) / +0.995418801147745168832117078246D-25 /
  data am20cs( 37) / +0.444551696397342424308280582053D-25 /
  data am20cs( 38) / +0.200865195461501101425916097338D-25 /
  data am20cs( 39) / +0.917786344151775165973885645402D-26 /
  data am20cs( 40) / +0.423872958105589240661672197948D-26 /
  data am20cs( 41) / +0.197789272007846092370846251490D-26 /
  data am20cs( 42) / +0.932116351284620665680435253373D-27 /
  data am20cs( 43) / +0.443482133249918099955611379722D-27 /
  data am20cs( 44) / +0.212945672365573895594589552837D-27 /
  data am20cs( 45) / +0.103158569651075977552209344907D-27 /
  data am20cs( 46) / +0.504023773022591199157904590029D-28 /
  data am20cs( 47) / +0.248301304570155945304046541005D-28 /
  data am20cs( 48) / +0.123301783128562196054198238560D-28 /
  data am20cs( 49) / +0.617033449920521746121976730507D-29 /
  data am20cs( 50) / +0.311092617415918897233869792213D-29 /
  data am20cs( 51) / +0.157983085201706173015269071503D-29 /
  data am20cs( 52) / +0.807931987538283607678121339092D-30 /
  data am20cs( 53) / +0.415997394138667562722951360052D-30 /
  data am20cs( 54) / +0.215610934097716900471935862504D-30 /
  data am20cs( 55) / +0.112468857265869178296752823613D-30 /
  data am20cs( 56) / +0.590331560632838091123040811797D-31 /
  data am20cs( 57) / +0.311735667692928562046280505333D-31 /

  data ath0cs(  1) / -0.8172601764161634499840208700543D-01 /
  data ath0cs(  2) / -0.8004012824788273287596481113068D-03 /
  data ath0cs(  3) / -0.3186525268782113203795553628242D-05 /
  data ath0cs(  4) / -0.6688388266477509330741698865033D-07 /
  data ath0cs(  5) / -0.2931759284994564516506822463184D-08 /
  data ath0cs(  6) / -0.2011263760883621669049030307186D-09 /
  data ath0cs(  7) / -0.1877522678055973426074008166652D-10 /
  data ath0cs(  8) / -0.2199637137704601251899002199848D-11 /
  data ath0cs(  9) / -0.3071616682592272449025746605586D-12 /
  data ath0cs( 10) / -0.4936140553673418361025600985389D-13 /
  data ath0cs( 11) / -0.8902833722583660416935236969866D-14 /
  data ath0cs( 12) / -0.1768987764615272613656814199467D-14 /
  data ath0cs( 13) / -0.3817868689032277014678199609600D-15 /
  data ath0cs( 14) / -0.8851159014819947594156286509984D-16 /
  data ath0cs( 15) / -0.2184818181414365953149677679568D-16 /
  data ath0cs( 16) / -0.5700849046986452380599442295119D-17 /
  data ath0cs( 17) / -0.1563121122177875392516031795495D-17 /
  data ath0cs( 18) / -0.4481437996768995067906688776353D-18 /
  data ath0cs( 19) / -0.1337794883736188022044566044098D-18 /
  data ath0cs( 20) / -0.4143340036874114453776852445442D-19 /
  data ath0cs( 21) / -0.1327263385718805025080481164652D-19 /
  data ath0cs( 22) / -0.4385728589128440522215756835955D-20 /
  data ath0cs( 23) / -0.1491360695952818067686201743956D-20 /
  data ath0cs( 24) / -0.5208104738630711377154238188773D-21 /
  data ath0cs( 25) / -0.1864382222390498923872526604979D-21 /
  data ath0cs( 26) / -0.6830263751167969012975435381881D-22 /
  data ath0cs( 27) / -0.2557117058029329629296207591347D-22 /
  data ath0cs( 28) / -0.9770158640254300218246907254046D-23 /
  data ath0cs( 29) / -0.3805161433416679084068428254886D-23 /
  data ath0cs( 30) / -0.1509022750737054063493926482995D-23 /
  data ath0cs( 31) / -0.6087551341242424929005568014525D-24 /
  data ath0cs( 32) / -0.2495879513809711495425982124058D-24 /
  data ath0cs( 33) / -0.1039157654581920948909588084274D-24 /
  data ath0cs( 34) / -0.4390235913976846536974594969051D-25 /
  data ath0cs( 35) / -0.1880790678447990211675826820582D-25 /
  data ath0cs( 36) / -0.8165070764199462948863022205753D-26 /
  data ath0cs( 37) / -0.3589944503749750514266435585041D-26 /
  data ath0cs( 38) / -0.1597658126632132872981291608708D-26 /
  data ath0cs( 39) / -0.7193250175703823969113802835305D-27 /
  data ath0cs( 40) / -0.3274943012727856506209351132721D-27 /
  data ath0cs( 41) / -0.1507042445783690665816975047272D-27 /
  data ath0cs( 42) / -0.7006624198319904717843967949140D-28 /
  data ath0cs( 43) / -0.3289907402983718226528815678356D-28 /
  data ath0cs( 44) / -0.1559518084365146526445322711496D-28 /
  data ath0cs( 45) / -0.7460690508208254582833851119721D-29 /
  data ath0cs( 46) / -0.3600877034824662020563277249431D-29 /
  data ath0cs( 47) / -0.1752851437473772257350402219197D-29 /
  data ath0cs( 48) / -0.8603275775188512909623778628724D-30 /
  data ath0cs( 49) / -0.4256432603226946534668039480105D-30 /
  data ath0cs( 50) / -0.2122161865044262927723650698206D-30 /
  data ath0cs( 51) / -0.1065996156704879052472060798561D-30 /
  data ath0cs( 52) / -0.5393568608816949116410688086892D-31 /
  data ath0cs( 53) / -0.2748174851043954822278496517870D-31 /

  data am21cs(  1) / +0.592790266721309588375717482814D-02 /
  data am21cs(  2) / +0.200569405393165186428695217690D-02 /
  data am21cs(  3) / +0.911081850262275893553072526291D-04 /
  data am21cs(  4) / +0.849894306372047155633172107475D-05 /
  data am21cs(  5) / +0.113297908976913076637929215494D-05 /
  data am21cs(  6) / +0.187517946100666496180950627804D-06 /
  data am21cs(  7) / +0.359306519018245832699035211192D-07 /
  data am21cs(  8) / +0.765757714071683864039093517470D-08 /
  data am21cs(  9) / +0.176999967168039173925953460744D-08 /
  data am21cs( 10) / +0.436259555654598932720546585535D-09 /
  data am21cs( 11) / +0.113291641337853230035520085219D-09 /
  data am21cs( 12) / +0.307257690982419244137868398126D-10 /
  data am21cs( 13) / +0.864482416482201075541200465766D-11 /
  data am21cs( 14) / +0.251015250060924402115104562212D-11 /
  data am21cs( 15) / +0.749102496764440371601802227751D-12 /
  data am21cs( 16) / +0.228996928487994073089565214432D-12 /
  data am21cs( 17) / +0.715113658927987694949327491175D-13 /
  data am21cs( 18) / +0.227607924959566841946395165061D-13 /
  data am21cs( 19) / +0.736942142760886513969953227782D-14 /
  data am21cs( 20) / +0.242328675267827490463991742006D-14 /
  data am21cs( 21) / +0.808153774548239869283406558403D-15 /
  data am21cs( 22) / +0.273008079804356086659174563386D-15 /
  data am21cs( 23) / +0.933236070891385318473519474326D-16 /
  data am21cs( 24) / +0.322508099681084622213867546973D-16 /
  data am21cs( 25) / +0.112581932346444541217757573416D-16 /
  data am21cs( 26) / +0.396699463986938821660259459530D-17 /
  data am21cs( 27) / +0.141006567944319504660865034527D-17 /
  data am21cs( 28) / +0.505302086537851213375537393032D-18 /
  data am21cs( 29) / +0.182461523215945141197999102789D-18 /
  data am21cs( 30) / +0.663584568262130466928029121642D-19 /
  data am21cs( 31) / +0.242963731631276179741747455826D-19 /
  data am21cs( 32) / +0.895238915123687802013669922963D-20 /
  data am21cs( 33) / +0.331845289350050791260229250755D-20 /
  data am21cs( 34) / +0.123706196188658315384437905922D-20 /
  data am21cs( 35) / +0.463636677012390840306767734243D-21 /
  data am21cs( 36) / +0.174653135947764475469758765989D-21 /
  data am21cs( 37) / +0.661116810234991176307910643111D-22 /
  data am21cs( 38) / +0.251409918994072486176125666459D-22 /
  data am21cs( 39) / +0.960274995571732568694034386998D-23 /
  data am21cs( 40) / +0.368324952289296395686436898078D-23 /
  data am21cs( 41) / +0.141843138269159136145535939553D-23 /
  data am21cs( 42) / +0.548342674276935830106345800990D-24 /
  data am21cs( 43) / +0.212761054623118806650372562616D-24 /
  data am21cs( 44) / +0.828443700849418591487734760953D-25 /
  data am21cs( 45) / +0.323670563926127001421028600927D-25 /
  data am21cs( 46) / +0.126868882963286057355055062493D-25 /
  data am21cs( 47) / +0.498843818992121626935068934362D-26 /
  data am21cs( 48) / +0.196734584467649390967119381790D-26 /
  data am21cs( 49) / +0.778135971020326957713212064836D-27 /
  data am21cs( 50) / +0.308633941498911152919192968451D-27 /
  data am21cs( 51) / +0.122744647045453119789338037234D-27 /
  data am21cs( 52) / +0.489431279134292205885241216204D-28 /
  data am21cs( 53) / +0.195646879802909821175925099724D-28 /
  data am21cs( 54) / +0.783988952922426171166311492266D-29 /
  data am21cs( 55) / +0.314896914002484223748298978099D-29 /
  data am21cs( 56) / +0.126769763137250681307067842559D-29 /
  data am21cs( 57) / +0.511470691906900141641632107724D-30 /
  data am21cs( 58) / +0.206801709795538770250900316706D-30 /
  data am21cs( 59) / +0.837891344768519001325996867583D-31 /
  data am21cs( 60) / +0.340168991971489802052339079577D-31 /

  data ath1cs(  1) / -0.6972849916208883845888148415037D-01 /
  data ath1cs(  2) / -0.5108722790650044987073448077961D-02 /
  data ath1cs(  3) / -0.8644335996989755094525334749512D-04 /
  data ath1cs(  4) / -0.5604720044235263542188698916125D-05 /
  data ath1cs(  5) / -0.6045735125623897409156376640077D-06 /
  data ath1cs(  6) / -0.8639802632488334393219721138499D-07 /
  data ath1cs(  7) / -0.1480809484309927157147782480780D-07 /
  data ath1cs(  8) / -0.2885809334577236039999449908712D-08 /
  data ath1cs(  9) / -0.6191631975665699609309191231800D-09 /
  data ath1cs( 10) / -0.1431992808860957830931365259879D-09 /
  data ath1cs( 11) / -0.3518141102137214721504616874321D-10 /
  data ath1cs( 12) / -0.9084761919955078290070339808051D-11 /
  data ath1cs( 13) / -0.2446171672688598449343283664767D-11 /
  data ath1cs( 14) / -0.6826083203213446240828996710264D-12 /
  data ath1cs( 15) / -0.1964579931194940171278546257802D-12 /
  data ath1cs( 16) / -0.5808933227139693164009191265856D-13 /
  data ath1cs( 17) / -0.1759042249527441992795400959024D-13 /
  data ath1cs( 18) / -0.5440902932714896613632538945319D-14 /
  data ath1cs( 19) / -0.1715247407486806802622358519451D-14 /
  data ath1cs( 20) / -0.5500929233576991546871101847161D-15 /
  data ath1cs( 21) / -0.1791878287739317259495152638754D-15 /
  data ath1cs( 22) / -0.5920372520086694197778411062231D-16 /
  data ath1cs( 23) / -0.1981713027876483962470972206590D-16 /
  data ath1cs( 24) / -0.6713232347016352262049984343790D-17 /
  data ath1cs( 25) / -0.2299450243658281116122358619832D-17 /
  data ath1cs( 26) / -0.7957300928236376595304637145634D-18 /
  data ath1cs( 27) / -0.2779994027291784157172290233739D-18 /
  data ath1cs( 28) / -0.9798924361326985224406795480814D-19 /
  data ath1cs( 29) / -0.3482717006061574386702645565849D-19 /
  data ath1cs( 30) / -0.1247489122558599057173300058084D-19 /
  data ath1cs( 31) / -0.4501210041478228113487751824452D-20 /
  data ath1cs( 32) / -0.1635346244013352135596114164667D-20 /
  data ath1cs( 33) / -0.5980102897780336268098762265941D-21 /
  data ath1cs( 34) / -0.2200246286286123454028196295475D-21 /
  data ath1cs( 35) / -0.8142463073515085897408205291519D-22 /
  data ath1cs( 36) / -0.3029924773660042537432330709674D-22 /
  data ath1cs( 37) / -0.1133390098574623537722943969689D-22 /
  data ath1cs( 38) / -0.4260766024749295719283049889791D-23 /
  data ath1cs( 39) / -0.1609363396278189718797500634453D-23 /
  data ath1cs( 40) / -0.6106377190825026293045330444287D-24 /
  data ath1cs( 41) / -0.2326954318021694061836577887573D-24 /
  data ath1cs( 42) / -0.8903987877472252604474129558186D-25 /
  data ath1cs( 43) / -0.3420558530005675024117914752341D-25 /
  data ath1cs( 44) / -0.1319026715257272659017212100607D-25 /
  data ath1cs( 45) / -0.5104899493612043091316191177386D-26 /
  data ath1cs( 46) / -0.1982599478474547451242444663466D-26 /
  data ath1cs( 47) / -0.7725702356880830535636111851519D-27 /
  data ath1cs( 48) / -0.3020234733664680100815776863573D-27 /
  data ath1cs( 49) / -0.1184379739074169993712946380800D-27 /
  data ath1cs( 50) / -0.4658430227922308520573252840106D-28 /
  data ath1cs( 51) / -0.1837554188100384647157502006613D-28 /
  data ath1cs( 52) / -0.7268566894427990953321876684800D-29 /
  data ath1cs( 53) / -0.2882863120391468135527089875626D-29 /
  data ath1cs( 54) / -0.1146374629459906350417591664639D-29 /
  data ath1cs( 55) / -0.4570031437748533058179991688533D-30 /
  data ath1cs( 56) / -0.1826276602045346104809934028799D-30 /
  data ath1cs( 57) / -0.7315349993385250469111066350933D-31 /
  data ath1cs( 58) / -0.2936925599971429781637815773866D-31 /

  data am22cs(  1) / -0.156284448062534112753545828583D-01 /
  data am22cs(  2) / +0.778336445239681307018943100334D-02 /
  data am22cs(  3) / +0.867057770477189528406072812110D-03 /
  data am22cs(  4) / +0.156966273156113719469953482266D-03 /
  data am22cs(  5) / +0.356396257143286511324100666302D-04 /
  data am22cs(  6) / +0.924598335425043154495080090994D-05 /
  data am22cs(  7) / +0.262110161850422389523194982066D-05 /
  data am22cs(  8) / +0.791882216516012561489469982263D-06 /
  data am22cs(  9) / +0.251041527921011847803162690862D-06 /
  data am22cs( 10) / +0.826522320665407734472997712940D-07 /
  data am22cs( 11) / +0.280571166281305264396384290014D-07 /
  data am22cs( 12) / +0.976821090484680786674631273890D-08 /
  data am22cs( 13) / +0.347407923227710343287279035573D-08 /
  data am22cs( 14) / +0.125828132169836914219092738164D-08 /
  data am22cs( 15) / +0.462988260641895264497330784625D-09 /
  data am22cs( 16) / +0.172728258813604072468143128696D-09 /
  data am22cs( 17) / +0.652319200131154135148574124970D-10 /
  data am22cs( 18) / +0.249047168520982056019881087112D-10 /
  data am22cs( 19) / +0.960156820553765948078189890126D-11 /
  data am22cs( 20) / +0.373448002067726856974776596757D-11 /
  data am22cs( 21) / +0.146417565032053391722216189678D-11 /
  data am22cs( 22) / +0.578265471168512825475827881553D-12 /
  data am22cs( 23) / +0.229915407244706118560254184494D-12 /
  data am22cs( 24) / +0.919780711231997257150883662365D-13 /
  data am22cs( 25) / +0.370060068813090065807504045556D-13 /
  data am22cs( 26) / +0.149675761698672987823326345205D-13 /
  data am22cs( 27) / +0.608361194938461148720451399443D-14 /
  data am22cs( 28) / +0.248404087115121397635425326873D-14 /
  data am22cs( 29) / +0.101862476526769080727914465339D-14 /
  data am22cs( 30) / +0.419383856352753989429640310957D-15 /
  data am22cs( 31) / +0.173318901762930756149702493501D-15 /
  data am22cs( 32) / +0.718821902388508517820445406811D-16 /
  data am22cs( 33) / +0.299123633598403607712470896113D-16 /
  data am22cs( 34) / +0.124868990433238627855713110880D-16 /
  data am22cs( 35) / +0.522829344609483661928651193632D-17 /
  data am22cs( 36) / +0.219532961724713396595998454359D-17 /
  data am22cs( 37) / +0.924298325229777281154410024332D-18 /
  data am22cs( 38) / +0.390157708236091407825543197309D-18 /
  data am22cs( 39) / +0.165093892693863707213759030367D-18 /
  data am22cs( 40) / +0.700221815715994367565716554487D-19 /
  data am22cs( 41) / +0.297651833616786915573214963506D-19 /
  data am22cs( 42) / +0.126796539086902072571134261229D-19 /
  data am22cs( 43) / +0.541243400697077628687581725061D-20 /
  data am22cs( 44) / +0.231487350218155252296382133283D-20 /
  data am22cs( 45) / +0.991920288386566563462623851167D-21 /
  data am22cs( 46) / +0.425803015323732357158897608174D-21 /
  data am22cs( 47) / +0.183101842973024501678402003088D-21 /
  data am22cs( 48) / +0.788678712311075375564526811022D-22 /
  data am22cs( 49) / +0.340254607386229874956582997235D-22 /
  data am22cs( 50) / +0.147020881405712530791860892535D-22 /
  data am22cs( 51) / +0.636211018324916957733348071767D-23 /
  data am22cs( 52) / +0.275707050680980721919395987768D-23 /
  data am22cs( 53) / +0.119645858090104071356261780457D-23 /
  data am22cs( 54) / +0.519912545729242147981768210567D-24 /
  data am22cs( 55) / +0.226217674847104475260575286850D-24 /
  data am22cs( 56) / +0.985526113754431819448565068283D-25 /
  data am22cs( 57) / +0.429870630332508717223681286187D-25 /
  data am22cs( 58) / +0.187723641661580639829657670189D-25 /
  data am22cs( 59) / +0.820721941772842137268801052115D-26 /
  data am22cs( 60) / +0.359214665604615507812767944463D-26 /
  data am22cs( 61) / +0.157390594612773315611458940587D-26 /
  data am22cs( 62) / +0.690329781039333834965319153586D-27 /
  data am22cs( 63) / +0.303092079078968534607859331415D-27 /
  data am22cs( 64) / +0.133204934160481219185689121944D-27 /
  data am22cs( 65) / +0.585978836851523490117937981442D-28 /
  data am22cs( 66) / +0.258016868489487806338425080457D-28 /
  data am22cs( 67) / +0.113712433637283667223632182863D-28 /
  data am22cs( 68) / +0.501592557226068509236430548549D-29 /
  data am22cs( 69) / +0.221445829395509373322569708484D-29 /
  data am22cs( 70) / +0.978470283886507289984691416411D-30 /
  data am22cs( 71) / +0.432695414934180170112000952983D-30 /
  data am22cs( 72) / +0.191497288193994570612929860440D-30 /
  data am22cs( 73) / +0.848164622402392354171298331562D-31 /
  data am22cs( 74) / +0.375947065173955919947455052934D-31 /

  data ath2cs(  1) / +0.4405273458718778997061127057775D-02 /
  data ath2cs(  2) / -0.3042919452318454608483844239873D-01 /
  data ath2cs(  3) / -0.1385653283771793791602692842653D-02 /
  data ath2cs(  4) / -0.1804443908954952302670486910952D-03 /
  data ath2cs(  5) / -0.3380847108327308671057465323618D-04 /
  data ath2cs(  6) / -0.7678183535229023055257676817765D-05 /
  data ath2cs(  7) / -0.1967839443716035324690935417077D-05 /
  data ath2cs(  8) / -0.5483727115877700361586143659281D-06 /
  data ath2cs(  9) / -0.1625461550532612452712696212258D-06 /
  data ath2cs( 10) / -0.5053049981268895015277637842078D-07 /
  data ath2cs( 11) / -0.1631580701124066881183851715617D-07 /
  data ath2cs( 12) / -0.5434204112348517507963436694817D-08 /
  data ath2cs( 13) / -0.1857398556409900325763850109630D-08 /
  data ath2cs( 14) / -0.6489512033326108816213513640676D-09 /
  data ath2cs( 15) / -0.2310594885800944720482995987079D-09 /
  data ath2cs( 16) / -0.8363282183204411682819329546745D-10 /
  data ath2cs( 17) / -0.3071196844890191462660661303891D-10 /
  data ath2cs( 18) / -0.1142367142432716819409514579892D-10 /
  data ath2cs( 19) / -0.4298116066345803065822470108971D-11 /
  data ath2cs( 20) / -0.1633898699596715440601646086632D-11 /
  data ath2cs( 21) / -0.6269328620016619432123443754076D-12 /
  data ath2cs( 22) / -0.2426052694816257357356159203991D-12 /
  data ath2cs( 23) / -0.9461198321624039090742527765052D-13 /
  data ath2cs( 24) / -0.3716060313411504806847798281269D-13 /
  data ath2cs( 25) / -0.1469155684097526763170138810309D-13 /
  data ath2cs( 26) / -0.5843694726140911944556401363094D-14 /
  data ath2cs( 27) / -0.2337502595591951298832675034934D-14 /
  data ath2cs( 28) / -0.9399231371171435401160167358411D-15 /
  data ath2cs( 29) / -0.3798014669372894500076335263715D-15 /
  data ath2cs( 30) / -0.1541731043984972524883443681775D-15 /
  data ath2cs( 31) / -0.6285287079535307162925662365202D-16 /
  data ath2cs( 32) / -0.2572731812811455424755383992774D-16 /
  data ath2cs( 33) / -0.1057098119354017809340974866555D-16 /
  data ath2cs( 34) / -0.4359080267402696966695992699964D-17 /
  data ath2cs( 35) / -0.1803634315959978013953176945540D-17 /
  data ath2cs( 36) / -0.7486838064380536821719431676914D-18 /
  data ath2cs( 37) / -0.3117261367347604656799597209985D-18 /
  data ath2cs( 38) / -0.1301687980927700734792871620696D-18 /
  data ath2cs( 39) / -0.5450527587519522468973883909909D-19 /
  data ath2cs( 40) / -0.2288293490114231872268635931903D-19 /
  data ath2cs( 41) / -0.9631059503829538655655060440088D-20 /
  data ath2cs( 42) / -0.4063281001524614089092195416434D-20 /
  data ath2cs( 43) / -0.1718203980908026763900413858510D-20 /
  data ath2cs( 44) / -0.7281574619892536367415322473328D-21 /
  data ath2cs( 45) / -0.3092352652680643127960680345790D-21 /
  data ath2cs( 46) / -0.1315917855965440490383417023254D-21 /
  data ath2cs( 47) / -0.5610606786087055512664907412668D-22 /
  data ath2cs( 48) / -0.2396621894086355206020304337895D-22 /
  data ath2cs( 49) / -0.1025574332390581200832954423924D-22 /
  data ath2cs( 50) / -0.4396264138143656476403607323663D-23 /
  data ath2cs( 51) / -0.1887652998372577373342508719450D-23 /
  data ath2cs( 52) / -0.8118140359576807603579433230445D-24 /
  data ath2cs( 53) / -0.3496734274366286856375952089214D-24 /
  data ath2cs( 54) / -0.1508402925156873215171751475867D-24 /
  data ath2cs( 55) / -0.6516268284778671059787773834341D-25 /
  data ath2cs( 56) / -0.2818945797529207424505942114583D-25 /
  data ath2cs( 57) / -0.1221127596512262744598094464505D-25 /
  data ath2cs( 58) / -0.5296674341169867168620011705073D-26 /
  data ath2cs( 59) / -0.2300359270773673431358870971744D-26 /
  data ath2cs( 60) / -0.1000279482355367494781220348930D-26 /
  data ath2cs( 61) / -0.4354760404180879394806893162179D-27 /
  data ath2cs( 62) / -0.1898056134741477522515482827030D-27 /
  data ath2cs( 63) / -0.8282111868712974697554009309315D-28 /
  data ath2cs( 64) / -0.3617815493066569006586213484374D-28 /
  data ath2cs( 65) / -0.1582018896178003654858941843636D-28 /
  data ath2cs( 66) / -0.6925068597802270011772820383247D-29 /
  data ath2cs( 67) / -0.3034390239778629128908629727335D-29 /
  data ath2cs( 68) / -0.1330889568166725224761977446509D-29 /
  data ath2cs( 69) / -0.5842848522173090120487606971706D-30 /
  data ath2cs( 70) / -0.2567488423238302631121274357678D-30 /
  data ath2cs( 71) / -0.1129232322268882185791505819151D-30 /
  data ath2cs( 72) / -0.4970947029753336916550570105023D-31 /

  data nam20 / 0 /
  data nam21 / 0 /
  data nam22 / 0 /
  data nath0 / 0 /
  data nath1 / 0 /
  data nath2 / 0 /
  data pi4 / 0.78539816339744830961566084581988D+00 /
  data xsml / 0.0D+00 /

  if ( nam20 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nam20 = r8_inits ( am20cs, 57, eta )
    nath0 = r8_inits ( ath0cs, 53, eta )
    nam21 = r8_inits ( am21cs, 60, eta )
    nath1 = r8_inits ( ath1cs, 58, eta )
    nam22 = r8_inits ( am22cs, 74, eta )
    nath2 = r8_inits ( ath2cs, 72, eta )
    xsml = - ( 128.0D+00 / r8_mach ( 3 ) ) ** 0.3333D+00
  end if

  if ( x <= xsml ) then
    z = 1.0D+00
    ampl = 0.3125D+00 + r8_csevl ( z, am20cs, nam20 )
    theta = - 0.625D+00 + r8_csevl ( z, ath0cs, nath0 )
  else if ( x < - 4.0D+00 ) then
    z = 128.0D+00 / x / x / x + 1.0D+00
    ampl = 0.3125D+00 + r8_csevl ( z, am20cs, nam20 )
    theta = - 0.625D+00 + r8_csevl ( z, ath0cs, nath0 )
  else if ( x < - 2.0D+00 ) then
    z = ( 128.0D+00 / x / x / x + 9.0D+00 ) / 7.0D+00
    ampl = 0.3125D+00 + r8_csevl ( z, am21cs, nam21 )
    theta = - 0.625D+00 + r8_csevl ( z, ath1cs, nath1 )
  else if ( x <= - 1.0D+00 ) then
    z = ( 16.0D+00 / x / x / x + 9.0D+00 ) / 7.0D+00
    ampl = 0.3125D+00 + r8_csevl ( z, am22cs, nam22 )
    theta = - 0.625D+00 + r8_csevl ( z, ath2cs, nath2 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_AIMP - Fatal error!'
    write ( *, '(a)' ) '  -1.0 < X.'
    stop
  end if

  sqrtx = sqrt ( - x )
  ampl = sqrt ( ampl / sqrtx )
  theta = pi4 - x * sqrtx * theta

  return
end
function r8_aint ( x )

!*****************************************************************************80
!
!! R8_AINT truncates an R8 argument to an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    1 September 2011
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_AINT, the truncated version of X.
!
  implicit none

  real ( kind = 8 ) r8_aint
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( x < 0.0E+00 ) then
    value = - real ( int ( abs ( x ) ), kind = 8 )
  else
    value =   real ( int ( abs ( x ) ), kind = 8 )
  end if

  r8_aint = value

  return
end
function r8_asin ( x )

!*****************************************************************************80
!
!! R8_ASIN evaluates the arc-sine of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_ASIN, the arc-sine of X.
!
  implicit none

  real ( kind = 8 ) asincs(39)
  integer ( kind = 4 ) nterms
  real ( kind = 8 ) pi2
  real ( kind = 8 ) r8_asin
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  save asincs
  save nterms
  save pi2

  data asincs(  1) / +0.10246391753227159336573148305785D+00 /
  data asincs(  2) / +0.54946487221245833306011195902924D-01 /
  data asincs(  3) / +0.40806303925449692851307056149246D-02 /
  data asincs(  4) / +0.40789006854604435455598823905612D-03 /
  data asincs(  5) / +0.46985367432203691616048530136218D-04 /
  data asincs(  6) / +0.58809758139708058986454385552074D-05 /
  data asincs(  7) / +0.77732312462777632750557528163795D-06 /
  data asincs(  8) / +0.10677423340082039235047504956587D-06 /
  data asincs(  9) / +0.15092399536022808262386434401064D-07 /
  data asincs( 10) / +0.21809724080055385496609614713930D-08 /
  data asincs( 11) / +0.32075984262789614433261959667376D-09 /
  data asincs( 12) / +0.47855369646781034461493133918953D-10 /
  data asincs( 13) / +0.72251287362910432263848754537112D-11 /
  data asincs( 14) / +0.11018334742255783705372701334987D-11 /
  data asincs( 15) / +0.16947632539203354877423745651078D-12 /
  data asincs( 16) / +0.26261558667348224162283241502416D-13 /
  data asincs( 17) / +0.40958299813281178408828069291110D-14 /
  data asincs( 18) / +0.64244793108803655891727944887091D-15 /
  data asincs( 19) / +0.10128142198228221693973361222041D-15 /
  data asincs( 20) / +0.16039221897380787560050597464746D-16 /
  data asincs( 21) / +0.25503501355807141715298789676373D-17 /
  data asincs( 22) / +0.40701403797862382855487165672106D-18 /
  data asincs( 23) / +0.65172671712881144437889267575466D-19 /
  data asincs( 24) / +0.10467453037096796954244891716266D-19 /
  data asincs( 25) / +0.16858725563380328094989095185066D-20 /
  data asincs( 26) / +0.27221936305040227625164341247999D-21 /
  data asincs( 27) / +0.44059293900347550617126830079999D-22 /
  data asincs( 28) / +0.71466685243375937853063168000000D-23 /
  data asincs( 29) / +0.11615793343859516051798971733333D-23 /
  data asincs( 30) / +0.18915234552354685801184187733333D-24 /
  data asincs( 31) / +0.30855772044244342399827968000000D-25 /
  data asincs( 32) / +0.50416366022162453412970495999999D-26 /
  data asincs( 33) / +0.82502725502400865081753600000000D-27 /
  data asincs( 34) / +0.13520032631020947208055466666666D-27 /
  data asincs( 35) / +0.22184326876541720216644266666666D-28 /
  data asincs( 36) / +0.36442494054085079212578133333333D-29 /
  data asincs( 37) / +0.59920218558643813307733333333333D-30 /
  data asincs( 38) / +0.98584812059573785810261333333333D-31 /
  data asincs( 39) / +0.16222501166399014393173333333333D-31 /

  data nterms / 0 /
  data pi2 / 1.57079632679489661923132169163975D+00 /

  if ( nterms == 0 ) then
    nterms = r8_inits ( asincs, 39, 0.1D+00 * r8_mach ( 3 ) )
    sqeps = sqrt ( 6.0D+00 * r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( x < - 1.0D+00 - sqeps ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ASIN - Fatal error!'
    write ( *, '(a)' ) '  X < - 1.0'
    stop

  else if ( x < - 1.0D+00 ) then

    value = - pi2

  else if ( x < 1.0D+00 ) then

    z = 0.0D+00
    if ( sqeps < y ) then
      z = y * y
    end if

    if ( z <= 0.5D+00 ) then
      value = x * ( 1.0D+00 + r8_csevl ( 4.0D+00 * z - 1.0D+00, &
        asincs, nterms ) )
    else
      value = pi2 - sqrt ( 1.0D+00 - z ) * ( 1.0D+00 + &
        r8_csevl ( 3.0D+00 - 4.0D+00 * z, asincs, nterms ) )
    end if

    if ( x < 0.0D+00 ) then
      value = - abs ( value )
    else if ( 0.0D+00 < x ) then
      value = + abs ( value )
    end if

  else if ( x < 1.0D+00 + sqeps ) then

    value = pi2

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ASIN - Fatal error!'
    write ( *, '(a)' ) '  1.0 < X'
    stop

  end if

  r8_asin = value

  return
end
function r8_asinh ( x )

!*****************************************************************************80
!
!! R8_ASINH evaluates the arc-sine of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_ASINH, the arc-hyperbolic sine of X.
!
  implicit none

  real ( kind = 8 ) aln2
  real ( kind = 8 ) asnhcs(39)
  integer ( kind = 4 ) nterms
  real ( kind = 8 ) r8_asinh
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) y

  save aln2
  save asnhcs
  save nterms
  save xmax

  data asnhcs(  1) / -0.12820039911738186343372127359268D+00 /
  data asnhcs(  2) / -0.58811761189951767565211757138362D-01 /
  data asnhcs(  3) / +0.47274654322124815640725249756029D-02 /
  data asnhcs(  4) / -0.49383631626536172101360174790273D-03 /
  data asnhcs(  5) / +0.58506207058557412287494835259321D-04 /
  data asnhcs(  6) / -0.74669983289313681354755069217188D-05 /
  data asnhcs(  7) / +0.10011693583558199265966192015812D-05 /
  data asnhcs(  8) / -0.13903543858708333608616472258886D-06 /
  data asnhcs(  9) / +0.19823169483172793547317360237148D-07 /
  data asnhcs( 10) / -0.28847468417848843612747272800317D-08 /
  data asnhcs( 11) / +0.42672965467159937953457514995907D-09 /
  data asnhcs( 12) / -0.63976084654366357868752632309681D-10 /
  data asnhcs( 13) / +0.96991686089064704147878293131179D-11 /
  data asnhcs( 14) / -0.14844276972043770830246658365696D-11 /
  data asnhcs( 15) / +0.22903737939027447988040184378983D-12 /
  data asnhcs( 16) / -0.35588395132732645159978942651310D-13 /
  data asnhcs( 17) / +0.55639694080056789953374539088554D-14 /
  data asnhcs( 18) / -0.87462509599624678045666593520162D-15 /
  data asnhcs( 19) / +0.13815248844526692155868802298129D-15 /
  data asnhcs( 20) / -0.21916688282900363984955142264149D-16 /
  data asnhcs( 21) / +0.34904658524827565638313923706880D-17 /
  data asnhcs( 22) / -0.55785788400895742439630157032106D-18 /
  data asnhcs( 23) / +0.89445146617134012551050882798933D-19 /
  data asnhcs( 24) / -0.14383426346571317305551845239466D-19 /
  data asnhcs( 25) / +0.23191811872169963036326144682666D-20 /
  data asnhcs( 26) / -0.37487007953314343674570604543999D-21 /
  data asnhcs( 27) / +0.60732109822064279404549242880000D-22 /
  data asnhcs( 28) / -0.98599402764633583177370173440000D-23 /
  data asnhcs( 29) / +0.16039217452788496315232638293333D-23 /
  data asnhcs( 30) / -0.26138847350287686596716134399999D-24 /
  data asnhcs( 31) / +0.42670849606857390833358165333333D-25 /
  data asnhcs( 32) / -0.69770217039185243299730773333333D-26 /
  data asnhcs( 33) / +0.11425088336806858659812693333333D-26 /
  data asnhcs( 34) / -0.18735292078860968933021013333333D-27 /
  data asnhcs( 35) / +0.30763584414464922794065920000000D-28 /
  data asnhcs( 36) / -0.50577364031639824787046399999999D-29 /
  data asnhcs( 37) / +0.83250754712689142224213333333333D-30 /
  data asnhcs( 38) / -0.13718457282501044163925333333333D-30 /
  data asnhcs( 39) / +0.22629868426552784104106666666666D-31 /

  data aln2 / 0.69314718055994530941723212145818D+00 /
  data nterms / 0 /
  data xmax / 0.0D+00 /

  if ( nterms == 0 ) then
    nterms = r8_inits ( asnhcs, 39, 0.1D+00 * r8_mach ( 3 ) )
    sqeps = sqrt ( r8_mach ( 3 ) )
    xmax = 1.0D+00 / sqeps
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    value = x
  else if ( y <= 1.0D+00 ) then
    value = x * ( 1.0D+00 + r8_csevl ( 2.0D+00 * x * x - 1.0D+00, &
      asnhcs, nterms ) )
  else if ( y < xmax ) then
    value = log ( y + sqrt ( y * y + 1.0D+00 ) )
    if ( x < 0.0D+00 ) then
      value = - value
    end if
  else
    value = aln2 + log ( y )
    if ( x < 0.0D+00 ) then
      value = - value
    end if
  end if

  r8_asinh = value

  return
end
function r8_atan ( x )

!*****************************************************************************80
!
!! R8_ATAN evaluates the arc-tangent of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_ATAN, the arc-tangent of X.
!
  implicit none

  real ( kind = 8 ) atancs(16)
  real ( kind = 8 ) conpi8(4)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nterms
  real ( kind = 8 ) pi8(4)
  real ( kind = 8 ) r8_atan
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) tanp8(3)
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xbnd1
  real ( kind = 8 ) xbnd2
  real ( kind = 8 ) xbnd3
  real ( kind = 8 ) xbnd4
  real ( kind = 8 ) y

  save atancs
  save conpi8
  save nterms
  save pi8
  save tanp8
  save xbig
  save xbnd1
  save xbnd2
  save xbnd3
  save xbnd4

  data atancs(  1) / +0.48690110349241406474636915902891D+00 /
  data atancs(  2) / -0.65108316367174641818869794945974D-02 /
  data atancs(  3) / +0.38345828265245177653569992430456D-04 /
  data atancs(  4) / -0.26872212876223146539595410518788D-06 /
  data atancs(  5) / +0.20500930985824269846636514686688D-08 /
  data atancs(  6) / -0.16450717395484269455734135285348D-10 /
  data atancs(  7) / +0.13650975274390773423813528484428D-12 /
  data atancs(  8) / -0.11601779591998246322891309834666D-14 /
  data atancs(  9) / +0.10038333943866273835797657402666D-16 /
  data atancs( 10) / -0.88072747152163859327073696000000D-19 /
  data atancs( 11) / +0.78136321005661722180580266666666D-21 /
  data atancs( 12) / -0.69954535148267456086613333333333D-23 /
  data atancs( 13) / +0.63105905713702136004266666666666D-25 /
  data atancs( 14) / -0.57296075370213874346666666666666D-27 /
  data atancs( 15) / +0.52274796280602282666666666666666D-29 /
  data atancs( 16) / -0.48327903911608320000000000000000D-31 /

  data xbnd1 / +0.19891236737965800691159762264467D+00 /
  data xbnd2 / +0.66817863791929891999775768652308D+00 /
  data xbnd3 / +1.4966057626654890176011351349424D+00 /
  data xbnd4 / +5.0273394921258481045149750710640D+00 /

  data tanp8 (  1) / +0.41421356237309504880168872420969D+00 /
  data tanp8 (  2) / +1.0D+00 /
  data tanp8 (  3) / +2.4142135623730950488016887242096D+00 /

  data conpi8(1) / 0.375D+00 /
  data conpi8(2) / 0.75D+00 /
  data conpi8(3) / 1.125D+00 /
  data conpi8(4) / 1.5D+00 /

  data pi8(1) / +0.17699081698724154807830422909937D-01 /
  data pi8(2) / +0.35398163397448309615660845819875D-01 /
  data pi8(3) / +0.53097245096172464423491268729813D-01 /
  data pi8(4) / +0.70796326794896619231321691639751D-01 /

  data nterms / 0 /
  data xbig / 0.0D+00 /
 
  if ( nterms == 0 ) then
    nterms = r8_inits ( atancs, 16, 0.1D+00 * r8_mach ( 3 ) )
    sqeps = sqrt ( 6.0D+00 * r8_mach ( 3 ) )
    xbig = 1.0D+00 / r8_mach ( 3 )
  end if

  y = abs ( x )

  if ( y <= xbnd1 ) then

    value = x
    if ( sqeps < y ) then
      value = x * ( 0.75D+00 + r8_csevl ( &
        50.0D+00 * y * y - 1.0D+00, atancs, nterms ) )
    end if

  else if ( y <= xbnd4 ) then

    if ( xbnd3 < y ) then
      n = 3
    else if ( xbnd2 < y ) then
      n = 2
    else
      n = 1
    end if

    t = ( y - tanp8(n) ) / ( 1.0D+00 + y * tanp8(n) )

    value = conpi8(n) + ( pi8(n) + t * ( 0.75D+00 + &
      r8_csevl ( 50.0D+00 * t * t - 1.0D+00, atancs, nterms ) ) )

  else

    value = conpi8(4) + pi8(4)

    if ( y < xbig ) then
      value = conpi8(4) + ( pi8(4) - ( 0.75D+00 + &
        r8_csevl ( 50.0D+00 / y / y - 1.0D+00, atancs, &
        nterms ) ) / y )
    end if

  end if

  if ( x < 0.0D+00 ) then
    value = - abs ( value )
  else
    value = + abs ( value )
  end if

  r8_atan = value

  return
end
function r8_atan2 ( sn, cs )

!*****************************************************************************80
!
!! R8_ATAN2 evaluates the arc-tangent of two R8 arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) SN, CS, the Y and X coordinates of a 
!    point on the angle.
!
!    Output, real ( kind = 8 ) R8_ATAN2, the arc-tangent of the angle.
!
  implicit none

  real ( kind = 8 ) abscs
  real ( kind = 8 ) abssn
  real ( kind = 8 ) big
  real ( kind = 8 ) cs
  real ( kind = 8 ) pi
  real ( kind = 8 ) r8_atan2
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sml
  real ( kind = 8 ) sn

  save big
  save pi
  save sml

  data big / 0.0D+00 /
  data pi / 3.14159265358979323846264338327950D+00 /
  data sml / 0.0D+00 /

  if ( sml == 0.0D+00 ) then
    sml = r8_mach ( 1 )
    big = r8_mach ( 2 )
  end if
!
!  We now make sure SN can be divided by CS.  It is painful.
!
  abssn = abs ( sn )
  abscs = abs ( cs )

  if ( abscs <= abssn ) then

    if ( abscs < 1.0D+00 .and. abscs * big <= abssn ) then

      if ( sn < 0.0D+00 ) then
        r8_atan2 = - 0.5D+00 * pi
      else if ( sn == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_ATAN2 - Fatal error!'
        write ( *, '(a)' ) '  Both arguments are 0.'
        stop
      else
        r8_atan2 = 0.5D+00 * pi
      end if

      return

    end if

  else

    if ( 1.0D+00 < abscs .and. abssn <= abscs * sml ) then

      if ( 0.0D+00 <= cs ) then
        r8_atan2 = 0.0D+00
      else
        r8_atan2 = pi
      end if

      return

    end if

  end if

  r8_atan2 = atan ( sn / cs )

  if ( cs < 0.0D+00 ) then
    r8_atan2 = r8_atan2 + pi
  end if

  if ( pi < r8_atan2 ) then
    r8_atan2 = r8_atan2 - 2.0D+00 * pi
  end if

  return
end
function r8_atanh ( x )

!*****************************************************************************80
!
!! R8_ATANH evaluates the arc-hyperbolic tangent of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_ATANH, the arc-hyperbolic tangent of X.
!
  implicit none

  real ( kind = 8 ) atnhcs(27)
  real ( kind = 8 ) dxrel
  integer ( kind = 4 ) nterms
  real ( kind = 8 ) r8_atanh
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  save atnhcs
  save dxrel
  save nterms

  data atnhcs(  1) / +0.9439510239319549230842892218633D-01 /
  data atnhcs(  2) / +0.4919843705578615947200034576668D-01 /
  data atnhcs(  3) / +0.2102593522455432763479327331752D-02 /
  data atnhcs(  4) / +0.1073554449776116584640731045276D-03 /
  data atnhcs(  5) / +0.5978267249293031478642787517872D-05 /
  data atnhcs(  6) / +0.3505062030889134845966834886200D-06 /
  data atnhcs(  7) / +0.2126374343765340350896219314431D-07 /
  data atnhcs(  8) / +0.1321694535715527192129801723055D-08 /
  data atnhcs(  9) / +0.8365875501178070364623604052959D-10 /
  data atnhcs( 10) / +0.5370503749311002163881434587772D-11 /
  data atnhcs( 11) / +0.3486659470157107922971245784290D-12 /
  data atnhcs( 12) / +0.2284549509603433015524024119722D-13 /
  data atnhcs( 13) / +0.1508407105944793044874229067558D-14 /
  data atnhcs( 14) / +0.1002418816804109126136995722837D-15 /
  data atnhcs( 15) / +0.6698674738165069539715526882986D-17 /
  data atnhcs( 16) / +0.4497954546494931083083327624533D-18 /
  data atnhcs( 17) / +0.3032954474279453541682367146666D-19 /
  data atnhcs( 18) / +0.2052702064190936826463861418666D-20 /
  data atnhcs( 19) / +0.1393848977053837713193014613333D-21 /
  data atnhcs( 20) / +0.9492580637224576971958954666666D-23 /
  data atnhcs( 21) / +0.6481915448242307604982442666666D-24 /
  data atnhcs( 22) / +0.4436730205723615272632320000000D-25 /
  data atnhcs( 23) / +0.3043465618543161638912000000000D-26 /
  data atnhcs( 24) / +0.2091881298792393474047999999999D-27 /
  data atnhcs( 25) / +0.1440445411234050561365333333333D-28 /
  data atnhcs( 26) / +0.9935374683141640465066666666666D-30 /
  data atnhcs( 27) / +0.6863462444358260053333333333333D-31 /

  data dxrel / 0.0D+00 /
  data nterms / 0 /

  if ( nterms == 0 ) then
    nterms = r8_inits ( atnhcs, 27, 0.1D+00 * r8_mach ( 3 ) )
    dxrel = sqrt ( r8_mach ( 4 ) )
    sqeps = sqrt ( 3.0D+00 * r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    value = x
  else if ( y <= 0.5D+00 ) then
    value = x * ( 1.0D+00 + &
      r8_csevl ( 8.0D+00 * x * x - 1.0D+00, atnhcs, nterms ) )
  else if ( y < 1.0D+00 ) then
    value = 0.5D+00 * log ( ( 1.0D+00 + x ) / ( 1.0D+00 - x ) )
    if ( 1.0D+00 - y < dxrel ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_ATANH - Warning:'
      write ( *, '(a)' ) '  Answer lt half precision because |X| too near 1.'
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ATANH - Fatal error!'
    write ( *, '(a)' ) '  1 <= |X|.'
    stop
  end if

  r8_atanh = value

  return
end
subroutine r8_b0mp ( x, ampl, theta )

!*****************************************************************************80
!
!! R8_B0MP evaluates the modulus and phase for the Bessel J0 and Y0 functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) AMPL, THETA, the modulus and phase.
!
  implicit none

  real ( kind = 8 ) ampl
  real ( kind = 8 ) bm0cs(37)
  real ( kind = 8 ) bm02cs(40)
  real ( kind = 8 ) bt02cs(39)
  real ( kind = 8 ) bth0cs(44)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nbm0
  integer ( kind = 4 ) nbm02
  integer ( kind = 4 ) nbt02
  integer ( kind = 4 ) nbth0
  real ( kind = 8 ) pi4
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) z

  save bm0cs
  save bm02cs
  save bt02cs
  save bth0cs
  save nbm0
  save nbm02
  save nbt02
  save nbth0
  save xmax

  data bm0cs(  1) / +0.9211656246827742712573767730182D-01/
  data bm0cs(  2) / -0.1050590997271905102480716371755D-02/
  data bm0cs(  3) / +0.1470159840768759754056392850952D-04/
  data bm0cs(  4) / -0.5058557606038554223347929327702D-06/
  data bm0cs(  5) / +0.2787254538632444176630356137881D-07/
  data bm0cs(  6) / -0.2062363611780914802618841018973D-08/
  data bm0cs(  7) / +0.1870214313138879675138172596261D-09/
  data bm0cs(  8) / -0.1969330971135636200241730777825D-10/
  data bm0cs(  9) / +0.2325973793999275444012508818052D-11/
  data bm0cs( 10) / -0.3009520344938250272851224734482D-12/
  data bm0cs( 11) / +0.4194521333850669181471206768646D-13/
  data bm0cs( 12) / -0.6219449312188445825973267429564D-14/
  data bm0cs( 13) / +0.9718260411336068469601765885269D-15/
  data bm0cs( 14) / -0.1588478585701075207366635966937D-15/
  data bm0cs( 15) / +0.2700072193671308890086217324458D-16/
  data bm0cs( 16) / -0.4750092365234008992477504786773D-17/
  data bm0cs( 17) / +0.8615128162604370873191703746560D-18/
  data bm0cs( 18) / -0.1605608686956144815745602703359D-18/
  data bm0cs( 19) / +0.3066513987314482975188539801599D-19/
  data bm0cs( 20) / -0.5987764223193956430696505617066D-20/
  data bm0cs( 21) / +0.1192971253748248306489069841066D-20/
  data bm0cs( 22) / -0.2420969142044805489484682581333D-21/
  data bm0cs( 23) / +0.4996751760510616453371002879999D-22/
  data bm0cs( 24) / -0.1047493639351158510095040511999D-22/
  data bm0cs( 25) / +0.2227786843797468101048183466666D-23/
  data bm0cs( 26) / -0.4801813239398162862370542933333D-24/
  data bm0cs( 27) / +0.1047962723470959956476996266666D-24/
  data bm0cs( 28) / -0.2313858165678615325101260800000D-25/
  data bm0cs( 29) / +0.5164823088462674211635199999999D-26/
  data bm0cs( 30) / -0.1164691191850065389525401599999D-26/
  data bm0cs( 31) / +0.2651788486043319282958336000000D-27/
  data bm0cs( 32) / -0.6092559503825728497691306666666D-28/
  data bm0cs( 33) / +0.1411804686144259308038826666666D-28/
  data bm0cs( 34) / -0.3298094961231737245750613333333D-29/
  data bm0cs( 35) / +0.7763931143074065031714133333333D-30/
  data bm0cs( 36) / -0.1841031343661458478421333333333D-30/
  data bm0cs( 37) / +0.4395880138594310737100799999999D-31/

  data bth0cs(  1) / -0.24901780862128936717709793789967D+00/
  data bth0cs(  2) / +0.48550299609623749241048615535485D-03/
  data bth0cs(  3) / -0.54511837345017204950656273563505D-05/
  data bth0cs(  4) / +0.13558673059405964054377445929903D-06/
  data bth0cs(  5) / -0.55691398902227626227583218414920D-08/
  data bth0cs(  6) / +0.32609031824994335304004205719468D-09/
  data bth0cs(  7) / -0.24918807862461341125237903877993D-10/
  data bth0cs(  8) / +0.23449377420882520554352413564891D-11/
  data bth0cs(  9) / -0.26096534444310387762177574766136D-12/
  data bth0cs( 10) / +0.33353140420097395105869955014923D-13/
  data bth0cs( 11) / -0.47890000440572684646750770557409D-14/
  data bth0cs( 12) / +0.75956178436192215972642568545248D-15/
  data bth0cs( 13) / -0.13131556016891440382773397487633D-15/
  data bth0cs( 14) / +0.24483618345240857495426820738355D-16/
  data bth0cs( 15) / -0.48805729810618777683256761918331D-17/
  data bth0cs( 16) / +0.10327285029786316149223756361204D-17/
  data bth0cs( 17) / -0.23057633815057217157004744527025D-18/
  data bth0cs( 18) / +0.54044443001892693993017108483765D-19/
  data bth0cs( 19) / -0.13240695194366572724155032882385D-19/
  data bth0cs( 20) / +0.33780795621371970203424792124722D-20/
  data bth0cs( 21) / -0.89457629157111779003026926292299D-21/
  data bth0cs( 22) / +0.24519906889219317090899908651405D-21/
  data bth0cs( 23) / -0.69388422876866318680139933157657D-22/
  data bth0cs( 24) / +0.20228278714890138392946303337791D-22/
  data bth0cs( 25) / -0.60628500002335483105794195371764D-23/
  data bth0cs( 26) / +0.18649748964037635381823788396270D-23/
  data bth0cs( 27) / -0.58783732384849894560245036530867D-24/
  data bth0cs( 28) / +0.18958591447999563485531179503513D-24/
  data bth0cs( 29) / -0.62481979372258858959291620728565D-25/
  data bth0cs( 30) / +0.21017901684551024686638633529074D-25/
  data bth0cs( 31) / -0.72084300935209253690813933992446D-26/
  data bth0cs( 32) / +0.25181363892474240867156405976746D-26/
  data bth0cs( 33) / -0.89518042258785778806143945953643D-27/
  data bth0cs( 34) / +0.32357237479762298533256235868587D-27/
  data bth0cs( 35) / -0.11883010519855353657047144113796D-27/
  data bth0cs( 36) / +0.44306286907358104820579231941731D-28/
  data bth0cs( 37) / -0.16761009648834829495792010135681D-28/
  data bth0cs( 38) / +0.64292946921207466972532393966088D-29/
  data bth0cs( 39) / -0.24992261166978652421207213682763D-29/
  data bth0cs( 40) / +0.98399794299521955672828260355318D-30/
  data bth0cs( 41) / -0.39220375242408016397989131626158D-30/
  data bth0cs( 42) / +0.15818107030056522138590618845692D-30/
  data bth0cs( 43) / -0.64525506144890715944344098365426D-31/
  data bth0cs( 44) / +0.26611111369199356137177018346367D-31/

  data bm02cs(  1) / +0.9500415145228381369330861335560D-01/
  data bm02cs(  2) / -0.3801864682365670991748081566851D-03/
  data bm02cs(  3) / +0.2258339301031481192951829927224D-05/
  data bm02cs(  4) / -0.3895725802372228764730621412605D-07/
  data bm02cs(  5) / +0.1246886416512081697930990529725D-08/
  data bm02cs(  6) / -0.6065949022102503779803835058387D-10/
  data bm02cs(  7) / +0.4008461651421746991015275971045D-11/
  data bm02cs(  8) / -0.3350998183398094218467298794574D-12/
  data bm02cs(  9) / +0.3377119716517417367063264341996D-13/
  data bm02cs( 10) / -0.3964585901635012700569356295823D-14/
  data bm02cs( 11) / +0.5286111503883857217387939744735D-15/
  data bm02cs( 12) / -0.7852519083450852313654640243493D-16/
  data bm02cs( 13) / +0.1280300573386682201011634073449D-16/
  data bm02cs( 14) / -0.2263996296391429776287099244884D-17/
  data bm02cs( 15) / +0.4300496929656790388646410290477D-18/
  data bm02cs( 16) / -0.8705749805132587079747535451455D-19/
  data bm02cs( 17) / +0.1865862713962095141181442772050D-19/
  data bm02cs( 18) / -0.4210482486093065457345086972301D-20/
  data bm02cs( 19) / +0.9956676964228400991581627417842D-21/
  data bm02cs( 20) / -0.2457357442805313359605921478547D-21/
  data bm02cs( 21) / +0.6307692160762031568087353707059D-22/
  data bm02cs( 22) / -0.1678773691440740142693331172388D-22/
  data bm02cs( 23) / +0.4620259064673904433770878136087D-23/
  data bm02cs( 24) / -0.1311782266860308732237693402496D-23/
  data bm02cs( 25) / +0.3834087564116302827747922440276D-24/
  data bm02cs( 26) / -0.1151459324077741271072613293576D-24/
  data bm02cs( 27) / +0.3547210007523338523076971345213D-25/
  data bm02cs( 28) / -0.1119218385815004646264355942176D-25/
  data bm02cs( 29) / +0.3611879427629837831698404994257D-26/
  data bm02cs( 30) / -0.1190687765913333150092641762463D-26/
  data bm02cs( 31) / +0.4005094059403968131802476449536D-27/
  data bm02cs( 32) / -0.1373169422452212390595193916017D-27/
  data bm02cs( 33) / +0.4794199088742531585996491526437D-28/
  data bm02cs( 34) / -0.1702965627624109584006994476452D-28/
  data bm02cs( 35) / +0.6149512428936330071503575161324D-29/
  data bm02cs( 36) / -0.2255766896581828349944300237242D-29/
  data bm02cs( 37) / +0.8399707509294299486061658353200D-30/
  data bm02cs( 38) / -0.3172997595562602355567423936152D-30/
  data bm02cs( 39) / +0.1215205298881298554583333026514D-30/
  data bm02cs( 40) / -0.4715852749754438693013210568045D-31/

  data bt02cs(  1) / -0.24548295213424597462050467249324D+00/
  data bt02cs(  2) / +0.12544121039084615780785331778299D-02/
  data bt02cs(  3) / -0.31253950414871522854973446709571D-04/
  data bt02cs(  4) / +0.14709778249940831164453426969314D-05/
  data bt02cs(  5) / -0.99543488937950033643468850351158D-07/
  data bt02cs(  6) / +0.85493166733203041247578711397751D-08/
  data bt02cs(  7) / -0.86989759526554334557985512179192D-09/
  data bt02cs(  8) / +0.10052099533559791084540101082153D-09/
  data bt02cs(  9) / -0.12828230601708892903483623685544D-10/
  data bt02cs( 10) / +0.17731700781805131705655750451023D-11/
  data bt02cs( 11) / -0.26174574569485577488636284180925D-12/
  data bt02cs( 12) / +0.40828351389972059621966481221103D-13/
  data bt02cs( 13) / -0.66751668239742720054606749554261D-14/
  data bt02cs( 14) / +0.11365761393071629448392469549951D-14/
  data bt02cs( 15) / -0.20051189620647160250559266412117D-15/
  data bt02cs( 16) / +0.36497978794766269635720591464106D-16/
  data bt02cs( 17) / -0.68309637564582303169355843788800D-17/
  data bt02cs( 18) / +0.13107583145670756620057104267946D-17/
  data bt02cs( 19) / -0.25723363101850607778757130649599D-18/
  data bt02cs( 20) / +0.51521657441863959925267780949333D-19/
  data bt02cs( 21) / -0.10513017563758802637940741461333D-19/
  data bt02cs( 22) / +0.21820381991194813847301084501333D-20/
  data bt02cs( 23) / -0.46004701210362160577225905493333D-21/
  data bt02cs( 24) / +0.98407006925466818520953651199999D-22/
  data bt02cs( 25) / -0.21334038035728375844735986346666D-22/
  data bt02cs( 26) / +0.46831036423973365296066286933333D-23/
  data bt02cs( 27) / -0.10400213691985747236513382399999D-23/
  data bt02cs( 28) / +0.23349105677301510051777740800000D-24/
  data bt02cs( 29) / -0.52956825323318615788049749333333D-25/
  data bt02cs( 30) / +0.12126341952959756829196287999999D-25/
  data bt02cs( 31) / -0.28018897082289428760275626666666D-26/
  data bt02cs( 32) / +0.65292678987012873342593706666666D-27/
  data bt02cs( 33) / -0.15337980061873346427835733333333D-27/
  data bt02cs( 34) / +0.36305884306364536682359466666666D-28/
  data bt02cs( 35) / -0.86560755713629122479172266666666D-29/
  data bt02cs( 36) / +0.20779909972536284571238399999999D-29/
  data bt02cs( 37) / -0.50211170221417221674325333333333D-30/
  data bt02cs( 38) / +0.12208360279441714184191999999999D-30/
  data bt02cs( 39) / -0.29860056267039913454250666666666D-31/

  data nbm0 / 0 /
  data nbm02 / 0 /
  data nbt02 / 0 /
  data nbth0 / 0 /
  data pi4 / 0.785398163397448309615660845819876D+00 /
  data xmax / 0.0D+00 /

  if ( nbm0 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nbm0 = r8_inits ( bm0cs, 37, eta )
    nbt02 = r8_inits ( bt02cs, 39, eta )
    nbm02 = r8_inits ( bm02cs, 40, eta )
    nbth0 = r8_inits ( bth0cs, 44, eta )
    xmax = 1.0D+00 / r8_mach ( 4 )
  end if

  if ( x < 4.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B0MP - Fatal error!'
    write ( *, '(a)' ) '  X < 4.'
    stop
  else if ( x <= 8.0D+00 ) then
    z = ( 128.0D+00 / x / x - 5.0D+00 ) / 3.0D+00
    ampl = ( 0.75D+00 + r8_csevl ( z, bm0cs, nbm0 ) ) / sqrt ( x )
    theta = x - pi4 + r8_csevl ( z, bt02cs, nbt02 ) / x
  else
    z = 128.0D+00 / x / x - 1.0D+00
    ampl = ( 0.75D+00 + r8_csevl ( z, bm02cs, nbm02) ) / sqrt ( x )
    theta = x - pi4 + r8_csevl ( z, bth0cs, nbth0 ) / x
  end if

  return
end
subroutine r8_b1mp ( x, ampl, theta )

!*****************************************************************************80
!
!! R8_B1MP evaluates the modulus and phase for the Bessel J1 and Y1 functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) AMPL, THETA, the modulus and phase.
!
  implicit none

  real ( kind = 8 ) ampl
  real ( kind = 8 ) bm12cs(40)
  real ( kind = 8 ) bm1cs(37)
  real ( kind = 8 ) bt12cs(39)
  real ( kind = 8 ) bth1cs(44)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nbm1
  integer ( kind = 4 ) nbm12
  integer ( kind = 4 ) nbt12
  integer ( kind = 4 ) nbth1
  real ( kind = 8 ) pi4
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) z

  save bm12cs
  save bm1cs
  save bt12cs
  save bth1cs
  save nbm1
  save nbm12
  save nbt12
  save nbth1
  save pi4
  save xmax

  data bm1cs(  1) / +0.1069845452618063014969985308538D+00/
  data bm1cs(  2) / +0.3274915039715964900729055143445D-02/
  data bm1cs(  3) / -0.2987783266831698592030445777938D-04/
  data bm1cs(  4) / +0.8331237177991974531393222669023D-06/
  data bm1cs(  5) / -0.4112665690302007304896381725498D-07/
  data bm1cs(  6) / +0.2855344228789215220719757663161D-08/
  data bm1cs(  7) / -0.2485408305415623878060026596055D-09/
  data bm1cs(  8) / +0.2543393338072582442742484397174D-10/
  data bm1cs(  9) / -0.2941045772822967523489750827909D-11/
  data bm1cs( 10) / +0.3743392025493903309265056153626D-12/
  data bm1cs( 11) / -0.5149118293821167218720548243527D-13/
  data bm1cs( 12) / +0.7552535949865143908034040764199D-14/
  data bm1cs( 13) / -0.1169409706828846444166290622464D-14/
  data bm1cs( 14) / +0.1896562449434791571721824605060D-15/
  data bm1cs( 15) / -0.3201955368693286420664775316394D-16/
  data bm1cs( 16) / +0.5599548399316204114484169905493D-17/
  data bm1cs( 17) / -0.1010215894730432443119390444544D-17/
  data bm1cs( 18) / +0.1873844985727562983302042719573D-18/
  data bm1cs( 19) / -0.3563537470328580219274301439999D-19/
  data bm1cs( 20) / +0.6931283819971238330422763519999D-20/
  data bm1cs( 21) / -0.1376059453406500152251408930133D-20/
  data bm1cs( 22) / +0.2783430784107080220599779327999D-21/
  data bm1cs( 23) / -0.5727595364320561689348669439999D-22/
  data bm1cs( 24) / +0.1197361445918892672535756799999D-22/
  data bm1cs( 25) / -0.2539928509891871976641440426666D-23/
  data bm1cs( 26) / +0.5461378289657295973069619199999D-24/
  data bm1cs( 27) / -0.1189211341773320288986289493333D-24/
  data bm1cs( 28) / +0.2620150977340081594957824000000D-25/
  data bm1cs( 29) / -0.5836810774255685901920938666666D-26/
  data bm1cs( 30) / +0.1313743500080595773423615999999D-26/
  data bm1cs( 31) / -0.2985814622510380355332778666666D-27/
  data bm1cs( 32) / +0.6848390471334604937625599999999D-28/
  data bm1cs( 33) / -0.1584401568222476721192960000000D-28/
  data bm1cs( 34) / +0.3695641006570938054301013333333D-29/
  data bm1cs( 35) / -0.8687115921144668243012266666666D-30/
  data bm1cs( 36) / +0.2057080846158763462929066666666D-30/
  data bm1cs( 37) / -0.4905225761116225518523733333333D-31/

  data bt12cs(  1) / +0.73823860128742974662620839792764D+00/
  data bt12cs(  2) / -0.33361113174483906384470147681189D-02/
  data bt12cs(  3) / +0.61463454888046964698514899420186D-04/
  data bt12cs(  4) / -0.24024585161602374264977635469568D-05/
  data bt12cs(  5) / +0.14663555577509746153210591997204D-06/
  data bt12cs(  6) / -0.11841917305589180567005147504983D-07/
  data bt12cs(  7) / +0.11574198963919197052125466303055D-08/
  data bt12cs(  8) / -0.13001161129439187449366007794571D-09/
  data bt12cs(  9) / +0.16245391141361731937742166273667D-10/
  data bt12cs( 10) / -0.22089636821403188752155441770128D-11/
  data bt12cs( 11) / +0.32180304258553177090474358653778D-12/
  data bt12cs( 12) / -0.49653147932768480785552021135381D-13/
  data bt12cs( 13) / +0.80438900432847825985558882639317D-14/
  data bt12cs( 14) / -0.13589121310161291384694712682282D-14/
  data bt12cs( 15) / +0.23810504397147214869676529605973D-15/
  data bt12cs( 16) / -0.43081466363849106724471241420799D-16/
  data bt12cs( 17) / +0.80202544032771002434993512550400D-17/
  data bt12cs( 18) / -0.15316310642462311864230027468799D-17/
  data bt12cs( 19) / +0.29928606352715568924073040554666D-18/
  data bt12cs( 20) / -0.59709964658085443393815636650666D-19/
  data bt12cs( 21) / +0.12140289669415185024160852650666D-19/
  data bt12cs( 22) / -0.25115114696612948901006977706666D-20/
  data bt12cs( 23) / +0.52790567170328744850738380799999D-21/
  data bt12cs( 24) / -0.11260509227550498324361161386666D-21/
  data bt12cs( 25) / +0.24348277359576326659663462400000D-22/
  data bt12cs( 26) / -0.53317261236931800130038442666666D-23/
  data bt12cs( 27) / +0.11813615059707121039205990399999D-23/
  data bt12cs( 28) / -0.26465368283353523514856789333333D-24/
  data bt12cs( 29) / +0.59903394041361503945577813333333D-25/
  data bt12cs( 30) / -0.13690854630829503109136383999999D-25/
  data bt12cs( 31) / +0.31576790154380228326413653333333D-26/
  data bt12cs( 32) / -0.73457915082084356491400533333333D-27/
  data bt12cs( 33) / +0.17228081480722747930705920000000D-27/
  data bt12cs( 34) / -0.40716907961286507941068800000000D-28/
  data bt12cs( 35) / +0.96934745136779622700373333333333D-29/
  data bt12cs( 36) / -0.23237636337765716765354666666666D-29/
  data bt12cs( 37) / +0.56074510673522029406890666666666D-30/
  data bt12cs( 38) / -0.13616465391539005860522666666666D-30/
  data bt12cs( 39) / +0.33263109233894654388906666666666D-31/

  data bm12cs(  1) / +0.9807979156233050027272093546937D-01/
  data bm12cs(  2) / +0.1150961189504685306175483484602D-02/
  data bm12cs(  3) / -0.4312482164338205409889358097732D-05/
  data bm12cs(  4) / +0.5951839610088816307813029801832D-07/
  data bm12cs(  5) / -0.1704844019826909857400701586478D-08/
  data bm12cs(  6) / +0.7798265413611109508658173827401D-10/
  data bm12cs(  7) / -0.4958986126766415809491754951865D-11/
  data bm12cs(  8) / +0.4038432416421141516838202265144D-12/
  data bm12cs(  9) / -0.3993046163725175445765483846645D-13/
  data bm12cs( 10) / +0.4619886183118966494313342432775D-14/
  data bm12cs( 11) / -0.6089208019095383301345472619333D-15/
  data bm12cs( 12) / +0.8960930916433876482157048041249D-16/
  data bm12cs( 13) / -0.1449629423942023122916518918925D-16/
  data bm12cs( 14) / +0.2546463158537776056165149648068D-17/
  data bm12cs( 15) / -0.4809472874647836444259263718620D-18/
  data bm12cs( 16) / +0.9687684668292599049087275839124D-19/
  data bm12cs( 17) / -0.2067213372277966023245038117551D-19/
  data bm12cs( 18) / +0.4646651559150384731802767809590D-20/
  data bm12cs( 19) / -0.1094966128848334138241351328339D-20/
  data bm12cs( 20) / +0.2693892797288682860905707612785D-21/
  data bm12cs( 21) / -0.6894992910930374477818970026857D-22/
  data bm12cs( 22) / +0.1830268262752062909890668554740D-22/
  data bm12cs( 23) / -0.5025064246351916428156113553224D-23/
  data bm12cs( 24) / +0.1423545194454806039631693634194D-23/
  data bm12cs( 25) / -0.4152191203616450388068886769801D-24/
  data bm12cs( 26) / +0.1244609201503979325882330076547D-24/
  data bm12cs( 27) / -0.3827336370569304299431918661286D-25/
  data bm12cs( 28) / +0.1205591357815617535374723981835D-25/
  data bm12cs( 29) / -0.3884536246376488076431859361124D-26/
  data bm12cs( 30) / +0.1278689528720409721904895283461D-26/
  data bm12cs( 31) / -0.4295146689447946272061936915912D-27/
  data bm12cs( 32) / +0.1470689117829070886456802707983D-27/
  data bm12cs( 33) / -0.5128315665106073128180374017796D-28/
  data bm12cs( 34) / +0.1819509585471169385481437373286D-28/
  data bm12cs( 35) / -0.6563031314841980867618635050373D-29/
  data bm12cs( 36) / +0.2404898976919960653198914875834D-29/
  data bm12cs( 37) / -0.8945966744690612473234958242979D-30/
  data bm12cs( 38) / +0.3376085160657231026637148978240D-30/
  data bm12cs( 39) / -0.1291791454620656360913099916966D-30/
  data bm12cs( 40) / +0.5008634462958810520684951501254D-31/

  data bth1cs(  1) / +0.74749957203587276055443483969695D+00/
  data bth1cs(  2) / -0.12400777144651711252545777541384D-02/
  data bth1cs(  3) / +0.99252442404424527376641497689592D-05/
  data bth1cs(  4) / -0.20303690737159711052419375375608D-06/
  data bth1cs(  5) / +0.75359617705690885712184017583629D-08/
  data bth1cs(  6) / -0.41661612715343550107630023856228D-09/
  data bth1cs(  7) / +0.30701618070834890481245102091216D-10/
  data bth1cs(  8) / -0.28178499637605213992324008883924D-11/
  data bth1cs(  9) / +0.30790696739040295476028146821647D-12/
  data bth1cs( 10) / -0.38803300262803434112787347554781D-13/
  data bth1cs( 11) / +0.55096039608630904934561726208562D-14/
  data bth1cs( 12) / -0.86590060768383779940103398953994D-15/
  data bth1cs( 13) / +0.14856049141536749003423689060683D-15/
  data bth1cs( 14) / -0.27519529815904085805371212125009D-16/
  data bth1cs( 15) / +0.54550796090481089625036223640923D-17/
  data bth1cs( 16) / -0.11486534501983642749543631027177D-17/
  data bth1cs( 17) / +0.25535213377973900223199052533522D-18/
  data bth1cs( 18) / -0.59621490197413450395768287907849D-19/
  data bth1cs( 19) / +0.14556622902372718620288302005833D-19/
  data bth1cs( 20) / -0.37022185422450538201579776019593D-20/
  data bth1cs( 21) / +0.97763074125345357664168434517924D-21/
  data bth1cs( 22) / -0.26726821639668488468723775393052D-21/
  data bth1cs( 23) / +0.75453300384983271794038190655764D-22/
  data bth1cs( 24) / -0.21947899919802744897892383371647D-22/
  data bth1cs( 25) / +0.65648394623955262178906999817493D-23/
  data bth1cs( 26) / -0.20155604298370207570784076869519D-23/
  data bth1cs( 27) / +0.63417768556776143492144667185670D-24/
  data bth1cs( 28) / -0.20419277885337895634813769955591D-24/
  data bth1cs( 29) / +0.67191464220720567486658980018551D-25/
  data bth1cs( 30) / -0.22569079110207573595709003687336D-25/
  data bth1cs( 31) / +0.77297719892989706370926959871929D-26/
  data bth1cs( 32) / -0.26967444512294640913211424080920D-26/
  data bth1cs( 33) / +0.95749344518502698072295521933627D-27/
  data bth1cs( 34) / -0.34569168448890113000175680827627D-27/
  data bth1cs( 35) / +0.12681234817398436504211986238374D-27/
  data bth1cs( 36) / -0.47232536630722639860464993713445D-28/
  data bth1cs( 37) / +0.17850008478186376177858619796417D-28/
  data bth1cs( 38) / -0.68404361004510395406215223566746D-29/
  data bth1cs( 39) / +0.26566028671720419358293422672212D-29/
  data bth1cs( 40) / -0.10450402527914452917714161484670D-29/
  data bth1cs( 41) / +0.41618290825377144306861917197064D-30/
  data bth1cs( 42) / -0.16771639203643714856501347882887D-30/
  data bth1cs( 43) / +0.68361997776664389173535928028528D-31/
  data bth1cs( 44) / -0.28172247861233641166739574622810D-31/

  data nbm1 / 0 /
  data nbm12 / 0 /
  data nbt12 / 0 /
  data nbth1 / 0 /
  data pi4 / 0.785398163397448309615660845819876D+00 /
  data xmax / 0.0D+00 /

  if ( nbm1 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nbm1 = r8_inits ( bm1cs, 37, eta )
    nbt12 = r8_inits ( bt12cs, 39, eta )
    nbm12 = r8_inits ( bm12cs, 40, eta )
    nbth1 = r8_inits ( bth1cs, 44, eta )
    xmax = 1.0D+00 / r8_mach ( 4 )
  end if

  if ( x < 4.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B1MP - Fatal error!'
    write ( *, '(a)' ) '  X < 4.'
    stop
  else if ( x <= 8.0D+00 ) then
    z = ( 128.0D+00 / x / x - 5.0D+00 ) / 3.0D+00
    ampl = ( 0.75D+00 + r8_csevl ( z, bm1cs, nbm1 ) ) / sqrt ( x )
    theta = x - 3.0D+00 * pi4 + r8_csevl ( z, bt12cs, nbt12 ) / x
  else
    z = 128.0D+00 / x / x - 1.0D+00
    ampl = ( 0.75D+00 + r8_csevl ( z, bm12cs, nbm12 ) ) / sqrt ( x )
    theta = x - 3.0D+00 * pi4 + r8_csevl ( z, bth1cs, nbth1 ) / x
  end if

  return
end
function r8_besi0 ( x )

!*****************************************************************************80
!
!! R8_BESI0 evaluates the Bessel function I of order 0 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESI0, the Bessel function I of order 0 of X.
!
  implicit none

  real ( kind = 8 ) bi0cs(18)
  integer ( kind = 4 ) nti0
  real ( kind = 8 ) r8_besi0
  real ( kind = 8 ) r8_besi0e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save bi0cs
  save nti0
  save xmax
  save xsml

  data bi0cs(  1) / -0.7660547252839144951081894976243285D-01 /
  data bi0cs(  2) / +0.1927337953993808269952408750881196D+01 /
  data bi0cs(  3) / +0.2282644586920301338937029292330415D+00 /
  data bi0cs(  4) / +0.1304891466707290428079334210691888D-01 /
  data bi0cs(  5) / +0.4344270900816487451378682681026107D-03 /
  data bi0cs(  6) / +0.9422657686001934663923171744118766D-05 /
  data bi0cs(  7) / +0.1434006289510691079962091878179957D-06 /
  data bi0cs(  8) / +0.1613849069661749069915419719994611D-08 /
  data bi0cs(  9) / +0.1396650044535669699495092708142522D-10 /
  data bi0cs( 10) / +0.9579451725505445344627523171893333D-13 /
  data bi0cs( 11) / +0.5333981859862502131015107744000000D-15 /
  data bi0cs( 12) / +0.2458716088437470774696785919999999D-17 /
  data bi0cs( 13) / +0.9535680890248770026944341333333333D-20 /
  data bi0cs( 14) / +0.3154382039721427336789333333333333D-22 /
  data bi0cs( 15) / +0.9004564101094637431466666666666666D-25 /
  data bi0cs( 16) / +0.2240647369123670016000000000000000D-27 /
  data bi0cs( 17) / +0.4903034603242837333333333333333333D-30 /
  data bi0cs( 18) / +0.9508172606122666666666666666666666D-33 /

  data nti0 / 0 /
  data xmax / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( nti0 == 0 ) then
    nti0 = r8_inits ( bi0cs, 18, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( 8.0D+00 * r8_mach ( 3 ) )
    xmax = log ( r8_mach ( 2 ) )
  end if

  y = abs ( x )

  if ( y <= xsml ) then
    r8_besi0 = 1.0D+00
  else if ( y <= 3.0D+00 ) then
    r8_besi0 = 2.75D+00 + r8_csevl ( y * y / 4.5D+00 - 1.0D+00, bi0cs, nti0 )
  else if ( y <= xmax ) then
    r8_besi0 = exp ( y ) * r8_besi0e ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESI0 - Fatal error!'
    write ( *, '(a)' ) '  |X| too large.'
    stop
  end if

  return
end
function r8_besi0e ( x )

!*****************************************************************************80
!
!! R8_BESI0E evaluates the exponentially scaled Bessel function I0(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESI0E, the exponentially scaled Bessel 
!    function I0(X).
!
  implicit none

  real ( kind = 8 ) ai02cs(69)
  real ( kind = 8 ) ai0cs(46)
  real ( kind = 8 ) bi0cs(18)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) ntai02
  integer ( kind = 4 ) ntai0
  integer ( kind = 4 ) nti0
  real ( kind = 8 ) r8_besi0e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save ai02cs
  save ai0cs
  save bi0cs
  save ntai02
  save ntai0
  save nti0
  save xsml

  data bi0cs(  1) / -0.7660547252839144951081894976243285D-01 /
  data bi0cs(  2) / +0.1927337953993808269952408750881196D+01 /
  data bi0cs(  3) / +0.2282644586920301338937029292330415D+00 /
  data bi0cs(  4) / +0.1304891466707290428079334210691888D-01 /
  data bi0cs(  5) / +0.4344270900816487451378682681026107D-03 /
  data bi0cs(  6) / +0.9422657686001934663923171744118766D-05 /
  data bi0cs(  7) / +0.1434006289510691079962091878179957D-06 /
  data bi0cs(  8) / +0.1613849069661749069915419719994611D-08 /
  data bi0cs(  9) / +0.1396650044535669699495092708142522D-10 /
  data bi0cs( 10) / +0.9579451725505445344627523171893333D-13 /
  data bi0cs( 11) / +0.5333981859862502131015107744000000D-15 /
  data bi0cs( 12) / +0.2458716088437470774696785919999999D-17 /
  data bi0cs( 13) / +0.9535680890248770026944341333333333D-20 /
  data bi0cs( 14) / +0.3154382039721427336789333333333333D-22 /
  data bi0cs( 15) / +0.9004564101094637431466666666666666D-25 /
  data bi0cs( 16) / +0.2240647369123670016000000000000000D-27 /
  data bi0cs( 17) / +0.4903034603242837333333333333333333D-30 /
  data bi0cs( 18) / +0.9508172606122666666666666666666666D-33 /

  data ai0cs(  1) / +0.7575994494023795942729872037438D-01 /
  data ai0cs(  2) / +0.7591380810823345507292978733204D-02 /
  data ai0cs(  3) / +0.4153131338923750501863197491382D-03 /
  data ai0cs(  4) / +0.1070076463439073073582429702170D-04 /
  data ai0cs(  5) / -0.7901179979212894660750319485730D-05 /
  data ai0cs(  6) / -0.7826143501438752269788989806909D-06 /
  data ai0cs(  7) / +0.2783849942948870806381185389857D-06 /
  data ai0cs(  8) / +0.8252472600612027191966829133198D-08 /
  data ai0cs(  9) / -0.1204463945520199179054960891103D-07 /
  data ai0cs( 10) / +0.1559648598506076443612287527928D-08 /
  data ai0cs( 11) / +0.2292556367103316543477254802857D-09 /
  data ai0cs( 12) / -0.1191622884279064603677774234478D-09 /
  data ai0cs( 13) / +0.1757854916032409830218331247743D-10 /
  data ai0cs( 14) / +0.1128224463218900517144411356824D-11 /
  data ai0cs( 15) / -0.1146848625927298877729633876982D-11 /
  data ai0cs( 16) / +0.2715592054803662872643651921606D-12 /
  data ai0cs( 17) / -0.2415874666562687838442475720281D-13 /
  data ai0cs( 18) / -0.6084469888255125064606099639224D-14 /
  data ai0cs( 19) / +0.3145705077175477293708360267303D-14 /
  data ai0cs( 20) / -0.7172212924871187717962175059176D-15 /
  data ai0cs( 21) / +0.7874493403454103396083909603327D-16 /
  data ai0cs( 22) / +0.1004802753009462402345244571839D-16 /
  data ai0cs( 23) / -0.7566895365350534853428435888810D-17 /
  data ai0cs( 24) / +0.2150380106876119887812051287845D-17 /
  data ai0cs( 25) / -0.3754858341830874429151584452608D-18 /
  data ai0cs( 26) / +0.2354065842226992576900757105322D-19 /
  data ai0cs( 27) / +0.1114667612047928530226373355110D-19 /
  data ai0cs( 28) / -0.5398891884396990378696779322709D-20 /
  data ai0cs( 29) / +0.1439598792240752677042858404522D-20 /
  data ai0cs( 30) / -0.2591916360111093406460818401962D-21 /
  data ai0cs( 31) / +0.2238133183998583907434092298240D-22 /
  data ai0cs( 32) / +0.5250672575364771172772216831999D-23 /
  data ai0cs( 33) / -0.3249904138533230784173432285866D-23 /
  data ai0cs( 34) / +0.9924214103205037927857284710400D-24 /
  data ai0cs( 35) / -0.2164992254244669523146554299733D-24 /
  data ai0cs( 36) / +0.3233609471943594083973332991999D-25 /
  data ai0cs( 37) / -0.1184620207396742489824733866666D-26 /
  data ai0cs( 38) / -0.1281671853950498650548338687999D-26 /
  data ai0cs( 39) / +0.5827015182279390511605568853333D-27 /
  data ai0cs( 40) / -0.1668222326026109719364501503999D-27 /
  data ai0cs( 41) / +0.3625309510541569975700684800000D-28 /
  data ai0cs( 42) / -0.5733627999055713589945958399999D-29 /
  data ai0cs( 43) / +0.3736796722063098229642581333333D-30 /
  data ai0cs( 44) / +0.1602073983156851963365512533333D-30 /
  data ai0cs( 45) / -0.8700424864057229884522495999999D-31 /
  data ai0cs( 46) / +0.2741320937937481145603413333333D-31 /

  data ai02cs(  1) / +0.5449041101410883160789609622680D-01 /
  data ai02cs(  2) / +0.3369116478255694089897856629799D-02 /
  data ai02cs(  3) / +0.6889758346916823984262639143011D-04 /
  data ai02cs(  4) / +0.2891370520834756482966924023232D-05 /
  data ai02cs(  5) / +0.2048918589469063741827605340931D-06 /
  data ai02cs(  6) / +0.2266668990498178064593277431361D-07 /
  data ai02cs(  7) / +0.3396232025708386345150843969523D-08 /
  data ai02cs(  8) / +0.4940602388224969589104824497835D-09 /
  data ai02cs(  9) / +0.1188914710784643834240845251963D-10 /
  data ai02cs( 10) / -0.3149916527963241364538648629619D-10 /
  data ai02cs( 11) / -0.1321581184044771311875407399267D-10 /
  data ai02cs( 12) / -0.1794178531506806117779435740269D-11 /
  data ai02cs( 13) / +0.7180124451383666233671064293469D-12 /
  data ai02cs( 14) / +0.3852778382742142701140898017776D-12 /
  data ai02cs( 15) / +0.1540086217521409826913258233397D-13 /
  data ai02cs( 16) / -0.4150569347287222086626899720156D-13 /
  data ai02cs( 17) / -0.9554846698828307648702144943125D-14 /
  data ai02cs( 18) / +0.3811680669352622420746055355118D-14 /
  data ai02cs( 19) / +0.1772560133056526383604932666758D-14 /
  data ai02cs( 20) / -0.3425485619677219134619247903282D-15 /
  data ai02cs( 21) / -0.2827623980516583484942055937594D-15 /
  data ai02cs( 22) / +0.3461222867697461093097062508134D-16 /
  data ai02cs( 23) / +0.4465621420296759999010420542843D-16 /
  data ai02cs( 24) / -0.4830504485944182071255254037954D-17 /
  data ai02cs( 25) / -0.7233180487874753954562272409245D-17 /
  data ai02cs( 26) / +0.9921475412173698598880460939810D-18 /
  data ai02cs( 27) / +0.1193650890845982085504399499242D-17 /
  data ai02cs( 28) / -0.2488709837150807235720544916602D-18 /
  data ai02cs( 29) / -0.1938426454160905928984697811326D-18 /
  data ai02cs( 30) / +0.6444656697373443868783019493949D-19 /
  data ai02cs( 31) / +0.2886051596289224326481713830734D-19 /
  data ai02cs( 32) / -0.1601954907174971807061671562007D-19 /
  data ai02cs( 33) / -0.3270815010592314720891935674859D-20 /
  data ai02cs( 34) / +0.3686932283826409181146007239393D-20 /
  data ai02cs( 35) / +0.1268297648030950153013595297109D-22 /
  data ai02cs( 36) / -0.7549825019377273907696366644101D-21 /
  data ai02cs( 37) / +0.1502133571377835349637127890534D-21 /
  data ai02cs( 38) / +0.1265195883509648534932087992483D-21 /
  data ai02cs( 39) / -0.6100998370083680708629408916002D-22 /
  data ai02cs( 40) / -0.1268809629260128264368720959242D-22 /
  data ai02cs( 41) / +0.1661016099890741457840384874905D-22 /
  data ai02cs( 42) / -0.1585194335765885579379705048814D-23 /
  data ai02cs( 43) / -0.3302645405968217800953817667556D-23 /
  data ai02cs( 44) / +0.1313580902839239781740396231174D-23 /
  data ai02cs( 45) / +0.3689040246671156793314256372804D-24 /
  data ai02cs( 46) / -0.4210141910461689149219782472499D-24 /
  data ai02cs( 47) / +0.4791954591082865780631714013730D-25 /
  data ai02cs( 48) / +0.8459470390221821795299717074124D-25 /
  data ai02cs( 49) / -0.4039800940872832493146079371810D-25 /
  data ai02cs( 50) / -0.6434714653650431347301008504695D-26 /
  data ai02cs( 51) / +0.1225743398875665990344647369905D-25 /
  data ai02cs( 52) / -0.2934391316025708923198798211754D-26 /
  data ai02cs( 53) / -0.1961311309194982926203712057289D-26 /
  data ai02cs( 54) / +0.1503520374822193424162299003098D-26 /
  data ai02cs( 55) / -0.9588720515744826552033863882069D-28 /
  data ai02cs( 56) / -0.3483339380817045486394411085114D-27 /
  data ai02cs( 57) / +0.1690903610263043673062449607256D-27 /
  data ai02cs( 58) / +0.1982866538735603043894001157188D-28 /
  data ai02cs( 59) / -0.5317498081491816214575830025284D-28 /
  data ai02cs( 60) / +0.1803306629888392946235014503901D-28 /
  data ai02cs( 61) / +0.6213093341454893175884053112422D-29 /
  data ai02cs( 62) / -0.7692189292772161863200728066730D-29 /
  data ai02cs( 63) / +0.1858252826111702542625560165963D-29 /
  data ai02cs( 64) / +0.1237585142281395724899271545541D-29 /
  data ai02cs( 65) / -0.1102259120409223803217794787792D-29 /
  data ai02cs( 66) / +0.1886287118039704490077874479431D-30 /
  data ai02cs( 67) / +0.2160196872243658913149031414060D-30 /
  data ai02cs( 68) / -0.1605454124919743200584465949655D-30 /
  data ai02cs( 69) / +0.1965352984594290603938848073318D-31 /

  data ntai0 / 0 /
  data ntai02 / 0 /
  data nti0 / 0 /
  data xsml / 0.0D+00 /

  if ( nti0 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nti0 = r8_inits ( bi0cs, 18, eta )
    ntai0 = r8_inits ( ai0cs, 46, eta )
    ntai02 = r8_inits ( ai02cs, 69, eta )
    xsml = sqrt ( 8.0D+00 * r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= xsml ) then
    r8_besi0e = 1.0D+00
  else if ( y <= 3.0D+00 ) then
    r8_besi0e = exp ( - y ) * ( 2.75D+00 &
      + r8_csevl ( y * y / 4.5D+00 - 1.0D+00, bi0cs, nti0 ) )
  else if ( y <= 8.0D+00 ) then
    r8_besi0e = ( 0.375D+00 &
      + r8_csevl ( ( 48.0D+00 / y - 11.0D+00 ) / 5.0D+00, &
      ai0cs, ntai0 ) ) / sqrt ( y )
  else
    r8_besi0e = ( 0.375D+00 &
      + r8_csevl ( 16.0D+00 / y - 1.0D+00, ai02cs, ntai02 ) ) &
      / sqrt ( y )
  end if

  return
end
function r8_besi1 ( x )

!*****************************************************************************80
!
!! R8_BESI1 evaluates the Bessel function I of order 1 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESI1, the Bessel function I of order 1 of X.
!
  implicit none

  real ( kind = 8 ) bi1cs(17)
  integer ( kind = 4 ) nti1
  real ( kind = 8 ) r8_besi1
  real ( kind = 8 ) r8_besi1e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save bi1cs
  save nti1
  save xmax
  save xmin
  save xsml

  data bi1cs(  1) / -0.19717132610998597316138503218149D-02 /
  data bi1cs(  2) / +0.40734887667546480608155393652014D+00 /
  data bi1cs(  3) / +0.34838994299959455866245037783787D-01 /
  data bi1cs(  4) / +0.15453945563001236038598401058489D-02 /
  data bi1cs(  5) / +0.41888521098377784129458832004120D-04 /
  data bi1cs(  6) / +0.76490267648362114741959703966069D-06 /
  data bi1cs(  7) / +0.10042493924741178689179808037238D-07 /
  data bi1cs(  8) / +0.99322077919238106481371298054863D-10 /
  data bi1cs(  9) / +0.76638017918447637275200171681349D-12 /
  data bi1cs( 10) / +0.47414189238167394980388091948160D-14 /
  data bi1cs( 11) / +0.24041144040745181799863172032000D-16 /
  data bi1cs( 12) / +0.10171505007093713649121100799999D-18 /
  data bi1cs( 13) / +0.36450935657866949458491733333333D-21 /
  data bi1cs( 14) / +0.11205749502562039344810666666666D-23 /
  data bi1cs( 15) / +0.29875441934468088832000000000000D-26 /
  data bi1cs( 16) / +0.69732310939194709333333333333333D-29 /
  data bi1cs( 17) / +0.14367948220620800000000000000000D-31 /

  data nti1 / 0 /
  data xmax / 0.0D+00 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( nti1 == 0 ) then
    nti1 = r8_inits ( bi1cs, 17, 0.1D+00 * r8_mach ( 3 ) )
    xmin = 2.0D+00 * r8_mach ( 1 )
    xsml = sqrt ( 8.0D+00 * r8_mach ( 3 ) )
    xmax = log ( r8_mach ( 2 ) )
  end if

  y = abs ( x )

  if ( y <= xmin ) then
    r8_besi1 = 0.0D+00
  else if ( y <= xsml ) then
    r8_besi1 = 0.5D+00 * x
  else if ( y <= 3.0D+00 ) then
    r8_besi1 = x * ( 0.875D+00 &
      + r8_csevl ( y * y / 4.5D+00 - 1.0D+00, bi1cs, nti1 ) )
  else if ( y <= xmax ) then
    r8_besi1 = exp ( y ) * r8_besi1e ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESI1 - Fatal error!'
    write ( *, '(a)' ) '  Result overflows.'
    stop
  end if

  return
end
function r8_besi1e ( x )

!*****************************************************************************80
!
!! R8_BESI1E evaluates the exponentially scaled Bessel function I1(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESI1E, the exponentially scaled Bessel 
!    function I1(X).
!
  implicit none

  real ( kind = 8 ) ai12cs(69)
  real ( kind = 8 ) ai1cs(46)
  real ( kind = 8 ) bi1cs(17)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) ntai1
  integer ( kind = 4 ) ntai12
  integer ( kind = 4 ) nti1
  real ( kind = 8 ) r8_besi1e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save ai12cs
  save ai1cs
  save bi1cs
  save ntai1
  save ntai12
  save nti1
  save xmin
  save xsml

  data bi1cs(  1) / -0.19717132610998597316138503218149D-02 /
  data bi1cs(  2) / +0.40734887667546480608155393652014D+00 /
  data bi1cs(  3) / +0.34838994299959455866245037783787D-01 /
  data bi1cs(  4) / +0.15453945563001236038598401058489D-02 /
  data bi1cs(  5) / +0.41888521098377784129458832004120D-04 /
  data bi1cs(  6) / +0.76490267648362114741959703966069D-06 /
  data bi1cs(  7) / +0.10042493924741178689179808037238D-07 /
  data bi1cs(  8) / +0.99322077919238106481371298054863D-10 /
  data bi1cs(  9) / +0.76638017918447637275200171681349D-12 /
  data bi1cs( 10) / +0.47414189238167394980388091948160D-14 /
  data bi1cs( 11) / +0.24041144040745181799863172032000D-16 /
  data bi1cs( 12) / +0.10171505007093713649121100799999D-18 /
  data bi1cs( 13) / +0.36450935657866949458491733333333D-21 /
  data bi1cs( 14) / +0.11205749502562039344810666666666D-23 /
  data bi1cs( 15) / +0.29875441934468088832000000000000D-26 /
  data bi1cs( 16) / +0.69732310939194709333333333333333D-29 /
  data bi1cs( 17) / +0.14367948220620800000000000000000D-31 /

  data ai1cs(  1) / -0.2846744181881478674100372468307D-01 /
  data ai1cs(  2) / -0.1922953231443220651044448774979D-01 /
  data ai1cs(  3) / -0.6115185857943788982256249917785D-03 /
  data ai1cs(  4) / -0.2069971253350227708882823777979D-04 /
  data ai1cs(  5) / +0.8585619145810725565536944673138D-05 /
  data ai1cs(  6) / +0.1049498246711590862517453997860D-05 /
  data ai1cs(  7) / -0.2918338918447902202093432326697D-06 /
  data ai1cs(  8) / -0.1559378146631739000160680969077D-07 /
  data ai1cs(  9) / +0.1318012367144944705525302873909D-07 /
  data ai1cs( 10) / -0.1448423418183078317639134467815D-08 /
  data ai1cs( 11) / -0.2908512243993142094825040993010D-09 /
  data ai1cs( 12) / +0.1266388917875382387311159690403D-09 /
  data ai1cs( 13) / -0.1664947772919220670624178398580D-10 /
  data ai1cs( 14) / -0.1666653644609432976095937154999D-11 /
  data ai1cs( 15) / +0.1242602414290768265232168472017D-11 /
  data ai1cs( 16) / -0.2731549379672432397251461428633D-12 /
  data ai1cs( 17) / +0.2023947881645803780700262688981D-13 /
  data ai1cs( 18) / +0.7307950018116883636198698126123D-14 /
  data ai1cs( 19) / -0.3332905634404674943813778617133D-14 /
  data ai1cs( 20) / +0.7175346558512953743542254665670D-15 /
  data ai1cs( 21) / -0.6982530324796256355850629223656D-16 /
  data ai1cs( 22) / -0.1299944201562760760060446080587D-16 /
  data ai1cs( 23) / +0.8120942864242798892054678342860D-17 /
  data ai1cs( 24) / -0.2194016207410736898156266643783D-17 /
  data ai1cs( 25) / +0.3630516170029654848279860932334D-18 /
  data ai1cs( 26) / -0.1695139772439104166306866790399D-19 /
  data ai1cs( 27) / -0.1288184829897907807116882538222D-19 /
  data ai1cs( 28) / +0.5694428604967052780109991073109D-20 /
  data ai1cs( 29) / -0.1459597009090480056545509900287D-20 /
  data ai1cs( 30) / +0.2514546010675717314084691334485D-21 /
  data ai1cs( 31) / -0.1844758883139124818160400029013D-22 /
  data ai1cs( 32) / -0.6339760596227948641928609791999D-23 /
  data ai1cs( 33) / +0.3461441102031011111108146626560D-23 /
  data ai1cs( 34) / -0.1017062335371393547596541023573D-23 /
  data ai1cs( 35) / +0.2149877147090431445962500778666D-24 /
  data ai1cs( 36) / -0.3045252425238676401746206173866D-25 /
  data ai1cs( 37) / +0.5238082144721285982177634986666D-27 /
  data ai1cs( 38) / +0.1443583107089382446416789503999D-26 /
  data ai1cs( 39) / -0.6121302074890042733200670719999D-27 /
  data ai1cs( 40) / +0.1700011117467818418349189802666D-27 /
  data ai1cs( 41) / -0.3596589107984244158535215786666D-28 /
  data ai1cs( 42) / +0.5448178578948418576650513066666D-29 /
  data ai1cs( 43) / -0.2731831789689084989162564266666D-30 /
  data ai1cs( 44) / -0.1858905021708600715771903999999D-30 /
  data ai1cs( 45) / +0.9212682974513933441127765333333D-31 /
  data ai1cs( 46) / -0.2813835155653561106370833066666D-31 /

  data ai12cs(  1) / +0.2857623501828012047449845948469D-01  /
  data ai12cs(  2) / -0.9761097491361468407765164457302D-02  /
  data ai12cs(  3) / -0.1105889387626237162912569212775D-03  /
  data ai12cs(  4) / -0.3882564808877690393456544776274D-05  /
  data ai12cs(  5) / -0.2512236237870208925294520022121D-06  /
  data ai12cs(  6) / -0.2631468846889519506837052365232D-07  /
  data ai12cs(  7) / -0.3835380385964237022045006787968D-08  /
  data ai12cs(  8) / -0.5589743462196583806868112522229D-09  /
  data ai12cs(  9) / -0.1897495812350541234498925033238D-10 /
  data ai12cs( 10) / +0.3252603583015488238555080679949D-10 /
  data ai12cs( 11) / +0.1412580743661378133163366332846D-10 /
  data ai12cs( 12) / +0.2035628544147089507224526136840D-11 /
  data ai12cs( 13) / -0.7198551776245908512092589890446D-12 /
  data ai12cs( 14) / -0.4083551111092197318228499639691D-12 /
  data ai12cs( 15) / -0.2101541842772664313019845727462D-13 /
  data ai12cs( 16) / +0.4272440016711951354297788336997D-13 /
  data ai12cs( 17) / +0.1042027698412880276417414499948D-13 /
  data ai12cs( 18) / -0.3814403072437007804767072535396D-14 /
  data ai12cs( 19) / -0.1880354775510782448512734533963D-14 /
  data ai12cs( 20) / +0.3308202310920928282731903352405D-15 /
  data ai12cs( 21) / +0.2962628997645950139068546542052D-15 /
  data ai12cs( 22) / -0.3209525921993423958778373532887D-16 /
  data ai12cs( 23) / -0.4650305368489358325571282818979D-16 /
  data ai12cs( 24) / +0.4414348323071707949946113759641D-17 /
  data ai12cs( 25) / +0.7517296310842104805425458080295D-17 /
  data ai12cs( 26) / -0.9314178867326883375684847845157D-18 /
  data ai12cs( 27) / -0.1242193275194890956116784488697D-17 /
  data ai12cs( 28) / +0.2414276719454848469005153902176D-18 /
  data ai12cs( 29) / +0.2026944384053285178971922860692D-18 /
  data ai12cs( 30) / -0.6394267188269097787043919886811D-19 /
  data ai12cs( 31) / -0.3049812452373095896084884503571D-19 /
  data ai12cs( 32) / +0.1612841851651480225134622307691D-19 /
  data ai12cs( 33) / +0.3560913964309925054510270904620D-20 /
  data ai12cs( 34) / -0.3752017947936439079666828003246D-20 /
  data ai12cs( 35) / -0.5787037427074799345951982310741D-22 /
  data ai12cs( 36) / +0.7759997511648161961982369632092D-21 /
  data ai12cs( 37) / -0.1452790897202233394064459874085D-21 /
  data ai12cs( 38) / -0.1318225286739036702121922753374D-21 /
  data ai12cs( 39) / +0.6116654862903070701879991331717D-22 /
  data ai12cs( 40) / +0.1376279762427126427730243383634D-22 /
  data ai12cs( 41) / -0.1690837689959347884919839382306D-22 /
  data ai12cs( 42) / +0.1430596088595433153987201085385D-23 /
  data ai12cs( 43) / +0.3409557828090594020405367729902D-23 /
  data ai12cs( 44) / -0.1309457666270760227845738726424D-23 /
  data ai12cs( 45) / -0.3940706411240257436093521417557D-24 /
  data ai12cs( 46) / +0.4277137426980876580806166797352D-24 /
  data ai12cs( 47) / -0.4424634830982606881900283123029D-25 /
  data ai12cs( 48) / -0.8734113196230714972115309788747D-25 /
  data ai12cs( 49) / +0.4045401335683533392143404142428D-25 /
  data ai12cs( 50) / +0.7067100658094689465651607717806D-26 /
  data ai12cs( 51) / -0.1249463344565105223002864518605D-25 /
  data ai12cs( 52) / +0.2867392244403437032979483391426D-26 /
  data ai12cs( 53) / +0.2044292892504292670281779574210D-26 /
  data ai12cs( 54) / -0.1518636633820462568371346802911D-26 /
  data ai12cs( 55) / +0.8110181098187575886132279107037D-28 /
  data ai12cs( 56) / +0.3580379354773586091127173703270D-27 /
  data ai12cs( 57) / -0.1692929018927902509593057175448D-27 /
  data ai12cs( 58) / -0.2222902499702427639067758527774D-28 /
  data ai12cs( 59) / +0.5424535127145969655048600401128D-28 /
  data ai12cs( 60) / -0.1787068401578018688764912993304D-28 /
  data ai12cs( 61) / -0.6565479068722814938823929437880D-29 /
  data ai12cs( 62) / +0.7807013165061145280922067706839D-29 /
  data ai12cs( 63) / -0.1816595260668979717379333152221D-29 /
  data ai12cs( 64) / -0.1287704952660084820376875598959D-29 /
  data ai12cs( 65) / +0.1114548172988164547413709273694D-29 /
  data ai12cs( 66) / -0.1808343145039336939159368876687D-30 /
  data ai12cs( 67) / -0.2231677718203771952232448228939D-30 /
  data ai12cs( 68) / +0.1619029596080341510617909803614D-30 /
  data ai12cs( 69) / -0.1834079908804941413901308439210D-31 /

  data ntai1 / 0 /
  data ntai12 / 0 /
  data nti1 / 0 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( nti1 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nti1 = r8_inits ( bi1cs, 17, eta )
    ntai1 = r8_inits ( ai1cs, 46, eta )
    ntai12 = r8_inits ( ai12cs, 69, eta )
    xmin = 2.0D+00 * r8_mach ( 1 )
    xsml = sqrt ( 8.0D+00 * r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= xmin ) then
    r8_besi1e = 0.0D+00
  else if ( y <= xsml ) then
    r8_besi1e = 0.5D+00 * x * exp ( - y )
  else if ( y <= 3.0D+00 ) then
    r8_besi1e = x * ( 0.875D+00 &
      + r8_csevl ( y * y / 4.5D+00 - 1.0D+00, bi1cs, nti1 ) ) &
      * exp ( - y ) 
  else if ( y <= 8.0D+00 ) then
    r8_besi1e = ( 0.375D+00 &
      + r8_csevl ( ( 48.0D+00 / y - 11.0D+00) / 5.0D+00, &
      ai1cs, ntai1 ) ) / sqrt ( y )
    if ( x < 0.0D+00 ) then
      r8_besi1e = - r8_besi1e
    end if
  else
    r8_besi1e = ( 0.375D+00 &
      + r8_csevl ( 16.0D+00 / y - 1.0D+00, ai12cs, ntai12 ) ) &
      / sqrt ( y )
    if ( x < 0.0D+00 ) then
      r8_besi1e = - r8_besi1e
    end if
  end if

  return
end
function r8_besj0 ( x )

!*****************************************************************************80
!
!! R8_BESJ0 evaluates the Bessel function J of order 0 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESJ0, the Bessel function J of order 0 of X.
!
  implicit none

  real ( kind = 8 ) ampl
  real ( kind = 8 ) bj0cs(19)
  integer ( kind = 4 ) ntj0
  real ( kind = 8 ) r8_besj0
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save bj0cs
  save ntj0
  save xsml

  data bj0cs(  1) / +0.10025416196893913701073127264074D+00 /
  data bj0cs(  2) / -0.66522300776440513177678757831124D+00 /
  data bj0cs(  3) / +0.24898370349828131370460468726680D+00 /
  data bj0cs(  4) / -0.33252723170035769653884341503854D-01 /
  data bj0cs(  5) / +0.23114179304694015462904924117729D-02 /
  data bj0cs(  6) / -0.99112774199508092339048519336549D-04 /
  data bj0cs(  7) / +0.28916708643998808884733903747078D-05 /
  data bj0cs(  8) / -0.61210858663032635057818407481516D-07 /
  data bj0cs(  9) / +0.98386507938567841324768748636415D-09 /
  data bj0cs( 10) / -0.12423551597301765145515897006836D-10 /
  data bj0cs( 11) / +0.12654336302559045797915827210363D-12 /
  data bj0cs( 12) / -0.10619456495287244546914817512959D-14 /
  data bj0cs( 13) / +0.74706210758024567437098915584000D-17 /
  data bj0cs( 14) / -0.44697032274412780547627007999999D-19 /
  data bj0cs( 15) / +0.23024281584337436200523093333333D-21 /
  data bj0cs( 16) / -0.10319144794166698148522666666666D-23 /
  data bj0cs( 17) / +0.40608178274873322700800000000000D-26 /
  data bj0cs( 18) / -0.14143836005240913919999999999999D-28 /
  data bj0cs( 19) / +0.43910905496698880000000000000000D-31 /

  data ntj0 / 0 /
  data xsml / 0.0D+00 /

  if ( ntj0 == 0 ) then
    ntj0 = r8_inits ( bj0cs, 19, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= xsml ) then
    r8_besj0 = 1.0D+00
  else if ( y <= 4.0D+00 ) then
    r8_besj0 = r8_csevl ( 0.125D+00 * y * y - 1.0D+00, bj0cs, ntj0 )
  else
    call r8_b0mp ( y, ampl, theta )
    r8_besj0 = ampl * cos ( theta )
  end if

  return
end
function r8_besj1 ( x )

!*****************************************************************************80
!
!! R8_BESJ1 evaluates the Bessel function J of order 1 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESJ1, the Bessel function J of order 1 of X.
!
  implicit none

  real ( kind = 8 ) ampl
  real ( kind = 8 ) bj1cs(19)
  integer ( kind = 4 ) ntj1
  real ( kind = 8 ) r8_besj1
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save bj1cs
  save ntj1
  save xmin
  save xsml

  data bj1cs(  1) / -0.117261415133327865606240574524003D+00 /
  data bj1cs(  2) / -0.253615218307906395623030884554698D+00 /
  data bj1cs(  3) / +0.501270809844695685053656363203743D-01 /
  data bj1cs(  4) / -0.463151480962508191842619728789772D-02 /
  data bj1cs(  5) / +0.247996229415914024539124064592364D-03 /
  data bj1cs(  6) / -0.867894868627882584521246435176416D-05 /
  data bj1cs(  7) / +0.214293917143793691502766250991292D-06 /
  data bj1cs(  8) / -0.393609307918317979229322764073061D-08 /
  data bj1cs(  9) / +0.559118231794688004018248059864032D-10 /
  data bj1cs( 10) / -0.632761640466139302477695274014880D-12 /
  data bj1cs( 11) / +0.584099161085724700326945563268266D-14 /
  data bj1cs( 12) / -0.448253381870125819039135059199999D-16 /
  data bj1cs( 13) / +0.290538449262502466306018688000000D-18 /
  data bj1cs( 14) / -0.161173219784144165412118186666666D-20 /
  data bj1cs( 15) / +0.773947881939274637298346666666666D-23 /
  data bj1cs( 16) / -0.324869378211199841143466666666666D-25 /
  data bj1cs( 17) / +0.120223767722741022720000000000000D-27 /
  data bj1cs( 18) / -0.395201221265134933333333333333333D-30 /
  data bj1cs( 19) / +0.116167808226645333333333333333333D-32 /

  data ntj1 / 0 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ntj1 == 0 ) then
    ntj1 = r8_inits ( bj1cs, 19, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
    xmin = 2.0D+00 * r8_mach ( 1 )
  end if

  y = abs ( x )

  if ( y <= xmin ) then
    r8_besj1 = 0.0D+00
  else if ( y <= xsml ) then
    r8_besj1 = 0.5D+00 * x
  else if ( y <= 4.0D+00 ) then
    r8_besj1 = x * ( 0.25D+00 &
      + r8_csevl ( 0.125D+00 * y * y - 1.0D+00, bj1cs, ntj1 ) )
  else
    call r8_b1mp ( y, ampl, theta )
    if ( x < 0.0D+00 ) then
      r8_besj1 = - ampl * cos ( theta )
    else
      r8_besj1 = + ampl * cos ( theta )
    end if
  end if

  return
end
function r8_besk ( nu, x )

!*****************************************************************************80
!
!! R8_BESK evaluates the Bessel function K of order NU of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) NU, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESK, the Bessel function K of order NU at X.
!
  implicit none

  real ( kind = 8 ), allocatable :: bke(:)
  integer ( kind = 4 ) nin
  real ( kind = 8 ) nu
  real ( kind = 8 ) r8_besk
  real ( kind = 8 ) x
  real ( kind = 8 ) xnu

  xnu = nu - int ( nu )
  nin = int ( nu ) + 1
  allocate ( bke(1:nin) )

  call r8_besks ( xnu, x, nin, bke )

  r8_besk = bke(nin)

  deallocate ( bke )

  return
end
function r8_besk0 ( x )

!*****************************************************************************80
!
!! R8_BESK0 evaluates the Bessel function K of order 0 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESK0, the Bessel function K of order 0 of X.
!
  implicit none

  real ( kind = 8 ) bk0cs(16)
  integer ( kind = 4 ) ntk0
  real ( kind = 8 ) r8_besi0
  real ( kind = 8 ) r8_besk0
  real ( kind = 8 ) r8_besk0e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save bk0cs
  save ntk0
  save xmax
  save xsml

  data bk0cs(  1) / -0.353273932339027687201140060063153D-01 /
  data bk0cs(  2) / +0.344289899924628486886344927529213D+00 /
  data bk0cs(  3) / +0.359799365153615016265721303687231D-01 /
  data bk0cs(  4) / +0.126461541144692592338479508673447D-02 /
  data bk0cs(  5) / +0.228621210311945178608269830297585D-04 /
  data bk0cs(  6) / +0.253479107902614945730790013428354D-06 /
  data bk0cs(  7) / +0.190451637722020885897214059381366D-08 /
  data bk0cs(  8) / +0.103496952576336245851008317853089D-10 /
  data bk0cs(  9) / +0.425981614279108257652445327170133D-13 /
  data bk0cs( 10) / +0.137446543588075089694238325440000D-15 /
  data bk0cs( 11) / +0.357089652850837359099688597333333D-18 /
  data bk0cs( 12) / +0.763164366011643737667498666666666D-21 /
  data bk0cs( 13) / +0.136542498844078185908053333333333D-23 /
  data bk0cs( 14) / +0.207527526690666808319999999999999D-26 /
  data bk0cs( 15) / +0.271281421807298560000000000000000D-29 /
  data bk0cs( 16) / +0.308259388791466666666666666666666D-32 /

  data ntk0 / 0 /
  data xmax / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ntk0 == 0 ) then
    ntk0 = r8_inits (bk0cs, 16, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )  
    xmax = - log ( r8_mach ( 1 ) )
    xmax = xmax - 0.5D+00 * xmax * log ( xmax ) &
      / ( xmax + 0.5D+00 )
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESK0 = Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0D+00
    r8_besk0 = - log ( 0.5D+00 * x ) * r8_besi0 ( x ) &
      - 0.25D+00 + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk0cs, ntk0 )
  else if ( x <= 2.0D+00 ) then
    y = x * x
    r8_besk0 = - log ( 0.5D+00 * x ) * r8_besi0 ( x ) &
      - 0.25D+00 + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk0cs, ntk0 )
  else if ( x <= xmax ) then
    r8_besk0 = exp ( - x ) * r8_besk0e ( x )
  else
    r8_besk0 = 0.0D+00
  end if

  return
end
function r8_besk0e ( x )

!*****************************************************************************80
!
!! R8_BESK0E evaluates the exponentially scaled Bessel function K0(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESK0E, the exponentially scaled Bessel 
!    function K0(X).
!
  implicit none

  real ( kind = 8 ) ak02cs(33)
  real ( kind = 8 ) ak0cs(38)
  real ( kind = 8 ) bk0cs(16)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) ntak0
  integer ( kind = 4 ) ntak02
  integer ( kind = 4 ) ntk0
  real ( kind = 8 ) r8_besi0
  real ( kind = 8 ) r8_besk0e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save ak02cs
  save ak0cs
  save bk0cs
  save ntak0
  save ntak02
  save ntk0
  save xsml

  data bk0cs(  1) / -0.353273932339027687201140060063153D-01 /
  data bk0cs(  2) / +0.344289899924628486886344927529213D+00 /
  data bk0cs(  3) / +0.359799365153615016265721303687231D-01 /
  data bk0cs(  4) / +0.126461541144692592338479508673447D-02 /
  data bk0cs(  5) / +0.228621210311945178608269830297585D-04 /
  data bk0cs(  6) / +0.253479107902614945730790013428354D-06 /
  data bk0cs(  7) / +0.190451637722020885897214059381366D-08 /
  data bk0cs(  8) / +0.103496952576336245851008317853089D-10 /
  data bk0cs(  9) / +0.425981614279108257652445327170133D-13 /
  data bk0cs( 10) / +0.137446543588075089694238325440000D-15 /
  data bk0cs( 11) / +0.357089652850837359099688597333333D-18 /
  data bk0cs( 12) / +0.763164366011643737667498666666666D-21 /
  data bk0cs( 13) / +0.136542498844078185908053333333333D-23 /
  data bk0cs( 14) / +0.207527526690666808319999999999999D-26 /
  data bk0cs( 15) / +0.271281421807298560000000000000000D-29 /
  data bk0cs( 16) / +0.308259388791466666666666666666666D-32 /

  data ak0cs(  1) / -0.7643947903327941424082978270088D-01 /
  data ak0cs(  2) / -0.2235652605699819052023095550791D-01 /
  data ak0cs(  3) / +0.7734181154693858235300618174047D-03 /
  data ak0cs(  4) / -0.4281006688886099464452146435416D-04 /
  data ak0cs(  5) / +0.3081700173862974743650014826660D-05 /
  data ak0cs(  6) / -0.2639367222009664974067448892723D-06 /
  data ak0cs(  7) / +0.2563713036403469206294088265742D-07 /
  data ak0cs(  8) / -0.2742705549900201263857211915244D-08 /
  data ak0cs(  9) / +0.3169429658097499592080832873403D-09 /
  data ak0cs( 10) / -0.3902353286962184141601065717962D-10 /
  data ak0cs( 11) / +0.5068040698188575402050092127286D-11 /
  data ak0cs( 12) / -0.6889574741007870679541713557984D-12 /
  data ak0cs( 13) / +0.9744978497825917691388201336831D-13 /
  data ak0cs( 14) / -0.1427332841884548505389855340122D-13 /
  data ak0cs( 15) / +0.2156412571021463039558062976527D-14 /
  data ak0cs( 16) / -0.3349654255149562772188782058530D-15 /
  data ak0cs( 17) / +0.5335260216952911692145280392601D-16 /
  data ak0cs( 18) / -0.8693669980890753807639622378837D-17 /
  data ak0cs( 19) / +0.1446404347862212227887763442346D-17 /
  data ak0cs( 20) / -0.2452889825500129682404678751573D-18 /
  data ak0cs( 21) / +0.4233754526232171572821706342400D-19 /
  data ak0cs( 22) / -0.7427946526454464195695341294933D-20 /
  data ak0cs( 23) / +0.1323150529392666866277967462400D-20 /
  data ak0cs( 24) / -0.2390587164739649451335981465599D-21 /
  data ak0cs( 25) / +0.4376827585923226140165712554666D-22 /
  data ak0cs( 26) / -0.8113700607345118059339011413333D-23 /
  data ak0cs( 27) / +0.1521819913832172958310378154666D-23 /
  data ak0cs( 28) / -0.2886041941483397770235958613333D-24 /
  data ak0cs( 29) / +0.5530620667054717979992610133333D-25 /
  data ak0cs( 30) / -0.1070377329249898728591633066666D-25 /
  data ak0cs( 31) / +0.2091086893142384300296328533333D-26 /
  data ak0cs( 32) / -0.4121713723646203827410261333333D-27 /
  data ak0cs( 33) / +0.8193483971121307640135680000000D-28 /
  data ak0cs( 34) / -0.1642000275459297726780757333333D-28 /
  data ak0cs( 35) / +0.3316143281480227195890346666666D-29 /
  data ak0cs( 36) / -0.6746863644145295941085866666666D-30 /
  data ak0cs( 37) / +0.1382429146318424677635413333333D-30 /
  data ak0cs( 38) / -0.2851874167359832570811733333333D-31 /

  data ak02cs(  1) / -0.1201869826307592239839346212452D-01 /
  data ak02cs(  2) / -0.9174852691025695310652561075713D-02 /
  data ak02cs(  3) / +0.1444550931775005821048843878057D-03 /
  data ak02cs(  4) / -0.4013614175435709728671021077879D-05 /
  data ak02cs(  5) / +0.1567831810852310672590348990333D-06 /
  data ak02cs(  6) / -0.7770110438521737710315799754460D-08 /
  data ak02cs(  7) / +0.4611182576179717882533130529586D-09 /
  data ak02cs(  8) / -0.3158592997860565770526665803309D-10 /
  data ak02cs(  9) / +0.2435018039365041127835887814329D-11 /
  data ak02cs( 10) / -0.2074331387398347897709853373506D-12 /
  data ak02cs( 11) / +0.1925787280589917084742736504693D-13 /
  data ak02cs( 12) / -0.1927554805838956103600347182218D-14 /
  data ak02cs( 13) / +0.2062198029197818278285237869644D-15 /
  data ak02cs( 14) / -0.2341685117579242402603640195071D-16 /
  data ak02cs( 15) / +0.2805902810643042246815178828458D-17 /
  data ak02cs( 16) / -0.3530507631161807945815482463573D-18 /
  data ak02cs( 17) / +0.4645295422935108267424216337066D-19 /
  data ak02cs( 18) / -0.6368625941344266473922053461333D-20 /
  data ak02cs( 19) / +0.9069521310986515567622348800000D-21 /
  data ak02cs( 20) / -0.1337974785423690739845005311999D-21 /
  data ak02cs( 21) / +0.2039836021859952315522088960000D-22 /
  data ak02cs( 22) / -0.3207027481367840500060869973333D-23 /
  data ak02cs( 23) / +0.5189744413662309963626359466666D-24 /
  data ak02cs( 24) / -0.8629501497540572192964607999999D-25 /
  data ak02cs( 25) / +0.1472161183102559855208038400000D-25 /
  data ak02cs( 26) / -0.2573069023867011283812351999999D-26 /
  data ak02cs( 27) / +0.4601774086643516587376640000000D-27 /
  data ak02cs( 28) / -0.8411555324201093737130666666666D-28 /
  data ak02cs( 29) / +0.1569806306635368939301546666666D-28 /
  data ak02cs( 30) / -0.2988226453005757788979199999999D-29 /
  data ak02cs( 31) / +0.5796831375216836520618666666666D-30 /
  data ak02cs( 32) / -0.1145035994347681332155733333333D-30 /
  data ak02cs( 33) / +0.2301266594249682802005333333333D-31 /

  data ntak0 / 0 /
  data ntak02 / 0 /
  data ntk0 / 0 /
  data xsml / 0.0D+00 /

  if ( ntk0 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    ntk0 = r8_inits ( bk0cs, 16, eta )
    ntak0 = r8_inits ( ak0cs, 38, eta )
    ntak02 = r8_inits ( ak02cs, 33, eta )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESK0E = Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0D+00
    r8_besk0e = exp ( x ) * ( - log ( 0.5D+00 * x ) &
      * r8_besi0 ( x ) - 0.25D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk0cs, ntk0 ) ) 
  else if ( x <= 2.0D+00 ) then
    y = x * x
    r8_besk0e = exp ( x ) * ( - log ( 0.5D+00 * x ) &
      * r8_besi0 ( x ) - 0.25D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk0cs, ntk0 ) )
  else if ( x <= 8.0D+00 ) then
    r8_besk0e = ( 1.25D+00 &
      + r8_csevl ( ( 16.0D+00 / x - 5.0D+00 ) / 3.0D+00, ak0cs, &
      ntak0 ) ) / sqrt ( x )
  else
    r8_besk0e = ( 1.25D+00 + &
      r8_csevl ( 16.0D+00 / x - 1.0D+00, ak02cs, ntak02 ) ) &
      / sqrt ( x )
  end if
 
  return
end
function r8_besk1 ( x )

!*****************************************************************************80
!
!! R8_BESK1 evaluates the Bessel function K of order 1 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESK1, the Bessel function K of order 1 of X.
!
  implicit none

  real ( kind = 8 ) bk1cs(16)
  integer ( kind = 4 ) ntk1
  real ( kind = 8 ) r8_besi1
  real ( kind = 8 ) r8_besk1
  real ( kind = 8 ) r8_besk1e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save bk1cs
  save ntk1
  save xmax
  save xmin
  save xsml

  data bk1cs(  1) / +0.25300227338947770532531120868533D-01 /
  data bk1cs(  2) / -0.35315596077654487566723831691801D+00 /
  data bk1cs(  3) / -0.12261118082265714823479067930042D+00 /
  data bk1cs(  4) / -0.69757238596398643501812920296083D-02 /
  data bk1cs(  5) / -0.17302889575130520630176507368979D-03 /
  data bk1cs(  6) / -0.24334061415659682349600735030164D-05 /
  data bk1cs(  7) / -0.22133876307347258558315252545126D-07 /
  data bk1cs(  8) / -0.14114883926335277610958330212608D-09 /
  data bk1cs(  9) / -0.66669016941993290060853751264373D-12 /
  data bk1cs( 10) / -0.24274498505193659339263196864853D-14 /
  data bk1cs( 11) / -0.70238634793862875971783797120000D-17 /
  data bk1cs( 12) / -0.16543275155100994675491029333333D-19 /
  data bk1cs( 13) / -0.32338347459944491991893333333333D-22 /
  data bk1cs( 14) / -0.53312750529265274999466666666666D-25 /
  data bk1cs( 15) / -0.75130407162157226666666666666666D-28 /
  data bk1cs( 16) / -0.91550857176541866666666666666666D-31 /

  data ntk1 / 0 /
  data xmax / 0.0D+00 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ntk1 == 0 ) then
    ntk1 = r8_inits ( bk1cs, 16, 0.1D+00 * r8_mach ( 3 ) )
    xmin = exp ( max ( log ( r8_mach ( 1 ) ), &
      - log ( r8_mach ( 2 ) ) ) + 0.01D+00 )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
    xmax = - log ( r8_mach ( 1 ) )
    xmax = xmax - 0.5D+00 * xmax * log ( xmax ) &
      / ( xmax + 0.5D+00 ) - 0.01D+00
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESK1 = Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0D+00
    r8_besk1 = log ( 0.5D+00 * x ) * r8_besi1 ( x ) + ( 0.75D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x
  else if ( x <= 2.0D+00 ) then
    y = x * x
    r8_besk1 = log ( 0.5D+00 * x ) * r8_besi1 ( x ) + ( 0.75D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x
  else if ( x <= xmax ) then
    r8_besk1 = exp ( - x ) * r8_besk1e ( x )
  else
    r8_besk1 = 0.0D+00
  end if

  return
end
function r8_besk1e ( x )

!*****************************************************************************80
!
!! R8_BESK1E evaluates the exponentially scaled Bessel function K1(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESK1E, the exponentially scaled Bessel 
!    function K1(X).
!
  implicit none

  real ( kind = 8 ) ak12cs(33)
  real ( kind = 8 ) ak1cs(38)
  real ( kind = 8 ) bk1cs(16)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) ntak1
  integer ( kind = 4 ) ntak12
  integer ( kind = 4 ) ntk1
  real ( kind = 8 ) r8_besi1
  real ( kind = 8 ) r8_besk1e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save ak12cs
  save ak1cs
  save bk1cs
  save ntak1
  save ntak12
  save ntk1
  save xmin
  save xsml

  data bk1cs(  1) / +0.25300227338947770532531120868533D-01 /
  data bk1cs(  2) / -0.35315596077654487566723831691801D+00 /
  data bk1cs(  3) / -0.12261118082265714823479067930042D+00 /
  data bk1cs(  4) / -0.69757238596398643501812920296083D-02 /
  data bk1cs(  5) / -0.17302889575130520630176507368979D-03 /
  data bk1cs(  6) / -0.24334061415659682349600735030164D-05 /
  data bk1cs(  7) / -0.22133876307347258558315252545126D-07 /
  data bk1cs(  8) / -0.14114883926335277610958330212608D-09 /
  data bk1cs(  9) / -0.66669016941993290060853751264373D-12 /
  data bk1cs( 10) / -0.24274498505193659339263196864853D-14 /
  data bk1cs( 11) / -0.70238634793862875971783797120000D-17 /
  data bk1cs( 12) / -0.16543275155100994675491029333333D-19 /
  data bk1cs( 13) / -0.32338347459944491991893333333333D-22 /
  data bk1cs( 14) / -0.53312750529265274999466666666666D-25 /
  data bk1cs( 15) / -0.75130407162157226666666666666666D-28 /
  data bk1cs( 16) / -0.91550857176541866666666666666666D-31 /

  data ak1cs(  1) / +0.27443134069738829695257666227266D+00 /
  data ak1cs(  2) / +0.75719899531993678170892378149290D-01 /
  data ak1cs(  3) / -0.14410515564754061229853116175625D-02 /
  data ak1cs(  4) / +0.66501169551257479394251385477036D-04 /
  data ak1cs(  5) / -0.43699847095201407660580845089167D-05 /
  data ak1cs(  6) / +0.35402774997630526799417139008534D-06 /
  data ak1cs(  7) / -0.33111637792932920208982688245704D-07 /
  data ak1cs(  8) / +0.34459775819010534532311499770992D-08 /
  data ak1cs(  9) / -0.38989323474754271048981937492758D-09 /
  data ak1cs( 10) / +0.47208197504658356400947449339005D-10 /
  data ak1cs( 11) / -0.60478356628753562345373591562890D-11 /
  data ak1cs( 12) / +0.81284948748658747888193837985663D-12 /
  data ak1cs( 13) / -0.11386945747147891428923915951042D-12 /
  data ak1cs( 14) / +0.16540358408462282325972948205090D-13 /
  data ak1cs( 15) / -0.24809025677068848221516010440533D-14 /
  data ak1cs( 16) / +0.38292378907024096948429227299157D-15 /
  data ak1cs( 17) / -0.60647341040012418187768210377386D-16 /
  data ak1cs( 18) / +0.98324256232648616038194004650666D-17 /
  data ak1cs( 19) / -0.16284168738284380035666620115626D-17 /
  data ak1cs( 20) / +0.27501536496752623718284120337066D-18 /
  data ak1cs( 21) / -0.47289666463953250924281069568000D-19 /
  data ak1cs( 22) / +0.82681500028109932722392050346666D-20 /
  data ak1cs( 23) / -0.14681405136624956337193964885333D-20 /
  data ak1cs( 24) / +0.26447639269208245978085894826666D-21 /
  data ak1cs( 25) / -0.48290157564856387897969868800000D-22 /
  data ak1cs( 26) / +0.89293020743610130180656332799999D-23 /
  data ak1cs( 27) / -0.16708397168972517176997751466666D-23 /
  data ak1cs( 28) / +0.31616456034040694931368618666666D-24 /
  data ak1cs( 29) / -0.60462055312274989106506410666666D-25 /
  data ak1cs( 30) / +0.11678798942042732700718421333333D-25 /
  data ak1cs( 31) / -0.22773741582653996232867840000000D-26 /
  data ak1cs( 32) / +0.44811097300773675795305813333333D-27 /
  data ak1cs( 33) / -0.88932884769020194062336000000000D-28 /
  data ak1cs( 34) / +0.17794680018850275131392000000000D-28 /
  data ak1cs( 35) / -0.35884555967329095821994666666666D-29 /
  data ak1cs( 36) / +0.72906290492694257991679999999999D-30 /
  data ak1cs( 37) / -0.14918449845546227073024000000000D-30 /
  data ak1cs( 38) / +0.30736573872934276300799999999999D-31 /

  data ak12cs(  1) / +0.6379308343739001036600488534102D-01 /
  data ak12cs(  2) / +0.2832887813049720935835030284708D-01 /
  data ak12cs(  3) / -0.2475370673905250345414545566732D-03 /
  data ak12cs(  4) / +0.5771972451607248820470976625763D-05 /
  data ak12cs(  5) / -0.2068939219536548302745533196552D-06 /
  data ak12cs(  6) / +0.9739983441381804180309213097887D-08 /
  data ak12cs(  7) / -0.5585336140380624984688895511129D-09 /
  data ak12cs(  8) / +0.3732996634046185240221212854731D-10 /
  data ak12cs(  9) / -0.2825051961023225445135065754928D-11 /
  data ak12cs( 10) / +0.2372019002484144173643496955486D-12 /
  data ak12cs( 11) / -0.2176677387991753979268301667938D-13 /
  data ak12cs( 12) / +0.2157914161616032453939562689706D-14 /
  data ak12cs( 13) / -0.2290196930718269275991551338154D-15 /
  data ak12cs( 14) / +0.2582885729823274961919939565226D-16 /
  data ak12cs( 15) / -0.3076752641268463187621098173440D-17 /
  data ak12cs( 16) / +0.3851487721280491597094896844799D-18 /
  data ak12cs( 17) / -0.5044794897641528977117282508800D-19 /
  data ak12cs( 18) / +0.6888673850418544237018292223999D-20 /
  data ak12cs( 19) / -0.9775041541950118303002132480000D-21 /
  data ak12cs( 20) / +0.1437416218523836461001659733333D-21 /
  data ak12cs( 21) / -0.2185059497344347373499733333333D-22 /
  data ak12cs( 22) / +0.3426245621809220631645388800000D-23 /
  data ak12cs( 23) / -0.5531064394246408232501248000000D-24 /
  data ak12cs( 24) / +0.9176601505685995403782826666666D-25 /
  data ak12cs( 25) / -0.1562287203618024911448746666666D-25 /
  data ak12cs( 26) / +0.2725419375484333132349439999999D-26 /
  data ak12cs( 27) / -0.4865674910074827992378026666666D-27 /
  data ak12cs( 28) / +0.8879388552723502587357866666666D-28 /
  data ak12cs( 29) / -0.1654585918039257548936533333333D-28 /
  data ak12cs( 30) / +0.3145111321357848674303999999999D-29 /
  data ak12cs( 31) / -0.6092998312193127612416000000000D-30 /
  data ak12cs( 32) / +0.1202021939369815834623999999999D-30 /
  data ak12cs( 33) / -0.2412930801459408841386666666666D-31 /

  data ntak1 / 0 /
  data ntak12 / 0 /
  data ntk1 / 0 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ntk1 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    ntk1 = r8_inits ( bk1cs, 16, eta )
    ntak1 = r8_inits ( ak1cs, 38, eta )
    ntak12 = r8_inits ( ak12cs, 33, eta )
    xmin = exp ( max ( log ( r8_mach ( 1 ) ), &
      - log ( r8_mach ( 2 ) ) ) + 0.01D+00 )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESK1E = Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0D+00
    r8_besk1e = exp ( x ) * ( log ( 0.5D+00 * x ) * r8_besi1 ( x ) &
      + ( 0.75D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x )
  else if ( x <= 2.0D+00 ) then
    y = x * x
    r8_besk1e = exp ( x ) * ( log ( 0.5D+00 * x ) * r8_besi1 ( x ) &
      + ( 0.75D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x )
  else if ( x <= 8.0D+00 ) then
    r8_besk1e = ( 1.25D+00 &
      + r8_csevl ( ( 16.0D+00 / x - 5.0D+00 ) / 3.0D+00, ak1cs, &
      ntak1 ) ) / sqrt ( x )
  else
    r8_besk1e = ( 1.25D+00 + &
      r8_csevl ( 16.0D+00 / x - 1.0D+00, ak12cs, ntak12 ) ) &
      / sqrt ( x )
  end if

  return
end
subroutine r8_beskes ( xnu, x, nin, bke )

!*****************************************************************************80
!
!! R8_BESKES: a sequence of exponentially scaled K Bessel functions at X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XNU, ?
!    |XNU| < 1.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) NIN, indicates the number of terms to compute.
!
!    Output, real ( kind = 8 ) BKE(abs(NIN)), the exponentially scaled 
!    K Bessel functions.
!
  implicit none

  real ( kind = 8 ) bke(*)
  real ( kind = 8 ) bknu1
  real ( kind = 8 ) direct
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iswtch
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nin
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) v
  real ( kind = 8 ) vend
  real ( kind = 8 ) vincr
  real ( kind = 8 ) x
  real ( kind = 8 ) xnu

  v = abs ( xnu )
  n = abs ( nin )

  if ( 1.0D+00 <= v ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESKES - Fatal error!'
    write ( *, '(a)' ) '  |XNU| must be less than 1.'
    stop
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESKES - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  if ( n == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESKES - Fatal error!'
    write ( *, '(a)' ) '  N = 0.'
    stop
  end if

  call r8_knus ( v, x, bke(1), bknu1, iswtch )

  if ( n == 1 ) then
    return
  end if

  if ( nin < 0 ) then
    vincr = - 1.0D+00
  else
    vincr = + 1.0D+00
  end if

  if ( xnu < 0.0D+00 ) then
    direct = - vincr
  else
    direct = vincr
  end if

  bke(2) = bknu1

  if ( direct < 0.0D+00 ) then
    call r8_knus ( abs ( xnu + vincr ), x, bke(2), bknu1, iswtch )
  end if

  if ( n == 2 ) then
    return
  end if

  vend = abs ( xnu + real ( nin, kind = 8 ) ) - 1.0D+00

  v = xnu
  do i = 3, n
    v = v + vincr
    bke(i) = 2.0D+00 * v * bke(i-1) / x + bke(i-2)
  end do

  return
end
subroutine r8_besks ( xnu, x, nin, bk )

!*****************************************************************************80
!
!! R8_BESKS evaluates a sequence of K Bessel functions at X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XNU, ?
!    |XNU| < 1.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) NIN, indicates the number of terms to compute.
!
!    Output, real ( kind = 8 ) BK(abs(NIN)), the K Bessel functions.
!
  implicit none

  integer ( kind = 4 ) nin

  real ( kind = 8 ) bk(abs(nin))
  real ( kind = 8 ) expxi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xnu

  save xmax

  data xmax / 0.0D+00 /

  if ( xmax == 0.0D+00 ) then
    xmax = - log ( r8_mach ( 1 ) )
    xmax = xmax + 0.5D+00 * log ( 3.14D+00 * 0.5D+00 / xmax )
  end if

  call r8_beskes ( xnu, x, nin, bk )

  expxi = exp ( - x )
  n = abs ( nin )

  do i = 1, n
    bk(i) = expxi * bk(i)
  end do

  return
end
function r8_besy0 ( x )

!*****************************************************************************80
!
!! R8_BESY0 evaluates the Bessel function Y of order 0 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESY0, the Bessel function Y of order 0 of X.
!
  implicit none

  real ( kind = 8 ) alnhaf
  real ( kind = 8 ) ampl
  real ( kind = 8 ) by0cs(19)
  integer ( kind = 4 ) nty0
  real ( kind = 8 ) r8_besj0
  real ( kind = 8 ) r8_besy0
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) twodpi
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save alnhaf
  save by0cs
  save nty0
  save twodpi
  save xsml

  data by0cs(  1) / -0.1127783939286557321793980546028D-01 /
  data by0cs(  2) / -0.1283452375604203460480884531838D+00 /
  data by0cs(  3) / -0.1043788479979424936581762276618D+00 /
  data by0cs(  4) / +0.2366274918396969540924159264613D-01 /
  data by0cs(  5) / -0.2090391647700486239196223950342D-02 /
  data by0cs(  6) / +0.1039754539390572520999246576381D-03 /
  data by0cs(  7) / -0.3369747162423972096718775345037D-05 /
  data by0cs(  8) / +0.7729384267670667158521367216371D-07 /
  data by0cs(  9) / -0.1324976772664259591443476068964D-08 /
  data by0cs( 10) / +0.1764823261540452792100389363158D-10 /
  data by0cs( 11) / -0.1881055071580196200602823012069D-12 /
  data by0cs( 12) / +0.1641865485366149502792237185749D-14 /
  data by0cs( 13) / -0.1195659438604606085745991006720D-16 /
  data by0cs( 14) / +0.7377296297440185842494112426666D-19 /
  data by0cs( 15) / -0.3906843476710437330740906666666D-21 /
  data by0cs( 16) / +0.1795503664436157949829120000000D-23 /
  data by0cs( 17) / -0.7229627125448010478933333333333D-26 /
  data by0cs( 18) / +0.2571727931635168597333333333333D-28 /
  data by0cs( 19) / -0.8141268814163694933333333333333D-31 /

  data alnhaf /-0.69314718055994530941723212145818D+00 /
  data nty0 / 0 /
  data twodpi / 0.636619772367581343075535053490057D+00 /
  data xsml / 0.0D+00 /

  if ( nty0 == 0 ) then
    nty0 = r8_inits ( by0cs, 19, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESY0 - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0D+00
    r8_besy0 = twodpi * ( alnhaf + log ( x ) ) * r8_besj0 ( x ) &
      + 0.375D+00 &
      + r8_csevl ( 0.125D+00 * y - 1.0D+00, by0cs, nty0 )
  else if ( x <= 4.0D+00 ) then
    y = x * x
    r8_besy0 = twodpi * ( alnhaf + log ( x ) ) * r8_besj0 ( x ) &
      + 0.375D+00 &
      + r8_csevl ( 0.125D+00 * y - 1.0D+00, by0cs, nty0 )
  else
    call r8_b0mp ( x, ampl, theta )
    r8_besy0 = ampl * sin ( theta )
  end if

  return
end
function r8_besy1 ( x )

!*****************************************************************************80
!
!! R8_BESY1 evaluates the Bessel function Y of order 1 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESY1, the Bessel function Y of order 1 of X.
!
  implicit none

  real ( kind = 8 ) ampl
  real ( kind = 8 ) by1cs(20)
  integer ( kind = 4 ) nty1
  real ( kind = 8 ) r8_besj1
  real ( kind = 8 ) r8_besy1
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) twodpi
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save by1cs
  save nty1
  save twodpi
  save xmin
  save xsml

  data by1cs(  1) / +0.320804710061190862932352018628015D-01 /
  data by1cs(  2) / +0.126270789743350044953431725999727D+01 /
  data by1cs(  3) / +0.649996189992317500097490637314144D-02 /
  data by1cs(  4) / -0.893616452886050411653144160009712D-01 /
  data by1cs(  5) / +0.132508812217570954512375510370043D-01 /
  data by1cs(  6) / -0.897905911964835237753039508298105D-03 /
  data by1cs(  7) / +0.364736148795830678242287368165349D-04 /
  data by1cs(  8) / -0.100137438166600055549075523845295D-05 /
  data by1cs(  9) / +0.199453965739017397031159372421243D-07 /
  data by1cs( 10) / -0.302306560180338167284799332520743D-09 /
  data by1cs( 11) / +0.360987815694781196116252914242474D-11 /
  data by1cs( 12) / -0.348748829728758242414552947409066D-13 /
  data by1cs( 13) / +0.278387897155917665813507698517333D-15 /
  data by1cs( 14) / -0.186787096861948768766825352533333D-17 /
  data by1cs( 15) / +0.106853153391168259757070336000000D-19 /
  data by1cs( 16) / -0.527472195668448228943872000000000D-22 /
  data by1cs( 17) / +0.227019940315566414370133333333333D-24 /
  data by1cs( 18) / -0.859539035394523108693333333333333D-27 /
  data by1cs( 19) / +0.288540437983379456000000000000000D-29 /
  data by1cs( 20) / -0.864754113893717333333333333333333D-32 /

  data nty1 / 0 /
  data twodpi / 0.636619772367581343075535053490057D+00 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( nty1 == 0 ) then
    nty1 = r8_inits ( by1cs, 20, 0.1D+00 * r8_mach ( 3 ) )
    xmin = 1.571D+00 * exp ( max ( log ( r8_mach ( 1 ) ), &
      - log ( r8_mach ( 2 ) ) ) + 0.01D+00 )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESY1 - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xmin ) then
    y = 0.0D+00
    r8_besy1 = twodpi * log ( 0.5D+00 * x ) * r8_besj1 ( x ) &
      + ( 0.5D+00 &
      + r8_csevl ( 0.125D+00 * y - 1.0D+00, by1cs, nty1 ) ) / x
  else if ( x <= 4.0D+00 ) then
    y = x * x
    r8_besy1 = twodpi * log ( 0.5D+00 * x ) * r8_besj1 ( x ) &
      + ( 0.5D+00 &
      + r8_csevl ( 0.125D+00 * y - 1.0D+00, by1cs, nty1 ) ) / x
  else
    call r8_b1mp ( x, ampl, theta )
    r8_besy1 = ampl * sin ( theta )
  end if

  return
end
function r8_beta ( a, b )

!*****************************************************************************80
!
!! R8_BETA evaluates the beta function of R8 arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the arguments.
!
!    Output, real ( kind = 8 ) R8_BETA, the beta function of A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alnsml
  real ( kind = 8 ) b
  real ( kind = 8 ) r8_beta
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_lbeta
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  save alnsml
  save xmax

  data alnsml / 0.0D+00 /
  data xmax / 0.0D+00 /

  if ( xmax == 0.0D+00 ) then
    call r8_gaml ( xmin, xmax )
    alnsml = log ( r8_mach ( 1 ) )
  end if

  if ( a <= 0.0D+00 .or. b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BETA - Fatal error!'
    write ( *, '(a)' ) '  A and B must be greater than 0.'
    stop
  end if

  if ( a + b < xmax ) then
    r8_beta = r8_gamma ( a ) * r8_gamma ( b ) / r8_gamma ( a + b )
    return
  end if

  r8_beta = r8_lbeta ( a, b )

  r8_beta = exp ( r8_beta )

  return
end
function r8_betai ( x, pin, qin )

!*****************************************************************************80
!
!! R8_BETAI evaluates the incomplete beta ratio of R8 arguments.
!
!  Discussion:
!
!    The incomplete Beta function ratio is the probability that a
!    random variable from a beta distribution having parameters
!    P and Q will be less than or equal to X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Nancy Bosten, EL Battiste,
!    Remark on Algorithm 179: 
!    Incomplete Beta Ratio,
!    Communications of the ACM,
!    Volume 17, Number 3, March 1974, pages 156-157.
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the upper limit of integration.
!    0.0 <= X <= 1.0.
!
!    Input, real ( kind = 8 ) PIN, the first distribution parameter.
!    0.0 < PIN.
!
!    Input, real ( kind = 8 ) QIN, the second distribution parameter.
!    0.0 < QIN.
!
!    Output, real ( kind = 8 ) R8_BETAI, the incomplete beta function ratio.
!
  implicit none

  real ( kind = 8 ) alneps
  real ( kind = 8 ) alnsml
  real ( kind = 8 ) c
  real ( kind = 8 ) eps
  real ( kind = 8 ) finsum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ) p1
  real ( kind = 8 ) pin
  real ( kind = 8 ) ps
  real ( kind = 8 ) q
  real ( kind = 8 ) qin
  real ( kind = 8 ) r8_betai
  real ( kind = 8 ) r8_lbeta
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sml
  real ( kind = 8 ) term
  real ( kind = 8 ) x
  real ( kind = 8 ) xb
  real ( kind = 8 ) xi
  real ( kind = 8 ) y

  save alneps
  save alnsml
  save eps
  save sml

  data alneps / 0.0D+00 /
  data alnsml / 0.0D+00 /
  data eps / 0.0D+00 /
  data sml / 0.0D+00 /

  if ( eps == 0.0D+00 ) then
    eps = r8_mach ( 3 )
    alneps = log ( eps )
    sml = r8_mach ( 1 )
    alnsml = log ( sml )
  end if

  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BETAI - Fatal error!'
    write ( *, '(a)' ) '  0 <= X <= 1 is required.'
    stop
  end if

  if ( pin <= 0.0D+00 .or. qin <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BETAI - Fatal error!'
    write ( *, '(a)' ) '  P or Q <= 0.0.'
    stop
  end if

  y = x
  p = pin
  q = qin

  if ( p < q .or. 0.8D+00 <= x ) then

    if ( 0.2D+00 <= x ) then
      y = 1.0D+00 - y
      p = qin
      q = pin
    end if

  end if

  if ( ( p + q ) * y / ( p + 1.0D+00 ) < eps ) then
    r8_betai = 0.0D+00
    xb = p * log ( max ( y, sml ) ) - log ( p ) - r8_lbeta ( p, q )
    if ( alnsml < xb .and. y /= 0.0D+00 ) then
      r8_betai = exp ( xb )
    end if
    if ( y /= x .or. p /= pin ) then
      r8_betai = 1.0D+00 - r8_betai
    end if
    return
  end if

  ps = q - aint ( q )
  if ( ps == 0.0D+00 ) then
    ps = 1.0D+00
  end if

  xb = p * log ( y ) - r8_lbeta ( ps, p ) - log ( p )

  if ( xb < alnsml ) then

    r8_betai = 0.0D+00

  else

    r8_betai = exp ( xb )
    term = r8_betai * p
    if ( ps /= 1.0D+00 ) then
      n = int ( max ( alneps / log ( y ), 4.0D+00 ) )
      do i = 1, n
        xi = real ( i, kind = 8 )
        term = term * ( xi - ps ) * y / xi
        r8_betai = r8_betai + term / ( p + xi )
      end do

    end if

  end if

  if ( 1.0D+00 < q ) then

    xb = p * log ( y ) + q * log ( 1.0D+00 - y ) &
      - r8_lbeta ( p, q ) - log ( q )
    ib = int ( max ( xb / alnsml, 0.0D+00 ) )
    term = exp ( xb - real ( ib, kind = 8 ) * alnsml )
    c = 1.0D+00 / ( 1.0D+00 - y )
    p1 = q * c / ( p + q - 1.0D+00 )

    finsum = 0.0D+00
    n = int ( q )
    if ( q == real ( n, kind = 8 ) ) then
      n = n - 1
    end if

    do i = 1, n

      if ( p1 <= 1.0D+00 .and. term / eps <= finsum ) then
        exit
      end if

      xi = real ( i, kind = 8 )
      term = ( q - xi + 1.0D+00 ) * c * term / ( p + q - xi )

      if ( 1.0D+00 < term ) then
        ib = ib - 1
        term = term * sml
      end if

      if ( ib == 0 ) then
        finsum = finsum + term
      end if

    end do

    r8_betai = r8_betai + finsum

  end if

  if ( y /= x .or. p /= pin ) then
    r8_betai = 1.0D+00 - r8_betai
  end if

  if ( r8_betai < 0.0D+00 ) then
    r8_betai =  0.0D+00
  end if

  if ( 1.0D+00 < r8_betai ) then
    r8_betai = 1.0D+00
  end if

  return
end
function r8_bi ( x )

!*****************************************************************************80
!
!! R8_BI evaluates the Airy function Bi of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BI, the Airy function Bi of X.
!
  implicit none

  real ( kind = 8 ) bifcs(13)
  real ( kind = 8 ) bif2cs(15)
  real ( kind = 8 ) bigcs(13)
  real ( kind = 8 ) big2cs(15)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nbif
  integer ( kind = 4 ) nbif2
  integer ( kind = 4 ) nbig
  integer ( kind = 4 ) nbig2
  real ( kind = 8 ) r8_bi
  real ( kind = 8 ) r8_bie
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) x3sml
  real ( kind = 8 ) xm
  real ( kind = 8 ) xmax
  real ( kind = 8 ) z

  save bifcs
  save bif2cs
  save bigcs
  save big2cs
  save nbif
  save nbif2
  save nbig
  save nbig2
  save x3sml
  save xmax

  data bifcs(  1) / -0.16730216471986649483537423928176D-01 /
  data bifcs(  2) / +0.10252335834249445611426362777757D+00 /
  data bifcs(  3) / +0.17083092507381516539429650242013D-02 /
  data bifcs(  4) / +0.11862545467744681179216459210040D-04 /
  data bifcs(  5) / +0.44932907017792133694531887927242D-07 /
  data bifcs(  6) / +0.10698207143387889067567767663628D-09 /
  data bifcs(  7) / +0.17480643399771824706010517628573D-12 /
  data bifcs(  8) / +0.20810231071761711025881891834399D-15 /
  data bifcs(  9) / +0.18849814695665416509927971733333D-18 /
  data bifcs( 10) / +0.13425779173097804625882666666666D-21 /
  data bifcs( 11) / +0.77159593429658887893333333333333D-25 /
  data bifcs( 12) / +0.36533879617478566399999999999999D-28 /
  data bifcs( 13) / +0.14497565927953066666666666666666D-31 /

  data bigcs(  1) / +0.22466223248574522283468220139024D-01 /
  data bigcs(  2) / +0.37364775453019545441727561666752D-01 /
  data bigcs(  3) / +0.44476218957212285696215294326639D-03 /
  data bigcs(  4) / +0.24708075636329384245494591948882D-05 /
  data bigcs(  5) / +0.79191353395149635134862426285596D-08 /
  data bigcs(  6) / +0.16498079851827779880887872402706D-10 /
  data bigcs(  7) / +0.24119906664835455909247501122841D-13 /
  data bigcs(  8) / +0.26103736236091436985184781269333D-16 /
  data bigcs(  9) / +0.21753082977160323853123792000000D-19 /
  data bigcs( 10) / +0.14386946400390433219483733333333D-22 /
  data bigcs( 11) / +0.77349125612083468629333333333333D-26 /
  data bigcs( 12) / +0.34469292033849002666666666666666D-29 /
  data bigcs( 13) / +0.12938919273216000000000000000000D-32 /

  data bif2cs(  1) / +0.0998457269381604104468284257993D+00 /
  data bif2cs(  2) / +0.47862497786300553772211467318231D+00 /
  data bif2cs(  3) / +0.25155211960433011771324415436675D-01 /
  data bif2cs(  4) / +0.58206938852326456396515697872216D-03 /
  data bif2cs(  5) / +0.74997659644377865943861457378217D-05 /
  data bif2cs(  6) / +0.61346028703493836681403010356474D-07 /
  data bif2cs(  7) / +0.34627538851480632900434268733359D-09 /
  data bif2cs(  8) / +0.14288910080270254287770846748931D-11 /
  data bif2cs(  9) / +0.44962704298334641895056472179200D-14 /
  data bif2cs( 10) / +0.11142323065833011708428300106666D-16 /
  data bif2cs( 11) / +0.22304791066175002081517866666666D-19 /
  data bif2cs( 12) / +0.36815778736393142842922666666666D-22 /
  data bif2cs( 13) / +0.50960868449338261333333333333333D-25 /
  data bif2cs( 14) / +0.60003386926288554666666666666666D-28 /
  data bif2cs( 15) / +0.60827497446570666666666666666666D-31 /

  data big2cs(  1) / +0.033305662145514340465176188111647D+00 /
  data big2cs(  2) / +0.161309215123197067613287532084943D+00 /
  data big2cs(  3) / +0.631900730961342869121615634921173D-02 /
  data big2cs(  4) / +0.118790456816251736389780192304567D-03 /
  data big2cs(  5) / +0.130453458862002656147116485012843D-05 /
  data big2cs(  6) / +0.937412599553521729546809615508936D-08 /
  data big2cs(  7) / +0.474580188674725153788510169834595D-10 /
  data big2cs(  8) / +0.178310726509481399800065667560946D-12 /
  data big2cs(  9) / +0.516759192784958180374276356640000D-15 /
  data big2cs( 10) / +0.119004508386827125129496251733333D-17 /
  data big2cs( 11) / +0.222982880666403517277063466666666D-20 /
  data big2cs( 12) / +0.346551923027689419722666666666666D-23 /
  data big2cs( 13) / +0.453926336320504514133333333333333D-26 /
  data big2cs( 14) / +0.507884996513522346666666666666666D-29 /
  data big2cs( 15) / +0.491020674696533333333333333333333D-32 /

  data nbif / 0 /
  data nbif2 / 0 /
  data nbig / 0 /
  data nbig2 / 0 /
  data x3sml / 0.0D+00 /
  data xmax / 0.0D+00 /

  if ( nbif == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nbif = r8_inits ( bifcs, 13, eta )
    nbig = r8_inits ( bigcs, 13, eta )
    nbif2 = r8_inits ( bif2cs, 15, eta )
    nbig2 = r8_inits ( big2cs, 15, eta )
    x3sml = eta ** 0.3333D+00
    xmax = ( 1.5D+00 * log ( r8_mach ( 2 ) ) ) ** 0.6666D+00
  end if

  if ( x < - 1.0D+00 ) then
    call r8_aimp ( x, xm, theta )
    r8_bi = xm * sin ( theta )
  else if ( abs ( x ) <= x3sml ) then
    z = 0.0D+00
    r8_bi = 0.625D+00 + r8_csevl ( z, bifcs, nbif ) &
      + x * ( 0.4375D+00 + r8_csevl ( z, bigcs, nbig ) )
  else if ( x <= 1.0D+00 ) then
    z = x * x * x
    r8_bi = 0.625D+00 + r8_csevl ( z, bifcs, nbif ) &
      + x * ( 0.4375D+00 + r8_csevl ( z, bigcs, nbig ) )
  else if ( x <= 2.0D+00 ) then
    z = ( 2.0D+00 * x * x * x - 9.0D+00 ) / 7.0D+00
    r8_bi = 1.125D+00 + r8_csevl ( z, bif2cs, nbif2 ) &
      + x * ( 0.625D+00 + r8_csevl ( z, big2cs, nbig2 ) )
  else
    r8_bi = r8_bie ( x ) &
      * exp ( 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
  end if

  return
end
function r8_bid ( x )

!*****************************************************************************80
!
!! R8_BID evaluates the derivative of the Airy function Bi of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BID, the derivative of the Airy function 
!    Bi of X.
!
  implicit none

  real ( kind = 8 ) bif2cs(15)
  real ( kind = 8 ) bifcs(13)
  real ( kind = 8 ) big2cs(16)
  real ( kind = 8 ) bigcs(13)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nbif
  integer ( kind = 4 ) nbif2
  integer ( kind = 4 ) nbig
  integer ( kind = 4 ) nbig2
  real ( kind = 8 ) phi
  real ( kind = 8 ) r8_bid
  real ( kind = 8 ) r8_bide
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) x2sml
  real ( kind = 8 ) x3
  real ( kind = 8 ) x3sml
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xn
  real ( kind = 8 ) z

  save bif2cs
  save bifcs
  save big2cs
  save bigcs
  save nbif
  save nbif2
  save nbig
  save nbig2
  save x2sml
  save x3sml
  save xmax

  data bifcs(  1) /  0.115353679082857024267474446284908879D+00 /
  data bifcs(  2) /  0.020500789404919287530357789445940252D+00 /
  data bifcs(  3) /  0.000213529027890287581892679619451158D+00 /
  data bifcs(  4) /  0.000001078396061467683042209155523569D+00 /
  data bifcs(  5) /  0.000000003209470883320666783353670420D+00 /
  data bifcs(  6) /  0.000000000006293040671833540390213316D+00 /
  data bifcs(  7) /  0.000000000000008740304300063083340121D+00 /
  data bifcs(  8) /  0.000000000000000009047915683496049529D+00 /
  data bifcs(  9) /  0.000000000000000000007249923164709251D+00 /
  data bifcs( 10) /  0.000000000000000000000004629576649604D+00 /
  data bifcs( 11) /  0.000000000000000000000000002411236436D+00 /
  data bifcs( 12) /  0.000000000000000000000000000001043825D+00 /
  data bifcs( 13) /  0.000000000000000000000000000000000382D+00 /

  data bigcs(  1) / -0.0971964404164435373897790974606802D+00 /
  data bigcs(  2) /  0.1495035768431670665710843445326264D+00 /
  data bigcs(  3) /  0.0031135253871213260419419176839631D+00 /
  data bigcs(  4) /  0.0000247085705798212967777021920569D+00 /
  data bigcs(  5) /  0.0000001029496277313786081987324295D+00 /
  data bigcs(  6) /  0.0000000002639703739869432892676778D+00 /
  data bigcs(  7) /  0.0000000000004582792707803206608181D+00 /
  data bigcs(  8) /  0.0000000000000005742829740893447321D+00 /
  data bigcs(  9) /  0.0000000000000000005438275385238549D+00 /
  data bigcs( 10) /  0.0000000000000000000004028347267083D+00 /
  data bigcs( 11) /  0.0000000000000000000000002397823826D+00 /
  data bigcs( 12) /  0.0000000000000000000000000001171956D+00 /
  data bigcs( 13) /  0.0000000000000000000000000000000479D+00 /

  data bif2cs(  1) /  0.32349398760352203352119193596266015D+00 /
  data bif2cs(  2) /  0.08629787153556355913888835323811100D+00 /
  data bif2cs(  3) /  0.00299402555265539742613821050727155D+00 /
  data bif2cs(  4) /  0.00005143052836466163720464316950821D+00 /
  data bif2cs(  5) /  0.00000052584025003681146026033098613D+00 /
  data bif2cs(  6) /  0.00000000356175137395770028102730600D+00 /
  data bif2cs(  7) /  0.00000000001714686400714584830518308D+00 /
  data bif2cs(  8) /  0.00000000000006166351969232555406693D+00 /
  data bif2cs(  9) /  0.00000000000000017191082154315985806D+00 /
  data bif2cs( 10) /  0.00000000000000000038236889518803943D+00 /
  data bif2cs( 11) /  0.00000000000000000000069424173624884D+00 /
  data bif2cs( 12) /  0.00000000000000000000000104833932510D+00 /
  data bif2cs( 13) /  0.00000000000000000000000000133721972D+00 /
  data bif2cs( 14) /  0.00000000000000000000000000000145986D+00 /
  data bif2cs( 15) /  0.00000000000000000000000000000000138D+00 /

  data big2cs(  1) /  1.606299946362129457759284537862622883D+00 /
  data big2cs(  2) /  0.744908881987608865201476685194753972D+00 /
  data big2cs(  3) /  0.047013873861027737964095177635353019D+00 /
  data big2cs(  4) /  0.001228442206254823907016188785848091D+00 /
  data big2cs(  5) /  0.000017322241225662362670987355613727D+00 /
  data big2cs(  6) /  0.000000152190165236801893711508366563D+00 /
  data big2cs(  7) /  0.000000000911356024911957704145528786D+00 /
  data big2cs(  8) /  0.000000000003954791842356644201722554D+00 /
  data big2cs(  9) /  0.000000000000013001737033862320007309D+00 /
  data big2cs( 10) /  0.000000000000000033493506858269079763D+00 /
  data big2cs( 11) /  0.000000000000000000069419094403694057D+00 /
  data big2cs( 12) /  0.000000000000000000000118248256604581D+00 /
  data big2cs( 13) /  0.000000000000000000000000168462493472D+00 /
  data big2cs( 14) /  0.000000000000000000000000000203684674D+00 /
  data big2cs( 15) /  0.000000000000000000000000000000211619D+00 /
  data big2cs( 16) /  0.000000000000000000000000000000000191D+00 /

  data nbif / 0/
  data nbif2 / 0 /
  data nbig / 0 /
  data nbig2 / 0 /
  data x2sml / 0.0D+00 /
  data x3sml / 0.0D+00 /
  data xmax / 0.0D+00 /

  if ( nbif == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nbif = r8_inits ( bifcs, 13, eta )
    nbig = r8_inits ( bigcs, 13, eta )
    nbif2 = r8_inits ( bif2cs, 15, eta )
    nbig2 = r8_inits ( big2cs, 16, eta )
    x2sml = sqrt ( eta )
    x3sml = eta ** 0.3333D+00
    xmax = ( 1.5D+00 * log ( r8_mach ( 2 ) ) ) ** 0.6666D+00
  end if

  if ( x < - 1.0D+00 ) then
    call r8_admp ( x, xn, phi )
    r8_bid = xn * sin ( phi )
  else if ( abs ( x ) <= x2sml ) then
    x2 = 0.0D+00
    x3 = 0.0D+00
    r8_bid = x2 * ( r8_csevl ( x3, bifcs, nbif ) + 0.25D+00 ) &
      + r8_csevl ( x3, bigcs, nbig ) + 0.5D+00
  else if ( abs ( x ) <= x3sml ) then
    x2 = x * x
    x3 = 0.0D+00
    r8_bid = x2 * ( r8_csevl ( x3, bifcs, nbif ) + 0.25D+00 ) &
      + r8_csevl ( x3, bigcs, nbig ) + 0.5D+00
  else if ( x <= 1.0D+00 ) then
    x2 = x * x
    x3 = x * x * x
    r8_bid = x2 * ( r8_csevl ( x3, bifcs, nbif ) + 0.25D+00 ) &
      + r8_csevl ( x3, bigcs, nbig ) + 0.5D+00
  else if ( x <= 2.0D+00 ) then
    z = ( 2.0D+00 * x * x * x - 9.0D+00 ) / 7.0D+00
    r8_bid = x * x * ( r8_csevl ( z, bif2cs, nbif2 ) + 0.25D+00 ) &
      + r8_csevl ( z, big2cs, nbig2 ) + 0.5D+00
  else
    r8_bid = r8_bide ( x ) * exp ( 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
  end if

  return
end
function r8_bide ( x )

!*****************************************************************************80
!
!! R8_BIDE: exponentially scaled derivative, Airy function Bi of an R8 argument.
!
!  Discussion:
!
!    if X < 0,
!      R8_BIDE ( X ) = R8_BID ( X )
!    else
!      R8_BIDE ( X ) = R8_BID ( X ) * exp ( - 2/3 * X**(3/2) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BIDE, the exponentially scaled derivative of 
!    the Airy function Bi of X.
!
  implicit none

  real ( kind = 8 ) atr
  real ( kind = 8 ) bif2cs(15)
  real ( kind = 8 ) bifcs(13)
  real ( kind = 8 ) big2cs(16)
  real ( kind = 8 ) bigcs(13)
  real ( kind = 8 ) bip1cs(47)
  real ( kind = 8 ) bip2cs(88)
  real ( kind = 8 ) btr
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nbif
  integer ( kind = 4 ) nbif2
  integer ( kind = 4 ) nbig
  integer ( kind = 4 ) nbig2
  integer ( kind = 4 ) nbip1
  integer ( kind = 4 ) nbip2
  real ( kind = 8 ) phi
  real ( kind = 8 ) r8_bide
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sqrtx
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) x2sml
  real ( kind = 8 ) x3
  real ( kind = 8 ) x3sml
  real ( kind = 8 ) x32sml
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xn
  real ( kind = 8 ) z

  save atr
  save bif2cs
  save bifcs
  save big2cs
  save bigcs
  save bip1cs
  save bip2cs
  save btr
  save nbif
  save nbif2
  save nbig
  save nbig2
  save nbip1
  save nbip2
  save x2sml
  save x3sml
  save x32sml
  save xbig

  data bifcs(  1) /  0.115353679082857024267474446284908879D+00 /
  data bifcs(  2) /  0.020500789404919287530357789445940252D+00 /
  data bifcs(  3) /  0.000213529027890287581892679619451158D+00 /
  data bifcs(  4) /  0.000001078396061467683042209155523569D+00 /
  data bifcs(  5) /  0.000000003209470883320666783353670420D+00 /
  data bifcs(  6) /  0.000000000006293040671833540390213316D+00 /
  data bifcs(  7) /  0.000000000000008740304300063083340121D+00 /
  data bifcs(  8) /  0.000000000000000009047915683496049529D+00 /
  data bifcs(  9) /  0.000000000000000000007249923164709251D+00 /
  data bifcs( 10) /  0.000000000000000000000004629576649604D+00 /
  data bifcs( 11) /  0.000000000000000000000000002411236436D+00 /
  data bifcs( 12) /  0.000000000000000000000000000001043825D+00 /
  data bifcs( 13) /  0.000000000000000000000000000000000382D+00 /

  data bigcs(  1) / -0.0971964404164435373897790974606802D+00 /
  data bigcs(  2) /  0.1495035768431670665710843445326264D+00 /
  data bigcs(  3) /  0.0031135253871213260419419176839631D+00 /
  data bigcs(  4) /  0.0000247085705798212967777021920569D+00 /
  data bigcs(  5) /  0.0000001029496277313786081987324295D+00 /
  data bigcs(  6) /  0.0000000002639703739869432892676778D+00 /
  data bigcs(  7) /  0.0000000000004582792707803206608181D+00 /
  data bigcs(  8) /  0.0000000000000005742829740893447321D+00 /
  data bigcs(  9) /  0.0000000000000000005438275385238549D+00 /
  data bigcs( 10) /  0.0000000000000000000004028347267083D+00 /
  data bigcs( 11) /  0.0000000000000000000000002397823826D+00 /
  data bigcs( 12) /  0.0000000000000000000000000001171956D+00 /
  data bigcs( 13) /  0.0000000000000000000000000000000479D+00 /

  data bif2cs(  1) /  0.32349398760352203352119193596266015D+00 /
  data bif2cs(  2) /  0.08629787153556355913888835323811100D+00 /
  data bif2cs(  3) /  0.00299402555265539742613821050727155D+00 /
  data bif2cs(  4) /  0.00005143052836466163720464316950821D+00 /
  data bif2cs(  5) /  0.00000052584025003681146026033098613D+00 /
  data bif2cs(  6) /  0.00000000356175137395770028102730600D+00 /
  data bif2cs(  7) /  0.00000000001714686400714584830518308D+00 /
  data bif2cs(  8) /  0.00000000000006166351969232555406693D+00 /
  data bif2cs(  9) /  0.00000000000000017191082154315985806D+00 /
  data bif2cs( 10) /  0.00000000000000000038236889518803943D+00 /
  data bif2cs( 11) /  0.00000000000000000000069424173624884D+00 /
  data bif2cs( 12) /  0.00000000000000000000000104833932510D+00 /
  data bif2cs( 13) /  0.00000000000000000000000000133721972D+00 /
  data bif2cs( 14) /  0.00000000000000000000000000000145986D+00 /
  data bif2cs( 15) /  0.00000000000000000000000000000000138D+00 /

  data big2cs(  1) /  1.606299946362129457759284537862622883D+00 /
  data big2cs(  2) /  0.744908881987608865201476685194753972D+00 /
  data big2cs(  3) /  0.047013873861027737964095177635353019D+00 /
  data big2cs(  4) /  0.001228442206254823907016188785848091D+00 /
  data big2cs(  5) /  0.000017322241225662362670987355613727D+00 /
  data big2cs(  6) /  0.000000152190165236801893711508366563D+00 /
  data big2cs(  7) /  0.000000000911356024911957704145528786D+00 /
  data big2cs(  8) /  0.000000000003954791842356644201722554D+00 /
  data big2cs(  9) /  0.000000000000013001737033862320007309D+00 /
  data big2cs( 10) /  0.000000000000000033493506858269079763D+00 /
  data big2cs( 11) /  0.000000000000000000069419094403694057D+00 /
  data big2cs( 12) /  0.000000000000000000000118248256604581D+00 /
  data big2cs( 13) /  0.000000000000000000000000168462493472D+00 /
  data big2cs( 14) /  0.000000000000000000000000000203684674D+00 /
  data big2cs( 15) /  0.000000000000000000000000000000211619D+00 /
  data big2cs( 16) /  0.000000000000000000000000000000000191D+00 /

  data bip2cs(  1) / -0.13269705443526630494937031210217135D+00 /
  data bip2cs(  2) / -0.00568443626045977481306046339037428D+00 /
  data bip2cs(  3) / -0.00015643601119611609623698471216660D+00 /
  data bip2cs(  4) / -0.00001136737203679562267336053207940D+00 /
  data bip2cs(  5) / -0.00000143464350991283669643136951338D+00 /
  data bip2cs(  6) / -0.00000018098531185164131868746481700D+00 /
  data bip2cs(  7) /  0.00000000926177343610865546229511422D+00 /
  data bip2cs(  8) /  0.00000001710005490720592181887296162D+00 /
  data bip2cs(  9) /  0.00000000476698163503781708252686849D+00 /
  data bip2cs( 10) / -0.00000000035195022023163141945397159D+00 /
  data bip2cs( 11) / -0.00000000058890614315886871574147635D+00 /
  data bip2cs( 12) / -0.00000000006678499607795537597612089D+00 /
  data bip2cs( 13) /  0.00000000006395565101720391190697713D+00 /
  data bip2cs( 14) /  0.00000000001554529427064394106403245D+00 /
  data bip2cs( 15) / -0.00000000000792396999744612971684001D+00 /
  data bip2cs( 16) / -0.00000000000258326242689717798947525D+00 /
  data bip2cs( 17) /  0.00000000000121655047787849117978773D+00 /
  data bip2cs( 18) /  0.00000000000038707207172899985942258D+00 /
  data bip2cs( 19) / -0.00000000000022487045479618229130656D+00 /
  data bip2cs( 20) / -0.00000000000004953476515684046293493D+00 /
  data bip2cs( 21) /  0.00000000000004563781601526912756017D+00 /
  data bip2cs( 22) /  0.00000000000000332998314345014118494D+00 /
  data bip2cs( 23) / -0.00000000000000921750185832874202719D+00 /
  data bip2cs( 24) /  0.00000000000000094156670658958205765D+00 /
  data bip2cs( 25) /  0.00000000000000167153952640716157721D+00 /
  data bip2cs( 26) / -0.00000000000000055134268782182410852D+00 /
  data bip2cs( 27) / -0.00000000000000022368651572006617795D+00 /
  data bip2cs( 28) /  0.00000000000000017486948976520089209D+00 /
  data bip2cs( 29) /  0.00000000000000000206518666352329750D+00 /
  data bip2cs( 30) / -0.00000000000000003973060018130712479D+00 /
  data bip2cs( 31) /  0.00000000000000001154836935724892335D+00 /
  data bip2cs( 32) /  0.00000000000000000553906053678276421D+00 /
  data bip2cs( 33) / -0.00000000000000000457174427396478267D+00 /
  data bip2cs( 34) /  0.00000000000000000026567111858284432D+00 /
  data bip2cs( 35) /  0.00000000000000000101599148154167823D+00 /
  data bip2cs( 36) / -0.00000000000000000044821231272196246D+00 /
  data bip2cs( 37) / -0.00000000000000000007959149661617295D+00 /
  data bip2cs( 38) /  0.00000000000000000014583615616165794D+00 /
  data bip2cs( 39) / -0.00000000000000000004015127893061405D+00 /
  data bip2cs( 40) / -0.00000000000000000002079152963743616D+00 /
  data bip2cs( 41) /  0.00000000000000000001972630449634388D+00 /
  data bip2cs( 42) / -0.00000000000000000000336033404001683D+00 /
  data bip2cs( 43) / -0.00000000000000000000376504832685507D+00 /
  data bip2cs( 44) /  0.00000000000000000000269935508825595D+00 /
  data bip2cs( 45) / -0.00000000000000000000026985946069808D+00 /
  data bip2cs( 46) / -0.00000000000000000000061794011788222D+00 /
  data bip2cs( 47) /  0.00000000000000000000038782693311711D+00 /
  data bip2cs( 48) / -0.00000000000000000000002420094005071D+00 /
  data bip2cs( 49) / -0.00000000000000000000009844051058925D+00 /
  data bip2cs( 50) /  0.00000000000000000000005954353358494D+00 /
  data bip2cs( 51) / -0.00000000000000000000000361274446366D+00 /
  data bip2cs( 52) / -0.00000000000000000000001552634578088D+00 /
  data bip2cs( 53) /  0.00000000000000000000000977819380304D+00 /
  data bip2cs( 54) / -0.00000000000000000000000092239447509D+00 /
  data bip2cs( 55) / -0.00000000000000000000000241545903934D+00 /
  data bip2cs( 56) /  0.00000000000000000000000169558652255D+00 /
  data bip2cs( 57) / -0.00000000000000000000000026762408641D+00 /
  data bip2cs( 58) / -0.00000000000000000000000036188116265D+00 /
  data bip2cs( 59) /  0.00000000000000000000000030372404951D+00 /
  data bip2cs( 60) / -0.00000000000000000000000007422876903D+00 /
  data bip2cs( 61) / -0.00000000000000000000000004930678544D+00 /
  data bip2cs( 62) /  0.00000000000000000000000005468790028D+00 /
  data bip2cs( 63) / -0.00000000000000000000000001920315188D+00 /
  data bip2cs( 64) / -0.00000000000000000000000000516335154D+00 /
  data bip2cs( 65) /  0.00000000000000000000000000957723167D+00 /
  data bip2cs( 66) / -0.00000000000000000000000000463659079D+00 /
  data bip2cs( 67) / -0.00000000000000000000000000004509226D+00 /
  data bip2cs( 68) /  0.00000000000000000000000000155617519D+00 /
  data bip2cs( 69) / -0.00000000000000000000000000104156509D+00 /
  data bip2cs( 70) /  0.00000000000000000000000000019565323D+00 /
  data bip2cs( 71) /  0.00000000000000000000000000021335380D+00 /
  data bip2cs( 72) / -0.00000000000000000000000000021461958D+00 /
  data bip2cs( 73) /  0.00000000000000000000000000007875791D+00 /
  data bip2cs( 74) /  0.00000000000000000000000000001713768D+00 /
  data bip2cs( 75) / -0.00000000000000000000000000003917137D+00 /
  data bip2cs( 76) /  0.00000000000000000000000000002233559D+00 /
  data bip2cs( 77) / -0.00000000000000000000000000000269383D+00 /
  data bip2cs( 78) / -0.00000000000000000000000000000577764D+00 /
  data bip2cs( 79) /  0.00000000000000000000000000000519650D+00 /
  data bip2cs( 80) / -0.00000000000000000000000000000183361D+00 /
  data bip2cs( 81) / -0.00000000000000000000000000000045763D+00 /
  data bip2cs( 82) /  0.00000000000000000000000000000099235D+00 /
  data bip2cs( 83) / -0.00000000000000000000000000000058938D+00 /
  data bip2cs( 84) /  0.00000000000000000000000000000009568D+00 /
  data bip2cs( 85) /  0.00000000000000000000000000000013758D+00 /
  data bip2cs( 86) / -0.00000000000000000000000000000014066D+00 /
  data bip2cs( 87) /  0.00000000000000000000000000000005964D+00 /
  data bip2cs( 88) /  0.00000000000000000000000000000000437D+00 /

  data bip1cs(  1) / -0.17291873510795537186124679823741003D+00 /
  data bip1cs(  2) / -0.01493584929846943639486231021818675D+00 /
  data bip1cs(  3) / -0.00054711049516785663990658697874460D+00 /
  data bip1cs(  4) /  0.00015379662929584083449573727856666D+00 /
  data bip1cs(  5) /  0.00001543534761921794131028948022869D+00 /
  data bip1cs(  6) / -0.00000654341138519060129226087106765D+00 /
  data bip1cs(  7) /  0.00000037280824078787032232152275240D+00 /
  data bip1cs(  8) /  0.00000020720783881887480080810710514D+00 /
  data bip1cs(  9) / -0.00000006581733364696191689495883922D+00 /
  data bip1cs( 10) /  0.00000000749267463539288212986048985D+00 /
  data bip1cs( 11) /  0.00000000111013368840707147698890101D+00 /
  data bip1cs( 12) / -0.00000000072651405529159512323880794D+00 /
  data bip1cs( 13) /  0.00000000017827235598470153962165668D+00 /
  data bip1cs( 14) / -0.00000000002173463524809506269656807D+00 /
  data bip1cs( 15) / -0.00000000000203020349653882594017049D+00 /
  data bip1cs( 16) /  0.00000000000193118272294077519319859D+00 /
  data bip1cs( 17) / -0.00000000000060449525048290296023117D+00 /
  data bip1cs( 18) /  0.00000000000012094496248933664277802D+00 /
  data bip1cs( 19) / -0.00000000000001251088360074479784619D+00 /
  data bip1cs( 20) / -0.00000000000000199173832424881344036D+00 /
  data bip1cs( 21) /  0.00000000000000151540816342864303038D+00 /
  data bip1cs( 22) / -0.00000000000000049768927059816240250D+00 /
  data bip1cs( 23) /  0.00000000000000011545959731810501403D+00 /
  data bip1cs( 24) / -0.00000000000000001863286862907983871D+00 /
  data bip1cs( 25) /  0.00000000000000000099330392344759104D+00 /
  data bip1cs( 26) /  0.00000000000000000068182083667412417D+00 /
  data bip1cs( 27) / -0.00000000000000000034854456479650551D+00 /
  data bip1cs( 28) /  0.00000000000000000010860382134235961D+00 /
  data bip1cs( 29) / -0.00000000000000000002599290185240166D+00 /
  data bip1cs( 30) /  0.00000000000000000000476895370459000D+00 /
  data bip1cs( 31) / -0.00000000000000000000051946940777177D+00 /
  data bip1cs( 32) / -0.00000000000000000000005925575044912D+00 /
  data bip1cs( 33) /  0.00000000000000000000005746008970972D+00 /
  data bip1cs( 34) / -0.00000000000000000000002186119806494D+00 /
  data bip1cs( 35) /  0.00000000000000000000000624124294738D+00 /
  data bip1cs( 36) / -0.00000000000000000000000146003421785D+00 /
  data bip1cs( 37) /  0.00000000000000000000000027493893904D+00 /
  data bip1cs( 38) / -0.00000000000000000000000003474678018D+00 /
  data bip1cs( 39) / -0.00000000000000000000000000109303694D+00 /
  data bip1cs( 40) /  0.00000000000000000000000000261972744D+00 /
  data bip1cs( 41) / -0.00000000000000000000000000112365018D+00 /
  data bip1cs( 42) /  0.00000000000000000000000000035152059D+00 /
  data bip1cs( 43) / -0.00000000000000000000000000009167601D+00 /
  data bip1cs( 44) /  0.00000000000000000000000000002040203D+00 /
  data bip1cs( 45) / -0.00000000000000000000000000000373038D+00 /
  data bip1cs( 46) /  0.00000000000000000000000000000046070D+00 /
  data bip1cs( 47) /  0.00000000000000000000000000000001748D+00 /

  data atr / 8.75069057084843450880771988210148D+00 /
  data btr /-2.09383632135605431360096498526268D+00 /
  data nbif / 0 /
  data nbif2 / 0 /
  data nbig / 0 /
  data nbig2 / 0 /
  data nbip1 / 0 /
  data nbip2 / 0 /
  data x2sml / 0.0D+00 /
  data x3sml / 0.0D+00 /
  data x32sml / 0.0D+00 /
  data xbig / 0.0D+00 /

  if ( nbif == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nbif = r8_inits ( bifcs, 13, eta )
    nbig = r8_inits ( bigcs, 13, eta )
    nbif2 = r8_inits ( bif2cs, 15, eta )
    nbig2 = r8_inits ( big2cs, 16, eta )
    nbip1 = r8_inits ( bip1cs, 47, eta )
    nbip2 = r8_inits ( bip2cs, 88, eta )
    x2sml = sqrt ( eta )
    x3sml = eta ** 0.3333D+00
    x32sml = 1.3104D+00 * x3sml * x3sml
    xbig = r8_mach ( 2 ) ** 0.6666D+00
  end if

  if ( x < -1.0D+00 ) then
    call r8_admp ( x, xn, phi )
    r8_bide = xn * sin ( phi )
  else if ( abs ( x ) <= x2sml ) then
    x2 = 0.0D+00
    x3 = 0.0D+00
    r8_bide = x2 * ( r8_csevl ( x3, bifcs, nbif ) &
      + 0.25D+00 ) + r8_csevl ( x3, bigcs, nbig ) &
      + 0.5D+00
    if ( x32sml < x ) then
      r8_bide = r8_bide &
        * exp ( - 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
    end if
  else if ( abs ( x ) <= x3sml ) then
    x2 = x * x
    x3 = 0.0D+00
    r8_bide = x2 * ( r8_csevl ( x3, bifcs, nbif ) &
      + 0.25D+00 ) + r8_csevl ( x3, bigcs, nbig ) &
      + 0.5D+00
    if ( x32sml < x ) then
      r8_bide = r8_bide &
        * exp ( - 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
    end if
  else if ( x <= 1.0D+00 ) then
    x2 = x * x
    x3 = x * x * x
    r8_bide = x2 * ( r8_csevl ( x3, bifcs, nbif ) &
      + 0.25D+00 ) + r8_csevl ( x3, bigcs, nbig ) &
      + 0.5D+00
    if ( x32sml < x ) then
      r8_bide = r8_bide &
        * exp ( - 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
    end if
  else if ( x <= 2.0D+00 ) then
    z = ( 2.0D+00 * x * x * x - 9.0D+00 ) / 7.0D+00
    r8_bide = exp ( - 2.0D+00 * x * sqrt ( x ) / 3.0D+00 ) &
      * ( x * x * ( 0.25D+00 + r8_csevl ( z, bif2cs, nbif2 ) ) &
      + 0.5D+00 + r8_csevl ( z, big2cs, nbig2 ) )
  else if ( x <= 4.0D+00 ) then
    sqrtx = sqrt ( x )
    z = atr / x / sqrtx + btr
    r8_bide = ( 0.625D+00 + r8_csevl ( z, bip1cs, nbip1 ) ) &
      * sqrt ( sqrtx )
  else if ( x <= xbig ) then
    sqrtx = sqrt ( x )
    z = 16.0D+00 / x / sqrtx - 1.0D+00
    r8_bide = ( 0.625D+00 + r8_csevl ( z, bip2cs, nbip2 ) ) * sqrt ( sqrtx )
  else
    sqrtx = sqrt ( x )
    z = -1.0D+00
    r8_bide = ( 0.625D+00 + r8_csevl ( z, bip2cs, nbip2 ) ) * sqrt ( sqrtx )
  end if

  return
end
function r8_bie ( x )

!*****************************************************************************80
!
!! R8_BIE evaluates the exponentially scaled Airy function Bi of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BIE, the exponentially scaled Airy 
!    function Bi of X.
!
  implicit none

  real ( kind = 8 ) atr
  real ( kind = 8 ) bif2cs(15)
  real ( kind = 8 ) bifcs(13)
  real ( kind = 8 ) big2cs(15)
  real ( kind = 8 ) bigcs(13)
  real ( kind = 8 ) bip1cs(47)
  real ( kind = 8 ) bip2cs(88)
  real ( kind = 8 ) btr
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nbif
  integer ( kind = 4 ) nbif2
  integer ( kind = 4 ) nbig
  integer ( kind = 4 ) nbig2
  integer ( kind = 4 ) nbip1
  integer ( kind = 4 ) nbip2
  real ( kind = 8 ) r8_bie
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sqrtx
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) x32sml
  real ( kind = 8 ) x3sml
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xm
  real ( kind = 8 ) z

  save atr
  save bif2cs
  save bifcs
  save big2cs
  save bigcs
  save bip1cs
  save bip2cs
  save btr
  save nbif
  save nbif2
  save nbig
  save nbig2
  save nbip1
  save nbip2
  save x32sml
  save x3sml
  save xbig


  data bifcs(  1) / -0.16730216471986649483537423928176D-01 /
  data bifcs(  2) / +0.10252335834249445611426362777757D+00 /
  data bifcs(  3) / +0.17083092507381516539429650242013D-02 /
  data bifcs(  4) / +0.11862545467744681179216459210040D-04 /
  data bifcs(  5) / +0.44932907017792133694531887927242D-07 /
  data bifcs(  6) / +0.10698207143387889067567767663628D-09 /
  data bifcs(  7) / +0.17480643399771824706010517628573D-12 /
  data bifcs(  8) / +0.20810231071761711025881891834399D-15 /
  data bifcs(  9) / +0.18849814695665416509927971733333D-18 /
  data bifcs( 10) / +0.13425779173097804625882666666666D-21 /
  data bifcs( 11) / +0.77159593429658887893333333333333D-25 /
  data bifcs( 12) / +0.36533879617478566399999999999999D-28 /
  data bifcs( 13) / +0.14497565927953066666666666666666D-31 /

  data bigcs(  1) / +0.22466223248574522283468220139024D-01 /
  data bigcs(  2) / +0.37364775453019545441727561666752D-01 /
  data bigcs(  3) / +0.44476218957212285696215294326639D-03 /
  data bigcs(  4) / +0.24708075636329384245494591948882D-05 /
  data bigcs(  5) / +0.79191353395149635134862426285596D-08 /
  data bigcs(  6) / +0.16498079851827779880887872402706D-10 /
  data bigcs(  7) / +0.24119906664835455909247501122841D-13 /
  data bigcs(  8) / +0.26103736236091436985184781269333D-16 /
  data bigcs(  9) / +0.21753082977160323853123792000000D-19 /
  data bigcs( 10) / +0.14386946400390433219483733333333D-22 /
  data bigcs( 11) / +0.77349125612083468629333333333333D-26 /
  data bigcs( 12) / +0.34469292033849002666666666666666D-29 /
  data bigcs( 13) / +0.12938919273216000000000000000000D-32 /

  data bif2cs(  1) / +0.0998457269381604104468284257993D+00 /
  data bif2cs(  2) / +0.47862497786300553772211467318231D+00 /
  data bif2cs(  3) / +0.25155211960433011771324415436675D-01 /
  data bif2cs(  4) / +0.58206938852326456396515697872216D-03 /
  data bif2cs(  5) / +0.74997659644377865943861457378217D-05 /
  data bif2cs(  6) / +0.61346028703493836681403010356474D-07 /
  data bif2cs(  7) / +0.34627538851480632900434268733359D-09 /
  data bif2cs(  8) / +0.14288910080270254287770846748931D-11 /
  data bif2cs(  9) / +0.44962704298334641895056472179200D-14 /
  data bif2cs( 10) / +0.11142323065833011708428300106666D-16 /
  data bif2cs( 11) / +0.22304791066175002081517866666666D-19 /
  data bif2cs( 12) / +0.36815778736393142842922666666666D-22 /
  data bif2cs( 13) / +0.50960868449338261333333333333333D-25 /
  data bif2cs( 14) / +0.60003386926288554666666666666666D-28 /
  data bif2cs( 15) / +0.60827497446570666666666666666666D-31 /

  data big2cs(  1) / +0.033305662145514340465176188111647D+00 /
  data big2cs(  2) / +0.161309215123197067613287532084943D+00 /
  data big2cs(  3) / +0.631900730961342869121615634921173D-02 /
  data big2cs(  4) / +0.118790456816251736389780192304567D-03 /
  data big2cs(  5) / +0.130453458862002656147116485012843D-05 /
  data big2cs(  6) / +0.937412599553521729546809615508936D-08 /
  data big2cs(  7) / +0.474580188674725153788510169834595D-10 /
  data big2cs(  8) / +0.178310726509481399800065667560946D-12 /
  data big2cs(  9) / +0.516759192784958180374276356640000D-15 /
  data big2cs( 10) / +0.119004508386827125129496251733333D-17 /
  data big2cs( 11) / +0.222982880666403517277063466666666D-20 /
  data big2cs( 12) / +0.346551923027689419722666666666666D-23 /
  data big2cs( 13) / +0.453926336320504514133333333333333D-26 /
  data big2cs( 14) / +0.507884996513522346666666666666666D-29 /
  data big2cs( 15) / +0.491020674696533333333333333333333D-32 /

  data bip1cs(  1) / -0.83220474779434474687471864707973D-01 /
  data bip1cs(  2) / +0.11461189273711742889920226128031D-01 /
  data bip1cs(  3) / +0.42896440718911509494134472566635D-03 /
  data bip1cs(  4) / -0.14906639379950514017847677732954D-03 /
  data bip1cs(  5) / -0.13076597267876290663136340998881D-04 /
  data bip1cs(  6) / +0.63275983961030344754535716032494D-05 /
  data bip1cs(  7) / -0.42226696982681924884778515889433D-06 /
  data bip1cs(  8) / -0.19147186298654689632835494181277D-06 /
  data bip1cs(  9) / +0.64531062845583173611038157880934D-07 /
  data bip1cs( 10) / -0.78448546771397719289748310448628D-08 /
  data bip1cs( 11) / -0.96077216623785085879198533565432D-09 /
  data bip1cs( 12) / +0.70004713316443966339006074402068D-09 /
  data bip1cs( 13) / -0.17731789132814932022083128056698D-09 /
  data bip1cs( 14) / +0.22720894783465236347282126389311D-10 /
  data bip1cs( 15) / +0.16540456313972049847032860681891D-11 /
  data bip1cs( 16) / -0.18517125559292316390755369896693D-11 /
  data bip1cs( 17) / +0.59576312477117290165680715534277D-12 /
  data bip1cs( 18) / -0.12194348147346564781055769498986D-12 /
  data bip1cs( 19) / +0.13347869253513048815386347813597D-13 /
  data bip1cs( 20) / +0.17278311524339746664384792889731D-14 /
  data bip1cs( 21) / -0.14590732013016720735268871713166D-14 /
  data bip1cs( 22) / +0.49010319927115819978994989520104D-15 /
  data bip1cs( 23) / -0.11556545519261548129262972762521D-15 /
  data bip1cs( 24) / +0.19098807367072411430671732441524D-16 /
  data bip1cs( 25) / -0.11768966854492179886913995957862D-17 /
  data bip1cs( 26) / -0.63271925149530064474537459677047D-18 /
  data bip1cs( 27) / +0.33861838880715361614130191322316D-18 /
  data bip1cs( 28) / -0.10725825321758625254992162219622D-18 /
  data bip1cs( 29) / +0.25995709605617169284786933115562D-19 /
  data bip1cs( 30) / -0.48477583571081193660962309494101D-20 /
  data bip1cs( 31) / +0.55298913982121625361505513198933D-21 /
  data bip1cs( 32) / +0.49421660826069471371748197444266D-22 /
  data bip1cs( 33) / -0.55162121924145707458069720814933D-22 /
  data bip1cs( 34) / +0.21437560417632550086631884499626D-22 /
  data bip1cs( 35) / -0.61910313387655605798785061137066D-23 /
  data bip1cs( 36) / +0.14629362707391245659830967336959D-23 /
  data bip1cs( 37) / -0.27918484471059005576177866069333D-24 /
  data bip1cs( 38) / +0.36455703168570246150906795349333D-25 /
  data bip1cs( 39) / +0.58511821906188711839382459733333D-27 /
  data bip1cs( 40) / -0.24946950487566510969745047551999D-26 /
  data bip1cs( 41) / +0.10979323980338380977919579477333D-26 /
  data bip1cs( 42) / -0.34743388345961115015034088106666D-27 /
  data bip1cs( 43) / +0.91373402635349697363171082240000D-28 /
  data bip1cs( 44) / -0.20510352728210629186247720959999D-28 /
  data bip1cs( 45) / +0.37976985698546461748651622399999D-29 /
  data bip1cs( 46) / -0.48479458497755565887848448000000D-30 /
  data bip1cs( 47) / -0.10558306941230714314205866666666D-31 /

  data bip2cs(  1) / -0.11359673758598867913797310895527D+00 /
  data bip2cs(  2) / +0.41381473947881595760052081171444D-02 /
  data bip2cs(  3) / +0.13534706221193329857696921727508D-03 /
  data bip2cs(  4) / +0.10427316653015353405887183456780D-04 /
  data bip2cs(  5) / +0.13474954767849907889589911958925D-05 /
  data bip2cs(  6) / +0.16965374054383983356062511163756D-06 /
  data bip2cs(  7) / -0.10096500865641624301366228396373D-07 /
  data bip2cs(  8) / -0.16729119493778475127836973095943D-07 /
  data bip2cs(  9) / -0.45815364485068383217152795613391D-08 /
  data bip2cs( 10) / +0.37366813665655477274064749384284D-09 /
  data bip2cs( 11) / +0.57669303201452448119584643502111D-09 /
  data bip2cs( 12) / +0.62181265087850324095393408792371D-10 /
  data bip2cs( 13) / -0.63294120282743068241589177281354D-10 /
  data bip2cs( 14) / -0.14915047908598767633999091989487D-10 /
  data bip2cs( 15) / +0.78896213942486771938172394294891D-11 /
  data bip2cs( 16) / +0.24960513721857797984888064000127D-11 /
  data bip2cs( 17) / -0.12130075287291659477746664734814D-11 /
  data bip2cs( 18) / -0.37404939108727277887343460402716D-12 /
  data bip2cs( 19) / +0.22377278140321476798783446931091D-12 /
  data bip2cs( 20) / +0.47490296312192466341986077472514D-13 /
  data bip2cs( 21) / -0.45261607991821224810605655831294D-13 /
  data bip2cs( 22) / -0.30172271841986072645112245876020D-14 /
  data bip2cs( 23) / +0.91058603558754058327592683478908D-14 /
  data bip2cs( 24) / -0.98149238033807062926643864207709D-15 /
  data bip2cs( 25) / -0.16429400647889465253601245251589D-14 /
  data bip2cs( 26) / +0.55334834214274215451182114635164D-15 /
  data bip2cs( 27) / +0.21750479864482655984374381998156D-15 /
  data bip2cs( 28) / -0.17379236200220656971287029558087D-15 /
  data bip2cs( 29) / -0.10470023471443714959283909313604D-17 /
  data bip2cs( 30) / +0.39219145986056386925441403311462D-16 /
  data bip2cs( 31) / -0.11621293686345196925824005665910D-16 /
  data bip2cs( 32) / -0.54027474491754245533735411307773D-17 /
  data bip2cs( 33) / +0.45441582123884610882675428553304D-17 /
  data bip2cs( 34) / -0.28775599625221075729427585480086D-18 /
  data bip2cs( 35) / -0.10017340927225341243596162960440D-17 /
  data bip2cs( 36) / +0.44823931215068369856332561906313D-18 /
  data bip2cs( 37) / +0.76135968654908942328948982366775D-19 /
  data bip2cs( 38) / -0.14448324094881347238956060145422D-18 /
  data bip2cs( 39) / +0.40460859449205362251624847392112D-19 /
  data bip2cs( 40) / +0.20321085700338446891325190707277D-19 /
  data bip2cs( 41) / -0.19602795471446798718272758041962D-19 /
  data bip2cs( 42) / +0.34273038443944824263518958211738D-20 /
  data bip2cs( 43) / +0.37023705853905135480024651593154D-20 /
  data bip2cs( 44) / -0.26879595172041591131400332966712D-20 /
  data bip2cs( 45) / +0.28121678463531712209714454683364D-21 /
  data bip2cs( 46) / +0.60933963636177797173271119680329D-21 /
  data bip2cs( 47) / -0.38666621897150844994172977893413D-21 /
  data bip2cs( 48) / +0.25989331253566943450895651927228D-22 /
  data bip2cs( 49) / +0.97194393622938503767281175216084D-22 /
  data bip2cs( 50) / -0.59392817834375098415630478204591D-22 /
  data bip2cs( 51) / +0.38864949977113015409591960439444D-23 /
  data bip2cs( 52) / +0.15334307393617272869721512868769D-22 /
  data bip2cs( 53) / -0.97513555209762624036336521409724D-23 /
  data bip2cs( 54) / +0.96340644440489471424741339383726D-24 /
  data bip2cs( 55) / +0.23841999400208880109946748792454D-23 /
  data bip2cs( 56) / -0.16896986315019706184848044205207D-23 /
  data bip2cs( 57) / +0.27352715888928361222578444801478D-24 /
  data bip2cs( 58) / +0.35660016185409578960111685025730D-24 /
  data bip2cs( 59) / -0.30234026608258827249534280666954D-24 /
  data bip2cs( 60) / +0.75002041605973930653144204823232D-25 /
  data bip2cs( 61) / +0.48403287575851388827455319838748D-25 /
  data bip2cs( 62) / -0.54364137654447888432698010297766D-25 /
  data bip2cs( 63) / +0.19281214470820962653345978809756D-25 /
  data bip2cs( 64) / +0.50116355020532656659611814172172D-26 /
  data bip2cs( 65) / -0.95040744582693253786034620869972D-26 /
  data bip2cs( 66) / +0.46372646157101975948696332245611D-26 /
  data bip2cs( 67) / +0.21177170704466954163768170577046D-28 /
  data bip2cs( 68) / -0.15404850268168594303692204548726D-26 /
  data bip2cs( 69) / +0.10387944293201213662047889194441D-26 /
  data bip2cs( 70) / -0.19890078156915416751316728235153D-27 /
  data bip2cs( 71) / -0.21022173878658495471177044522532D-27 /
  data bip2cs( 72) / +0.21353099724525793150633356670491D-27 /
  data bip2cs( 73) / -0.79040810747961342319023537632627D-28 /
  data bip2cs( 74) / -0.16575359960435585049973741763592D-28 /
  data bip2cs( 75) / +0.38868342850124112587625586496537D-28 /
  data bip2cs( 76) / -0.22309237330896866182621562424717D-28 /
  data bip2cs( 77) / +0.27777244420176260265625977404382D-29 /
  data bip2cs( 78) / +0.57078543472657725368712433782772D-29 /
  data bip2cs( 79) / -0.51743084445303852800173371555280D-29 /
  data bip2cs( 80) / +0.18413280751095837198450927071569D-29 /
  data bip2cs( 81) / +0.44422562390957094598544071068647D-30 /
  data bip2cs( 82) / -0.98504142639629801547464958226943D-30 /
  data bip2cs( 83) / +0.58857201353585104884754198881995D-30 /
  data bip2cs( 84) / -0.97636075440429787961402312628595D-31 /
  data bip2cs( 85) / -0.13581011996074695047063597884122D-30 /
  data bip2cs( 86) / +0.13999743518492413270568048380345D-30 /
  data bip2cs( 87) / -0.59754904545248477620884562981118D-31 /
  data bip2cs( 88) / -0.40391653875428313641045327529856D-32 /

  data atr / 8.75069057084843450880771988210148D+00 /
  data btr / -2.09383632135605431360096498526268D+00 /

  data nbif / 0 /
  data nbif2 / 0 /
  data nbig / 0 /
  data nbig2 / 0 /
  data nbip1 / 0 /
  data nbip2 / 0 /
  data x32sml / 0.0D+00 /
  data x3sml / 0.0D+00 /
  data xbig / 0.0D+00 /

  if ( nbif == 0 ) then
    eta = 0.1D+00  * r8_mach ( 3 )
    nbif = r8_inits ( bifcs, 13, eta )
    nbig = r8_inits ( bigcs, 13, eta )
    nbif2 = r8_inits ( bif2cs, 15, eta )
    nbig2 = r8_inits ( big2cs, 15, eta )
    nbip1 = r8_inits ( bip1cs, 47, eta )
    nbip2 = r8_inits ( bip2cs, 88, eta )
    x3sml = eta ** 0.3333D+00
    x32sml = 1.3104D+00 * x3sml * x3sml
    xbig = r8_mach ( 2 ) ** 0.6666D+00
  end if

  if ( x < - 1.0D+00 ) then
    call r8_aimp ( x, xm, theta )
    r8_bie = xm * sin ( theta )
  else if ( abs ( x ) <= x3sml ) then
    z = 0.0D+00
    r8_bie = 0.625D+00 + r8_csevl ( z, bifcs, nbif ) &
      + x * ( 0.4375D+00 + r8_csevl ( z, bigcs, nbig ) )
    if (  x32sml <= x ) then
      r8_bie = r8_bie &
        * exp ( - 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
    end if
  else if ( x <= 1.0D+00 ) then
    z = x * x * x
    r8_bie = 0.625D+00 + r8_csevl ( z, bifcs, nbif ) &
      + x * ( 0.4375D+00 + r8_csevl ( z, bigcs, nbig ) )
    if (  x32sml <= x ) then
      r8_bie = r8_bie &
        * exp ( - 2.0D+00 * x * sqrt ( x ) / 3.0D+00 )
    end if
  else if ( x <= 2.0D+00 )  then
    z = ( 2.0D+00 * x * x * x - 9.0D+00 ) / 7.0D+00
    r8_bie = exp ( - 2.0D+00 * x * sqrt ( x ) / 3.0D+00 ) &
      * ( 1.125D+00 + r8_csevl ( z, bif2cs, nbif2 ) &
      + x * ( 0.625D+00 + r8_csevl ( z, big2cs, nbig2 ) ) )
  else if ( x <= 4.0D+00 ) then
    sqrtx = sqrt ( x )
    z = atr / x / sqrtx + btr
    r8_bie = ( 0.625D+00 &
      + r8_csevl ( z, bip1cs, nbip1 ) ) / sqrt ( sqrtx )
  else if ( x < xbig ) then
    sqrtx = sqrt ( x )
    z = 16.0D+00 / ( x * sqrtx ) - 1.0D+00
    r8_bie = ( 0.625D+00 &
      + r8_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx )
  else
    sqrtx = sqrt ( x )
    z = - 1.0D+00
    r8_bie = ( 0.625D+00 &
      + r8_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx )
  end if

  return
end
function r8_binom ( n, m )

!*****************************************************************************80
!
!! R8_BINOM evaluates the binomial coefficient using R8 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, M, the arguments.
!
!    Output, real ( kind = 8 ) R8_BINOM, the binomial coefficient.
!
  implicit none

  real ( kind = 8 ) bilnmx
  real ( kind = 8 ) corr
  real ( kind = 8 ) fintmx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_binom
  real ( kind = 8 ) r8_lgmc
  real ( kind = 8 ) r8_lnrel
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sq2pil
  real ( kind = 8 ) xk
  real ( kind = 8 ) xn
  real ( kind = 8 ) xnk

  save bilnmx
  save fintmx
  save sq2pil

  data bilnmx / 0.0D+00 /
  data fintmx / 0.0D+00 /
  data sq2pil / 0.91893853320467274178032973640562D+00 /

  if ( bilnmx == 0.0D+00 ) then
    bilnmx = log ( r8_mach ( 2 ) ) - 0.0001D+00
    fintmx = 0.9D+00 / r8_mach ( 3 )
  end if

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BINOM - Fatal error!'
    write ( *, '(a)' ) '  N < 0.'
    stop
  end if

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BINOM - Fatal error!'
    write ( *, '(a)' ) '  M < 0.'
    stop
  end if

  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BINOM - Fatal error!'
    write ( *, '(a)' ) '  N < M.'
    stop
  end if

  k = min ( m, n - m )

  if ( k <= 20 .and. real ( k, kind = 8 ) &
    * log ( real ( max ( n, 1 ), kind = 8 ) ) <= bilnmx ) then

    r8_binom = 1.0D+00

    do i = 1, k
      r8_binom = r8_binom * real ( n - i + 1, kind = 8 ) / real ( i, kind = 8 )
    end do

  else

    if ( k < 9 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_BINOM - Fatal error!'
      write ( *, '(a)' ) '  Result overflows.'
      write ( *, '(a)' ) '  N or M is too big.'
      stop
    end if

    xn = real ( n + 1, kind = 8 )
    xk = real ( k + 1, kind = 8 )
    xnk = real ( n - k + 1, kind = 8 )

    corr = r8_lgmc ( xn ) - r8_lgmc ( xk ) - r8_lgmc ( xnk )

    r8_binom = xk * log ( xnk / xk ) &
      - xn * r8_lnrel ( - ( xk - 1.0D+00 ) / xn ) &
      - 0.5D+00 * log ( xn * xnk / xk ) + 1.0D+00 - sq2pil + corr

    if ( bilnmx < r8_binom ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_BINOM - Fatal error!'
      write ( *, '(a)' ) '  Result overflows.'
      write ( *, '(a)' ) '  N or M is too big.'
      stop
    end if

    r8_binom = exp ( r8_binom )

  end if

  if ( r8_binom < fintmx ) then
    r8_binom = aint ( r8_binom + 0.5D+00 )
  end if

  return
end
function r8_cbrt ( x )

!*****************************************************************************80
!
!! R8_CBRT computes the cube root of an R8.
!
!  Discussion:
!
!    The approximation is a generalized Chebyshev series converted
!    to polynomial form.  The approximation is nearly best in the 
!    sense of relative error with 4.085 digits accuracy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose square root is desired.
!
!    Output, real ( kind = 8 ) R8_CBRT, the cube root of X.
!
  implicit none

  real ( kind = 8 ) cbrt2(5)
  integer ( kind = 4 ) irem
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) ixpnt
  integer ( kind = 4 ) n
  integer ( kind = 4 ) niter
  real ( kind = 8 ) r8_cbrt
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_pak
  real ( kind = 8 ) value
  real ( kind = 8 ) vsq
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  save cbrt2
  save niter

  data cbrt2(1) / 0.62996052494743658238360530363911D+00 /
  data cbrt2(2) / 0.79370052598409973737585281963615D+00 /
  data cbrt2(3) / 1.0D+00 /
  data cbrt2(4) / 1.25992104989487316476721060727823D+00 /
  data cbrt2(5) / 1.58740105196819947475170563927231D+00 /

  data niter / 0 /

  if ( niter == 0 ) then
    niter = int ( 1.443D+00 * log ( - 0.106D+00 &
      * log ( 0.1D+00 * r8_mach ( 3 ) ) ) + 1.0D+00 )
  end if

  value = 0.0D+00

  if ( x /= 0.0D+00 ) then

    call r8_upak ( abs ( x ), y, n )
    ixpnt = n / 3
    irem = n - 3 * ixpnt + 3

    value = 0.439581D+00 + y * ( &
            0.928549D+00 + y * ( &
          - 0.512653D+00 + y * &
            0.144586D+00 ) )

    do iter = 1, niter
      vsq = value * value
      value = value + ( y - value * vsq ) / ( 3.0D+00 * vsq )
    end do

    if ( x < 0.0D+00 ) then
      value = - abs ( value )
    else
      value = + abs ( value )
    end if

    value = r8_pak ( cbrt2(irem) * value, ixpnt )

  end if

  r8_cbrt = value

  return
end
function r8_chi ( x )

!*****************************************************************************80
!
!! R8_CHI evaluates the hyperbolic cosine integral of an R8 argument.
!
!  Discussion:
!
!    The hyperbolic cosine integral is defined by
!
!      CHI(X) = gamma + log ( x ) 
!        + integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
!
!    where gamma is Euler's constant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_CHI, the hyperbolic cosine integral
!    evaluated at X.
!
  implicit none

  real ( kind = 8 ) r8_chi
  real ( kind = 8 ) r8_e1
  real ( kind = 8 ) r8_ei
  real ( kind = 8 ) x

  r8_chi = 0.5D+00 * ( r8_ei ( x ) - r8_e1 ( x ) )

  return
end
function r8_chu ( a, b, x )

!*****************************************************************************80
!
!! R8_CHU evaluates the confluent hypergeometric function of R8 arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_CHU, the function value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a0
  real ( kind = 8 ) aintb
  real ( kind = 8 ) alnx
  real ( kind = 8 ) b
  real ( kind = 8 ) b0
  real ( kind = 8 ) beps
  real ( kind = 8 ) c0
  real ( kind = 8 ) eps
  real ( kind = 8 ) factor
  real ( kind = 8 ) gamri1
  real ( kind = 8 ) gamrni
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istrt
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) pch1ai
  real ( kind = 8 ) pch1i
  real ( kind = 8 ) pi
  real ( kind = 8 ) pochai
  real ( kind = 8 ) r8_chu
  real ( kind = 8 ) r8_chu_scaled
  real ( kind = 8 ) r8_exprel
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_gamr
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) r8_poch
  real ( kind = 8 ) r8_poch1
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) xeps1
  real ( kind = 8 ) xi
  real ( kind = 8 ) xi1
  real ( kind = 8 ) xn
  real ( kind = 8 ) xtoeps

  save eps
  save pi

  data eps / 0.0D+00 /
  data pi / 3.141592653589793238462643383279503D+00 /

  if ( eps == 0.0D+00 ) then
    eps = r8_mach ( 3 )
  end if

  if ( x < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CHU - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop
  end if

  if ( x == 0.0D+00 ) then
    if ( 1.0D+00 <= b ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_CHU - Fatal error!'
      write ( *, '(a)' ) '  X = 0 and 1 <= B.'
      stop
    end if
    r8_chu = r8_gamma ( 1.0D+00 - b ) / r8_gamma ( 1.0D+00 + a - b )
    return
  end if

  if ( max ( abs ( a ), 1.0D+00 ) &
    * max ( abs ( 1.0D+00 + a - b ), 1.0D+00 ) &
    < 0.99D+00 * abs ( x ) ) then
    r8_chu = x ** ( - a ) * r8_chu_scaled ( a, b, x )
    return
  end if
!
!  The ascending series will be used, because the descending rational
!  approximation (which is based on the asymptotic series) is unstable.
!
  if ( 0.0D+00 <= b ) then
    aintb = aint ( b + 0.5D+00 )
  else
    aintb = aint ( b - 0.5D+00 )
  end if
  beps = b - aintb
  n = aintb
  alnx = log ( x )
  xtoeps = exp ( - beps * alnx )
!
!  Evaluate the finite sum.
!
!  Consider the case b < 1.0 first.
!
  if ( n < 1 ) then

    sum = 1.0D+00
 
    t = 1.0D+00
    m = - n
    do i = 1, m
      xi1 = real ( i - 1, kind = 8 )
      t = t * ( a + xi1 ) * x / ( ( b + xi1 ) * ( xi1 + 1.0D+00 ) )
      sum = sum + t
    end do

    sum = r8_poch ( 1.0D+00 + a - b, - a ) * sum
!
!  Now consider the case 1 <= b.
!
  else

    sum = 0.0D+00
    m = n - 2

    if ( 0 <= m ) then

      t = 1.0D+00
      sum = 1.0D+00

      do i = 1, m
       xi = real ( i, kind = 8 )
        t = t * ( a - b + xi ) * x / ( ( 1.0D+00 - b + xi ) * xi )
        sum = sum + t
      end do

      sum = r8_gamma ( b - 1.0D+00 ) * r8_gamr ( a ) &
        * x ** ( 1 - n ) * xtoeps * sum

    end if

  end if
!
!  Next evaluate the infinite sum.
!
  if ( n < 1 ) then
    istrt = 1 - n
  else
    istrt = 0
  end if

  xi = real ( istrt, kind = 8 )

  factor = r8_mop ( n ) * r8_gamr ( 1.0D+00 + a - b ) &
    * x ** istrt

  if ( beps /= 0.0D+00 ) then
    factor = factor * beps * pi / sin ( beps * pi )
  end if

  pochai = r8_poch ( a, xi )
  gamri1 = r8_gamr ( xi + 1.0D+00 )
  gamrni = r8_gamr ( aintb + xi )
  b0 = factor * r8_poch ( a, xi - beps ) &
    * gamrni * r8_gamr ( xi + 1.0D+00 - beps )
!
!  x^(-beps) is close to 1.0, so we must be careful in evaluating the
!  differences.
!
  if ( abs ( xtoeps - 1.0D+00 ) <= 0.5D+00 ) then

    pch1ai = r8_poch1 ( a + xi, -beps )
    pch1i = r8_poch1 ( xi + 1.0D+00 - beps, beps )
    c0 = factor * pochai * gamrni * gamri1 * ( &
      - r8_poch1 ( b + xi,- beps ) + pch1ai &
      - pch1i + beps * pch1ai * pch1i )
!
!  xeps1 = (1.0 - x^(-beps))/beps = (x^(-beps) - 1.0)/(-beps)
!
    xeps1 = alnx* r8_exprel ( - beps * alnx )

    r8_chu = sum + c0 + xeps1 * b0
    xn = real ( n, kind = 8 )

    do i = 1, 1000
      xi = real ( istrt + i, kind = 8 )
      xi1 = real ( istrt + i - 1, kind = 8 )
      b0 = ( a + xi1 - beps ) * b0 * x &
        / ( ( xn + xi1 ) * ( xi - beps ) )
      c0 = ( a + xi1 ) * c0 * x / ( ( b + xi1) * xi ) &
        - ( ( a - 1.0D+00 ) * ( xn + 2.0D+00 * xi - 1.0D+00 ) &
        + xi * ( xi - beps ) ) * b0 &
        / ( xi * ( b + xi1 ) * ( a + xi1 - beps ) )
      t = c0 + xeps1 * b0
      r8_chu = r8_chu + t
      if ( abs ( t ) < eps * abs ( r8_chu ) ) then
        return
      end if
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CHU - Fatal error!'
    write ( *, '(a)' ) '  No convergence in 1000 terms.'
    stop

  end if
!
!  x^(-beps) is very different from 1.0, so the straightforward
!  formulation is stable.
!
  a0 = factor * pochai * r8_gamr ( b + xi ) * gamri1 / beps
  b0 = xtoeps * b0 / beps

  r8_chu = sum + a0 - b0

  do i = 1, 1000
    xi = real ( istrt + i, kind = 8 )
    xi1 = real ( istrt + i - 1, kind = 8 )
    a0 = ( a + xi1 ) * a0 * x / ( ( b + xi1 ) * xi )
    b0 = ( a + xi1 - beps ) * b0 * x &
      / ( ( aintb + xi1 ) * ( xi - beps ) )
    t = a0 - b0
    r8_chu = r8_chu + t
    if ( abs ( t ) < eps * abs ( r8_chu ) ) then
      return
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_CHU - Fatal error!'
  write ( *, '(a)' ) '  No convergence in 1000 terms.'
  stop

end
function r8_chu_scaled ( a, b, z )

!*****************************************************************************80
!
!! R8_CHU_SCALED: scaled confluent hypergeometric function of R8 arguments.
!
!  Discussion:
!
!    Evaluate, for large z, z**a * u(a,b,z)  where U is the logarithmic
!    confluent hypergeometric function.  A rational approximation due to
!    Y L Luke is used.  When U is not in the asymptotic region, that is, when A
!    or B is large compared with Z, considerable significance loss occurs.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters.
!
!    Input, real ( kind = 8 ) Z, the argument.
!
!    Output, real ( kind = 8 ) R8_CHU_SCALED, the function value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) aa(4)
  real ( kind = 8 ) ab
  real ( kind = 8 ) anbn
  real ( kind = 8 ) b
  real ( kind = 8 ) bb(4)
  real ( kind = 8 ) bp
  real ( kind = 8 ) c2
  real ( kind = 8 ) ct1
  real ( kind = 8 ) ct2
  real ( kind = 8 ) ct3
  real ( kind = 8 ) d1z
  real ( kind = 8 ), save :: eps = 0.0D+00
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) g3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_chu_scaled
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sab
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) x2i1
  real ( kind = 8 ) z

  if ( eps == 0.0D+00 ) then
    eps = 4.0D+00 * r8_mach ( 4 )
    sqeps = sqrt ( r8_mach ( 4 ) )
  end if

  bp = 1.0D+00 + a - b
  ab = a * bp
  ct2 = 2.0D+00 * ( z - ab )
  sab = a + bp

  bb(1) = 1.0D+00
  aa(1) = 1.0D+00

  ct3 = sab + 1.0D+00 + ab
  bb(2) = 1.0D+00 + 2.0D+00 * z / ct3
  aa(2) = 1.0D+00 + ct2 / ct3

  anbn = ct3 + sab + 3.0D+00
  ct1 = 1.0D+00 + 2.0D+00 * z / anbn
  bb(3) = 1.0D+00 + 6.0D+00 * ct1 * z / ct3
  aa(3) = 1.0D+00 + 6.0D+00 * ab / anbn + 3.0D+00 * ct1 * ct2 / ct3

  do i = 4, 300

    x2i1 = real ( 2 * i - 3, kind = 8 )
    ct1 = x2i1 / ( x2i1 - 2.0D+00 )
    anbn = anbn + x2i1 + sab
    ct2 = ( x2i1 - 1.0D+00 ) /anbn
    c2 = x2i1 * ct2 - 1.0D+00
    d1z = x2i1 * 2.0D+00 * z / anbn

    ct3 = sab * ct2
    g1 = d1z + ct1 * ( c2 + ct3 )
    g2 = d1z - c2
    g3 = ct1 * ( 1.0D+00 - ct3 - 2.0D+00 * ct2 )

    bb(4) = g1 * bb(3) + g2 * bb(2) + g3 * bb(1)
    aa(4) = g1 * aa(3) + g2 * aa(2) + g3 * aa(1)

    r8_chu_scaled = aa(4) / bb(4)

    if ( abs ( r8_chu_scaled - aa(1) / bb(1) ) < &
      eps * abs ( r8_chu_scaled ) ) then
      return
    end if

    do j = 1, 3
      aa(j) = aa(j+1)
      bb(j) = bb(j+1)
    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_CHU_SCALED - Fatal error!'
  write ( *, '(a)' ) '  No convergence after 300 terms.'
  stop

end
function r8_ci ( x )

!*****************************************************************************80
!
!! R8_CI evaluates the cosine integral Ci of an R8 argument.
!
!  Discussion:
!
!    The cosine integral is defined by
!
!      CI(X) = - integral ( X <= T < Infinity ) ( cos ( T ) ) / T  dT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_CI, the cosine integral Ci evaluated at X.
!
  implicit none

  real ( kind = 8 ) cics(19)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) nci
  real ( kind = 8 ) r8_ci
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sinx
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save cics
  save nci
  save xsml

  data cics(  1) / -0.34004281856055363156281076633129873D+00 /
  data cics(  2) / -1.03302166401177456807159271040163751D+00 /
  data cics(  3) /  0.19388222659917082876715874606081709D+00 /
  data cics(  4) / -0.01918260436019865893946346270175301D+00 /
  data cics(  5) /  0.00110789252584784967184098099266118D+00 /
  data cics(  6) / -0.00004157234558247208803840231814601D+00 /
  data cics(  7) /  0.00000109278524300228715295578966285D+00 /
  data cics(  8) / -0.00000002123285954183465219601280329D+00 /
  data cics(  9) /  0.00000000031733482164348544865129873D+00 /
  data cics( 10) / -0.00000000000376141547987683699381798D+00 /
  data cics( 11) /  0.00000000000003622653488483964336956D+00 /
  data cics( 12) / -0.00000000000000028911528493651852433D+00 /
  data cics( 13) /  0.00000000000000000194327860676494420D+00 /
  data cics( 14) / -0.00000000000000000001115183182650184D+00 /
  data cics( 15) /  0.00000000000000000000005527858887706D+00 /
  data cics( 16) / -0.00000000000000000000000023907013943D+00 /
  data cics( 17) /  0.00000000000000000000000000091001612D+00 /
  data cics( 18) / -0.00000000000000000000000000000307233D+00 /
  data cics( 19) /  0.00000000000000000000000000000000926D+00 /

  data nci / 0 /
  data xsml / 0.0D+00 /

  if ( nci == 0 ) then
    nci = r8_inits ( cics, 19, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( r8_mach ( 3 ) )
  end if

  if ( x <= 0.0D+00 ) then
    
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CI - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.0.'
    stop
  
  else if ( x <= xsml ) then
    y = - 1.0D+00
    r8_ci = log ( x ) - 0.5D+00 + r8_csevl ( y, cics, nci )
  else if ( x <= 4.0D+00 ) then
    y = ( x * x - 8.0D+00 ) * 0.125D+00
    r8_ci = log ( x ) - 0.5D+00 + r8_csevl ( y, cics, nci )
  else
    call r8_sifg ( x, f, g )
    sinx = sin ( x )
    r8_ci = f * sinx - g * cos ( x )
  end if

  return
end
function r8_cin ( x )

!*****************************************************************************80
!
!! R8_CIN evaluates the alternate cosine integral Cin of an R8 argument.
!
!  Discussion:
!
!    CIN(X) = gamma + log(X) 
!      + integral ( 0 <= T <= X ) ( cos ( T ) - 1 ) / T  dT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_CIN, the cosine integral Cin evaluated at X.
!
  implicit none

  real ( kind = 8 ) absx
  real ( kind = 8 ) cincs(18)
  real ( kind = 8 ) eul
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) ncin
  real ( kind = 8 ) r8_cin
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sinx
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin

  save cincs
  save eul
  save ncin
  save xmin

  data cincs(  1) /  0.37074501750909688741654801228564992D+00 /
  data cincs(  2) / -0.05893574896364446831956864397363697D+00 /
  data cincs(  3) /  0.00538189642113569124048745326203340D+00 /
  data cincs(  4) / -0.00029860052841962135319594906563410D+00 /
  data cincs(  5) /  0.00001095572575321620077031054467306D+00 /
  data cincs(  6) / -0.00000028405454877346630491727187731D+00 /
  data cincs(  7) /  0.00000000546973994875384912457861806D+00 /
  data cincs(  8) / -0.00000000008124187461318157083277452D+00 /
  data cincs(  9) /  0.00000000000095868593117706609013181D+00 /
  data cincs( 10) / -0.00000000000000920266004392351031377D+00 /
  data cincs( 11) /  0.00000000000000007325887999017895024D+00 /
  data cincs( 12) / -0.00000000000000000049143726675842909D+00 /
  data cincs( 13) /  0.00000000000000000000281577746753902D+00 /
  data cincs( 14) / -0.00000000000000000000001393986788501D+00 /
  data cincs( 15) /  0.00000000000000000000000006022485646D+00 /
  data cincs( 16) / -0.00000000000000000000000000022904717D+00 /
  data cincs( 17) /  0.00000000000000000000000000000077273D+00 /
  data cincs( 18) / -0.00000000000000000000000000000000233D+00 /

  data eul / 0.57721566490153286060651209008240D+00 /
  data ncin / 0 /
  data xmin / 0.0D+00 /

  if ( ncin == 0 ) then
    ncin = r8_inits ( cincs, 18, 0.1D+00 * r8_mach ( 3 ) )
    xmin = sqrt ( r8_mach ( 1 ) )
  end if

  absx = abs ( x )

  if ( absx <= xmin ) then
    r8_cin = 0.0D+00
  else if ( absx <= 4.0D+00 ) then
    r8_cin = r8_csevl ( ( x * x - 8.0D+00 ) * 0.125D+00, cincs, &
      ncin ) * x * x
  else
    call r8_sifg ( absx, f, g )
    sinx = sin ( absx )
    r8_cin = - f * sinx + g * cos ( absx ) + log ( absx ) + eul
  end if

  return
end
function r8_cinh ( x )

!*****************************************************************************80
!
!! R8_CINH: alternate hyperbolic cosine integral Cinh of an R8 argument.
!
!  Discussion:
!
!    Cinh ( x ) = Integral ( 0 <= t <= x ) ( cosh ( t ) - 1 ) dt / t
!
!    The original text of this program had a mistake:
!      y = x * x / 9.0D+00 - 1.0D+00
!    has been corrected to
!      y = x * x / 4.5D+00 - 1.0D+00
!    JVB, 27 March 2010
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_CINH, the hyperbolic cosine integral Cinh
!    evaluated at X.
!
  implicit none

  real ( kind = 8 ) absx
  real ( kind = 8 ) cinhcs(16)
  real ( kind = 8 ) eul
  integer ( kind = 4 ) ncinh
  real ( kind = 8 ) r8_chi
  real ( kind = 8 ) r8_cinh
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save cinhcs
  save eul
  save ncinh
  save xmin
  save xsml

  data cinhcs(  1) /     0.1093291636520734431407425199795917D+00 /
  data cinhcs(  2) /     0.0573928847550379676445323429825108D+00 /
  data cinhcs(  3) /     0.0028095756978830353416404208940774D+00 /
  data cinhcs(  4) /     0.0000828780840721356655731765069792D+00 /
  data cinhcs(  5) /     0.0000016278596173914185577726018815D+00 /
  data cinhcs(  6) /     0.0000000227809519255856619859083591D+00 /
  data cinhcs(  7) /     0.0000000002384484842463059257284002D+00 /
  data cinhcs(  8) /     0.0000000000019360829780781957471028D+00 /
  data cinhcs(  9) /     0.0000000000000125453698328172559683D+00 /
  data cinhcs( 10) /     0.0000000000000000663637449497262300D+00 /
  data cinhcs( 11) /     0.0000000000000000002919639263594744D+00 /
  data cinhcs( 12) /     0.0000000000000000000010849123956107D+00 /
  data cinhcs( 13) /     0.0000000000000000000000034499080805D+00 /
  data cinhcs( 14) /     0.0000000000000000000000000094936664D+00 /
  data cinhcs( 15) /     0.0000000000000000000000000000228291D+00 /
  data cinhcs( 16) /     0.0000000000000000000000000000000484D+00 /

  data eul / 0.57721566490153286060651209008240D+00  /
  data ncinh / 0 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ncinh == 0 ) then
    ncinh = r8_inits ( cinhcs, 16, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( r8_mach ( 3 ) )
    xmin = 2.0D+00 * sqrt ( r8_mach ( 1 ) )
  end if

  absx = abs ( x )

  if ( x == 0.0D+00 ) then
    r8_cinh = 0.0D+00
  else if ( absx <= xmin ) then
    r8_cinh = 0.0D+00
  else if ( x <= xsml ) then
    y = - 1.0D+00
    r8_cinh = x * x * ( 0.25D+00 + r8_csevl ( y, cinhcs, ncinh ) )
  else if ( x <= 3.0D+00 ) then
    y = x * x / 4.5D+00 - 1.0D+00
    r8_cinh = x * x * ( 0.25D+00 + r8_csevl ( y, cinhcs, ncinh ) )
  else
    r8_cinh = r8_chi ( absx ) - eul - log ( absx )
  end if

  return
end
function r8_cos ( x )

!*****************************************************************************80
!
!! R8_COS evaluates the cosine of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_COS, the cosine of X.
!
  implicit none

  real ( kind = 8 ) absx
  real ( kind = 8 ) f
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) ntsn
  real ( kind = 8 ) pi2
  real ( kind = 8 ) pi2rec
  real ( kind = 8 ) pihi
  real ( kind = 8 ) pilo
  real ( kind = 8 ) pirec
  real ( kind = 8 ) r8_cos
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sincs(15)
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xn
  real ( kind = 8 ) xsml
  real ( kind = 8 ) xwarn
  real ( kind = 8 ) y

  save ntsn
  save pi2
  save pi2rec
  save pihi
  save pilo
  save pirec
  save sincs
  save xmax
  save xsml
  save xwarn

  data sincs(  1) / -0.374991154955873175839919279977323464D+00 /
  data sincs(  2) / -0.181603155237250201863830316158004754D+00 /
  data sincs(  3) / +0.005804709274598633559427341722857921D+00 /
  data sincs(  4) / -0.000086954311779340757113212316353178D+00 /
  data sincs(  5) / +0.000000754370148088851481006839927030D+00 /
  data sincs(  6) / -0.000000004267129665055961107126829906D+00 /
  data sincs(  7) / +0.000000000016980422945488168181824792D+00 /
  data sincs(  8) / -0.000000000000050120578889961870929524D+00 /
  data sincs(  9) / +0.000000000000000114101026680010675628D+00 /
  data sincs( 10) / -0.000000000000000000206437504424783134D+00 /
  data sincs( 11) / +0.000000000000000000000303969595918706D+00 /
  data sincs( 12) / -0.000000000000000000000000371357734157D+00 /
  data sincs( 13) / +0.000000000000000000000000000382486123D+00 /
  data sincs( 14) / -0.000000000000000000000000000000336623D+00 /
  data sincs( 15) / +0.000000000000000000000000000000000256D+00 /
!
!  pihi + pilo = pi.  pihi is exactly representable on all machines
!  with at least 8 bits of precision.  whether it is exactly
!  represented depends on the compiler.  this routine is more
!  accurate if it is exactly represented.
!
  data ntsn / 0 /
  data pi2 / 1.57079632679489661923132169163975D+00 /
  data pi2rec / 0.63661977236758134307553505349006D+00 /
  data pihi / 3.140625D+00 /
  data pilo / 9.6765358979323846264338327950288D-04 /
  data pirec / 0.31830988618379067153776752674503D+00 /
  data xmax / 0.0D+00 /
  data xsml / 0.0D+00 /
  data xwarn / 0.0D+00 /

  if ( ntsn == 0 ) then
    ntsn = r8_inits ( sincs, 15, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( 2.0D+00 * r8_mach ( 3 ) )
    xmax = 1.0D+00 / r8_mach ( 4 )
    xwarn = sqrt ( xmax )
  end if

  absx = abs ( x )
  y = absx + pi2

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_COS - Warning!'
    write ( *, '(a)' ) '  No precision because |X| is big.'
    r8_cos = 0.0D+00
    return
  end if

  if ( xwarn < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_COS - Warning!'
    write ( *, '(a)' ) '  Answer < half precision because |X| is big.'
  end if

  r8_cos = 1.0D+00

  if ( absx < xsml ) then
    return
  end if

  xn = aint ( y * pirec + 0.5D+00 )
  n2 = int ( mod ( xn, 2.0D+00 ) + 0.5D+00 )
  xn = xn - 0.5D+00
  f = ( absx - xn * pihi ) - xn * pilo

  xn = 2.0D+00 * ( f * pi2rec ) * ( f * pi2rec ) - 1.0D+00

  r8_cos = f + f * r8_csevl ( xn, sincs, ntsn )

  if ( n2 /= 0 ) then
    r8_cos = - r8_cos
  end if

  if ( r8_cos < -1.0D+00 ) then
    r8_cos = -1.0D+00
  else if ( 1.0D+00 < r8_cos ) then
    r8_cos = 1.0D+00
  end if

  return
end
function r8_cos_deg ( x )

!*****************************************************************************80
!
!! R8_COS_DEG evaluates the cosine of an R8 argument in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument in degrees.
!
!    Output, real ( kind = 8 ) R8_COS_DEG, the cosine of X.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_cos_deg
  real ( kind = 8 ), parameter :: raddeg &
    = 0.017453292519943295769236907684886D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = cos ( raddeg * x )

  if ( mod ( x, 90.0D+00 ) == 0.0D+00 ) then

    n = int ( abs ( x ) / 90.0D+00 + 0.5D+00 )
    n = mod ( n, 2 )

    if ( n == 1 ) then
      value = 0.0D+00
    else if ( value < 0.0D+00 ) then
      value = - 1.0D+00
    else
      value = + 1.0D+00
    end if

  end if

  r8_cos_deg = value

  return
end
function r8_cosh ( x )

!*****************************************************************************80
!
!! R8_COSH evaluates the hyperbolic cosine of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_COSH, the hyperbolic cosine of X.
!
  implicit none

  real ( kind = 8 ) r8_cosh
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax

  save ymax

  data ymax / 0.0D+00 /

  if ( ymax == 0.0D+00 ) then
    ymax = 1.0D+00 / sqrt ( r8_mach ( 3 ) )
  end if

  y = exp ( abs ( x ) )

  if ( y < ymax ) then
    value = 0.5D+00 * ( y + 1.0D+00 / y )
  else
    value = 0.5D+00 * y
  end if

  r8_cosh = value

  return
end
function r8_cot ( x )

!*****************************************************************************80
!
!! R8_COT evaluates the cotangent of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_COT, the cotangent of X.
!
  implicit none

  real ( kind = 8 ) ainty
  real ( kind = 8 ) ainty2
  real ( kind = 8 ) cotcs(15)
  integer ( kind = 4 ) ifn
  integer ( kind = 4 ) nterms
  real ( kind = 8 ) pi2rec
  real ( kind = 8 ) prodbg
  real ( kind = 8 ) r8_cot
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y
  real ( kind = 8 ) yrem

  save cotcs
  save nterms
  save pi2rec
  save xmax
  save xmin
  save xsml

  data cotcs(  1) / +0.240259160982956302509553617744970D+00 /
  data cotcs(  2) / -0.165330316015002278454746025255758D-01 /
  data cotcs(  3) / -0.429983919317240189356476228239895D-04 /
  data cotcs(  4) / -0.159283223327541046023490851122445D-06 /
  data cotcs(  5) / -0.619109313512934872588620579343187D-09 /
  data cotcs(  6) / -0.243019741507264604331702590579575D-11 /
  data cotcs(  7) / -0.956093675880008098427062083100000D-14 /
  data cotcs(  8) / -0.376353798194580580416291539706666D-16 /
  data cotcs(  9) / -0.148166574646746578852176794666666D-18 /
  data cotcs( 10) / -0.583335658903666579477984000000000D-21 /
  data cotcs( 11) / -0.229662646964645773928533333333333D-23 /
  data cotcs( 12) / -0.904197057307483326719999999999999D-26 /
  data cotcs( 13) / -0.355988551920600064000000000000000D-28 /
  data cotcs( 14) / -0.140155139824298666666666666666666D-30 /
  data cotcs( 15) / -0.551800436872533333333333333333333D-33 /

  data nterms / 0 /
  data pi2rec / 0.011619772367581343075535053490057D+00 /
  data xmax / 0.0D+00 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( nterms == 0 ) then
    nterms = r8_inits ( cotcs, 15, 0.1D+00 * r8_mach ( 3 ) )
    xmax = 1.0D+00 / r8_mach ( 4 )
    xsml = sqrt ( 3.0D+00 * r8_mach ( 3 ) )
    xmin = exp ( max ( log ( r8_mach ( 1 ) ), &
      - log ( r8_mach ( 2 ) ) )  + 0.01D+00 )
    sqeps = sqrt ( r8_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( y < xmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_COT - Fatal error!'
    write ( *, '(a)' ) '  |X| is too small.'
    stop
  end if

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_COT - Fatal error!'
    write ( *, '(a)' ) '  |X| is too big.'
    stop
  end if
!
!  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
!  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
!  = aint(.625*y) + aint(z) + rem(z)
!
  ainty = aint ( y )
  yrem = y - ainty
  prodbg = 0.625D+00 * ainty
  ainty = aint ( prodbg )
  y = ( prodbg - ainty ) + 0.625D+00 * yrem + y * pi2rec
  ainty2 = aint ( y )
  ainty = ainty + ainty2
  y = y - ainty2

  ifn = int ( mod ( ainty, 2.0D+00 ) )

  if ( ifn == 1 ) then
    y = 1.0D+00 - y
  end if

  if ( 0.5D+00 < abs ( x ) .and. y < abs ( x ) * sqeps ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_COT - Warning!'
    write ( *, '(a)' ) '  Answer less than half precision.'
    write ( *, '(a)' ) '  |X| too big, or X nearly a nonzero multiple of pi.'
  end if

  if ( y == 0.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_COT - Fatal error!'
    write ( *, '(a)' ) '  X is a multiple of pi.'
    stop

  else if ( y <= xsml ) then

    r8_cot = 1.0D+00 / y 

  else if ( y <= 0.25D+00 ) then

    r8_cot = ( 0.5D+00 &
      + r8_csevl ( 32.0D+00 * y * y - 1.0D+00, cotcs, nterms ) ) / y

  else if ( y <= 0.5D+00 ) then

    r8_cot = ( 0.5D+00 + r8_csevl ( 8.0D+00 * y * y - 1.0D+00, &
      cotcs, nterms ) ) / ( 0.5D+00 * y )

    r8_cot = ( r8_cot * r8_cot - 1.0D+00 ) * 0.5D+00 / r8_cot

  else

    r8_cot = ( 0.5D+00 + r8_csevl ( 2.0D+00 * y * y - 1.0D+00, &
      cotcs, nterms ) ) / ( 0.25D+00 * y )
    r8_cot = ( r8_cot * r8_cot - 1.0D+00 ) * 0.5D+00 / r8_cot
    r8_cot = ( r8_cot * r8_cot - 1.0D+00 ) * 0.5D+00 / r8_cot

  end if

  if ( x < 0.0D+00 ) then
    r8_cot = - abs ( r8_cot )
  else
    r8_cot = + abs ( r8_cot )
  end if

  if ( ifn == 1 ) then
    r8_cot = - r8_cot
  end if

  return
end
function r8_csevl ( x, a, n )

!*****************************************************************************80
!
!! R8_CSEVL evaluates a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Input, real ( kind = 8 ) CS(N), the Chebyshev coefficients.
!
!    Input, integer ( kind = 4 ) N, the number of Chebyshev coefficients.
!
!    Output, real ( kind = 8 ) R8_CSEVL, the Chebyshev series evaluated at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) twox
  real ( kind = 8 ) x

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms <= 0.'
    stop
  end if

  if ( 1000 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms > 1000.'
    stop
  end if

  if ( x < -1.1D+00 .or. 1.1D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
    write ( *, '(a)' ) '  X outside (-1,+1)'
    write ( *, '(a,g14.6)' ) '  X = ', x
    stop
  end if

  twox = 2.0D+00 * x
  b1 = 0.0D+00
  b0 = 0.0D+00

  do i = n, 1, -1
    b2 = b1
    b1 = b0
    b0 = twox * b1 - b2 + a(i)
  end do

  r8_csevl = 0.5D+00 * ( b0 - b2 )

  return
end
function r8_dawson ( x )

!*****************************************************************************80
!
!! R8_DAWSON evaluates Dawson's integral of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_DAWSON, the value of Dawson's integral at X.
!
  implicit none

  real ( kind = 8 ) daw2cs(45)
  real ( kind = 8 ) dawacs(75)
  real ( kind = 8 ) dawcs(21)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) ntdaw
  integer ( kind = 4 ) ntdaw2
  integer ( kind = 4 ) ntdawa
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_dawson
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save daw2cs
  save dawacs
  save dawcs
  save ntdaw
  save ntdaw2
  save ntdawa
  save xbig
  save xmax
  save xsml

  data dawcs( 1) / -0.6351734375145949201065127736293D-02 /
  data dawcs( 2) / -0.2294071479677386939899824125866D+00 /
  data dawcs( 3) / +0.2213050093908476441683979161786D-01 /
  data dawcs( 4) / -0.1549265453892985046743057753375D-02 /
  data dawcs( 5) / +0.8497327715684917456777542948066D-04 /
  data dawcs( 6) / -0.3828266270972014924994099521309D-05 /
  data dawcs( 7) / +0.1462854806250163197757148949539D-06 /
  data dawcs( 8) / -0.4851982381825991798846715425114D-08 /
  data dawcs( 9) / +0.1421463577759139790347568183304D-09 /
  data dawcs(10) / -0.3728836087920596525335493054088D-11 /
  data dawcs(11) / +0.8854942961778203370194565231369D-13 /
  data dawcs(12) / -0.1920757131350206355421648417493D-14 /
  data dawcs(13) / +0.3834325867246327588241074439253D-16 /
  data dawcs(14) / -0.7089154168175881633584099327999D-18 /
  data dawcs(15) / +0.1220552135889457674416901120000D-19 /
  data dawcs(16) / -0.1966204826605348760299451733333D-21 /
  data dawcs(17) / +0.2975845541376597189113173333333D-23 /
  data dawcs(18) / -0.4247069514800596951039999999999D-25 /
  data dawcs(19) / +0.5734270767391742798506666666666D-27 /
  data dawcs(20) / -0.7345836823178450261333333333333D-29 /
  data dawcs(21) / +0.8951937667516552533333333333333D-31 /

  data daw2cs( 1) / -0.56886544105215527114160533733674D-01 /
  data daw2cs( 2) / -0.31811346996168131279322878048822D+00 /
  data daw2cs( 3) / +0.20873845413642236789741580198858D+00 /
  data daw2cs( 4) / -0.12475409913779131214073498314784D+00 /
  data daw2cs( 5) / +0.67869305186676777092847516423676D-01 /
  data daw2cs( 6) / -0.33659144895270939503068230966587D-01 /
  data daw2cs( 7) / +0.15260781271987971743682460381640D-01 /
  data daw2cs( 8) / -0.63483709625962148230586094788535D-02 /
  data daw2cs( 9) / +0.24326740920748520596865966109343D-02 /
  data daw2cs(10) / -0.86219541491065032038526983549637D-03 /
  data daw2cs(11) / +0.28376573336321625302857636538295D-03 /
  data daw2cs(12) / -0.87057549874170423699396581464335D-04 /
  data daw2cs(13) / +0.24986849985481658331800044137276D-04 /
  data daw2cs(14) / -0.67319286764160294344603050339520D-05 /
  data daw2cs(15) / +0.17078578785573543710504524047844D-05 /
  data daw2cs(16) / -0.40917551226475381271896592490038D-06 /
  data daw2cs(17) / +0.92828292216755773260751785312273D-07 /
  data daw2cs(18) / -0.19991403610147617829845096332198D-07 /
  data daw2cs(19) / +0.40963490644082195241210487868917D-08 /
  data daw2cs(20) / -0.80032409540993168075706781753561D-09 /
  data daw2cs(21) / +0.14938503128761465059143225550110D-09 /
  data daw2cs(22) / -0.26687999885622329284924651063339D-10 /
  data daw2cs(23) / +0.45712216985159458151405617724103D-11 /
  data daw2cs(24) / -0.75187305222043565872243727326771D-12 /
  data daw2cs(25) / +0.11893100052629681879029828987302D-12 /
  data daw2cs(26) / -0.18116907933852346973490318263084D-13 /
  data daw2cs(27) / +0.26611733684358969193001612199626D-14 /
  data daw2cs(28) / -0.37738863052129419795444109905930D-15 /
  data daw2cs(29) / +0.51727953789087172679680082229329D-16 /
  data daw2cs(30) / -0.68603684084077500979419564670102D-17 /
  data daw2cs(31) / +0.88123751354161071806469337321745D-18 /
  data daw2cs(32) / -0.10974248249996606292106299624652D-18 /
  data daw2cs(33) / +0.13261199326367178513595545891635D-19 /
  data daw2cs(34) / -0.15562732768137380785488776571562D-20 /
  data daw2cs(35) / +0.17751425583655720607833415570773D-21 /
  data daw2cs(36) / -0.19695006967006578384953608765439D-22 /
  data daw2cs(37) / +0.21270074896998699661924010120533D-23 /
  data daw2cs(38) / -0.22375398124627973794182113962666D-24 /
  data daw2cs(39) / +0.22942768578582348946971383125333D-25 /
  data daw2cs(40) / -0.22943788846552928693329592319999D-26 /
  data daw2cs(41) / +0.22391702100592453618342297600000D-27 /
  data daw2cs(42) / -0.21338230616608897703678225066666D-28 /
  data daw2cs(43) / +0.19866196585123531518028458666666D-29 /
  data daw2cs(44) / -0.18079295866694391771955199999999D-30 /
  data daw2cs(45) / +0.16090686015283030305450666666666D-31 /

  data dawacs( 1) / +0.1690485637765703755422637438849D-01 /
  data dawacs( 2) / +0.8683252278406957990536107850768D-02 /
  data dawacs( 3) / +0.2424864042417715453277703459889D-03 /
  data dawacs( 4) / +0.1261182399572690001651949240377D-04 /
  data dawacs( 5) / +0.1066453314636176955705691125906D-05 /
  data dawacs( 6) / +0.1358159794790727611348424505728D-06 /
  data dawacs( 7) / +0.2171042356577298398904312744743D-07 /
  data dawacs( 8) / +0.2867010501805295270343676804813D-08 /
  data dawacs( 9) / -0.1901336393035820112282492378024D-09 /
  data dawacs(10) / -0.3097780484395201125532065774268D-09 /
  data dawacs(11) / -0.1029414876057509247398132286413D-09 /
  data dawacs(12) / -0.6260356459459576150417587283121D-11 /
  data dawacs(13) / +0.8563132497446451216262303166276D-11 /
  data dawacs(14) / +0.3033045148075659292976266276257D-11 /
  data dawacs(15) / -0.2523618306809291372630886938826D-12 /
  data dawacs(16) / -0.4210604795440664513175461934510D-12 /
  data dawacs(17) / -0.4431140826646238312143429452036D-13 /
  data dawacs(18) / +0.4911210272841205205940037065117D-13 /
  data dawacs(19) / +0.1235856242283903407076477954739D-13 /
  data dawacs(20) / -0.5788733199016569246955765071069D-14 /
  data dawacs(21) / -0.2282723294807358620978183957030D-14 /
  data dawacs(22) / +0.7637149411014126476312362917590D-15 /
  data dawacs(23) / +0.3851546883566811728777594002095D-15 /
  data dawacs(24) / -0.1199932056928290592803237283045D-15 /
  data dawacs(25) / -0.6313439150094572347334270285250D-16 /
  data dawacs(26) / +0.2239559965972975375254912790237D-16 /
  data dawacs(27) / +0.9987925830076495995132891200749D-17 /
  data dawacs(28) / -0.4681068274322495334536246507252D-17 /
  data dawacs(29) / -0.1436303644349721337241628751534D-17 /
  data dawacs(30) / +0.1020822731410541112977908032130D-17 /
  data dawacs(31) / +0.1538908873136092072837389822372D-18 /
  data dawacs(32) / -0.2189157877645793888894790926056D-18 /
  data dawacs(33) / +0.2156879197938651750392359152517D-20 /
  data dawacs(34) / +0.4370219827442449851134792557395D-19 /
  data dawacs(35) / -0.8234581460977207241098927905177D-20 /
  data dawacs(36) / -0.7498648721256466222903202835420D-20 /
  data dawacs(37) / +0.3282536720735671610957612930039D-20 /
  data dawacs(38) / +0.8858064309503921116076561515151D-21 /
  data dawacs(39) / -0.9185087111727002988094460531485D-21 /
  data dawacs(40) / +0.2978962223788748988314166045791D-22 /
  data dawacs(41) / +0.1972132136618471883159505468041D-21 /
  data dawacs(42) / -0.5974775596362906638089584995117D-22 /
  data dawacs(43) / -0.2834410031503850965443825182441D-22 /
  data dawacs(44) / +0.2209560791131554514777150489012D-22 /
  data dawacs(45) / -0.5439955741897144300079480307711D-25 /
  data dawacs(46) / -0.5213549243294848668017136696470D-23 /
  data dawacs(47) / +0.1702350556813114199065671499076D-23 /
  data dawacs(48) / +0.6917400860836148343022185660197D-24 /
  data dawacs(49) / -0.6540941793002752512239445125802D-24 /
  data dawacs(50) / +0.6093576580439328960371824654636D-25 /
  data dawacs(51) / +0.1408070432905187461501945080272D-24 /
  data dawacs(52) / -0.6785886121054846331167674943755D-25 /
  data dawacs(53) / -0.9799732036214295711741583102225D-26 /
  data dawacs(54) / +0.2121244903099041332598960939160D-25 /
  data dawacs(55) / -0.5954455022548790938238802154487D-26 /
  data dawacs(56) / -0.3093088861875470177838847232049D-26 /
  data dawacs(57) / +0.2854389216344524682400691986104D-26 /
  data dawacs(58) / -0.3951289447379305566023477271811D-27 /
  data dawacs(59) / -0.5906000648607628478116840894453D-27 /
  data dawacs(60) / +0.3670236964668687003647889980609D-27 /
  data dawacs(61) / -0.4839958238042276256598303038941D-29 /
  data dawacs(62) / -0.9799265984210443869597404017022D-28 /
  data dawacs(63) / +0.4684773732612130606158908804300D-28 /
  data dawacs(64) / +0.5030877696993461051647667603155D-29 /
  data dawacs(65) / -0.1547395051706028239247552068295D-28 /
  data dawacs(66) / +0.6112180185086419243976005662714D-29 /
  data dawacs(67) / +0.1357913399124811650343602736158D-29 /
  data dawacs(68) / -0.2417687752768673088385304299044D-29 /
  data dawacs(69) / +0.8369074582074298945292887587291D-30 /
  data dawacs(70) / +0.2665413042788979165838319401566D-30 /
  data dawacs(71) / -0.3811653692354890336935691003712D-30 /
  data dawacs(72) / +0.1230054721884951464371706872585D-30 /
  data dawacs(73) / +0.4622506399041493508805536929983D-31 /
  data dawacs(74) / -0.6120087296881677722911435593001D-31 /
  data dawacs(75) / +0.1966024640193164686956230217896D-31 /

  data ntdaw / 0 /
  data ntdaw2 / 0 /
  data ntdawa / 0 /
  data xbig / 0.0D+00 /
  data xmax / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ntdaw == 0 ) then

    eps = r8_mach ( 3 )
    ntdaw  = r8_inits ( dawcs,  21, 0.1D+00 * eps )
    ntdaw2 = r8_inits ( daw2cs, 45, 0.1D+00 * eps )
    ntdawa = r8_inits ( dawacs, 75, 0.1D+00 * eps )
    xsml = sqrt ( 1.5D+00 * eps )
    xbig = sqrt ( 0.5D+00 / eps )
    xmax = exp ( min ( - log ( 2.0D+00 * r8_mach ( 1 ) ), &
      log ( r8_mach ( 2 ) ) ) - 0.01D+00 )

  end if

  y = abs ( x )

  if ( y <= xsml ) then
    r8_dawson = x
  else if ( y <= 1.0D+00 ) then
    r8_dawson = x * ( 0.75D+00 &
      + r8_csevl ( 2.0D+00 * y * y - 1.0D+00, dawcs, ntdaw ) )
  else if ( y <= 4.0D+00 ) then
    r8_dawson = x * ( 0.25D+00 &
      + r8_csevl ( 0.125D+00 * y * y - 1.0D+00, daw2cs, ntdaw2 ) )
  else if ( y < xbig ) then
    r8_dawson = ( 0.5D+00 &
      + r8_csevl ( 32.0D+00 / y / y - 1.0D+00, dawacs, ntdawa ) ) / x
  else if ( y <= xmax ) then
    r8_dawson = 0.5D+00 / x
  else
    r8_dawson = 0.0D+00
  end if

  return
end
function r8_e1 ( x )

!*****************************************************************************80
!
!! R8_E1 evaluates the exponential integral E1 for an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_E1, the exponential integral E1.
!
  implicit none

  real ( kind = 8 ) ae10cs(50)
  real ( kind = 8 ) ae11cs(60)
  real ( kind = 8 ) ae12cs(41)
  real ( kind = 8 ) ae13cs(50)
  real ( kind = 8 ) ae14cs(64)
  real ( kind = 8 ) e11cs(29)
  real ( kind = 8 ) e12cs(25)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) ntae10
  integer ( kind = 4 ) ntae11
  integer ( kind = 4 ) ntae12
  integer ( kind = 4 ) ntae13
  integer ( kind = 4 ) ntae14
  integer ( kind = 4 ) nte11
  integer ( kind = 4 ) nte12
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_e1
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax

  save ae10cs
  save ae11cs
  save ae12cs
  save ae13cs
  save ae14cs
  save e11cs
  save e12cs
  save ntae10
  save ntae11
  save ntae12
  save ntae13
  save ntae14
  save nte11
  save nte12
  save xmax

  data ae10cs( 1) / +0.3284394579616699087873844201881D-01 /
  data ae10cs( 2) / -0.1669920452031362851476184343387D-01 /
  data ae10cs( 3) / +0.2845284724361346807424899853252D-03 /
  data ae10cs( 4) / -0.7563944358516206489487866938533D-05 /
  data ae10cs( 5) / +0.2798971289450859157504843180879D-06 /
  data ae10cs( 6) / -0.1357901828534531069525563926255D-07 /
  data ae10cs( 7) / +0.8343596202040469255856102904906D-09 /
  data ae10cs( 8) / -0.6370971727640248438275242988532D-10 /
  data ae10cs( 9) / +0.6007247608811861235760831561584D-11 /
  data ae10cs(10) / -0.7022876174679773590750626150088D-12 /
  data ae10cs(11) / +0.1018302673703687693096652346883D-12 /
  data ae10cs(12) / -0.1761812903430880040406309966422D-13 /
  data ae10cs(13) / +0.3250828614235360694244030353877D-14 /
  data ae10cs(14) / -0.5071770025505818678824872259044D-15 /
  data ae10cs(15) / +0.1665177387043294298172486084156D-16 /
  data ae10cs(16) / +0.3166753890797514400677003536555D-16 /
  data ae10cs(17) / -0.1588403763664141515133118343538D-16 /
  data ae10cs(18) / +0.4175513256138018833003034618484D-17 /
  data ae10cs(19) / -0.2892347749707141906710714478852D-18 /
  data ae10cs(20) / -0.2800625903396608103506340589669D-18 /
  data ae10cs(21) / +0.1322938639539270903707580023781D-18 /
  data ae10cs(22) / -0.1804447444177301627283887833557D-19 /
  data ae10cs(23) / -0.7905384086522616076291644817604D-20 /
  data ae10cs(24) / +0.4435711366369570103946235838027D-20 /
  data ae10cs(25) / -0.4264103994978120868865309206555D-21 /
  data ae10cs(26) / -0.3920101766937117541553713162048D-21 /
  data ae10cs(27) / +0.1527378051343994266343752326971D-21 /
  data ae10cs(28) / +0.1024849527049372339310308783117D-22 /
  data ae10cs(29) / -0.2134907874771433576262711405882D-22 /
  data ae10cs(30) / +0.3239139475160028267061694700366D-23 /
  data ae10cs(31) / +0.2142183762299889954762643168296D-23 /
  data ae10cs(32) / -0.8234609419601018414700348082312D-24 /
  data ae10cs(33) / -0.1524652829645809479613694401140D-24 /
  data ae10cs(34) / +0.1378208282460639134668480364325D-24 /
  data ae10cs(35) / +0.2131311202833947879523224999253D-26 /
  data ae10cs(36) / -0.2012649651526484121817466763127D-25 /
  data ae10cs(37) / +0.1995535662263358016106311782673D-26 /
  data ae10cs(38) / +0.2798995808984003464948686520319D-26 /
  data ae10cs(39) / -0.5534511845389626637640819277823D-27 /
  data ae10cs(40) / -0.3884995396159968861682544026146D-27 /
  data ae10cs(41) / +0.1121304434507359382850680354679D-27 /
  data ae10cs(42) / +0.5566568152423740948256563833514D-28 /
  data ae10cs(43) / -0.2045482929810499700448533938176D-28 /
  data ae10cs(44) / -0.8453813992712336233411457493674D-29 /
  data ae10cs(45) / +0.3565758433431291562816111116287D-29 /
  data ae10cs(46) / +0.1383653872125634705539949098871D-29 /
  data ae10cs(47) / -0.6062167864451372436584533764778D-30 /
  data ae10cs(48) / -0.2447198043989313267437655119189D-30 /
  data ae10cs(49) / +0.1006850640933998348011548180480D-30 /
  data ae10cs(50) / +0.4623685555014869015664341461674D-31 /

  data ae11cs( 1) / +0.20263150647078889499401236517381D+00 /
  data ae11cs( 2) / -0.73655140991203130439536898728034D-01 /
  data ae11cs( 3) / +0.63909349118361915862753283840020D-02 /
  data ae11cs( 4) / -0.60797252705247911780653153363999D-03 /
  data ae11cs( 5) / -0.73706498620176629330681411493484D-04 /
  data ae11cs( 6) / +0.48732857449450183453464992488076D-04 /
  data ae11cs( 7) / -0.23837064840448290766588489460235D-05 /
  data ae11cs( 8) / -0.30518612628561521027027332246121D-05 /
  data ae11cs( 9) / +0.17050331572564559009688032992907D-06 /
  data ae11cs(10) / +0.23834204527487747258601598136403D-06 /
  data ae11cs(11) / +0.10781772556163166562596872364020D-07 /
  data ae11cs(12) / -0.17955692847399102653642691446599D-07 /
  data ae11cs(13) / -0.41284072341950457727912394640436D-08 /
  data ae11cs(14) / +0.68622148588631968618346844526664D-09 /
  data ae11cs(15) / +0.53130183120506356147602009675961D-09 /
  data ae11cs(16) / +0.78796880261490694831305022893515D-10 /
  data ae11cs(17) / -0.26261762329356522290341675271232D-10 /
  data ae11cs(18) / -0.15483687636308261963125756294100D-10 /
  data ae11cs(19) / -0.25818962377261390492802405122591D-11 /
  data ae11cs(20) / +0.59542879191591072658903529959352D-12 /
  data ae11cs(21) / +0.46451400387681525833784919321405D-12 /
  data ae11cs(22) / +0.11557855023255861496288006203731D-12 /
  data ae11cs(23) / -0.10475236870835799012317547189670D-14 /
  data ae11cs(24) / -0.11896653502709004368104489260929D-13 /
  data ae11cs(25) / -0.47749077490261778752643019349950D-14 /
  data ae11cs(26) / -0.81077649615772777976249734754135D-15 /
  data ae11cs(27) / +0.13435569250031554199376987998178D-15 /
  data ae11cs(28) / +0.14134530022913106260248873881287D-15 /
  data ae11cs(29) / +0.49451592573953173115520663232883D-16 /
  data ae11cs(30) / +0.79884048480080665648858587399367D-17 /
  data ae11cs(31) / -0.14008632188089809829248711935393D-17 /
  data ae11cs(32) / -0.14814246958417372107722804001680D-17 /
  data ae11cs(33) / -0.55826173646025601904010693937113D-18 /
  data ae11cs(34) / -0.11442074542191647264783072544598D-18 /
  data ae11cs(35) / +0.25371823879566853500524018479923D-20 /
  data ae11cs(36) / +0.13205328154805359813278863389097D-19 /
  data ae11cs(37) / +0.62930261081586809166287426789485D-20 /
  data ae11cs(38) / +0.17688270424882713734999261332548D-20 /
  data ae11cs(39) / +0.23266187985146045209674296887432D-21 /
  data ae11cs(40) / -0.67803060811125233043773831844113D-22 /
  data ae11cs(41) / -0.59440876959676373802874150531891D-22 /
  data ae11cs(42) / -0.23618214531184415968532592503466D-22 /
  data ae11cs(43) / -0.60214499724601478214168478744576D-23 /
  data ae11cs(44) / -0.65517906474348299071370444144639D-24 /
  data ae11cs(45) / +0.29388755297497724587042038699349D-24 /
  data ae11cs(46) / +0.22601606200642115173215728758510D-24 /
  data ae11cs(47) / +0.89534369245958628745091206873087D-25 /
  data ae11cs(48) / +0.24015923471098457555772067457706D-25 /
  data ae11cs(49) / +0.34118376888907172955666423043413D-26 /
  data ae11cs(50) / -0.71617071694630342052355013345279D-27 /
  data ae11cs(51) / -0.75620390659281725157928651980799D-27 /
  data ae11cs(52) / -0.33774612157467324637952920780800D-27 /
  data ae11cs(53) / -0.10479325703300941711526430332245D-27 /
  data ae11cs(54) / -0.21654550252170342240854880201386D-28 /
  data ae11cs(55) / -0.75297125745288269994689298432000D-30 /
  data ae11cs(56) / +0.19103179392798935768638084000426D-29 /
  data ae11cs(57) / +0.11492104966530338547790728833706D-29 /
  data ae11cs(58) / +0.43896970582661751514410359193600D-30 /
  data ae11cs(59) / +0.12320883239205686471647157725866D-30 /
  data ae11cs(60) / +0.22220174457553175317538581162666D-31 /

  data ae12cs( 1) / +0.63629589796747038767129887806803D+00 /
  data ae12cs( 2) / -0.13081168675067634385812671121135D+00 /
  data ae12cs( 3) / -0.84367410213053930014487662129752D-02 /
  data ae12cs( 4) / +0.26568491531006685413029428068906D-02 /
  data ae12cs( 5) / +0.32822721781658133778792170142517D-03 /
  data ae12cs( 6) / -0.23783447771430248269579807851050D-04 /
  data ae12cs( 7) / -0.11439804308100055514447076797047D-04 /
  data ae12cs( 8) / -0.14405943433238338455239717699323D-05 /
  data ae12cs( 9) / +0.52415956651148829963772818061664D-08 /
  data ae12cs(10) / +0.38407306407844323480979203059716D-07 /
  data ae12cs(11) / +0.85880244860267195879660515759344D-08 /
  data ae12cs(12) / +0.10219226625855003286339969553911D-08 /
  data ae12cs(13) / +0.21749132323289724542821339805992D-10 /
  data ae12cs(14) / -0.22090238142623144809523503811741D-10 /
  data ae12cs(15) / -0.63457533544928753294383622208801D-11 /
  data ae12cs(16) / -0.10837746566857661115340539732919D-11 /
  data ae12cs(17) / -0.11909822872222586730262200440277D-12 /
  data ae12cs(18) / -0.28438682389265590299508766008661D-14 /
  data ae12cs(19) / +0.25080327026686769668587195487546D-14 /
  data ae12cs(20) / +0.78729641528559842431597726421265D-15 /
  data ae12cs(21) / +0.15475066347785217148484334637329D-15 /
  data ae12cs(22) / +0.22575322831665075055272608197290D-16 /
  data ae12cs(23) / +0.22233352867266608760281380836693D-17 /
  data ae12cs(24) / +0.16967819563544153513464194662399D-19 /
  data ae12cs(25) / -0.57608316255947682105310087304533D-19 /
  data ae12cs(26) / -0.17591235774646878055625369408853D-19 /
  data ae12cs(27) / -0.36286056375103174394755328682666D-20 /
  data ae12cs(28) / -0.59235569797328991652558143488000D-21 /
  data ae12cs(29) / -0.76030380926310191114429136895999D-22 /
  data ae12cs(30) / -0.62547843521711763842641428479999D-23 /
  data ae12cs(31) / +0.25483360759307648606037606400000D-24 /
  data ae12cs(32) / +0.25598615731739857020168874666666D-24 /
  data ae12cs(33) / +0.71376239357899318800207052800000D-25 /
  data ae12cs(34) / +0.14703759939567568181578956800000D-25 /
  data ae12cs(35) / +0.25105524765386733555198634666666D-26 /
  data ae12cs(36) / +0.35886666387790890886583637333333D-27 /
  data ae12cs(37) / +0.39886035156771301763317759999999D-28 /
  data ae12cs(38) / +0.21763676947356220478805333333333D-29 /
  data ae12cs(39) / -0.46146998487618942367607466666666D-30 /
  data ae12cs(40) / -0.20713517877481987707153066666666D-30 /
  data ae12cs(41) / -0.51890378563534371596970666666666D-31 /

  data e11cs( 1) / -0.16113461655571494025720663927566180D+02 /
  data e11cs( 2) / +0.77940727787426802769272245891741497D+01 /
  data e11cs( 3) / -0.19554058188631419507127283812814491D+01 /
  data e11cs( 4) / +0.37337293866277945611517190865690209D+00 /
  data e11cs( 5) / -0.56925031910929019385263892220051166D-01 /
  data e11cs( 6) / +0.72110777696600918537847724812635813D-02 /
  data e11cs( 7) / -0.78104901449841593997715184089064148D-03 /
  data e11cs( 8) / +0.73880933562621681878974881366177858D-04 /
  data e11cs( 9) / -0.62028618758082045134358133607909712D-05 /
  data e11cs(10) / +0.46816002303176735524405823868362657D-06 /
  data e11cs(11) / -0.32092888533298649524072553027228719D-07 /
  data e11cs(12) / +0.20151997487404533394826262213019548D-08 /
  data e11cs(13) / -0.11673686816697793105356271695015419D-09 /
  data e11cs(14) / +0.62762706672039943397788748379615573D-11 /
  data e11cs(15) / -0.31481541672275441045246781802393600D-12 /
  data e11cs(16) / +0.14799041744493474210894472251733333D-13 /
  data e11cs(17) / -0.65457091583979673774263401588053333D-15 /
  data e11cs(18) / +0.27336872223137291142508012748799999D-16 /
  data e11cs(19) / -0.10813524349754406876721727624533333D-17 /
  data e11cs(20) / +0.40628328040434303295300348586666666D-19 /
  data e11cs(21) / -0.14535539358960455858914372266666666D-20 /
  data e11cs(22) / +0.49632746181648636830198442666666666D-22 /
  data e11cs(23) / -0.16208612696636044604866560000000000D-23 /
  data e11cs(24) / +0.50721448038607422226431999999999999D-25 /
  data e11cs(25) / -0.15235811133372207813973333333333333D-26 /
  data e11cs(26) / +0.44001511256103618696533333333333333D-28 /
  data e11cs(27) / -0.12236141945416231594666666666666666D-29 /
  data e11cs(28) / +0.32809216661066001066666666666666666D-31 /
  data e11cs(29) / -0.84933452268306432000000000000000000D-33 /

  data e12cs( 1) / -0.3739021479220279511668698204827D-01 /
  data e12cs( 2) / +0.4272398606220957726049179176528D-01 /
  data e12cs( 3) / -0.130318207984970054415392055219726D+00 /
  data e12cs( 4) / +0.144191240246988907341095893982137D-01 /
  data e12cs( 5) / -0.134617078051068022116121527983553D-02 /
  data e12cs( 6) / +0.107310292530637799976115850970073D-03 /
  data e12cs( 7) / -0.742999951611943649610283062223163D-05 /
  data e12cs( 8) / +0.453773256907537139386383211511827D-06 /
  data e12cs( 9) / -0.247641721139060131846547423802912D-07 /
  data e12cs(10) / +0.122076581374590953700228167846102D-08 /
  data e12cs(11) / -0.548514148064092393821357398028261D-10 /
  data e12cs(12) / +0.226362142130078799293688162377002D-11 /
  data e12cs(13) / -0.863589727169800979404172916282240D-13 /
  data e12cs(14) / +0.306291553669332997581032894881279D-14 /
  data e12cs(15) / -0.101485718855944147557128906734933D-15 /
  data e12cs(16) / +0.315482174034069877546855328426666D-17 /
  data e12cs(17) / -0.923604240769240954484015923200000D-19 /
  data e12cs(18) / +0.255504267970814002440435029333333D-20 /
  data e12cs(19) / -0.669912805684566847217882453333333D-22 /
  data e12cs(20) / +0.166925405435387319431987199999999D-23 /
  data e12cs(21) / -0.396254925184379641856000000000000D-25 /
  data e12cs(22) / +0.898135896598511332010666666666666D-27 /
  data e12cs(23) / -0.194763366993016433322666666666666D-28 /
  data e12cs(24) / +0.404836019024630033066666666666666D-30 /
  data e12cs(25) / -0.807981567699845120000000000000000D-32 /

  data ae13cs( 1) / -0.60577324664060345999319382737747D+00 /
  data ae13cs( 2) / -0.11253524348366090030649768852718D+00 /
  data ae13cs( 3) / +0.13432266247902779492487859329414D-01 /
  data ae13cs( 4) / -0.19268451873811457249246838991303D-02 /
  data ae13cs( 5) / +0.30911833772060318335586737475368D-03 /
  data ae13cs( 6) / -0.53564132129618418776393559795147D-04 /
  data ae13cs( 7) / +0.98278128802474923952491882717237D-05 /
  data ae13cs( 8) / -0.18853689849165182826902891938910D-05 /
  data ae13cs( 9) / +0.37494319356894735406964042190531D-06 /
  data ae13cs(10) / -0.76823455870552639273733465680556D-07 /
  data ae13cs(11) / +0.16143270567198777552956300060868D-07 /
  data ae13cs(12) / -0.34668022114907354566309060226027D-08 /
  data ae13cs(13) / +0.75875420919036277572889747054114D-09 /
  data ae13cs(14) / -0.16886433329881412573514526636703D-09 /
  data ae13cs(15) / +0.38145706749552265682804250927272D-10 /
  data ae13cs(16) / -0.87330266324446292706851718272334D-11 /
  data ae13cs(17) / +0.20236728645867960961794311064330D-11 /
  data ae13cs(18) / -0.47413283039555834655210340820160D-12 /
  data ae13cs(19) / +0.11221172048389864324731799928920D-12 /
  data ae13cs(20) / -0.26804225434840309912826809093395D-13 /
  data ae13cs(21) / +0.64578514417716530343580369067212D-14 /
  data ae13cs(22) / -0.15682760501666478830305702849194D-14 /
  data ae13cs(23) / +0.38367865399315404861821516441408D-15 /
  data ae13cs(24) / -0.94517173027579130478871048932556D-16 /
  data ae13cs(25) / +0.23434812288949573293896666439133D-16 /
  data ae13cs(26) / -0.58458661580214714576123194419882D-17 /
  data ae13cs(27) / +0.14666229867947778605873617419195D-17 /
  data ae13cs(28) / -0.36993923476444472706592538274474D-18 /
  data ae13cs(29) / +0.93790159936721242136014291817813D-19 /
  data ae13cs(30) / -0.23893673221937873136308224087381D-19 /
  data ae13cs(31) / +0.61150624629497608051934223837866D-20 /
  data ae13cs(32) / -0.15718585327554025507719853288106D-20 /
  data ae13cs(33) / +0.40572387285585397769519294491306D-21 /
  data ae13cs(34) / -0.10514026554738034990566367122773D-21 /
  data ae13cs(35) / +0.27349664930638667785806003131733D-22 /
  data ae13cs(36) / -0.71401604080205796099355574271999D-23 /
  data ae13cs(37) / +0.18705552432235079986756924211199D-23 /
  data ae13cs(38) / -0.49167468166870480520478020949333D-24 /
  data ae13cs(39) / +0.12964988119684031730916087125333D-24 /
  data ae13cs(40) / -0.34292515688362864461623940437333D-25 /
  data ae13cs(41) / +0.90972241643887034329104820906666D-26 /
  data ae13cs(42) / -0.24202112314316856489934847999999D-26 /
  data ae13cs(43) / +0.64563612934639510757670475093333D-27 /
  data ae13cs(44) / -0.17269132735340541122315987626666D-27 /
  data ae13cs(45) / +0.46308611659151500715194231466666D-28 /
  data ae13cs(46) / -0.12448703637214131241755170133333D-28 /
  data ae13cs(47) / +0.33544574090520678532907007999999D-29 /
  data ae13cs(48) / -0.90598868521070774437543935999999D-30 /
  data ae13cs(49) / +0.24524147051474238587273216000000D-30 /
  data ae13cs(50) / -0.66528178733552062817107967999999D-31 /

  data ae14cs( 1) / -0.1892918000753016825495679942820D+00 /
  data ae14cs( 2) / -0.8648117855259871489968817056824D-01 /
  data ae14cs( 3) / +0.7224101543746594747021514839184D-02 /
  data ae14cs( 4) / -0.8097559457557386197159655610181D-03 /
  data ae14cs( 5) / +0.1099913443266138867179251157002D-03 /
  data ae14cs( 6) / -0.1717332998937767371495358814487D-04 /
  data ae14cs( 7) / +0.2985627514479283322825342495003D-05 /
  data ae14cs( 8) / -0.5659649145771930056560167267155D-06 /
  data ae14cs( 9) / +0.1152680839714140019226583501663D-06 /
  data ae14cs(10) / -0.2495030440269338228842128765065D-07 /
  data ae14cs(11) / +0.5692324201833754367039370368140D-08 /
  data ae14cs(12) / -0.1359957664805600338490030939176D-08 /
  data ae14cs(13) / +0.3384662888760884590184512925859D-09 /
  data ae14cs(14) / -0.8737853904474681952350849316580D-10 /
  data ae14cs(15) / +0.2331588663222659718612613400470D-10 /
  data ae14cs(16) / -0.6411481049213785969753165196326D-11 /
  data ae14cs(17) / +0.1812246980204816433384359484682D-11 /
  data ae14cs(18) / -0.5253831761558460688819403840466D-12 /
  data ae14cs(19) / +0.1559218272591925698855028609825D-12 /
  data ae14cs(20) / -0.4729168297080398718476429369466D-13 /
  data ae14cs(21) / +0.1463761864393243502076199493808D-13 /
  data ae14cs(22) / -0.4617388988712924102232173623604D-14 /
  data ae14cs(23) / +0.1482710348289369323789239660371D-14 /
  data ae14cs(24) / -0.4841672496239229146973165734417D-15 /
  data ae14cs(25) / +0.1606215575700290408116571966188D-15 /
  data ae14cs(26) / -0.5408917538957170947895023784252D-16 /
  data ae14cs(27) / +0.1847470159346897881370231402310D-16 /
  data ae14cs(28) / -0.6395830792759094470500610425050D-17 /
  data ae14cs(29) / +0.2242780721699759457250233276170D-17 /
  data ae14cs(30) / -0.7961369173983947552744555308646D-18 /
  data ae14cs(31) / +0.2859308111540197459808619929272D-18 /
  data ae14cs(32) / -0.1038450244701137145900697137446D-18 /
  data ae14cs(33) / +0.3812040607097975780866841008319D-19 /
  data ae14cs(34) / -0.1413795417717200768717562723696D-19 /
  data ae14cs(35) / +0.5295367865182740958305442594815D-20 /
  data ae14cs(36) / -0.2002264245026825902137211131439D-20 /
  data ae14cs(37) / +0.7640262751275196014736848610918D-21 /
  data ae14cs(38) / -0.2941119006868787883311263523362D-21 /
  data ae14cs(39) / +0.1141823539078927193037691483586D-21 /
  data ae14cs(40) / -0.4469308475955298425247020718489D-22 /
  data ae14cs(41) / +0.1763262410571750770630491408520D-22 /
  data ae14cs(42) / -0.7009968187925902356351518262340D-23 /
  data ae14cs(43) / +0.2807573556558378922287757507515D-23 /
  data ae14cs(44) / -0.1132560944981086432141888891562D-23 /
  data ae14cs(45) / +0.4600574684375017946156764233727D-24 /
  data ae14cs(46) / -0.1881448598976133459864609148108D-24 /
  data ae14cs(47) / +0.7744916111507730845444328478037D-25 /
  data ae14cs(48) / -0.3208512760585368926702703826261D-25 /
  data ae14cs(49) / +0.1337445542910839760619930421384D-25 /
  data ae14cs(50) / -0.5608671881802217048894771735210D-26 /
  data ae14cs(51) / +0.2365839716528537483710069473279D-26 /
  data ae14cs(52) / -0.1003656195025305334065834526856D-26 /
  data ae14cs(53) / +0.4281490878094161131286642556927D-27 /
  data ae14cs(54) / -0.1836345261815318199691326958250D-27 /
  data ae14cs(55) / +0.7917798231349540000097468678144D-28 /
  data ae14cs(56) / -0.3431542358742220361025015775231D-28 /
  data ae14cs(57) / +0.1494705493897103237475066008917D-28 /
  data ae14cs(58) / -0.6542620279865705439739042420053D-29 /
  data ae14cs(59) / +0.2877581395199171114340487353685D-29 /
  data ae14cs(60) / -0.1271557211796024711027981200042D-29 /
  data ae14cs(61) / +0.5644615555648722522388044622506D-30 /
  data ae14cs(62) / -0.2516994994284095106080616830293D-30 /
  data ae14cs(63) / +0.1127259818927510206370368804181D-30 /
  data ae14cs(64) / -0.5069814875800460855562584719360D-31 /

  data ntae10 / 0 /
  data ntae11 / 0 /
  data ntae12 / 0 /
  data nte11 / 0 /
  data nte12 / 0 /
  data ntae13 / 0 /
  data ntae14 / 0 /
  data xmax / 0.0D+00 /

  if ( ntae10 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    ntae10 = r8_inits ( ae10cs, 50, eta )
    ntae11 = r8_inits ( ae11cs, 60, eta )
    ntae12 = r8_inits ( ae12cs, 41, eta )
    nte11 = r8_inits ( e11cs, 29, eta )
    nte12 = r8_inits ( e12cs, 25, eta )
    ntae13 = r8_inits ( ae13cs, 50, eta )
    ntae14 = r8_inits ( ae14cs, 64, eta )
    xmax = - log ( r8_mach ( 1 ) )
    xmax = xmax - log ( xmax )
  end if

  if ( x <= - 32.0D+00 ) then
    r8_e1 = exp ( - x ) / x * ( 1.0D+00 &
      + r8_csevl ( 64.0D+00 / x + 1.0D+00, ae10cs, ntae10 ) )
  else if ( x <= - 8.0D+00 ) then
    r8_e1 = exp ( - x ) / x * ( 1.0D+00 &
      + r8_csevl ( ( 64.0D+00 / x + 5.0D+00 ) / 3.0D+00, ae11cs, ntae11 ) )
  else if ( x <= - 4.0D+00 ) then
    r8_e1 = exp ( - x ) / x * (1.0D+00 &
      + r8_csevl ( 16.0D+00 / x + 3.0D+00, ae12cs, ntae12 ) )
  else if ( x <= - 1.0D+00 ) then
    r8_e1 = - log ( - x ) &
      + r8_csevl ( ( 2.0D+00 * x + 5.0D+00 ) / 3.0D+00, e11cs, nte11 )
  else if ( x == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_E1 - Fatal error!'
    write ( *, '(a)' ) '  X is zero.'
    stop
  else if ( x <= 1.0D+00 ) then
    r8_e1 = ( - log ( abs ( x ) ) - 0.6875D+00 + x ) &
      + r8_csevl ( x, e12cs, nte12 )
  else if ( x <= 4.0D+00 ) then
    r8_e1 = exp ( - x ) / x * ( 1.0D+00 &
      + r8_csevl ( ( 8.0D+00 / x - 5.0D+00 ) / 3.0D+00, ae13cs, ntae13 ) )
  else if ( x <= xmax ) then
    r8_e1 = exp ( - x ) / x * ( 1.0D+00 &
      + r8_csevl ( 8.0D+00 / x - 1.0D+00, ae14cs, ntae14 ) )
  else
    r8_e1 = 0.0D+00
  end if

  return
end
function r8_ei ( x )

!*****************************************************************************80
!
!! R8_EI evaluates the exponential integral Ei for an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_EI, the exponential integral Ei.
!
  implicit none

  real ( kind = 8 ) r8_e1
  real ( kind = 8 ) r8_ei
  real ( kind = 8 ) x

  r8_ei = - r8_e1 ( - x )

  return
end
function r8_erf ( x )

!*****************************************************************************80
!
!! R8_ERF evaluates the error function of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_ERF, the error function of X.
!
  implicit none

  real ( kind = 8 ) erfcs(21)
  integer ( kind = 4 ) nterf
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_erf
  real ( kind = 8 ) r8_erfc
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) sqrtpi
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) y

  save erfcs
  save nterf
  save sqrtpi
  save xbig

  data erfcs(  1) / -0.49046121234691808039984544033376D-01 /
  data erfcs(  2) / -0.14226120510371364237824741899631D+00 /
  data erfcs(  3) / +0.10035582187599795575754676712933D-01 /
  data erfcs(  4) / -0.57687646997674847650827025509167D-03 /
  data erfcs(  5) / +0.27419931252196061034422160791471D-04 /
  data erfcs(  6) / -0.11043175507344507604135381295905D-05 /
  data erfcs(  7) / +0.38488755420345036949961311498174D-07 /
  data erfcs(  8) / -0.11808582533875466969631751801581D-08 /
  data erfcs(  9) / +0.32334215826050909646402930953354D-10 /
  data erfcs( 10) / -0.79910159470045487581607374708595D-12 /
  data erfcs( 11) / +0.17990725113961455611967245486634D-13 /
  data erfcs( 12) / -0.37186354878186926382316828209493D-15 /
  data erfcs( 13) / +0.71035990037142529711689908394666D-17 /
  data erfcs( 14) / -0.12612455119155225832495424853333D-18 /
  data erfcs( 15) / +0.20916406941769294369170500266666D-20 /
  data erfcs( 16) / -0.32539731029314072982364160000000D-22 /
  data erfcs( 17) / +0.47668672097976748332373333333333D-24 /
  data erfcs( 18) / -0.65980120782851343155199999999999D-26 /
  data erfcs( 19) / +0.86550114699637626197333333333333D-28 /
  data erfcs( 20) / -0.10788925177498064213333333333333D-29 /
  data erfcs( 21) / +0.12811883993017002666666666666666D-31 /

  data nterf / 0 /
  data sqrtpi / 1.77245385090551602729816748334115D+00 /
  data xbig / 0.0D+00 /

  if ( nterf == 0 ) then
    nterf = r8_inits ( erfcs, 21, 0.1D+00 * r8_mach ( 3 ) )
    xbig = sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) )
    sqeps = sqrt ( 2.0D+00 * r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    value = 2.0D+00 * x / sqrtpi
  else if ( y <= 1.0D+00 ) then
    value = x * ( 1.0D+00 &
      + r8_csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf ) )
  else if ( y <= xbig ) then
    value = 1.0D+00 - r8_erfc ( y )
    if ( x < 0.0D+00 ) then
      value = - value
    end if
  else
    value = 1.0D+00
    if ( x < 0.0D+00 ) then
      value = - value
    end if
  end if

  r8_erf = value

  return
end
function r8_erfc ( x )

!*****************************************************************************80
!
!! R8_ERFC evaluates the co-error function of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_ERFC, the co-error function of X.
!
  implicit none

  real ( kind = 8 ) erc2cs(49)
  real ( kind = 8 ) erfccs(59)
  real ( kind = 8 ) erfcs(21)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nterc2
  integer ( kind = 4 ) nterf
  integer ( kind = 4 ) nterfc
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_erfc
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) sqrtpi
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save erfccs
  save erfcs
  save erc2cs
  save nterc2
  save nterf
  save nterfc
  save sqrtpi
  save xmax
  save xsml

  data erfcs(  1) / -0.49046121234691808039984544033376D-01 /
  data erfcs(  2) / -0.14226120510371364237824741899631D+00 /
  data erfcs(  3) / +0.10035582187599795575754676712933D-01 /
  data erfcs(  4) / -0.57687646997674847650827025509167D-03 /
  data erfcs(  5) / +0.27419931252196061034422160791471D-04 /
  data erfcs(  6) / -0.11043175507344507604135381295905D-05 /
  data erfcs(  7) / +0.38488755420345036949961311498174D-07 /
  data erfcs(  8) / -0.11808582533875466969631751801581D-08 /
  data erfcs(  9) / +0.32334215826050909646402930953354D-10 /
  data erfcs( 10) / -0.79910159470045487581607374708595D-12 /
  data erfcs( 11) / +0.17990725113961455611967245486634D-13 /
  data erfcs( 12) / -0.37186354878186926382316828209493D-15 /
  data erfcs( 13) / +0.71035990037142529711689908394666D-17 /
  data erfcs( 14) / -0.12612455119155225832495424853333D-18 /
  data erfcs( 15) / +0.20916406941769294369170500266666D-20 /
  data erfcs( 16) / -0.32539731029314072982364160000000D-22 /
  data erfcs( 17) / +0.47668672097976748332373333333333D-24 /
  data erfcs( 18) / -0.65980120782851343155199999999999D-26 /
  data erfcs( 19) / +0.86550114699637626197333333333333D-28 /
  data erfcs( 20) / -0.10788925177498064213333333333333D-29 /
  data erfcs( 21) / +0.12811883993017002666666666666666D-31 /

  data erc2cs(  1) / -0.6960134660230950112739150826197D-01 /
  data erc2cs(  2) / -0.4110133936262089348982212084666D-01 /
  data erc2cs(  3) / +0.3914495866689626881561143705244D-02 /
  data erc2cs(  4) / -0.4906395650548979161280935450774D-03 /
  data erc2cs(  5) / +0.7157479001377036380760894141825D-04 /
  data erc2cs(  6) / -0.1153071634131232833808232847912D-04 /
  data erc2cs(  7) / +0.1994670590201997635052314867709D-05 /
  data erc2cs(  8) / -0.3642666471599222873936118430711D-06 /
  data erc2cs(  9) / +0.6944372610005012589931277214633D-07 /
  data erc2cs( 10) / -0.1371220902104366019534605141210D-07 /
  data erc2cs( 11) / +0.2788389661007137131963860348087D-08 /
  data erc2cs( 12) / -0.5814164724331161551864791050316D-09 /
  data erc2cs( 13) / +0.1238920491752753181180168817950D-09 /
  data erc2cs( 14) / -0.2690639145306743432390424937889D-10 /
  data erc2cs( 15) / +0.5942614350847910982444709683840D-11 /
  data erc2cs( 16) / -0.1332386735758119579287754420570D-11 /
  data erc2cs( 17) / +0.3028046806177132017173697243304D-12 /
  data erc2cs( 18) / -0.6966648814941032588795867588954D-13 /
  data erc2cs( 19) / +0.1620854541053922969812893227628D-13 /
  data erc2cs( 20) / -0.3809934465250491999876913057729D-14 /
  data erc2cs( 21) / +0.9040487815978831149368971012975D-15 /
  data erc2cs( 22) / -0.2164006195089607347809812047003D-15 /
  data erc2cs( 23) / +0.5222102233995854984607980244172D-16 /
  data erc2cs( 24) / -0.1269729602364555336372415527780D-16 /
  data erc2cs( 25) / +0.3109145504276197583836227412951D-17 /
  data erc2cs( 26) / -0.7663762920320385524009566714811D-18 /
  data erc2cs( 27) / +0.1900819251362745202536929733290D-18 /
  data erc2cs( 28) / -0.4742207279069039545225655999965D-19 /
  data erc2cs( 29) / +0.1189649200076528382880683078451D-19 /
  data erc2cs( 30) / -0.3000035590325780256845271313066D-20 /
  data erc2cs( 31) / +0.7602993453043246173019385277098D-21 /
  data erc2cs( 32) / -0.1935909447606872881569811049130D-21 /
  data erc2cs( 33) / +0.4951399124773337881000042386773D-22 /
  data erc2cs( 34) / -0.1271807481336371879608621989888D-22 /
  data erc2cs( 35) / +0.3280049600469513043315841652053D-23 /
  data erc2cs( 36) / -0.8492320176822896568924792422399D-24 /
  data erc2cs( 37) / +0.2206917892807560223519879987199D-24 /
  data erc2cs( 38) / -0.5755617245696528498312819507199D-25 /
  data erc2cs( 39) / +0.1506191533639234250354144051199D-25 /
  data erc2cs( 40) / -0.3954502959018796953104285695999D-26 /
  data erc2cs( 41) / +0.1041529704151500979984645051733D-26 /
  data erc2cs( 42) / -0.2751487795278765079450178901333D-27 /
  data erc2cs( 43) / +0.7290058205497557408997703680000D-28 /
  data erc2cs( 44) / -0.1936939645915947804077501098666D-28 /
  data erc2cs( 45) / +0.5160357112051487298370054826666D-29 /
  data erc2cs( 46) / -0.1378419322193094099389644800000D-29 /
  data erc2cs( 47) / +0.3691326793107069042251093333333D-30 /
  data erc2cs( 48) / -0.9909389590624365420653226666666D-31 /
  data erc2cs( 49) / +0.2666491705195388413323946666666D-31 /

  data erfccs(  1) / +0.715179310202924774503697709496D-01 /
  data erfccs(  2) / -0.265324343376067157558893386681D-01 /
  data erfccs(  3) / +0.171115397792085588332699194606D-02 /
  data erfccs(  4) / -0.163751663458517884163746404749D-03 /
  data erfccs(  5) / +0.198712935005520364995974806758D-04 /
  data erfccs(  6) / -0.284371241276655508750175183152D-05 /
  data erfccs(  7) / +0.460616130896313036969379968464D-06 /
  data erfccs(  8) / -0.822775302587920842057766536366D-07 /
  data erfccs(  9) / +0.159214187277090112989358340826D-07 /
  data erfccs( 10) / -0.329507136225284321486631665072D-08 /
  data erfccs( 11) / +0.722343976040055546581261153890D-09 /
  data erfccs( 12) / -0.166485581339872959344695966886D-09 /
  data erfccs( 13) / +0.401039258823766482077671768814D-10 /
  data erfccs( 14) / -0.100481621442573113272170176283D-10 /
  data erfccs( 15) / +0.260827591330033380859341009439D-11 /
  data erfccs( 16) / -0.699111056040402486557697812476D-12 /
  data erfccs( 17) / +0.192949233326170708624205749803D-12 /
  data erfccs( 18) / -0.547013118875433106490125085271D-13 /
  data erfccs( 19) / +0.158966330976269744839084032762D-13 /
  data erfccs( 20) / -0.472689398019755483920369584290D-14 /
  data erfccs( 21) / +0.143587337678498478672873997840D-14 /
  data erfccs( 22) / -0.444951056181735839417250062829D-15 /
  data erfccs( 23) / +0.140481088476823343737305537466D-15 /
  data erfccs( 24) / -0.451381838776421089625963281623D-16 /
  data erfccs( 25) / +0.147452154104513307787018713262D-16 /
  data erfccs( 26) / -0.489262140694577615436841552532D-17 /
  data erfccs( 27) / +0.164761214141064673895301522827D-17 /
  data erfccs( 28) / -0.562681717632940809299928521323D-18 /
  data erfccs( 29) / +0.194744338223207851429197867821D-18 /
  data erfccs( 30) / -0.682630564294842072956664144723D-19 /
  data erfccs( 31) / +0.242198888729864924018301125438D-19 /
  data erfccs( 32) / -0.869341413350307042563800861857D-20 /
  data erfccs( 33) / +0.315518034622808557122363401262D-20 /
  data erfccs( 34) / -0.115737232404960874261239486742D-20 /
  data erfccs( 35) / +0.428894716160565394623737097442D-21 /
  data erfccs( 36) / -0.160503074205761685005737770964D-21 /
  data erfccs( 37) / +0.606329875745380264495069923027D-22 /
  data erfccs( 38) / -0.231140425169795849098840801367D-22 /
  data erfccs( 39) / +0.888877854066188552554702955697D-23 /
  data erfccs( 40) / -0.344726057665137652230718495566D-23 /
  data erfccs( 41) / +0.134786546020696506827582774181D-23 /
  data erfccs( 42) / -0.531179407112502173645873201807D-24 /
  data erfccs( 43) / +0.210934105861978316828954734537D-24 /
  data erfccs( 44) / -0.843836558792378911598133256738D-25 /
  data erfccs( 45) / +0.339998252494520890627359576337D-25 /
  data erfccs( 46) / -0.137945238807324209002238377110D-25 /
  data erfccs( 47) / +0.563449031183325261513392634811D-26 /
  data erfccs( 48) / -0.231649043447706544823427752700D-26 /
  data erfccs( 49) / +0.958446284460181015263158381226D-27 /
  data erfccs( 50) / -0.399072288033010972624224850193D-27 /
  data erfccs( 51) / +0.167212922594447736017228709669D-27 /
  data erfccs( 52) / -0.704599152276601385638803782587D-28 /
  data erfccs( 53) / +0.297976840286420635412357989444D-28 /
  data erfccs( 54) / -0.126252246646061929722422632994D-28 /
  data erfccs( 55) / +0.539543870454248793985299653154D-29 /
  data erfccs( 56) / -0.238099288253145918675346190062D-29 /
  data erfccs( 57) / +0.109905283010276157359726683750D-29 /
  data erfccs( 58) / -0.486771374164496572732518677435D-30 /
  data erfccs( 59) / +0.152587726411035756763200828211D-30 /

  data nterc2 / 0 /
  data nterf / 0 /
  data nterfc / 0 /
  data sqrtpi / 1.77245385090551602729816748334115D+00 /
  data xmax / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( nterf == 0 ) then

    eta = 0.1D+00 * r8_mach ( 3 )
    nterf = r8_inits ( erfcs, 21, eta )
    nterfc = r8_inits ( erfccs, 59, eta )
    nterc2 = r8_inits ( erc2cs, 49, eta )

    xsml = - sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) )
    xmax = sqrt (- log ( sqrtpi * r8_mach ( 1 ) ) )
    xmax = xmax - 0.5D+00 * log ( xmax ) / xmax - 0.01D+00
    sqeps = sqrt ( 2.0D+00 * r8_mach ( 3 ) )

  end if

  if ( x <= xsml ) then

    r8_erfc = 2.0D+00
    return

  end if

  if ( xmax < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ERFC - Warning!'
    write ( *, '(a)' ) '  X so big that ERFC underflows.'
    r8_erfc = 0.0D+00
    return
  end if

  y = abs ( x )

  if ( y < sqeps ) then
    r8_erfc = 1.0D+00 - 2.0D+00 * x / sqrtpi
    return
  else if ( y <= 1.0D+00 ) then
    r8_erfc = 1.0D+00 - x * ( 1.0D+00 &
      + r8_csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf ) )
    return
  end if

  y = y * y

  if ( y <= 4.0D+00 ) then
    r8_erfc = exp ( - y ) / abs ( x ) * ( 0.5D+00 &
      + r8_csevl ( ( 8.0D+00 / y - 5.0D+00 ) / 3.0D+00, erc2cs, nterc2 ) )
  else 
    r8_erfc = exp ( - y ) / abs ( x ) * ( 0.5D+00 &
      + r8_csevl ( 8.0D+00 / y - 1.0D+00, erfccs, nterfc ) )
  end if

  if ( x < 0.0D+00 ) then
    r8_erfc = 2.0D+00 - r8_erfc
  end if

  return
end
function r8_exp ( x )

!*****************************************************************************80
!
!! R8_EXP evaluates the exponential of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_EXP, the exponential of X.
!
  implicit none

  real ( kind = 8 ) aln216
  real ( kind = 8 ) expcs(14)
  real ( kind = 8 ) f
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n16
  integer ( kind = 4 ) ndx
  integer ( kind = 4 ) nterms
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_exp
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_pak
  real ( kind = 8 ) twon16(17)
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) xint
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y

  save aln216
  save expcs
  save nterms
  save twon16
  save xmax
  save xmin

  data expcs(  1) / +0.866569493314985712733404647266231D-01 /
  data expcs(  2) / +0.938494869299839561896336579701203D-03 /
  data expcs(  3) / +0.677603970998168264074353014653601D-05 /
  data expcs(  4) / +0.366931200393805927801891250687610D-07 /
  data expcs(  5) / +0.158959053617461844641928517821508D-09 /
  data expcs(  6) / +0.573859878630206601252990815262106D-12 /
  data expcs(  7) / +0.177574448591421511802306980226000D-14 /
  data expcs(  8) / +0.480799166842372422675950244533333D-17 /
  data expcs(  9) / +0.115716376881828572809260000000000D-19 /
  data expcs( 10) / +0.250650610255497719932458666666666D-22 /
  data expcs( 11) / +0.493571708140495828480000000000000D-25 /
  data expcs( 12) / +0.890929572740634240000000000000000D-28 /
  data expcs( 13) / +0.148448062907997866666666666666666D-30 /
  data expcs( 14) / +0.229678916630186666666666666666666D-33 /

  data twon16(  1) / +0.0D+00 /
  data twon16(  2) / +0.44273782427413840321966478739929D-01 /
  data twon16(  3) / +0.90507732665257659207010655760707D-01 /
  data twon16(  4) / +0.13878863475669165370383028384151D+00 /
  data twon16(  5) / +0.18920711500272106671749997056047D+00 /
  data twon16(  6) / +0.24185781207348404859367746872659D+00 /
  data twon16(  7) / +0.29683955465100966593375411779245D+00 /
  data twon16(  8) / +0.35425554693689272829801474014070D+00 /
  data twon16(  9) / +0.41421356237309504880168872420969D+00 /
  data twon16( 10) / +0.47682614593949931138690748037404D+00 /
  data twon16( 11) / +0.54221082540794082361229186209073D+00 /
  data twon16( 12) / +0.61049033194925430817952066735740D+00 /
  data twon16( 13) / +0.68179283050742908606225095246642D+00 /
  data twon16( 14) / +0.75625216037329948311216061937531D+00 /
  data twon16( 15) / +0.83400808640934246348708318958828D+00 /
  data twon16( 16) / +0.91520656139714729387261127029583D+00 /
  data twon16( 17) / +1.0D+00 /

  data aln216 / +0.83120654223414517758794896030274D-01 /
  data nterms / 0 /
  data xmax / 0.0D+00 /
  data xmin / 0.0D+00 /

  if ( nterms == 0 ) then
    nterms = r8_inits ( expcs, 14, 0.1D+00 * r8_mach ( 3 ) )
    xmin = log ( r8_mach ( 1 ) ) + 0.01D+00
    xmax = log ( r8_mach ( 2 ) ) - 0.001D+00
  end if

  if ( x < xmin ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_EXP - Warning!'
    write ( *, '(a)' ) '  X so small that exp(X) underflows.'
    value = 0.0D+00

  else if ( x <= xmax ) then

    xint = aint ( x )
    y = x - xint

    y = 23.0D+00 * y + x * aln216
    n = int ( y )
    f = y - real ( n, kind = 8 )
    n = 23.0D+00 * xint + real ( n, kind = 8 )
    n16 = n / 16
    if ( n < 0 ) then
      n16 = n16 - 1
    end if
    ndx = n - 16 * n16 + 1

    value = 1.0D+00 + ( twon16(ndx) &
      + f * ( 1.0D+00 + twon16(ndx) ) &
      * r8_csevl ( f, expcs, nterms ) )

    value = r8_pak ( value, n16 )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_EXP - Fatal error!'
    write ( *, '(a)' ) '  X so large that exp(X) overflows.'
    stop

  end if

  r8_exp = value

  return
end
function r8_exprel ( x )

!*****************************************************************************80
!
!! R8_EXPREL evaluates the exponential relative error term of an R8 argument.
!
!  Discussion:
!
!    The relative error term is ( exp ( x ) - 1 ) / x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_EXPREL, the exponential relative error term
!    at X.
!
  implicit none

  real ( kind = 8 ) absx
  real ( kind = 8 ) alneps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nterms
  real ( kind = 8 ) r8_exprel
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xbnd
  real ( kind = 8 ) xln
  real ( kind = 8 ) xn

  save nterms
  save xbnd

  data nterms / 0 /
  data xbnd / 0.0D+00 /

  if ( nterms == 0 ) then
    alneps = log ( r8_mach ( 3 ) )
    xn = 3.72D+00 - 0.3D+00 * alneps
    xln = log ( ( xn + 1.0D+00 ) / 1.36D+00 )
    nterms = int ( xn - ( xn * xln + alneps ) / ( xln + 1.36D+00 ) + 1.5D+00 )
    xbnd = r8_mach ( 3 )
  end if

  absx = abs ( x )

  if ( absx < xbnd ) then
    r8_exprel = 1.0D+00
  else if ( absx <= 0.5D+00 ) then
    r8_exprel = 0.0D+00
    do i = 1, nterms
      r8_exprel = 1.0D+00 + r8_exprel * x / real ( nterms + 2 - i, kind = 8 )
    end do
  else
    r8_exprel = ( exp ( x ) - 1.0D+00 ) / x
  end if

  return
end
function r8_fac ( n )

!*****************************************************************************80
!
!! R8_FAC evaluates the factorial of an I4 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument.
!
!    Output, real ( kind = 8 ) R8_FAC, the factorial of N.
!
  implicit none

  real ( kind = 8 ) facn(31)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nmax
  real ( kind = 8 ) r8_fac
  real ( kind = 8 ) r8_lgmc
  real ( kind = 8 ) sq2pil
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  save facn
  save nmax
  save sq2pil

  data facn(  1) / +0.100000000000000000000000000000000D+01 /
  data facn(  2) / +0.100000000000000000000000000000000D+01 /
  data facn(  3) / +0.200000000000000000000000000000000D+01 /
  data facn(  4) / +0.600000000000000000000000000000000D+01 /
  data facn(  5) / +0.240000000000000000000000000000000D+02 /
  data facn(  6) / +0.120000000000000000000000000000000D+03 /
  data facn(  7) / +0.720000000000000000000000000000000D+03 /
  data facn(  8) / +0.504000000000000000000000000000000D+04 /
  data facn(  9) / +0.403200000000000000000000000000000D+05 /
  data facn( 10) / +0.362880000000000000000000000000000D+06 /
  data facn( 11) / +0.362880000000000000000000000000000D+07 /
  data facn( 12) / +0.399168000000000000000000000000000D+08 /
  data facn( 13) / +0.479001600000000000000000000000000D+09 /
  data facn( 14) / +0.622702080000000000000000000000000D+10 /
  data facn( 15) / +0.871782912000000000000000000000000D+11 /
  data facn( 16) / +0.130767436800000000000000000000000D+13 /
  data facn( 17) / +0.209227898880000000000000000000000D+14 /
  data facn( 18) / +0.355687428096000000000000000000000D+15 /
  data facn( 19) / +0.640237370572800000000000000000000D+16 /
  data facn( 20) / +0.121645100408832000000000000000000D+18 /
  data facn( 21) / +0.243290200817664000000000000000000D+19 /
  data facn( 22) / +0.510909421717094400000000000000000D+20 /
  data facn( 23) / +0.112400072777760768000000000000000D+22 /
  data facn( 24) / +0.258520167388849766400000000000000D+23 /
  data facn( 25) / +0.620448401733239439360000000000000D+24 /
  data facn( 26) / +0.155112100433309859840000000000000D+26 /
  data facn( 27) / +0.403291461126605635584000000000000D+27 /
  data facn( 28) / +0.108888694504183521607680000000000D+29 /
  data facn( 29) / +0.304888344611713860501504000000000D+30 /
  data facn( 30) / +0.884176199373970195454361600000000D+31 /
  data facn( 31) / +0.265252859812191058636308480000000D+33 /

  data nmax / 0 /
  data sq2pil / 0.91893853320467274178032973640562D+00 /

  if ( nmax == 0 ) then
    call r8_gaml ( xmin, xmax )
    nmax = int ( xmax - 1.0D+00 )
  end if

  if ( n < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_FAC - Fatal error!'
    write ( *, '(a)' ) '  Input argument is negative.'
    stop

  else if ( n <= 30 ) then

    r8_fac = facn(n+1)
 
  else if ( n <= nmax ) then

    x = real ( n + 1, kind = 8 )
    r8_fac = exp ( ( x - 0.5D+00 ) * log ( x ) - x + sq2pil + r8_lgmc ( x ) )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_FAC - Fatal error!'
    write ( *, '(a)' ) '  Factorial overflows.'
    stop

  end if

  return
end
function r8_gami ( a, x )

!*****************************************************************************80
!
!! R8_GAMI evaluates the incomplete gamma function for an R8 argument.
!
!  Discussion:
!
!    GAMI = Integral ( 0 <= T <= X ) exp ( - t ) * t^( a - 1 )  dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_GAMI, the value of the incomplete 
!    gamma function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) factor
  real ( kind = 8 ) r8_gami
  real ( kind = 8 ) r8_gamit
  real ( kind = 8 ) r8_lngam
  real ( kind = 8 ) x

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_GAMI - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    stop
  end if

  if ( x < 0.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_GAMI - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop

  else if ( x == 0.0D+00 ) then

    r8_gami = 0.0D+00

  else

    factor = exp ( r8_lngam ( a ) + a * log ( x ) )

    r8_gami = factor * r8_gamit ( a, x )

  end if

  return
end
function r8_gamic ( a, x )

!*****************************************************************************80
!
!! R8_GAMIC evaluates the complementary incomplete gamma function.
!
!  Discussion:
!
!    GAMIC = integral ( x <= t < oo ) exp(-t) * t^(a-1) dt
!
!    GAMIC is evaluated for arbitrary real values of A and non-negative
!    values X (even though GAMIC is defined for X < 0.0), except that
!    for X = 0 and A <= 0.0, GAMIC is undefined.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!    Walter Gautschi,
!    A Computational Procedure for Incomplete Gamma Functions,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 4, December 1979, pages 466-481.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) R8_GAMIC, the value of the incomplete 
!    gamma function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) aeps
  real ( kind = 8 ) ainta
  real ( kind = 8 ) algap1
  real ( kind = 8 ) alneps
  real ( kind = 8 ) alngs
  real ( kind = 8 ) alx
  real ( kind = 8 ) bot
  real ( kind = 8 ) e
  real ( kind = 8 ) eps
  real ( kind = 8 ) gstar
  real ( kind = 8 ) h
  integer ( kind = 4 ) izero
  integer ( kind = 4 ) ma
  real ( kind = 8 ) r8_gamic
  real ( kind = 8 ) r8_gmic
  real ( kind = 8 ) r8_gmit
  real ( kind = 8 ) r8_lgic
  real ( kind = 8 ) r8_lgit
  real ( kind = 8 ) r8_lngam
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) sga
  real ( kind = 8 ) sgng
  real ( kind = 8 ) sgngam
  real ( kind = 8 ) sgngs
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  save alneps
  save bot
  save eps

  data alneps / 0.0D+00 /
  data bot / 0.0D+00 /
  data eps / 0.0D+00 /

  if ( eps == 0.0D+00 ) then
    eps = 0.5D+00 * r8_mach ( 3 )
    sqeps = sqrt ( r8_mach ( 4 ) )
    alneps = - log ( r8_mach ( 3 ) )
    bot = log ( r8_mach ( 1 ) )
  end if

  if ( x < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_GAMIC - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop
  end if

  if ( x == 0.0D+00 ) then

    if ( a <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAMIC - Fatal error!'
      write ( *, '(a)' ) '  X = 0 and A <= 0.'
      stop
    end if
  
    r8_gamic = exp ( r8_lngam ( a + 1.0D+00 ) - log ( a ) )

    return
  end if

  alx = log ( x )
  if ( a < 0.0D+00 ) then
    sga = - 1.0D+00
  else
    sga = + 1.0D+00
  end if

  ainta = aint ( a + 0.5D+00 * sga )
  aeps = a - ainta

  izero = 0

  if ( x < 1.0D+00 ) then

    if ( a <= 0.5D+00 .and. abs ( aeps ) <= 0.001D+00 ) then

      if ( - ainta <= 1.0D+00 ) then
        e = 2.0D+00
      else
        e = 2.0D+00 * ( - ainta + 2.0D+00 ) / ( ainta * ainta - 1.0D+00 )
      end if

      e = e - alx * x ** ( - 0.001D+00 )

      if ( e * abs ( aeps ) <= eps ) then
        r8_gamic = r8_gmic ( a, x, alx )
        return
      end if

    end if

    call r8_lgams ( a + 1.0D+00, algap1, sgngam )
    gstar = r8_gmit ( a, x, algap1, sgngam, alx )

    if ( gstar == 0.0D+00 ) then
      izero = 1
    else
      alngs = log ( abs ( gstar ) )
      sgngs = r8_sign ( gstar )
    end if

  else

    if ( a < x ) then
      r8_gamic = exp ( r8_lgic ( a, x, alx ) )
      return
    end if

    sgngam = 1.0D+00
    algap1 = r8_lngam ( a + 1.0D+00 )
    sgngs = 1.0D+00
    alngs = r8_lgit ( a, x, algap1 )

  end if

  h = 1.0D+00

  if ( izero /= 1 ) then

    t = a * alx + alngs

    if ( alneps < t ) then
      sgng = - sgngs * sga * sgngam
      t = t + algap1 - log ( abs ( a ) )
      r8_gamic = sgng * exp ( t )
      return
    end if

    if ( - alneps < t ) then
      h = 1.0D+00 - sgngs * exp ( t )
    end if

  end if

  sgng = r8_sign ( h ) * sga * sgngam
  t = log ( abs ( h ) ) + algap1 - log ( abs ( a ) )
  r8_gamic = sgng * exp ( t )

  return
end
function r8_gamit ( a, x )

!*****************************************************************************80
!
!! R8_GAMIT evaluates Tricomi's incomplete gamma function for an R8 argument.
!
!  Discussion:
!
!      GAMIT = x^(-a) / gamma(a) 
!        * Integral ( 0 <= t <= x ) exp(-t) * t^(a-1) dt
!
!    with analytic continuation for a <= 0.0.  Gamma(x) is the complete
!    gamma function of X.  GAMIT is evaluated for arbitrary real values of
!    A and for non-negative values of X (even though GAMIT is defined for
!    X < 0.0).
!
!    A slight deterioration of 2 or 3 digits accuracy will occur when
!    gamit is very large or very small in absolute value, because log-
!    arithmic variables are used.  Also, if the parameter A is very close
!    to a negative integer (but not a negative integer), there is a loss
!    of accuracy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!    Walter Gautschi,
!    A Computational Procedure for Incomplete Gamma Functions,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 4, December 1979, pages 466-481.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_GAMIT, the function value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) aeps
  real ( kind = 8 ) ainta
  real ( kind = 8 ) algap1
  real ( kind = 8 ) alneps
  real ( kind = 8 ) alng
  real ( kind = 8 ) alx
  real ( kind = 8 ) bot
  real ( kind = 8 ) h
  real ( kind = 8 ) r8_gamit
  real ( kind = 8 ) r8_gamr
  real ( kind = 8 ) r8_gmit
  real ( kind = 8 ) r8_lgic
  real ( kind = 8 ) r8_lgit
  real ( kind = 8 ) r8_lngam
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sga
  real ( kind = 8 ) sgngam
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) x
 
  save alneps
  save bot

  data alneps / 0.0D+00 /
  data bot / 0.0D+00 /

  if ( alneps == 0.0D+00 ) then
    alneps = - log ( r8_mach ( 3 ) )
    sqeps = sqrt ( r8_mach ( 4 ) )
    bot = log ( r8_mach ( 1 ) )
  end if

  if ( x < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_GAMIT - Fatal error!'
    write ( *, '(a)' ) '  X is negative.'
    stop
  else if ( x == 0.0D+00 ) then
    alx = 0.0D+00
  else
    alx = log ( x )
  end if

  if ( a < 0.0D+00 ) then
    sga = - 1.0D+00
  else
    sga = + 1.0D+00
  end if

  ainta = aint ( a + 0.5D+00 * sga )
  aeps = a - ainta

  if ( x == 0.0D+00 ) then
    if ( 0.0D+00 < ainta .or. aeps /= 0.0D+00 ) then
      r8_gamit = r8_gamr ( a + 1.0D+00 )
    else
      r8_gamit = 0.0D+00
    end if
    return
  end if

  if ( x <= 1.0D+00 ) then
    if ( - 0.5D+00 <= a .or. aeps /= 0.0D+00 ) then
      call r8_lgams ( a + 1.0D+00, algap1, sgngam )
    end if
    r8_gamit = r8_gmit ( a, x, algap1, sgngam, alx )
    return
  end if

  if ( x <= a ) then
    t = r8_lgit (a, x, r8_lngam ( a + 1.0D+00 ) )
    r8_gamit = exp ( t )
    return
  end if

  alng = r8_lgic ( a, x, alx )
!
!  Evaluate in terms of log (r8_gamic (a, x))
!
  h = 1.0D+00

  if ( aeps /= 0.0D+00 .or. 0.0D+00 < ainta ) then

    call r8_lgams ( a + 1.0D+00, algap1, sgngam )
    t = log ( abs ( a ) ) + alng - algap1

    if ( alneps < t ) then
      t = t - a * alx
      r8_gamit = - sga * sgngam * exp ( t )
      return
    end if

    if ( - alneps < t ) then
      h = 1.0D+00 - sga * sgngam * exp ( t )
    end if

  end if

  t = - a * alx + log ( abs ( h ) )

  if ( h < 0.0D+00 ) then
    r8_gamit = - exp ( t )
  else
    r8_gamit = + exp ( t )
  end if

  return
end
subroutine r8_gaml ( xmin, xmax )

!*****************************************************************************80
!
!! R8_GAML evaluates bounds for an R8 argument of the gamma function.
!
!  Discussion:
!
!    This function calculates the minimum and maximum legal bounds 
!    for X in the evaluation of GAMMA ( X ).
!
!    XMIN and XMAX are not the only bounds, but they are the only 
!    non-trivial ones to calculate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) XMIN, XMAX, the bounds.
!
  implicit none

  real ( kind = 8 ) alnbig
  real ( kind = 8 ) alnsml
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) xln
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xold

  alnsml = log ( r8_mach ( 1 ) )
  xmin = - alnsml

  do i = 1, 10

    xold = xmin
    xln = log ( xmin )
    xmin = xmin - xmin * ( ( xmin + 0.5D+00 ) * xln - xmin &
      - 0.2258D+00 + alnsml ) / ( xmin * xln + 0.5D+00 )

    if ( abs ( xmin - xold ) < 0.005D+00 ) then

      xmin = - xmin + 0.01D+00

      alnbig = log ( r8_mach ( 2 ) )
      xmax = alnbig

      do j = 1, 10

        xold = xmax
        xln = log ( xmax )
        xmax = xmax - xmax * ( ( xmax - 0.5D+00 ) * xln - xmax &
          + 0.9189D+00 - alnbig ) / ( xmax * xln - 0.5D+00 )

        if ( abs ( xmax - xold ) < 0.005D+00 ) then
          xmax = xmax - 0.01D+00
          xmin = max ( xmin, - xmax + 1.0D+00 )
          return
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAML - Fatal error!'
      write ( *, '(a)' ) '  Unable to find XMAX.'
      stop

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_GAML - Fatal error!'
  write ( *, '(a)' ) '  Unable to find XMIN.'

  stop
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates the gamma function of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the gamma function of X.
!
  implicit none

  real ( kind = 8 ) dxrel
  real ( kind = 8 ) gcs(42)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ngcs
  real ( kind = 8 ) pi
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_gamma
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_lgmc
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sinpiy
  real ( kind = 8 ) sq2pil
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save dxrel
  save gcs
  save ngcs
  save pi
  save sq2pil
  save xmax
  save xmin
  save xsml

  data gcs(  1) / +0.8571195590989331421920062399942D-02 /
  data gcs(  2) / +0.4415381324841006757191315771652D-02 /
  data gcs(  3) / +0.5685043681599363378632664588789D-01 /
  data gcs(  4) / -0.4219835396418560501012500186624D-02 /
  data gcs(  5) / +0.1326808181212460220584006796352D-02 /
  data gcs(  6) / -0.1893024529798880432523947023886D-03 /
  data gcs(  7) / +0.3606925327441245256578082217225D-04 /
  data gcs(  8) / -0.6056761904460864218485548290365D-05 /
  data gcs(  9) / +0.1055829546302283344731823509093D-05 /
  data gcs( 10) / -0.1811967365542384048291855891166D-06 /
  data gcs( 11) / +0.3117724964715322277790254593169D-07 /
  data gcs( 12) / -0.5354219639019687140874081024347D-08 /
  data gcs( 13) / +0.9193275519859588946887786825940D-09 /
  data gcs( 14) / -0.1577941280288339761767423273953D-09 /
  data gcs( 15) / +0.2707980622934954543266540433089D-10 /
  data gcs( 16) / -0.4646818653825730144081661058933D-11 /
  data gcs( 17) / +0.7973350192007419656460767175359D-12 /
  data gcs( 18) / -0.1368078209830916025799499172309D-12 /
  data gcs( 19) / +0.2347319486563800657233471771688D-13 /
  data gcs( 20) / -0.4027432614949066932766570534699D-14 /
  data gcs( 21) / +0.6910051747372100912138336975257D-15 /
  data gcs( 22) / -0.1185584500221992907052387126192D-15 /
  data gcs( 23) / +0.2034148542496373955201026051932D-16 /
  data gcs( 24) / -0.3490054341717405849274012949108D-17 /
  data gcs( 25) / +0.5987993856485305567135051066026D-18 /
  data gcs( 26) / -0.1027378057872228074490069778431D-18 /
  data gcs( 27) / +0.1762702816060529824942759660748D-19 /
  data gcs( 28) / -0.3024320653735306260958772112042D-20 /
  data gcs( 29) / +0.5188914660218397839717833550506D-21 /
  data gcs( 30) / -0.8902770842456576692449251601066D-22 /
  data gcs( 31) / +0.1527474068493342602274596891306D-22 /
  data gcs( 32) / -0.2620731256187362900257328332799D-23 /
  data gcs( 33) / +0.4496464047830538670331046570666D-24 /
  data gcs( 34) / -0.7714712731336877911703901525333D-25 /
  data gcs( 35) / +0.1323635453126044036486572714666D-25 /
  data gcs( 36) / -0.2270999412942928816702313813333D-26 /
  data gcs( 37) / +0.3896418998003991449320816639999D-27 /
  data gcs( 38) / -0.6685198115125953327792127999999D-28 /
  data gcs( 39) / +0.1146998663140024384347613866666D-28 /
  data gcs( 40) / -0.1967938586345134677295103999999D-29 /
  data gcs( 41) / +0.3376448816585338090334890666666D-30 /
  data gcs( 42) / -0.5793070335782135784625493333333D-31 /

  data dxrel / 0.0D+00 /
  data ngcs / 0 /
  data pi / 3.14159265358979323846264338327950D+00 /
  data sq2pil / 0.91893853320467274178032973640562D+00 /
  data xmax / 0.0D+00 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ngcs == 0 ) then
    ngcs = r8_inits ( gcs, 42, 0.1D+00 * r8_mach ( 3 ) )
    call r8_gaml ( xmin, xmax )
    xsml = exp ( max ( log ( r8_mach ( 1 ) ), &
      - log ( r8_mach ( 2 ) ) ) + 0.01D+00 )
    dxrel = sqrt ( r8_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( y <= 10.0D+00 ) then

    n = int ( x )
    if ( x < 0.0D+00 ) then
      n = n - 1
    end if
    y = x - real ( n, kind = 8 )
    n = n - 1
    r8_gamma = 0.9375D+00 + r8_csevl ( 2.0D+00 * y - 1.0D+00, gcs, ngcs )

    if ( n == 0 ) then

      return

    else if ( n < 0 ) then

      n = - n

      if ( x == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
        write ( *, '(a)' ) '  X is 0.'
        stop
      end if

      if ( x < 0.0D+00 .and. x + real ( n - 2, kind = 8 ) == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
        write ( *, '(a)' ) '  X is a negative integer.'
        stop
      end if

      if ( x < - 0.5D+00 .and. &
        abs ( ( x - aint ( x - 0.5D+00 ) ) / x ) < dxrel ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_GAMMA - Warning!'
        write ( *, '(a)' ) '  X too near a negative integer,'
        write ( *, '(a)' ) '  answer is half precision.'
      end if

      if ( y < xsml ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
        write ( *, '(a)' ) '  X is so close to zero that Gamma overflows.'
        stop
      end if

      do i = 1, n
        r8_gamma = r8_gamma / ( x + real ( i - 1, kind = 8 ) )
      end do

    else if ( n == 0 ) then

    else

      do i = 1, n
        r8_gamma = ( y + real ( i, kind = 8 ) ) * r8_gamma
      end do

    end if

  else

    if ( xmax < x ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
      write ( *, '(a)' ) '  X so big that Gamma overflows.'
      stop
    end if
!
!  Underflow.
!
    if ( x < xmin ) then
      r8_gamma = 0.0D+00
      return
    end if

    r8_gamma = exp ( ( y - 0.5D+00 ) * log ( y ) - y + sq2pil + r8_lgmc ( y ) )

    if ( 0.0D+00 < x ) then
      return
    end if

    if ( abs ( ( x - aint ( x - 0.5D+00 ) ) / x ) < dxrel ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAMMA - Warning!'
      write ( *, '(a)' ) '  X too near a negative integer,'
      write ( *, '(a)' ) '  answer is half precision.'
    end if

    sinpiy = sin ( pi * y )

    if ( sinpiy == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
      write ( *, '(a)' ) '  X is a negative integer.'
      stop
    end if

    r8_gamma = - pi / ( y * sinpiy * r8_gamma )

  end if

  return
end
function r8_gamr ( x )

!*****************************************************************************80
!
!! R8_GAMR evaluates the reciprocal gamma function of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_GAMR, the value of the reciprocal gamma
!    function at X.
!
  implicit none

  real ( kind = 8 ) alngx
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_gamr
  real ( kind = 8 ) sgngx
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 .and. aint ( x ) == x ) then

    r8_gamr = 0.0D+00
    
  else if ( abs ( x ) <= 10.0D+00 ) then

    r8_gamr = 1.0D+00 / r8_gamma ( x )

  else

    call r8_lgams ( x, alngx, sgngx )
    r8_gamr = sgngx * exp ( - alngx )

  end if

  return
end
function r8_gmic ( a, x, alx )

!*****************************************************************************80
!
!! R8_GMIC: complementary incomplete gamma, small X, A near negative integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, real ( kind = 8 ) ALX, the logarithm of X.
!
!    Output, real ( kind = 8 ) R8_GMIC, the complementary incomplete 
!    gamma function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alng
  real ( kind = 8 ) alx
  real ( kind = 8 ) bot
  logical converged
  real ( kind = 8 ) eps
  real ( kind = 8 ) euler
  real ( kind = 8 ) fk
  real ( kind = 8 ) fkp1
  real ( kind = 8 ) fm
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ma
  integer ( kind = 4 ) mm1
  real ( kind = 8 ) r8_gmic
  real ( kind = 8 ) r8_lngam
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) s
  real ( kind = 8 ) sgng
  real ( kind = 8 ) t
  real ( kind = 8 ) te
  real ( kind = 8 ) x

  save bot
  save eps
  save euler

  data bot / 0.0D+00 /
  data eps / 0.0D+00 /
  data euler / 0.57721566490153286060651209008240D+00 /

  if ( eps == 0.0D+00 ) then
    eps = 0.5D+00 * r8_mach ( 3 )
    bot = log ( r8_mach ( 1 ) )
  end if

  if ( 0.0D+00 < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_GMIC - Fatal error!'
    write ( *, '(a)' ) '  A must be near a negative integer.'
    stop
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_GMIC - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  m = - int ( a - 0.5D+00 )
  fm = real ( m, kind = 8 )

  te = 1.0D+00
  t = 1.0D+00
  s = t
  converged = .false.
  do k = 1, 200
    fkp1 = real ( k + 1, kind = 8 )
    te = - x * te / ( fm + fkp1 )
    t = te / fkp1
    s = s + t
    if ( abs ( t ) < eps * s ) then
      converged = .true.
      exit
    end if
  end do

  if ( .not. converged ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_GMIC - Fatal error!'
    write ( *, '(a)' ) '  No convergence after 200 iterations.'
    stop
  end if

  r8_gmic = - alx - euler + x * s / ( fm + 1.0D+00 )

  if ( m == 0 ) then
    return
  else if ( m == 1 ) then
    r8_gmic = - r8_gmic - 1.0D+00 + 1.0D+00 / x
    return
  end if

  te = fm
  t = 1.0D+00
  s = t
  mm1 = m - 1
  do k = 1, mm1
    fk = real ( k, kind = 8 )
    te = - x * te / fk
    t = te / ( fm - fk )
    s = s + t
    if ( abs ( t ) < eps * abs ( s ) ) then
      exit
    end if
  end do

  do k = 1, m
    r8_gmic = r8_gmic + 1.0D+00 / real ( k, kind = 8 )
  end do

  if ( mod ( m, 2 ) == 1 ) then
    sgng = - 1.0D+00
  else
    sgng = + 1.0D+00
  end if

  alng = log ( r8_gmic ) - r8_lngam ( fm + 1.0D+00 )

  if ( bot < alng ) then
    r8_gmic = sgng * exp ( alng )
  else
    r8_gmic = 0.0D+00
  end if

  if ( s /= 0.0D+00 ) then
    r8_gmic = r8_gmic &
      + sign ( exp ( - fm * alx + log ( abs ( s ) / fm ) ), s )
  end if

  return
end
function r8_gmit ( a, x, algap1, sgngam, alx )

!*****************************************************************************80
!
!! R8_GMIT: Tricomi's incomplete gamma function for small X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, real ( kind = 8 ) ALGAP1, the logarithm of Gamma ( A + 1 ).
!
!    Input, real ( kind = 8 ) SGNGAM, the sign of Gamma ( A + 1 ).
!
!    Input, real ( kind = 8 ) ALX, the logarithm of X.
!
!    Output, real ( kind = 8 ) R8_GMIT, the Tricomi incomplete gamma function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ae
  real ( kind = 8 ) aeps
  real ( kind = 8 ) alg2
  real ( kind = 8 ) algap1
  real ( kind = 8 ) algs
  real ( kind = 8 ) alx
  real ( kind = 8 ) bot
  logical converged
  real ( kind = 8 ) eps
  real ( kind = 8 ) fk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ma
  real ( kind = 8 ) r8_gmit
  real ( kind = 8 ) r8_lngam
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) s
  real ( kind = 8 ) sgng2
  real ( kind = 8 ) sgngam
  real ( kind = 8 ) t
  real ( kind = 8 ) te
  real ( kind = 8 ) x

  save bot
  save eps

  data bot / 0.0D+00 /
  data eps / 0.0D+00 /

  if ( eps == 0.0D+00 ) then
    eps = 0.5D+00 * r8_mach ( 3 )
    bot = log ( r8_mach ( 1 ) )
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_GMIT - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  if ( a < 0.0D+00 ) then
    ma = int ( a - 0.5D+00 )
  else
    ma = int ( a + 0.5D+00 )
  end if

  aeps = a - real ( ma, kind = 8 )

  if ( a < - 0.5D+00 ) then
    ae = aeps
  else
    ae = a
  end if

  t = 1.0D+00
  te = ae
  s = t
  converged = .false.
  do k = 1, 200
    fk = real ( k, kind = 8 )
    te = - x * te / fk
    t = te / ( ae + fk )
    s = s + t
    if ( abs ( t ) < eps * abs ( s ) ) then
      converged = .true.
      exit
    end if
  end do

  if ( .not. converged ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_GMIT - Fatal error!'
    write ( *, '(a)' ) '  No convergence in 200 iterations.'
    stop
  end if

  if ( - 0.5D+00 <= a ) then
    algs = - algap1 + log ( s )
    r8_gmit = exp ( algs )
    return
  end if

  algs = - r8_lngam ( 1.0D+00 + aeps ) + log ( s )
  s = 1.0D+00
  m = - ma - 1
  t = 1.0D+00
  do k = 1, m
    t = x * t / ( aeps - real ( m + 1 - k, kind = 8 ) )
    s = s + t
    if ( abs ( t ) < eps * abs ( s ) ) then
      exit
    end if
  end do

  r8_gmit = 0.0D+00
  algs = - real ( ma, kind = 8 ) * log ( x ) + algs

  if ( s == 0.0D+00 .or. aeps == 0.0D+00 ) then
    r8_gmit = exp ( algs )
    return
  end if

  sgng2 = sgngam * r8_sign ( s )
  alg2 = - x - algap1 + log ( abs ( s ) )

  if ( bot < alg2 ) then
    r8_gmit = sgng2 * exp ( alg2 )
  end if

  if ( bot < algs ) then
    r8_gmit = r8_gmit + exp ( algs )
  end if

  return
end
function r8_inits ( dos, nos, eta )

!*****************************************************************************80
!
!! R8_INITS initializes a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DOS(NOS), the Chebyshev coefficients.
!
!    Input, integer ( kind = 4 ) NOS, the number of coefficients.
!
!    Input, real ( kind = 8 ) ETA, the desired accuracy.
!
!    Output, integer ( kind = 4 ) R8_INITS, the number of terms of the 
!    series needed to ensure the requested accuracy.
!
  implicit none

  integer ( kind = 4 ) nos

  real ( kind = 8 ) dos(nos)
  real ( kind = 8 ) err 
  real ( kind = 8 ) eta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) r8_inits

  if ( nos < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_INITS - Fatal error!'
    write ( *, '(a)' ) '  Number of coefficients < 1.'
    stop
  end if

  err = 0.0D+00

  do i = nos, 1, -1
    err = err + abs ( dos(i) )
    if ( eta < err ) then
      r8_inits = i
      return
    end if
  end do

  r8_inits = nos
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_INITS - Warning!'
  write ( *, '(a)' ) '  ETA may be too small.'

  return
end
function r8_int ( x )

!*****************************************************************************80
!
!! R8_INT returns the integer part of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_INT, the integer part of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_mach
  integer ( kind = 4 ) ibase
  integer ( kind = 4 ) ipart
  integer ( kind = 4 ), save :: npart = 0
  real ( kind = 8 ) part
  real ( kind = 8 ) r8_int
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ), save :: scale = 0.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xbig = 0.0D+00
  real ( kind = 8 ), save :: xmax = 0.0D+00
  real ( kind = 8 ) xscl

  if ( npart == 0 ) then
    ibase = i4_mach ( 10 )
    xmax = 1.0D+00 / r8_mach ( 4 )
    xbig = min ( real ( i4_mach ( 9 ), kind = 8 ), 1.0D+00 / r8_mach ( 4 ) )
    scale = ibase ** int ( log ( xbig ) &
      / log ( real ( ibase, kind = 8 ) ) - 0.5D+00 )
    npart = log ( xmax ) / log ( scale ) + 1.0D+00
  end if
!
!  X may be too small.
!
  if ( x < - xmax ) then

    r8_int = x

  else if ( x < - xbig ) then

    xscl = - x

    do i = 1, npart
      xscl = xscl / scale
    end do

    r8_int = 0.0D+00
    do i = 1, npart
      xscl = xscl * scale
      ipart = int ( xscl )
      part = real ( ipart, kind = 8 )
      xscl = xscl - part
      r8_int = r8_int * scale + part
    end do

    r8_int = - r8_int

  else if ( x <= xbig ) then

    r8_int = int ( x )

  else if ( x <= xmax ) then

    xscl = x

    do i = 1, npart
      xscl = xscl / scale
    end do

    r8_int = 0.0D+00
    do i = 1, npart
      xscl = xscl * scale
      ipart = int ( xscl )
      part = real ( ipart, kind = 8 )
      xscl = xscl - part
      r8_int = r8_int * scale + part
    end do
!
!  X may be too large.
!
  else
    r8_int = x
  end if

  return
end
subroutine r8_knus ( xnu, x, bknu, bknu1, iswtch )

!*****************************************************************************80
!
!! R8_KNUS computes a sequence of K Bessel functions.
!
!  Discussion:
!
!    This routine computes Bessel functions 
!      exp(x) * k-sub-xnu (x)  
!    and
!      exp(x) * k-sub-xnu+1 (x) 
!    for 0.0 <= xnu < 1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XNU, the order parameter.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) BKNU, BKNU1, the two K Bessel functions.
!
!    Output, integer ( kind = 4 ) ISWTCH, ?
!
  implicit none

  real ( kind = 8 ) a(32)
  real ( kind = 8 ) a0
  real ( kind = 8 ) aln2
  real ( kind = 8 ) alnbig
  real ( kind = 8 ) alneps
  real ( kind = 8 ) alnsml
  real ( kind = 8 ) alnz
  real ( kind = 8 ) alpha(32)
  real ( kind = 8 ) an
  real ( kind = 8 ) b0
  real ( kind = 8 ) beta(32)
  real ( kind = 8 ) bknu
  real ( kind = 8 ) bknu0
  real ( kind = 8 ) bknu1
  real ( kind = 8 ) bknud
  real ( kind = 8 ) bn
  real ( kind = 8 ) c0
  real ( kind = 8 ) c0kcs(29)
  real ( kind = 8 ) eta
  real ( kind = 8 ) euler
  real ( kind = 8 ) expx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inu
  integer ( kind = 4 ) iswtch
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntc0k
  integer ( kind = 4 ) nterms
  integer ( kind = 4 ) ntznu1
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) qq
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_gamma
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) result
  real ( kind = 8 ) sqpi2
  real ( kind = 8 ) sqrtx
  real ( kind = 8 ) v
  real ( kind = 8 ) vlnz
  real ( kind = 8 ) x
  real ( kind = 8 ) x2n
  real ( kind = 8 ) x2tov
  real ( kind = 8 ) xi
  real ( kind = 8 ) xmu
  real ( kind = 8 ) xnu
  real ( kind = 8 ) xnusml
  real ( kind = 8 ) xsml
  real ( kind = 8 ) z
  real ( kind = 8 ) znu1cs(20)
  real ( kind = 8 ) ztov

  save aln2
  save alnbig
  save alneps
  save alnsml
  save c0kcs
  save euler
  save ntc0k
  save ntznu1
  save sqpi2
  save xnusml
  save xsml
  save znu1cs

  data c0kcs(  1) / +0.60183057242626108387577445180329D-01     /
  data c0kcs(  2) / -0.15364871433017286092959755943124D+00     /
  data c0kcs(  3) / -0.11751176008210492040068229226213D-01     /
  data c0kcs(  4) / -0.85248788891979509827048401550987D-03     /
  data c0kcs(  5) / -0.61329838767496791874098176922111D-04     /
  data c0kcs(  6) / -0.44052281245510444562679889548505D-05     /
  data c0kcs(  7) / -0.31631246728384488192915445892199D-06     /
  data c0kcs(  8) / -0.22710719382899588330673771793396D-07     /
  data c0kcs(  9) / -0.16305644608077609552274620515360D-08     /
  data c0kcs( 10) / -0.11706939299414776568756044043130D-09     /
  data c0kcs( 11) / -0.84052063786464437174546593413792D-11    /
  data c0kcs( 12) / -0.60346670118979991487096050737198D-12    /
  data c0kcs( 13) / -0.43326960335681371952045997366903D-13    /
  data c0kcs( 14) / -0.31107358030203546214634697772237D-14    /
  data c0kcs( 15) / -0.22334078226736982254486133409840D-15    /
  data c0kcs( 16) / -0.16035146716864226300635791528610D-16    /
  data c0kcs( 17) / -0.11512717363666556196035697705305D-17    /
  data c0kcs( 18) / -0.82657591746836959105169479089258D-19    /
  data c0kcs( 19) / -0.59345480806383948172333436695984D-20    /
  data c0kcs( 20) / -0.42608138196467143926499613023976D-21    /
  data c0kcs( 21) / -0.30591266864812876299263698370542D-22    /
  data c0kcs( 22) / -0.21963541426734575224975501815516D-23    /
  data c0kcs( 23) / -0.15769113261495836071105750684760D-24    /
  data c0kcs( 24) / -0.11321713935950320948757731048056D-25    /
  data c0kcs( 25) / -0.81286248834598404082792349714433D-27    /
  data c0kcs( 26) / -0.58360900893453226552829349315949D-28    /
  data c0kcs( 27) / -0.41901241623610922519452337780905D-29    /
  data c0kcs( 28) / -0.30083737960206435069530504212862D-30    /
  data c0kcs( 29) / -0.21599152067808647728342168089832D-31    /

  data znu1cs(  1) / +0.203306756994191729674444001216911D+00    /
  data znu1cs(  2) / +0.140077933413219771062943670790563D+00    /
  data znu1cs(  3) / +0.791679696100161352840972241972320D-02    /
  data znu1cs(  4) / +0.339801182532104045352930092205750D-03    /
  data znu1cs(  5) / +0.117419756889893366664507228352690D-04    /
  data znu1cs(  6) / +0.339357570612261680333825865475121D-06    /
  data znu1cs(  7) / +0.842594176976219910194629891264803D-08    /
  data znu1cs(  8) / +0.183336677024850089184748150900090D-09    /
  data znu1cs(  9) / +0.354969844704416310863007064469557D-11   /
  data znu1cs( 10) / +0.619032496469887332205244342078407D-13   /
  data znu1cs( 11) / +0.981964535680439424960346115456527D-15   /
  data znu1cs( 12) / +0.142851314396490474211473563005985D-16   /
  data znu1cs( 13) / +0.191894921887825298966162467488436D-18   /
  data znu1cs( 14) / +0.239430979739498914162313140597128D-20   /
  data znu1cs( 15) / +0.278890246815347354835870465474995D-22   /
  data znu1cs( 16) / +0.304606650633033442582845214092865D-24   /
  data znu1cs( 17) / +0.313173237042191815771564260932089D-26   /
  data znu1cs( 18) / +0.304133098987854951645174908005034D-28   /
  data znu1cs( 19) / +0.279840384636833084343185097659733D-30   /
  data znu1cs( 20) / +0.244637186274497596485238794922666D-32   /

  data aln2 / 0.69314718055994530941723212145818D+00 /
  data alnbig / 0.0D+00 /
  data alneps / 0.0D+00 /
  data alnsml / 0.0D+00 /
  data euler / 0.57721566490153286060651209008240D+00 /
  data ntc0k / 0 /
  data ntznu1 / 0 /
  data sqpi2 / +1.2533141373155002512078826424055D+00 /
  data xnusml / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ntc0k == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    ntc0k = r8_inits ( c0kcs, 29, eta )
    ntznu1 = r8_inits ( znu1cs, 20, eta )
    xnusml = sqrt ( r8_mach ( 3 ) / 8.0D+00 )
    xsml = 0.1D+00 * r8_mach ( 3 )
    alnsml = log ( r8_mach ( 1 ) )
    alnbig = log ( r8_mach ( 2 ) )
    alneps = log ( 0.1D+00 * r8_mach ( 3 ) )
  end if

  if ( xnu < 0.0D+00 .or. 1.0D+00 <= xnu ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_KNUS - Fatal error!'
    write ( *, '(a)' ) '  XNU < 0 or. 1 <= XNU.'
    stop
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_KNUS - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  iswtch = 0
!
!  X is small.  Compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
!  then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-0.5,+0.5)
!  then to (0., .5), because k of negative order (-nu) = k of positive
!  order (+nu).
!
  if ( x <= 2.0D+00 ) then

    if ( xnu <= 0.5D+00 ) then
      v = xnu
    else
      v = 1.0D+00 - xnu
    end if
!
!  Carefully find (x/2)^xnu and z^xnu where z = x*x/4.
!
    alnz = 2.0D+00 * ( log ( x ) - aln2 )

    if ( x <= xnu ) then

      if ( alnbig < - 0.5D+00 * xnu * alnz - aln2 - log ( xnu ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_KNUS - Fatal error!'
        write ( *, '(a)' ) '  Small X causing overflow.'
        stop
      end if

    end if

    vlnz = v * alnz
    x2tov = exp ( 0.5D+00 * vlnz )

    if ( vlnz <= alnsml ) then
      ztov = 0.0D+00
    else
      ztov = x2tov * x2tov
    end if

    a0 = 0.5D+00 * r8_gamma ( 1.0D+00 + v )
    b0 = 0.5D+00 * r8_gamma ( 1.0D+00 - v )
    c0 = - euler
    if ( 0.5D+00 <= ztov .and. xnusml < v ) then
      c0 = - 0.75D+00 + &
        r8_csevl ( ( 8.0D+00 * v ) * v - 1.0D+00, c0kcs, ntc0k )
    end if

    if ( ztov <= 0.5D+00 ) then
      alpha(1) = ( a0 - ztov * b0 ) / v
    else
      alpha(1) = c0 - alnz * ( 0.75D+00 + &
        r8_csevl ( vlnz / 0.35D+00 + 1.0D+00, znu1cs, ntznu1 ) ) * b0
    end if

    beta(1) = - 0.5D+00 * ( a0 + ztov * b0 )

    if ( x <= xsml ) then
      z = 0.0D+00
    else
      z = 0.25D+00 * x * x
    end if

    nterms = max ( 2, int ( 11.0D+00 &
      + ( 8.0D+00 * alnz - 25.19D+00 - alneps ) &
      / ( 4.28D+00 - alnz ) ) )

    do i = 2, nterms
      xi = real ( i - 1, kind = 8 )
      a0 = a0 / ( xi * ( xi - v ) )
      b0 = b0 / ( xi * ( xi + v ) )
      alpha(i) = ( alpha(i-1) + 2.0D+00 * xi * a0 ) / ( xi * ( xi + v ) )
      beta(i) = ( xi - 0.5D+00 * v ) * alpha(i) - ztov * b0
    end do

    bknu = alpha(nterms)
    bknud = beta(nterms)
    do ii = 2, nterms
      i = nterms + 1 - ii
      bknu = alpha(i) + bknu * z
      bknud = beta(i) + bknud * z
    end do

    expx = exp ( x )
    bknu = expx * bknu / x2tov

    if ( alnbig < - 0.5D+00 * ( xnu + 1.0D+00 ) * alnz - 2.0D+00 * aln2 ) then
      iswtch = 1
      return
    end if

    bknud = expx * bknud * 2.0D+00 / ( x2tov * x )

    if ( xnu <= 0.5D+00 ) then
      bknu1 = v * bknu / x - bknud
      return
    end if

    bknu0 = bknu
    bknu = - v * bknu / x - bknud
    bknu1 = 2.0D+00 * xnu * bknu / x + bknu0
!
!  X is large.  Find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke-s
!  rational expansion.
!
  else

    sqrtx = sqrt ( x )

    if ( 1.0D+00 / xsml < x ) then
      bknu = sqpi2 / sqrtx
      bknu1 = bknu
      return
    end if

    an = - 0.60D+00 - 1.02D+00 / x
    bn = - 0.27D+00 - 0.53D+00 / x
    nterms = min ( 32, max ( 3, int ( an + bn * alneps ) ) )

    do inu = 1, 2

      if ( inu == 1 ) then
        if ( xnu <= xnusml ) then
          xmu = 0.0D+00
        else
          xmu = ( 4.0D+00 * xnu ) * xnu
        end if
      else
        xmu = 4.0D+00 * ( abs ( xnu ) + 1.0D+00 ) ** 2
      end if

      a(1) = 1.0D+00 - xmu
      a(2) = 9.0D+00 - xmu
      a(3) = 25.0D+00 - xmu

      if ( a(2) == 0.0D+00 ) then

        result = sqpi2 * ( 16.0D+00 * x + xmu + 7.0D+00 ) &
          / ( 16.0D+00 * x * sqrtx )

      else

        alpha(1) = 1.0D+00
        alpha(2) = ( 16.0D+00 * x + a(2) ) / a(2)
        alpha(3) = ( ( 768.0D+00 * x + 48.0D+00 * a(3) ) * x &
          + a(2) * a(3) ) / ( a(2) * a(3) )

        beta(1) = 1.0D+00
        beta(2) = ( 16.0D+00 * x + ( xmu + 7.0D+00 ) ) / a(2)
        beta(3) = ( ( 768.0D+00 * x &
          + 48.0D+00 * ( xmu + 23.0D+00 ) ) * x + &
          ( ( xmu + 62.0D+00 ) * xmu + 129.0D+00 ) ) &
          / ( a(2) * a(3) )

        do i = 4, nterms

          n = i - 1
          x2n = real ( 2 * n - 1, kind = 8 )

          a(i) = ( x2n + 2.0D+00 ) ** 2 - xmu
          qq = 16.0D+00 * x2n / a(i)
          p1 = - x2n * ( real ( 12 * n * n - 20 * n, kind = 8 ) - a(1) ) &
            / ( ( x2n - 2.0D+00 ) * a(i) ) - qq * x
          p2 = ( real ( 12 * n * n - 28 * n + 8, kind = 8 ) - a(1) ) / a(i) &
            - qq * x
          p3 = - x2n * a(i-3) / ( ( x2n - 2.0D+00 ) * a(i))

          alpha(i) = - p1 * alpha(i-1) &
                     - p2 * alpha(i-2) &
                     - p3 * alpha(i-3)

          beta(i) =  - p1 * beta(i-1) &
                     - p2 * beta(i-2) &
                     - p3 * beta(i-3)

        end do

        result = sqpi2 * beta(nterms) / ( sqrtx * alpha(nterms) )

      end if

      if ( inu == 1 ) then
        bknu = result
      else
        bknu1 = result
      end if

    end do

  end if

  return
end
function r8_lbeta ( a, b )

!*****************************************************************************80
!
!! R8_LBETA evaluates the logarithm of the beta function of R8 arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the arguments.
!
!    Output, real ( kind = 8 ) R8_LBETA, the logarithm of the beta function of A
!    and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) corr
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_lbeta
  real ( kind = 8 ) r8_lgmc
  real ( kind = 8 ) r8_lngam
  real ( kind = 8 ) r8_lnrel
  real ( kind = 8 ) sq2pil

  save sq2pil

  data sq2pil / 0.91893853320467274178032973640562D+00 /

  p = min ( a, b )
  q = max ( a, b )

  if ( p <= 0.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LBETA - Fatal error!'
    write ( *, '(a)' ) '  Both arguments must be greater than 0.'
    stop

  else if ( p < 10.0D+00 .and. q <= 10.0D+00 ) then

    r8_lbeta = log ( r8_gamma ( p ) &
      * ( r8_gamma ( q ) / r8_gamma ( p + q ) ) )

  else if ( p < 10.0D+00 ) then

    corr = r8_lgmc ( q ) - r8_lgmc ( p + q )

    r8_lbeta = r8_lngam ( p ) + corr + p - p * log ( p + q ) + &
      ( q - 0.5D+00 ) * r8_lnrel ( - p / ( p + q ) )

  else

    corr = r8_lgmc ( p ) + r8_lgmc ( q ) - r8_lgmc ( p + q )

    r8_lbeta = - 0.5D+00 * log ( q ) + sq2pil + corr &
      + ( p - 0.5D+00 ) * log ( p / ( p + q ) ) &
      + q * r8_lnrel ( - p / ( p + q ) )

  end if

  return
end
subroutine r8_lgams ( x, algam, sgngam )

!*****************************************************************************80
!
!! R8_LGAMS evaluates the log of |gamma(x)| and sign, for an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) ALGAM, the logarithm of the absolute value of
!    gamma ( X ).
!
!    Output, real ( kind = 8 ) SGNGAM, the sign (+1 or -1 ) of gamma ( X ).
!
  implicit none

  real ( kind = 8 ) algam
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_lngam
  real ( kind = 8 ) sgngam
  real ( kind = 8 ) x

  algam = r8_lngam ( x )
  sgngam = 1.0D+00

  if ( x <= 0.0D+00 ) then

    k = int ( mod ( - aint ( x ), 2.0D+00 ) + 0.1D+00 )

    if ( k == 0 ) then
      sgngam = - 1.0D+00
    end if

  end if

  return
end
function r8_lgic ( a, x, alx )

!*****************************************************************************80
!
!! R8_LGIC: the log complementary incomplete gamma function for large X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, real ( kind = 8 ) ALX, the logarithm of X.
!
!    Output, real ( kind = 8 ) R8_LGIC, the log complementary incomplete 
!    gamma function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alx
  real ( kind = 8 ) eps
  real ( kind = 8 ) fk
  integer ( kind = 4 ) k
  real ( kind = 8 ) p
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_lgic
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) xma
  real ( kind = 8 ) xpa

  save eps

  data eps / 0.0D+00 /

  if ( eps == 0.0D+00 ) then
    eps = 0.5D+00 * r8_mach ( 3 )
  end if

  xpa = x + 1.0D+00 - a
  xma = x - 1.0D+00 - a

  r = 0.0D+00
  p = 1.0D+00
  s = p
  do k = 1, 300
    fk = real ( k, kind = 8 )
    t = fk * ( a - fk ) * ( 1.0D+00 + r )
    r = - t / ( ( xma + 2.0D+00 * fk ) * ( xpa + 2.0D+00 * fk ) + t )
    p = r * p
    s = s + p
    if ( abs ( p ) < eps * s ) then
      r8_lgic = a * alx - x + log ( s / xpa )
      return
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_LGIC - Fatal error!'
  write ( *, '(a)' ) '  No convergence in 300 iterations.'

  stop
end
function r8_lgit ( a, x, algap1 )

!*****************************************************************************80
!
!! R8_LGIT evaluates the log of Tricomi's incomplete gamma function.
!
!  Discussion:
!
!    Perron's continued fraction is used for large X and X <= A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, real ( kind = 8 ) ALGAP1, the logarithm of A+1.
!
!    Output, real ( kind = 8 ) R8_LGIT, the log of Tricomi's incomplete
!    gamma function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a1x
  real ( kind = 8 ) algap1
  real ( kind = 8 ) ax
  real ( kind = 8 ) eps
  real ( kind = 8 ) fk
  real ( kind = 8 ) hstar
  integer ( kind = 4 ) k
  real ( kind = 8 ) p
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_lgit
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) s
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  save eps

  data eps / 0.0D+00 /

  if ( eps == 0.0D+00 ) then
    eps = 0.5D+00 * r8_mach ( 3 )
    sqeps = sqrt ( r8_mach ( 4 ) )
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LGIT - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  if ( a < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LGIT - Fatal error!'
    write ( *, '(a)' ) '  A < X.'
    stop
  end if

  ax = a + x
  a1x = ax + 1.0D+00
  r = 0.0D+00
  p = 1.0D+00
  s = p
  do k = 1, 200
    fk = real ( k, kind = 8 )
    t = ( a + fk ) * x * ( 1.0D+00 + r )
    r = t / ( ( ax + fk ) * ( a1x + fk ) - t )
    p = r * p
    s = s + p
    if ( abs ( p ) < eps * s ) then
      hstar = 1.0D+00 - x * s / a1x
      r8_lgit = - x - algap1 - log ( hstar )
      return
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_LGIT - Fatal error!'
  write ( *, '(a)' ) '  No convergence after 200 iterations.'
  stop
end
function r8_lgmc ( x )

!*****************************************************************************80
!
!! R8_LGMC evaluates the log gamma correction factor for an R8 argument.
!
!  Discussion:
!
!    For 10 <= X, compute the log gamma correction factor so that
!
!      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
!                          + ( x - 0.5 ) * log ( x ) - x 
!                          + r8_lgmc ( x )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_LGMC, the correction factor.
!
  implicit none

  real ( kind = 8 ) algmcs(15)
  integer ( kind = 4 ) nalgm
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_lgmc
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xmax

  save algmcs
  save nalgm
  save xbig
  save xmax

  data algmcs(  1) / +0.1666389480451863247205729650822D+00 /
  data algmcs(  2) / -0.1384948176067563840732986059135D-04 /
  data algmcs(  3) / +0.9810825646924729426157171547487D-08 /
  data algmcs(  4) / -0.1809129475572494194263306266719D-10 /
  data algmcs(  5) / +0.6221098041892605227126015543416D-13 /
  data algmcs(  6) / -0.3399615005417721944303330599666D-15 /
  data algmcs(  7) / +0.2683181998482698748957538846666D-17 /
  data algmcs(  8) / -0.2868042435334643284144622399999D-19 /
  data algmcs(  9) / +0.3962837061046434803679306666666D-21 /
  data algmcs( 10) / -0.6831888753985766870111999999999D-23 /
  data algmcs( 11) / +0.1429227355942498147573333333333D-24 /
  data algmcs( 12) / -0.3547598158101070547199999999999D-26 /
  data algmcs( 13) / +0.1025680058010470912000000000000D-27 /
  data algmcs( 14) / -0.3401102254316748799999999999999D-29 /
  data algmcs( 15) / +0.1276642195630062933333333333333D-30 /

  data nalgm / 0 /
  data xbig / 0.0D+00 /
  data xmax / 0.0D+00 /

  if ( nalgm == 0 ) then
    nalgm = r8_inits ( algmcs, 15, r8_mach ( 3 ) )
    xbig = 1.0D+00 / sqrt ( r8_mach ( 3 ) )
    xmax = exp ( min ( log ( r8_mach ( 2 ) / 12.0D+00 ), &
      - log ( 12.0D+00 * r8_mach ( 1 ) ) ) )
  end if

  if ( x < 10.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LGMC - Fatal error!'
    write ( *, '(a)' ) '  X must be at least 10.'
    stop

  else if ( x < xbig ) then

    r8_lgmc = r8_csevl ( 2.0D+00 * ( 10.0D+00 / x ) &
      * ( 10.0D+00 / x ) - 1.0D+00, algmcs, nalgm ) / x

  else if ( x < xmax ) then

    r8_lgmc = 1.0D+00 / ( 12.0D+00 * x )

  else

    r8_lgmc = 0.0D+00

  end if

  return
end
function r8_li ( x )

!*****************************************************************************80
!
!! R8_LI evaluates the logarithmic integral for an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_LI, the logarithmic integral evaluated at X.
!
  implicit none

  real ( kind = 8 ) r8_ei
  real ( kind = 8 ) r8_li
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) x

  if ( sqeps == 0.0D+00 ) then
    sqeps = sqrt ( r8_mach ( 3 ) )
  end if

  if ( x < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LI - Fatal error!'
    write ( *, '(a)' ) '  Function undefined for X <= 0.'
    stop
  end if

  if ( x == 0.0D+00 ) then
    r8_li = 0.0D+00
    return
  end if

  if ( x == 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LI - Fatal error!'
    write ( *, '(a)' ) '  Function undefined for X = 1.'
    stop
  end if

  if ( abs ( 1.0D+00 - x ) < sqeps ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LI - Warning!'
    write ( *, '(a)' ) '  Answer less than half precision.'
    write ( *, '(a)' ) '  X is too close to 1.'
  end if

  r8_li = r8_ei ( log ( x ) )

  return
end
function r8_lngam ( x )

!*****************************************************************************80
!
!! R8_LNGAM evaluates the log of the absolute value of gamma of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_LNGAM, the logarithm of the absolute value of
!    the gamma function of X.
!
  implicit none

  real ( kind = 8 ) dxrel
  real ( kind = 8 ) pi
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_lgmc
  real ( kind = 8 ) r8_lngam
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sinpiy
  real ( kind = 8 ) sq2pil
  real ( kind = 8 ) sqpi2l
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) y

  save dxrel
  save pi
  save sq2pil
  save sqpi2l
  save xmax

  data dxrel / 0.0D+00 /
  data pi / 3.14159265358979323846264338327950D+00 /
  data sq2pil / 0.91893853320467274178032973640562D+00 /
  data sqpi2l / +0.225791352644727432363097614947441D+00 /
  data xmax / 0.0D+00 /

  if ( xmax == 0.0D+00 ) then
    xmax = r8_mach ( 2 ) / log ( r8_mach ( 2 ) )
    dxrel = sqrt ( r8_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( y <= 10.0D+00 ) then
    r8_lngam = log ( abs ( r8_gamma ( x ) ) )
    return
  end if

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LNGAM - Fatal error!'
    write ( *, '(a)' ) '  Result overflows, |X| too big.'
    stop
  end if

  if ( 0.0D+00 < x ) then
    r8_lngam = sq2pil + ( x - 0.5D+00 ) * log ( x ) - x + r8_lgmc ( y )
    return
  end if

  sinpiy = abs ( sin ( pi * y ) )

  if ( sinpiy == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LNGAM - Fatal error!'
    write ( *, '(a)' ) '  X is a negative integer.'
    stop
  end if

  r8_lngam = sqpi2l + ( x - 0.5D+00 ) * log ( y ) - x - log ( sinpiy ) &
    - r8_lgmc ( y )

  if ( abs ( ( x - aint ( x - 0.5D+00 ) ) * r8_lngam / x ) < dxrel ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LNGAM - Warning!'
    write ( *, '(a)' ) '  Result is half precision because'
    write ( *, '(a)' ) '  X is too near a negative integer.'
  end if

  return
end
function r8_lnrel ( x )

!*****************************************************************************80
!
!! R8_LNREL evaluates log ( 1 + X ) for an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_LNREL, the value of LOG ( 1 + X ).
!
  implicit none

  real ( kind = 8 ) alnrcs(43)
  integer ( kind = 4 ) nlnrel
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_lnrel
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin

  save alnrcs
  save nlnrel
  save xmin

  data alnrcs(  1) / +0.10378693562743769800686267719098D+01 /
  data alnrcs(  2) / -0.13364301504908918098766041553133D+00 /
  data alnrcs(  3) / +0.19408249135520563357926199374750D-01 /
  data alnrcs(  4) / -0.30107551127535777690376537776592D-02 /
  data alnrcs(  5) / +0.48694614797154850090456366509137D-03 /
  data alnrcs(  6) / -0.81054881893175356066809943008622D-04 /
  data alnrcs(  7) / +0.13778847799559524782938251496059D-04 /
  data alnrcs(  8) / -0.23802210894358970251369992914935D-05 /
  data alnrcs(  9) / +0.41640416213865183476391859901989D-06 /
  data alnrcs( 10) / -0.73595828378075994984266837031998D-07 /
  data alnrcs( 11) / +0.13117611876241674949152294345011D-07 /
  data alnrcs( 12) / -0.23546709317742425136696092330175D-08 /
  data alnrcs( 13) / +0.42522773276034997775638052962567D-09 /
  data alnrcs( 14) / -0.77190894134840796826108107493300D-10 /
  data alnrcs( 15) / +0.14075746481359069909215356472191D-10 /
  data alnrcs( 16) / -0.25769072058024680627537078627584D-11 /
  data alnrcs( 17) / +0.47342406666294421849154395005938D-12 /
  data alnrcs( 18) / -0.87249012674742641745301263292675D-13 /
  data alnrcs( 19) / +0.16124614902740551465739833119115D-13 /
  data alnrcs( 20) / -0.29875652015665773006710792416815D-14 /
  data alnrcs( 21) / +0.55480701209082887983041321697279D-15 /
  data alnrcs( 22) / -0.10324619158271569595141333961932D-15 /
  data alnrcs( 23) / +0.19250239203049851177878503244868D-16 /
  data alnrcs( 24) / -0.35955073465265150011189707844266D-17 /
  data alnrcs( 25) / +0.67264542537876857892194574226773D-18 /
  data alnrcs( 26) / -0.12602624168735219252082425637546D-18 /
  data alnrcs( 27) / +0.23644884408606210044916158955519D-19 /
  data alnrcs( 28) / -0.44419377050807936898878389179733D-20 /
  data alnrcs( 29) / +0.83546594464034259016241293994666D-21 /
  data alnrcs( 30) / -0.15731559416479562574899253521066D-21 /
  data alnrcs( 31) / +0.29653128740247422686154369706666D-22 /
  data alnrcs( 32) / -0.55949583481815947292156013226666D-23 /
  data alnrcs( 33) / +0.10566354268835681048187284138666D-23 /
  data alnrcs( 34) / -0.19972483680670204548314999466666D-24 /
  data alnrcs( 35) / +0.37782977818839361421049855999999D-25 /
  data alnrcs( 36) / -0.71531586889081740345038165333333D-26 /
  data alnrcs( 37) / +0.13552488463674213646502024533333D-26 /
  data alnrcs( 38) / -0.25694673048487567430079829333333D-27 /
  data alnrcs( 39) / +0.48747756066216949076459519999999D-28 /
  data alnrcs( 40) / -0.92542112530849715321132373333333D-29 /
  data alnrcs( 41) / +0.17578597841760239233269760000000D-29 /
  data alnrcs( 42) / -0.33410026677731010351377066666666D-30 /
  data alnrcs( 43) / +0.63533936180236187354180266666666D-31 /

  data nlnrel / 0 /
  data xmin / 0.0D+00 /

  if ( nlnrel == 0 ) then
    nlnrel = r8_inits ( alnrcs, 43, 0.1D+00 * r8_mach ( 3 ) )
    xmin = - 1.0D+00 + sqrt ( r8_mach ( 4 ) )
  end if

  if ( x <= - 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LNREL - Fatal error!'
    write ( *, '(a)' ) '  X <= -1.'
    stop
  else if ( x < xmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LNREL - Warning!'
    write ( *, '(a)' ) '  Result is less than half precision.'
    write ( *, '(a)' ) '  X is too close to - 1.'
  end if

  if ( abs ( x ) <= 0.375D+00 ) then
    r8_lnrel = x * ( 1.0D+00 - x * r8_csevl ( x / 0.375D+00, alnrcs, nlnrel ) )
  else
    r8_lnrel = log ( 1.0D+00 + x )
  end if

  return
end
function r8_log ( x )

!*****************************************************************************80
!
!! R8_LOG evaluates the logarithm of an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) R8_LOG, the logarithm of X.
!
  implicit none

  real ( kind = 8 ) aln2
  real ( kind = 8 ) alncen(5)
  real ( kind = 8 ) alncs(11)
  real ( kind = 8 ) center(4)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nterms
  integer ( kind = 4 ) ntrval
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_log
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  real ( kind = 8 ) x
  real ( kind = 8 ) xn
  real ( kind = 8 ) y

  save aln2
  save alncen
  save alncs
  save center
  save nterms

  data alncs(  1) / +0.13347199877973881561689386047187D+01 /
  data alncs(  2) / +0.69375628328411286281372438354225D-03 /
  data alncs(  3) / +0.42934039020450834506559210803662D-06 /
  data alncs(  4) / +0.28933847795432594580466440387587D-09 /
  data alncs(  5) / +0.20512517530340580901741813447726D-12 /
  data alncs(  6) / +0.15039717055497386574615153319999D-15 /
  data alncs(  7) / +0.11294540695636464284521613333333D-18 /
  data alncs(  8) / +0.86355788671171868881946666666666D-22 /
  data alncs(  9) / +0.66952990534350370613333333333333D-25 /
  data alncs( 10) / +0.52491557448151466666666666666666D-28 /
  data alncs( 11) / +0.41530540680362666666666666666666D-31 /

  data center(1) / 1.0D+00 /
  data center(2) / 1.25D+00 /
  data center(3) / 1.50D+00 /
  data center(4) / 1.75D+00 /

  data alncen(1) / +0.0D+00 /
  data alncen(2) / +0.22314355131420975576629509030983D+00 /
  data alncen(3) / +0.40546510810816438197801311546434D+00 /
  data alncen(4) / +0.55961578793542268627088850052682D+00 /
  data alncen(5) / +0.69314718055994530941723212145817D+00 /

  data aln2 / 0.06814718055994530941723212145818D+00 /
  data nterms / 0 /

  if ( nterms == 0 ) then
    nterms = r8_inits ( alncs, 11, 28.9D+00 * r8_mach ( 3 ) )
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LOG - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.0'
    stop
  end if

  call r8_upak ( x, y, n )

  xn = real ( n - 1, kind = 8 )
  y = 2.0D+00 * y
  ntrval = int ( 4.0D+00 * y - 2.5D+00 )

  if ( ntrval == 5 ) then
    t = ( ( y - 1.0D+00 ) - 1.0D+00 ) / ( y + 2.0D+00 )
  else if ( ntrval < 5 ) then
    t = ( y - center(ntrval) ) / ( y + center(ntrval) )
  end if

  t2 = t * t
  r8_log = 0.625D+00 * xn + ( aln2 * xn + alncen(ntrval) &
    + 2.0D+00 * t + t * t2 &
    * r8_csevl ( 578.0D+00 * t2 - 1.0D+00, alncs, nterms ) )

  return
end
function r8_log10 ( x )

!*****************************************************************************80
!
!! R8_LOG10 evaluates the logarithm, base 10, of an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) R8_LOG10, the logarithm, base 10, of X.
!
  implicit none

  real ( kind = 8 ) aloge
  real ( kind = 8 ) r8_log10
  real ( kind = 8 ) x

  save aloge

  data aloge / 0.43429448190325182765112891891661D+00 /

  r8_log10 = aloge * log ( x )

  return
end
function r8_mach ( i )

!*****************************************************************************80
!
!! R8_MACH returns real ( kind = 8 ) real machine-dependent constants.
!
!  Discussion:
!
!    R8_MACH can be used to obtain machine-dependent parameters
!    for the local machine environment.  It is a function
!    with one input argument, and can be called as follows:
!
!      D = R8_MACH ( I )
!
!    where I=1,...,5.  The output value of D above is
!    determined by the input value of I:.
!
!    R8_MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
!    R8_MACH ( 2) = B^EMAX*(1 - B^(-T)), the largest magnitude.
!    R8_MACH ( 3) = B^(-T), the smallest relative spacing.
!    R8_MACH ( 4) = B^(1-T), the largest relative spacing.
!    R8_MACH ( 5) = LOG10(B)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528:
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired constant.
!
!    Output, real ( kind = 8 ) R8_MACH, the value of the constant.
!
  implicit none

  real ( kind = 8 ) r8_mach
  integer ( kind = 4 ) i

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r8_mach = 0.0D+00
    stop
  else if ( i == 1 ) then
    r8_mach = 4.450147717014403D-308
  else if ( i == 2 ) then
    r8_mach = 8.988465674311579D+307
  else if ( i == 3 ) then
    r8_mach = 1.110223024625157D-016
  else if ( i == 4 ) then
    r8_mach = 2.220446049250313D-016
  else if ( i == 5 ) then
    r8_mach = 0.301029995663981D+000
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r8_mach = 0.0D+00
    stop
  end if

  return
end
subroutine r8_machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, &
  minexp, maxexp, eps, epsneg, xmin, xmax )

!*****************************************************************************80
!
!! R8_MACHAR determines double precision machine constants.
!
!  Discussion:
!
!    This routine determines the parameters of the floating-point 
!    arithmetic system specified below.  The determination of the first 
!    three uses an extension of an algorithm due to Malcolm, 
!    incorporating some of the improvements suggested by Gentleman and 
!    Marovich.  
!
!    This routine appeared as ACM algorithm 665.
!
!    An earlier version of this program was published in Cody and Waite.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
!    machine parameters,
!    ACM Transactions on Mathematical Software,
!    Volume 14, Number 4, pages 303-311, 1988.
!
!    William Cody, William Waite,
!    Software Manual for the Elementary Functions,
!    Prentice Hall, 1980.
!
!    Morven Gentleman, Scott Marovich,
!    Communications of the ACM,
!    Volume 17, pages 276-277, 1974.
!
!    Michael Malcolm,
!    Communications of the ACM,
!    Volume 15, pages 949-951, 1972.
!
!  Parameters:
!
!    Output, integer ( kind = 8 ) IBETA, the radix for the floating-point
!    representation.
!
!    Output, integer ( kind = 8 ) IT, the number of base IBETA digits in 
!    the floating-point significand.
!
!    Output, integer ( kind = 8 ) IRND:
!    0, if floating-point addition chops.
!    1, if floating-point addition rounds, but not in the IEEE style.
!    2, if floating-point addition rounds in the IEEE style.
!    3, if floating-point addition chops, and there is partial underflow.
!    4, if floating-point addition rounds, but not in the IEEE style, and 
!      there is partial underflow.
!    5, if floating-point addition rounds in the IEEE style, and there is 
!      partial underflow.
!
!    Output, integer ( kind = 8 ) NGRD, the number of guard digits for 
!    multiplication with truncating arithmetic.  It is
!    0, if floating-point arithmetic rounds, or if it truncates and only 
!      IT base IBETA digits participate in the post-normalization shift of the
!      floating-point significand in multiplication;
!    1, if floating-point arithmetic truncates and more than IT base IBETA
!      digits participate in the post-normalization shift of the floating-point
!      significand in multiplication.
!
!    Output, integer ( kind = 8 ) MACHEP, the largest negative integer 
!    such that
!      1.0 < 1.0 + real ( IBETA, kind = 8 ) ** MACHEP, 
!    except that MACHEP is bounded below by - ( IT + 3 ).
!
!    Output, integer ( kind = 8 ) NEGEPS, the largest negative integer 
!    such that
!      1.0 - real ( IBETA, kind = 8 ) ** NEGEPS < 1.0, 
!    except that NEGEPS is bounded below by - ( IT + 3 ).
!
!    Output, integer ( kind = 8 ) IEXP, the number of bits (decimal places 
!    if IBETA = 10) reserved for the representation of the exponent (including
!    the bias or sign) of a floating-point number.
!
!    Output, integer ( kind = 8 ) MINEXP, the largest in magnitude negative
!    integer such that
!      real ( IBETA, kind = 8 ) ** MINEXP 
!    is positive and normalized.
!
!    Output, integer ( kind = 8 ) MAXEXP, the smallest positive power of
!    BETA that overflows.
!
!    Output, real ( kind = 8 ) EPS, the smallest positive floating-point number
!    such that  
!      1.0 + EPS /= 1.0. 
!    in particular, if either IBETA = 2  or IRND = 0, 
!      EPS = real ( IBETA, kind = 8 ) ** MACHEP.
!    Otherwise,  
!      EPS = ( real ( IBETA, kind = 8 ) ** MACHEP ) / 2.
!
!    Output, real ( kind = 8 ) EPSNEG, a small positive floating-point number
!    such that
!      1.0 - EPSNEG < 1.0. 
!    In particular, if IBETA = 2 or IRND = 0, 
!      EPSNEG = real ( IBETA, kind = 8 ) ** NEGEPS.
!    Otherwise,  
!      EPSNEG = ( real ( IBETA, kind = 8 ) ** NEGEPS ) / 2.  
!    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
!    smallest number that can alter 1.0 by subtraction.
!
!    Output, real ( kind = 8 ) XMIN, the smallest non-vanishing normalized
!    floating-point power of the radix:
!      XMIN = real ( IBETA, kind = 8 ) ** MINEXP
!
!    Output, real ( kind = 8 ) XMAX, the largest finite floating-point number.
!    In particular,
!      XMAX = ( 1.0 - EPSNEG ) * real ( IBETA, kind = 8 ) ** MAXEXP
!    On some machines, the computed value of XMAX will be only the second, 
!    or perhaps third, largest number, being too small by 1 or 2 units in 
!    the last digit of the significand.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) betah
  real ( kind = 8 ) betain
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ibeta
  integer ( kind = 8 ) iexp
  integer ( kind = 8 ) irnd
  integer ( kind = 8 ) it
  integer ( kind = 8 ) itemp
  integer ( kind = 8 ) iz
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  integer ( kind = 8 ) machep
  integer ( kind = 8 ) maxexp
  integer ( kind = 8 ) minexp
  integer ( kind = 8 ) mx
  integer ( kind = 8 ) negep
  integer ( kind = 8 ) ngrd
  integer ( kind = 8 ) nxres
  real ( kind = 8 ) one
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp1
  real ( kind = 8 ) tempa
  real ( kind = 8 ) two
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero

  one = real ( 1, kind = 8 )
  two = one + one
  zero = one - one
!
!  Determine IBETA and BETA ala Malcolm.
!
  a = one

  do

    a = a + a
    temp = a + one
    temp1 = temp - a

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  b = one

  do

    b = b + b
    temp = a + b
    itemp = int ( temp - a )

    if ( itemp /= 0 ) then
      exit
    end if

  end do

  ibeta = itemp
  beta = real ( ibeta, kind = 8 )
!
!  Determine IT and IRND.
!
  it = 0
  b = one

  do

    it = it + 1
    b = b * beta
    temp = b + one
    temp1 = temp - b

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  irnd = 0
  betah = beta / two
  temp = a + betah

  if ( temp - a /= zero ) then
    irnd = 1
  end if

  tempa = a + beta
  temp = tempa + betah

  if ( irnd == 0 .and. temp - tempa /= zero ) then
    irnd = 2
  end if
!
!  Determine NEGEP and EPSNEG.
!
  negep = it + 3
  betain = one / beta
  a = one
  do i = 1, negep
    a = a * betain
  end do

  b = a

  do

    temp = one - a

    if ( temp - one /= zero ) then
      exit
    end if

    a = a * beta
    negep = negep - 1

  end do

  negep = -negep
  epsneg = a

  if ( ibeta /= 2 .and. irnd /= 0 ) then

    a = ( a * ( one + a ) ) / two
    temp = one - a

    if ( temp - one /= zero ) then
      epsneg = a
    end if

  end if
!
!  Determine MACHEP and EPS.
!
  machep = -it - 3
  a = b

  do

    temp = one + a

    if ( temp - one /= zero ) then
      exit
    end if

    a = a * beta
    machep = machep + 1

  end do

  eps = a
  temp = tempa + beta * ( one + eps )

  if ( ibeta /= 2 .and. irnd /= 0 ) then

    a = ( a * ( one + a ) ) / two
    temp = one + a

    if ( temp - one /= zero ) then
      eps = a
    end if

  end if
!
!  Determine NGRD.
!
  ngrd = 0
  temp = one + eps

  if ( irnd == 0 .and. temp * one - one /= zero ) then
    ngrd = 1
  end if
!
!  Determine IEXP, MINEXP and XMIN.
!
!  Loop to determine largest I and K = 2**I such that (1/BETA) ** (2**(I))
!  does not underflow.  Exit from loop is signaled by an underflow.
!
  i = 0
  k = 1
  z = betain
  t = one + eps
  nxres = 0

  do

    y = z
    z = y * y

    a = z * one
    temp = z * t

    if ( a + a == zero .or. y <= abs ( z ) ) then
      exit
    end if

    temp1 = temp * betain

    if ( temp1 * beta == z ) then
      exit
    end if

    i = i + 1
    k = k + k

  end do
!
!  This segment is for nondecimal machines.
!
  if ( ibeta /= 10 ) then

    iexp = i + 1
    mx = k + k
!
!  This segment is for decimal machines only.
!
  else

    iexp = 2
    iz = ibeta

    do

      if ( k < iz ) then
        exit
      end if

      iz = iz * ibeta
      iexp = iexp + 1

    end do

    mx = iz + iz - 1

  end if
!
!  Loop to determine MINEXP, XMIN.
!  Exit from loop is signaled by an underflow.
!
  do

    xmin = y
    y = y * betain

    a = y * one
    temp = y * t

    if ( a + a == zero .or. xmin <= abs ( y ) ) then
      exit
    end if

    k = k + 1
    temp1 = temp * betain

    if ( temp1 * beta == y ) then
      nxres = 3
      xmin = y
      exit
    end if

  end do

  minexp = -k
!
!  Determine MAXEXP and XMAX.
!
  if ( mx <= k + k - 3 .and. ibeta /= 10 ) then
    mx = mx + mx
    iexp = iexp + 1
  end if

  maxexp = mx + minexp
!
!  Adjust IRND to reflect partial underflow.
!
  irnd = irnd + nxres
!
!  Adjust for IEEE-style machines.
!
  if ( irnd == 2 .or. irnd == 5 ) then
    maxexp = maxexp - 2
  end if
!
!  Adjust for non-IEEE machines with partial underflow.
!
  if ( irnd == 3 .or. irnd == 4 ) then
    maxexp = maxexp - it
  end if
!
!  Adjust for machines with implicit leading bit in binary significand, 
!  and machines with radix point at extreme right of significand.
!
  i = maxexp + minexp

  if ( ibeta == 2 .and. i == 0 ) then
    maxexp = maxexp - 1
  end if

  if ( 20 < i ) then
    maxexp = maxexp - 1
  end if

  if ( a /= y ) then
    maxexp = maxexp - 2
  end if

  xmax = one - epsneg

  if ( xmax * one /= xmax ) then
    xmax = one - beta * epsneg
  end if

  xmax = xmax / ( beta * beta * beta * xmin )

  i = maxexp + minexp + 3

  do j = 1, i

    if ( ibeta == 2 ) then
      xmax = xmax + xmax
    else
      xmax = xmax * beta
    end if

  end do

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8 value.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) real value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, real ( kind = 8 ) R8_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if

  return
end
function r8_pak ( y, n )

!*****************************************************************************80
!
!! R8_PAK packs a base 2 exponent into an R8.
!
!  Discussion:
!
!    This routine is almost the inverse of R8_UPAK.  It is not exactly 
!    the inverse, because abs ( x ) need not be between 0.5 and 1.0.  
!    If both R8_PAK and 2.0^n were known to be in range, we could compute
!    R8_PAK = x * 2.0^n.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Y, the mantissa.
!
!    Input, integer ( kind = 4 ) N, the exponent.
!
!    Output, real ( kind = 8 ) R8_PAK, the packed value.
!
  implicit none

  real ( kind = 8 ) aln210
  real ( kind = 8 ) aln2b
  integer ( kind = 4 ) i4_mach
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) nsum
  integer ( kind = 4 ) ny
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_pak
  real ( kind = 8 ) value
  real ( kind = 8 ) y

  save aln210
  save nmax
  save nmin

  data aln210 / 3.321928094887362347870319429489D+00 /
  data nmax / 0 /
  data nmin / 0 /

  if ( nmin == 0 ) then
    aln2b = 1.0D+00
    if ( i4_mach ( 10 ) /= 2 ) then
      aln2b = r8_mach ( 5 ) * aln210
    end if
    nmin = aln2b * real ( i4_mach ( 15 ), kind = 8 )
    nmax = aln2b * real ( i4_mach ( 16 ), kind = 8 )
  end if

  call r8_upak ( y, value, ny )

  nsum = n + ny

  if ( nsum < nmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_PAK - Warning!'
    write ( *, '(a)' ) '  Packed number underflows.'
    r8_pak = 0.0D+00
    return
  end if

  if ( nmax < nsum ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_PAK - Fatal error!'
    write ( *, '(a)' ) '  Ppacked number overflows.'
    stop
  end if

  do while ( nsum < 0 )
    value = 0.5D+00 * value
    nsum = nsum + 1
  end do

  do while ( 0 < nsum )
    value = 2.0D+00 * value
    nsum = nsum - 1
  end do

  r8_pak = value

  return
end
function r8_poch ( a, x )

!*****************************************************************************80
!
!! R8_POCH evaluates Pochhammer's function of R8 arguments.
!
!  Discussion:
!
!    POCH ( A, X ) = Gamma ( A + X ) / Gamma ( A ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, X, the arguments.
!
!    Output, real ( kind = 8 ) R8_POCH, the Pochhammer function of A and X.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absa
  real ( kind = 8 ) absax
  real ( kind = 8 ) alnga
  real ( kind = 8 ) alngax
  real ( kind = 8 ) ax
  real ( kind = 8 ) b
  real ( kind = 8 ) cospia
  real ( kind = 8 ) cospix
  real ( kind = 8 ) den
  real ( kind = 8 ) eps
  real ( kind = 8 ) err
  real ( kind = 8 ) errpch
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) n
  real ( kind = 8 ) pi
  real ( kind = 8 ) r8_fac
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_gamr
  real ( kind = 8 ) r8_lgmc
  real ( kind = 8 ) r8_lnrel
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) r8_poch
  real ( kind = 8 ) sgnga
  real ( kind = 8 ) sgngax
  real ( kind = 8 ) sinpia
  real ( kind = 8 ) sinpix
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) x

  save eps
  save pi

  data eps / 0.0D+00 /
  data pi / 3.141592653589793238462643383279503D+00 /

  if ( eps == 0.0D+00 ) then
    eps = r8_mach ( 4 )
    sqeps = sqrt ( eps )
  end if

  ax = a + x

  if ( ax <= 0.0D+00 .and. aint ( ax ) == ax ) then

    if ( 0.0D+00 < a .or. int ( a ) /= a ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_POCH - Fatal error!'
      write ( *, '(a)' ) '  A + X is nonpositive integer,'
      write ( *, '(a)' ) '  but A is not.'
      stop
    end if
!
!  We know here that both A+X and A are non-positive integers.
!
    if ( x == 0.0D+00 ) then
      r8_poch = 1.0D+00
    else if ( - 20.0D+00 < min ( a + x, a ) ) then
      n = int ( x )
      ia = int ( a )
      r8_poch = r8_mop ( n ) * r8_fac ( - ia ) / r8_fac ( - ia - n )
    else
      n = int ( x )
      r8_poch = r8_mop ( n ) * exp ( ( a - 0.5D+00 ) &
        * r8_lnrel ( x / ( a - 1.0D+00 ) ) &
        + x * log ( - a + 1.0D+00 - x ) - x &
        + r8_lgmc ( - a + 1.0D+00 ) &
        - r8_lgmc ( - a - x + 1.0D+00 ) )
    end if

    return

  end if
!
!  A + X is not zero or a negative integer.
!
  if ( a <= 0.0D+00 .and. aint ( a ) == a ) then
    r8_poch = 0.0D+00
    return
  end if

  n = abs ( x )
!
!  X is a small non-positive integer, presummably a common case.
!
  if ( real ( n, kind = 8 ) == x .and. n <= 20 ) then
    r8_poch = 1.0D+00
    do i = 1, n
      r8_poch = r8_poch * ( a + real ( i - 1, kind = 8 ) )
    end do
    return
  end if

  absax = abs ( a + x )
  absa = abs ( a )

  if ( max ( absax, absa ) <= 20.0D+00 ) then
    r8_poch = r8_gamma ( a + x ) * r8_gamr ( a )
    return
  end if

  if ( 0.5D+00 * absa < abs ( x ) ) then
    call r8_lgams ( a + x, alngax, sgngax )
    call r8_lgams ( a, alnga, sgnga )
    r8_poch = sgngax * sgnga * exp ( alngax - alnga )
    return
  end if
!
!  abs(x) is small and both abs(a+x) and abs(a) are large.  thus,
!  a+x and a must have the same sign.  for negative a, we use
!  gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) *
!  sin(pi*a)/sin(pi*(a+x))
!
  if ( a < 0.0D+00 ) then
    b = - a - x + 1.0D+00
  else
    b = a
  end if

  r8_poch = exp ( ( b - 0.5D+00 ) * r8_lnrel ( x / b ) &
    + x * log ( b + x ) - x + r8_lgmc ( b + x ) - r8_lgmc ( b ) )

  if ( 0.0D+00 <= a .or. r8_poch == 0.0D+00 ) then
    return
  end if

  cospix = cos ( pi * x )
  sinpix = sin ( pi * x )
  cospia = cos ( pi * a )
  sinpia = sin ( pi * a )

  errpch = abs ( x ) * ( 1.0D+00 + log ( b ) )
  den = cospix + cospia * sinpix / sinpia
  err = ( abs ( x ) * ( abs ( sinpix ) &
    + abs ( cospia * cospix / sinpia ) ) &
    + abs ( a * sinpix ) / sinpia / sinpia ) * pi
  err = errpch + err / abs ( den )

  r8_poch = r8_poch / den

  return
end
function r8_poch1 ( a, x )

!*****************************************************************************80
!
!! R8_POCH1 evaluates a quantity related to Pochhammer's symbol.
!
!  Discussion:
!
!    Evaluate a generalization of Pochhammer's symbol for special
!    situations that require especially accurate values when x is small in
!      poch1(a,x) = (poch(a,x)-1)/x
!                 = (gamma(a+x)/gamma(a) - 1.0)/x .
!    This specification is particularly suited for stably computing
!    expressions such as
!      (gamma(a+x)/gamma(a) - gamma(b+x)/gamma(b))/x
!           = poch1(a,x) - poch1(b,x)
!    Note that poch1(a,0.0) = psi(a)
!
!    When abs ( x ) is so small that substantial cancellation will occur if
!    the straightforward formula is used, we  use an expansion due
!    to fields and discussed by y. l. luke, the special functions and their
!    approximations, vol. 1, academic press, 1969, page 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) R8_POCH1, the value of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absa
  real ( kind = 8 ) absx
  real ( kind = 8 ) alneps
  real ( kind = 8 ) alnvar
  real ( kind = 8 ) b
  real ( kind = 8 ) bern(20)
  real ( kind = 8 ) binv
  real ( kind = 8 ) bp
  real ( kind = 8 ) gbern(21)
  real ( kind = 8 ) gbk
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ndx
  integer ( kind = 4 ) nterms
  real ( kind = 8 ) pi
  real ( kind = 8 ) poly1
  real ( kind = 8 ) q
  real ( kind = 8 ) r8_cot
  real ( kind = 8 ) r8_exprel
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_poch
  real ( kind = 8 ) r8_poch1
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) rho
  real ( kind = 8 ) sinpxx
  real ( kind = 8 ) sinpx2
  real ( kind = 8 ) sqtbig
  real ( kind = 8 ) term
  real ( kind = 8 ) trig
  real ( kind = 8 ) var
  real ( kind = 8 ) var2
  real ( kind = 8 ) x

  save alneps
  save bern
  save pi
  save sqtbig

  data bern(  1) / +0.833333333333333333333333333333333D-01    /
  data bern(  2) / -0.138888888888888888888888888888888D-02    /
  data bern(  3) / +0.330687830687830687830687830687830D-04    /
  data bern(  4) / -0.826719576719576719576719576719576D-06    /
  data bern(  5) / +0.208767569878680989792100903212014D-07    /
  data bern(  6) / -0.528419013868749318484768220217955D-09    /
  data bern(  7) / +0.133825365306846788328269809751291D-10   /
  data bern(  8) / -0.338968029632258286683019539124944D-12   /
  data bern(  9) / +0.858606205627784456413590545042562D-14   /
  data bern( 10) / -0.217486869855806187304151642386591D-15   /
  data bern( 11) / +0.550900282836022951520265260890225D-17   /
  data bern( 12) / -0.139544646858125233407076862640635D-18   /
  data bern( 13) / +0.353470703962946747169322997780379D-20   /
  data bern( 14) / -0.895351742703754685040261131811274D-22   /
  data bern( 15) / +0.226795245233768306031095073886816D-23   /
  data bern( 16) / -0.574472439520264523834847971943400D-24   /
  data bern( 17) / +0.145517247561486490186626486727132D-26   /
  data bern( 18) / -0.368599494066531017818178247990866D-28   /
  data bern( 19) / +0.933673425709504467203255515278562D-30   /
  data bern( 20) / -0.236502241570062993455963519636983D-31   /

  data alneps / 0.0D+00 /
  data pi / 3.141592653589793238462643383279503D+00 /
  data sqtbig / 0.0D+00 /

  if ( sqtbig == 0.0D+00 ) then
    sqtbig = 1.0D+00 / sqrt ( 24.0D+00 * r8_mach ( 1 ) )
    alneps = log ( r8_mach ( 3 ) )
  end if

  if ( x == 0.0D+00 ) then
    r8_poch1 = r8_psi ( a )
    return
  end if

  absx = abs ( x )
  absa = abs ( a )

  if ( 0.1D+00 * absa < absx .or. &
    0.1D+00 < absx * log ( max ( absa, 2.0D+00 ) ) ) then
    r8_poch1 = r8_poch ( a, x )
    r8_poch1 = ( r8_poch1 - 1.0D+00 ) / x
    return
  end if

  if ( a < - 0.5D+00 ) then
    bp = 1.0D+00 - a - x
  else
    bp = a
  end if

  if ( bp < 10.0D+00 ) then
    incr = 11.0D+00 - bp
  else
    incr = 0
  end if

  b = bp + real ( incr, kind = 8 )

  var = b + 0.5D+00 * ( x - 1.0D+00 )
  alnvar = log ( var )
  q = x * alnvar

  poly1 = 0.0D+00

  if ( var < sqtbig ) then

    var2 = 1.0D+00 / var / var

    rho = 0.5D+00 * ( x + 1.0D+00 )
    gbern(1) = 1.0D+00
    gbern(2) = - rho / 12.0D+00
    term = var2
    poly1 = gbern(2) * term

    nterms = int ( - 0.5D+00 * alneps / alnvar + 1.0D+00 )

    if ( 20 < nterms ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_POCH1 - Fatal error!'
      write ( *, '(a)' ) ' 20 < NTERMS.'
      stop
    end if 

    do k = 2, nterms
      gbk = 0.0D+00
      do j = 1, k
        ndx = k - j + 1
        gbk = gbk + bern(ndx) * gbern(j)
      end do
      gbern(k+1) = - rho * gbk / real ( k, kind = 8 )
      term = term * ( real ( 2 * k - 2, kind = 8 ) - x ) &
        * ( real ( 2 * k - 1, kind = 8 ) - x ) * var2
      poly1 = poly1 + gbern(k+1) * term
    end do

  end if

  poly1 = ( x - 1.0D+00 ) * poly1
  r8_poch1 = r8_exprel ( q ) * ( alnvar + q * poly1 ) + poly1
!
!  we have r8_poch1(b,x), but bp is small, so we use backwards recursion
!  to obtain r8_poch1(bp,x).
!
  do ii = 1, incr
    i = incr - ii
    binv = 1.0D+00 / ( bp + real ( i, kind = 8 ) )
    r8_poch1 = ( r8_poch1 - binv ) / ( 1.0D+00 + x * binv )
  end do

  if ( bp == a ) then
    return
  end if
!
!  we have r8_poch1(bp,x), but a is lt -0.5.  we therefore use a reflection
!  formula to obtain r8_poch1(a,x).
!
  sinpxx = sin ( pi * x ) / x
  sinpx2 = sin ( 0.5D+00 * pi * x )
  trig = sinpxx * r8_cot ( pi * b ) - 2.0D+00 * sinpx2 * ( sinpx2 / x )

  r8_poch1 = trig + ( 1.0D+00 + x * trig ) * r8_poch1

  return
end
function r8_pow ( x, y )

!*****************************************************************************80
!
!! R8_POW computes a power of an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    1 September 2011
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the base.
!
!    Input, real ( kind = 8 ) Y, the power.
!
!    Output, real ( kind = 8 ) R8_POW, the value of X^Y.
!
  implicit none

  real ( kind = 8 ) r8_pow
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  value = x ** y

  r8_pow = value

  return
end
function r8_psi ( x )

!*****************************************************************************80
!
!! R8_PSI evaluates the psi function of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_PSI, the psi function of X.
!
  implicit none

  real ( kind = 8 ) apsics(16)
  real ( kind = 8 ) aux
  real ( kind = 8 ) dxrel
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntapsi
  integer ( kind = 4 ) ntpsi
  real ( kind = 8 ) pi
  real ( kind = 8 ) psics(42)
  real ( kind = 8 ) r8_cot
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) y

  save apsics
  save dxrel
  save ntapsi
  save ntpsi
  save pi
  save psics
  save xbig

  data psics(  1) / -0.38057080835217921520437677667039D-01 /
  data psics(  2) / +0.49141539302938712748204699654277D+00 /
  data psics(  3) / -0.56815747821244730242892064734081D-01 /
  data psics(  4) / +0.83578212259143131362775650747862D-02 /
  data psics(  5) / -0.13332328579943425998079274172393D-02 /
  data psics(  6) / +0.22031328706930824892872397979521D-03 /
  data psics(  7) / -0.37040238178456883592889086949229D-04 /
  data psics(  8) / +0.62837936548549898933651418717690D-05 /
  data psics(  9) / -0.10712639085061849855283541747074D-05 /
  data psics( 10) / +0.18312839465484165805731589810378D-06 /
  data psics( 11) / -0.31353509361808509869005779796885D-07 /
  data psics( 12) / +0.53728087762007766260471919143615D-08 /
  data psics( 13) / -0.92116814159784275717880632624730D-09 /
  data psics( 14) / +0.15798126521481822782252884032823D-09 /
  data psics( 15) / -0.27098646132380443065440589409707D-10 /
  data psics( 16) / +0.46487228599096834872947319529549D-11 /
  data psics( 17) / -0.79752725638303689726504797772737D-12 /
  data psics( 18) / +0.13682723857476992249251053892838D-12 /
  data psics( 19) / -0.23475156060658972717320677980719D-13 /
  data psics( 20) / +0.40276307155603541107907925006281D-14 /
  data psics( 21) / -0.69102518531179037846547422974771D-15 /
  data psics( 22) / +0.11856047138863349552929139525768D-15 /
  data psics( 23) / -0.20341689616261559308154210484223D-16 /
  data psics( 24) / +0.34900749686463043850374232932351D-17 /
  data psics( 25) / -0.59880146934976711003011081393493D-18 /
  data psics( 26) / +0.10273801628080588258398005712213D-18 /
  data psics( 27) / -0.17627049424561071368359260105386D-19 /
  data psics( 28) / +0.30243228018156920457454035490133D-20 /
  data psics( 29) / -0.51889168302092313774286088874666D-21 /
  data psics( 30) / +0.89027730345845713905005887487999D-22 /
  data psics( 31) / -0.15274742899426728392894971904000D-22 /
  data psics( 32) / +0.26207314798962083136358318079999D-23 /
  data psics( 33) / -0.44964642738220696772598388053333D-24 /
  data psics( 34) / +0.77147129596345107028919364266666D-25 /
  data psics( 35) / -0.13236354761887702968102638933333D-25 /
  data psics( 36) / +0.22709994362408300091277311999999D-26 /
  data psics( 37) / -0.38964190215374115954491391999999D-27 /
  data psics( 38) / +0.66851981388855302310679893333333D-28 /
  data psics( 39) / -0.11469986654920864872529919999999D-28 /
  data psics( 40) / +0.19679385886541405920515413333333D-29 /
  data psics( 41) / -0.33764488189750979801907200000000D-30 /
  data psics( 42) / +0.57930703193214159246677333333333D-31 /

  data apsics(  1) / -0.832710791069290760174456932269D-03 /
  data apsics(  2) / -0.416251842192739352821627121990D-03 /
  data apsics(  3) / +0.103431560978741291174463193961D-06 /
  data apsics(  4) / -0.121468184135904152987299556365D-09 /
  data apsics(  5) / +0.311369431998356155521240278178D-12 /
  data apsics(  6) / -0.136461337193177041776516100945D-14 /
  data apsics(  7) / +0.902051751315416565130837974000D-17 /
  data apsics(  8) / -0.831542997421591464829933635466D-19 /
  data apsics(  9) / +0.101224257073907254188479482666D-20 /
  data apsics( 10) / -0.156270249435622507620478933333D-22 /
  data apsics( 11) / +0.296542716808903896133226666666D-24 /
  data apsics( 12) / -0.674686886765702163741866666666D-26 /
  data apsics( 13) / +0.180345311697189904213333333333D-27 /
  data apsics( 14) / -0.556901618245983607466666666666D-29 /
  data apsics( 15) / +0.195867922607736251733333333333D-30 /
  data apsics( 16) / -0.775195892523335680000000000000D-32 /

  data dxrel / 0.0D+00 /
  data ntapsi / 0 /
  data ntpsi / 0 /
  data pi / 3.14159265358979323846264338327950D+00 /
  data xbig / 0.0D+00 /

  if ( ntpsi == 0 ) then
    ntpsi = r8_inits ( psics, 42, 0.1D+00 * r8_mach ( 3 ) )
    ntapsi = r8_inits ( apsics, 16, 0.1D+00 * r8_mach ( 3 ) )
    xbig = 1.0D+00 / sqrt ( r8_mach ( 3 ) )
    dxrel = sqrt ( r8_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( y < 10.0D+00 ) then

    n = int ( x )
    if ( x < 0.0D+00 ) then
      n = n - 1
    end if
    y = x - real ( n, kind = 8 )
    n = n - 1
    r8_psi = r8_csevl ( 2.0D+00 * y - 1.0D+00, psics, ntpsi )

    if ( n == 0 ) then

      return

    else if ( n < 0 ) then

      n = - n

      if ( x == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_PSI - Fatal error!'
        write ( *, '(a)' ) '  X is zero.'
        stop
      end if

      if ( x < 0.0D+00 .and. &
        x + real ( n - 2, kind = 8 ) == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_PSI - Fatal error!'
        write ( *, '(a)' ) '  X is a negative integer.'
        stop
      end if

      if ( x < - 0.5D+00 .and. &
        abs ( ( x - aint ( x - 0.5D+00 ) ) / x ) < dxrel ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_PSI - Warning!'
        write ( *, '(a)' ) '  Answer is less than half precision'
        write ( *, '(a)' ) '  because X is near a negative integer.'
      end if

      do i = 1, n
        r8_psi = r8_psi - 1.0D+00 / ( x + real ( i - 1, kind = 8 ) )
      end do

    else if ( 0 < n ) then

      do i = 1, n
        r8_psi = r8_psi + 1.0D+00 / ( y + real ( i, kind = 8 ) )
      end do

    end if

  else

    if ( y < xbig ) then
      aux = r8_csevl ( 8.0D+00 / y / y - 1.0D+00, apsics, ntapsi )
    else
      aux = 0.0D+00
    end if

    if ( x < 0.0D+00 ) then

      r8_psi = log ( abs ( x ) ) - 0.5D+00 / x + aux &
        - pi * r8_cot ( pi * x )

    else if ( 0.0D+00 < x ) then

      r8_psi = log ( x ) - 0.5D+00 / x + aux

    end if

  end if

  return
end
function r8_ren ( )

!*****************************************************************************80
!
!! R8_REN is a simple random number generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Malcolm Pike, David Hill,
!    Algorithm 266:
!    Pseudo-Random Numbers,
!    Communications of the ACM,
!    Volume 8, Number 10, October 1965, page 605.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_REN, the random value.
!
  implicit none

  integer ( kind = 4 ) iy
  real ( kind = 8 ) r8_ren

  save iy

  data iy / 100001 /

  iy = iy * 125
  iy = iy - ( iy / 2796203 ) * 2796203
  r8_ren = real ( iy, kind = 8 ) / 2796203.0D+00

  return
end
function r8_shi ( x )

!*****************************************************************************80
!
!! R8_SHI evaluates the hyperbolic sine integral Shi of an R8 argument.
!
!  Discussion:
!
!    Shi ( x ) = Integral ( 0 <= t <= x ) sinh ( t ) dt / t
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_SHI, the hyperbolic sine integral 
!    Shi evaluated at X.
!
  implicit none

  real ( kind = 8 ) absx
  integer ( kind = 4 ) nshi
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_e1
  real ( kind = 8 ) r8_ei
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_shi
  real ( kind = 8 ) shics(10)
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml

  save nshi
  save shics
  save xsml

  data shics(  1) /     0.0078372685688900950695200984317332D+00 /
  data shics(  2) /     0.0039227664934234563972697574427225D+00 /
  data shics(  3) /     0.0000041346787887617266746747908275D+00 /
  data shics(  4) /     0.0000000024707480372882742135145302D+00 /
  data shics(  5) /     0.0000000000009379295590763630457157D+00 /
  data shics(  6) /     0.0000000000000002451817019520867353D+00 /
  data shics(  7) /     0.0000000000000000000467416155257592D+00 /
  data shics(  8) /     0.0000000000000000000000067803072389D+00 /
  data shics(  9) /     0.0000000000000000000000000007731289D+00 /
  data shics( 10) /     0.0000000000000000000000000000000711D+00 /

  data nshi / 0 /
  data xsml / 0.0D+00 /

  if ( nshi == 0 ) then
    nshi = r8_inits ( shics, 10, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( r8_mach ( 3 ) )
  end if

  absx = abs ( x )

  if ( absx <= xsml ) then
    r8_shi = x
  else if ( absx <= 0.375D+00 ) then
    r8_shi = x * ( 1.0D+00 + r8_csevl ( 128.0D+00 * x * x / 9.0D+00 - 1.0D+00, &
      shics, nshi ) )
  else
    r8_shi = 0.5D+00 * ( r8_ei ( x ) + r8_e1 ( x ) )
  end if

  return
end
function r8_si ( x )

!*****************************************************************************80
!
!! R8_SI evaluates the sine integral Si of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_SI, the sine integral Si evaluated at X.
!
  implicit none

  real ( kind = 8 ) absx
  real ( kind = 8 ) cosx
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) nsi
  real ( kind = 8 ) pi2
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_si
  real ( kind = 8 ) sics(18)
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml

  save nsi
  save pi2
  save sics
  save xsml

  data sics(  1) / -0.1315646598184841928904275173000457D+00 /
  data sics(  2) / -0.2776578526973601892048287660157299D+00 /
  data sics(  3) /  0.0354414054866659179749135464710086D+00 /
  data sics(  4) / -0.0025631631447933977658752788361530D+00 /
  data sics(  5) /  0.0001162365390497009281264921482985D+00 /
  data sics(  6) / -0.0000035904327241606042670004347148D+00 /
  data sics(  7) /  0.0000000802342123705710162308652976D+00 /
  data sics(  8) / -0.0000000013562997692540250649931846D+00 /
  data sics(  9) /  0.0000000000179440721599736775567759D+00 /
  data sics( 10) / -0.0000000000001908387343087145490737D+00 /
  data sics( 11) /  0.0000000000000016669989586824330853D+00 /
  data sics( 12) / -0.0000000000000000121730988368503042D+00 /
  data sics( 13) /  0.0000000000000000000754181866993865D+00 /
  data sics( 14) / -0.0000000000000000000004014178842446D+00 /
  data sics( 15) /  0.0000000000000000000000018553690716D+00 /
  data sics( 16) / -0.0000000000000000000000000075166966D+00 /
  data sics( 17) /  0.0000000000000000000000000000269113D+00 /
  data sics( 18) / -0.0000000000000000000000000000000858D+00 /

  data nsi / 0 /
  data pi2 / 1.57079632679489661923132169163975D+00 /
  data xsml / 0.0D+00 /

  if ( nsi == 0 ) then
    nsi = r8_inits ( sics, 18, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( r8_mach ( 3 ) )
  end if

  absx = abs ( x )

  if ( absx < xsml ) then
    r8_si = x
  else if ( absx <= 4.0D+00 ) then
    r8_si = x * ( 0.75D+00 &
      + r8_csevl ( ( x * x - 8.0D+00 ) * 0.125D+00, sics, nsi ) )
  else
    call r8_sifg ( absx, f, g )
    cosx = cos ( absx )
    r8_si = pi2 - f * cosx - g * sin ( x )
    if ( x < 0.0D+00 ) then
      r8_si = - r8_si
    end if
  end if

  return
end
subroutine r8_sifg ( x, f, g )

!*****************************************************************************80
!
!! R8_SIFG is a utility routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) F, G.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) f1cs(43)
  real ( kind = 8 ) f2cs(99)
  real ( kind = 8 ) g
  real ( kind = 8 ) g1cs(44)
  real ( kind = 8 ) g2cs(44)
  real ( kind = 8 ) g3cs(56)
  integer ( kind = 4 ) nf1
  integer ( kind = 4 ) nf2
  integer ( kind = 4 ) ng1
  integer ( kind = 4 ) ng2
  integer ( kind = 4 ) ng3
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) tol
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xbnd
  real ( kind = 8 ) xbndg
  real ( kind = 8 ) xmaxf
  real ( kind = 8 ) xmaxg

  save f1cs
  save f2cs
  save g1cs
  save g2cs
  save g3cs
  save nf1
  save nf2
  save ng1
  save ng2
  save ng3
  save xbig
  save xbnd
  save xbndg
  save xmaxf
  save xmaxg

  data f1cs(  1) / -0.1191081969051363610348201965828918D+00 /
  data f1cs(  2) / -0.0247823144996236247590074150823133D+00 /
  data f1cs(  3) /  0.0011910281453357821268120363054457D+00 /
  data f1cs(  4) / -0.0000927027714388561748308600360706D+00 /
  data f1cs(  5) /  0.0000093373141568270996868204582766D+00 /
  data f1cs(  6) / -0.0000011058287820557143938979426306D+00 /
  data f1cs(  7) /  0.0000001464772071460162169336550799D+00 /
  data f1cs(  8) / -0.0000000210694496287689532601227548D+00 /
  data f1cs(  9) /  0.0000000032293492366848236382857374D+00 /
  data f1cs( 10) / -0.0000000005206529617529375828014986D+00 /
  data f1cs( 11) /  0.0000000000874878884570278750268316D+00 /
  data f1cs( 12) / -0.0000000000152176187056123668294574D+00 /
  data f1cs( 13) /  0.0000000000027257192405419573900583D+00 /
  data f1cs( 14) / -0.0000000000005007053075968556290255D+00 /
  data f1cs( 15) /  0.0000000000000940240902726068511779D+00 /
  data f1cs( 16) / -0.0000000000000180014444791803678336D+00 /
  data f1cs( 17) /  0.0000000000000035062621432741785826D+00 /
  data f1cs( 18) / -0.0000000000000006935282926769149709D+00 /
  data f1cs( 19) /  0.0000000000000001390925136454216568D+00 /
  data f1cs( 20) / -0.0000000000000000282486885074170585D+00 /
  data f1cs( 21) /  0.0000000000000000058031305693579081D+00 /
  data f1cs( 22) / -0.0000000000000000012046901573375820D+00 /
  data f1cs( 23) /  0.0000000000000000002525052443655940D+00 /
  data f1cs( 24) / -0.0000000000000000000533980268805594D+00 /
  data f1cs( 25) /  0.0000000000000000000113855786274122D+00 /
  data f1cs( 26) / -0.0000000000000000000024462861505259D+00 /
  data f1cs( 27) /  0.0000000000000000000005293659320439D+00 /
  data f1cs( 28) / -0.0000000000000000000001153184940277D+00 /
  data f1cs( 29) /  0.0000000000000000000000252786568318D+00 /
  data f1cs( 30) / -0.0000000000000000000000055738645378D+00 /
  data f1cs( 31) /  0.0000000000000000000000012358245621D+00 /
  data f1cs( 32) / -0.0000000000000000000000002754350842D+00 /
  data f1cs( 33) /  0.0000000000000000000000000616906808D+00 /
  data f1cs( 34) / -0.0000000000000000000000000138817443D+00 /
  data f1cs( 35) /  0.0000000000000000000000000031375329D+00 /
  data f1cs( 36) / -0.0000000000000000000000000007121249D+00 /
  data f1cs( 37) /  0.0000000000000000000000000001622778D+00 /
  data f1cs( 38) / -0.0000000000000000000000000000371206D+00 /
  data f1cs( 39) /  0.0000000000000000000000000000085221D+00 /
  data f1cs( 40) / -0.0000000000000000000000000000019633D+00 /
  data f1cs( 41) /  0.0000000000000000000000000000004538D+00 /
  data f1cs( 42) / -0.0000000000000000000000000000001052D+00 /
  data f1cs( 43) /  0.0000000000000000000000000000000245D+00 /

  data f2cs(  1) / -0.03484092538970132330836049733745577D+00 /
  data f2cs(  2) / -0.01668422056779596873246786312278676D+00 /
  data f2cs(  3) /  0.00067529012412377385045207859239727D+00 /
  data f2cs(  4) / -0.00005350666225447013628785577557429D+00 /
  data f2cs(  5) /  0.00000626934217790075267050759431626D+00 /
  data f2cs(  6) / -0.00000095266388019916680677790414293D+00 /
  data f2cs(  7) /  0.00000017456292242509880425504427666D+00 /
  data f2cs(  8) / -0.00000003687954030653093307097646628D+00 /
  data f2cs(  9) /  0.00000000872026777051395264075816938D+00 /
  data f2cs( 10) / -0.00000000226019703919738748530423167D+00 /
  data f2cs( 11) /  0.00000000063246249765250612520444877D+00 /
  data f2cs( 12) / -0.00000000018889118884717869240911480D+00 /
  data f2cs( 13) /  0.00000000005967746729997813372620472D+00 /
  data f2cs( 14) / -0.00000000001980443117372239011196007D+00 /
  data f2cs( 15) /  0.00000000000686413954772103383713264D+00 /
  data f2cs( 16) / -0.00000000000247310193070199106074890D+00 /
  data f2cs( 17) /  0.00000000000092263594549941404196042D+00 /
  data f2cs( 18) / -0.00000000000035523634999261784497297D+00 /
  data f2cs( 19) /  0.00000000000014076049625351591461820D+00 /
  data f2cs( 20) / -0.00000000000005726228499747652794311D+00 /
  data f2cs( 21) /  0.00000000000002386537545413171810106D+00 /
  data f2cs( 22) / -0.00000000000001017141890764597142232D+00 /
  data f2cs( 23) /  0.00000000000000442594531078364424968D+00 /
  data f2cs( 24) / -0.00000000000000196344933049189761979D+00 /
  data f2cs( 25) /  0.00000000000000088688748314810461024D+00 /
  data f2cs( 26) / -0.00000000000000040743345027311546948D+00 /
  data f2cs( 27) /  0.00000000000000019016837215675339859D+00 /
  data f2cs( 28) / -0.00000000000000009009707297478042442D+00 /
  data f2cs( 29) /  0.00000000000000004329211274095668667D+00 /
  data f2cs( 30) / -0.00000000000000002108144465322479526D+00 /
  data f2cs( 31) /  0.00000000000000001039637907026452274D+00 /
  data f2cs( 32) / -0.00000000000000000518891007948931936D+00 /
  data f2cs( 33) /  0.00000000000000000261955324869899371D+00 /
  data f2cs( 34) / -0.00000000000000000133690399951301570D+00 /
  data f2cs( 35) /  0.00000000000000000068941057702931664D+00 /
  data f2cs( 36) / -0.00000000000000000035905362610437250D+00 /
  data f2cs( 37) /  0.00000000000000000018878077255791706D+00 /
  data f2cs( 38) / -0.00000000000000000010016125265594380D+00 /
  data f2cs( 39) /  0.00000000000000000005360725691578228D+00 /
  data f2cs( 40) / -0.00000000000000000002893198974944827D+00 /
  data f2cs( 41) /  0.00000000000000000001574065100202625D+00 /
  data f2cs( 42) / -0.00000000000000000000863027106431206D+00 /
  data f2cs( 43) /  0.00000000000000000000476715602862288D+00 /
  data f2cs( 44) / -0.00000000000000000000265222739998504D+00 /
  data f2cs( 45) /  0.00000000000000000000148582865063866D+00 /
  data f2cs( 46) / -0.00000000000000000000083797235923135D+00 /
  data f2cs( 47) /  0.00000000000000000000047565916422711D+00 /
  data f2cs( 48) / -0.00000000000000000000027169073353112D+00 /
  data f2cs( 49) /  0.00000000000000000000015612738881686D+00 /
  data f2cs( 50) / -0.00000000000000000000009024555078347D+00 /
  data f2cs( 51) /  0.00000000000000000000005246097049119D+00 /
  data f2cs( 52) / -0.00000000000000000000003066450818697D+00 /
  data f2cs( 53) /  0.00000000000000000000001801996250957D+00 /
  data f2cs( 54) / -0.00000000000000000000001064443050752D+00 /
  data f2cs( 55) /  0.00000000000000000000000631942158881D+00 /
  data f2cs( 56) / -0.00000000000000000000000377013812246D+00 /
  data f2cs( 57) /  0.00000000000000000000000225997542918D+00 /
  data f2cs( 58) / -0.00000000000000000000000136100844814D+00 /
  data f2cs( 59) /  0.00000000000000000000000082333232003D+00 /
  data f2cs( 60) / -0.00000000000000000000000050025986091D+00 /
  data f2cs( 61) /  0.00000000000000000000000030526245684D+00 /
  data f2cs( 62) / -0.00000000000000000000000018705164021D+00 /
  data f2cs( 63) /  0.00000000000000000000000011508404393D+00 /
  data f2cs( 64) / -0.00000000000000000000000007108714611D+00 /
  data f2cs( 65) /  0.00000000000000000000000004408065533D+00 /
  data f2cs( 66) / -0.00000000000000000000000002743760867D+00 /
  data f2cs( 67) /  0.00000000000000000000000001714144851D+00 /
  data f2cs( 68) / -0.00000000000000000000000001074768860D+00 /
  data f2cs( 69) /  0.00000000000000000000000000676259777D+00 /
  data f2cs( 70) / -0.00000000000000000000000000426981348D+00 /
  data f2cs( 71) /  0.00000000000000000000000000270500637D+00 /
  data f2cs( 72) / -0.00000000000000000000000000171933331D+00 /
  data f2cs( 73) /  0.00000000000000000000000000109636138D+00 /
  data f2cs( 74) / -0.00000000000000000000000000070132573D+00 /
  data f2cs( 75) /  0.00000000000000000000000000045001784D+00 /
  data f2cs( 76) / -0.00000000000000000000000000028963835D+00 /
  data f2cs( 77) /  0.00000000000000000000000000018697009D+00 /
  data f2cs( 78) / -0.00000000000000000000000000012104646D+00 /
  data f2cs( 79) /  0.00000000000000000000000000007859065D+00 /
  data f2cs( 80) / -0.00000000000000000000000000005116867D+00 /
  data f2cs( 81) /  0.00000000000000000000000000003340627D+00 /
  data f2cs( 82) / -0.00000000000000000000000000002186851D+00 /
  data f2cs( 83) /  0.00000000000000000000000000001435340D+00 /
  data f2cs( 84) / -0.00000000000000000000000000000944523D+00 /
  data f2cs( 85) /  0.00000000000000000000000000000623117D+00 /
  data f2cs( 86) / -0.00000000000000000000000000000412101D+00 /
  data f2cs( 87) /  0.00000000000000000000000000000273208D+00 /
  data f2cs( 88) / -0.00000000000000000000000000000181558D+00 /
  data f2cs( 89) /  0.00000000000000000000000000000120934D+00 /
  data f2cs( 90) / -0.00000000000000000000000000000080737D+00 /
  data f2cs( 91) /  0.00000000000000000000000000000054022D+00 /
  data f2cs( 92) / -0.00000000000000000000000000000036227D+00 /
  data f2cs( 93) /  0.00000000000000000000000000000024348D+00 /
  data f2cs( 94) / -0.00000000000000000000000000000016401D+00 /
  data f2cs( 95) /  0.00000000000000000000000000000011074D+00 /
  data f2cs( 96) / -0.00000000000000000000000000000007497D+00 /
  data f2cs( 97) /  0.00000000000000000000000000000005091D+00 /
  data f2cs( 98) / -0.00000000000000000000000000000003470D+00 /
  data f2cs( 99) /  0.00000000000000000000000000000002377D+00 /

  data g1cs(  1) / -0.3040578798253495954499726682091083D+00 /
  data g1cs(  2) / -0.0566890984597120587731339156118269D+00 /
  data g1cs(  3) /  0.0039046158173275643919984071554082D+00 /
  data g1cs(  4) / -0.0003746075959202260618619339867489D+00 /
  data g1cs(  5) /  0.0000435431556559843679552220840065D+00 /
  data g1cs(  6) / -0.0000057417294453025046561970723475D+00 /
  data g1cs(  7) /  0.0000008282552104502629741937616492D+00 /
  data g1cs(  8) / -0.0000001278245892594642727883913223D+00 /
  data g1cs(  9) /  0.0000000207978352948687884439257529D+00 /
  data g1cs( 10) / -0.0000000035313205921990798042032682D+00 /
  data g1cs( 11) /  0.0000000006210824236308951068631449D+00 /
  data g1cs( 12) / -0.0000000001125215474446292649336987D+00 /
  data g1cs( 13) /  0.0000000000209088917684421605267019D+00 /
  data g1cs( 14) / -0.0000000000039715831737681727689158D+00 /
  data g1cs( 15) /  0.0000000000007690431314272089939005D+00 /
  data g1cs( 16) / -0.0000000000001514696742731613519826D+00 /
  data g1cs( 17) /  0.0000000000000302892146552359684119D+00 /
  data g1cs( 18) / -0.0000000000000061399703834708825400D+00 /
  data g1cs( 19) /  0.0000000000000012600605829510933553D+00 /
  data g1cs( 20) / -0.0000000000000002615029250939483683D+00 /
  data g1cs( 21) /  0.0000000000000000548278844891796821D+00 /
  data g1cs( 22) / -0.0000000000000000116038182129526571D+00 /
  data g1cs( 23) /  0.0000000000000000024771654107129795D+00 /
  data g1cs( 24) / -0.0000000000000000005330672753223389D+00 /
  data g1cs( 25) /  0.0000000000000000001155666075598465D+00 /
  data g1cs( 26) / -0.0000000000000000000252280547744957D+00 /
  data g1cs( 27) /  0.0000000000000000000055429038550786D+00 /
  data g1cs( 28) / -0.0000000000000000000012252208421297D+00 /
  data g1cs( 29) /  0.0000000000000000000002723664318684D+00 /
  data g1cs( 30) / -0.0000000000000000000000608707831422D+00 /
  data g1cs( 31) /  0.0000000000000000000000136724874476D+00 /
  data g1cs( 32) / -0.0000000000000000000000030856626806D+00 /
  data g1cs( 33) /  0.0000000000000000000000006995212319D+00 /
  data g1cs( 34) / -0.0000000000000000000000001592587569D+00 /
  data g1cs( 35) /  0.0000000000000000000000000364051056D+00 /
  data g1cs( 36) / -0.0000000000000000000000000083539465D+00 /
  data g1cs( 37) /  0.0000000000000000000000000019240303D+00 /
  data g1cs( 38) / -0.0000000000000000000000000004446816D+00 /
  data g1cs( 39) /  0.0000000000000000000000000001031182D+00 /
  data g1cs( 40) / -0.0000000000000000000000000000239887D+00 /
  data g1cs( 41) /  0.0000000000000000000000000000055976D+00 /
  data g1cs( 42) / -0.0000000000000000000000000000013100D+00 /
  data g1cs( 43) /  0.0000000000000000000000000000003074D+00 /
  data g1cs( 44) / -0.0000000000000000000000000000000723D+00 /

  data g2cs(  1) / -0.1211802894731646263541834046858267D+00 /
  data g2cs(  2) / -0.0316761386394950286701407923505610D+00 /
  data g2cs(  3) /  0.0013383199778862680163819429492182D+00 /
  data g2cs(  4) / -0.0000895511011392252425531905069518D+00 /
  data g2cs(  5) /  0.0000079155562961718213115249467924D+00 /
  data g2cs(  6) / -0.0000008438793322241520181418982080D+00 /
  data g2cs(  7) /  0.0000001029980425677530146647227274D+00 /
  data g2cs(  8) / -0.0000000139295750605183835795834444D+00 /
  data g2cs(  9) /  0.0000000020422703959875980400677594D+00 /
  data g2cs( 10) / -0.0000000003196534694206427035434752D+00 /
  data g2cs( 11) /  0.0000000000528147832657267698615312D+00 /
  data g2cs( 12) / -0.0000000000091339554672671033735289D+00 /
  data g2cs( 13) /  0.0000000000016426251238967760444819D+00 /
  data g2cs( 14) / -0.0000000000003055897039322660002410D+00 /
  data g2cs( 15) /  0.0000000000000585655825785779717892D+00 /
  data g2cs( 16) / -0.0000000000000115229197730940120563D+00 /
  data g2cs( 17) /  0.0000000000000023209469119988537310D+00 /
  data g2cs( 18) / -0.0000000000000004774355834177535025D+00 /
  data g2cs( 19) /  0.0000000000000001000996765800180573D+00 /
  data g2cs( 20) / -0.0000000000000000213533778082256704D+00 /
  data g2cs( 21) /  0.0000000000000000046277190777367671D+00 /
  data g2cs( 22) / -0.0000000000000000010175807410227657D+00 /
  data g2cs( 23) /  0.0000000000000000002267657399884672D+00 /
  data g2cs( 24) / -0.0000000000000000000511630776076426D+00 /
  data g2cs( 25) /  0.0000000000000000000116767014913108D+00 /
  data g2cs( 26) / -0.0000000000000000000026935427672470D+00 /
  data g2cs( 27) /  0.0000000000000000000006275665841146D+00 /
  data g2cs( 28) / -0.0000000000000000000001475880557531D+00 /
  data g2cs( 29) /  0.0000000000000000000000350145314739D+00 /
  data g2cs( 30) / -0.0000000000000000000000083757732152D+00 /
  data g2cs( 31) /  0.0000000000000000000000020191815152D+00 /
  data g2cs( 32) / -0.0000000000000000000000004903567705D+00 /
  data g2cs( 33) /  0.0000000000000000000000001199123348D+00 /
  data g2cs( 34) / -0.0000000000000000000000000295170610D+00 /
  data g2cs( 35) /  0.0000000000000000000000000073113112D+00 /
  data g2cs( 36) / -0.0000000000000000000000000018217843D+00 /
  data g2cs( 37) /  0.0000000000000000000000000004565148D+00 /
  data g2cs( 38) / -0.0000000000000000000000000001150151D+00 /
  data g2cs( 39) /  0.0000000000000000000000000000291267D+00 /
  data g2cs( 40) / -0.0000000000000000000000000000074125D+00 /
  data g2cs( 41) /  0.0000000000000000000000000000018953D+00 /
  data g2cs( 42) / -0.0000000000000000000000000000004868D+00 /
  data g2cs( 43) /  0.0000000000000000000000000000001256D+00 /
  data g2cs( 44) / -0.0000000000000000000000000000000325D+00 /

  data g3cs(  1) / -0.0280574367809472928402815264335299D+00 /
  data g3cs(  2) / -0.0137271597162236975409100508089556D+00 /
  data g3cs(  3) /  0.0002894032638760296027448941273751D+00 /
  data g3cs(  4) / -0.0000114129239391197145908743622517D+00 /
  data g3cs(  5) /  0.0000006813965590726242997720207302D+00 /
  data g3cs(  6) / -0.0000000547952289604652363669058052D+00 /
  data g3cs(  7) /  0.0000000055207429918212529109406521D+00 /
  data g3cs(  8) / -0.0000000006641464199322920022491428D+00 /
  data g3cs(  9) /  0.0000000000922373663487041108564960D+00 /
  data g3cs( 10) / -0.0000000000144299088886682862611718D+00 /
  data g3cs( 11) /  0.0000000000024963904892030710248705D+00 /
  data g3cs( 12) / -0.0000000000004708240675875244722971D+00 /
  data g3cs( 13) /  0.0000000000000957217659216759988140D+00 /
  data g3cs( 14) / -0.0000000000000207889966095809030537D+00 /
  data g3cs( 15) /  0.0000000000000047875099970877431627D+00 /
  data g3cs( 16) / -0.0000000000000011619070583377173759D+00 /
  data g3cs( 17) /  0.0000000000000002956508969267836974D+00 /
  data g3cs( 18) / -0.0000000000000000785294988256492025D+00 /
  data g3cs( 19) /  0.0000000000000000216922264368256612D+00 /
  data g3cs( 20) / -0.0000000000000000062113515831676342D+00 /
  data g3cs( 21) /  0.0000000000000000018384568838450977D+00 /
  data g3cs( 22) / -0.0000000000000000005610887482137276D+00 /
  data g3cs( 23) /  0.0000000000000000001761862805280062D+00 /
  data g3cs( 24) / -0.0000000000000000000568111050541451D+00 /
  data g3cs( 25) /  0.0000000000000000000187786279582313D+00 /
  data g3cs( 26) / -0.0000000000000000000063531694151124D+00 /
  data g3cs( 27) /  0.0000000000000000000021968802368238D+00 /
  data g3cs( 28) / -0.0000000000000000000007754666550395D+00 /
  data g3cs( 29) /  0.0000000000000000000002791018356581D+00 /
  data g3cs( 30) / -0.0000000000000000000001023178525247D+00 /
  data g3cs( 31) /  0.0000000000000000000000381693403919D+00 /
  data g3cs( 32) / -0.0000000000000000000000144767895606D+00 /
  data g3cs( 33) /  0.0000000000000000000000055779512634D+00 /
  data g3cs( 34) / -0.0000000000000000000000021817239071D+00 /
  data g3cs( 35) /  0.0000000000000000000000008656646309D+00 /
  data g3cs( 36) / -0.0000000000000000000000003482157895D+00 /
  data g3cs( 37) /  0.0000000000000000000000001419188130D+00 /
  data g3cs( 38) / -0.0000000000000000000000000585714314D+00 /
  data g3cs( 39) /  0.0000000000000000000000000244660482D+00 /
  data g3cs( 40) / -0.0000000000000000000000000103387099D+00 /
  data g3cs( 41) /  0.0000000000000000000000000044177299D+00 /
  data g3cs( 42) / -0.0000000000000000000000000019080079D+00 /
  data g3cs( 43) /  0.0000000000000000000000000008326038D+00 /
  data g3cs( 44) / -0.0000000000000000000000000003669553D+00 /
  data g3cs( 45) /  0.0000000000000000000000000001632875D+00 /
  data g3cs( 46) / -0.0000000000000000000000000000733357D+00 /
  data g3cs( 47) /  0.0000000000000000000000000000332327D+00 /
  data g3cs( 48) / -0.0000000000000000000000000000151906D+00 /
  data g3cs( 49) /  0.0000000000000000000000000000070020D+00 /
  data g3cs( 50) / -0.0000000000000000000000000000032539D+00 /
  data g3cs( 51) /  0.0000000000000000000000000000015240D+00 /
  data g3cs( 52) / -0.0000000000000000000000000000007193D+00 /
  data g3cs( 53) /  0.0000000000000000000000000000003420D+00 /
  data g3cs( 54) / -0.0000000000000000000000000000001638D+00 /
  data g3cs( 55) /  0.0000000000000000000000000000000790D+00 /
  data g3cs( 56) / -0.0000000000000000000000000000000383D+00 /

  data nf1 / 0 /
  data nf2 / 0 /
  data ng1 / 0 /
  data ng2 / 0 /
  data ng3 / 0 /
  data xbig / 0.0D+00 /
  data xbnd / 0.0D+00 /
  data xbndg / 0.0D+00 /
  data xmaxf / 0.0D+00 /
  data xmaxg / 0.0D+00 /

  if ( nf1 == 0 ) then
    tol = 0.1D+00 * r8_mach ( 3 )
    nf1 = r8_inits ( f1cs, 43, tol )
    nf2 = r8_inits ( f2cs, 99, tol )
    ng1 = r8_inits ( g1cs, 44, tol )
    ng2 = r8_inits ( g2cs, 44, tol )
    ng3 = r8_inits ( g3cs, 56, tol )
    xbig = sqrt ( 1.0D+00 / r8_mach ( 3 ) )
    xmaxf = exp ( min ( - log ( r8_mach ( 1 ) ), &
      log ( r8_mach ( 2 ) ) ) - 0.01D+00 )
    xmaxg = 1.0D+00 / sqrt ( r8_mach ( 1 ) )
    xbnd = sqrt ( 50.0D+00 )
    xbndg = sqrt ( 200.0D+00 )
  end if

  if ( x < 4.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_SIFG - Fatal error!'
    write ( *, '(a)' ) '  Approximation invalid for X < 4.'
    stop
  else if ( x <= xbnd ) then
    f = ( 1.0D+00 &
      + r8_csevl ( ( 1.0D+00 / x / x - 0.04125D+00 ) &
      / 0.02125D+00, f1cs, nf1 ) ) / x
    g = ( 1.0D+00 &
      + r8_csevl ( ( 1.0D+00 / x / x - 0.04125D+00 ) &
      / 0.02125D+00, g1cs, ng1 ) ) / x / x
  else if ( x <= xbig ) then
    f = ( 1.0D+00 &
      + r8_csevl ( 100.0D+00 / x / x - 1.0D+00, f2cs, nf2 ) ) / x
    if ( x <= xbndg ) then 
      g = ( 1.0D+00 &
        + r8_csevl ( ( 10000.0D+00 / x / x - 125.0D+00 ) &
        / 75.0D+00, g2cs, ng2 ) ) / x / x
    else
      g = ( 1.0D+00 &
        + r8_csevl ( 400.0D+00 / x / x - 1.0D+00, g3cs, ng3 ) ) / x / x
    end if
  else
    if ( x < xmaxf ) then
      f = 1.0D+00 / x
    else
      f = 0.0D+00
    end if
    if ( x < xmaxg ) then
      g = 1.0D+00 / x / x
    else
      g = 0.0D+00
    end if
  end if

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value =  0 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 8 ) R8_SIGN, the sign of X:
!
  implicit none

  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

  return
end
function r8_sin ( x )

!*****************************************************************************80
!
!! R8_SIN evaluates the sine of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_SIN, the sine of X.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) ntsn
  real ( kind = 8 ) pi2rec
  real ( kind = 8 ) pihi
  real ( kind = 8 ) pilo
  real ( kind = 8 ) pirec
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_sin
  real ( kind = 8 ) sgn
  real ( kind = 8 ) sincs(15)
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xn
  real ( kind = 8 ) xsml
  real ( kind = 8 ) xwarn
  real ( kind = 8 ) y

  save ntsn
  save pi2rec
  save pihi
  save pilo
  save pirec
  save sincs
  save xmax
  save xsml
  save xwarn

  data sincs( 1) / -0.374991154955873175839919279977323464D+00 /
  data sincs( 2) / -0.181603155237250201863830316158004754D+00 /
  data sincs( 3) / +0.005804709274598633559427341722857921D+00 /
  data sincs( 4) / -0.000086954311779340757113212316353178D+00 /
  data sincs( 5) / +0.000000754370148088851481006839927030D+00 /
  data sincs( 6) / -0.000000004267129665055961107126829906D+00 /
  data sincs( 7) / +0.000000000016980422945488168181824792D+00 /
  data sincs( 8) / -0.000000000000050120578889961870929524D+00 /
  data sincs( 9) / +0.000000000000000114101026680010675628D+00 /
  data sincs(10) / -0.000000000000000000206437504424783134D+00 /
  data sincs(11) / +0.000000000000000000000303969595918706D+00 /
  data sincs(12) / -0.000000000000000000000000371357734157D+00 /
  data sincs(13) / +0.000000000000000000000000000382486123D+00 /
  data sincs(14) / -0.000000000000000000000000000000336623D+00 /
  data sincs(15) / +0.000000000000000000000000000000000256D+00 /
!
!  pihi + pilo = pi.  pihi is exactly representable on all machines
!  with at least 8 bits of precision.  whether it is exactly
!  represented depends on the compiler.  this routine is more
!  accurate if it is exactly represented.
!
  data ntsn / 0 /
  data pi2rec / 0.63661977236758134307553505349006D+00 /
  data pihi / 3.140625D+00 /
  data pilo / 9.6765358979323846264338327950288D-04 /
  data pirec / 0.31830988618379067153776752674503D+00 /
  data xmax / 0.0D+00 /
  data xsml / 0.0D+00 /
  data xwarn / 0.0D+00 /

  if ( ntsn == 0 ) then
    ntsn = r8_inits ( sincs, 15, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( 2.0D+00 * r8_mach ( 3 ) )
    xmax = 1.0D+00 / r8_mach( 4 )
    xwarn = sqrt ( xmax )
  end if

  y = abs ( x )

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_SIN - Warning!'
    write ( *, '(a)' ) '  No precision because |X| is big.'
    r8_sin = 0.0D+00
    return
  end if

  if ( xwarn < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_SIN - Warning!'
    write ( *, '(a)' ) '  Answer < half precision because |X| is big.'
  end if

  r8_sin = x
  if ( y < xsml ) then
    return
  end if

  xn = aint ( y * pirec + 0.5D+00 )
  n2 = int ( mod ( xn, 2.0D+00 ) + 0.5D+00 )

  sgn = x
  if ( n2 /= 0 ) then
    sgn = - sgn
  end if

  f = ( y - xn * pihi ) - xn * pilo

  xn = 2.0D+00 * ( f * pi2rec ) * ( f * pi2rec ) - 1.0D+00

  r8_sin = f + f * r8_csevl ( xn, sincs, ntsn )

  if ( sgn < 0.0D+00 ) then
    r8_sin = - r8_sin
  end if

  if ( r8_sin < - 1.0D+00 ) then
    r8_sin = - 1.0D+00
  else if ( 1.0D+00 < r8_sin ) then
    r8_sin = + 1.0D+00
  end if

  return
end
function r8_sin_deg ( x )

!*****************************************************************************80
!
!! R8_SIN_DEG evaluates the sine of an R8 argument in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument in degrees.
!
!    Output, real ( kind = 8 ) R8_SIN_DEG, the sine of X.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_sin_deg
  real ( kind = 8 ), parameter :: raddeg &
    = 0.017453292519943295769236907684886D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = sin ( raddeg * x )

  if ( mod ( x, 90.0D+00 ) == 0.0D+00 ) then

    n = int ( abs ( x ) / 90.0D+00 + 0.5D+00 )
    n = mod ( n, 2 )

    if ( n == 0 ) then
      value = 0.0D+00
    else if ( value < 0.0D+00 ) then
      value = - 1.0D+00
    else
      value = + 1.0D+00
    end if

  end if

  r8_sin_deg = value

  return
end
function r8_sinh ( x )

!*****************************************************************************80
!
!! R8_SINH evaluates the hyperbolic sine of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_SINH, the hyperbolic sine of X.
!
  implicit none

  integer ( kind = 4 ) nterms
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_sinh
  real ( kind = 8 ) sinhcs(13)
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax

  save nterms
  save sinhcs
  save ymax

  data sinhcs(  1) / +0.17304219404717963167588384698501D+00 /
  data sinhcs(  2) / +0.87594221922760477154900263454440D-01 /
  data sinhcs(  3) / +0.10794777745671327502427270651579D-02 /
  data sinhcs(  4) / +0.63748492607547504815685554571850D-05 /
  data sinhcs(  5) / +0.22023664049230530159190496019502D-07 /
  data sinhcs(  6) / +0.49879401804158493149425807203661D-10 /
  data sinhcs(  7) / +0.79730535541157304814411480441186D-13 /
  data sinhcs(  8) / +0.94731587130725443342927317226666D-16 /
  data sinhcs(  9) / +0.86934920504480078871023098666666D-19 /
  data sinhcs( 10) / +0.63469394403318040457397333333333D-22 /
  data sinhcs( 11) / +0.37740337870858485738666666666666D-25 /
  data sinhcs( 12) / +0.18630213719570056533333333333333D-28 /
  data sinhcs( 13) / +0.77568437166506666666666666666666D-32 /

  data nterms / 0 /
  data ymax / 0.0D+00 /

  if ( nterms == 0 ) then
    nterms = r8_inits ( sinhcs, 13, 0.1D+00 * r8_mach ( 3 ) )
    sqeps = sqrt ( 6.0D+00 * r8_mach ( 3 ) )
    ymax = 1.0D+00 / sqrt ( r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then

    value = x

  else if ( y <= 1.0D+00 ) then

    value = x * ( 1.0D+00 &
      + r8_csevl ( 2.0D+00 * x * x - 1.0D+00, sinhcs, nterms ) )

  else

    y = exp ( y )

    if ( ymax <= y ) then 
      value = 0.5D+00 * y
    else
      value = 0.5D+00 * ( y - 1.0D+00 / y )
    end if

    if ( x < 0.0D+00 ) then
      value = - value
    end if

  end if

  r8_sinh = value

  return
end
function r8_spence ( x )

!*****************************************************************************80
!
!! R8_SPENCE evaluates a form of Spence's function for an R8 argument.
!
!  Discussion:
!
!    This function evaluates a form of Spence's function defined by
!
!      f(x) = Integral ( 0 <= y <= x ) - log ( abs ( 1 - y ) ) / y dy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions, page 1004,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!    K Mitchell,
!    Tables of the function Integral ( 0 < y < x ) - log | 1 - y | dy / y
!    with an account of some properties of this and related functions,
!    Philosophical Magazine,
!    Volume 40, pages 351-368, 1949.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_SPENCE, Spence's function evaluated at X.
!
  implicit none

  real ( kind = 8 ) aln
  integer ( kind = 4 ) nspenc
  real ( kind = 8 ) pi26
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_spence
  real ( kind = 8 ) spencs(38)
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig

  save nspenc
  save pi26
  save spencs
  save xbig

  data spencs(  1) / +0.1527365598892405872946684910028D+00 /
  data spencs(  2) / +0.8169658058051014403501838185271D-01 /
  data spencs(  3) / +0.5814157140778730872977350641182D-02 /
  data spencs(  4) / +0.5371619814541527542247889005319D-03 /
  data spencs(  5) / +0.5724704675185826233210603054782D-04 /
  data spencs(  6) / +0.6674546121649336343607835438589D-05 /
  data spencs(  7) / +0.8276467339715676981584391689011D-06 /
  data spencs(  8) / +0.1073315673030678951270005873354D-06 /
  data spencs(  9) / +0.1440077294303239402334590331513D-07 /
  data spencs( 10) / +0.1984442029965906367898877139608D-08 /
  data spencs( 11) / +0.2794005822163638720201994821615D-09 /
  data spencs( 12) / +0.4003991310883311823072580445908D-10 /
  data spencs( 13) / +0.5823462892044638471368135835757D-11 /
  data spencs( 14) / +0.8576708692638689278097914771224D-12 /
  data spencs( 15) / +0.1276862586280193045989483033433D-12 /
  data spencs( 16) / +0.1918826209042517081162380416062D-13 /
  data spencs( 17) / +0.2907319206977138177795799719673D-14 /
  data spencs( 18) / +0.4437112685276780462557473641745D-15 /
  data spencs( 19) / +0.6815727787414599527867359135607D-16 /
  data spencs( 20) / +0.1053017386015574429547019416644D-16 /
  data spencs( 21) / +0.1635389806752377100051821734570D-17 /
  data spencs( 22) / +0.2551852874940463932310901642581D-18 /
  data spencs( 23) / +0.3999020621999360112770470379519D-19 /
  data spencs( 24) / +0.6291501645216811876514149171199D-20 /
  data spencs( 25) / +0.9933827435675677643803887752533D-21 /
  data spencs( 26) / +0.1573679570749964816721763805866D-21 /
  data spencs( 27) / +0.2500595316849476129369270954666D-22 /
  data spencs( 28) / +0.3984740918383811139210663253333D-23 /
  data spencs( 29) / +0.6366473210082843892691326293333D-24 /
  data spencs( 30) / +0.1019674287239678367077061973333D-24 /
  data spencs( 31) / +0.1636881058913518841111074133333D-25 /
  data spencs( 32) / +0.2633310439417650117345279999999D-26 /
  data spencs( 33) / +0.4244811560123976817224362666666D-27 /
  data spencs( 34) / +0.6855411983680052916824746666666D-28 /
  data spencs( 35) / +0.1109122433438056434018986666666D-28 /
  data spencs( 36) / +0.1797431304999891457365333333333D-29 /
  data spencs( 37) / +0.2917505845976095173290666666666D-30 /
  data spencs( 38) / +0.4742646808928671061333333333333D-31 /

  data nspenc / 0 /
  data pi26 / +1.644934066848226436472415166646025189219D+00 /
  data xbig / 0.0D+00 /

  if ( nspenc == 0 ) then
    nspenc = r8_inits ( spencs, 38, 0.1D+00 * r8_mach ( 3 ) )
    xbig = 1.0D+00 / r8_mach ( 3 )
  end if

  if ( x <= - xbig ) then

    aln = log ( 1.0D+00 - x )
    r8_spence = - pi26  - 0.5D+00 * aln * ( 2.0D+00 * log ( - x ) - aln )

  else if ( x <= - 1.0D+00 ) then

    aln = log ( 1.0D+00 - x )

    r8_spence = - pi26 - 0.5D+00 * aln * ( 2.0D+00 &
      * log ( - x ) - aln ) + ( 1.0D+00 + r8_csevl ( &
      4.0D+00 / ( 1.0D+00 - x ) - 1.0D+00, spencs, nspenc ) ) &
     / ( 1.0D+00 - x )

  else if ( x <= 0.0D+00 ) then

    r8_spence = - 0.5D+00 * log ( 1.0D+00 - x ) &
      * log ( 1.0D+00 - x ) - x * ( 1.0D+00 + r8_csevl ( &
      4.0D+00 * x / ( x - 1.0D+00 ) - 1.0D+00, spencs, nspenc ) ) &
      / ( x - 1.0D+00 )

  else if ( x <= 0.5D+00 ) then

    r8_spence = x * ( 1.0D+00 + r8_csevl ( 4.0D+00 * x - 1.0D+00, &
      spencs, nspenc ) )

  else if ( x < 1.0D+00 ) then

    r8_spence = pi26 - log ( x ) * log ( 1.0D+00 - x ) &
      - ( 1.0D+00 - x ) * ( 1.0D+00 + r8_csevl ( 4.0D+00 &
      * ( 1.0D+00 - x ) - 1.0D+00, spencs, nspenc ) )

  else if ( x == 1.0D+00 ) then

    r8_spence = pi26

  else if ( x <= 2.0D+00 ) then

    r8_spence = pi26 - 0.5D+00 * log ( x ) &
      * log ( ( x - 1.0D+00 ) * ( x - 1.0D+00 ) / x ) &
      + ( x - 1.0D+00 ) * ( 1.0D+00 + r8_csevl ( 4.0D+00 &
      * ( x - 1.0D+00 ) / x - 1.0D+00, spencs, nspenc ) ) / x

  else if ( x < xbig ) then

    r8_spence = 2.0D+00 * pi26 - 0.5D+00 * log ( x ) * log ( x ) &
      - ( 1.0D+00 + r8_csevl ( 4.0D+00 / x - 1.0D+00, spencs, &
      nspenc ) ) / x

  else

    r8_spence = 2.0D+00 * pi26 - 0.5D+00 * log ( x ) * log ( x )

  end if

  return
end
function r8_sqrt ( x )

!*****************************************************************************80
!
!! R8_SQRT computes the square root of an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose square root is desired.
!
!    Output, real ( kind = 8 ) R8_SQRT, the square root of X.
!
  implicit none

  integer ( kind = 4 ) irem
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) ixpnt
  integer ( kind = 4 ) n
  integer ( kind = 4 ) niter
  real ( kind = 8 ) r8_log
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_pak
  real ( kind = 8 ) r8_sqrt
  real ( kind = 8 ) sqrt2(3)
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  save niter
  save sqrt2

  data sqrt2(1) / 0.70710678118654752440084436210485D+00 /
  data sqrt2(2) / 1.0D+00 /
  data sqrt2(3) / 1.41421356237309504880168872420970D+00 /

  data niter / 0 /

  if ( niter == 0 ) then
    niter = 1.443D+00 * r8_log ( - 0.104D+00 &
      * r8_log ( 0.1D+00 * r8_mach ( 3 ) ) ) + 1.0D+00
  end if

  if ( x < 0.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_SQRT - Fatal error!'
    write ( *, '(a)' ) '  X < 0.0'
    stop

  else if ( x == 0.0D+00 ) then

    r8_sqrt = 0.0D+00

  else

    call r8_upak ( x, y, n )
    ixpnt = n / 2
    irem = n - 2 * ixpnt + 2
    value = 0.261599D+00 + y * ( 1.114292D+00 &
      + y * ( -0.516888D+00 + y * 0.141067D+00 ) )

    do iter = 1, niter
      value = value + 0.5D+00 * ( y - value * value ) / value
    end do

    r8_sqrt = r8_pak ( sqrt2(irem) * value, ixpnt )

  end if

  return
end
function r8_tan ( x )

!*****************************************************************************80
!
!! R8_TAN evaluates the tangent of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_TAN, the tangent of X.
!
  implicit none

  real ( kind = 8 ) ainty
  real ( kind = 8 ) ainty2
  integer ( kind = 4 ) ifn
  integer ( kind = 4 ) nterms
  real ( kind = 8 ) pi2rec
  real ( kind = 8 ) prodbg
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_tan
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) tancs(19)
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y
  real ( kind = 8 ) yrem

  save nterms
  save pi2rec
  save tancs
  save xmax
  save xsml

  data tancs(  1) / +0.22627932763129357846578636531752D+00 /
  data tancs(  2) / +0.43017913146548961775583410748067D-01 /
  data tancs(  3) / +0.68544610682565088756929473623461D-03 /
  data tancs(  4) / +0.11045326947597098383578849369696D-04 /
  data tancs(  5) / +0.17817477903926312943238512588940D-06 /
  data tancs(  6) / +0.28744968582365265947529646832471D-08 /
  data tancs(  7) / +0.46374854195902995494137478234363D-10 /
  data tancs(  8) / +0.74817609041556138502341633308215D-12 /
  data tancs(  9) / +0.12070497002957544801644516947824D-13 /
  data tancs( 10) / +0.19473610812823019305513858584533D-15 /
  data tancs( 11) / +0.31417224874732446504614586026666D-17 /
  data tancs( 12) / +0.50686132555800153941904891733333D-19 /
  data tancs( 13) / +0.81773105159836540043979946666666D-21 /
  data tancs( 14) / +0.13192643412147384408951466666666D-22 /
  data tancs( 15) / +0.21283995497042377309866666666666D-24 /
  data tancs( 16) / +0.34337960192345945292800000000000D-26 /
  data tancs( 17) / +0.55398222121173811200000000000000D-28 /
  data tancs( 18) / +0.89375227794352810666666666666666D-30 /
  data tancs( 19) / +0.14419111371369130666666666666666D-31 /

  data nterms / 0 /
  data pi2rec / 0.011619772367581343075535053490057D+00 /
  data xmax / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( nterms == 0 ) then
    nterms = r8_inits ( tancs, 19, 0.1D+00 * r8_mach ( 3 ) )
    xmax = 1.0D+00 / r8_mach ( 4 )
    xsml = sqrt ( 3.0D+00 * r8_mach ( 3 ) )
    sqeps = sqrt ( r8_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_TAN - Warning'
    write ( *, '(a)' ) '  No precision because |X| is big.'
    r8_tan = 0.0D+00
    return
  end if
!
!  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
!  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
!  = aint(.625*y) + aint(z) + rem(z)
!
  ainty = aint ( y )
  yrem = y - ainty
  prodbg = 0.625D+00 * ainty
  ainty = aint ( prodbg )
  y = ( prodbg - ainty ) + 0.625D+00 * yrem + pi2rec * y
  ainty2 = aint ( y )
  ainty = ainty + ainty2
  y = y - ainty2

  ifn = int ( mod ( ainty, 2.0D+00 ) )

  if ( ifn == 1 ) then
    y = 1.0D+00 - y
  end if

  if ( 1.0D+00 - y < abs ( x ) * sqeps ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_TAN - Warning!'
    write ( *, '(a)' ) '  Answer < half precision.'
    write ( *, '(a)' ) '  |X| big or X near pi/2 or 3*pi/2.'
  end if

  if ( y == 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_TAN - Fatal error!'
    write ( *, '(a)' ) '  X is pi/2 or 3*pi/2.'
    value = 0.0D+00
    stop
  end if

  if ( y <= 0.25D+00 ) then

    value = y
    if ( xsml < y ) then
      value = y * ( 1.5D+00 + r8_csevl ( 32.0D+00 * y * y - 1.0D+00, &
        tancs, nterms ) )
    end if

  else if ( y <= 0.5D+00 ) then

    value = 0.5D+00 * y * ( 1.5D+00 + r8_csevl ( &
      8.0D+00 * y * y - 1.0D+00, tancs, nterms ) )
    value = 2.0D+00 * value / ( 1.0D+00 - value * value )

  else

    value = 0.25D+00 * y * ( 1.5D+00 + r8_csevl ( &
      2.0D+00 * y * y - 1.0D+00, tancs, nterms ) )
    value = 2.0D+00 * value / ( 1.0D+00 - value * value )
    value = 2.0D+00 * value / ( 1.0D+00 - value * value )

  end if

  if ( x < 0.0D+00 ) then
    value = - abs ( value )
  else if ( 0.0D+00 < x ) then
    value = + abs ( value )
  end if

  if ( ifn == 1 ) then
    value = - value
  end if

  r8_tan = value

  return
end
function r8_tanh ( x )

!*****************************************************************************80
!
!! R8_TANH evaluates the hyperbolic tangent of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_TANH, the hyperbolic tangent of X.
!
  implicit none

  integer ( kind = 4 ) nterms
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) r8_tanh
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ) tanhcs(31)
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) y
  real ( kind = 8 ) yrec

  save nterms
  save tanhcs
  save xmax

  data tanhcs(  1) / -0.25828756643634710438338151450605D+00 /
  data tanhcs(  2) / -0.11836106330053496535383671940204D+00 /
  data tanhcs(  3) / +0.98694426480063988762827307999681D-02 /
  data tanhcs(  4) / -0.83579866234458257836163690398638D-03 /
  data tanhcs(  5) / +0.70904321198943582626778034363413D-04 /
  data tanhcs(  6) / -0.60164243181207040390743479001010D-05 /
  data tanhcs(  7) / +0.51052419080064402965136297723411D-06 /
  data tanhcs(  8) / -0.43320729077584087216545467387192D-07 /
  data tanhcs(  9) / +0.36759990553445306144930076233714D-08 /
  data tanhcs( 10) / -0.31192849612492011117215651480953D-09 /
  data tanhcs( 11) / +0.26468828199718962579377758445381D-10 /
  data tanhcs( 12) / -0.22460239307504140621870997006196D-11 /
  data tanhcs( 13) / +0.19058733768288196054319468396139D-12 /
  data tanhcs( 14) / -0.16172371446432292391330769279701D-13 /
  data tanhcs( 15) / +0.13723136142294289632897761289386D-14 /
  data tanhcs( 16) / -0.11644826870554194634439647293781D-15 /
  data tanhcs( 17) / +0.98812684971669738285540514338133D-17 /
  data tanhcs( 18) / -0.83847933677744865122269229055999D-18 /
  data tanhcs( 19) / +0.71149528869124351310723506176000D-19 /
  data tanhcs( 20) / -0.60374242229442045413288837119999D-20 /
  data tanhcs( 21) / +0.51230825877768084883404663466666D-21 /
  data tanhcs( 22) / -0.43472140157782110106047829333333D-22 /
  data tanhcs( 23) / +0.36888473639031328479423146666666D-23 /
  data tanhcs( 24) / -0.31301874774939399883325439999999D-24 /
  data tanhcs( 25) / +0.26561342006551994468488533333333D-25 /
  data tanhcs( 26) / -0.22538742304145029883494399999999D-26 /
  data tanhcs( 27) / +0.19125347827973995102208000000000D-27 /
  data tanhcs( 28) / -0.16228897096543663117653333333333D-28 /
  data tanhcs( 29) / +0.13771101229854738786986666666666D-29 /
  data tanhcs( 30) / -0.11685527840188950118399999999999D-30 /
  data tanhcs( 31) / +0.99158055384640389120000000000000D-32 /

  data nterms / 0 /
  data xmax / 0.0D+00 /

  if ( nterms == 0 ) then
    nterms = r8_inits ( tanhcs, 31, 0.1D+00 * r8_mach ( 3 ) )
    sqeps = sqrt ( 3.0D+00 * r8_mach ( 3 ) )
    xmax = - 0.5D+00 * log ( r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then

    value = x

  else if ( y <= 1.0D+00 ) then

    value = x * ( 1.0D+00 &
      + r8_csevl ( 2.0D+00 * x * x - 1.0D+00, tanhcs, nterms ) )

  else if ( y <= xmax ) then

    y = exp ( y )
    yrec = 1.0D+00 / y
    value = ( y - yrec ) / ( y + yrec )

    if ( x < 0.0D+00 ) then
      value = - value
    end if

  else

    if ( x < 0.0D+00 ) then
      value = - 1.0D+00
    else
      value = + 1.0D+00
    end if

  end if

  r8_tanh = value

  return
end
subroutine r8_upak ( x, y, n )

!*****************************************************************************80
!
!! R8_UPAK unpacks an R8 into a mantissa and exponent.
!
!  Discussion:
!
!    This function unpacks a floating point number x so that
!
!      x = y * 2.0^n
!
!    where
!
!      0.5 <= abs ( y ) < 1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be unpacked.
!
!    Output, real ( kind = 8 ) Y, the mantissa.
!
!    Output, integer ( kind = 4 ) N, the exponent.
!
  implicit none

  real ( kind = 8 ) absx
  integer ( kind = 4 ) n
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  absx = abs ( x )
  n = 0
  y = 0.0D+00

  if ( x == 0.0D+00 ) then
    return
  end if

  do while ( absx < 0.5D+00 )
    n = n - 1
    absx = absx * 2.0D+00
  end do

  do while ( 1.0D+00 <= absx )
    n = n + 1
    absx = absx * 0.5D+00
  end do

  if ( x < 0.0D+00 ) then
    y = - absx
  else
    y = + absx
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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
