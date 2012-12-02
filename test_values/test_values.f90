subroutine abram0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ABRAM0_VALUES returns some values of the Abramowitz0 function.
!
!  Discussion:
!
!    The function is defined by:
!
!      ABRAM0(X) = integral ( 0 <= T < +oo ) exp ( -T * T - X / T ) dT
!
!    The data was reported by McLeod.
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.87377726306985360531D+00, &
    0.84721859650456925922D+00, &
    0.77288934483988301615D+00, &
    0.59684345853450151603D+00, &
    0.29871735283675888392D+00, &
    0.15004596450516388138D+00, &
    0.11114662419157955096D+00, &
    0.83909567153151897766D-01, &
    0.56552321717943417515D-01, &
    0.49876496603033790206D-01, &
    0.44100889219762791328D-01, &
    0.19738535180254062496D-01, &
    0.86193088287161479900D-02, &
    0.40224788162540127227D-02, &
    0.19718658458164884826D-02, &
    0.10045868340133538505D-02, &
    0.15726917263304498649D-03, &
    0.10352666912350263437D-04, &
    0.91229759190956745069D-06, &
    0.25628287737952698742D-09 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0019531250D+00, &
    0.0078125000D+00, &
    0.0312500000D+00, &
    0.1250000000D+00, &
    0.5000000000D+00, &
    1.0000000000D+00, &
    1.2500000000D+00, &
    1.5000000000D+00, &
    1.8750000000D+00, &
    2.0000000000D+00, &
    2.1250000000D+00, &
    3.0000000000D+00, &
    4.0000000000D+00, &
    5.0000000000D+00, &
    6.0000000000D+00, &
    7.0000000000D+00, &
    10.0000000000D+00, &
    15.0000000000D+00, &
    20.0000000000D+00, &
    40.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine abram1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ABRAM1_VALUES returns some values of the Abramowitz1 function.
!
!  Discussion:
!
!    The function is defined by:
!
!      ABRAM1(x) = integral ( 0 <= t < +oo ) t * exp ( -t^2 - x / t ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.49828219848799921792D+00, &
    0.49324391773047288556D+00, &
    0.47431612784691234649D+00, &
    0.41095983258760410149D+00, &
    0.25317617388227035867D+00, &
    0.14656338138597777543D+00, &
    0.11421547056018366587D+00, &
    0.90026307383483764795D-01, &
    0.64088214170742303375D-01, &
    0.57446614314166191085D-01, &
    0.51581624564800730959D-01, &
    0.25263719555776416016D-01, &
    0.11930803330196594536D-01, &
    0.59270542280915272465D-02, &
    0.30609215358017829567D-02, &
    0.16307382136979552833D-02, &
    0.28371851916959455295D-03, &
    0.21122150121323238154D-04, &
    0.20344578892601627337D-05, &
    0.71116517236209642290D-09 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0019531250D+00, &
    0.0078125000D+00, &
    0.0312500000D+00, &
    0.1250000000D+00, &
    0.5000000000D+00, &
    1.0000000000D+00, &
    1.2500000000D+00, &
    1.5000000000D+00, &
    1.8750000000D+00, &
    2.0000000000D+00, &
    2.1250000000D+00, &
    3.0000000000D+00, &
    4.0000000000D+00, &
    5.0000000000D+00, &
    6.0000000000D+00, &
    7.0000000000D+00, &
    10.0000000000D+00, &
    15.0000000000D+00, &
    20.0000000000D+00, &
    40.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine abram2_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ABRAM2_VALUES returns some values of the Abramowitz2 function.
!
!  Discussion:
!
!    The function is defined by:
!
!      ABRAM2(x) = integral ( 0 <= t < +oo ) t^2 * exp( -t^2 - x / t ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.44213858162107913430D+00, &
    0.43923379545684026308D+00, &
    0.42789857297092602234D+00, &
    0.38652825661854504406D+00, &
    0.26538204413231368110D+00, &
    0.16848734838334595000D+00, &
    0.13609200032513227112D+00, &
    0.11070330027727917352D+00, &
    0.82126019995530382267D-01, &
    0.74538781999594581763D-01, &
    0.67732034377612811390D-01, &
    0.35641808698811851022D-01, &
    0.17956589956618269083D-01, &
    0.94058737143575370625D-02, &
    0.50809356204299213556D-02, &
    0.28149565414209719359D-02, &
    0.53808696422559303431D-03, &
    0.44821756380146327259D-04, &
    0.46890678427324100410D-05, &
    0.20161544850996420504D-08 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0019531250D+00, &
    0.0078125000D+00, &
    0.0312500000D+00, &
    0.1250000000D+00, &
    0.5000000000D+00, &
    1.0000000000D+00, &
    1.2500000000D+00, &
    1.5000000000D+00, &
    1.8750000000D+00, &
    2.0000000000D+00, &
    2.1250000000D+00, &
    3.0000000000D+00, &
    4.0000000000D+00, &
    5.0000000000D+00, &
    6.0000000000D+00, &
    7.0000000000D+00, &
    10.0000000000D+00, &
    15.0000000000D+00, &
    20.0000000000D+00, &
    40.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine agm_values ( n_data, a, b, fx )

!*****************************************************************************80
!
!! AGM_VALUES returns some values of the AGM.
!
!  Discussion:
!
!    The AGM is defined for nonnegative A and B.
!
!    The AGM of numbers A and B is defined by setting
!
!      A(0) = A,
!      B(0) = B
!
!      A(N+1) = ( A(N) + B(N) ) / 2
!      B(N+1) = sqrt ( A(N) * B(N) )
!
!    The two sequences both converge to AGM(A,B).
!
!    In Mathematica, the AGM can be evaluated by
!
!      ArithmeticGeometricMean [ a, b ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the argument ofs the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     22.0D+00, &
     83.0D+00, &
     42.0D+00, &
     26.0D+00, &
      4.0D+00, &
      6.0D+00, &
     40.0D+00, &
     80.0D+00, &
     90.0D+00, &
      9.0D+00, &
     53.0D+00, &
      1.0D+00, &
      1.0D+00, &
      1.0D+00, &
      1.5D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     96.0D+00, &
     56.0D+00, &
      7.0D+00, &
     11.0D+00, &
     63.0D+00, &
     45.0D+00, &
     75.0D+00, &
      0.0D+00, &
     35.0D+00, &
      1.0D+00, &
     53.0D+00, &
      2.0D+00, &
      4.0D+00, &
      8.0D+00, &
      8.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     52.274641198704240049D+00, &
     68.836530059858524345D+00, &
     20.659301196734009322D+00, &
     17.696854873743648823D+00, &
     23.867049721753300163D+00, &
     20.717015982805991662D+00, &
     56.127842255616681863D+00, &
      0.000000000000000000D+00, &
     59.269565081229636528D+00, &
     3.9362355036495554780D+00, &
     53.000000000000000000D+00, &
     1.4567910310469068692D+00, &
     2.2430285802876025701D+00, &
     3.6157561775973627487D+00, &
     4.0816924080221632670D+00 /)
  integer ( kind = 4 ) n_data

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine airy_ai_values ( n_data, x, ai )

!*****************************************************************************80
!
!! AIRY_AI_VALUES returns some values of the Airy Ai(x) function.
!
!  Discussion:
!
!    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
!    solutions of the differential equation:
!
!      W'' - X * W = 0
!
!    In Mathematica, the function can be evaluated by:
!
!      AiryAi[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) AI, the value of the Airy AI function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) ai
  real ( kind = 8 ), save, dimension ( n_max ) :: ai_vec = (/ &
    0.3550280538878172D+00, &
    0.3292031299435381D+00, &
    0.3037031542863820D+00, &
    0.2788064819550049D+00, &
    0.2547423542956763D+00, &
    0.2316936064808335D+00, &
    0.2098000616663795D+00, &
    0.1891624003981501D+00, &
    0.1698463174443649D+00, &
    0.1518868036405444D+00, &
    0.1352924163128814D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    ai = 0.0D+00
  else
    x = x_vec(n_data)
    ai = ai_vec(n_data)
  end if

  return
end
subroutine airy_ai_int_values ( n_data, x, fx )

!*****************************************************************************80
!
!! AIRY_AI_INT_VALUES returns some values of the integral of the Airy function.
!
!  Discussion:
!
!    The function is defined by:
!
!      AIRY_AI_INT(x) = integral ( 0 <= t <= x ) Ai(t) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.75228838916610124300D+00, &
    -0.57348350185854889466D+00, &
    -0.76569840313421291743D+00, &
    -0.65181015505382467421D+00, &
    -0.55881974894471876922D+00, &
    -0.56902352870716815309D+00, &
    -0.47800749642926168100D+00, &
    -0.46567398346706861416D+00, &
    -0.96783140945618013679D-01, &
    -0.34683049857035607494D-03, &
     0.34658366917927930790D-03, &
     0.27657581846051227124D-02, &
     0.14595330491185717833D+00, &
     0.23631734191710977960D+00, &
     0.33289264538612212697D+00, &
     0.33318759129779422976D+00, &
     0.33332945170523851439D+00, &
     0.33333331724248357420D+00, &
     0.33333333329916901594D+00, &
     0.33333333333329380187D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -12.0000000000D+00, &
    -11.0000000000D+00, &
    -10.0000000000D+00, &
     -9.5000000000D+00, &
     -9.0000000000D+00, &
     -6.5000000000D+00, &
     -4.0000000000D+00, &
     -1.0000000000D+00, &
     -0.2500000000D+00, &
     -0.0009765625D+00, &
      0.0009765625D+00, &
      0.0078125000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      4.0000000000D+00, &
      4.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     12.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine airy_ai_prime_values ( n_data, x, aip )

!*****************************************************************************80
!
!! AIRY_AI_PRIME_VALUES returns some values of the Airy function Ai'(x).
!
!  Discussion:
!
!    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
!    solutions of the differential equation:
!
!      W'' - X * W = 0
!
!    In Mathematica, the function can be evaluated by:
!
!      AiryAiPrime[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) AIP, the derivative of the Airy AI function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) aip
  real ( kind = 8 ), save, dimension ( n_max ) :: aip_vec = (/ &
    -0.2588194037928068D+00, &
    -0.2571304219075862D+00, &
    -0.2524054702856195D+00, &
    -0.2451463642190548D+00, &
    -0.2358320344192082D+00, &
    -0.2249105326646839D+00, &
    -0.2127932593891585D+00, &
    -0.1998511915822805D+00, &
    -0.1864128638072717D+00, &
    -0.1727638434616347D+00, &
    -0.1591474412967932D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    aip = 0.0D+00
  else
    x = x_vec(n_data)
    aip = aip_vec(n_data)
  end if

  return
end
subroutine airy_bi_values ( n_data, x, bi )

!*****************************************************************************80
!
!! AIRY_BI_VALUES returns some values of the Airy Bi(x) function.
!
!  Discussion:
!
!    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
!    solutions of the differential equation:
!
!      W'' - X * W = 0
!
!    In Mathematica, the function can be evaluated by:
!
!      AiryBi[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BI, the value of the Airy BI function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) bi
  real ( kind = 8 ), save, dimension ( n_max ) :: bi_vec = (/ &
    0.6149266274460007D+00, &
    0.6598616901941892D+00, &
    0.7054642029186612D+00, &
    0.7524855850873156D+00, &
    0.8017730000135972D+00, &
    0.8542770431031555D+00, &
    0.9110633416949405D+00, &
    0.9733286558781659D+00, &
    0.1042422171231561D+01, &
    0.1119872813134447D+01, &
    0.1207423594952871D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    bi = 0.0D+00
  else
    x = x_vec(n_data)
    bi = bi_vec(n_data)
  end if

  return
end
subroutine airy_bi_int_values ( n_data, x, fx )

!*****************************************************************************80
!
!! AIRY_BI_INT_VALUES returns some values of the integral of the Airy function.
!
!  Discussion:
!
!    The function is defined by:
!
!      AIRY_BI_INT(x) = integral ( 0 <= t <= x ) Bi(t) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.17660819031554631869D-01, &
    -0.15040424806140020451D-01, &
     0.14756446293227661920D-01, &
    -0.11847304264848446271D+00, &
    -0.64916741266165856037D-01, &
     0.97260832464381044540D-01, &
     0.50760058495287539119D-01, &
    -0.37300500963429492179D+00, &
    -0.13962988442666578531D+00, &
    -0.12001735266723296160D-02, &
     0.12018836117890354598D-02, &
     0.36533846550952011043D+00, &
     0.87276911673800812196D+00, &
     0.48219475263803429675D+02, &
     0.44006525804904178439D+06, &
     0.17608153976228301458D+07, &
     0.73779211705220007228D+07, &
     0.14780980310740671617D+09, &
     0.97037614223613433849D+11, &
     0.11632737638809878460D+15 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -12.0000000000D+00, &
    -10.0000000000D+00, &
     -8.0000000000D+00, &
     -7.5000000000D+00, &
     -7.0000000000D+00, &
     -6.5000000000D+00, &
     -4.0000000000D+00, &
     -1.0000000000D+00, &
     -0.2500000000D+00, &
     -0.0019531250D+00, &
      0.0019531250D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      4.0000000000D+00, &
      8.0000000000D+00, &
      8.5000000000D+00, &
      9.0000000000D+00, &
     10.0000000000D+00, &
     12.0000000000D+00, &
     14.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine airy_bi_prime_values ( n_data, x, bip )

!*****************************************************************************80
!
!! AIRY_BI_PRIME_VALUES returns some values of the Airy function Bi'(x).
!
!  Discussion:
!
!    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
!    solutions of the differential equation:
!
!      W'' - X * W = 0
!
!    In Mathematica, the function can be evaluated by:
!
!      AiryBiPrime[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BIP, the derivative of the Airy BI function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) bip
  real ( kind = 8 ), save, dimension ( n_max ) :: bip_vec = (/ &
    0.4482883573538264D+00, &
    0.4515126311496465D+00, &
    0.4617892843621509D+00, &
    0.4800490287524480D+00, &
    0.5072816760506224D+00, &
    0.5445725641405923D+00, &
    0.5931444786342857D+00, &
    0.6544059191721400D+00, &
    0.7300069016152518D+00, &
    0.8219038903072090D+00, &
    0.9324359333927756D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    bip = 0.0D+00
  else
    x = x_vec(n_data)
    bip = bip_vec(n_data)
  end if

  return
end
subroutine airy_cai_values ( n_data, x, cai )

!*****************************************************************************80
!
!! AIRY_CAI_VALUES returns some values of the Airy Ai(x) with complex argument.
!
!  Discussion:
!
!    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
!    solutions of the differential equation:
!
!      W'' - X * W = 0
!
!    In Mathematica, the function can be evaluated by:
!
!      AiryAi[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, complex ( kind = 8 ) X, the argument of the function.
!
!    Output, complex ( kind = 8 ) CAI, the value of the Airy AI function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  complex ( kind = 8 ) cai
  complex ( kind = 8 ), save, dimension ( n_max ) :: cai_vec = (/ &
    (0.1352924163128814D+00,   0.0000000000000000D+00),  &
    (0.1433824486882056D+00,  -0.1092193342707378D+00),  &
    (0.2215404472324631D+00,  -0.2588711788891803D+00),  &
    (0.4763929771766866D+00,  -0.3036484220291284D+00),  &
    (0.5983692170633874D+00,  -0.08154602160771214D+00), &
    (0.5355608832923521D+00,   0.00000000000000000D+00), &
    (0.5983692170633874D+00,   0.08154602160771214D+00), &
    (0.4763929771766866D+00,   0.3036484220291284D+00),  &
    (0.2215404472324631D+00,   0.2588711788891803D+00),  &
    (0.1433824486882056D+00,   0.1092193342707378D+00) /)
  integer ( kind = 4 ) n_data
  complex ( kind = 8 ) x
  complex ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    ( 1.000000000000000D+00,    0.0000000000000000D+00), &
    ( 0.8090169943749474D+00,   0.5877852522924731D+00), &
    ( 0.3090169943749474D+00,   0.9510565162951536D+00), &
    (-0.3090169943749474D+00,   0.9510565162951536D+00), &
    (-0.8090169943749474D+00,   0.5877852522924731D+00), &
    (-1.0000000000000000D+00,   0.0000000000000000D+00), &
    (-0.8090169943749474D+00,  -0.5877852522924731D+00), &
    (-0.3090169943749474D+00,  -0.9510565162951536D+00), &
    ( 0.3090169943749474D+00,  -0.9510565162951536D+00), &
    ( 0.8090169943749474D+00,  -0.5877852522924731D+00) /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = ( 0.0D+00, 0.0D+00 )
    cai = ( 0.0D+00, 0.0D+00 )
  else
    x = x_vec(n_data)
    cai = cai_vec(n_data)
  end if

  return
end
subroutine airy_cbi_values ( n_data, x, cbi )

!*****************************************************************************80
!
!! AIRY_CBI_VALUES returns some values of the Airy Bi(x) with complex argument.
!
!  Discussion:
!
!    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
!    solutions of the differential equation:
!
!      W'' - X * W = 0
!
!    In Mathematica, the function can be evaluated by:
!
!      AiryBi[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, complex ( kind = 8 ) X, the argument of the function.
!
!    Output, complex ( kind = 8 ) CBI, the value of the Airy BI function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  complex ( kind = 8 ) cbi
  complex ( kind = 8 ), save, dimension ( n_max ) :: cbi_vec = (/ &
    ( 1.207423594952871D+00,    0.0000000000000000D+00), &
    ( 0.9127160108293936D+00,   0.3800456133135556D+00), &
    ( 0.6824453575635721D+00,   0.3343047153635002D+00), &
    ( 0.5726265660086474D+00,   0.3988641086982559D+00), &
    ( 0.2511841251049547D+00,   0.3401447690712719D+00), &
    ( 0.1039973894969446D+00,   0.0000000000000000D+00), &
    ( 0.2511841251049547D+00,  -0.3401447690712719D+00), &
    ( 0.5726265660086474D+00,  -0.3988641086982559D+00), &
    ( 0.6824453575635721D+00,  -0.3343047153635002D+00), &
    ( 0.9127160108293936D+00,  -0.3800456133135556D+00) /)
  integer ( kind = 4 ) n_data
  complex ( kind = 8 ) x
  complex ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    ( 1.000000000000000D+00,    0.0000000000000000D+00), &
    ( 0.8090169943749474D+00,   0.5877852522924731D+00), &
    ( 0.3090169943749474D+00,   0.9510565162951536D+00), &
    (-0.3090169943749474D+00,   0.9510565162951536D+00), &
    (-0.8090169943749474D+00,   0.5877852522924731D+00), &
    (-1.0000000000000000D+00,   0.0000000000000000D+00), &
    (-0.8090169943749474D+00,  -0.5877852522924731D+00), &
    (-0.3090169943749474D+00,  -0.9510565162951536D+00), &
    ( 0.3090169943749474D+00,  -0.9510565162951536D+00), &
    ( 0.8090169943749474D+00,  -0.5877852522924731D+00) /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = ( 0.0D+00, 0.0D+00 )
    cbi = ( 0.0D+00, 0.0D+00 )
  else
    x = x_vec(n_data)
    cbi = cbi_vec(n_data)
  end if

  return
end
subroutine airy_gi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! AIRY_GI_VALUES returns some values of the Airy Gi function.
!
!  Discussion:
!
!    The function is defined by:
!
!      AIRY_GI(x) = integral ( 0 <= t < +oo ) sin ( x*t+t^3/3) dt / pi
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.20468308070040542435D+00, &
     0.18374662832557904078D+00, &
    -0.11667221729601528265D+00, &
     0.31466934902729557596D+00, &
    -0.37089040722426257729D+00, &
    -0.25293059772424019694D+00, &
     0.28967410658692701936D+00, &
    -0.34644836492634090590D+00, &
     0.28076035913873049496D+00, &
     0.21814994508094865815D+00, &
     0.20526679000810503329D+00, &
     0.22123695363784773258D+00, &
     0.23521843981043793760D+00, &
     0.82834303363768729338D-01, &
     0.45757385490989281893D-01, &
     0.44150012014605159922D-01, &
     0.39951133719508907541D-01, &
     0.35467706833949671483D-01, &
     0.31896005100679587981D-01, &
     0.26556892713512410405D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     -0.0019531250D+00, &
     -0.1250000000D+00, &
     -1.0000000000D+00, &
     -4.0000000000D+00, &
     -8.0000000000D+00, &
     -8.2500000000D+00, &
     -9.0000000000D+00, &
    -10.0000000000D+00, &
    -11.0000000000D+00, &
    -13.0000000000D+00, &
      0.0019531250D+00, &
      0.1250000000D+00, &
      1.0000000000D+00, &
      4.0000000000D+00, &
      7.0000000000D+00, &
      7.2500000000D+00, &
      8.0000000000D+00, &
      9.0000000000D+00, &
     10.0000000000D+00, &
     12.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine airy_hi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! AIRY_HI_VALUES returns some values of the Airy Hi function.
!
!  Discussion:
!
!    The function is defined by:
!
!      AIRY_HI(x) = integral ( 0 <= t < +oo ) exp(x*t-t^3/3) dt / pi
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.40936798278458884024D+00, &
    0.37495291608048868619D+00, &
    0.22066960679295989454D+00, &
    0.77565356679703713590D-01, &
    0.39638826473124717315D-01, &
    0.38450072575004151871D-01, &
    0.35273216868317898556D-01, &
    0.31768535282502272742D-01, &
    0.28894408288051391369D-01, &
    0.24463284011678541180D-01, &
    0.41053540139998941517D+00, &
    0.44993502381204990817D+00, &
    0.97220515514243332184D+00, &
    0.83764237105104371193D+02, &
    0.80327744952044756016D+05, &
    0.15514138847749108298D+06, &
    0.11995859641733262114D+07, &
    0.21472868855967642259D+08, &
    0.45564115351632913590D+09, &
    0.32980722582904761929D+12 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     -0.0019531250D+00, &
     -0.1250000000D+00, &
     -1.0000000000D+00, &
     -4.0000000000D+00, &
     -8.0000000000D+00, &
     -8.2500000000D+00, &
     -9.0000000000D+00, &
    -10.0000000000D+00, &
    -11.0000000000D+00, &
    -13.0000000000D+00, &
      0.0019531250D+00, &
      0.1250000000D+00, &
      1.0000000000D+00, &
      4.0000000000D+00, &
      7.0000000000D+00, &
      7.2500000000D+00, &
      8.0000000000D+00, &
      9.0000000000D+00, &
     10.0000000000D+00, &
     12.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine arccos_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ARCCOS_VALUES returns some values of the arc cosine function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ArcCos[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.6709637479564564156D+00, &
    1.5707963267948966192D+00, &
    1.4706289056333368229D+00, &
    1.3694384060045658278D+00, &
    1.2661036727794991113D+00, &
    1.1592794807274085998D+00, &
    1.0471975511965977462D+00, &
    0.92729521800161223243D+00, &
    0.79539883018414355549D+00, &
    0.64350110879328438680D+00, &
    0.45102681179626243254D+00, &
    0.00000000000000000000D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.1D+00, &
     0.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.3D+00, &
     0.4D+00, &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine arccosh_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ARCCOSH_VALUES returns some values of the hyperbolic arc cosine function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ArcCosh[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000000D+00, &
    0.14130376948564857735D+00, &
    0.44356825438511518913D+00, &
    0.62236250371477866781D+00, &
    0.75643291085695958624D+00, &
    0.86701472649056510395D+00, &
    0.96242365011920689500D+00, &
    1.3169578969248167086D+00, &
    1.7627471740390860505D+00, &
    1.8115262724608531070D+00, &
    2.0634370688955605467D+00, &
    2.2924316695611776878D+00, &
    2.9932228461263808979D+00, &
    5.2982923656104845907D+00, &
    7.6009022095419886114D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
       1.0D+00, &
       1.01D+00, &
       1.1D+00, &
       1.2D+00, &
       1.3D+00, &
       1.4D+00, &
       1.5D+00, &
       2.0D+00, &
       3.0D+00, &
       3.1415926535897932385D+00, &
       4.0D+00, &
       5.0D+00, &
      10.0D+00, &
     100.0D+00, &
    1000.0D+00  /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine arcsin_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ARCSIN_VALUES returns some values of the arc sine function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ArcSin[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.10016742116155979635D+00, &
     0.00000000000000000000D+00, &
     0.10016742116155979635D+00, &
     0.20135792079033079146D+00, &
     0.30469265401539750797D+00, &
     0.41151684606748801938D+00, &
     0.52359877559829887308D+00, &
     0.64350110879328438680D+00, &
     0.77539749661075306374D+00, &
     0.92729521800161223243D+00, &
     1.1197695149986341867D+00, &
     1.5707963267948966192D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.1D+00, &
     0.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.3D+00, &
     0.4D+00, &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine arcsinh_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ARCSINH_VALUES returns some values of the hyperbolic arc sine function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ArcSinh[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -2.3124383412727526203D+00, &
    -0.88137358701954302523D+00, &
     0.00000000000000000000D+00, &
     0.099834078899207563327D+00, &
     0.19869011034924140647D+00, &
     0.29567304756342243910D+00, &
     0.39003531977071527608D+00, &
     0.48121182505960344750D+00, &
     0.56882489873224753010D+00, &
     0.65266656608235578681D+00, &
     0.73266825604541086415D+00, &
     0.80886693565278246251D+00, &
     0.88137358701954302523D+00, &
     1.4436354751788103425D+00, &
     1.8184464592320668235D+00, &
     2.0947125472611012942D+00, &
     2.3124383412727526203D+00, &
     2.9982229502979697388D+00, &
     5.2983423656105887574D+00, &
     7.6009027095419886115D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
       -5.0D+00, &
       -1.0D+00, &
        0.0D+00, &
        0.1D+00, &
        0.2D+00, &
        0.3D+00, &
        0.4D+00, &
        0.5D+00, &
        0.6D+00, &
        0.7D+00, &
        0.8D+00, &
        0.9D+00, &
        1.0D+00, &
        2.0D+00, &
        3.0D+00, &
        4.0D+00, &
        5.0D+00, &
       10.0D+00, &
      100.0D+00, &
     1000.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine arctan_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ARCTAN_VALUES returns some values of the arc tangent function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ArcTan[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.00000000000000000000D+00, &
    0.24497866312686415417D+00, &
    0.32175055439664219340D+00, &
    0.46364760900080611621D+00, &
    0.78539816339744830962D+00, &
    1.1071487177940905030D+00, &
    1.2490457723982544258D+00, &
    1.3258176636680324651D+00, &
    1.3734007669450158609D+00, &
    1.4711276743037345919D+00, &
    1.5208379310729538578D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00000000000000000000D+00, &
    0.25000000000000000000D+00, &
    0.33333333333333333333D+00, &
    0.50000000000000000000D+00, &
    1.0000000000000000000D+00, &
    2.0000000000000000000D+00, &
    3.0000000000000000000D+00, &
    4.0000000000000000000D+00, &
    5.0000000000000000000D+00, &
    10.000000000000000000D+00, &
    20.000000000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine arctan_int_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ARCTAN_INT_VALUES returns some values of the inverse tangent integral.
!
!  Discussion:
!
!    The function is defined by:
!
!      ARCTAN_INT(x) = integral ( 0 <= t <= x ) arctan ( t ) / t dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.19531241721588483191D-02, &
    -0.39062433772980711281D-02, &
     0.78124470192576499535D-02, &
     0.15624576181996527280D-01, &
    -0.31246610349485401551D-01, &
     0.62472911335014397321D-01, &
     0.12478419717389654039D+00, &
    -0.24830175098230686908D+00, &
     0.48722235829452235711D+00, &
     0.91596559417721901505D+00, &
     0.12749694484943800618D+01, &
    -0.15760154034463234224D+01, &
     0.24258878412859089996D+01, &
     0.33911633326292997361D+01, &
     0.44176450919422186583D+01, &
    -0.47556713749547247774D+01, &
     0.50961912150934111303D+01, &
     0.53759175735714876256D+01, &
    -0.61649904785027487422D+01, &
     0.72437843013083534973D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
     -0.0039062500D+00, &
      0.0078125000D+00, &
      0.0156250000D+00, &
     -0.0312500000D+00, &
      0.0625000000D+00, &
      0.1250000000D+00, &
     -0.2500000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
     -2.0000000000D+00, &
      4.0000000000D+00, &
      8.0000000000D+00, &
     16.0000000000D+00, &
    -20.0000000000D+00, &
     25.0000000000D+00, &
     30.0000000000D+00, &
    -50.0000000000D+00, &
    100.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine arctan2_values ( n_data, x, y, f )

!*****************************************************************************80
!
!! ARCTAN2_VALUES: arc tangent function of two arguments.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ArcTan[x,y]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, Y, the arguments of the function.
!
!    Output, real ( kind = 8 ) F, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 19

  real ( kind = 8 ) f
  real ( kind = 8 ), save, dimension ( n_max ) :: f_vec = (/ &
   -1.5707963267948966192D+00, &
   -1.0471975511965977462D+00, &
   -0.52359877559829887308D+00, &
    0.00000000000000000000D+00, &
    0.52359877559829887308D+00, &
    1.0471975511965977462D+00, &
    1.5707963267948966192D+00, &
    2.0943951023931954923D+00, &
    2.6179938779914943654D+00, &
    3.1415926535897932385D+00, &
   -2.6179938779914943654D+00, &
   -2.0943951023931954923D+00, &
   -1.5707963267948966192D+00, &
   -1.0471975511965977462D+00, &
   -0.52359877559829887308D+00, &
    0.00000000000000000000D+00, &
    0.52359877559829887308D+00, &
    1.0471975511965977462D+00, &
    1.5707963267948966192D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.86602540378443864676D+00, &
    1.00000000000000000000D+00, &
    0.86602540378443864676D+00, &
    0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
   -0.50000000000000000000D+00, &
   -0.86602540378443864676D+00, &
   -1.00000000000000000000D+00, &
   -0.86602540378443864676D+00, &
   -0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.86602540378443864676D+00, &
    1.00000000000000000000D+00, &
    0.86602540378443864676D+00, &
    0.50000000000000000000D+00, &
    0.00000000000000000000D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( n_max ) :: y_vec = (/ &
   -1.00000000000000000000D+00, &
   -0.86602540378443864676D+00, &
   -0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.86602540378443864676D+00, &
    1.00000000000000000000D+00, &
    0.86602540378443864676D+00, &
    0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
   -0.50000000000000000000D+00, &
   -0.86602540378443864676D+00, &
   -1.00000000000000000000D+00, &
   -0.86602540378443864676D+00, &
   -0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.86602540378443864676D+00, &
    1.00000000000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    y = 0.0D+00
    f = 0.0D+00
  else
    x = x_vec(n_data)
    y = y_vec(n_data)
    f = f_vec(n_data)
  end if

  return
end
subroutine arctanh_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ARCTANH_VALUES returns some values of the hyperbolic arc tangent function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ArcTanh[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.54930614433405484570D+00, &
     0.00000000000000000000D+00, &
     0.0010000003333335333335D+00, &
     0.10033534773107558064D+00, &
     0.20273255405408219099D+00, &
     0.30951960420311171547D+00, &
     0.42364893019360180686D+00, &
     0.54930614433405484570D+00, &
     0.69314718055994530942D+00, &
     0.86730052769405319443D+00, &
     1.0986122886681096914D+00, &
     1.4722194895832202300D+00, &
     2.6466524123622461977D+00, &
     3.8002011672502000318D+00, &
     7.2543286192620472067D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.5D+00, &
     0.0D+00, &
     0.001D+00, &
     0.1D+00, &
     0.2D+00, &
     0.3D+00, &
     0.4D+00, &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     0.99D+00, &
     0.999D+00, &
     0.999999D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bei0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BEI0_VALUES returns some values of the Kelvin BEI function of order NU = 0.
!
!  Discussion:
!
!    The function is defined by:
!
!      BER(NU,X) + i * BEI(NU,X) = exp(NU*Pi*I) * J(NU,X*exp(-PI*I/4))
!
!    where J(NU,X) is the J Bessel function.
!
!    In Mathematica, BEI(NU,X) can be defined by:
!
!      Im [ Exp [ NU * Pi * I ] * BesselJ [ NU, X * Exp[ -Pi * I / 4 ] ] ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.06249321838219946D+00, &
    0.2495660400366597D+00, &
    0.5575600623030867D+00, &
    0.9722916273066612D+00, &
    1.457182044159804D+00, &
   1.937586785266043D+00, &
   2.283249966853915D+00, &
   2.292690322699300D+00, &
   1.686017203632139D+00, &
   0.1160343815502004D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    2.5D+00, &
    3.0D+00, &
    3.5D+00, &
    4.0D+00, &
    4.5D+00, &
    5.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bei1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BEI1_VALUES returns some values of the Kelvin BEI function of order NU = 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      BER(NU,X) + i * BEI(NU,X) = exp(NU*Pi*I) * J(NU,X*exp(-PI*I/4))
!
!    where J(NU,X) is the J Bessel function.
!
!    In Mathematica, BEI(NU,X) can be defined by:
!
!      Im [ Exp [ NU * Pi * I ] * BesselJ [ NU, X * Exp[ -Pi * I / 4 ] ] ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.1711951797170153D+00, &
    0.3075566313755366D+00, &
    0.3678649890020899D+00, &
    0.2997754370020335D+00, &
    0.03866844396595048D+00, &
   -0.4874541770160708D+00, &
   -1.344042373111174D+00, &
   -2.563821688561078D+00, &
   -4.105685408400878D+00, &
   -5.797907901792625D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    2.5D+00, &
    3.0D+00, &
    3.5D+00, &
    4.0D+00, &
    4.5D+00, &
    5.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bell_values ( n_data, n, c )

!*****************************************************************************80
!
!! BELL_VALUES returns some values of the Bell numbers.
!
!  Discussion:
!
!    The Bell number B(N) is the number of restricted growth functions on N.
!
!    Note that the Stirling numbers of the second kind, S^m_n, count the
!    number of partitions of N objects into M classes, and so it is
!    true that
!
!      B(N) = S^1_N + S^2_N + ... + S^N_N.
!
!    The Bell numbers were named for Eric Temple Bell.
!
!    In Mathematica, the function can be evaluated by
!
!      Sum[StirlingS2[n,m],{m,1,n}]
!
!  Definition:
!
!    The Bell number B(N) is defined as the number of partitions (of
!    any size) of a set of N distinguishable objects.
!
!    A partition of a set is a division of the objects of the set into
!    subsets.
!
!  Examples:
!
!    There are 15 partitions of a set of 4 objects:
!
!      (1234),
!      (123) (4),
!      (124) (3),
!      (12) (34),
!      (12) (3) (4),
!      (134) (2),
!      (13) (24),
!      (13) (2) (4),
!      (14) (23),
!      (1) (234),
!      (1) (23) (4),
!      (14) (2) (3),
!      (1) (24) (3),
!      (1) (2) (34),
!      (1) (2) (3) (4).
!
!    and so B(4) = 15.
!
!  First values:
!
!     N         B(N)
!     0           1
!     1           1
!     2           2
!     3           5
!     4          15
!     5          52
!     6         203
!     7         877
!     8        4140
!     9       21147
!    10      115975
!
!  Recursion:
!
!    B(I) = sum ( 1 <= J <=I ) Binomial ( I-1, J-1 ) * B(I-J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the Bell number.
!
!    Output, integer ( kind = 4 ) C, the value of the Bell number.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
    1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine ber0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BER0_VALUES returns some values of the Kelvin BER function of order NU = 0.
!
!  Discussion:
!
!    The function is defined by:
!
!      BER(NU,X) + i * BEI(NU,X) = exp(NU*Pi*I) * J(NU,X*exp(-PI*I/4))
!
!    where J(NU,X) is the J Bessel function.
!
!    In Mathematica, BER(NU,X) can be defined by:
!
!      Re [ Exp [ NU * Pi * I ] * BesselJ [ NU, X * Exp[ -Pi * I / 4 ] ] ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.0000000000000000D+00, &
    0.9990234639908383D+00, &
    0.9843817812130869D+00, &
    0.9210721835462558D+00, &
    0.7517341827138082D+00, &
    0.3999684171295313D+00, &
   -0.2213802495986939D+00, &
   -1.193598179589928D+00, &
   -2.563416557258580D+00, &
   -4.299086551599756D+00, &
   -6.230082478666358D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    2.5D+00, &
    3.0D+00, &
    3.5D+00, &
    4.0D+00, &
    4.5D+00, &
    5.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine ber1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BER1_VALUES returns some values of the Kelvin BER function of order NU = 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      BER(NU,X) + i * BEI(NU,X) = exp(NU*Pi*I) * J(NU,X*exp(-PI*I/4))
!
!    where J(NU,X) is the J Bessel function.
!
!    In Mathematica, BER(NU,X) can be defined by:
!
!      Re [ Exp [ NU * Pi * I ] * BesselJ [ NU, X * Exp[ -Pi * I / 4 ] ] ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
    -0.1822431237551121D+00, &
    -0.3958682610197114D+00, &
    -0.6648654179597691D+00, &
    -0.9970776519264285D+00, &
    -1.373096897645111D+00, &
   -1.732644221128481D+00, &
   -1.959644131289749D+00, &
   -1.869248459031899D+00, &
   -1.202821631480086D+00, &
    0.3597766667766728D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    2.5D+00, &
    3.0D+00, &
    3.5D+00, &
    4.0D+00, &
    4.5D+00, &
    5.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bernoulli_number_values ( n_data, n, c )

!*****************************************************************************80
!
!! BERNOULLI_NUMBER_VALUES returns some values of the Bernoulli numbers.
!
!  Discussion:
!
!    The Bernoulli numbers are rational.
!
!    If we define the sum of the M-th powers of the first N integers as:
!
!      SIGMA(M,N) = sum ( 0 <= I <= N ) I^M
!
!    and let C(I,J) be the combinatorial coefficient:
!
!      C(I,J) = I! / ( ( I - J )! * J! )
!
!    then the Bernoulli numbers B(J) satisfy:
!
!      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)^(M+1-J)
!
!    In Mathematica, the function can be evaluated by:
!
!      BernoulliB[n]
!
!  First values:
!
!   B0  1                   =         1.00000000000
!   B1 -1/2                 =        -0.50000000000
!   B2  1/6                 =         1.66666666666
!   B3  0                   =         0
!   B4 -1/30                =        -0.03333333333
!   B5  0                   =         0
!   B6  1/42                =         0.02380952380
!   B7  0                   =         0
!   B8 -1/30                =        -0.03333333333
!   B9  0                   =         0
!  B10  5/66                =         0.07575757575
!  B11  0                   =         0
!  B12 -691/2730            =        -0.25311355311
!  B13  0                   =         0
!  B14  7/6                 =         1.16666666666
!  B15  0                   =         0
!  B16 -3617/510            =        -7.09215686274
!  B17  0                   =         0
!  B18  43867/798           =        54.97117794486
!  B19  0                   =         0
!  B20 -174611/330          =      -529.12424242424
!  B21  0                   =         0
!  B22  854,513/138         =      6192.123
!  B23  0                   =         0
!  B24 -236364091/2730      =    -86580.257
!  B25  0                   =         0
!  B26  8553103/6           =   1425517.16666
!  B27  0                   =         0
!  B28 -23749461029/870     = -27298231.0678
!  B29  0                   =         0
!  B30  8615841276005/14322 = 601580873.901
!
!  Recursion:
!
!    With C(N+1,K) denoting the standard binomial coefficient,
!
!    B(0) = 1.0
!    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
!
!  Special Values:
!
!    Except for B(1), all Bernoulli numbers of odd index are 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the Bernoulli number.
!
!    Output, real ( kind = 8 ) C, the value of the Bernoulli number.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  real ( kind = 8 ) c
  real ( kind = 8 ), save, dimension ( n_max ) :: c_vec = (/ &
     0.1000000000000000D+01, &
    -0.5000000000000000D+00, &
     0.1666666666666667D+00, &
     0.0000000000000000D+00, &
    -0.3333333333333333D-01, &
    -0.2380952380952380D-01, &
    -0.3333333333333333D-01, &
     0.7575757575757575D-01, &
    -0.5291242424242424D+03, &
     0.6015808739006424D+09 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 6,  8, 10, 20, 30 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0.0D+00
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine bernoulli_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! BERNOULLI_POLY_VALUES returns some values of the Bernoulli polynomials.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BernoulliB[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the Bernoulli polynomial.
!
!    Output, real ( kind = 8 ) X, the argument of the Bernoulli polynomial.
!
!    Output, real ( kind = 8 ) FX, the value of the Bernoulli polynomial.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 27

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
    -0.3000000000000000D+00, &
     0.6666666666666667D-02, &
     0.4800000000000000D-01, &
    -0.7733333333333333D-02, &
    -0.2368000000000000D-01, &
     0.6913523809523810D-02, &
     0.2490880000000000D-01, &
    -0.1014997333333333D-01, &
    -0.4527820800000000D-01, &
     0.2332631815757576D-01, &
    -0.3125000000000000D+00, &
    -0.1142400000000000D+00, &
    -0.0176800000000000D+00, &
     0.0156800000000000D+00, &
     0.0147400000000000D+00, &
     0.0000000000000000D+00, &
    -0.1524000000000000D-01, &
    -0.2368000000000000D-01, &
    -0.2282000000000000D-01, &
    -0.1376000000000000D-01, &
     0.0000000000000000D+01, &
     0.1376000000000000D-01, &
     0.2282000000000000D-01, &
     0.2368000000000000D-01, &
     0.1524000000000000D-01, &
     0.0000000000000000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, &
     1, &
     2, &
     3, &
     4, &
     5, &
     6, &
     7, &
     8, &
     9, &
    10, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
    -0.5D+00, &
    -0.4D+00, &
    -0.3D+00, &
    -0.2D+00, &
    -0.1D+00, &
     0.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.3D+00, &
     0.4D+00, &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bernstein_poly_values ( n_data, n, k, x, b )

!*****************************************************************************80
!
!! BERNSTEIN_POLY_VALUES returns some values of the Bernstein polynomials.
!
!  Discussion:
!
!    The Bernstein polynomials are assumed to be based on [0,1].
!
!    The formula for the Bernstein polynomials is
!
!      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
!
!    In Mathematica, the function can be evaluated by:
!
!      Binomial[n,i] * (1-x)^(n-i) * x^i
!
!  First values:
!
!    B(0,0)(X) = 1
!
!    B(1,0)(X) =      1-X
!    B(1,1)(X) =               X
!
!    B(2,0)(X) =     (1-X)^2
!    B(2,1)(X) = 2 * (1-X)   * X
!    B(2,2)(X) =               X^2
!
!    B(3,0)(X) =     (1-X)^3
!    B(3,1)(X) = 3 * (1-X)^2 * X
!    B(3,2)(X) = 3 * (1-X)   * X^2
!    B(3,3)(X) =               X^3
!
!    B(4,0)(X) =     (1-X)^4
!    B(4,1)(X) = 4 * (1-X)^3 * X
!    B(4,2)(X) = 6 * (1-X)^2 * X^2
!    B(4,3)(X) = 4 * (1-X)   * X^3
!    B(4,4)(X) =               X^4
!
!  Special values:
!
!    B(N,I)(X) has a unique maximum value at X = I/N.
!
!    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
!
!    B(N,I)(1/2) = C(N,K) / 2^N
!
!    For a fixed X and N, the polynomials add up to 1:
!
!      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Output, integer ( kind = 4 ) K, the index of the polynomial.
!
!    Output, real ( kind = 8 ) X, the argument of the polynomial.
!
!    Output, real ( kind = 8 ) B, the value of the polynomial B(N,K)(X).
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    0.1000000000000000D+01, &
    0.7500000000000000D+00, &
    0.2500000000000000D+00, &
    0.5625000000000000D+00, &
    0.3750000000000000D+00, &
    0.6250000000000000D-01, &
    0.4218750000000000D+00, &
    0.4218750000000000D+00, &
    0.1406250000000000D+00, &
    0.1562500000000000D-01, &
    0.3164062500000000D+00, &
    0.4218750000000000D+00, &
    0.2109375000000000D+00, &
    0.4687500000000000D-01, &
    0.3906250000000000D-02 /)
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save, dimension ( n_max ) :: k_vec = (/ &
    0, &
    0, 1, &
    0, 1, 2, &
    0, 1, 2, 3, &
    0, 1, 2, 3, 4 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0, &
    1, 1, &
    2, 2, 2, &
    3, 3, 3, 3, &
    4, 4, 4, 4, 4 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    k = 0
    x = 0.0D+00
    b = 0.0D+00
  else
    n = n_vec(n_data)
    k = k_vec(n_data)
    x = x_vec(n_data)
    b = b_vec(n_data)
  end if

  return
end
subroutine bessel_i0_int_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I0_INT_VALUES returns some values of the Bessel I0 integral.
!
!  Discussion:
!
!    The function is defined by:
!
!      I0_INT(x) = integral ( 0 <= t <= x ) I0(t) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.19531256208818052282D-02, &
    -0.39062549670565734544D-02, &
     0.62520348032546565850D-01, &
     0.12516285581366971819D+00, &
    -0.51051480879740303760D+00, &
     0.10865210970235898158D+01, &
     0.27750019054282535299D+01, &
    -0.13775208868039716639D+02, &
     0.46424372058106108576D+03, &
     0.64111867658021584522D+07, &
    -0.10414860803175857953D+08, &
     0.44758598913855743089D+08, &
    -0.11852985311558287888D+09, &
     0.31430078220715992752D+09, &
    -0.83440212900794309620D+09, &
     0.22175367579074298261D+10, &
     0.58991731842803636487D+10, &
    -0.41857073244691522147D+11, &
     0.79553885818472357663D+12, &
     0.15089715082719201025D+17 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
     -0.0039062500D+00, &
      0.0625000000D+00, &
      0.1250000000D+00, &
     -0.5000000000D+00, &
      1.0000000000D+00, &
      2.0000000000D+00, &
     -4.0000000000D+00, &
      8.0000000000D+00, &
     18.0000000000D+00, &
    -18.5000000000D+00, &
     20.0000000000D+00, &
    -21.0000000000D+00, &
     22.0000000000D+00, &
    -23.0000000000D+00, &
     24.0000000000D+00, &
     25.0000000000D+00, &
    -27.0000000000D+00, &
     30.0000000000D+00, &
     40.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_i0_spherical_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I0_SPHERICAL_VALUES: some values of the Spherical Bessel function i0.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi/(2*x)] * BesselI[1/2,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.001667500198440D+00, &
    1.006680012705470D+00, &
    1.026880814507039D+00, &
    1.061089303580402D+00, &
    1.110132477734529D+00, &
    1.175201193643801D+00, &
    1.257884462843477D+00, &
    1.360215358179667D+00, &
    1.484729970750144D+00, &
    1.634541271164267D+00, &
    1.813430203923509D+00, &
    2.025956895698133D+00, &
    2.277595505698373D+00, &
    2.574897010920645D+00, &
    2.925685126512827D+00, &
    3.339291642469967D+00, &
    3.826838748926716D+00, &
    4.401577467270101D+00, &
    5.079293155726485D+00, &
    5.878791279137455D+00, &
    6.822479299281938D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00, &
    3.2D+00, &
    3.4D+00, &
    3.6D+00, &
    3.8D+00, &
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_i0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I0_VALUES returns some values of the I0 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    The modified Bessel function I0(Z) corresponds to N = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1000000000000000D+01, &
    0.1010025027795146D+01, &
    0.1040401782229341D+01, &
    0.1092045364317340D+01, &
    0.1166514922869803D+01, &
    0.1266065877752008D+01, &
    0.1393725584134064D+01, &
    0.1553395099731217D+01, &
    0.1749980639738909D+01, &
    0.1989559356618051D+01, &
    0.2279585302336067D+01, &
    0.3289839144050123D+01, &
    0.4880792585865024D+01, &
    0.7378203432225480D+01, &
    0.1130192195213633D+02, &
    0.1748117185560928D+02, &
    0.2723987182360445D+02, &
    0.6723440697647798D+02, &
    0.4275641157218048D+03, &
    0.2815716628466254D+04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.20D+00, &
    0.40D+00, &
    0.60D+00, &
    0.80D+00, &
    0.10D+01, &
    0.12D+01, &
    0.14D+01, &
    0.16D+01, &
    0.18D+01, &
    0.20D+01, &
    0.25D+01, &
    0.30D+01, &
    0.35D+01, &
    0.40D+01, &
    0.45D+01, &
    0.50D+01, &
    0.60D+01, &
    0.80D+01, &
    0.10D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_i1_spherical_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I1_SPHERICAL_VALUES: some values of the Spherical Bessel function i1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi/(2*x)] * BesselI[3/2,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.03336667857363341D+00, &
    0.06693371456802954D+00, &
    0.1354788933285401D+00, &
    0.2072931911031093D+00, &
    0.2841280857128948D+00, &
    0.3678794411714423D+00, &
    0.4606425870674146D+00, &
    0.5647736480096238D+00, &
    0.6829590627779635D+00, &
    0.8182955028627777D+00, &
    0.9743827435800610D+00, &
    1.155432469636406D+00, &
    1.366396525527973D+00, &
    1.613118767572064D+00, &
    1.902515460838681D+00, &
    2.242790117769266D+00, &
    2.643689828630357D+00, &
    3.116811526884873D+00, &
    3.675968313148932D+00, &
    4.337627987747642D+00, &
    5.121438384183637D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00, &
    3.2D+00, &
    3.4D+00, &
    3.6D+00, &
    3.8D+00, &
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_i1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I1_VALUES returns some values of the I1 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.1005008340281251D+00, &
    0.2040267557335706D+00, &
    0.3137040256049221D+00, &
    0.4328648026206398D+00, &
    0.5651591039924850D+00, &
    0.7146779415526431D+00, &
    0.8860919814143274D+00, &
    0.1084810635129880D+01, &
    0.1317167230391899D+01, &
    0.1590636854637329D+01, &
    0.2516716245288698D+01, &
    0.3953370217402609D+01, &
    0.6205834922258365D+01, &
    0.9759465153704450D+01, &
    0.1538922275373592D+02, &
    0.2433564214245053D+02, &
    0.6134193677764024D+02, &
    0.3998731367825601D+03, &
    0.2670988303701255D+04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.20D+00, &
    0.40D+00, &
    0.60D+00, &
    0.80D+00, &
    0.10D+01, &
    0.12D+01, &
    0.14D+01, &
    0.16D+01, &
    0.18D+01, &
    0.20D+01, &
    0.25D+01, &
    0.30D+01, &
    0.35D+01, &
    0.40D+01, &
    0.45D+01, &
    0.50D+01, &
    0.60D+01, &
    0.80D+01, &
    0.10D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_in_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_IN_VALUES returns some values of the In Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5016687513894678D-02, &
    0.1357476697670383D+00, &
    0.6889484476987382D+00, &
    0.1276466147819164D+01, &
    0.2245212440929951D+01, &
    0.1750561496662424D+02, &
    0.2281518967726004D+04, &
    0.3931278522104076D+08, &
    0.2216842492433190D-01, &
    0.2127399592398527D+00, &
    0.1033115016915114D+02, &
    0.1758380716610853D+04, &
    0.2677764138883941D+21, &
    0.2714631559569719D-03, &
    0.9825679323131702D-02, &
    0.2157974547322546D+01, &
    0.7771882864032600D+03, &
    0.2278548307911282D+21, &
    0.2752948039836874D-09, &
    0.3016963879350684D-06, &
    0.4580044419176051D-02, &
    0.2189170616372337D+02, &
    0.1071597159477637D+21, &
    0.3966835985819020D-24, &
    0.4310560576109548D-18, &
    0.5024239357971806D-10, &
    0.1250799735644948D-03, &
    0.5442008402752998D+19 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) nu
  integer ( kind = 4 ), save, dimension ( n_max ) :: nu_vec = (/ &
     2,  2,  2,  2, &
     2,  2,  2,  2, &
     3,  3,  3,  3, &
     3,  5,  5,  5, &
     5,  5, 10, 10, &
    10, 10, 10, 20, &
    20, 20, 20, 20 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_ix_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_IX_VALUES returns some values of the Ix Bessel function.
!
!  Discussion:
!
!    This set of data considers the less common case in which the
!    index of the Bessel function In is actually not an integer.
!    We may suggest this case by occasionally replacing the symbol
!    "In" by "Ix".
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.3592084175833614D+00,  &
    0.9376748882454876D+00,  &
    2.046236863089055D+00,   &
    3.053093538196718D+00,   &
    4.614822903407601D+00,   &
    26.47754749755907D+00,   &
    2778.784603874571D+00,   &
    4.327974627242893D+07,   &
    0.2935253263474798D+00,  &
    1.099473188633110D+00,   &
    21.18444226479414D+00,   &
    2500.906154942118D+00,   &
    2.866653715931464D+20,   &
    0.05709890920304825D+00, &
    0.3970270801393905D+00,  &
    13.76688213868258D+00,   &
    2028.512757391936D+00,   &
    2.753157630035402D+20,   &
    0.4139416015642352D+00,  &
    1.340196758982897D+00,   &
    22.85715510364670D+00,   &
    2593.006763432002D+00,   &
    2.886630075077766D+20,   &
    0.03590910483251082D+00, &
    0.2931108636266483D+00,  &
    11.99397010023068D+00,   &
    1894.575731562383D+00,   &
    2.716911375760483D+20 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 8 ), save, dimension ( n_max ) :: nu_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_j0_int_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_J0_INT_VALUES returns some values of the Bessel J0 integral.
!
!  Discussion:
!
!    The function is defined by:
!
!      J0_INT(x) = integral ( 0 <= t <= x ) J0(t) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.97656242238978822427D-03, &
     0.39062450329491108875D-02, &
    -0.62479657927917933620D-01, &
     0.12483733492120479139D+00, &
    -0.48968050664604505505D+00, &
     0.91973041008976023931D+00, &
    -0.14257702931970265690D+01, &
     0.10247341594606064818D+01, &
    -0.12107468348304501655D+01, &
     0.11008652032736190799D+01, &
    -0.10060334829904124192D+01, &
     0.81330572662485953519D+00, &
    -0.10583788214211277585D+01, &
     0.87101492116545875169D+00, &
    -0.88424908882547488420D+00, &
     0.11257761503599914603D+01, &
    -0.90141212258183461184D+00, &
     0.91441344369647797803D+00, &
    -0.94482281938334394886D+00, &
     0.92266255696016607257D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0009765625D+00, &
      0.0039062500D+00, &
     -0.0625000000D+00, &
      0.1250000000D+00, &
     -0.5000000000D+00, &
      1.0000000000D+00, &
     -2.0000000000D+00, &
      4.0000000000D+00, &
     -8.0000000000D+00, &
     16.0000000000D+00, &
    -16.5000000000D+00, &
     18.0000000000D+00, &
    -20.0000000000D+00, &
     25.0000000000D+00, &
    -30.0000000000D+00, &
     40.0000000000D+00, &
    -50.0000000000D+00, &
     75.0000000000D+00, &
    -80.0000000000D+00, &
    100.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_j0_spherical_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_J0_SPHERICAL_VALUES: some values of the Spherical Bessel function j0.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi/(2*x)] * BesselJ[1/2,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.9983341664682815D+00, &
     0.9933466539753061D+00, &
     0.9735458557716262D+00, &
     0.9410707889917256D+00, &
     0.8966951136244035D+00, &
     0.8414709848078965D+00, &
     0.7766992383060220D+00, &
     0.7038926642774716D+00, &
     0.6247335019009407D+00, &
     0.5410264615989973D+00, &
     0.4546487134128408D+00, &
     0.3674983653725410D+00, &
     0.2814429918963129D+00, &
     0.1982697583928709D+00, &
     0.1196386250556803D+00, &
     0.4704000268662241D-01, &
    -0.1824191982111872D-01, &
    -0.7515914765495039D-01, &
    -0.1229223453596812D+00, &
    -0.1610152344586103D+00, &
    -0.1892006238269821D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00, &
    3.2D+00, &
    3.4D+00, &
    3.6D+00, &
    3.8D+00, &
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_j0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_J0_VALUES returns some values of the J0 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1775967713143383D+00, &
    -0.3971498098638474D+00, &
    -0.2600519549019334D+00, &
     0.2238907791412357D+00, &
     0.7651976865579666D+00, &
     0.1000000000000000D+01, &
     0.7651976865579666D+00, &
     0.2238907791412357D+00, &
    -0.2600519549019334D+00, &
    -0.3971498098638474D+00, &
    -0.1775967713143383D+00, &
     0.1506452572509969D+00, &
     0.3000792705195556D+00, &
     0.1716508071375539D+00, &
    -0.9033361118287613D-01, &
    -0.2459357644513483D+00, &
    -0.1711903004071961D+00, &
     0.4768931079683354D-01, &
     0.2069261023770678D+00, &
     0.1710734761104587D+00, &
    -0.1422447282678077D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -5.0D+00, &
    -4.0D+00, &
    -3.0D+00, &
    -2.0D+00, &
    -1.0D+00, &
     0.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_j1_spherical_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_J1_SPHERICAL_VALUES: some values of the Spherical Bessel function j1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi/(2*x)] * BesselJ[3/2,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.3330001190255757D-01, &
    0.6640038067032223D-01, &
    0.1312121544218529D+00, &
    0.1928919568034122D+00, &
    0.2499855053465475D+00, &
    0.3011686789397568D+00, &
    0.3452845698577903D+00, &
    0.3813753724123076D+00, &
    0.4087081401263934D+00, &
    0.4267936423844913D+00, &
    0.4353977749799916D+00, &
    0.4345452193763121D+00, &
    0.4245152947656493D+00, &
    0.4058301968314685D+00, &
    0.3792360591872637D+00, &
    0.3456774997623560D+00, &
    0.3062665174917607D+00, &
    0.2622467779189737D+00, &
    0.2149544641595738D+00, &
    0.1657769677515280D+00, &
    0.1161107492591575D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00, &
    3.2D+00, &
    3.4D+00, &
    3.6D+00, &
    3.8D+00, &
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_j1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_J1_VALUES returns some values of the J1 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.3275791375914652D+00, &
     0.6604332802354914D-01, &
    -0.3390589585259365D+00, &
    -0.5767248077568734D+00, &
    -0.4400505857449335D+00, &
     0.0000000000000000D+00, &
     0.4400505857449335D+00, &
     0.5767248077568734D+00, &
     0.3390589585259365D+00, &
    -0.6604332802354914D-01, &
    -0.3275791375914652D+00, &
    -0.2766838581275656D+00, &
    -0.4682823482345833D-02, &
     0.2346363468539146D+00, &
     0.2453117865733253D+00, &
     0.4347274616886144D-01, &
    -0.1767852989567215D+00, &
    -0.2234471044906276D+00, &
    -0.7031805212177837D-01, &
     0.1333751546987933D+00, &
     0.2051040386135228D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -5.0D+00, &
    -4.0D+00, &
    -3.0D+00, &
    -2.0D+00, &
    -1.0D+00, &
     0.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_jn_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_JN_VALUES returns some values of the Jn Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1149034849319005D+00, &
     0.3528340286156377D+00, &
     0.4656511627775222D-01, &
     0.2546303136851206D+00, &
    -0.5971280079425882D-01, &
     0.2497577302112344D-03, &
     0.7039629755871685D-02, &
     0.2611405461201701D+00, &
    -0.2340615281867936D+00, &
    -0.8140024769656964D-01, &
     0.2630615123687453D-09, &
     0.2515386282716737D-06, &
     0.1467802647310474D-02, &
     0.2074861066333589D+00, &
    -0.1138478491494694D+00, &
     0.3873503008524658D-24, &
     0.3918972805090754D-18, &
     0.2770330052128942D-10, &
     0.1151336924781340D-04, &
    -0.1167043527595797D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) nu
  integer ( kind = 4 ), save, dimension ( n_max ) :: nu_vec = (/ &
     2,  2,  2,  2, &
     2,  5,  5,  5, &
     5,  5, 10, 10, &
    10, 10, 10, 20, &
    20, 20, 20, 20 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_jx_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_JX_VALUES returns some values of the Jx Bessel function.
!
!  Discussion:
!
!    This set of data considers the less common case in which the
!    index of the Bessel function Jn is actually not an integer.
!    We may suggest this case by occasionally replacing the symbol
!    "Jn" by "Jx".
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.3544507442114011D+00, &
        0.6713967071418031D+00, &
        0.5130161365618278D+00, &
        0.3020049060623657D+00, &
        0.06500818287737578D+00, &
       -0.3421679847981618D+00, &
       -0.1372637357550505D+00, &
        0.1628807638550299D+00, &
        0.2402978391234270D+00, &
        0.4912937786871623D+00, &
       -0.1696513061447408D+00, &
        0.1979824927558931D+00, &
       -0.1094768729883180D+00, &
        0.04949681022847794D+00, &
        0.2239245314689158D+00, &
        0.2403772011113174D+00, &
        0.1966584835818184D+00, &
        0.02303721950962553D+00, &
        0.3314145508558904D+00, &
        0.5461734240402840D+00, &
       -0.2616584152094124D+00, &
        0.1296035513791289D+00, &
       -0.1117432171933552D+00, &
        0.03142623570527935D+00, &
        0.1717922192746527D+00, &
        0.3126634069544786D+00, &
        0.1340289119304364D+00, &
        0.06235967135106445D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 8 ), save, dimension ( n_max ) :: nu_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_k0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_K0_VALUES returns some values of the K0 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    The modified Bessel function K0(Z) corresponds to N = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselK[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.2427069024702017D+01, &
    0.1752703855528146D+01, &
    0.1114529134524434D+01, &
    0.7775220919047293D+00, &
    0.5653471052658957D+00, &
    0.4210244382407083D+00, &
    0.3185082202865936D+00, &
    0.2436550611815419D+00, &
    0.1879547519693323D+00, &
    0.1459314004898280D+00, &
    0.1138938727495334D+00, &
    0.6234755320036619D-01, &
    0.3473950438627925D-01, &
    0.1959889717036849D-01, &
    0.1115967608585302D-01, &
    0.6399857243233975D-02, &
    0.3691098334042594D-02, &
    0.1243994328013123D-02, &
    0.1464707052228154D-03, &
    0.1778006231616765D-04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     0.2D+00, &
     0.4D+00, &
     0.6D+00, &
     0.8D+00, &
     1.0D+00, &
     1.2D+00, &
     1.4D+00, &
     1.6D+00, &
     1.8D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00, &
     4.5D+00, &
     5.0D+00, &
     6.0D+00, &
     8.0D+00, &
    10.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_k0_int_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_K0_INT_VALUES returns some values of the Bessel K0 integral.
!
!  Discussion:
!
!    The function is defined by:
!
!      K0_INT(x) = integral ( 0 <= t <= x ) K0(t) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.78587929563466784589D-02, &
    0.26019991617330578111D-01, &
    0.24311842237541167904D+00, &
    0.39999633750480508861D+00, &
    0.92710252093114907345D+00, &
    0.12425098486237782662D+01, &
    0.14736757343168286825D+01, &
    0.15606495706051741364D+01, &
    0.15673873907283660493D+01, &
    0.15696345532693743714D+01, &
    0.15701153443250786355D+01, &
    0.15706574852894436220D+01, &
    0.15707793116159788598D+01, &
    0.15707942066465767196D+01, &
    0.15707962315469192247D+01, &
    0.15707963262340149876D+01, &
    0.15707963267948756308D+01, &
    0.15707963267948966192D+01, &
    0.15707963267948966192D+01, &
    0.15707963267948966192D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0009765625D+00, &
      0.0039062500D+00, &
      0.0625000000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      2.0000000000D+00, &
      4.0000000000D+00, &
      5.0000000000D+00, &
      6.0000000000D+00, &
      6.5000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     12.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00, &
     80.0000000000D+00, &
    100.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_k1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_K1_VALUES returns some values of the K1 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    The modified Bessel function K1(Z) corresponds to N = 1.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselK[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.9853844780870606D+01, &
    0.4775972543220472D+01, &
    0.2184354424732687D+01, &
    0.1302834939763502D+01, &
    0.8617816344721803D+00, &
    0.6019072301972346D+00, &
    0.4345923910607150D+00, &
    0.3208359022298758D+00, &
    0.2406339113576119D+00, &
    0.1826230998017470D+00, &
    0.1398658818165224D+00, &
    0.7389081634774706D-01, &
    0.4015643112819418D-01, &
    0.2223939292592383D-01, &
    0.1248349888726843D-01, &
    0.7078094908968090D-02, &
    0.4044613445452164D-02, &
    0.1343919717735509D-02, &
    0.1553692118050011D-03, &
    0.1864877345382558D-04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     0.2D+00, &
     0.4D+00, &
     0.6D+00, &
     0.8D+00, &
     1.0D+00, &
     1.2D+00, &
     1.4D+00, &
     1.6D+00, &
     1.8D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00, &
     4.5D+00, &
     5.0D+00, &
     6.0D+00, &
     8.0D+00, &
    10.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_kn_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_KN_VALUES returns some values of the Kn Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 * W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselK[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.4951242928773287D+02, &
    0.1624838898635177D+01, &
    0.2537597545660559D+00, &
    0.1214602062785638D+00, &
    0.6151045847174204D-01, &
    0.5308943712223460D-02, &
    0.2150981700693277D-04, &
    0.6329543612292228D-09, &
    0.7101262824737945D+01, &
    0.6473853909486342D+00, &
    0.8291768415230932D-02, &
    0.2725270025659869D-04, &
    0.3727936773826211D-22, &
    0.3609605896012407D+03, &
    0.9431049100596467D+01, &
    0.3270627371203186D-01, &
    0.5754184998531228D-04, &
    0.4367182254100986D-22, &
    0.1807132899010295D+09, &
    0.1624824039795591D+06, &
    0.9758562829177810D+01, &
    0.1614255300390670D-02, &
    0.9150988209987996D-22, &
    0.6294369360424535D+23, &
    0.5770856852700241D+17, &
    0.4827000520621485D+09, &
    0.1787442782077055D+03, &
    0.1706148379722035D-20 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) nu
  integer ( kind = 4 ), save, dimension ( n_max ) :: nu_vec = (/ &
     2,  2,  2,  2, &
     2,  2,  2,  2, &
     3,  3,  3,  3, &
     3,  5,  5,  5, &
     5,  5, 10, 10, &
    10, 10, 10, 20, &
    20, 20, 20, 20 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_kx_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_KX_VALUES returns some values of the Kx Bessel function.
!
!  Discussion:
!
!    This set of data considers the less common case in which the
!    index of the Bessel function Kn is actually not an integer.
!    We may suggest this case by occasionally replacing the symbol
!    "Kn" by "Kx".
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselK[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
       2.294489339798475D+00, &
       0.4610685044478946D+00, &
       0.1199377719680614D+00, &
       0.06506594315400999D+00, &
       0.03602598513176459D+00, &
       0.003776613374642883D+00, &
       0.00001799347809370518D+00, &
       5.776373974707445D-10, &
       0.9221370088957891D+00, &
       0.1799066579520922D+00, &
       0.004531936049571459D+00, &
       0.00001979282590307570D+00, &
       3.486992497366216D-23, &
       3.227479531135262D+00, &
       0.3897977588961997D+00, &
       0.006495775004385758D+00, &
       0.00002393132586462789D+00, &
       3.627839645299048D-23, &
       0.7311451879202114D+00, &
       0.1567475478393932D+00, &
       0.004257389528177461D+00, &
       0.00001915541065869563D+00, &
       3.463337593569306D-23, &
       4.731184839919541D+00, &
       0.4976876225514758D+00, &
       0.007300864610941163D+00, &
       0.00002546421294106458D+00, &
       3.675275677913656D-23 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 8 ), save, dimension ( n_max ) :: nu_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_y0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_Y0_VALUES returns some values of the Y0 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1534238651350367D+01, &
     0.8825696421567696D-01, &
     0.5103756726497451D+00, &
     0.3768500100127904D+00, &
    -0.1694073932506499D-01, &
    -0.3085176252490338D+00, &
    -0.2881946839815792D+00, &
    -0.2594974396720926D-01, &
     0.2235214893875662D+00, &
     0.2499366982850247D+00, &
     0.5567116728359939D-01, &
    -0.1688473238920795D+00, &
    -0.2252373126343614D+00, &
    -0.7820786452787591D-01, &
     0.1271925685821837D+00, &
     0.2054642960389183D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_y0_int_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_Y0_INT_VALUES returns some values of the Bessel Y0 integral.
!
!  Discussion:
!
!    The function is defined by:
!
!      Y0_INT(x) = integral ( 0 <= t <= x ) Y0(t) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.91442642860172110926D-02, &
    -0.29682047390397591290D-01, &
    -0.25391431276585388961D+00, &
    -0.56179545591464028187D+00, &
    -0.63706937660742309754D+00, &
    -0.28219285008510084123D+00, &
     0.38366964785312561103D+00, &
    -0.12595061285798929390D+00, &
     0.24129031832266684828D+00, &
     0.17138069757627037938D+00, &
     0.18958142627134083732D+00, &
     0.17203846136449706946D+00, &
    -0.16821597677215029611D+00, &
    -0.93607927351428988679D-01, &
     0.88229711948036648408D-01, &
    -0.89324662736274161841D-02, &
    -0.54814071000063488284D-01, &
    -0.94958246003466381588D-01, &
    -0.19598064853404969850D-01, &
    -0.83084772357154773468D-02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0078125000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      2.0000000000D+00, &
      4.0000000000D+00, &
      6.0000000000D+00, &
     10.0000000000D+00, &
     16.0000000000D+00, &
     16.2500000000D+00, &
     17.0000000000D+00, &
     20.0000000000D+00, &
     25.0000000000D+00, &
     30.0000000000D+00, &
     40.0000000000D+00, &
     50.0000000000D+00, &
     70.0000000000D+00, &
    100.0000000000D+00, &
    125.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_y0_spherical_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_Y0_SPHERICAL_VALUES: some values of the Spherical Bessel function y0.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi/(2*x)] * BesselY[1/2,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.9950041652780258D+01, &
    -0.4900332889206208D+01, &
    -0.2302652485007213D+01, &
    -0.1375559358182797D+01, &
    -0.8708833866839568D+00, &
    -0.5403023058681397D+00, &
    -0.3019647953972280D+00, &
    -0.1214051020716007D+00, &
     0.1824970143830545D-01, &
     0.1262233859406039D+00, &
     0.2080734182735712D+00, &
     0.2675005078433390D+00, &
     0.3072473814755190D+00, &
     0.3295725974495951D+00, &
     0.3365079788102351D+00, &
     0.3299974988668152D+00, &
     0.3119671174358603D+00, &
     0.2843524095821944D+00, &
     0.2490995600928186D+00, &
     0.2081493978722149D+00, &
     0.1634109052159030D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00, &
    3.2D+00, &
    3.4D+00, &
    3.6D+00, &
    3.8D+00, &
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_y1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_Y1_VALUES returns some values of the Y1 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.6458951094702027D+01, &
    -0.7812128213002887D+00, &
    -0.1070324315409375D+00, &
     0.3246744247918000D+00, &
     0.3979257105571000D+00, &
     0.1478631433912268D+00, &
    -0.1750103443003983D+00, &
    -0.3026672370241849D+00, &
    -0.1580604617312475D+00, &
     0.1043145751967159D+00, &
     0.2490154242069539D+00, &
     0.1637055374149429D+00, &
    -0.5709921826089652D-01, &
    -0.2100814084206935D+00, &
    -0.1666448418561723D+00, &
     0.2107362803687351D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_y1_spherical_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_Y1_SPHERICAL_VALUES: some values of the Spherical Bessel function y1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi/(2*x)] * BesselY[3/2,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1004987506942709D+03, &
    -0.2549501110000635D+02, &
    -0.6730177068289658D+01, &
    -0.3233669719296388D+01, &
    -0.1985299346979349D+01, &
    -0.1381773290676036D+01, &
    -0.1028336567803712D+01, &
    -0.7906105943286149D+00, &
    -0.6133274385019998D+00, &
    -0.4709023582986618D+00, &
    -0.3506120042760553D+00, &
    -0.2459072254437506D+00, &
    -0.1534232496148467D+00, &
    -0.7151106706610352D-01, &
     0.5427959479750482D-03, &
     0.6295916360231598D-01, &
     0.1157316440198251D+00, &
     0.1587922092967723D+00, &
     0.1921166676076864D+00, &
     0.2157913917934037D+00, &
     0.2300533501309578D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00, &
    3.2D+00, &
    3.4D+00, &
    3.6D+00, &
    3.8D+00, &
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_yn_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_YN_VALUES returns some values of the Yn Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1650682606816254D+01, &
    -0.6174081041906827D+00, &
     0.3676628826055245D+00, &
    -0.5868082442208615D-02, &
     0.9579316872759649D-01, &
    -0.2604058666258122D+03, &
    -0.9935989128481975D+01, &
    -0.4536948224911019D+00, &
     0.1354030476893623D+00, &
    -0.7854841391308165D-01, &
    -0.1216180142786892D+09, &
    -0.1291845422080393D+06, &
    -0.2512911009561010D+02, &
    -0.3598141521834027D+00, &
     0.5723897182053514D-02, &
    -0.4113970314835505D+23, &
    -0.4081651388998367D+17, &
    -0.5933965296914321D+09, &
    -0.1597483848269626D+04, &
     0.1644263394811578D-01 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) nu
  integer ( kind = 4 ), save, dimension ( n_max ) :: nu_vec = (/ &
     2,  2,  2,  2, &
     2,  5,  5,  5, &
     5,  5, 10, 10, &
    10, 10, 10, 20, &
    20, 20, 20, 20 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_yx_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_YX_VALUES returns some values of the Yx Bessel function.
!
!  Discussion:
!
!    This set of data considers the less common case in which the
!    index of the Bessel function Yn is actually not an integer.
!    We may suggest this case by occasionally replacing the symbol
!    "Yn" by "Yx".
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
       -1.748560416961876D+00, &
       -0.4310988680183761D+00, &
        0.2347857104062485D+00, &
        0.4042783022390569D+00, &
        0.4560488207946332D+00, &
       -0.1012177091851084D+00, &
        0.2117088663313982D+00, &
       -0.07280690478506185D+00, &
       -1.102495575160179D+00, &
       -0.3956232813587035D+00, &
        0.3219244429611401D+00, &
        0.1584346223881903D+00, &
        0.02742813676191382D+00, &
       -2.876387857462161D+00, &
       -0.8282206324443037D+00, &
        0.2943723749617925D+00, &
       -0.1641784796149411D+00, &
        0.1105304445562544D+00, &
       -0.9319659251969881D+00, &
       -0.2609445010948933D+00, &
        0.2492796362185881D+00, &
        0.2174410301416733D+00, &
       -0.01578576650557229D+00, &
       -4.023453301501028D+00, &
       -0.9588998694752389D+00, &
        0.2264260361047367D+00, &
       -0.2193617736566760D+00, &
        0.09413988344515077D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 8 ), save, dimension ( n_max ) :: nu_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine beta_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BETA_CDF_VALUES returns some values of the Beta CDF.
!
!  Discussion:
!
!    The incomplete Beta function may be written
!
!      BETA_INC(A,B,X) = integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
!                      / integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
!
!    Thus,
!
!      BETA_INC(A,B,0.0) = 0.0
!      BETA_INC(A,B,1.0) = 1.0
!
!    The incomplete Beta function is also sometimes called the
!    "modified" Beta function, or the "normalized" Beta function
!    or the Beta CDF (cumulative density function).
!
!    In Mathematica, the function can be evaluated by:
!
!      BETA[X,A,B] / BETA[A,B]
!
!    The function can also be evaluated by using the Statistics package:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = BetaDistribution [ a, b ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968,
!    LC: QA351.P38.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 42

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.5D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    30.0D+00, &
    30.0D+00, &
    40.0D+00, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.2D+01, &
     0.3D+01, &
     0.4D+01, &
     0.5D+01 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.0D+00, &
     0.5D+00, &
     5.0D+00, &
     5.0D+00, &
    10.0D+00, &
     5.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.2D+01, &
     0.3D+01, &
     0.4D+01, &
     0.5D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6376856085851985D-01, &
    0.2048327646991335D+00, &
    0.1000000000000000D+01, &
    0.0000000000000000D+00, &
    0.5012562893380045D-02, &
    0.5131670194948620D-01, &
    0.2928932188134525D+00, &
    0.5000000000000000D+00, &
    0.2800000000000000D-01, &
    0.1040000000000000D+00, &
    0.2160000000000000D+00, &
    0.3520000000000000D+00, &
    0.5000000000000000D+00, &
    0.6480000000000000D+00, &
    0.7840000000000000D+00, &
    0.8960000000000000D+00, &
    0.9720000000000000D+00, &
    0.4361908850559777D+00, &
    0.1516409096347099D+00, &
    0.8978271484375000D-01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.4598773297575791D+00, &
    0.2146816102371739D+00, &
    0.9507364826957875D+00, &
    0.5000000000000000D+00, &
    0.8979413687105918D+00, &
    0.2241297491808366D+00, &
    0.7586405487192086D+00, &
    0.7001783247477069D+00, &
    0.5131670194948620D-01, &
    0.1055728090000841D+00, &
    0.1633399734659245D+00, &
    0.2254033307585166D+00, &
    0.3600000000000000D+00, &
    0.4880000000000000D+00, &
    0.5904000000000000D+00, &
    0.6723200000000000D+00, &
    0.2160000000000000D+00, &
    0.8370000000000000D-01, &
    0.3078000000000000D-01, &
    0.1093500000000000D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.10D+00, &
    1.00D+00, &
    0.00D+00, &
    0.01D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    0.50D+00, &
    0.90D+00, &
    0.50D+00, &
    1.00D+00, &
    0.50D+00, &
    0.80D+00, &
    0.60D+00, &
    0.80D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.70D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine beta_inc_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BETA_INC_VALUES returns some values of the incomplete Beta function.
!
!  Discussion:
!
!    The incomplete Beta function may be written
!
!      BETA_INC(A,B,X) = integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
!                      / integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
!
!    Thus,
!
!      BETA_INC(A,B,0.0) = 0.0
!      BETA_INC(A,B,1.0) = 1.0
!
!    The incomplete Beta function is also sometimes called the
!    "modified" Beta function, or the "normalized" Beta function
!    or the Beta CDF (cumulative density function).
!
!    In Mathematica, the function can be evaluated by:
!
!      BETA[X,A,B] / BETA[A,B]
!
!    The function can also be evaluated by using the Statistics package:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = BetaDistribution [ a, b ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968,
!    LC: QA351.P38.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 42

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.5D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    30.0D+00, &
    30.0D+00, &
    40.0D+00, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.2D+01, &
     0.3D+01, &
     0.4D+01, &
     0.5D+01 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.0D+00, &
     0.5D+00, &
     5.0D+00, &
     5.0D+00, &
    10.0D+00, &
     5.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.2D+01, &
     0.3D+01, &
     0.4D+01, &
     0.5D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6376856085851985D-01, &
    0.2048327646991335D+00, &
    0.1000000000000000D+01, &
    0.0000000000000000D+00, &
    0.5012562893380045D-02, &
    0.5131670194948620D-01, &
    0.2928932188134525D+00, &
    0.5000000000000000D+00, &
    0.2800000000000000D-01, &
    0.1040000000000000D+00, &
    0.2160000000000000D+00, &
    0.3520000000000000D+00, &
    0.5000000000000000D+00, &
    0.6480000000000000D+00, &
    0.7840000000000000D+00, &
    0.8960000000000000D+00, &
    0.9720000000000000D+00, &
    0.4361908850559777D+00, &
    0.1516409096347099D+00, &
    0.8978271484375000D-01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.4598773297575791D+00, &
    0.2146816102371739D+00, &
    0.9507364826957875D+00, &
    0.5000000000000000D+00, &
    0.8979413687105918D+00, &
    0.2241297491808366D+00, &
    0.7586405487192086D+00, &
    0.7001783247477069D+00, &
    0.5131670194948620D-01, &
    0.1055728090000841D+00, &
    0.1633399734659245D+00, &
    0.2254033307585166D+00, &
    0.3600000000000000D+00, &
    0.4880000000000000D+00, &
    0.5904000000000000D+00, &
    0.6723200000000000D+00, &
    0.2160000000000000D+00, &
    0.8370000000000000D-01, &
    0.3078000000000000D-01, &
    0.1093500000000000D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.10D+00, &
    1.00D+00, &
    0.00D+00, &
    0.01D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    0.50D+00, &
    0.90D+00, &
    0.50D+00, &
    1.00D+00, &
    0.50D+00, &
    0.80D+00, &
    0.60D+00, &
    0.80D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.70D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine beta_log_values ( n_data, x, y, fxy )

!*****************************************************************************80
!
!! BETA_LOG_VALUES returns some values of the logarithm of the Beta function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Log[Beta[x]]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, Y, the arguments of the function.
!
!    Output, real ( kind = 8 ) FXY, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fxy
  real ( kind = 8 ), save, dimension ( n_max ) :: fxy_vec = (/ &
     0.1609437912434100D+01, &
     0.9162907318741551D+00, &
     0.5108256237659907D+00, &
     0.2231435513142098D+00, &
     0.1609437912434100D+01, &
     0.9162907318741551D+00, &
     0.0000000000000000D+00, &
    -0.1791759469228055D+01, &
    -0.3401197381662155D+01, &
    -0.4941642422609304D+01, &
    -0.6445719819385578D+01, &
    -0.3737669618283368D+01, &
    -0.5123963979403259D+01, &
    -0.6222576268071369D+01, &
    -0.7138866999945524D+01, &
    -0.7927324360309794D+01, &
    -0.9393661429103221D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    7.0D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( n_max ) :: y_vec = (/ &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    0.2D+00, &
    0.4D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    6.0D+00, &
    7.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    y = 0.0D+00
    fxy = 0.0D+00
  else
    x = x_vec(n_data)
    y = y_vec(n_data)
    fxy = fxy_vec(n_data)
  end if

  return
end
subroutine beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

!*****************************************************************************80
!
!! BETA_NONCENTRAL_CDF_VALUES returns some values of the noncentral Beta CDF.
!
!  Discussion:
!
!    The values presented here are taken from the reference, where they
!    were given to a limited number of decimal places.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R Chattamvelli, R Shanmugam,
!    Algorithm AS 310:
!    Computing the Non-central Beta Distribution Function,
!    Applied Statistics,
!    Volume 46, Number 1, 1997, pages 146-156.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the shape parameters.
!
!    Output, real ( kind = 8 ) LAMBDA, the noncentrality parameter.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
        5.0D+00, &
        5.0D+00, &
        5.0D+00, &
       10.0D+00, &
       10.0D+00, &
       10.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       10.0D+00, &
       10.0D+00, &
       15.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       30.0D+00, &
       30.0D+00, &
       10.0D+00, &
       10.0D+00, &
       10.0D+00, &
       15.0D+00, &
       10.0D+00, &
       12.0D+00, &
       30.0D+00, &
       35.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
        5.0D+00, &
        5.0D+00, &
        5.0D+00, &
       10.0D+00, &
       10.0D+00, &
       10.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       10.0D+00, &
        5.0D+00, &
       10.0D+00, &
       30.0D+00, &
       50.0D+00, &
       20.0D+00, &
       40.0D+00, &
        5.0D+00, &
       10.0D+00, &
       30.0D+00, &
       20.0D+00, &
        5.0D+00, &
       17.0D+00, &
       30.0D+00, &
       30.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
       0.4563021D+00, &
       0.1041337D+00, &
       0.6022353D+00, &
       0.9187770D+00, &
       0.6008106D+00, &
       0.0902850D+00, &
       0.9998655D+00, &
       0.9925997D+00, &
       0.9641112D+00, &
       0.9376626573D+00, &
       0.7306817858D+00, &
       0.1604256918D+00, &
       0.1867485313D+00, &
       0.6559386874D+00, &
       0.9796881486D+00, &
       0.1162386423D+00, &
       0.9930430054D+00, &
       0.0506899273D+00, &
       0.1030959706D+00, &
       0.9978417832D+00, &
       0.2555552369D+00, &
       0.0668307064D+00, &
       0.0113601067D+00, &
       0.7813366615D+00, &
       0.8867126477D+00 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
        54.0D+00, &
       140.0D+00, &
       170.0D+00, &
        54.0D+00, &
       140.0D+00, &
       250.0D+00, &
        54.0D+00, &
       140.0D+00, &
       250.0D+00, &
       150.0D+00, &
       120.0D+00, &
        80.0D+00, &
       110.0D+00, &
        65.0D+00, &
       130.0D+00, &
        80.0D+00, &
       130.0D+00, &
        20.0D+00, &
        54.0D+00, &
        80.0D+00, &
       120.0D+00, &
        55.0D+00, &
        64.0D+00, &
       140.0D+00, &
        20.0D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
       0.8640D+00, &
       0.9000D+00, &
       0.9560D+00, &
       0.8686D+00, &
       0.9000D+00, &
       0.9000D+00, &
       0.8787D+00, &
       0.9000D+00, &
       0.9220D+00, &
       0.868D+00, &
       0.900D+00, &
       0.880D+00, &
       0.850D+00, &
       0.660D+00, &
       0.720D+00, &
       0.720D+00, &
       0.800D+00, &
       0.644D+00, &
       0.700D+00, &
       0.780D+00, &
       0.760D+00, &
       0.795D+00, &
       0.560D+00, &
       0.800D+00, &
       0.670D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine beta_values ( n_data, x, y, fxy )

!*****************************************************************************80
!
!! BETA_VALUES returns some values of the Beta function.
!
!  Discussion:
!
!    Beta(X,Y) = ( Gamma(X) * Gamma(Y) ) / Gamma(X+Y)
!
!    Both X and Y must be greater than 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      Beta[X,Y]
!
!  Properties:
!
!    Beta(X,Y) = Beta(Y,X).
!    Beta(X,Y) = integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
!    Beta(X,Y) = Gamma(X) * Gamma(Y) / Gamma(X+Y)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, Y, the arguments of the function.
!
!    Output, real ( kind = 8 ) FXY, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    0.5000000000000000D+01, &
    0.2500000000000000D+01, &
    0.1666666666666667D+01, &
    0.1250000000000000D+01, &
    0.5000000000000000D+01, &
    0.2500000000000000D+01, &
    0.1000000000000000D+01, &
    0.1666666666666667D+00, &
    0.3333333333333333D-01, &
    0.7142857142857143D-02, &
    0.1587301587301587D-02, &
    0.2380952380952381D-01, &
    0.5952380952380952D-02, &
    0.1984126984126984D-02, &
    0.7936507936507937D-03, &
    0.3607503607503608D-03, &
    0.8325008325008325D-04 /)
  real ( kind = 8 ) fxy
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    7.0D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( n_max ) :: y_vec = (/ &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    0.2D+00, &
    0.4D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    6.0D+00, &
    7.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    y = 0.0D+00
    fxy = 0.0D+00
  else
    x = x_vec(n_data)
    y = y_vec(n_data)
    fxy = b_vec(n_data)
  end if

  return
end
subroutine binomial_values ( n_data, a, b, fx )

!*****************************************************************************80
!
!! BINOMIAL_VALUES returns some values of the binomial coefficients.
!
!  Discussion:
!
!    The formula for the binomial coefficient is
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!    In Mathematica, the function can be evaluated by:
!
!      Binomial[n,k]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) A, B, the arguments of the function.
!
!    Output, integer ( kind = 4 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) a
  integer ( kind = 4 ), save, dimension ( n_max ) :: a_vec = (/ &
     1,  6,  6,  6, 15, &
    15, 15, 15, 15, 15, &
    15, 25, 25, 25, 25, &
    25, 25, 25, 25, 25 /)
  integer ( kind = 4 ) b
  integer ( kind = 4 ), save, dimension ( n_max ) :: b_vec = (/ &
     0,  1,  3,  5,  1, &
     3,  5,  7,  9, 11, &
    13,  1,  3,  5,  7, &
     9, 11, 13, 15, 17 /)
  integer ( kind = 4 ) fx
  integer ( kind = 4 ), save, dimension ( n_max ) :: fx_vec = (/ &
           1, &
           6, &
          20, &
           6, &
          15, &
         455, &
        3003, &
        6435, &
        5005, &
        1365, &
         105, &
          25, &
        2300, &
       53130, &
      480700, &
     2042975, &
     4457400, &
     5200300, &
     3268760, &
     1081575 /)
  integer ( kind = 4 ) n_data

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0
    b = 0
    fx = 0
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine binomial_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
!
!  Discussion:
!
!    CDF(X)(A,B) is the probability of at most X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = BinomialDistribution [ n, p ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) A, a parameter of the function.
!
!    Output, real ( kind = 8 ) B, a parameter of the function.
!
!    Output, integer ( kind = 4 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  integer ( kind = 4 ) a
  integer ( kind = 4 ), save, dimension ( n_max ) :: a_vec = (/ &
     2,  2,  2,  2, &
     2,  4,  4,  4, &
     4, 10, 10, 10, &
    10, 10, 10, 10, &
    10 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    0.05D+00, &
    0.05D+00, &
    0.05D+00, &
    0.50D+00, &
    0.50D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.05D+00, &
    0.10D+00, &
    0.15D+00, &
    0.20D+00, &
    0.25D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.9025000000000000D+00, &
    0.9975000000000000D+00, &
    0.1000000000000000D+01, &
    0.2500000000000000D+00, &
    0.7500000000000000D+00, &
    0.3164062500000000D+00, &
    0.7382812500000000D+00, &
    0.9492187500000000D+00, &
    0.9960937500000000D+00, &
    0.9999363101685547D+00, &
    0.9983650626000000D+00, &
    0.9901259090013672D+00, &
    0.9672065024000000D+00, &
    0.9218730926513672D+00, &
    0.8497316674000000D+00, &
    0.6331032576000000D+00, &
    0.3769531250000000D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) x
  integer ( kind = 4 ), save, dimension ( n_max ) :: x_vec = (/ &
     0, 1, 2, 0, &
     1, 0, 1, 2, &
     3, 4, 4, 4, &
     4, 4, 4, 4, &
     4 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0
    b = 0.0D+00
    x = 0
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bivariate_normal_cdf_values ( n_data, x, y, r, fxy )

!*****************************************************************************80
!
!! BIVARIATE_NORMAL_CDF_VALUES returns some values of the bivariate normal CDF.
!
!  Discussion:
!
!    FXY is the probability that two variables A and B, which are
!    related by a bivariate normal distribution with correlation R,
!    respectively satisfy A <= X and B <= Y.
!
!    Mathematica can evaluate the bivariate normal CDF via the commands:
!
!      <<MultivariateStatistics`
!      cdf = CDF[MultinormalDistribution[{0,0}{{1,r},{r,1}}],{x,y}]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    National Bureau of Standards,
!    Tables of the Bivariate Normal Distribution and Related Functions,
!    Applied Mathematics Series, Number 50, 1959.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, Y, the parameters of the function.
!
!    Output, real ( kind = 8 ) R, the correlation value.
!
!    Output, real ( kind = 8 ) FXY, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 41

  real ( kind = 8 ) fxy
  real ( kind = 8 ), save, dimension ( n_max ) :: fxy_vec = (/ &
  0.02260327218569867D+00, &
  0.1548729518584100D+00, &
  0.4687428083352184D+00, &
  0.7452035868929476D+00, &
  0.8318608306874188D+00, &
  0.8410314261134202D+00, &
  0.1377019384919464D+00, &
  0.1621749501739030D+00, &
  0.1827411243233119D+00, &
  0.2010067421506235D+00, &
  0.2177751155265290D+00, &
  0.2335088436446962D+00, &
  0.2485057781834286D+00, &
  0.2629747825154868D+00, &
  0.2770729823404738D+00, &
  0.2909261168683812D+00, &
  0.3046406378726738D+00, &
  0.3183113449213638D+00, &
  0.3320262544108028D+00, &
  0.3458686754647614D+00, &
  0.3599150462310668D+00, &
  0.3742210899871168D+00, &
  0.3887706405282320D+00, &
  0.4032765198361344D+00, &
  0.4162100291953678D+00, &
  0.6508271498838664D+00, &
  0.8318608306874188D+00, &
  0.0000000000000000D+00, &
  0.1666666666539970D+00, &
  0.2500000000000000D+00, &
  0.3333333333328906D+00, &
  0.5000000000000000D+00, &
  0.7452035868929476D+00, &
  0.1548729518584100D+00, &
  0.1548729518584100D+00, &
  0.06251409470431653D+00, &
  0.7452035868929476D+00, &
  0.1548729518584100D+00, &
  0.1548729518584100D+00, &
  0.06251409470431653D+00, &
  0.6337020457912916D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r
  real ( kind = 8 ), save, dimension ( n_max ) :: r_vec = (/ &
     0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00, &
     0.500D+00, -0.900D+00, -0.800D+00, -0.700D+00, -0.600D+00, &
    -0.500D+00, -0.400D+00, -0.300D+00, -0.200D+00, -0.100D+00, &
     0.000D+00,  0.100D+00,  0.200D+00,  0.300D+00,  0.400D+00, &
     0.500D+00,  0.600D+00,  0.700D+00,  0.800D+00,  0.900D+00, &
     0.673D+00,  0.500D+00, -1.000D+00, -0.500D+00,  0.000D+00, &
     0.500D+00,  1.000D+00,  0.500D+00,  0.500D+00,  0.500D+00, &
     0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00, &
     0.500D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -2.0D+00, -1.0D+00,  0.0D+00,  1.0D+00,  2.0D+00, &
     3.0D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
    -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
    -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
    -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
     1.0D+00,  2.0D+00,  0.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00, -1.0D+00, &
    -1.0D+00,  1.0D+00,  1.0D+00, -1.0D+00, -1.0D+00, &
     0.7071067811865475D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( n_max ) :: y_vec = (/ &
     1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
     1.0D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  1.0D+00,  0.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00, -1.0D+00,  1.0D+00, &
    -1.0D+00,  1.0D+00, -1.0D+00,  1.0D+00, -1.0D+00, &
     0.7071067811865475D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    r = 0.0D+00
    x = 0.0D+00
    y = 0.0D+00
    fxy = 0.0D+00
  else
    r = r_vec(n_data)
    x = x_vec(n_data)
    y = y_vec(n_data)
    fxy = fxy_vec(n_data)
  end if

  return
end
subroutine catalan_values ( n_data, n, c )

!*****************************************************************************80
!
!! CATALAN_VALUES returns some values of the Catalan numbers.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Binomial[2*n,n] / ( n + 1 )
!
!  First values:
!
!     C(0)     1
!     C(1)     1
!     C(2)     2
!     C(3)     5
!     C(4)    14
!     C(5)    42
!     C(6)   132
!     C(7)   429
!     C(8)  1430
!     C(9)  4862
!    C(10) 16796
!
!  Formula:
!
!    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
!         = 1 / (N+1) * COMB ( 2N, N )
!         = 1 / (2N+1) * COMB ( 2N+1, N+1).
!
!  Recursion:
!
!    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
!    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
!
!  Discussion:
!
!    The Catalan number C(N) counts:
!
!    1) the number of binary trees on N vertices;
!    2) the number of ordered trees on N+1 vertices;
!    3) the number of full binary trees on 2N+1 vertices;
!    4) the number of well formed sequences of 2N parentheses;
!    5) the number of ways 2N ballots can be counted, in order,
!       with N positive and N negative, so that the running sum
!       is never negative;
!    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
!    7) the number of monotone functions from [1..N} to [1..N} which
!       satisfy f(i) <= i for all i;
!    8) the number of ways to triangulate a polygon with N+2 vertices.
!
!  Example:
!
!    N = 3
!
!    ()()()
!    ()(())
!    (()())
!    (())()
!    ((()))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the Catalan number.
!
!    Output, integer ( kind = 4 ) C, the value of the Catalan number.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
    1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine cauchy_cdf_values ( n_data, mu, sigma, x, fx )

!*****************************************************************************80
!
!! CAUCHY_CDF_VALUES returns some values of the Cauchy CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = CauchyDistribution [ mu, sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) SIGMA, the variance of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.8524163823495667D+00, &
    0.9220208696226307D+00, &
    0.9474315432887466D+00, &
    0.6475836176504333D+00, &
    0.6024163823495667D+00, &
    0.5779791303773693D+00, &
    0.5628329581890012D+00, &
    0.6475836176504333D+00, &
    0.5000000000000000D+00, &
    0.3524163823495667D+00, &
    0.2500000000000000D+00 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cbrt_values ( n_data, x, fx )

!*****************************************************************************80
!
!! CBRT_VALUES returns some values of the cube root function.
!
!  Discussion:
!
!    CBRT(X) = real number Y such that Y * Y * Y = X.
!
!    In Mathematica, the function can be evaluated by:
!
!      Sign[x] * ( Abs[x] )^(1/3)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 14

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
   -0.0020082988563383484484D+00, &
    0.44814047465571647087D+00, &
   -0.46415888336127788924D+00, &
    0.73680629972807732116D+00, &
   -1.0000000000000000000D+00, &
    1.2599210498948731648D+00, &
   -1.4422495703074083823D+00, &
    1.4645918875615232630D+00, &
   -2.6684016487219448673D+00, &
    3.0723168256858472933D+00, &
   -4.1408177494228532500D+00, &
    4.5947008922070398061D+00, &
   -497.93385921817447440D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
   -0.8100000073710001D-08, &
    0.9000000000000000D-01, &
   -0.1000000000000000D+00, &
    0.4000000000000000D+00, &
   -0.1000000000000000D+01, &
    0.2000000000000000D+01, &
   -0.3000000000000000D+01, &
    0.31415926535897932385D+01, &
   -0.1900000000000000D+02, &
    0.2900000000000000D+02, &
   -0.7100000000000000D+02, &
    0.9700000000000000D+02, &
   -0.1234567890000000D+09 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cheby_t_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! CHEBY_T_POLY_VALUES returns values of Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ChebyshevT[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.8000000000000000D+00, &
     0.2800000000000000D+00, &
    -0.3520000000000000D+00, &
    -0.8432000000000000D+00, &
    -0.9971200000000000D+00, &
    -0.7521920000000000D+00, &
    -0.2063872000000000D+00, &
     0.4219724800000000D+00, &
     0.8815431680000000D+00, &
     0.9884965888000000D+00, &
     0.7000513740800000D+00, &
     0.1315856097280000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cheby_u_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! CHEBY_U_POLY_VALUES returns values of Chebyshev polynomials U(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ChebyshevU[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.1600000000000000D+01, &
     0.1560000000000000D+01, &
     0.8960000000000000D+00, &
    -0.1264000000000000D+00, &
    -0.1098240000000000D+01, &
    -0.1630784000000000D+01, &
    -0.1511014400000000D+01, &
    -0.7868390400000000D+00, &
     0.2520719360000000D+00, &
     0.1190154137600000D+01, &
     0.1652174684160000D+01, &
     0.1453325357056000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cheby_v_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! CHEBY_V_POLY_VALUES returns values of Chebyshev polynomials V(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      u = Sqrt[(x+1)/2],
!      ChebyshevT[2*n+1,u] / u
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     1.0000000000000000D+00, &
     0.6000000000000000D+00, &
    -0.0400000000000000D+00, &
    -0.6640000000000000D+00, &
    -1.0224000000000000D+00, &
    -0.9718400000000000D+00, &
    -0.5325440000000000D+00, &
     0.1197696000000000D+00, &
     0.7241753600000000D+00, &
     1.0389109760000000D+00, &
     0.9380822016000000D+00, &
     0.4620205465600000D+00, &
    -0.1988493271040000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cheby_w_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! CHEBY_W_POLY_VALUES returns values of Chebyshev polynomials W(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      u = Sqrt[(x+1)/2],
!      ChebyshevU[2*n,u]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     1.000000000000000D+00, &
     2.600000000000000D+00, &
     3.160000000000000D+00, &
     2.456000000000000D+00, &
     0.769600000000000D+00, &
    -1.224640000000000D+00, &
    -2.729024000000000D+00, &
    -3.141798400000000D+00, &
    -2.297853440000000D+00, &
    -0.534767104000000D+00, &
     1.442226073600000D+00, &
     2.842328821760000D+00, &
     3.105500041216000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine chi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! CHI_VALUES returns some values of the hyperbolic cosine integral function.
!
!  Discussion:
!
!    The hyperbolic cosine integral is defined by
!
!      chi ( x ) = gamma + log ( x )
!        + integral ( 0 <= t < x ) ( cosh ( t ) - 1 ) / t  dt
!
!    where gamma is Euler's constant.
!
!    In Mathematica, the function can be evaluated by:
!
!      CoshIntegral[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.05277684495649362D+00, &
     0.1577508933739787D+00, &
     0.3455691756953907D+00, &
     0.5183999848333915D+00, &
     0.6813138871854339D+00, &
     0.8378669409802082D+00, &
     1.141841924170595D+00, &
     1.445494075789644D+00, &
     1.759505807660965D+00, &
     2.092577214062032D+00, &
     2.452666922646915D+00, &
     3.524425488354165D+00, &
     4.960392094765610D+00, &
     6.959191927647393D+00, &
     9.813547558823186D+00, &
    13.96581164859243D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     1.0D+00, &
     1.2D+00, &
     1.4D+00, &
     1.6D+00, &
     1.8D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00, &
     4.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine chi_square_cdf_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = ChiSquareDistribution [ df ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  integer ( kind = 4 ) a
  integer ( kind = 4 ), save, dimension ( n_max ) :: a_vec = (/ &
     1,  2,  1,  2, &
     1,  2,  3,  4, &
     1,  2,  3,  4, &
     5,  3,  3,  3, &
     3,  3, 10, 10, &
    10 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7965567455405796D-01, &
    0.4987520807317687D-02, &
    0.1124629160182849D+00, &
    0.9950166250831946D-02, &
    0.4729107431344619D+00, &
    0.1812692469220181D+00, &
    0.5975750516063926D-01, &
    0.1752309630642177D-01, &
    0.6826894921370859D+00, &
    0.3934693402873666D+00, &
    0.1987480430987992D+00, &
    0.9020401043104986D-01, &
    0.3743422675270363D-01, &
    0.4275932955291202D+00, &
    0.6083748237289110D+00, &
    0.7385358700508894D+00, &
    0.8282028557032669D+00, &
    0.8883897749052874D+00, &
    0.1721156299558408D-03, &
    0.3659846827343712D-02, &
    0.1857593622214067D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.01D+00, &
    0.02D+00, &
    0.02D+00, &
    0.40D+00, &
    0.40D+00, &
    0.40D+00, &
    0.40D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    3.00D+00, &
    4.00D+00, &
    5.00D+00, &
    6.00D+00, &
    1.00D+00, &
    2.00D+00, &
    3.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine chi_square_noncentral_cdf_values ( n_data, df, lambda, x, cdf )

!*****************************************************************************80
!
!! CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NoncentralChiSquareDistribution [ df, lambda ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) DF, the number of degrees of freedom.
!
!    Output, real ( kind = 8 ) LAMBDA, the noncentrality parameter.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) CDF, the noncentral chi CDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) cdf
  real ( kind = 8 ), save, dimension ( n_max ) :: cdf_vec = (/ &
    0.8399444269398261D+00, &
    0.6959060300435139D+00, &
    0.5350879697078847D+00, &
    0.7647841496310313D+00, &
    0.6206436532195436D+00, &
    0.4691667375373180D+00, &
    0.3070884345937569D+00, &
    0.2203818092990903D+00, &
    0.1500251895581519D+00, &
    0.3071163194335791D-02, &
    0.1763982670131894D-02, &
    0.9816792594625022D-03, &
    0.1651753140866208D-01, &
    0.2023419573950451D-03, &
    0.4984476352854074D-06, &
    0.1513252400654827D-01, &
    0.2090414910614367D-02, &
    0.2465021206048452D-03, &
    0.2636835050342939D-01, &
    0.1857983220079215D-01, &
    0.1305736595486640D-01, &
    0.5838039534819351D-01, &
    0.4249784402463712D-01, &
    0.3082137716021596D-01, &
    0.1057878223400849D+00, &
    0.7940842984598509D-01, &
    0.5932010895599639D-01, &
    0.2110395656918684D+00 /)
  integer ( kind = 4 ) df
  integer ( kind = 4 ), save, dimension ( n_max ) :: df_vec = (/ &
      1,   2,   3, &
      1,   2,   3, &
      1,   2,   3, &
      1,   2,   3, &
     60,  80, 100, &
      1,   2,   3, &
     10,  10,  10, &
     10,  10,  10, &
     10,  10,  10, &
      8 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
      0.5D+00, &
      0.5D+00, &
      0.5D+00, &
      1.0D+00, &
      1.0D+00, &
      1.0D+00, &
      5.0D+00, &
      5.0D+00, &
      5.0D+00, &
     20.0D+00, &
     20.0D+00, &
     20.0D+00, &
     30.0D+00, &
     30.0D+00, &
     30.0D+00, &
      5.0D+00, &
      5.0D+00, &
      5.0D+00, &
      2.0D+00, &
      3.0D+00, &
      4.0D+00, &
      2.0D+00, &
      3.0D+00, &
      4.0D+00, &
      2.0D+00, &
      3.0D+00, &
      4.0D+00, &
      0.5D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
    60.000D+00, &
    60.000D+00, &
    60.000D+00, &
     0.050D+00, &
     0.050D+00, &
     0.050D+00, &
     4.000D+00, &
     4.000D+00, &
     4.000D+00, &
     5.000D+00, &
     5.000D+00, &
     5.000D+00, &
     6.000D+00, &
     6.000D+00, &
     6.000D+00, &
     5.000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    lambda = 0.0D+00
    df = 0
    cdf = 0.0D+00
  else
    x = x_vec(n_data)
    lambda = lambda_vec(n_data)
    df = df_vec(n_data)
    cdf = cdf_vec(n_data)
  end if

  return
end
subroutine ci_values ( n_data, x, fx )

!*****************************************************************************80
!
!! CI_VALUES returns some values of the cosine integral function.
!
!  Discussion:
!
!    The cosine integral is defined by
!
!      ci ( x ) = - integral ( x <= t < +oo ) ( cos ( t ) ) / t  dt
!
!    In Mathematica, the function can be evaluated by:
!
!      CosIntegral[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1777840788066129D+00, &
    -0.2227070695927976D-01, &
     0.1005147070088978D+00, &
     0.1982786159524672D+00, &
     0.2760678304677729D+00, &
     0.3374039229009681D+00, &
     0.4204591828942405D+00, &
     0.4620065850946773D+00, &
     0.4717325169318778D+00, &
     0.4568111294183369D+00, &
     0.4229808287748650D+00, &
     0.2858711963653835D+00, &
     0.1196297860080003D+00, &
    -0.3212854851248112D-01, &
    -0.1409816978869304D+00, &
    -0.1934911221017388D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     1.0D+00, &
     1.2D+00, &
     1.4D+00, &
     1.6D+00, &
     1.8D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00, &
     4.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cin_values ( n_data, x, fx )

!*****************************************************************************80
!
!! CIN_VALUES returns some values of the alternate cosine integral function.
!
!  Discussion:
!
!    The alternate cosine integral is defined by
!
!      cin ( x ) = gamma + log ( x )
!        + integral ( 0 <= t <= x ) ( cos ( t ) - 1 ) / t  dt
!
!    where gamma is Euler's constant.
!
!    In Mathematica, the function can be evaluated by:
!
!      EulerGamma + Log[x] - CosIntegral[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6185256314820045D-01, &
    0.8866074809482194D-01, &
    0.1200260139539026D+00, &
    0.1557934976348559D+00, &
    0.1957873187759337D+00, &
    0.2398117420005647D+00, &
    0.3390780388012470D+00, &
    0.4516813164280685D+00, &
    0.5754867772153906D+00, &
    0.7081912003853150D+00, &
    0.8473820166866132D+00, &
    0.1207635200410304D+01, &
    0.1556198167561642D+01, &
    0.1862107181909382D+01, &
    0.2104491723908354D+01, &
    0.2274784183779546D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.5D+00, &
    3.0D+00, &
    3.5D+00, &
    4.0D+00, &
    4.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cinh_values ( n_data, x, fx )

!*****************************************************************************80
!
!! CINH_VALUES returns some values of the alternate hyperbolic cosine integral.
!
!  Discussion:
!
!    The alternate hyperbolic cosine integral is defined by
!
!      CINH(X) =integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
!
!    In Mathematica, the function can be evaluated by:
!
!      Integrate [ ( Cosh[t] - 1 ) / t, { t, 0, x } ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.00000000000000000D+00, &
     0.06315467070191883D+00, &
     0.09136085223843649D+00, &
     0.1250284547325902D+00, &
     0.1643278712460683D+00, &
     0.2094587379417273D+00, &
     0.2606512760786754D+00, &
     0.3823047024751071D+00, &
     0.5318061742668980D+00, &
     0.7122865135136963D+00, &
     0.9275748842583805D+00, &
     1.182304077185436D+00, &
     2.030919091578478D+00, &
     3.284564141195967D+00, &
     5.129213294250493D+00, &
     7.850037532801762D+00, &
    11.88451858691463D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.0D+00, &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     1.0D+00, &
     1.2D+00, &
     1.4D+00, &
     1.6D+00, &
     1.8D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00, &
     4.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine clausen_values ( n_data, x, fx )

!*****************************************************************************80
!
!! CLAUSEN_VALUES returns some values of the Clausen's integral.
!
!  Discussion:
!
!    The function is defined by:
!
!      CLAUSEN(x) = integral ( 0 <= t <= x ) -ln ( 2 * sin ( t / 2 ) ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.14137352886760576684D-01, &
     0.13955467081981281934D+00, &
    -0.38495732156574238507D+00, &
     0.84831187770367927099D+00, &
     0.10139591323607685043D+01, &
    -0.93921859275409211003D+00, &
     0.72714605086327924743D+00, &
     0.43359820323553277936D+00, &
    -0.98026209391301421161D-01, &
    -0.56814394442986978080D+00, &
    -0.70969701784448921625D+00, &
     0.99282013254695671871D+00, &
    -0.98127747477447367875D+00, &
    -0.64078266570172320959D+00, &
     0.86027963733231192456D+00, &
     0.39071647608680211043D+00, &
     0.47574793926539191502D+00, &
     0.10105014481412878253D+01, &
     0.96332089044363075154D+00, &
    -0.61782699481929311757D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
     -0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
     -1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
     -3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
     -5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
    -10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
    -30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine clebsch_gordan_values ( n_data, j1, j2, j3, m1, m2, m3, fx )

!*****************************************************************************80
!
!! CLEBSCH_GORDAN_VALUES returns some values of the Clebsch-Gordan function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ClebschGordan[{j1,m1},{j2,m2},{j3,m3}]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) J1, J2, J3, M1, M2, M3, the arguments
!    of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.7071067811865475D+00, &
     1.000000000000000D+00, &
     0.5773502691896258D+00, &
    -0.2581988897471611D+00, &
    -0.6324555320336759D+00, &
    -0.7745966692414834D+00, &
     0.4082482904638630D+00, &
     0.8164965809277260D+00, &
     0.5345224838248488D+00, &
     0.2672612419124244D+00, &
     0.8944271909999159D+00, &
     0.3380617018914066D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) j1
  real ( kind = 8 ), save, dimension ( n_max ) :: j1_vec = (/ &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    2.0D+00, &
    1.5D+00, &
    1.5D+00 /)
  real ( kind = 8 ) j2
  real ( kind = 8 ), save, dimension ( n_max ) :: j2_vec = (/ &
    0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    1.5D+00, &
    1.5D+00, &
    1.5D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00 /)
  real ( kind = 8 ) j3
  real ( kind = 8 ), save, dimension ( n_max ) :: j3_vec = (/ &
    1.0D+00, &
    1.0D+00, &
    1.5D+00, &
    1.5D+00, &
    1.5D+00, &
    1.5D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.5D+00, &
    3.5D+00 /)
  real ( kind = 8 ) m1
  real ( kind = 8 ), save, dimension ( n_max ) :: m1_vec = (/ &
    0.5D+00, &
    0.5D+00, &
   -0.5D+00, &
    0.0D+00, &
   -1.0D+00, &
    0.0D+00, &
    1.0D+00, &
    0.0D+00, &
    2.0D+00, &
    1.0D+00, &
    0.5D+00, &
    1.5D+00 /)
  real ( kind = 8 ) m2
  real ( kind = 8 ), save, dimension ( n_max ) :: m2_vec = (/ &
   -0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    0.5D+00, &
    1.5D+00, &
    1.5D+00, &
   -1.0D+00, &
    0.0D+00, &
   -2.0D+00, &
   -1.0D+00, &
    1.0D+00, &
   -1.0D+00 /)
  real ( kind = 8 ) m3
  real ( kind = 8 ), save, dimension ( n_max ) :: m3_vec = (/ &
    0.0D+00, &
    1.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    1.5D+00, &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    1.5D+00, &
    0.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    j1 = 0.0D+00
    j2 = 0.0D+00
    j3 = 0.0D+00
    m1 = 0.0D+00
    m2 = 0.0D+00
    m3 = 0.0D+00
    fx = 0.0D+00
  else
    j1 = j1_vec(n_data)
    j2 = j2_vec(n_data)
    j3 = j3_vec(n_data)
    m1 = m1_vec(n_data)
    m2 = m2_vec(n_data)
    m3 = m3_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine collatz_count_values ( n_data, n, count )

!*****************************************************************************80
!
!! COLLATZ_COUNT_VALUES returns some values of the Collatz count function.
!
!  Discussion:
!
!    The rules for generation of the Collatz sequence are recursive.
!    If T is the current entry of the sequence, (T is
!    assumed to be a positive integer), then the next
!    entry, U is determined as follows:
!
!      if T is 1 (or less)
!        terminate the sequence;
!      else if T is even
!        U = T/2.
!      else (if T is odd and not 1)
!        U = 3*T+1;
!
!    The Collatz count is the length of the Collatz sequence for a given
!    starting value.  By convention, we include the initial value in the
!    count, so the minimum value of the count is 1.
!
!     N  Sequence                                                 Count
!
!     1                                                               1
!     2   1                                                           2
!     3  10,  5, 16,  8,  4,  2,  1                                   8
!     4   2   1                                                       3
!     5  16,  8,  4,  2,  1                                           6
!     6   3, 10,  5, 16,  8,  4,  2,  1                               9
!     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
!     8   4,  2,  1                                                   4
!     9  28, 14,  7, ...                                             20
!    10   5, 16,  8,  4,  2,  1                                       7
!    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
!    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the initial value of a Collatz sequence.
!
!    Output, integer ( kind = 4 ) COUNT, the length of the Collatz sequence
!    starting with N.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) count
  integer ( kind = 4 ), save, dimension ( n_max ) :: count_vec = (/ &
      1,   2,   8,   3,   6,   9,   17,   4,  20,   7, &
    112,  25,  26,  27,  17,  28,  111,  18,  83,  29 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     27,  50, 100, 200, 300, 400, 500, 600, 700, 800 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    count = 0
  else
    n = n_vec(n_data)
    count = count_vec(n_data)
  end if

  return
end
subroutine cos_values ( n_data, x, fx )

!*****************************************************************************80
!
!! COS_VALUES returns some values of the cosine function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Cos[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     1.0000000000000000000D+00, &
     0.96592582628906828675D+00, &
     0.87758256189037271612D+00, &
     0.86602540378443864676D+00, &
     0.70710678118654752440D+00, &
     0.54030230586813971740D+00, &
     0.50000000000000000000D+00, &
     0.00000000000000000000D+00, &
    -0.41614683654714238700D+00, &
    -0.98999249660044545727D+00, &
    -1.0000000000000000000D+00, &
    -0.65364362086361191464D+00, &
     0.28366218546322626447D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000000D+00, &
    0.26179938779914943654D+00, &
    0.50000000000000000000D+00, &
    0.52359877559829887308D+00, &
    0.78539816339744830962D+00, &
    1.0000000000000000000D+00, &
    1.0471975511965977462D+00, &
    1.5707963267948966192D+00, &
    2.0000000000000000000D+00, &
    3.0000000000000000000D+00, &
    3.1415926535897932385D+00, &
    4.0000000000000000000D+00, &
    5.0000000000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cos_degree_values ( n_data, x, fx )

!*****************************************************************************80
!
!c COS_DEGREE_VALUES: the cosine function with argument in degrees.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Cos[x Degree]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.99619469809174553230D+00, &
     1.0000000000000000000D+00, &
     0.99984769515639123916D+00, &
     0.99939082701909573001D+00, &
     0.99862953475457387378D+00, &
     0.99756405025982424761D+00, &
     0.99619469809174553230D+00, &
     0.98480775301220805937D+00, &
     0.96592582628906828675D+00, &
     0.86602540378443864676D+00, &
     0.70710678118654752440D+00, &
     0.50000000000000000000D+00, &
     0.25881904510252076235D+00, &
     0.087155742747658173558D+00, &
     0.069756473744125300776D+00, &
     0.052335956242943832722D+00, &
     0.034899496702500971646D+00, &
     0.017452406437283512819D+00, &
     0.000000000000000000000D+00, &
    -0.017452406437283512819D+00, &
    -0.25881904510252076235D+00, &
    -1.0000000000000000000D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     -5.0D+00, &
      0.0D+00, &
      1.0D+00, &
      2.0D+00, &
      3.0D+00, &
      4.0D+00, &
      5.0D+00, &
     10.0D+00, &
     15.0D+00, &
     30.0D+00, &
     45.0D+00, &
     60.0D+00, &
     75.0D+00, &
     85.0D+00, &
     86.0D+00, &
     87.0D+00, &
     88.0D+00, &
     89.0D+00, &
     90.0D+00, &
     91.0D+00, &
    105.0D+00, &
    180.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cos_power_int_values ( n_data, a, b, n, fx )

!*****************************************************************************80
!
!! COS_POWER_INT_VALUES returns some values of the cosine power integral.
!
!  Discussion:
!
!    The function has the form
!
!      COS_POWER_INT(A,B,N) = integral ( A <= T <= B ) ( cos(T) )^N dt
!
!    In Mathematica, the function can be evaluated by:
!
!      Integrate [ ( Cos[x] )^n, { x, a, b } ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, integer ( kind = 4 ) N, the power.
!
!    Output, real ( kind = 8 ) FX, the function value.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     3.141592653589793D+00, &
     0.0D+00, &
     1.570796326794897D+00, &
     0.0D+00, &
     1.178097245096172D+00, &
     0.0D+00, &
     0.9817477042468104D+00, &
     0.0D+00, &
     0.8590292412159591D+00, &
     0.0D+00, &
     0.7731263170943632D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, &
     1, &
     2, &
     3, &
     4, &
     5, &
     6, &
     7, &
     8, &
     9, &
    10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    n = 0
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cosh_values ( n_data, x, fx )

!*****************************************************************************80
!
!! COSH_VALUES returns some values of the hyperbolic cosine function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Cosh[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
      74.209948524787844444D+00, &
       1.5430806348152437785D+00, &
       1.0000000000000000000D+00, &
       1.0050041680558035990D+00, &
       1.0200667556190758463D+00, &
       1.0453385141288604850D+00, &
       1.0810723718384548093D+00, &
       1.1276259652063807852D+00, &
       1.1854652182422677038D+00, &
       1.2551690056309430182D+00, &
       1.3374349463048445980D+00, &
       1.4330863854487743878D+00, &
       1.5430806348152437785D+00, &
       3.7621956910836314596D+00, &
      10.067661995777765842D+00, &
      27.308232836016486629D+00, &
      74.209948524787844444D+00, &
   11013.232920103323140D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
   -5.0D+00, &
   -1.0D+00, &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
   10.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cot_values ( n_data, x, fx )

!*****************************************************************************80
!
!! COT_VALUES returns some values of the cotangent function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Cot[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
   11.972209353628661620D+00, &
    3.7320508075688772935D+00, &
    1.8304877217124519193D+00, &
    1.7320508075688772935D+00, &
    1.0000000000000000000D+00, &
    0.64209261593433070301D+00, &
    0.57735026918962576451D+00, &
    0.26794919243112270647D+00, &
    0.00000000000000000000D+00, &
    0.13165249758739585347D+00, &
    0.065543462815238228565D+00, &
   -0.45765755436028576375D+00, &
   -7.0152525514345334694D+00, &
    0.86369115445061661395D+00, &
   -0.29581291553274554043D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
   0.083333333333333333333D+00, &
   0.26179938779914943654D+00, &
   0.50000000000000000000D+00, &
   0.52359877559829887308D+00, &
   0.78539816339744830962D+00, &
   1.0000000000000000000D+00, &
   1.0471975511965977462D+00, &
   1.3089969389957471827D+00, &
   1.5707963267948966192D+00, &
   1.4398966328953219010D+00, &
   1.5053464798451092601D+00, &
   2.0000000000000000000D+00, &
   3.0000000000000000000D+00, &
   4.0000000000000000000D+00, &
   5.0000000000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cp_values ( n_data, tc, p, cp )

!*****************************************************************************80
!
!! CP_VALUES returns some values of the specific heat at constant pressure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) CP, the specific heat at constant pressure,
!    in KJ/(kg K).
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 24

  real ( kind = 8 ) cp
  real ( kind = 8 ), save, dimension ( n_max ) :: cp_vec = (/ &
    4.228D+00, &
    2.042D+00, &
    1.975D+00, &
    2.013D+00, &
    2.040D+00, &
    2.070D+00, &
    2.135D+00, &
    2.203D+00, &
    2.378D+00, &
    2.541D+00, &
    2.792D+00, &
    2.931D+00, &
    4.226D+00, &
    4.223D+00, &
    4.202D+00, &
    4.177D+00, &
    4.130D+00, &
    4.089D+00, &
    4.053D+00, &
    4.021D+00, &
    3.909D+00, &
    3.844D+00, &
    3.786D+00, &
    2.890D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       5.0D+00, &
      10.0D+00, &
      50.0D+00, &
     100.0D+00, &
     200.0D+00, &
     300.0D+00, &
     400.0D+00, &
     500.0D+00, &
    1000.0D+00, &
    1500.0D+00, &
    2000.0D+00, &
    5000.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
       0.0D+00, &
     100.0D+00, &
     200.0D+00, &
     300.0D+00, &
     350.0D+00, &
     400.0D+00, &
     500.0D+00, &
     600.0D+00, &
     850.0D+00, &
    1100.0D+00, &
    1600.0D+00, &
    2000.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    tc = 0.0D+00
    p = 0.0D+00
    cp = 0.0D+00
  else
    tc = tc_vec(n_data)
    p = p_vec(n_data)
    cp = cp_vec(n_data)
  end if

  return
end
subroutine dawson_values ( n_data, x, fx )

!*****************************************************************************80
!
!! DAWSON_VALUES returns some values of Dawson's integral.
!
!  Discussion:
!
!    The definition of Dawson's integral is
!
!      D(X) = exp ( -X * X ) * integral ( 0 <= Y <= X ) exp ( Y * Y ) dY
!
!    Dawson's integral has a maximum at roughly
!
!      X = 0.9241388730
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi] * Exp[-x^2] * I * Erf[I*x] / 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.9933599239785286D-01, &
    0.1947510333680280D+00, &
    0.2826316650213119D+00, &
    0.3599434819348881D+00, &
    0.4244363835020223D+00, &
    0.4747632036629779D+00, &
    0.5105040575592318D+00, &
    0.5321017070563654D+00, &
    0.5407243187262987D+00, &
    0.5380795069127684D+00, &
    0.5262066799705525D+00, &
    0.5072734964077396D+00, &
    0.4833975173848241D+00, &
    0.4565072375268973D+00, &
    0.4282490710853986D+00, &
    0.3999398943230814D+00, &
    0.3725593489740788D+00, &
    0.3467727691148722D+00, &
    0.3229743193228178D+00, &
    0.3013403889237920D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine debye1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! DEBYE1_VALUES returns some values of Debye's function of order 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE1(x) = 1 / x * integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.99951182471380889183D+00, &
     0.99221462647120597836D+00, &
     0.96918395997895308324D+00, &
     0.88192715679060552968D+00, &
     0.77750463411224827642D+00, &
     0.68614531078940204342D+00, &
     0.60694728460981007205D+00, &
     0.53878956907785587703D+00, &
     0.48043521957304283829D+00, &
     0.38814802129793784501D+00, &
     0.36930802829242526815D+00, &
     0.32087619770014612104D+00, &
     0.29423996623154246701D+00, &
     0.27126046678502189985D+00, &
     0.20523930310221503723D+00, &
     0.16444346567994602563D+00, &
     0.10966194482735821276D+00, &
     0.82246701178200016086D-01, &
     0.54831135561510852445D-01, &
     0.32898681336964528729D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine debye2_values ( n_data, x, fx )

!*****************************************************************************80
!
!! DEBYE2_VALUES returns some values of Debye's function of order 2.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE2(x) = 2 / x^2 * integral ( 0 <= t <= x ) t^2 / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.99934911727904599738D+00, &
     0.98962402299599181205D+00, &
     0.95898426200345986743D+00, &
     0.84372119334725358934D+00, &
     0.70787847562782928288D+00, &
     0.59149637225671282917D+00, &
     0.49308264399053185014D+00, &
     0.41079413579749669069D+00, &
     0.34261396060786351671D+00, &
     0.24055368752127897660D+00, &
     0.22082770061202308232D+00, &
     0.17232915939014138975D+00, &
     0.14724346738730182894D+00, &
     0.12666919046715789982D+00, &
     0.74268805954862854626D-01, &
     0.47971498020121871622D-01, &
     0.21369201683658373846D-01, &
     0.12020564476446432799D-01, &
     0.53424751249537071952D-02, &
     0.19232910450553508562D-02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine debye3_values ( n_data, x, fx )

!*****************************************************************************80
!
!! DEBYE3_VALUES returns some values of Debye's function of order 3.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE3(x) = 3 / x^3 * integral ( 0 <= t <= x ) t^3 / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.99926776885985461940D+00, &
     0.98833007755734698212D+00, &
     0.95390610472023510237D+00, &
     0.82496296897623372315D+00, &
     0.67441556407781468010D+00, &
     0.54710665141286285468D+00, &
     0.44112847372762418113D+00, &
     0.35413603481042394211D+00, &
     0.28357982814342246206D+00, &
     0.18173691382177474795D+00, &
     0.16277924385112436877D+00, &
     0.11759741179993396450D+00, &
     0.95240802723158889887D-01, &
     0.77581324733763020269D-01, &
     0.36560295673194845002D-01, &
     0.19295765690345489563D-01, &
     0.57712632276188798621D-02, &
     0.24352200674805479827D-02, &
     0.72154882216335666096D-03, &
     0.15585454565440389896D-03 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine debye4_values ( n_data, x, fx )

!*****************************************************************************80
!
!! DEBYE4_VALUES returns some values of Debye's function of order 4.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE4(x) = 4 / x^4 * integral ( 0 <= t <= x ) t^4 / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.99921896192761576256D+00, &
     0.98755425280996071022D+00, &
     0.95086788606389739976D+00, &
     0.81384569172034042516D+00, &
     0.65487406888673697092D+00, &
     0.52162830964878715188D+00, &
     0.41189273671788528876D+00, &
     0.32295434858707304628D+00, &
     0.25187863642883314410D+00, &
     0.15185461258672022043D+00, &
     0.13372661145921413299D+00, &
     0.91471377664481164749D-01, &
     0.71227828197462523663D-01, &
     0.55676547822738862783D-01, &
     0.21967566525574960096D-01, &
     0.96736755602711590082D-02, &
     0.19646978158351837850D-02, &
     0.62214648623965450200D-03, &
     0.12289514092077854510D-03, &
     0.15927210319002161231D-04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine dedekind_sum_values ( n_data, p, q, n, d )

!*****************************************************************************80
!
!! DEDEKIND_SUM_VALUES returns some values of the Dedekind sum.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hans Rademacher, Emil Grosswald,
!    Dedekind Sums,
!    Mathematics Association of America, 1972,
!    LC: QA241.R2.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA
!    by 1, and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) P, Q, the arguments of the function.
!
!    Output, integer ( kind = 4 ) N, D, the numerator and denominator of the
!    function value.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 95

  integer ( kind = 4 ) d
  integer ( kind = 4 ), save, dimension ( n_max ) :: d_vec = (/ &
     1,  1, 18,  8,  5, 18, 14, 16, 27,  5, &
    22, 72, 13, 14, 90, 32, 17, 27, 38, 40, &
     1, 18,  1, 14, 27, 22, 13, 18, 17, 38, &
     1,  1,  8,  1, 14, 16,  1, 22, 13, 14, &
    32, 17, 38,  8,  1, 18,  5, 14, 27, 22, &
    13, 90,  1, 38,  1,  1, 18,  8, 18, 14, &
    16, 27, 22, 72,  1, 14, 32, 17, 27, 38, &
     1,  5, 14, 22, 13, 17, 38,  1,  1, 18, &
     8,  1, 18, 16, 27,  1, 22, 72, 13, 18, &
    32, 17, 27, 38,  8 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  0,  1,  1,  1,  5,  5,  7, 14,  3, &
    15, 55, 11, 13, 91, 35, 20, 34, 51, 57, &
     0, -1,  0,  1,  4,  5,  4,  7,  8, 21, &
     0,  0, -1,  0, -1,  1,  0,  3,  1,  3, &
     5,  5,  9,  3,  0,  1, -1,  1, -4,  3, &
    -1, 19,  0, 11,  0,  0, -1,  1, -5, -1, &
    -1,  4, -5, -1,  0,  3, -5,  1,  2, 11, &
     0,  1, -5,  5, -4,  5, -9,  0,  0,  1, &
    -1,  0,  5, -7, -4,  0, -3,  1,  4, -7, &
    -3,  1, -2,  3,  3 /)
  integer ( kind = 4 ) p
  integer ( kind = 4 ), save, dimension ( n_max ) :: p_vec = (/ &
     1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
     1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
     2,  2,  2,  2,  2,  2,  2,  2,  2,  2, &
     3,  3,  3,  3,  3,  3,  3,  3,  3,  3, &
     3,  3,  3,  3,  4,  4,  4,  4,  4,  4, &
     4,  4,  4,  4,  5,  5,  5,  5,  5,  5, &
     5,  5,  5,  5,  5,  5,  5,  5,  5,  5, &
     6,  6,  6,  6,  6,  6,  6,  7,  7,  7, &
     7,  7,  7,  7,  7,  7,  7,  7,  7,  7, &
     7,  7,  7,  7,  7 /)
  integer ( kind = 4 ) q
  integer ( kind = 4 ), save, dimension ( n_max ) :: q_vec = (/ &
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
     1,  3,  5,  7,  9, 11, 13, 15, 17, 19, &
     1,  2,  4,  5,  7,  8, 10, 11, 13, 14, &
    16, 17, 19, 20,  1,  3,  5,  7,  9, 11, &
    13, 15, 17, 19,  1,  2,  3,  4,  6,  7, &
     8,  9, 11, 12, 13, 14, 16, 17, 18, 19, &
     1,  5,  7, 11, 13, 17, 19,  1,  2,  3, &
     4,  5,  6,  8,  9, 10, 11, 12, 13, 15, &
    16, 17, 18, 19, 20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    p = 0
    q = 0
    n = 0
    d = 0
  else
    p = p_vec(n_data)
    q = q_vec(n_data)
    n = n_vec(n_data)
    d = d_vec(n_data)
  end if

  return
end
subroutine dielectric_values ( n_data, tc, p, eps )

!*****************************************************************************80
!
!! DIELECTRIC_VALUES returns some values of the static dielectric constant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher, George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) EPS, the dielectric constant, dimensionless.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) eps
  real ( kind = 8 ), save, dimension ( n_max ) :: eps_vec = (/ &
     88.29D+00, &
     90.07D+00, &
     92.02D+00, &
     95.14D+00, &
    100.77D+00, &
     78.85D+00, &
     70.27D+00, &
     62.60D+00, &
     55.78D+00, &
     44.31D+00, &
     35.11D+00, &
     20.40D+00, &
      1.17D+00, &
      1.11D+00, &
      1.08D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
     100.0D+00, &
     500.0D+00, &
    1000.0D+00, &
    2000.0D+00, &
    5000.0D+00, &
     100.0D+00, &
     100.0D+00, &
     100.0D+00, &
     100.0D+00, &
     100.0D+00, &
     100.0D+00, &
     100.0D+00, &
     100.0D+00, &
     100.0D+00, &
     100.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
     25.0D+00, &
     50.0D+00, &
     75.0D+00, &
    100.0D+00, &
    150.0D+00, &
    200.0D+00, &
    300.0D+00, &
    400.0D+00, &
    500.0D+00, &
    600.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    tc = 0.0D+00
    p = 0.0D+00
    eps = 0.0D+00
  else
    tc = tc_vec(n_data)
    p = p_vec(n_data)
    eps = eps_vec(n_data)
  end if

  return
end
subroutine dilogarithm_values ( n_data, x, fx )

!*****************************************************************************80
!
!! DILOGARITHM_VALUES returns some values of the dilogarithm function.
!
!  Discussion:
!
!    The dilogarithm is defined as
!
!      Li_2(X) = - integral ( 0 <= T <= X ) ln ( 1 - T ) / T dT
!
!    The dilogarithm is also known as Spence's integral.
!
!    In Abramowitz and Stegun form of the function is different,
!    and is equivalent to evaluated Li_2(1-X).
!
!    The dilogarithm is the special case, with N = 2, of the
!    polylogarithm Li_N(X).
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyLog[2,X]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.5063929246449603D-01, &
    0.1026177910993911D+00, &
    0.1560350339454831D+00, &
    0.2110037754397048D+00, &
    0.2676526390827326D+00, &
    0.3261295100754761D+00, &
    0.3866059411605865D+00, &
    0.4492829744712817D+00, &
    0.5143989891542119D+00, &
    0.5822405264650125D+00, &
    0.6531576315069018D+00, &
    0.7275863077163334D+00, &
    0.8060826895177240D+00, &
    0.8893776242860387D+00, &
    0.9784693929303061D+00, &
    0.1074794600008248D+01, &
    0.1180581123830255D+01, &
    0.1299714723004959D+01, &
    0.1440633796970039D+01, &
    0.1644934066848226D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.05D+00, &
    0.10D+00, &
    0.15D+00, &
    0.20D+00, &
    0.25D+00, &
    0.30D+00, &
    0.35D+00, &
    0.40D+00, &
    0.45D+00, &
    0.50D+00, &
    0.55D+00, &
    0.60D+00, &
    0.65D+00, &
    0.70D+00, &
    0.75D+00, &
    0.80D+00, &
    0.85D+00, &
    0.90D+00, &
    0.95D+00, &
    0.10D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine e1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! E1_VALUES returns some values of the exponential integral function E1(X).
!
!  Discussion:
!
!    The exponential integral E1(X) is defined by the formula:
!
!      E1(X) = integral ( 1 <= T <= +oo ) exp ( -X*T ) / T dT
!
!    In Mathematica, the function can be evaluated by:
!
!      ExpIntegralE[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5597735947761608D+00, &
    0.4543795031894021D+00, &
    0.3737688432335091D+00, &
    0.3105965785455430D+00, &
    0.2601839393259996D+00, &
    0.2193839343955203D+00, &
    0.1859909045360402D+00, &
    0.1584084368514626D+00, &
    0.1354509578491291D+00, &
    0.1162193125713579D+00, &
    0.1000195824066327D+00, &
    0.8630833369753979D-01, &
    0.7465464440125305D-01, &
    0.6471312936386886D-01, &
    0.5620437817453485D-01, &
    0.4890051070806112D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine ei_values ( n_data, x, fx )

!*****************************************************************************80
!
!! EI_VALUES returns some values of the exponential integral function EI(X).
!
!  Discussion:
!
!    The exponential integral EI(X) has the formula:
!
!      EI(X) = - integral ( -X <= T <= +oo ) exp ( -T ) / T dT
!
!    In Mathematica, the function can be evaluated by:
!
!      ExpIntegralEi[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.4542199048631736D+00, &
    0.7698812899373594D+00, &
    0.1064907194624291D+01, &
    0.1347396548212326D+01, &
    0.1622811713696867D+01, &
    0.1895117816355937D+01, &
    0.2167378279563403D+01, &
    0.2442092285192652D+01, &
    0.2721398880232024D+01, &
    0.3007207464150646D+01, &
    0.3301285449129798D+01, &
    0.3605319949019469D+01, &
    0.3920963201354904D+01, &
    0.4249867557487934D+01, &
    0.4593713686953585D+01, &
    0.4954234356001890D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine elliptic_ea_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ELLIPTIC_EA_VALUES returns values of the complete elliptic integral E(ALPHA).
!
!  Discussion:
!
!    This is one form of what is sometimes called the complete elliptic
!    integral of the second kind.
!
!    The function is defined by the formula:
!
!      E(ALPHA) = integral ( 0 <= T <= PI/2 )
!        sqrt ( 1 - sin ( ALPHA )^2 * sin ( T )^2 ) dT
!
!    In Mathematica, the function can be evaluated by:
!
!      EllipticE[(Sin[Pi*alpha/180])^2]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function, measured
!    in degrees.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 19

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.570796326794897D+00, &
    1.567809073977622D+00, &
    1.558887196601596D+00, &
    1.544150496914673D+00, &
    1.523799205259774D+00, &
    1.498114928422116D+00, &
    1.467462209339427D+00, &
    1.432290969306756D+00, &
    1.393140248523812D+00, &
    1.350643881047676D+00, &
    1.305539094297794D+00, &
    1.258679624779997D+00, &
    1.211056027568459D+00, &
    1.163827964493139D+00, &
    1.118377737969864D+00, &
    1.076405113076403D+00, &
    1.040114395706010D+00, &
    1.012663506234396D+00, &
    1.000000000000000D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.0D+00, &
     5.0D+00, &
    10.0D+00, &
    15.0D+00, &
    20.0D+00, &
    25.0D+00, &
    30.0D+00, &
    35.0D+00, &
    40.0D+00, &
    45.0D+00, &
    50.0D+00, &
    55.0D+00, &
    60.0D+00, &
    65.0D+00, &
    70.0D+00, &
    75.0D+00, &
    80.0D+00, &
    85.0D+00, &
    90.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine elliptic_em_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ELLIPTIC_EM_VALUES returns values of the complete elliptic integral E(M).
!
!  Discussion:
!
!    This is one form of what is sometimes called the complete elliptic
!    integral of the second kind.
!
!    The function is defined by the formula:
!
!      E(M) = integral ( 0 <= T <= PI/2 )
!        sqrt ( 1 - M * sin ( T )^2 ) dT
!
!    In Mathematica, the function can be evaluated by:
!
!      EllipticE[m]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.570796326794897D+00, &
    1.550973351780472D+00, &
    1.530757636897763D+00, &
    1.510121832092819D+00, &
    1.489035058095853D+00, &
    1.467462209339427D+00, &
    1.445363064412665D+00, &
    1.422691133490879D+00, &
    1.399392138897432D+00, &
    1.375401971871116D+00, &
    1.350643881047676D+00, &
    1.325024497958230D+00, &
    1.298428035046913D+00, &
    1.270707479650149D+00, &
    1.241670567945823D+00, &
    1.211056027568459D+00, &
    1.178489924327839D+00, &
    1.143395791883166D+00, &
    1.104774732704073D+00, &
    1.060473727766278D+00, &
    1.000000000000000D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.05D+00, &
    0.10D+00, &
    0.15D+00, &
    0.20D+00, &
    0.25D+00, &
    0.30D+00, &
    0.35D+00, &
    0.40D+00, &
    0.45D+00, &
    0.50D+00, &
    0.55D+00, &
    0.60D+00, &
    0.65D+00, &
    0.70D+00, &
    0.75D+00, &
    0.80D+00, &
    0.85D+00, &
    0.90D+00, &
    0.95D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine elliptic_ka_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ELLIPTIC_KA_VALUES returns values of the complete elliptic integral K(ALPHA).
!
!  Discussion:
!
!    This is one form of what is sometimes called the complete elliptic integral
!    of the first kind.
!
!    The function is defined by the formula:
!
!      K(ALPHA) = integral ( 0 <= T <= PI/2 )
!        dT / sqrt ( 1 - sin ( ALPHA )^2 * sin ( T )^2 )
!
!    In Mathematica, the function can be evaluated by:
!
!      EllipticK[(Sin[alpha*Pi/180])^2]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function, measured
!    in degrees.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1570796326794897D+01, &
    0.1573792130924768D+01, &
    0.1582842804338351D+01, &
    0.1598142002112540D+01, &
    0.1620025899124204D+01, &
    0.1648995218478530D+01, &
    0.1685750354812596D+01, &
    0.1731245175657058D+01, &
    0.1786769134885021D+01, &
    0.1854074677301372D+01, &
    0.1935581096004722D+01, &
    0.2034715312185791D+01, &
    0.2156515647499643D+01, &
    0.2308786798167196D+01, &
    0.2504550079001634D+01, &
    0.2768063145368768D+01, &
    0.3153385251887839D+01, &
    0.3831741999784146D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.0D+00, &
     5.0D+00, &
    10.0D+00, &
    15.0D+00, &
    20.0D+00, &
    25.0D+00, &
    30.0D+00, &
    35.0D+00, &
    40.0D+00, &
    45.0D+00, &
    50.0D+00, &
    55.0D+00, &
    60.0D+00, &
    65.0D+00, &
    70.0D+00, &
    75.0D+00, &
    80.0D+00, &
    85.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine elliptic_km_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ELLIPTIC_KM_VALUES returns values of the complete elliptic integral K(M).
!
!  Discussion:
!
!    This is one form of what is sometimes called the complete elliptic
!    integral of the first kind.
!
!    The function is defined by the formula:
!
!      K(M) = integral ( 0 <= T <= PI/2 )
!        dT / sqrt ( 1 - M * sin ( T )^2 )
!
!    In Mathematica, the function can be evaluated by:
!
!      EllipticK[m]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.570796326794897D+00, &
    1.591003453790792D+00, &
    1.612441348720219D+00, &
    1.635256732264580D+00, &
    1.659623598610528D+00, &
    1.685750354812596D+00, &
    1.713889448178791D+00, &
    1.744350597225613D+00, &
    1.777519371491253D+00, &
    1.813883936816983D+00, &
    1.854074677301372D+00, &
    1.898924910271554D+00, &
    1.949567749806026D+00, &
    2.007598398424376D+00, &
    2.075363135292469D+00, &
    2.156515647499643D+00, &
    2.257205326820854D+00, &
    2.389016486325580D+00, &
    2.578092113348173D+00, &
    2.908337248444552D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.00D+00, &
     0.05D+00, &
     0.10D+00, &
     0.15D+00, &
     0.20D+00, &
     0.25D+00, &
     0.30D+00, &
     0.35D+00, &
     0.40D+00, &
     0.45D+00, &
     0.50D+00, &
     0.55D+00, &
     0.60D+00, &
     0.65D+00, &
     0.70D+00, &
     0.75D+00, &
     0.80D+00, &
     0.85D+00, &
     0.90D+00, &
     0.95D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine erf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ERF_VALUES returns some values of the ERF or "error" function.
!
!  Discussion:
!
!    The error function is defined by:
!
!      ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
!
!    In Mathematica, the function can be evaluated by:
!
!      Erf[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.1124629160182849D+00, &
    0.2227025892104785D+00, &
    0.3286267594591274D+00, &
    0.4283923550466685D+00, &
    0.5204998778130465D+00, &
    0.6038560908479259D+00, &
    0.6778011938374185D+00, &
    0.7421009647076605D+00, &
    0.7969082124228321D+00, &
    0.8427007929497149D+00, &
    0.8802050695740817D+00, &
    0.9103139782296354D+00, &
    0.9340079449406524D+00, &
    0.9522851197626488D+00, &
    0.9661051464753107D+00, &
    0.9763483833446440D+00, &
    0.9837904585907746D+00, &
    0.9890905016357307D+00, &
    0.9927904292352575D+00, &
    0.9953222650189527D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine erfc_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ERFC_VALUES returns some values of the complementary error function.
!
!  Discussion:
!
!    The complementary error function is defined by:
!
!      ERFC(X) = 1 -
!        ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
!
!    In Mathematica, the function can be evaluated by:
!
!      Erfc[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.000000000000000D+00, &
    0.7772974107895215D+00, &
    0.5716076449533315D+00, &
    0.3961439091520741D+00, &
    0.2578990352923395D+00, &
    0.1572992070502851D+00, &
    0.08968602177036462D+00, &
    0.04771488023735119D+00, &
    0.02365161665535599D+00, &
    0.01090949836426929D+00, &
    0.004677734981047266D+00, &
    0.001862846297981891D+00, &
    0.0006885138966450786D+00, &
    0.0002360344165293492D+00, &
    0.00007501319466545902D+00, &
    0.00002209049699858544D+00, &
    6.025761151762095D-06, &
    1.521993362862285D-06, &
    3.558629930076853D-07, &
    7.700392745696413D-08, &
    1.541725790028002D-08 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00, &
    3.2D+00, &
    3.4D+00, &
    3.6D+00, &
    3.8D+00, &
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine euler_number_values ( n_data, n, c )

!*****************************************************************************80
!
!! EULER_NUMBER_VALUES returns some values of the Euler numbers.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      EulerE[n]
!
!    These numbers rapidly get too big to store in an ordinary integer!
!
!    The terms of odd index are 0.
!
!    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
!
!  First terms:
!
!    E0  = 1
!    E1  = 0
!    E2  = -1
!    E3  = 0
!    E4  = 5
!    E5  = 0
!    E6  = -61
!    E7  = 0
!    E8  = 1385
!    E9  = 0
!    E10 = -50521
!    E11 = 0
!    E12 = 2702765
!    E13 = 0
!    E14 = -199360981
!    E15 = 0
!    E16 = 19391512145
!    E17 = 0
!    E18 = -2404879675441
!    E19 = 0
!    E20 = 370371188237525
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the Euler number.
!
!    Output, integer ( kind = 4 ) C, the value of the Euler number.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 8

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
    1, 0, -1, 5, 61, 1385, -50521, 2702765 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, 1, 2, 4, 6, 8, 10, 12 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine euler_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! EULER_POLY_VALUES returns some values of the Euler polynomials.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      EulerE[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the Euler polynomial.
!
!    Output, real ( kind = 8 ) X, the argument of the Euler polynomial.
!
!    Output, real ( kind = 8 ) FX, the value of the Euler polynomial.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 27

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.100000000000D+01, &
    -0.300000000000D+00, &
    -0.160000000000D+00, &
     0.198000000000D+00, &
     0.185600000000D+00, &
    -0.403680000000D+00, &
    -0.560896000000D+00, &
     0.171878880000D+01, &
     0.318043136000D+01, &
    -0.125394670080D+02, &
    -0.289999384576D+02, &
    -0.625000000000D-01, &
    -0.174240000000D+00, &
    -0.297680000000D+00, &
    -0.404320000000D+00, &
    -0.475260000000D+00, &
    -0.500000000000D+00, &
    -0.475240000000D+00, &
    -0.403680000000D+00, &
    -0.292820000000D+00, &
    -0.153760000000D+00, &
     0.000000000000D+00, &
     0.153760000000D+00, &
     0.292820000000D+00, &
     0.403680000000D+00, &
     0.475240000000D+00, &
     0.500000000000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, &
     1, &
     2, &
     3, &
     4, &
     5, &
     6, &
     7, &
     8, &
     9, &
    10, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
     0.2D+00, &
    -0.5D+00, &
    -0.4D+00, &
    -0.3D+00, &
    -0.2D+00, &
    -0.1D+00, &
     0.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.3D+00, &
     0.4D+00, &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine exp_values ( n_data, x, fx )

!*****************************************************************************80
!
!! EXP_VALUES returns some values of the exponential function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Exp[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 24

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.000045399929762484851536D+00, &
    0.0067379469990854670966D+00, &
    0.36787944117144232160D+00, &
    1.0000000000000000000D+00, &
    1.0000000100000000500D+00, &
    1.0001000050001666708D+00, &
    1.0010005001667083417D+00, &
    1.0100501670841680575D+00, &
    1.1051709180756476248D+00, &
    1.2214027581601698339D+00, &
    1.3498588075760031040D+00, &
    1.4918246976412703178D+00, &
    1.6487212707001281468D+00, &
    1.8221188003905089749D+00, &
    2.0137527074704765216D+00, &
    2.2255409284924676046D+00, &
    2.4596031111569496638D+00, &
    2.7182818284590452354D+00, &
    7.3890560989306502272D+00, &
    23.140692632779269006D+00, &
    148.41315910257660342D+00, &
    22026.465794806716517D+00, &
    4.8516519540979027797D+08, &
    2.3538526683701998541D+17 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     -10.0D+00, &
      -5.0D+00, &
      -1.0D+00, &
       0.0D+00, &
       0.00000001D+00, &
       0.0001D+00, &
       0.001D+00, &
       0.01D+00, &
       0.1D+00, &
       0.2D+00, &
       0.3D+00, &
       0.4D+00, &
       0.5D+00, &
       0.6D+00, &
       0.7D+00, &
       0.8D+00, &
       0.9D+00, &
       1.0D+00, &
       2.0D+00, &
       3.1415926535897932385D+00, &
       5.0D+00, &
      10.0D+00, &
      20.0D+00, &
      40.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine exp3_int_values ( n_data, x, fx )

!*****************************************************************************80
!
!! EXP3_INT_VALUES returns some values of the EXP3 integral function.
!
!  Discussion:
!
!    The function is defined by:
!
!      EXP3_INT(x) = integral ( 0 <= t <= x ) exp ( -t^3 ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.19531249963620212007D-02, &
     0.78124990686775522671D-02, &
     0.31249761583499728667D-01, &
     0.12493899888803079984D+00, &
     0.48491714311363971332D+00, &
     0.80751118213967145286D+00, &
     0.86889265412623270696D+00, &
     0.88861722235357162648D+00, &
     0.89286018500218176869D+00, &
     0.89295351429387631138D+00, &
     0.89297479112737843939D+00, &
     0.89297880579798112220D+00, &
     0.89297950317496621294D+00, &
     0.89297951152951902903D+00, &
     0.89297951156918122102D+00, &
     0.89297951156924734716D+00, &
     0.89297951156924917298D+00, &
     0.89297951156924921121D+00, &
     0.89297951156924921122D+00, &
     0.89297951156924921122D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.0019531250D+00, &
     0.0078125000D+00, &
     0.0312500000D+00, &
     0.1250000000D+00, &
     0.5000000000D+00, &
     1.0000000000D+00, &
     1.2500000000D+00, &
     1.5000000000D+00, &
     1.8750000000D+00, &
     2.0000000000D+00, &
     2.1250000000D+00, &
     2.2500000000D+00, &
     2.5000000000D+00, &
     2.7500000000D+00, &
     3.0000000000D+00, &
     3.1250000000D+00, &
     3.2500000000D+00, &
     3.5000000000D+00, &
     3.7500000000D+00, &
     4.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine exponential_cdf_values ( n_data, lambda, x, fx )

!*****************************************************************************80
!
!! EXPONENTIAL_CDF_VALUES returns some values of the Exponential CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = ExponentialDistribution [ lambda ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) LAMBDA, the parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 9

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.3934693402873666D+00, &
    0.6321205588285577D+00, &
    0.7768698398515702D+00, &
    0.8646647167633873D+00, &
    0.8646647167633873D+00, &
    0.9816843611112658D+00, &
    0.9975212478233336D+00, &
    0.9996645373720975D+00, &
    0.9999546000702375D+00 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine extreme_values_cdf_values ( n_data, alpha, beta, x, fx )

!*****************************************************************************80
!
!! EXTREME_VALUES_CDF_VALUES returns some values of the Extreme Values CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = ExtremeValuesDistribution [ alpha, beta ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) ALPHA, the first parameter of the distribution.
!
!    Output, real ( kind = 8 ) BETA, the second parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) alpha
  real ( kind = 8 ), save, dimension ( n_max ) :: alpha_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  real ( kind = 8 ) beta
  real ( kind = 8 ), save, dimension ( n_max ) :: beta_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.3678794411714423D+00, &
    0.8734230184931166D+00, &
    0.9818510730616665D+00, &
    0.9975243173927525D+00, &
    0.5452392118926051D+00, &
    0.4884435800065159D+00, &
    0.4589560693076638D+00, &
    0.4409910259429826D+00, &
    0.5452392118926051D+00, &
    0.3678794411714423D+00, &
    0.1922956455479649D+00, &
    0.6598803584531254D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    alpha = 0.0D+00
    beta = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    alpha = alpha_vec(n_data)
    beta = beta_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine f_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! F_CDF_VALUES returns some values of the F CDF test function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = FRatioDistribution [ dfn, dfd ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) A, integer B, the parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) a
  integer ( kind = 4 ), save, dimension ( n_max ) :: a_vec = (/ &
    1, &
    1, &
    5, &
    1, &
    2, &
    4, &
    1, &
    6, &
    8, &
    1, &
    3, &
    6, &
    1, &
    1, &
    1, &
    1, &
    2, &
    3, &
    4, &
    5 /)
  integer ( kind = 4 ) b
  integer ( kind = 4 ), save, dimension ( n_max ) :: b_vec = (/ &
     1, &
     5, &
     1, &
     5, &
    10, &
    20, &
     5, &
     6, &
    16, &
     5, &
    10, &
    12, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.4999714850534485D+00, &
    0.4996034370170990D+00, &
    0.7496993658293228D+00, &
    0.7504656462757382D+00, &
    0.7514156325324275D+00, &
    0.8999867031372156D+00, &
    0.8997127554259699D+00, &
    0.9002845660853669D+00, &
    0.9500248817817622D+00, &
    0.9500574946122442D+00, &
    0.9501926400000000D+00, &
    0.9750133887312993D+00, &
    0.9900022327445249D+00, &
    0.9949977837872073D+00, &
    0.9989999621122122D+00, &
    0.5687988496283079D+00, &
    0.5351452100063650D+00, &
    0.5143428032407864D+00, &
    0.5000000000000000D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     1.00D+00, &
     0.528D+00, &
     1.89D+00, &
     1.69D+00, &
     1.60D+00, &
     1.47D+00, &
     4.06D+00, &
     3.05D+00, &
     2.09D+00, &
     6.61D+00, &
     3.71D+00, &
     3.00D+00, &
    10.01D+00, &
    16.26D+00, &
    22.78D+00, &
    47.18D+00, &
     1.00D+00, &
     1.00D+00, &
     1.00D+00, &
     1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0
    b = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine f_noncentral_cdf_values ( n_data, n1, n2, lambda, x, fx )

!*****************************************************************************80
!
!! F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NoncentralFRatioDistribution [ n1, n2, lambda ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N1, integer N2, the numerator and denominator
!    degrees of freedom.
!
!    Output, real ( kind = 8 ) LAMBDA, the noncentrality parameter.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.6367825323508774D+00, &
    0.5840916116305482D+00, &
    0.3234431872392788D+00, &
    0.4501187879813550D+00, &
    0.6078881441188312D+00, &
    0.7059275551414605D+00, &
    0.7721782003263727D+00, &
    0.8191049017635072D+00, &
    0.3170348430749965D+00, &
    0.4327218008454471D+00, &
    0.4502696915707327D+00, &
    0.4261881186594096D+00, &
    0.6753687206341544D+00, &
    0.4229108778879005D+00, &
    0.6927667261228938D+00, &
    0.3632174676491226D+00, &
    0.4210054012695865D+00, &
    0.4266672258818927D+00, &
    0.4464016600524644D+00, &
    0.8445888579504827D+00, &
    0.4339300273343604D+00 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
    0.00D+00, &
    0.00D+00, &
    0.25D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    1.00D+00, &
    1.00D+00, &
    0.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n1
  integer ( kind = 4 ), save, dimension ( n_max ) :: n1_vec = (/ &
     1,  1,  1,  1, &
     1,  1,  1,  1, &
     1,  1,  2,  2, &
     3,  3,  4,  4, &
     5,  5,  6,  6, &
     8, 16 /)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ), save, dimension ( n_max ) :: n2_vec = (/ &
     1,  5,  5,  5, &
     5,  5,  5,  5, &
     5,  5,  5, 10, &
     5,  5,  5,  5, &
     1,  5,  6, 12, &
    16,  8 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    0.50D+00, &
    1.00D+00, &
    2.00D+00, &
    3.00D+00, &
    4.00D+00, &
    5.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    2.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n1 = 0
    n2 = 0
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n1 = n1_vec(n_data)
    n2 = n2_vec(n_data)
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine factorial_rising_values ( n_data, m, n, fmn )

!*****************************************************************************80
!
!! FACTORIAL_RISING_VALUES returns values of the integer Pochhammer function.
!
!  Discussion:
!
!    The integer Pochhammer function is sometimes symbolized by (m)_n.
!
!    The definition of the integer Pochhammer function is
!
!      (m)_n = (m-1+n)! / (m-1)!
!            = ( m ) * ( m + 1 ) * ( m + 2 ) ... * ( m - 1 + n )
!            = Gamma ( m + n ) / Gamma ( m )
!
!    We assume 0 <= N <= M.
!
!    In Mathematica, the function can be evaluated by:
!
!      Pochhammer[m,n]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) M, N, the arguments of the function.
!
!    Output, integer ( kind = 4 ) FMN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 8

  integer ( kind = 4 ) fmn
  integer ( kind = 4 ), save, dimension ( n_max ) :: fmn_vec = (/ &
     1, 10, 4000, 110, 6840, 840, 970200, 5040 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    50, 10, 4000, 10, 18, 4, 98, 1 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,   1,   2,  3, 4,  3, 7 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    m = 0
    n = 0
    fmn = 0
  else
    m = m_vec(n_data)
    n = n_vec(n_data)
    fmn = fmn_vec(n_data)
  end if

  return
end
subroutine factorial_values ( n_data, n, fn )

!*****************************************************************************80
!
!! FACTORIAL_VALUES returns values of the factorial function.
!
!  Discussion:
!
!    0! = 1
!    I! = Product ( 1 <= J <= I ) I
!
!    In Mathematica, the function can be evaluated by:
!
!      n!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, integer ( kind = 4 ) FN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  integer ( kind = 4 ), save, dimension ( n_max ) :: fn_vec = (/ &
            1, &
            1, &
            2, &
            6, &
           24, &
          120, &
          720, &
         5040, &
        40320, &
       362880, &
      3628800, &
     39916800, &
    479001600 /)
  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2,  3, &
     4,  5,  6,  7, &
     8,  9, 10, 11, &
    12 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    fn = 0
  else
    n = n_vec(n_data)
    fn = fn_vec(n_data)
  end if

  return
end
subroutine factorial2_values ( n_data, n, fn )

!*****************************************************************************80
!
!! FACTORIAL2_VALUES returns values of the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!    In Mathematica, the function can be evaluated by:
!
!      n!!
!
!  Example:
!
!     N    N!!
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, integer ( kind = 4 ) FN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  integer ( kind = 4 ), save, dimension ( n_max ) :: fn_vec = (/ &
          1, &
          1, &
          2, &
          3, &
          8, &
         15, &
         48, &
        105, &
        384, &
        945, &
       3840, &
      10395, &
      46080, &
     135135, &
     645120, &
    2027025 /)
  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, &
     1,  2,  3,  4,  5, &
     6,  7,  8,  9, 10, &
    11, 12, 13, 14, 15 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    fn = 0
  else
    n = n_vec(n_data)
    fn = fn_vec(n_data)
  end if

  return
end
subroutine fresnel_cos_values ( n_data, x, fx )

!*****************************************************************************80
!
!! FRESNEL_COS_VALUES returns values of the Fresnel cosine integral function.
!
!  Discussion:
!
!    The Fresnel cosine integral is defined by:
!
!      C(X) = integral ( 0 <= T <= X ) cos ( PI * T^2 / 2 ) / T dT
!
!    In Mathematica, the function can be evaluated by:
!
!      FresnelC[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.1999210575944531D+00, &
    0.3974807591723594D+00, &
    0.5810954469916523D+00, &
    0.7228441718963561D+00, &
    0.7798934003768228D+00, &
    0.7154377229230734D+00, &
    0.5430957835462564D+00, &
    0.3654616834404877D+00, &
    0.3336329272215571D+00, &
    0.4882534060753408D+00, &
    0.6362860449033195D+00, &
    0.5549614058564281D+00, &
    0.3889374961919690D+00, &
    0.4674916516989059D+00, &
    0.6057207892976856D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine fresnel_sin_values ( n_data, x, fx )

!*****************************************************************************80
!
!! FRESNEL_SIN_VALUES returns some values of the Fresnel sine integral function.
!
!  Discussion:
!
!    The Fresnel sine integral is defined by
!
!      S(X) = integral ( 0 <= T <= X ) sin ( pi * T^2 / 2 ) / T dT
!
!    In Mathematica, the function can be evaluated by:
!
!      FresnelS[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.4187609161656762D-02, &
    0.3335943266061318D-01, &
    0.1105402073593870D+00, &
    0.2493413930539178D+00, &
    0.4382591473903548D+00, &
    0.6234009185462497D+00, &
    0.7135250773634121D+00, &
    0.6388876835093809D+00, &
    0.4509387692675831D+00, &
    0.3434156783636982D+00, &
    0.4557046121246569D+00, &
    0.6196899649456836D+00, &
    0.5499893231527195D+00, &
    0.3915284435431718D+00, &
    0.4963129989673750D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine frobenius_number_data_values ( n_data, order, c, f )

!*****************************************************************************80
!
!! FROBENIUS_NUMBER_DATA_VALUES returns data for the Frobenius problem.
!
!  Discussion:
!
!    The user should first call FROBENIUS_NUMBER_ORDER_VALUES to get the
!    order or size of the "C" vector that will be returned by this routine.
!
!    The Frobenius number of order N and data C is the solution of the
!    Frobenius coin sum problem for N coin denominations C(1) through C(N).
!
!    The Frobenius coin sum problem assumes the existence of
!    N coin denominations, and asks for the largest value that cannot
!    be formed by any combination of coins of these denominations.
!
!    The coin denominations are assumed to be distinct positive integers.
!
!    For general N, this problem is fairly difficult to handle.
!
!    For N = 2, it is known that:
!
!    * if C1 and C2 are not relatively prime, then
!      there are infinitely large values that cannot be formed.
!
!    * otherwise, the largest value that cannot be formed is
!      C1 * C2 - C1 - C2, and that exactly half the values between
!      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
!
!    As a simple example, if C1 = 2 and C2 = 7, then the largest
!    unrepresentable value is 5, and there are (5+1)/2 = 3
!    unrepresentable values, namely 1, 3, and 5.
!
!    For a general N, and a set of coin denominations C1, C2, ..., CN,
!    the Frobenius number F(N, C(1:N) ) is defined as the largest value
!    B for which the equation
!
!      C1*X1 + C2*X2 + ... + CN*XN = B
!
!    has no nonnegative integer solution X(1:N).
!
!    In Mathematica, the Frobenius number can be determined by the command:
!
!      FrobeniusNumber[ {C1,...,CN} ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerard Cornuejols, Regina Urbaniak, Robert Weismantel, Laurence Wolsey,
!    Decomposition of Integer Programs and of Generating Sets,
!    in Algorithms - ESA '97,
!    Lecture Notes in Computer Science 1284,
!    edited by R Burkard, G Woeginger,
!    Springer, 1997, pages 92-103.
!
!    Robert Owens,
!    An Algorithm to Solve the Frobenius Problem,
!    Mathematics Magazine,
!    Volume 76, Number 4, October 2003, 264-275.
!
!    James Sylvester,
!    Question 7382,
!    Mathematical Questions with their Solutions,
!    Educational Times,
!    Volume 41, page 21, 1884.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_DATA.  Unlike most other routines in this
!    library, this routine assumes that N_DATA has already been processed by a
!    call to FROBENIUS_NUMBER_ORDER_VALUE.  Therefore, this routine will return
!    the next set of data as long as N_DATA is in the expected range.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the problem.
!
!    Output, integer ( kind = 4 ) C(ORDER), the denominations of the problem.
!
!    Output, integer ( kind = 4 ) F, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: cvec_max = 77
  integer ( kind = 4 ), parameter :: n_max = 19
  integer ( kind = 4 ) order

  integer ( kind = 4 ) c(order)
  integer ( kind = 4 ), save, dimension ( cvec_max ) :: c_vec = (/ &
        2,     5, &
        3,    17, &
        4,    19, &
        5,    13, &
       12,    11, &
       99,   100, &
        6,     9,    20, &
        5,    17,    23, &
      137,   251,   256, &
       31,    41,    47,    61, &
      271,   277,   281,   283, &
       10,    18,    26,    33,    35, &
       34,    37,    38,    40,    43, &
    12223, 12224, 36674, 61119, 85569, &
     1000,  1476,  3764,  4864,  4871,  7773, &
    12228, 36679, 36682, 46908, 61139, 73365, &
    12137, 36405, 24269, 36407, 84545, 60683, &
    13211, 13212, 39638, 66060, 52864, 79268, 92482, &
    13429, 26850, 26855, 40280, 40281, 53711, 53714, 67141 &
  /)
  integer ( kind = 4 ) f
  integer ( kind = 4 ), save, dimension ( n_max ) :: f_vec = (/ &
          3, &
         31, &
         53, &
         47, &
        109, &
       9701, &
         43, &
         41, &
       4948, &
        240, &
      13022, &
         67, &
        175, &
   89643481, &
      47350, &
   89716838, &
   58925134, &
  104723595, &
   45094583 /)
  integer ( kind = 4 ), save :: v_data = 0
  integer ( kind = 4 ) n_data

  if ( n_data < 1 .or. n_max < n_data ) then
    n_data = 0
    v_data = 0
    c(1:order) = 0
    f = 0
  else
    c(1:order) = c_vec(v_data+1:v_data+order)
    v_data = v_data + order
    if ( n_data == n_max ) then
      v_data = 0
    end if
    f = f_vec(n_data)
  end if

  return
end
subroutine frobenius_number_order_values ( n_data, order )

!*****************************************************************************80
!
!! FROBENIUS_NUMBER_ORDER_VALUES returns orders of the Frobenius problem.
!
!  Discussion:
!
!    This routine returns the order or size of a Frobenius problem.
!    To get the corresponding data, call FROBENIUS_NUMBER_DATA_VALUES.
!
!    The Frobenius number of order N and data C is the solution of a Frobenius
!    coin sum problem for N coin denominations C(1) through C(N).
!
!    The Frobenius coin sum problem assumes the existence of
!    N coin denominations, and asks for the largest value that cannot
!    be formed by any combination of coins of these denominations.
!
!    The coin denominations are assumed to be distinct positive integers.
!
!    For general order N, this problem is fairly difficult to handle.
!
!    For order N = 2, it is known that:
!
!    * if C1 and C2 are not relatively prime, then
!      there are infinitely large values that cannot be formed.
!
!    * otherwise, the largest value that cannot be formed is
!      C1 * C2 - C1 - C2, and that exactly half the values between
!      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
!
!    As a simple example, if C1 = 2 and C2 = 7, then the largest
!    unrepresentable value is 5, and there are (5+1)/2 = 3
!    unrepresentable values, namely 1, 3, and 5.
!
!    For a general N, and a set of coin denominations C1, C2, ..., CN,
!    the Frobenius number F(N, C(1:N) ) is defined as the largest value
!    B for which the equation
!
!      C1*X1 + C2*X2 + ... + CN*XN = B
!
!    has no nonnegative integer solution X(1:N).
!
!    In Mathematica, the Frobenius number can be determined by the command:
!
!      FrobeniusNumber[ {C1,...,CN} ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerard Cornuejols, Regina Urbaniak, Robert Weismantel, Laurence Wolsey,
!    Decomposition of Integer Programs and of Generating Sets,
!    in Algorithms - ESA '97,
!    Lecture Notes in Computer Science 1284,
!    edited by R Burkard, G Woeginger,
!    Springer, 1997, pages 92-103,
!    LC: QA76.9.A43.E83.
!
!    Robert Owens,
!    An Algorithm to Solve the Frobenius Problem,
!    Mathematics Magazine,
!    Volume 76, Number 4, October 2003, 264-275.
!
!    James Sylvester,
!    Question 7382,
!    Mathematical Questions with their Solutions,
!    Educational Times,
!    Volume 41, page 21, 1884.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) ORDER, the order of a Frobenius problem.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 19

  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) order
  integer ( kind = 4 ), save, dimension ( n_max ) :: order_vec = (/ &
    2, &
    2, &
    2, &
    2, &
    2, &
    2, &
    3, &
    3, &
    3, &
    4, &
    4, &
    5, &
    5, &
    5, &
    6, &
    6, &
    6, &
    7, &
    8  &
  /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    order = 0
  else
    order = order_vec(n_data)
  end if

  return
end
subroutine frobenius_number_order2_values ( n_data, c1, c2, f )

!*****************************************************************************80
!
!! FROBENIUS_NUMBER_ORDER2_VALUES: some values of the order 2 Frobenius number.
!
!  Discussion:
!
!    The Frobenius number of order N is the solution of the Frobenius
!    coin sum problem for N coin denominations.
!
!    The Frobenius coin sum problem assumes the existence of
!    N coin denominations, and asks for the largest value that cannot
!    be formed by any combination of coins of these denominations.
!
!    The coin denominations are assumed to be distinct positive integers.
!
!    For general N, this problem is fairly difficult to handle.
!
!    For N = 2, it is known that:
!
!    * if C1 and C2 are not relatively prime, then
!      there are infinitely large values that cannot be formed.
!
!    * otherwise, the largest value that cannot be formed is
!      C1 * C2 - C1 - C2, and that exactly half the values between
!      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
!
!    As a simple example, if C1 = 2 and C2 = 7, then the largest
!    unrepresentable value is 5, and there are (5+1)/2 = 3
!    unrepresentable values, namely 1, 3, and 5.
!
!    For a general N, and a set of coin denominations C1, C2, ..., CN,
!    the Frobenius number F(N, C(1:N) ) is defined as the largest value
!    B for which the equation
!
!      C1*X1 + C2*X2 + ... + CN*XN = B
!
!    has no nonnegative integer solution X(1:N).
!
!    In Mathematica, the Frobenius number can be determined by the command:
!
!      FrobeniusNumber[ {C1,...,CN} ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerard Cornuejols, Regina Urbaniak, Robert Weismantel, Laurence Wolsey,
!    Decomposition of Integer Programs and of Generating Sets,
!    in Algorithms - ESA '97,
!    Lecture Notes in Computer Science 1284,
!    edited by R Burkard, G Woeginger,
!    Springer, 1997, pages 92-103,
!    LC: QA76.9.A43.E83.
!
!    Robert Owens,
!    An Algorithm to Solve the Frobenius Problem,
!    Mathematics Magazine,
!    Volume 76, Number 4, October 2003, 264-275.
!
!    James Sylvester,
!    Question 7382,
!    Mathematical Questions with their Solutions,
!    Educational Times,
!    Volume 41, page 21, 1884.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) C1, C2, the parameters of the function.
!
!    Output, integer ( kind = 4 ) F, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 6

  integer ( kind = 4 ) c1
  integer ( kind = 4 ), save, dimension ( n_max ) :: c1_vec = (/ &
     2, &
     3, &
     4, &
     5, &
    12, &
    99 /)
  integer ( kind = 4 ) c2
  integer ( kind = 4 ), save, dimension ( n_max ) :: c2_vec = (/ &
      5, &
     17, &
     19, &
     13, &
     11, &
    100 /)
  integer ( kind = 4 ) f
  integer ( kind = 4 ), save, dimension ( n_max ) :: f_vec = (/ &
    3, &
    31, &
    53, &
    47, &
    109, &
    9701 /)
  integer ( kind = 4 ) n_data

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    c1 = 0
    c2 = 0
    f = 0
  else
    c1 = c1_vec(n_data)
    c2 = c2_vec(n_data)
    f = f_vec(n_data)
  end if

  return
end
subroutine gamma_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_VALUES returns some values of the Gamma function.
!
!  Discussion:
!
!    The Gamma function is defined as:
!
!      Gamma(Z) = integral ( 0 <= T < +oo) T^(Z-1) exp(-T) dT
!
!    It satisfies the recursion:
!
!      Gamma(X+1) = X * Gamma(X)
!
!    Gamma is undefined for nonpositive integral X.
!    Gamma(0.5) = sqrt(PI)
!    For N a positive integer, Gamma(N+1) = N!, the standard factorial.
!
!    In Mathematica, the function can be evaluated by:
!
!      Gamma[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.3544907701811032D+01, &
    -0.1005871979644108D+03, &
     0.9943258511915060D+02, &
     0.9513507698668732D+01, &
     0.4590843711998803D+01, &
     0.2218159543757688D+01, &
     0.1772453850905516D+01, &
     0.1489192248812817D+01, &
     0.1164229713725303D+01, &
     0.1000000000000000D+01, &
     0.9513507698668732D+00, &
     0.9181687423997606D+00, &
     0.8974706963062772D+00, &
     0.8872638175030753D+00, &
     0.8862269254527580D+00, &
     0.8935153492876903D+00, &
     0.9086387328532904D+00, &
     0.9313837709802427D+00, &
     0.9617658319073874D+00, &
     0.1000000000000000D+01, &
     0.2000000000000000D+01, &
     0.6000000000000000D+01, &
     0.3628800000000000D+06, &
     0.1216451004088320D+18, &
     0.8841761993739702D+31 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.50D+00, &
    -0.01D+00, &
     0.01D+00, &
     0.10D+00, &
     0.20D+00, &
     0.40D+00, &
     0.50D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_cdf_values ( n_data, mu, sigma, x, fx )

!*****************************************************************************80
!
!! GAMMA_CDF_VALUES returns some values of the Gamma CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = GammaDistribution [ mu, sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) SIGMA, the variance of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8646647167633873D+00, &
    0.9816843611112658D+00, &
    0.9975212478233336D+00, &
    0.9996645373720975D+00, &
    0.6321205588285577D+00, &
    0.4865828809674080D+00, &
    0.3934693402873666D+00, &
    0.3296799539643607D+00, &
    0.4421745996289254D+00, &
    0.1911531694619419D+00, &
    0.6564245437845009D-01, &
    0.1857593622214067D-01 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_inc_p_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! GAMMA_INC_P_VALUES: values of normalized incomplete Gamma function P(A,X).
!
!  Discussion:
!
!    The (normalized) incomplete Gamma function is defined as:
!
!      P(A,X) = 1/Gamma(A) * integral ( 0 <= T <= X ) T^(A-1) * exp(-T) dT.
!
!    With this definition, for all A and X,
!
!      0 <= P(A,X) <= 1
!
!    and
!
!      P(A,oo) = 1.0
!
!    In Mathematica, the function can be evaluated by:
!
!      1 - GammaRegularized[A,X]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+01, &
    0.10D+01, &
    0.10D+01, &
    0.11D+01, &
    0.11D+01, &
    0.11D+01, &
    0.20D+01, &
    0.20D+01, &
    0.20D+01, &
    0.60D+01, &
    0.60D+01, &
    0.11D+02, &
    0.26D+02, &
    0.41D+02 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7382350532339351D+00, &
    0.9083579897300343D+00, &
    0.9886559833621947D+00, &
    0.3014646416966613D+00, &
    0.7793286380801532D+00, &
    0.9918490284064973D+00, &
    0.9516258196404043D-01, &
    0.6321205588285577D+00, &
    0.9932620530009145D+00, &
    0.7205974576054322D-01, &
    0.5891809618706485D+00, &
    0.9915368159845525D+00, &
    0.1018582711118352D-01, &
    0.4421745996289254D+00, &
    0.9927049442755639D+00, &
    0.4202103819530612D-01, &
    0.9796589705830716D+00, &
    0.9226039842296429D+00, &
    0.4470785799755852D+00, &
    0.7444549220718699D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.30D-01, &
    0.30D+00, &
    0.15D+01, &
    0.75D-01, &
    0.75D+00, &
    0.35D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.15D+00, &
    0.15D+01, &
    0.70D+01, &
    0.25D+01, &
    0.12D+02, &
    0.16D+02, &
    0.25D+02, &
    0.45D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_inc_q_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! GAMMA_INC_Q_VALUES: values the normalized incomplete Gamma function Q(A,X).
!
!  Discussion:
!
!    The (normalized) incomplete Gamma function is defined as:
!
!      Q(A,X) = 1/Gamma(A) * integral ( X <= T < oo ) T^(A-1) * exp(-T) dT.
!
!    With this definition, for all A and X,
!
!      0 <= Q(A,X) <= 1
!
!    and
!
!      Q(A,oo) = 0.0
!
!    In Mathematica, the function can be evaluated by:
!
!      GammaRegularized[A,X]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+01, &
    0.10D+01, &
    0.10D+01, &
    0.11D+01, &
    0.11D+01, &
    0.11D+01, &
    0.20D+01, &
    0.20D+01, &
    0.20D+01, &
    0.60D+01, &
    0.60D+01, &
    0.11D+02, &
    0.26D+02, &
    0.41D+02 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.2617649467660649D+00, &
    0.09164201026996572D+00, &
    0.01134401663780527D+00, &
    0.6985353583033387D+00, &
    0.2206713619198468D+00, &
    0.008150971593502700D+00, &
    0.9048374180359596D+00, &
    0.3678794411714423D+00, &
    0.006737946999085467D+00, &
    0.9279402542394568D+00, &
    0.4108190381293515D+00, &
    0.008463184015447498D+00, &
    0.9898141728888165D+00, &
    0.5578254003710746D+00, &
    0.007295055724436130D+00, &
    0.9579789618046939D+00, &
    0.02034102941692837D+00, &
    0.07739601577035708D+00, &
    0.5529214200244148D+00, &
    0.2555450779281301D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.30D-01, &
    0.30D+00, &
    0.15D+01, &
    0.75D-01, &
    0.75D+00, &
    0.35D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.15D+00, &
    0.15D+01, &
    0.70D+01, &
    0.25D+01, &
    0.12D+02, &
    0.16D+02, &
    0.25D+02, &
    0.45D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_inc_tricomi_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! GAMMA_INC_TRICOMI_VALUES: Tricomi's incomplete Gamma function.
!
!  Discussion:
!
!    Tricomi's incomplete Gamma function is defined as:
!
!      1/Gamma(A) * 1/X^A * integral ( 0 <= T <= X ) T^(A-1) * exp(-T) dT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+01, &
    0.10D+01, &
    0.10D+01, &
    0.11D+01, &
    0.11D+01, &
    0.11D+01, &
    0.20D+01, &
    0.20D+01, &
    0.20D+01, &
    0.60D+01, &
    0.60D+01, &
    0.11D+02, &
    0.26D+02, &
    0.41D+02 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.048292641463504D+00, &
    1.024577737369574D+00, &
    0.9493712443185374D+00, &
    1.100793230316492D+00, &
    0.8998911979655218D+00, &
    0.5301656062431039D+00, &
    0.9516258196404043D+00, &
    0.6321205588285577D+00, &
    0.1986524106001829D+00, &
    0.9071784510537487D+00, &
    0.5891809618706485D+00, &
    0.1688269752193589D+00, &
    0.4527034271637121D+00, &
    0.1965220442795224D+00, &
    0.02025928457705232D+00, &
    0.0001721181724479739D+00, &
    3.280858070850586D-07, &
    5.244396471821590D-14, &
    2.013462926183376D-37, &
    1.230623887499875D-68 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.30D-01, &
    0.30D+00, &
    0.15D+01, &
    0.75D-01, &
    0.75D+00, &
    0.35D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.15D+00, &
    0.15D+01, &
    0.70D+01, &
    0.25D+01, &
    0.12D+02, &
    0.16D+02, &
    0.25D+02, &
    0.45D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_inc_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! GAMMA_INC_VALUES: values of the incomplete Gamma function.
!
!  Discussion:
!
!    The incomplete Gamma function is defined as:
!
!      integral ( X <= T < +oo ) T^(A-1) * exp(-T) dT.
!
!    In Mathematica, the function can be evaluated by:
!
!      Gamma[A,X]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+01, &
    0.10D+01, &
    0.10D+01, &
    0.11D+01, &
    0.11D+01, &
    0.11D+01, &
    0.20D+01, &
    0.20D+01, &
    0.20D+01, &
    0.60D+01, &
    0.60D+01, &
    0.11D+02, &
    0.26D+02, &
    0.41D+02 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    2.490302836300570D+00, &
    0.8718369702247978D+00, &
    0.1079213896175866D+00, &
    1.238121685818417D+00, &
    0.3911298052193973D+00, &
    0.01444722098952533D+00, &
    0.9048374180359596D+00, &
    0.3678794411714423D+00, &
    0.006737946999085467D+00, &
    0.8827966752611692D+00, &
    0.3908330082003269D+00, &
    0.008051456628620993D+00, &
    0.9898141728888165D+00, &
    0.5578254003710746D+00, &
    0.007295055724436130D+00, &
    114.9574754165633D+00, &
    2.440923530031405D+00, &
    280854.6620274718D+00, &
    8.576480283455533D+24, &
    2.085031346403364D+47 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.30D-01, &
    0.30D+00, &
    0.15D+01, &
    0.75D-01, &
    0.75D+00, &
    0.35D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.15D+00, &
    0.15D+01, &
    0.70D+01, &
    0.25D+01, &
    0.12D+02, &
    0.16D+02, &
    0.25D+02, &
    0.45D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_LOG_VALUES returns some values of the Log Gamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Log[Gamma[x]]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1524063822430784D+01, &
     0.7966778177017837D+00, &
     0.3982338580692348D+00, &
     0.1520596783998375D+00, &
     0.0000000000000000D+00, &
    -0.4987244125983972D-01, &
    -0.8537409000331584D-01, &
    -0.1081748095078604D+00, &
    -0.1196129141723712D+00, &
    -0.1207822376352452D+00, &
    -0.1125917656967557D+00, &
    -0.9580769740706586D-01, &
    -0.7108387291437216D-01, &
    -0.3898427592308333D-01, &
    0.00000000000000000D+00, &
    0.69314718055994530D+00, &
    0.17917594692280550D+01, &
    0.12801827480081469D+02, &
    0.39339884187199494D+02, &
    0.71257038967168009D+02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.20D+00, &
     0.40D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gegenbauer_poly_values ( n_data, n, a, x, fx )

!*****************************************************************************80
!
!! GEGENBAUER_POLY_VALUES returns some values of the Gegenbauer polynomials.
!
!  Discussion:
!
!    The Gegenbauer polynomials are also known as the "spherical
!    polynomials" or "ultraspherical polynomials".
!
!    In Mathematica, the function can be evaluated by:
!
!      GegenbauerC[n,m,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order parameter of the function.
!
!    Output, real ( kind = 8 ) A, the real parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 38

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00, &
     3.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
      1.0000000000D+00, &
      0.2000000000D+00, &
     -0.4400000000D+00, &
     -0.2800000000D+00, &
      0.2320000000D+00, &
      0.3075200000D+00, &
     -0.0805760000D+00, &
     -0.2935168000D+00, &
     -0.0395648000D+00, &
      0.2459712000D+00, &
      0.1290720256D+00, &
      0.0000000000D+00, &
     -0.3600000000D+00, &
     -0.0800000000D+00, &
      0.8400000000D+00, &
      2.4000000000D+00, &
      4.6000000000D+00, &
      7.4400000000D+00, &
     10.9200000000D+00, &
     15.0400000000D+00, &
     19.8000000000D+00, &
     25.2000000000D+00, &
     -9.0000000000D+00, &
     -0.1612800000D+00, &
     -6.6729600000D+00, &
     -8.3750400000D+00, &
     -5.5267200000D+00, &
      0.0000000000D+00, &
      5.5267200000D+00, &
      8.3750400000D+00, &
      6.6729600000D+00, &
      0.1612800000D+00, &
     -9.0000000000D+00, &
    -15.4252800000D+00, &
     -9.6969600000D+00, &
     22.4409600000D+00, &
    100.8892800000D+00, &
    252.0000000000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10,  2, &
     2,  2,  2, &
     2,  2,  2, &
     2,  2,  2, &
     2,  5,  5, &
     5,  5,  5, &
     5,  5,  5, &
     5,  5,  5, &
     5,  5,  5, &
     5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.20D+00, &
     0.20D+00, &
     0.20D+00, &
     0.20D+00, &
     0.20D+00, &
     0.20D+00, &
     0.20D+00, &
     0.20D+00, &
     0.20D+00, &
     0.20D+00, &
     0.20D+00, &
     0.40D+00, &
     0.40D+00, &
     0.40D+00, &
     0.40D+00, &
     0.40D+00, &
     0.40D+00, &
     0.40D+00, &
     0.40D+00, &
     0.40D+00, &
     0.40D+00, &
     0.40D+00, &
    -0.50D+00, &
    -0.40D+00, &
    -0.30D+00, &
    -0.20D+00, &
    -0.10D+00, &
     0.00D+00, &
     0.10D+00, &
     0.20D+00, &
     0.30D+00, &
     0.40D+00, &
     0.50D+00, &
     0.60D+00, &
     0.70D+00, &
     0.80D+00, &
     0.90D+00, &
     1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine geometric_cdf_values ( n_data, x, p, cdf )

!*****************************************************************************80
!
!! GEOMETRIC_CDF_VALUES returns values of the geometric CDF.
!
!  Discussion:
!
!    The geometric or Pascal probability density function gives the
!    probability that the first success will happen on the X-th Bernoulli
!    trial, given that the probability of a success on a single trial is P.
!
!    The value of CDF ( X, P ) is the probability that the first success
!    will happen on or before the X-th trial.
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = GeometricDistribution [ p ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, Steven Kokoska,
!    Standard Probability and Statistical Tables,
!    CRC Press, 2000,
!    ISBN: 1-58488-059-7,
!    LC: QA273.3.Z95.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) X, the number of trials.
!
!    Output, real ( kind = 8 ) P, the probability of success
!    on one trial.
!
!    Output, real ( kind = 8 ) CDF, the cumulative density function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 14

  real ( kind = 8 ) cdf
  real ( kind = 8 ), save, dimension ( n_max ) :: cdf_vec = (/ &
    0.1900000000000000D+00, &
    0.2710000000000000D+00, &
    0.3439000000000000D+00, &
    0.6861894039100000D+00, &
    0.3600000000000000D+00, &
    0.4880000000000000D+00, &
    0.5904000000000000D+00, &
    0.9141006540800000D+00, &
    0.7599000000000000D+00, &
    0.8704000000000000D+00, &
    0.9375000000000000D+00, &
    0.9843750000000000D+00, &
    0.9995117187500000D+00, &
    0.9999000000000000D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
    0.1D+00, &
    0.1D+00, &
    0.1D+00, &
    0.1D+00, &
    0.2D+00, &
    0.2D+00, &
    0.2D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    0.9D+00 /)
  integer ( kind = 4 ) x
  integer ( kind = 4 ), save, dimension ( n_max ) :: x_vec = (/ &
    1,  2,  3, 10, 1, &
    2,  3, 10,  3, 3, &
    3,  5, 10,  3 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0
    p = 0.0D+00
    cdf = 0.0D+00
  else
    x = x_vec(n_data)
    p = p_vec(n_data)
    cdf = cdf_vec(n_data)
  end if

  return
end
subroutine goodwin_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GOODWIN_VALUES returns some values of the Goodwin and Staton function.
!
!  Discussion:
!
!    The function is defined by:
!
!      GOODWIN(x) = integral ( 0 <= t < +oo ) exp ( -t^2 ) / ( t + x ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.59531540040441651584D+01, &
     0.45769601268624494109D+01, &
     0.32288921331902217638D+01, &
     0.19746110873568719362D+01, &
     0.96356046208697728563D+00, &
     0.60513365250334458174D+00, &
     0.51305506459532198016D+00, &
     0.44598602820946133091D+00, &
     0.37344458206879749357D+00, &
     0.35433592884953063055D+00, &
     0.33712156518881920994D+00, &
     0.29436170729362979176D+00, &
     0.25193499644897222840D+00, &
     0.22028778222123939276D+00, &
     0.19575258237698917033D+00, &
     0.17616303166670699424D+00, &
     0.16015469479664778673D+00, &
     0.14096116876193391066D+00, &
     0.13554987191049066274D+00, &
     0.11751605060085098084D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.0019531250D+00, &
     0.0078125000D+00, &
     0.0312500000D+00, &
     0.1250000000D+00, &
     0.5000000000D+00, &
     1.0000000000D+00, &
     1.2500000000D+00, &
     1.5000000000D+00, &
     1.8750000000D+00, &
     2.0000000000D+00, &
     2.1250000000D+00, &
     2.5000000000D+00, &
     3.0000000000D+00, &
     3.5000000000D+00, &
     4.0000000000D+00, &
     4.5000000000D+00, &
     5.0000000000D+00, &
     5.7500000000D+00, &
     6.0000000000D+00, &
     7.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gud_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GUD_VALUES returns some values of the Gudermannian function.
!
!  Discussion:
!
!    The Gudermannian function relates the hyperbolic and trigonomentric
!    functions.  For any argument X, there is a corresponding value
!    GD so that
!
!      SINH(X) = TAN(GD).
!
!    This value GD is called the Gudermannian of X and symbolized
!    GD(X).  The inverse Gudermannian function is given as input a value
!    GD and computes the corresponding value X.
!
!    GD(X) = 2 * arctan ( exp ( X ) ) - PI / 2
!
!    In Mathematica, the function can be evaluated by:
!
!      2 * Atan[Exp[x]] - Pi/2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1301760336046015D+01, &
    -0.8657694832396586D+00, &
     0.0000000000000000D+00, &
     0.9983374879348662D-01, &
     0.1986798470079397D+00, &
     0.4803810791337294D+00, &
     0.8657694832396586D+00, &
     0.1131728345250509D+01, &
     0.1301760336046015D+01, &
     0.1406993568936154D+01, &
     0.1471304341117193D+01, &
     0.1510419907545700D+01, &
     0.1534169144334733D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -2.0D+00, &
    -1.0D+00, &
     0.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     1.5D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hermite_function_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! HERMITE_FUNCTION_VALUES: values of the Hermite function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Hf(n,x) = HermiteH[n,x] 
!        * Exp [ -1/2 * x^2] / Sqrt [ 2^n * n! * Sqrt[Pi] ]
!
!    The Hermite functions are orthonormal:
!
!      integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 23

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7511255444649425D+00,  0.0000000000000000D+00, -0.5311259660135985D+00, &
    0.0000000000000000D+00,  0.4599685791773266D+00,  0.0000000000000000D+00, &
    0.4555806720113325D+00,  0.6442883651134752D+00,  0.3221441825567376D+00, &
   -0.2630296236233334D+00, -0.4649750762925110D+00, -0.5881521185179581D-01, &
    0.3905052515434106D+00,  0.2631861423064045D+00, -0.2336911435996523D+00, &
   -0.3582973361472840D+00,  0.6146344487883041D-01,  0.3678312067984882D+00, &
    0.9131969309166278D-01,  0.4385750950032321D+00, -0.2624689527931006D-01, &
    0.5138426125477819D+00,  0.9355563118061758D-01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  2,  &
    3,  4,  5,  &
    0,  1,  2,  &
    3,  4,  5,  &
    6,  7,  8,  &
    9, 10, 11,  &
   12,  5,  5,  &
    5,  5  /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 0.5D+00, 2.0D+00, &
    3.0D+00, 4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hermite_poly_phys_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! HERMITE_POLY_PHYS_VALUES: values of the physicist's Hermite polynomial.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      HermiteH[n,x]
!
!  First terms:
!
!      1
!      2 X
!      4 X^2     -  2
!      8 X^3     - 12 X
!     16 X^4     - 48 X^2     + 12
!     32 X^5    - 160 X^3    + 120 X
!     64 X^6    - 480 X^4    + 720 X^2    - 120
!    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
!    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
!    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
!   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
!    = sqrt ( PI ) * 2^N * N!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
      0.1000000000000000D+01, &
      0.1000000000000000D+02, &
      0.9800000000000000D+02, &
      0.9400000000000000D+03, &
      0.8812000000000000D+04, &
      0.8060000000000000D+05, &
      0.7178800000000000D+06, &
      0.6211600000000000D+07, &
      0.5206568000000000D+08, &
      0.4212712000000000D+09, &
      0.3275529760000000D+10, &
      0.2432987360000000D+11, &
      0.1712370812800000D+12, &
      0.0000000000000000D+00, &
      0.4100000000000000D+02, &
     -0.8000000000000000D+01, &
      0.3816000000000000D+04, &
      0.3041200000000000D+07 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5, 5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    3.0D+00, &
    1.0D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hermite_poly_prob_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! HERMITE_POLY_PROB_VALUES: values of the probabilist's Hermite polynomial.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      He(n,x) = HermiteH[n,x/Sqrt[2]] / Sqrt [ 2^n ] 
!
!  First terms:
!
!   1
!   X
!   X^2  -  1
!   X^3  -  3 X
!   X^4  -  6 X^2 +   3
!   X^5  - 10 X^3 +  15 X
!   X^6  - 15 X^4 +  45 X^2 -   15
!   X^7  - 21 X^5 + 105 X^3 -  105 X
!   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
!   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
!   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
!
!  Recursion:
!
!    He(0,X) = 1,
!    He(1,X) = X,
!    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
!
!  Norm:
!
!    integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
!    = sqrt ( 2 * pi ) * N! * delta ( M, N )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.000000000000000D+00, &
    5.000000000000000D+00, &
    24.00000000000000D+00, &
    110.0000000000000D+00, &
    478.0000000000000D+00, &
    1950.000000000000D+00, &
    7360.000000000000D+00, &
    25100.00000000000D+00, &
    73980.00000000000D+00, &
    169100.0000000000D+00, &
    179680.0000000000D+00, &
   -792600.0000000000D+00, &
   -5939480.000000000D+00, &
    0.000000000000000D+00, &
    6.281250000000000D+00, &
    6.000000000000000D+00, &
    18.00000000000000D+00, &
    90150.00000000000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    3.0D+00, &
    1.0D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hyper_1f1_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! HYPER_1F1_VALUES returns some values of the hypergeometric function 1F1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = Hypergeometric1F1 [ a, b, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, X, the parameters.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 24

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) ::  a_vec = (/ &
    -2.500D+00, &
    -0.500D+00, &
     0.500D+00, &
     2.500D+00, &
    -2.500D+00, &
    -0.500D+00, &
     0.500D+00, &
     2.500D+00, &
    -2.500D+00, &
    -0.500D+00, &
     0.500D+00, &
     2.500D+00, &
     0.825D+00, &
     1.100D+00, &
     1.650D+00, &
     3.300D+00, &
     0.825D+00, &
     1.100D+00, &
     1.650D+00, &
     3.300D+00, &
     0.825D+00, &
     1.100D+00, &
     1.650D+00,&
     3.300D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) ::  b_vec = (/ &
     3.3D+00, &
     1.1D+00, &
     1.1D+00, &
     3.3D+00, &
     3.3D+00, &
     1.1D+00, &
     1.1D+00, &
     3.3D+00, &
     3.3D+00, &
     1.1D+00, &
     1.1D+00, &
     3.3D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) ::  fx_vec = (/ &
     0.81879926689265186854D+00, &
     0.88283984828032972070D+00, &
     1.1245023764952626690D+00, &
     1.2101049301639599598D+00, &
     0.12723045536781567174D+00, &
     0.12326016871544045107D+00, &
     2.3297954665128293051D+00, &
     3.3890020264468009733D+00, &
    -0.18819510282516768874D+00, &
    -1.0764203806547022727D+00, &
     5.7521824680907968433D+00, &
     9.9998567403304086593D+00, &
     1.0317208964319891384D+00, &
     1.0424867029249952040D+00, &
     1.0643112000949092012D+00, &
     1.1321844369742336326D+00, &
     1.2328402688568452181D+00, &
     1.3200654482027340732D+00, &
     1.5104811522310825217D+00, &
     2.2307520785940524365D+00, &
     1.5197286298183137741D+00, &
     1.7364938170250847619D+00, &
     2.2492330307668135926D+00, &
     4.6377737119178965298D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) ::  x_vec = (/ &
     0.25D+00, &
     0.25D+00, &
     0.25D+00, &
     0.25D+00, &
     1.55D+00, &
     1.55D+00, &
     1.55D+00, &
     1.55D+00, &
     2.85D+00, &
     2.85D+00, &
     2.85D+00, &
     2.85D+00, &
     0.25D+00, &
     0.25D+00, &
     0.25D+00, &
     0.25D+00, &
     1.55D+00, &
     1.55D+00, &
     1.55D+00, &
     1.55D+00, &
     2.85D+00, &
     2.85D+00, &
     2.85D+00, &
     2.85D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hyper_2f1_values ( n_data, a, b, c, x, fx )

!*****************************************************************************80
!
!! HYPER_2F1_VALUES returns some values of the hypergeometric function 2F1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = Hypergeometric2F1 [ a, b, c, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, C, X, the parameters.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 24

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00 /)
  real ( kind = 8 ) c
  real ( kind = 8 ), save, dimension ( n_max ) :: c_vec = (/ &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.72356129348997784913D+00, &
    0.97911109345277961340D+00, &
    1.0216578140088564160D+00, &
    1.4051563200112126405D+00, &
    0.46961431639821611095D+00, &
    0.95296194977446325454D+00, &
    1.0512814213947987916D+00, &
    2.3999062904777858999D+00, &
    0.29106095928414718320D+00, &
    0.92536967910373175753D+00, &
    1.0865504094806997287D+00, &
    5.7381565526189046578D+00, &
    15090.669748704606754D+00, &
   -104.31170067364349677D+00, &
    21.175050707768812938D+00, &
    4.1946915819031922850D+00, &
    1.0170777974048815592D+10, &
   -24708.635322489155868D+00, &
    1372.2304548384989560D+00, &
    58.092728706394652211D+00, &
    5.8682087615124176162D+18, &
   -4.4635010147295996680D+08, &
    5.3835057561295731310D+06, &
    20396.913776019659426D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    c = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    c = c_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hypergeometric_cdf_values ( n_data, sam, suc, pop, n, fx )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
!
!  Discussion:
!
!    CDF(X)(A,B) is the probability of at most X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!    In Mathematica, the function can be evaluated by:
!
!      dist = HypergeometricDistribution [ sam, suc, pop ]
!      CDF [ dist, n ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) SAM, SUC, POP, the sample size,
!    success size, and population parameters of the function.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6001858177500578D-01, &
    0.2615284665839845D+00, &
    0.6695237889132748D+00, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5332595856827856D+00, &
    0.1819495964117640D+00, &
    0.4448047017527730D-01, &
    0.9999991751316731D+00, &
    0.9926860896560750D+00, &
    0.8410799901444538D+00, &
    0.3459800113391901D+00, &
    0.0000000000000000D+00, &
    0.2088888139634505D-02, &
    0.3876752992448843D+00, &
    0.9135215248834896D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     7,  8,  9, 10, &
     6,  6,  6,  6, &
     6,  6,  6,  6, &
     0,  0,  0,  0 /)
  integer ( kind = 4 ) pop
  integer ( kind = 4 ), save, dimension ( n_max ) :: pop_vec = (/ &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    90,  200, 1000, 10000 /)
  integer ( kind = 4 ) sam
  integer ( kind = 4 ), save, dimension ( n_max ) :: sam_vec = (/ &
    10, 10, 10, 10, &
     6,  7,  8,  9, &
    10, 10, 10, 10, &
    10, 10, 10, 10 /)
  integer ( kind = 4 ) suc
  integer ( kind = 4 ), save, dimension ( n_max ) :: suc_vec = (/ &
    90, 90, 90, 90, &
    90, 90, 90, 90, &
    10, 30, 50, 70, &
    90, 90, 90, 90 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    sam = 0
    suc = 0
    pop = 0
    n = 0
    fx = 0.0D+00
  else
    sam = sam_vec(n_data)
    suc = suc_vec(n_data)
    pop = pop_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hypergeometric_pdf_values ( n_data, sam, suc, pop, n, fx )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_PDF_VALUES returns some values of the hypergeometric PDF.
!
!  Discussion:
!
!    PDF(X)(A,B) is the probability of X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!    In Mathematica, the function can be evaluated by:
!
!      dist = HypergeometricDistribution [ sam, suc, pop ]
!      PDF [ dist, n ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) SAM, SUC, POP, the sample size,
!    success size, and population parameters of the function.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.05179370533242827D+00, &
    0.2015098848089788D+00, &
    0.4079953223292903D+00, &
    0.3304762110867252D+00, &
    0.5223047493549780D+00, &
    0.3889503452643453D+00, &
    0.1505614239732950D+00, &
    0.03927689321042477D+00, &
    0.00003099828465518108D+00, &
    0.03145116093938197D+00, &
    0.2114132170316862D+00, &
    0.2075776621999210D+00, &
    0.0000000000000000D+00, &
    0.002088888139634505D+00, &
    0.3876752992448843D+00, &
    0.9135215248834896D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     7,  8,  9, 10, &
     6,  6,  6,  6, &
     6,  6,  6,  6, &
     0,  0,  0,  0 /)
  integer ( kind = 4 ) pop
  integer ( kind = 4 ), save, dimension ( n_max ) :: pop_vec = (/ &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    90,  200, 1000, 10000 /)
  integer ( kind = 4 ) sam
  integer ( kind = 4 ), save, dimension ( n_max ) :: sam_vec = (/ &
    10, 10, 10, 10, &
     6,  7,  8,  9, &
    10, 10, 10, 10, &
    10, 10, 10, 10 /)
  integer ( kind = 4 ) suc
  integer ( kind = 4 ), save, dimension ( n_max ) :: suc_vec = (/ &
    90, 90, 90, 90, &
    90, 90, 90, 90, &
    10, 30, 50, 70, &
    90, 90, 90, 90 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    sam = 0
    suc = 0
    pop = 0
    n = 0
    fx = 0.0D+00
  else
    sam = sam_vec(n_data)
    suc = suc_vec(n_data)
    pop = pop_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hypergeometric_u_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_U_VALUES: some values of the hypergeometric function U(a,b,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = HypergeometricU [ a, b, x ]
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
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, X, the parameters.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 24

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) ::  a_vec = (/ &
    -2.500D+00, &
    -0.500D+00, &
     0.500D+00, &
     2.500D+00, &
    -2.500D+00, &
    -0.500D+00, &
     0.500D+00, &
     2.500D+00, &
    -2.500D+00, &
    -0.500D+00, &
     0.500D+00, &
     2.500D+00, &
     0.825D+00, &
     1.100D+00, &
     1.650D+00, &
     3.300D+00, &
     0.825D+00, &
     1.100D+00, &
     1.650D+00, &
     3.300D+00, &
     0.825D+00, &
     1.100D+00, &
     1.650D+00,&
     3.300D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) ::  b_vec = (/ &
     3.3D+00, &
     1.1D+00, &
     1.1D+00, &
     3.3D+00, &
     3.3D+00, &
     1.1D+00, &
     1.1D+00, &
     3.3D+00, &
     3.3D+00, &
     1.1D+00, &
     1.1D+00, &
     3.3D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00, &
     6.7D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) ::  fx_vec = (/ &
         -68.693628728078601389D+00, &
          -0.0029710551374761070801D+00, &
           1.5008631742177797301D+00, &
          20.614688244200596134D+00, &
           7.4563815469305551938D+00, &
           1.0155793767749293733D+00, &
           0.73446538936622668912D+00, &
           0.28046404941879399225D+00, &
           3.4508153741446547607D+00, &
           1.5156637368753063495D+00, &
           0.56042118587934993510D+00, &
           0.064897147735134223341D+00, &
      223432.02356977463356D+00, &
      263079.25980740811495D+00, &
      269802.90319351274132D+00, &
       82809.311335606553425D+00, &
          26.465684783131844524D+00, &
          28.093506172516056560D+00, &
          23.889164624518872504D+00, &
           4.5338847857070388229D+00, &
           3.0224469362694842535D+00, &
           2.8040650913713359934D+00, &
           1.9262578111480172682D+00, &
           0.23020518115860909098D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) ::  x_vec = (/ &
     0.25D+00, &
     0.25D+00, &
     0.25D+00, &
     0.25D+00, &
     1.55D+00, &
     1.55D+00, &
     1.55D+00, &
     1.55D+00, &
     2.85D+00, &
     2.85D+00, &
     2.85D+00, &
     2.85D+00, &
     0.25D+00, &
     0.25D+00, &
     0.25D+00, &
     0.25D+00, &
     1.55D+00, &
     1.55D+00, &
     1.55D+00, &
     1.55D+00, &
     2.85D+00, &
     2.85D+00, &
     2.85D+00, &
     2.85D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine i0ml0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! I0ML0_VALUES returns some values of the I0ML0 function.
!
!  Discussion:
!
!    The function is defined by:
!
!      I0ML0(x) = I0(x) - L0(x)
!
!    I0(x) is the modified Bessel function of the first kind of order 0,
!    L0(x) is the modified Struve function of order 0.
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.99875755515461749793D+00, &
    0.99011358230706643807D+00, &
    0.92419435310023947018D+00, &
    0.73624267134714273902D+00, &
    0.55582269181411744686D+00, &
    0.34215154434462160628D+00, &
    0.17087174888774706539D+00, &
    0.81081008709219208918D-01, &
    0.53449421441089580702D-01, &
    0.39950321008923244846D-01, &
    0.39330637437584921392D-01, &
    0.37582274342808670750D-01, &
    0.31912486554480390343D-01, &
    0.25506146883504738403D-01, &
    0.21244480317825292412D-01, &
    0.15925498348551684335D-01, &
    0.12737506927242585015D-01, &
    0.84897750814784916847D-02, &
    0.63668349178454469153D-02, &
    0.50932843163122551114D-02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0156250000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      2.0000000000D+00, &
      4.0000000000D+00, &
      8.0000000000D+00, &
     12.0000000000D+00, &
     16.0000000000D+00, &
     16.2500000000D+00, &
     17.0000000000D+00, &
     20.0000000000D+00, &
     25.0000000000D+00, &
     30.0000000000D+00, &
     40.0000000000D+00, &
     50.0000000000D+00, &
     75.0000000000D+00, &
    100.0000000000D+00, &
    125.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine i1ml1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! I1ML1_VALUES returns some values of the I1ML1 function.
!
!  Discussion:
!
!    The function is defined by:
!
!      I1ML1(x) = I1(x) - L1(x)
!
!    I1(x) is the modified Bessel function of the first kind of order 1,
!    L1(x) is the modified Struve function of order 1.
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.97575346155386267134D-03, &
    0.77609293280609272733D-02, &
    0.59302966404545373770D-01, &
    0.20395212276737365307D+00, &
    0.33839472293667639038D+00, &
    0.48787706726961324579D+00, &
    0.59018734196576517506D+00, &
    0.62604539530312149476D+00, &
    0.63209315274909764698D+00, &
    0.63410179313235359215D+00, &
    0.63417966797578128188D+00, &
    0.63439268632392089434D+00, &
    0.63501579073257770690D+00, &
    0.63559616677359459337D+00, &
    0.63591001826697110312D+00, &
    0.63622113181751073643D+00, &
    0.63636481702133606597D+00, &
    0.63650653499619902120D+00, &
    0.63655609126300261851D+00, &
    0.63657902087183929223D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0156250000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      2.0000000000D+00, &
      4.0000000000D+00, &
      8.0000000000D+00, &
     12.0000000000D+00, &
     16.0000000000D+00, &
     16.2500000000D+00, &
     17.0000000000D+00, &
     20.0000000000D+00, &
     25.0000000000D+00, &
     30.0000000000D+00, &
     40.0000000000D+00, &
     50.0000000000D+00, &
     75.0000000000D+00, &
    100.0000000000D+00, &
    125.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine int_values ( n_data, x, fx )

!*****************************************************************************80
!
!! INT_VALUES returns some values of the "integer part" function.
!
!  Discussion:
!
!    INT(X) = returns the integer part of a real number.
!
!    The result is returned as a real number.
!
!    The result is computed by rounding the absolute value of the
!    input towards 0, and then restoring the sign.
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
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     -2.00D+00, &
     -1.00D+00, &
     -1.00D+00, &
     -1.00D+00, &    
     -1.00D+00, &    
     -1.00D+00, &      
      0.00D+00, &
      0.00D+00, &    
      0.00D+00, &
      0.00D+00, &
      0.00D+00, &
      0.00D+00, &
      0.00D+00, &
      0.00D+00, &
      0.00D+00, &
      0.00D+00, &
      0.00D+00, &
      0.00D+00, &
      0.00D+00, &
      1.00D+00, &
      1.00D+00, &
      1.00D+00, &
      1.00D+00, &
      1.00D+00, &
      2.00D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     -2.01D+00, &
     -1.99D+00, &
     -1.50D+00, &
     -1.10D+00, &    
     -1.01D+00, &    
     -1.00D+00, &      
     -0.99D+00, &
     -0.90D+00, &    
     -0.51D+00, &
     -0.50D+00, &
     -0.49D+00, &
     -0.01D+00, &
      0.00D+00, &
      0.01D+00, &
      0.49D+00, &
      0.50D+00, &
      0.51D+00, &
      0.90D+00, &
      0.99D+00, &
      1.00D+00, &
      1.01D+00, &
      1.10D+00, &
      1.50D+00, &
      1.99D+00, &
      2.01D+00  /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine jacobi_cn_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! JACOBI_CN_VALUES returns some values of the Jacobi elliptic function CN(A,X).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      JacobiCN[ x, a ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.9950041652780258D+00, &
     0.9800665778412416D+00, &
     0.8775825618903727D+00, &
     0.5403023058681397D+00, &
    -0.4161468365471424D+00, &
     0.9950124626090582D+00, &
     0.9801976276784098D+00, &
     0.8822663948904403D+00, &
     0.5959765676721407D+00, &
    -0.1031836155277618D+00, &
     0.9950207489532265D+00, &
     0.9803279976447253D+00, &
     0.8868188839700739D+00, &
     0.6480542736638854D+00, &
     0.2658022288340797D+00, &
     0.3661899347368653D-01, &
     0.9803279976447253D+00, &
     0.8868188839700739D+00, &
     0.6480542736638854D+00, &
     0.2658022288340797D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     4.0D+00, &
    -0.2D+00, &
    -0.5D+00, &
    -1.0D+00, &
    -2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine jacobi_dn_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! JACOBI_DN_VALUES returns some values of the Jacobi elliptic function DN(A,X).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      JacobiDN[ x, a ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.9975093485144243D+00, &
    0.9901483195224800D+00, &
    0.9429724257773857D+00, &
    0.8231610016315963D+00, &
    0.7108610477840873D+00, &
    0.9950207489532265D+00, &
    0.9803279976447253D+00, &
    0.8868188839700739D+00, &
    0.6480542736638854D+00, &
    0.2658022288340797D+00, &
    0.3661899347368653D-01, &
    0.9803279976447253D+00, &
    0.8868188839700739D+00, &
    0.6480542736638854D+00, &
    0.2658022288340797D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     4.0D+00, &
    -0.2D+00, &
    -0.5D+00, &
    -1.0D+00, &
    -2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine jacobi_poly_values ( n_data, n, a, b, x, fx )

!*****************************************************************************80
!
!! JACOBI_POLY_VALUES returns some values of the Jacobi polynomial.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      JacobiP[ n, a, b, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Output, real ( kind = 8 ) A, B, parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 26

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, &
     3.0D+00, 4.0D+00, 5.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, &
    3.0D+00, 4.0D+00, 5.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.2500000000000000D+00, &
    -0.3750000000000000D+00, &
    -0.4843750000000000D+00, &
    -0.1328125000000000D+00, &
     0.2753906250000000D+00, &
    -0.1640625000000000D+00, &
    -0.1174804687500000D+01, &
    -0.2361328125000000D+01, &
    -0.2616210937500000D+01, &
     0.1171875000000000D+00, &
     0.4218750000000000D+00, &
     0.5048828125000000D+00, &
     0.5097656250000000D+00, &
     0.4306640625000000D+00, &
    -0.6000000000000000D+01, &
     0.3862000000000000D-01, &
     0.8118400000000000D+00, &
     0.3666000000000000D-01, &
    -0.4851200000000000D+00, &
    -0.3125000000000000D+00, &
     0.1891200000000000D+00, &
     0.4023400000000000D+00, &
     0.1216000000000000D-01, &
    -0.4396200000000000D+00, &
     0.1000000000000000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, 1, 2, 3, &
     4, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
    -1.0D+00, &
    -0.8D+00, &
    -0.6D+00, &
    -0.4D+00, &
    -0.2D+00, &
     0.0D+00, &
     0.2D+00, &
     0.4D+00, &
     0.6D+00, &
     0.8D+00, &
     1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine jacobi_sn_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! JACOBI_SN_VALUES returns some values of the Jacobi elliptic function SN(A,X).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      JacobiSN[ x, a ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.9983341664682815D-01, &
     0.1986693307950612D+00, &
     0.4794255386042030D+00, &
     0.8414709848078965D+00, &
     0.9092974268256817D+00, &
     0.9975068547462484D-01, &
     0.1980217429819704D+00, &
     0.4707504736556573D+00, &
     0.8030018248956439D+00, &
     0.9946623253580177D+00, &
     0.9966799462495582D-01, &
     0.1973753202249040D+00, &
     0.4621171572600098D+00, &
     0.7615941559557649D+00, &
     0.9640275800758169D+00, &
     0.9993292997390670D+00, &
    -0.1973753202249040D+00, &
    -0.4621171572600098D+00, &
    -0.7615941559557649D+00, &
    -0.9640275800758169D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     0.1D+00, &
     0.2D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     4.0D+00, &
    -0.2D+00, &
    -0.5D+00, &
    -1.0D+00, &
    -2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine jed_ce_values ( n_data, jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_CE_VALUES returns the Common Era dates for Julian Ephemeris Dates.
!
!  Discussion:
!
!    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
!    starting from noon on 1 January 4713 BCE.
!
!    The CE or Common Era is the day, month and year under the
!    hybrid Julian/Gregorian Calendar, with a transition from Julian
!    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
!
!    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
!    years BC/BCE are indicated by a negative year value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Reingold, Nachum Dershowitz,
!    Calendrical Calculations: The Millennium Edition,
!    Cambridge University Press, 2001,
!    ISBN: 0 521 77752 6
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, the Common Era date.
!
!    Output, real ( kind = 8 ) F, the fractional part of the day.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 51

  integer ( kind = 4 ) d
  integer ( kind = 4 ), save, dimension ( n_max ) :: d_vec = (/ &
    01, &
    02, &
    26, &
    08, &
    06, &
    18, &
    08, &
    09, &
    01, &
    26, &
    26, &
    01, &
    01, &
    29, &
    31, &
    01, &
    03, &
    03, &
    29, &
    24, &
    24, &
    29, &
    03, &
    11, &
    12, &
    24, &
    19, &
    15, &
    16, &
    16, &
    21, &
    17, &
    09, &
    04, &
    15, &
    04, &
    13, &
    14, &
    18, &
    22, &
    21, &
    24, &
    17, &
    31, &
    01, &
    06, &
    25, &
    01, &
    09, &
    23, &
    01 /)
  real ( kind = 8 ) f
  real ( kind = 8 ), save, dimension ( n_max ) :: f_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.00D+00, &
    0.00D+00, &
    0.25D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.50D+00, &
    0.50D+00, &
    0.00D+00, &
    0.50D+00, &
    0.50D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.81D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.33D+00, &
    0.00D+00, &
    0.50D+00 /)
  real ( kind = 8 ) jed
  real ( kind = 8 ), save, dimension ( n_max ) :: jed_vec = (/ &
           0.00D+00, &
           1.00D+00, &
      259261.00D+00, &
      347998.50D+00, &
      584282.50D+00, &
      588465.75D+00, &
      758325.50D+00, &
     1438178.50D+00, &
     1446389.50D+00, &
     1448637.50D+00, &
     1448637.50D+00, &
     1607708.50D+00, &
     1607738.50D+00, &
     1713262.50D+00, &
     1721422.50D+00, &
     1721423.50D+00, &
     1721425.50D+00, &
     1721425.50D+00, &
     1724220.50D+00, &
     1741959.50D+00, &
     1749994.50D+00, &
     1825029.50D+00, &
     1862836.50D+00, &
     1922867.50D+00, &
     1936747.50D+00, &
     1940351.50D+00, &
     1948320.50D+00, &
     1948438.50D+00, &
     1948439.50D+00, &
     1952062.50D+00, &
     1952067.50D+00, &
     2114872.50D+00, &
     2289425.50D+00, &
     2299160.00D+00, &
     2299161.00D+00, &
     2333269.50D+00, &
     2361221.00D+00, &
     2361222.00D+00, &
     2372547.50D+00, &
     2375839.50D+00, &
     2394646.50D+00, &
     2394710.50D+00, &
     2400000.50D+00, &
     2415020.31D+00, &
     2440587.50D+00, &
     2444244.50D+00, &
     2450138.50D+00, &
     2451544.50D+00, &
     2453073.83D+00, &
     2456284.50D+00, &
     2913943.00D+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
     01, &
     01, &
     10, &
     10, &
     09, &
     02, &
     03, &
     07, &
     01, &
     02, &
     02, &
     09, &
     10, &
     08, &
     12, &
     01, &
     01, &
     01, &
     08, &
     03, &
     03, &
     08, &
     03, &
     07, &
     07, &
     05, &
     03, &
     07, &
     07, &
     06, &
     06, &
     03, &
     02, &
     10, &
     10, &
     03, &
     09, &
     09, &
     09, &
     09, &
     03, &
     05, &
     11, &
     12, &
     01, &
     01, &
     02, &
     01, &
     03, &
     12, &
     01 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) y
  integer ( kind = 4 ), save, dimension ( n_max ) :: y_vec = (/ &
    -4713, &
    -4713, &
    -4004, &
    -3761, &
    -3114, &
    -3102, &
    -2637, &
     -776, &
     -753, &
     -747, &
     -747, &
     -312, &
     -312, &
      -23, &
       -1, &
        1, &
        1, &
        1, &
        8, &
       57, &
       79, &
      284, &
      388, &
      552, &
      590, &
      600, &
      622, &
      622, &
      622, &
      632, &
      632, &
     1078, &
     1556, &
     1582, &
     1582, &
     1676, &
     1752, &
     1752, &
     1783, &
     1792, &
     1844, &
     1844, &
     1858, &
     1899, &
     1970, &
     1980, &
     1996, &
     2000, &
     2004, &
     2012, &
     3266 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    jed = 0.0D+00
    y = 0
    m = 0
    d = 0
    f = 0.0D+00
  else
    jed = jed_vec(n_data)
    y = y_vec(n_data)
    m = m_vec(n_data)
    d = d_vec(n_data)
    f = f_vec(n_data)
  end if

  return
end
subroutine jed_mjd_values ( n_data, jed, mjd )

!*****************************************************************************80
!
!! JED_MJD_VALUES returns the MJD for Julian Ephemeris Dates.
!
!  Discussion:
!
!    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
!    starting from noon on 1 January 4713 BCE.
!
!    The MJD (Modified Julian Day) counts days starting from midnight,
!    17 November 1858.  This essentially subtracts 2400000.5 days from the JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Reingold, Nachum Dershowitz,
!    Calendrical Calculations: The Millennium Edition,
!    Cambridge University Press, 2001,
!    ISBN: 0 521 77752 6
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) MJD, the Modified Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 33

  real ( kind = 8 ) jed
  real ( kind = 8 ), save, dimension ( n_max ) :: jed_vec = (/ &
    1507231.5D+00, &
    1660037.5D+00, &
    1746893.5D+00, &
    1770641.5D+00, &
    1892731.5D+00, &
    1931579.5D+00, &
    1974851.5D+00, &
    2091164.5D+00, &
    2121509.5D+00, &
    2155779.5D+00, &
    2174029.5D+00, &
    2191584.5D+00, &
    2195261.5D+00, &
    2229274.5D+00, &
    2245580.5D+00, &
    2266100.5D+00, &
    2288542.5D+00, &
    2290901.5D+00, &
    2323140.5D+00, &
    2334848.5D+00, &
    2348020.5D+00, &
    2366978.5D+00, &
    2385648.5D+00, &
    2392825.5D+00, &
    2416223.5D+00, &
    2425848.5D+00, &
    2430266.5D+00, &
    2430833.5D+00, &
    2431004.5D+00, &
    2448698.5D+00, &
    2450138.5D+00, &
    2465737.5D+00, &
    2486076.5D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) mjd
  real ( kind = 8 ), save, dimension ( n_max ) :: mjd_vec = (/ &
    -892769.0D+00, &
    -739963.0D+00, &
    -653107.0D+00, &
    -629359.0D+00, &
    -507269.0D+00, &
    -468421.0D+00, &
    -425149.0D+00, &
    -308836.0D+00, &
    -278491.0D+00, &
    -244221.0D+00, &
    -225971.0D+00, &
    -208416.0D+00, &
    -204739.0D+00, &
    -170726.0D+00, &
    -154420.0D+00, &
    -133900.0D+00, &
    -111458.0D+00, &
    -109099.0D+00, &
     -76860.0D+00, &
     -65152.0D+00, &
     -51980.0D+00, &
     -33022.0D+00, &
     -14352.0D+00, &
      -7175.0D+00, &
      16223.0D+00, &
      25848.0D+00, &
      30266.0D+00, &
      30833.0D+00, &
      31004.0D+00, &
      48698.0D+00, &
      50138.0D+00, &
      65737.0D+00, &
      86076.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    jed = 0.0D+00
    mjd = 0.0D+00
  else
    jed = jed_vec(n_data)
    mjd = mjd_vec(n_data)
  end if

  return
end
subroutine jed_rd_values ( n_data, jed, rd )

!*****************************************************************************80
!
!! JED_RD_VALUES returns the RD for Julian Ephemeris Dates.
!
!  Discussion:
!
!    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
!    starting from noon on 1 January 4713 BCE.
!
!    The RD is the Reingold Dershowitz Date, which counts days from
!    midnight, 1 January year 1 in the Gregorian calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Reingold, Nachum Dershowitz,
!    Calendrical Calculations: The Millennium Edition,
!    Cambridge University Press, 2001,
!    ISBN: 0 521 77752 6
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) RD, the Reingold Dershowitz Date.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 33

  real ( kind = 8 ) jed
  real ( kind = 8 ), save, dimension ( n_max ) :: jed_vec = (/ &
    1507231.5D+00, &
    1660037.5D+00, &
    1746893.5D+00, &
    1770641.5D+00, &
    1892731.5D+00, &
    1931579.5D+00, &
    1974851.5D+00, &
    2091164.5D+00, &
    2121509.5D+00, &
    2155779.5D+00, &
    2174029.5D+00, &
    2191584.5D+00, &
    2195261.5D+00, &
    2229274.5D+00, &
    2245580.5D+00, &
    2266100.5D+00, &
    2288542.5D+00, &
    2290901.5D+00, &
    2323140.5D+00, &
    2334848.5D+00, &
    2348020.5D+00, &
    2366978.5D+00, &
    2385648.5D+00, &
    2392825.5D+00, &
    2416223.5D+00, &
    2425848.5D+00, &
    2430266.5D+00, &
    2430833.5D+00, &
    2431004.5D+00, &
    2448698.5D+00, &
    2450138.5D+00, &
    2465737.5D+00, &
    2486076.5D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) rd
  real ( kind = 8 ), save, dimension ( n_max ) :: rd_vec = (/ &
    -214193.0D+00, &
     -61387.0D+00, &
      25469.0D+00, &
      49217.0D+00, &
     171307.0D+00, &
     210155.0D+00, &
     253427.0D+00, &
     369740.0D+00, &
     400085.0D+00, &
     434355.0D+00, &
     452605.0D+00, &
     470160.0D+00, &
     473837.0D+00, &
     507850.0D+00, &
     524156.0D+00, &
     544676.0D+00, &
     567118.0D+00, &
     569477.0D+00, &
     601716.0D+00, &
     613424.0D+00, &
     626596.0D+00, &
     645554.0D+00, &
     664224.0D+00, &
     671401.0D+00, &
     694799.0D+00, &
     704424.0D+00, &
     708842.0D+00, &
     709409.0D+00, &
     709580.0D+00, &
     727274.0D+00, &
     728714.0D+00, &
     744313.0D+00, &
     764652.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    jed = 0.0D+00
    rd = 0.0D+00
  else
    jed = jed_vec(n_data)
    rd = rd_vec(n_data)
  end if

  return
end
subroutine jed_weekday_values ( n_data, jed, weekday )

!*****************************************************************************80
!
!! JED_WEEKDAY_VALUES returns the day of the week for Julian Ephemeris Dates.
!
!  Discussion:
!
!    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
!    starting from noon on 1 January 4713 BCE.
!
!    Weekdays are numbered as follows:
!
!    1  Sunday
!    2  Monday
!    3  Tuesday
!    4  Wednesday
!    5  Thursday
!    6  Friday
!    7  Saturday
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Reingold, Nachum Dershowitz,
!    Calendrical Calculations: The Millennium Edition,
!    Cambridge University Press, 2001,
!    ISBN: 0 521 77752 6
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) WEEKDAY, the day of the week.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 33

  real ( kind = 8 ) jed
  real ( kind = 8 ), save, dimension ( n_max ) :: jed_vec = (/ &
    1507231.5D+00, &
    1660037.5D+00, &
    1746893.5D+00, &
    1770641.5D+00, &
    1892731.5D+00, &
    1931579.5D+00, &
    1974851.5D+00, &
    2091164.5D+00, &
    2121509.5D+00, &
    2155779.5D+00, &
    2174029.5D+00, &
    2191584.5D+00, &
    2195261.5D+00, &
    2229274.5D+00, &
    2245580.5D+00, &
    2266100.5D+00, &
    2288542.5D+00, &
    2290901.5D+00, &
    2323140.5D+00, &
    2334848.5D+00, &
    2348020.5D+00, &
    2366978.5D+00, &
    2385648.5D+00, &
    2392825.5D+00, &
    2416223.5D+00, &
    2425848.5D+00, &
    2430266.5D+00, &
    2430833.5D+00, &
    2431004.5D+00, &
    2448698.5D+00, &
    2450138.5D+00, &
    2465737.5D+00, &
    2486076.5D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) weekday
  integer ( kind = 4 ), save, dimension ( n_max ) :: weekday_vec = (/ &
    1, 4, 4, 1, 4, &
    2, 7, 1, 1, 6, &
    7, 6, 1, 1, 4, &
    7, 7, 7, 4, 1, &
    6, 1, 2, 4, 1, &
    1, 2, 2, 5, 3, &
    1, 4, 1 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    jed = 0.0D+00
    weekday = 0
  else
    jed = jed_vec(n_data)
    weekday = weekday_vec(n_data)
  end if

  return
end
subroutine kei0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! KEI0_VALUES returns some values of the Kelvin KEI function of order NU = 0.
!
!  Discussion:
!
!    The function is defined by:
!
!      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
!
!    where K(NU,X) is the K Bessel function.
!
!    In Mathematica, KEI(NU,X) can be defined by:
!
!      Im [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    LC: QA76.95.W65,
!    ISBN: 0-521-64314-7.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    -0.6715816950943676D+00, &
    -0.4949946365187199D+00, &
    -0.3313955623385585D+00, &
    -0.2024000677647043D+00, &
    -0.1106960991556749D+00, &
   -0.05112188404598678D+00, &
   -0.01600256851827124D+00, &
   0.002198399294972520D+00, &
   0.009720918540151990D+00, &
   0.01118758650986964D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    2.5D+00, &
    3.0D+00, &
    3.5D+00, &
    4.0D+00, &
    4.5D+00, &
    5.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine kei1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! KEI1_VALUES returns some values of the Kelvin KEI function of order NU = 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
!
!    where K(NU,X) is the K Bessel function.
!
!    In Mathematica, KEI(NU,X) can be defined by:
!
!      Im [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    LC: QA76.95.W65,
!    ISBN: 0-521-64314-7.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -1.051182085412523D+00, &
    -0.2419959664297382D+00, &
    0.001008680985009855D+00, &
    0.08004939780706674D+00, &
    0.09331378813535750D+00, &
   0.08027022252392219D+00, &
   0.05937625647622691D+00, &
   0.03916601076917133D+00, &
   0.02300216024690250D+00, &
   0.01157775439325247D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.5D+00, &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    2.5D+00, &
    3.0D+00, &
    3.5D+00, &
    4.0D+00, &
    4.5D+00, &
    5.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine ker0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! KER0_VALUES returns some values of the Kelvin KER function of order NU = 0.
!
!  Discussion:
!
!    The function is defined by:
!
!      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
!
!    where K(NU,X) is the K Bessel function.
!
!    In Mathematica, KER(NU,X) can be defined by:
!
!      Re [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    LC: QA76.95.W65,
!    ISBN: 0-521-64314-7.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8559058721186342D+00, &
    0.2867062087283160D+00, &
    0.05293491548771044D+00, &
   -0.04166451399150953D+00, &
   -0.06968797258904534D+00, &
   -0.06702923330379870D+00, &
   -0.05263927724224119D+00, &
   -0.03617884789954761D+00, &
   -0.02199987504667382D+00, &
   -0.01151172719949066D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.5D+00, &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    2.5D+00, &
    3.0D+00, &
    3.5D+00, &
    4.0D+00, &
    4.5D+00, &
    5.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine ker1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! KER1_VALUES returns some values of the Kelvin KER function of order NU = 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
!
!    where K(NU,X) is the K Bessel function.
!
!    In Mathematica, KER(NU,X) can be defined by:
!
!      Re [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    LC: QA76.95.W65,
!    ISBN: 0-521-64314-7.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
   -1.522403406532090D+00, &
   -0.7403222768419827D+00, &
   -0.4170442851662574D+00, &
   -0.2308059295181230D+00, &
   -0.1172561358598705D+00, &
   -0.04989830778751491D+00, &
   -0.01272324936181659D+00, &
    0.005351296460277448D+00, &
    0.01209090413515866D+00, &
    0.01273739048421857D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.5D+00, &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    2.5D+00, &
    3.0D+00, &
    3.5D+00, &
    4.0D+00, &
    4.5D+00, &
    5.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine laguerre_associated_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! LAGUERRE_ASSOCIATED_VALUES returns values of associated Laguerre polynomials.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,m,x]
!
!    The associated Laguerre polynomials may be generalized so that the
!    parameter M is allowed to take on arbitrary noninteger values.
!    The resulting function is known as the generalized Laguerre function.
!
!    The polynomials satisfy the differential equation:
!
!      X * Y'' + (M+1-X) * Y' + (N-M) * Y = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, integer ( kind = 4 ) M, the parameter.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.1625000000000000D+01, &
    0.1479166666666667D+01, &
    0.1148437500000000D+01, &
    0.4586666666666667D+00, &
    0.2878666666666667D+01, &
    0.8098666666666667D+01, &
    0.1711866666666667D+02, &
    0.1045328776041667D+02, &
    0.1329019368489583D+02, &
    0.5622453647189670D+02, &
    0.7484729341779436D+02, &
    0.3238912982762806D+03, &
    0.4426100000097533D+03, &
    0.1936876572288250D+04 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 0, 0, &
    0, 1, 1, 1, &
    1, 0, 1, 2, &
    3, 2, 2, 3, &
    3, 4, 4, 5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    1,  2,  3,  4, &
    5,  1,  2,  3, &
    4,  3,  3,  3, &
    3,  4,  5,  6, &
    7,  8,  9, 10 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine laguerre_general_values ( n_data, n, a, x, fx )

!*****************************************************************************80
!
!! LAGUERRE_GENERAL_VALUES returns values of the generalized Laguerre function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,a,x]
!
!    The functions satisfy the following differential equation:
!
!      X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0
!
!    Function values can be generated by the recursion:
!
!      L(0,ALPHA)(X) = 1
!      L(1,ALPHA)(X) = 1+ALPHA-X
!
!      L(N,ALPHA)(X) = ( (2*N-1+ALPHA-X) * L(N-1,ALPHA)(X)
!                     - (N-1+ALPHA) * L(N-2,ALPHA)(X) ) / N
!
!    The parameter ALPHA is required to be greater than -1.
!
!    For ALPHA = 0, the generalized Laguerre function L(N,ALPHA)(X)
!    is equal to the Laguerre polynomial L(N)(X).
!
!    For ALPHA integral, the generalized Laguerre function
!    L(N,ALPHA)(X) equals the associated Laguerre polynomial L(N,ALPHA)(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) A, the parameter.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.00D+00, &
    0.25D+00, &
    0.50D+00, &
    0.75D+00, &
    1.50D+00, &
    2.50D+00, &
    5.00D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.3726399739583333D-01, &
     0.3494791666666667D+00, &
     0.8710042317708333D+00, &
     0.1672395833333333D+01, &
     0.6657625325520833D+01, &
     0.2395726725260417D+02, &
     0.2031344319661458D+03, &
     0.1284193996800000D+02, &
     0.5359924801587302D+01, &
     0.9204589064126984D+00, &
    -0.1341585114857143D+01, &
    -0.2119726307555556D+01, &
    -0.1959193658349206D+01, &
     0.1000000000000000D+01, &
     0.5450000000000000D+01, &
     0.1720125000000000D+02, &
     0.4110393750000000D+02, &
     0.8239745859375000D+02, &
     0.1460179186171875D+03, &
     0.2359204608298828D+03 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     8, &
     8, &
     8, &
     8, &
     8, &
     8, &
     0, &
     1, &
     2, &
     3, &
     4, &
     5, &
     6 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.00D+00, &
    0.20D+00, &
    0.40D+00, &
    0.60D+00, &
    0.80D+00, &
    1.00D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine laguerre_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_VALUES: some values of the Laguerre polynomial L(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,x]
!
!  Differential equation:
!
!    X * Y'' + (1-X) * Y' + N * Y = 0
!
!  First terms:
!
!      1
!     -X    +  1
!   (  X^2 -  4 X     +  2 ) / 2
!   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
!   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
!   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
!   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
!   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3 + 52920 X^2 - 35280 X 
!     + 5040 ) / 5040
!
!  Recursion:
!
!    L(0)(X) = 1,
!    L(1)(X) = 1-X,
!    N * L(N)(X) = (2*N-1-X) * L(N-1)(X) - (N-1) * L(N-2)(X)
!
!  Orthogonality:
!
!    integral ( 0 <= X < oo ) exp ( - X ) * L(N)(X) * L(M)(X) dX
!    = 0 if N /= M
!    = 1 if N == M
!
!  Special values:
!
!    L(N)(0) = 1.
!
!  Relations:
!
!    L(N)(X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.0000000000000000D+00, &
    -0.5000000000000000D+00, &
    -0.6666666666666667D+00, &
    -0.6250000000000000D+00, &
    -0.4666666666666667D+00, &
    -0.2569444444444444D+00, &
    -0.4047619047619048D-01, &
     0.1539930555555556D+00, &
     0.3097442680776014D+00, &
     0.4189459325396825D+00, &
     0.4801341790925124D+00, &
     0.4962122235082305D+00, &
    -0.4455729166666667D+00, &
     0.8500000000000000D+00, &
    -0.3166666666666667D+01, &
     0.3433333333333333D+02 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    0.5D+00, &
    3.0D+00, &
    5.0D+00, &
    1.0D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine lambert_w_values ( n_data, x, fx )

!*****************************************************************************80
!
!! LAMBERT_W_VALUES returns some values of the Lambert W function.
!
!  Discussion:
!
!    The function W(X) is defined implicitly by:
!
!      W(X) * e^W(X) = X
!
!    The function is also known as the "Omega" function.
!
!    In Mathematica, the function can be evaluated by:
!
!      W = ProductLog [ X ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Corless, Gaston Gonnet, David Hare, David Jeffrey,
!    Donald Knuth,
!    On the Lambert W Function,
!    Advances in Computational Mathematics,
!    Volume 5, 1996, pages 329-359.
!
!    Brian Hayes,
!    "Why W?",
!    The American Scientist,
!    Volume 93, March-April 2005, pages 104-108.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.3517337112491958D+00, &
    0.5671432904097839D+00, &
    0.7258613577662263D+00, &
    0.8526055020137255D+00, &
    0.9585863567287029D+00, &
    0.1000000000000000D+01, &
    0.1049908894964040D+01, &
    0.1130289326974136D+01, &
    0.1202167873197043D+01, &
    0.1267237814307435D+01, &
    0.1326724665242200D+01, &
    0.1381545379445041D+01, &
    0.1432404775898300D+01, &
    0.1479856830173851D+01, &
    0.1524345204984144D+01, &
    0.1566230953782388D+01, &
    0.1605811996320178D+01, &
    0.1745528002740699D+01, &
    0.3385630140290050D+01, &
    0.5249602852401596D+01, &
    0.1138335808614005D+02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.2000000000000000D+01, &
    0.2500000000000000D+01, &
    0.2718281828459045D+01, &
    0.3000000000000000D+01, &
    0.3500000000000000D+01, &
    0.4000000000000000D+01, &
    0.4500000000000000D+01, &
    0.5000000000000000D+01, &
    0.5500000000000000D+01, &
    0.6000000000000000D+01, &
    0.6500000000000000D+01, &
    0.7000000000000000D+01, &
    0.7500000000000000D+01, &
    0.8000000000000000D+01, &
    0.1000000000000000D+02, &
    0.1000000000000000D+03, &
    0.1000000000000000D+04, &
    0.1000000000000000D+07 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine laplace_cdf_values ( n_data, mu, beta, x, fx )

!*****************************************************************************80
!
!! LAPLACE_CDF_VALUES returns some values of the Laplace CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = LaplaceDistribution [ mu, beta ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) BETA, the shape parameter.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) beta
  real ( kind = 8 ), save, dimension ( n_max ) :: beta_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.8160602794142788D+00, &
    0.9323323583816937D+00, &
    0.9751064658160680D+00, &
    0.6967346701436833D+00, &
    0.6417343447131054D+00, &
    0.6105996084642976D+00, &
    0.5906346234610091D+00, &
    0.5000000000000000D+00, &
    0.3032653298563167D+00, &
    0.1839397205857212D+00, &
    0.1115650800742149D+00 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    beta = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    beta = beta_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_associated_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! LEGENDRE_ASSOCIATED_VALUES returns values of associated Legendre functions.
!
!  Discussion:
!
!    The function considered is the associated Legendre polynomial 
!    P^M_N(X) or Pm(n,m,x).
!
!    In Mathematica, the function can be evaluated by:
!
!      LegendreP [ n, m, x ]
!
!  Differential equation:
!
!    (1-X^2) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X^2)) * Y = 0
!
!  First terms:
!
!    M = 0  ( = Legendre polynomials of first kind P(N,X) )
!
!    Pm(0,0) =    1
!    Pm(1,0) =    1 X
!    Pm(2,0) = (  3 X^2 -   1)/2
!    Pm(3,0) = (  5 X^3 -   3 X)/2
!    Pm(4,0) = ( 35 X^4 -  30 X^2 +   3)/8
!    Pm(5,0) = ( 63 X^5 -  70 X^3 +  15 X)/8
!    Pm(6,0) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
!    Pm(7,0) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
!
!    M = 1
!
!    Pm(0,1) =   0
!    Pm(1,1) =   1 * SQRT(1-X^2)
!    Pm(2,1) =   3 * SQRT(1-X^2) * X
!    Pm(3,1) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
!    Pm(4,1) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)
!
!    M = 2
!
!    Pm(0,2) =   0
!    Pm(1,2) =   0
!    Pm(2,2) =   3 * (1-X^2)
!    Pm(3,2) =  15 * (1-X^2) * X
!    Pm(4,2) = 7.5 * (1-X^2) * (7*X^2-1)
!
!    M = 3
!
!    Pm(0,3) =   0
!    Pm(1,3) =   0
!    Pm(2,3) =   0
!    Pm(3,3) =  15 * (1-X^2)^1.5
!    Pm(4,3) = 105 * (1-X^2)^1.5 * X
!
!    M = 4
!
!    Pm(0,4) =   0
!    Pm(1,4) =   0
!    Pm(2,4) =   0
!    Pm(3,4) =   0
!    Pm(4,4) = 105 * (1-X^2)^2
!
!  Recursion:
!
!    if N < M:
!      Pm(N,M) = 0
!    if N = M:
!      Pm(N,M) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
!      all the odd integers less than or equal to N.
!    if N = M+1:
!      Pm(N,M) = X*(2*M+1)*Pm(M,M)
!    if M+1 < N:
!      Pm(N,M) = ( X*(2*N-1)*Pm(N-1,M) - (N+M-1)*Pm(N-2,M) )/(N-M)
!
!  Restrictions:
!
!    -1 <= X <= 1
!     0 <= M <= N
!
!  Special values:
!
!    Pm(N,0)(X) = P(N,X), that is, for M=0, the associated Legendre
!    polynomial of the first kind equals the Legendre polynomial of the
!    first kind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) X, the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
    -0.5000000000000000D+00, &
     0.0000000000000000D+00, &
     0.3750000000000000D+00, &
     0.0000000000000000D+00, &
    -0.8660254037844386D+00, &
    -0.1299038105676658D+01, &
    -0.3247595264191645D+00, &
     0.1353164693413185D+01, &
    -0.2800000000000000D+00, &
     0.1175755076535925D+01, &
     0.2880000000000000D+01, &
    -0.1410906091843111D+02, &
    -0.3955078125000000D+01, &
    -0.9997558593750000D+01, &
     0.8265311444100484D+02, &
     0.2024442836815152D+02, &
    -0.4237997531890869D+03, &
     0.1638320624828339D+04, &
    -0.2025687389227225D+05 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 0, 0, &
    0, 1, 1, 1, &
    1, 0, 1, 2, &
    3, 2, 2, 3, &
    3, 4, 4, 5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    1,  2,  3,  4, &
    5,  1,  2,  3, &
    4,  3,  3,  3, &
    3,  4,  5,  6, &
    7,  8,  9, 10 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_associated_normalized_sphere_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES: normalized associated Legendre.
!
!  Discussion:
!
!    The function considered is the associated Legendre polynomial P^M_N(X).
!
!    In Mathematica, the function can be evaluated by:
!
!      LegendreP [ n, m, x ]
!
!    The function is normalized by dividing by
!
!      sqrt ( 4 * pi * ( n + m )! / ( 4 * pi * n + 1 ) / ( n - m )! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) X, the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.2820947917738781D+00, &
     0.2443012559514600D+00, &
    -0.2992067103010745D+00, &
    -0.07884789131313000D+00, &
    -0.3345232717786446D+00, &
     0.2897056515173922D+00, &
    -0.3265292910163510D+00, &
    -0.06997056236064664D+00, &
     0.3832445536624809D+00, &
    -0.2709948227475519D+00, &
    -0.2446290772414100D+00, &
     0.2560660384200185D+00, &
     0.1881693403754876D+00, &
    -0.4064922341213279D+00, &
     0.2489246395003027D+00, &
     0.08405804426339821D+00, &
     0.3293793022891428D+00, &
    -0.1588847984307093D+00, &
    -0.2808712959945307D+00, &
     0.4127948151484925D+00, &
    -0.2260970318780046D+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 1, 0, &
    1, 2, 0, 1, &
    2, 3, 0, 1, &
    2, 3, 4, 0, &
    1, 2, 3, 4, &
    5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  1,  2, &
    2,  2,  3,  3, &
    3,  3,  4,  4, &
    4,  4,  4,  5, &
    5,  5,  5,  5, &
    5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_associated_normalized_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! LEGENDRE_ASSOCIATED_NORMALIZED_VALUES: normalized associated Legendre.
!
!  Discussion:
!
!    The function considered is the associated Legendre polynomial P^M_N(X).
!
!    In Mathematica, the function can be evaluated by:
!
!      LegendreP [ n, m, x ]
!
!    The function is normalized by dividing by
!
!      sqrt ( 4 * pi * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) X, the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7071067811865475D+00, &
    0.6123724356957945D+00, &
   -0.7500000000000000D+00, &
   -0.1976423537605237D+00, &
   -0.8385254915624211D+00, &
    0.7261843774138907D+00, &
   -0.8184875533567997D+00, &
   -0.1753901900050285D+00, &
    0.9606516343087123D+00, &
   -0.6792832849776299D+00, &
   -0.6131941618102092D+00, &
    0.6418623720763665D+00, &
    0.4716705890038619D+00, &
   -0.1018924927466445D+01, &
    0.6239615396237876D+00, &
    0.2107022704608181D+00, &
    0.8256314721961969D+00, &
   -0.3982651281554632D+00, &
   -0.7040399320721435D+00, &
    0.1034723155272289D+01, &
   -0.5667412129155530D+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 1, 0, &
    1, 2, 0, 1, &
    2, 3, 0, 1, &
    2, 3, 4, 0, &
    1, 2, 3, 4, &
    5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  1,  2, &
    2,  2,  3,  3, &
    3,  3,  4,  4, &
    4,  4,  4,  5, &
    5,  5,  5,  5, &
    5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_function_q_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! LEGENDRE_FUNCTION_Q_VALUES returns values of the Legendre Q function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LegendreQ[n,x]
!
!  Differential equation:
!
!    (1-X*X) Y'' - 2 X Y' + N (N+1) = 0
!
!  First terms:
!
!    Qf(0,X) = 0.5 * log((1+X)/(1-X))
!    Qf(1,X) = Qf(0,X) * X - 1
!    Qf(2,X) = Qf(0,X) * (3*X^2-1)/4 - 1.5*X
!    Qf(3,X) = Qf(0,X) * (5*X^3-3*X)/4 - 2.5*X^2 + 2/3
!    Qf(4,X) = Qf(0,X) * (35*X^4-30*X62+3)/16 - 35/8 * X^3 + 55/24 * X
!    Qf(5,X) = Qf(0,X) * (63*X^5-70*X^3+15*X)/16 - 63/8*X^4 + 49/8*X^2 - 8/15
!
!  Recursion:
!
!    Qf(0,X) = 0.5 * log ( (1+X) / (1-X) )
!    Qf(1,X) = 0.5 * X * log ( (1+X) / (1-X) ) - 1.0
!
!    Qf(N,X) = ( (2*N-1) * X * Qf(N-1,X) - (N-1) * Qf(N-2,X) ) / N
!
!  Restrictions:
!
!    -1 < X < 1
!
!  Special values:
!
!    Note that the Legendre function Qf(N,X) is equal to the
!    associated Legendre function of the second kind,
!    Q(N,M)(X) with M = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.2554128118829953D+00, &
    -0.9361467970292512D+00, &
    -0.4787614548274669D+00, &
     0.4246139251747229D+00, &
     0.5448396833845414D+00, &
    -0.9451328261673470D-01, &
    -0.4973516573531213D+00, &
    -0.1499018843853194D+00, &
     0.3649161918783626D+00, &
     0.3055676545072885D+00, &
    -0.1832799367995643D+00, &
     0.6666666666666667D+00, &
     0.6268672028763330D+00, &
     0.5099015515315237D+00, &
     0.3232754180589764D+00, &
     0.8026113738148187D-01, &
    -0.1986547714794823D+00, &
    -0.4828663183349136D+00, &
    -0.7252886849144386D+00, &
    -0.8454443502398846D+00, &
    -0.6627096245052618D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3,  3,  3 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! LEGENDRE_POLY_VALUES returns values of the Legendre polynomials.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LegendreP [ n, x ]
!
!  Differential equation:
!
!    (1-X*X) * P(N,X)'' - 2 * X * P(N,X)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0,X) =       1
!    P( 1,X) =       1 X
!    P( 2,X) =  (    3 X^2 -       1)/2
!    P( 3,X) =  (    5 X^3 -     3 X)/2
!    P( 4,X) =  (   35 X^4 -    30 X^2 +     3)/8
!    P( 5,X) =  (   63 X^5 -    70 X^3 +    15 X)/8
!    P( 6,X) =  (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
!    P( 7,X) =  (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
!    P( 8,X) =  ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
!    P( 9,X) =  (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
!    P(10,X) =  (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2
!                 -63 ) /256
!
!  Recursion:
!
!    P(0,X) = 1
!    P(1,X) = X
!    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!
!    P'(0,X) = 0
!    P'(1,X) = 1
!    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!
!  Formula:
!
!    P(N)(X) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
!
!  Orthogonality:
!
!    integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
!      = 0 if I =/= J
!      = 2 / ( 2*I+1 ) if I = J.
!
!  Approximation:
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!
!      C0*P(0,X) + C1*P(1,X) + ... + CN*P(N,X)
!
!    where
!
!      C(I) = (2*I+1)/(2) * integral ( -1 <= X <= 1 ) F(X) P(I,X) dx.
!
!  Special values:
!
!    P(N,1) = 1.
!    P(N,-1) = (-1)^N.
!    | P(N,X) | <= 1 in [-1,1].
!
!    Pm(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
!    function of the first kind and order N equals the Legendre polynomial
!    of the first kind and order N.
!
!    The N zeroes of P(N,X) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.2500000000000000D+00, &
    -0.4062500000000000D+00, &
    -0.3359375000000000D+00, &
     0.1577148437500000D+00, &
     0.3397216796875000D+00, &
     0.2427673339843750D-01, &
    -0.2799186706542969D+00, &
    -0.1524540185928345D+00, &
     0.1768244206905365D+00, &
     0.2212002165615559D+00, &
     0.0000000000000000D+00, &
    -0.1475000000000000D+00, &
    -0.2800000000000000D+00, &
    -0.3825000000000000D+00, &
    -0.4400000000000000D+00, &
    -0.4375000000000000D+00, &
    -0.3600000000000000D+00, &
    -0.1925000000000000D+00, &
     0.8000000000000000D-01, &
     0.4725000000000000D+00, &
     0.1000000000000000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine lerch_values ( n_data, z, s, a, fx )

!*****************************************************************************80
!
!! LERCH_VALUES returns some values of the Lerch transcendent function.
!
!  Discussion:
!
!    The Lerch function is defined as
!
!      Phi(z,s,a) = Sum ( 0 <= k < +oo ) z^k / ( a + k )^s
!
!    omitting any terms with ( a + k ) = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      LerchPhi[z,s,a]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) Z, the parameters of the function.
!
!    Output, integer ( kind = 4 ) S, the parameters of the function.
!
!    Output, real ( kind = 8 ) A, the parameters of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    3.0D+00, &
    3.0D+00, &
    3.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1644934066848226D+01, &
    0.1202056903159594D+01, &
    0.1000994575127818D+01, &
    0.1164481052930025D+01, &
    0.1074426387216080D+01, &
    0.1000492641212014D+01, &
    0.2959190697935714D+00, &
    0.1394507503935608D+00, &
    0.9823175058446061D-03, &
    0.1177910993911311D+00, &
    0.3868447922298962D-01, &
    0.1703149614186634D-04 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) s
  integer ( kind = 4 ), save, dimension ( n_max ) :: s_vec = (/ &
     2, 3, 10, &
     2, 3, 10, &
     2, 3, 10, &
     2, 3, 10 /)
  real ( kind = 8 ) z
  real ( kind = 8 ), save, dimension ( n_max ) :: z_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.3333333333333333D+00, &
    0.3333333333333333D+00, &
    0.3333333333333333D+00, &
    0.1000000000000000D+00, &
    0.1000000000000000D+00, &
    0.1000000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    z = 0.0D+00
    s = 0
    a = 0.0D+00
    fx = 0.0D+00
  else
    z = z_vec(n_data)
    s = s_vec(n_data)
    a = a_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine linear_system_values ( n_data, nrow, ncol, nsys, a, x, b  )

!*****************************************************************************80
!
!! LINEAR_SYSTEM_VALUES returns some linear systems.
!
!  Discussion:
!
!    Each call to this routine returns scalars NROW, NCOL and NSYS,
!    which give the dimensions of the linear system
!
!      A(NROW,NCOL) * X(NCOL,NSYS) = B(NROW,NSYS)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns.
!
!    Output, integer ( kind = 4 ) NSYS, the number of systems.
!
!    Output, real ( kind = 8 ) A(NROW,NCOL), the matrix.
!
!    Output, real ( kind = 8 ) X(NCOL,NSYS), the solutions of the linear system.
!
!    Output, real ( kind = 8 ) B(NROW,NSYS), the right hand sides.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nsys

  real ( kind = 8 ) a(nrow,ncol)
  real ( kind = 8 ) b(nrow,nsys)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), parameter :: n_max = 4
  real ( kind = 8 ) x(ncol,nsys)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then

    n_data = 0
    nrow = 0
    ncol = 0
    nsys = 0

  else if ( n_data == 1 ) then

    nrow = 3
    ncol = 3
    nsys = 2

    a(1:nrow,1:ncol) = reshape ( (/ &
      1.0, 0.0, 0.0, &
      0.0, 2.0, 0.0, &
      0.0, 0.0, 3.0 /), (/ nrow, ncol /) )
    x(1:ncol,1:nsys) = reshape ( (/ &
      1.0, 0.0, 0.0, &
      1.0, 1.0, 1.0 /), (/ ncol, nsys /) )
    b(1:nrow,1:nsys) = reshape ( (/ &
      1.0, 0.0, 0.0, &
      1.0, 2.0, 3.0 /), (/ nrow, nsys /) )

  else if ( n_data == 2 ) then

    nrow = 3
    ncol = 3
    nsys = 2

    a(1:nrow,1:ncol) = reshape ( (/ &
      1.0, 2.0, 3.0, &
      2.0, 2.0, 3.0, &
      3.0, 3.0, 3.0 /), (/ nrow, ncol /) )
    x(1:ncol,1:nsys) = reshape ( (/ &
      1.0, 1.0, 1.0, &
      1.0, 2.0, 3.0 /), (/ ncol, nsys /) )
    b(1:nrow,1:nsys) = reshape ( (/ &
       6.0,  7.0,  9.0, &
      14.0, 15.0, 18.0 /), (/ nrow, nsys /) )

  else if ( n_data == 3 ) then

    nrow = 5
    ncol = 5
    nsys = 2

    a(1:nrow,1:ncol) = reshape ( (/ &
      1.0, 2.0, 3.0, 4.0, 5.0, &
      2.0, 3.0, 4.0, 5.0, 1.0, &
      3.0, 4.0, 5.0, 1.0, 2.0, &
      4.0, 5.0, 1.0, 2.0, 3.0, &
      5.0, 1.0, 2.0, 3.0, 4.0 /), (/ nrow, ncol /) )
    x(1:nrow,1:nsys) = reshape ( (/ &
      0.066667,  0.066667,  0.066667, 0.066667, 0.066667, &
      1.0,       0.0,       0.0,      0.0,      0.0 /), (/ nrow, nsys /) )
    b(1:ncol,1:nsys) = reshape ( (/ &
      1.0, 1.0, 1.0, 1.0, 1.0, &
      1.0, 2.0, 3.0, 4.0, 5.0 /), (/ ncol, nsys /) )

  else if ( n_data == 4 ) then

    nrow = 5
    ncol = 5
    nsys = 2

    a(1:nrow,1:ncol) = reshape ( (/ &
      1.4, 1.6, 3.8, 4.6, 2.6, &
      2.1, 1.5, 8.0, 8.2, 2.9, &
      2.1, 1.1, 9.6, 8.4, 0.1, &
      7.4, 0.7, 5.4, 0.4, 9.6, &
      9.6, 5.0, 8.8, 8.0, 7.7 /), (/ nrow, ncol /) )
    x(1:nrow,1:nsys) = reshape ( (/ &
      -5.313077,   5.735670, -2.507606, -1.058741,  0.999381, &
      31.601006, -28.594793, 13.389395,  2.780322, -3.008797 /), &
      (/ nrow, nsys /) )
    b(1:ncol,1:nsys) = reshape ( (/ &
      1.1, 1.6, 4.7, 9.1, 0.1, &
      4.0, 9.3, 8.4, 0.4, 4.1 /), (/ ncol, nsys /) )

  end if

  return
end
subroutine lobachevsky_values ( n_data, x, fx )

!*****************************************************************************80
!
!! LOBACHEVSKY_VALUES returns some values of the Lobachevsky function.
!
!  Discussion:
!
!    The function is defined by:
!
!      LOBACHEVSKY(x) = integral ( 0 <= t <= x ) -ln ( abs ( cos ( t ) ) dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.12417639065161393857D-08, &
    0.79473344770001088225D-07, &
    0.50867598186208834198D-05, &
    0.32603097901207200319D-03, &
    0.21380536815408214419D-01, &
    0.18753816902083824050D+00, &
    0.83051199971883645115D+00, &
    0.18854362426679034904D+01, &
    0.21315988986516411053D+01, &
    0.21771120185613427221D+01, &
    0.22921027921896650849D+01, &
    0.39137195028784495586D+01, &
    0.43513563983836427904D+01, &
    0.44200644968478185898D+01, &
    0.65656013133623829156D+01, &
    0.10825504661504599479D+02, &
    0.13365512855474227325D+02, &
    0.21131002685639959927D+02, &
    0.34838236589449117389D+02, &
    0.69657062437837394278D+02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0078125000D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      5.0000000000D+00, &
      6.0000000000D+00, &
      7.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00, &
    100.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! LOG_VALUES returns some values of the natural logarithm function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Log[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -11.512925464970228420D+00, &
     -4.6051701859880913680D+00, &
     -2.3025850929940456840D+00, &
     -1.6094379124341003746D+00, &
     -1.2039728043259359926D+00, &
     -0.91629073187415506518D+00, &
     -0.69314718055994530942D+00, &
     -0.51082562376599068321D+00, &
     -0.35667494393873237891D+00, &
     -0.22314355131420975577D+00, &
     -0.10536051565782630123D+00, &
      0.00000000000000000000D+00, &
      0.69314718055994530942D+00, &
      1.0986122886681096914D+00, &
      1.1447298858494001741D+00, &
      1.6094379124341003746D+00, &
      2.3025850929940456840D+00, &
      2.9957322735539909934D+00, &
      4.6051701859880913680D+00, &
      18.631401766168018033D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D-05, &
    1.0D-02, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    3.1415926535897932385D+00, &
    5.0D+00, &
    10.0D+00, &
    20.0D+00, &
    100.0D+00, &
    123456789.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine log_normal_cdf_values ( n_data, mu, sigma, x, fx )

!*****************************************************************************80
!
!! LOG_NORMAL_CDF_VALUES returns some values of the Log Normal CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = LogNormalDistribution [ mu, sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) SIGMA, the shape parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.2275013194817921D-01, &
    0.2697049307349095D+00, &
    0.5781741008028732D+00, &
    0.7801170895122241D+00, &
    0.4390310097476894D+00, &
    0.4592655190218048D+00, &
    0.4694258497695908D+00, &
    0.4755320473858733D+00, &
    0.3261051056816658D+00, &
    0.1708799040927608D+00, &
    0.7343256357952060D-01, &
    0.2554673736161761D-01 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine log_series_cdf_values ( n_data, t, n, fx )

!*****************************************************************************80
!
!! LOG_SERIES_CDF_VALUES returns some values of the log series CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = LogSeriesDistribution [ t ]
!      CDF [ dist, n ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) T, the parameter of the function.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 29

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.9491221581029903D+00, &
    0.9433541128559735D+00, &
    0.9361094611773272D+00, &
    0.9267370278044118D+00, &
    0.9141358246245129D+00, &
    0.8962840235449100D+00, &
    0.8690148741955517D+00, &
    0.8221011541254772D+00, &
    0.7213475204444817D+00, &
    0.6068261510845583D+00, &
    0.5410106403333613D+00, &
    0.4970679476476894D+00, &
    0.4650921887927060D+00, &
    0.4404842934597863D+00, &
    0.4207860535926143D+00, &
    0.4045507673897055D+00, &
    0.3908650337129266D+00, &
    0.2149757685421097D+00, &
    0.0000000000000000D+00, &
    0.2149757685421097D+00, &
    0.3213887739704539D+00, &
    0.3916213575531612D+00, &
    0.4437690508633213D+00, &
    0.4850700239649681D+00, &
    0.5191433267738267D+00, &
    0.5480569580144867D+00, &
    0.5731033910767085D+00, &
    0.5951442521714636D+00, &
    0.6147826594068904D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     1, 1, 1, 1, 1, &
     1, 1, 1, 1, 1, &
     1, 1, 1, 1, 1, &
     1, 1, 1, 0, 1, &
     2, 3, 4, 5, 6, &
     7, 8, 9, 10 /)
  real ( kind = 8 ) t
  real ( kind = 8 ), save, dimension ( n_max ) :: t_vec = (/ &
    0.1000000000000000D+00, &
    0.1111111111111111D+00, &
    0.1250000000000000D+00, &
    0.1428571428571429D+00, &
    0.1666666666666667D+00, &
    0.2000000000000000D+00, &
    0.2500000000000000D+00, &
    0.3333333333333333D+00, &
    0.5000000000000000D+00, &
    0.6666666666666667D+00, &
    0.7500000000000000D+00, &
    0.8000000000000000D+00, &
    0.8333333333333333D+00, &
    0.8571485714857149D+00, &
    0.8750000000000000D+00, &
    0.8888888888888889D+00, &
    0.9000000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    t = 0.0D+00
    n = 0
    fx = 0.0D+00
  else
    t = t_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine log10_values ( n_data, x, fx )

!*****************************************************************************80
!
!! LOG10_VALUES returns some values of the logarithm 10 function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Log10[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer, parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
   -5.0000000000000000000D+00, &
   -2.0000000000000000000D+00, &
   -1.0000000000000000000D+00, &
   -0.69897000433601880479D+00, &
   -0.52287874528033756270D+00, &
   -0.39794000867203760957D+00, &
   -0.30102999566398119521D+00, &
   -0.22184874961635636749D+00, &
   -0.15490195998574316929D+00, &
   -0.096910013008056414359D+00, &
   -0.045757490560675125410D+00, &
    0.000000000000000000000D+00, &
    0.30102999566398119521D+00, &
    0.47712125471966243730D+00, &
    0.49714987269413385435D+00, &
    0.69897000433601880479D+00, &
    1.0000000000000000000D+00, &
    1.3010299956639811952D+00, &
    2.0000000000000000000D+00, &
    8.0915149771692704475D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
   1.0D-05,  &
   1.0D-02,  &
   0.1D+00,  &
   0.2D+00,  &
   0.3D+00,  &
   0.4D+00,  &
   0.5D+00,  &
   0.6D+00,  &
   0.7D+00,  &
   0.8D+00,  &
   0.9D+00,  &
   1.0D+00,  &
   2.0D+00,  &
   3.0D+00,  &
   3.1415926535897932385D+00,  &
   5.0D+00,  &
   10.0D+00,  &
   20.0D+00, &
   100.0D+00, &
   123456789.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine logarithmic_integral_values ( n_data, x, fx )

!*****************************************************************************80
!
!! LOGARITHMIC_INTEGRAL_VALUES returns values of the logarithmic integral LI(X).
!
!  Discussion:
!
!    The logarithmic integral is defined as:
!
!      LI(X) = integral ( 0 <= T <= Z ) dT / log ( T )
!
!    The principal value of the integral is taken.  There is a
!    branch cut discontinuity in the complex plane from -oo to +1.
!
!    Abramowitz and Stegun assume 1 < X.
!
!    In Mathematica, the function can be evaluated by:
!
!      LogIntegral[x]
!
!    There is a simple relationship with the exponential integral EI:
!
!      LI(X) = EI(LN(X))
!
!    The function LI(X) provides a good approximation to PI(X),
!    the number of primes less than or equal to X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
    -0.3238978959329102D-01, &
    -0.8512648672879405D-01, &
    -0.1574149028946895D+00, &
    -0.2529494192126213D+00, &
    -0.3786710430610880D+00, &
    -0.5468514142104170D+00, &
    -0.7809468775455607D+00, &
    -0.1134011957382327D+01, &
    -0.1775800683423525D+01, &
    -0.2443622553873225D+01, &
    -0.3124190050507211D+01, &
    -0.2872935510329120D+01, &
    -0.2164282524138207D+01, &
    -0.1440351296279408D+01, &
    -0.6864884538258716D+00, &
     0.1250649863152964D+00, &
     0.1045163780117493D+01, &
     0.2967585095039051D+01, &
     0.5253718299558931D+01, &
     0.8519716463711059D+01, &
     0.1360509217709172D+02, &
     0.2193466832805100D+02, &
     0.3604254831722944D+02, &
     0.6051306533791733D+02, &
     0.1037211171690373D+03, &
     0.1810780396816945D+03, &
     0.3211144156746837D+03 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.000000D+00, &
    0.100000D+00, &
    0.200000D+00, &
    0.300000D+00, &
    0.400000D+00, &
    0.500000D+00, &
    0.600000D+00, &
    0.700000D+00, &
    0.800000D+00, &
    0.900000D+00, &
    0.950000D+00, &
    0.975000D+00, &
    0.103125D+01, &
    0.106250D+01, &
    0.112500D+01, &
    0.125000D+01, &
    0.150000D+01, &
    0.200000D+01, &
    0.400000D+01, &
    0.800000D+01, &
    0.160000D+02, &
    0.320000D+02, &
    0.640000D+02, &
    0.128000D+03, &
    0.256000D+03, &
    0.512000D+03, &
    0.102400D+04, &
    0.204800D+04 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine logistic_cdf_values ( n_data, mu, beta, x, fx )

!*****************************************************************************80
!
!! LOGISTIC_CDF_VALUES returns some values of the Logistic CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = LogisticDistribution [ mu, beta ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) BETA, the shape parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) beta
  real ( kind = 8 ), save, dimension ( n_max ) :: beta_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.8807970779778824D+00, &
    0.9820137900379084D+00, &
    0.9975273768433652D+00, &
    0.6224593312018546D+00, &
    0.5825702064623147D+00, &
    0.5621765008857981D+00, &
    0.5498339973124779D+00, &
    0.6224593312018546D+00, &
    0.5000000000000000D+00, &
    0.3775406687981454D+00, &
    0.2689414213699951D+00 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    beta = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    beta = beta_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine mathieu_even_values ( n_data, r, q, a )

!*****************************************************************************80
!
!! MATHIEU_EVEN_VALUES returns eigenvalues of even Mathieu solutions.
!
!  Discussion:
!
!    Mathieu's differential equation can be written
!
!      d2y/dx2 + ( a - 2 * q * cos ( 2 * x ) ) * y = 0
!
!    For each integer Q, there are sets of eigenvalues and
!    associated periodic eigensolutions.  We denote by A(R,Q)
!    the R-th eigenvalue associated with an even periodic
!    solution for Q, and B(R,Q) the R-th eigenvalue associated
!    with an odd periodic solution for Q.
!
!    In Mathematica, the eigenvalues for even functions can
!    be evaluated by
!
!      MathieuCharacteristicA [ r, q ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) R, the index of the eigenvalue.
!
!    Output, integer ( kind = 4 ) Q, the value of the parameter Q.
!
!    Output, real ( kind = 8 ) A, the eigenvalue of the even solution
!    of Mathieu's equation, A(Q,R).
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
      0.0D+00, &
      1.0D+00, &
      4.0D+00, &
    225.0D+00, &
     25.0D+00, &
     25.54997174998161D+00,  &
     27.70376873393928D+00,  &
     31.95782125217288D+00,  &
     36.64498973413284D+00,  &
     40.05019098580771D+00,  &
     -5.800046020851508D+00, &
    -40.25677954656679D+00,  &
    -14.49130142517482D+00,  &
      5.077983197543472D+00, &
    100.5067700246813D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) q
  integer ( kind = 4 ), save, dimension ( n_max ) :: q_vec = (/ &
     0, &
     0, &
     0, &
     0, &
     0, &
     5, &
    10, &
    15, &
    20, &
    25, &
     5, &
    25, &
    20, &
    15, &
    10 /)
  integer ( kind = 4 ) r
  integer ( kind = 4 ), save, dimension ( n_max ) :: r_vec = (/ &
     0, &
     1, &
     2, &
    15, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     0, &
     0, &
     1, &
     2, &
    10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    r = 0
    q = 0
    a = 0.0D+00
  else
    r = r_vec(n_data)
    q = q_vec(n_data)
    a = a_vec(n_data)
  end if

  return
end
subroutine mathieu_odd_values ( n_data, r, q, b )

!*****************************************************************************80
!
!! MATHIEU_ODD_VALUES returns eigenvalues of odd Mathieu solutions.
!
!  Discussion:
!
!    Mathieu's differential equation can be written
!
!      d2y/dx2 + ( a - 2 * q * cos ( 2 * x ) ) * y = 0
!
!    For each integer Q, there are sets of eigenvalues and
!    associated periodic eigensolutions.  We denote by A(R,Q)
!    the R-th eigenvalue associated with an even periodic
!    solution for Q, and B(R,Q) the R-th eigenvalue associated
!    with an odd periodic solution for Q.
!
!    In Mathematica, the eigenvalues for odd functions can
!    be evaluated by
!
!      MathieuCharacteristicB [ r, q ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) R, the index of the eigenvalue.
!
!    Output, integer ( kind = 4 ) Q, the value of the parameter Q.
!
!    Output, real ( kind = 8 ) B, the eigenvalue of the odd solution
!    of Mathieu's equation, B(Q,R).
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
      1.0D+00,               &
    -40.25677898468416D+00,  &
      4.0D+00,               &
      2.099460445486665D+00, &
     -2.382158235956956D+00, &
     -8.099346798895896D+00, &
    -14.49106325598072D+00,  &
    -21.31486062224985D+00,  &
     27.96788059671747D+00,  &
    100.0D+00,               &
    100.5067694628784D+00,   &
    225.0D+00,               &
    225.8951534161768D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) q
  integer ( kind = 4 ), save, dimension ( n_max ) :: q_vec = (/ &
     0, &
    25, &
     0, &
     5, &
    10, &
    15, &
    20, &
    25, &
    15, &
     0, &
    10, &
     0, &
    20 /)
  integer ( kind = 4 ) r
  integer ( kind = 4 ), save, dimension ( n_max ) :: r_vec = (/ &
    1, &
    1, &
    2, &
    2, &
    2, &
    2, &
    2, &
    2, &
    5, &
   10, &
   10, &
   15, &
   15 /)
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    r = 0
    q = 0
    b = 0.0D+00
  else
    r = r_vec(n_data)
    q = q_vec(n_data)
    b = b_vec(n_data)
  end if

  return
end
subroutine mertens_values ( n_data, n, c )

!*****************************************************************************80
!
!! MERTENS_VALUES returns some values of the Mertens function.
!
!  Discussion:
!
!    The Mertens function M(N) is the sum from 1 to N of the Moebius
!    function MU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc Deleglise, Joel Rivat,
!    Computing the Summation of the Moebius Function,
!    Experimental Mathematics,
!    Volume 5, 1996, pages 291-295.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the Mertens function.
!
!    Output, integer ( kind = 4 ) C, the value of the Mertens function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 15

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
      1,   0,  -1,   -1,  -2,  -1,  -2,  -2,   -2,  -1, &
     -2,  -2,   1,    2, -23 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     11,  12,  100, 1000, 10000 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine moebius_values ( n_data, n, c )

!*****************************************************************************80
!
!! MOEBIUS_VALUES returns some values of the Moebius function.
!
!  Discussion:
!
!    MU(N) is defined as follows:
!
!      MU(N) = 1 if N = 1;
!              0 if N is divisible by the square of a prime;
!              (-1)**K, if N is the product of K distinct primes.
!
!    In Mathematica, the function can be evaluated by:
!
!      MoebiusMu[n]
!
!  First values:
!
!     N  MU(N)
!
!     1    1
!     2   -1
!     3   -1
!     4    0
!     5   -1
!     6    1
!     7   -1
!     8    0
!     9    0
!    10    1
!    11   -1
!    12    0
!    13   -1
!    14    1
!    15    1
!    16    0
!    17   -1
!    18    0
!    19   -1
!    20    0
!
!  Note:
!
!    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
!    if N is a square, cube, etc.
!
!  Formula:
!
!    The Moebius function is related to Euler's totient function:
!
!      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the Moebius function.
!
!    Output, integer ( kind = 4 ) C, the value of the Moebius function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
      1,  -1,  -1,   0,  -1,   1,  -1,   0,   0,   1, &
     -1,   0,  -1,   1,   1,   0,  -1,   0,  -1,   0 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     11,  12,  13,  14,  15,  16,  17,  18,  19,  20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine negative_binomial_cdf_values ( n_data, f, s, p, cdf )

!*****************************************************************************80
!
!! NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
!
!  Discussion:
!
!    Assume that a coin has a probability P of coming up heads on
!    any one trial.  Suppose that we plan to flip the coin until we
!    achieve a total of S heads.  If we let F represent the number of
!    tails that occur in this process, then the value of F satisfies
!    a negative binomial PDF:
!
!      PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
!
!    The negative binomial CDF is the probability that there are F or
!    fewer failures upon the attainment of the S-th success.  Thus,
!
!      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = NegativeBinomialDistribution [ s, p ]
!      CDF [ dist, f ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Frank Powell,
!    Statistical Tables for Sociology, Biology and Physical Sciences,
!    Cambridge University Press, 1982,
!    ISBN: 0521284732,
!    LC: QA276.25.S73.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) F, the maximum number of failures.
!
!    Output, integer ( kind = 4 ) S, the number of successes.
!
!    Output, real ( kind = 8 ) P, the probability of a success on one trial.
!
!    Output, real ( kind = 8 ) CDF, the probability of at most F failures
!    before the S-th success.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 27

  real ( kind = 8 ) cdf
  real ( kind = 8 ), save, dimension ( n_max ) :: cdf_vec = (/ &
    0.6367187500000000D+00, &
    0.3632812500000000D+00, &
    0.1445312500000000D+00, &
    0.5000000000000000D+00, &
    0.2265625000000000D+00, &
    0.6250000000000000D-01, &
    0.3437500000000000D+00, &
    0.1093750000000000D+00, &
    0.1562500000000000D-01, &
    0.1792000000000000D+00, &
    0.4096000000000000D-01, &
    0.4096000000000000D-02, &
    0.7047000000000000D-01, &
    0.1093500000000000D-01, &
    0.7290000000000000D-03, &
    0.9861587127990000D+00, &
    0.9149749500510000D+00, &
    0.7471846521450000D+00, &
    0.8499053647030009D+00, &
    0.5497160941090026D+00, &
    0.2662040052146710D+00, &
    0.6513215599000000D+00, &
    0.2639010709000000D+00, &
    0.7019082640000000D-01, &
    0.1000000000000000D+01, &
    0.1990000000000000D-01, &
    0.1000000000000000D-03 /)
  integer ( kind = 4 ) f
  integer ( kind = 4 ), save, dimension ( n_max ) :: f_vec = (/ &
     4,  3,  2, &
     3,  2,  1, &
     2,  1,  0, &
     2,  1,  0, &
     2,  1,  0, &
    11, 10,  9, &
    17, 16, 15, &
     9,  8,  7, &
     2,  1,  0 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.40D+00, &
    0.40D+00, &
    0.40D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D-01, &
    0.10D-01, &
    0.10D-01 /)
  integer ( kind = 4 ) s
  integer ( kind = 4 ), save, dimension ( n_max ) :: s_vec = (/ &
    4, 5, 6, &
    4, 5, 6, &
    4, 5, 6, &
    4, 5, 6, &
    4, 5, 6, &
    1, 2, 3, &
    1, 2, 3, &
    1, 2, 3, &
    0, 1, 2 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    f = 0
    s = 0
    p = 0.0D+00
    cdf = 0.0D+00
  else
    f = f_vec(n_data)
    s = s_vec(n_data)
    p = p_vec(n_data)
    cdf = cdf_vec(n_data)
  end if

  return
end
subroutine nine_j_values ( n_data, j1, j2, j3, j4, j5, j6, j7, j8, j9, fx )

!*****************************************************************************80
!
!! NINE_J_VALUES returns some values of the Wigner 9J function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) J1, J2, J3, J4, J5, J6, J7, J8, J9,
!    the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 9

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0004270039294528318D+00, &
    -0.001228915451058514D+00, &
    -0.0001944260688400887D+00, &
     0.003338419923885592D+00, &
    -0.0007958936865080434D+00, &
    -0.004338208690251972D+00, &
     0.05379143536399187D+00, &
     0.006211299937499411D+00, &
     0.03042903097250921D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) j1
  real ( kind = 8 ), save, dimension ( n_max ) :: j1_vec = (/ &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    1.0D+00, &
    1.5D+00, &
    2.0D+00, &
    0.5D+00, &
    1.0D+00, &
    1.5D+00  /)
  real ( kind = 8 ) j2
  real ( kind = 8 ), save, dimension ( n_max ) :: j2_vec = (/ &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    3.0D+00, &
    3.0D+00, &
    3.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00  /)
  real ( kind = 8 ) j3
  real ( kind = 8 ), save, dimension ( n_max ) :: j3_vec = (/ &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00 /)
  real ( kind = 8 ) j4
  real ( kind = 8 ), save, dimension ( n_max ) :: j4_vec = (/ &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    4.0D+00, &
    4.0D+00, &
    4.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00 /)
  real ( kind = 8 ) j5
  real ( kind = 8 ), save, dimension ( n_max ) :: j5_vec = (/ &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    1.5D+00, &
    1.5D+00, &
    1.5D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00 /)
  real ( kind = 8 ) j6
  real ( kind = 8 ), save, dimension ( n_max ) :: j6_vec = (/ &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    3.0D+00, &
    3.0D+00, &
    3.0D+00, &
    1.5D+00, &
    1.5D+00, &
    1.5D+00 /)
  real ( kind = 8 ) j7
  real ( kind = 8 ), save, dimension ( n_max ) :: j7_vec = (/ &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    3.5D+00, &
    3.5D+00, &
    3.5D+00, &
    1.5D+00, &
    1.5D+00, &
    1.5D+00 /)
  real ( kind = 8 ) j8
  real ( kind = 8 ), save, dimension ( n_max ) :: j8_vec = (/ &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00 /)
  real ( kind = 8 ) j9
  real ( kind = 8 ), save, dimension ( n_max ) :: j9_vec = (/ &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    1.5D+00, &
    1.5D+00, &
    1.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    j1 = 0.0D+00
    j2 = 0.0D+00
    j3 = 0.0D+00
    j4 = 0.0D+00
    j5 = 0.0D+00
    j6 = 0.0D+00
    j7 = 0.0D+00
    j8 = 0.0D+00
    j9 = 0.0D+00
    fx = 0.0D+00
  else
    j1 = j1_vec(n_data)
    j2 = j2_vec(n_data)
    j3 = j3_vec(n_data)
    j4 = j4_vec(n_data)
    j5 = j5_vec(n_data)
    j6 = j6_vec(n_data)
    j7 = j7_vec(n_data)
    j8 = j8_vec(n_data)
    j9 = j9_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine normal_cdf_values ( n_data, mu, sigma, x, fx )

!*****************************************************************************80
!
!! NORMAL_CDF_VALUES returns some values of the Normal CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NormalDistribution [ mu, sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) SIGMA, the variance of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.9772498680518208D+00, &
    0.9999683287581669D+00, &
    0.9999999990134124D+00, &
    0.6914624612740131D+00, &
    0.6305586598182364D+00, &
    0.5987063256829237D+00, &
    0.5792597094391030D+00, &
    0.6914624612740131D+00, &
    0.5000000000000000D+00, &
    0.3085375387259869D+00, &
    0.1586552539314571D+00 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine normal_01_cdf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NormalDistribution [ 0, 1 ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.5398278372770290D+00, &
    0.5792597094391030D+00, &
    0.6179114221889526D+00, &
    0.6554217416103242D+00, &
    0.6914624612740131D+00, &
    0.7257468822499270D+00, &
    0.7580363477769270D+00, &
    0.7881446014166033D+00, &
    0.8159398746532405D+00, &
    0.8413447460685429D+00, &
    0.9331927987311419D+00, &
    0.9772498680518208D+00, &
    0.9937903346742239D+00, &
    0.9986501019683699D+00, &
    0.9997673709209645D+00, &
    0.9999683287581669D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
    0.1000000000000000D+00, &
    0.2000000000000000D+00, &
    0.3000000000000000D+00, &
    0.4000000000000000D+00, &
    0.5000000000000000D+00, &
    0.6000000000000000D+00, &
    0.7000000000000000D+00, &
    0.8000000000000000D+00, &
    0.9000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.2000000000000000D+01, &
    0.2500000000000000D+01, &
    0.3000000000000000D+01, &
    0.3500000000000000D+01, &
    0.4000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine omega_values ( n_data, n, c )

!*****************************************************************************80
!
!! OMEGA_VALUES returns some values of the OMEGA function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by
!
!      Length [ FactorInteger [ n ] ]
!
!  First values:
!
!     N   OMEGA(N)
!
!     1    1
!     2    1
!     3    1
!     4    1
!     5    1
!     6    2
!     7    1
!     8    1
!     9    1
!    10    2
!    11    1
!    12    2
!    13    1
!    14    2
!    15    2
!    16    1
!    17    1
!    18    2
!    19    1
!    20    2
!
!  Formula:
!
!    If N = 1, then
!
!      OMEGA(N) = 1
!
!    else if the prime factorization of N is
!
!      N = P1**E1 * P2**E2 * ... * PM**EM,
!
!    then
!
!      OMEGA(N) = M
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the OMEGA function.
!
!    Output, integer ( kind = 4 ) C, the value of the OMEGA function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 23

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
      1,   1,   1,   1,   1, &
      2,   1,   1,   1,   2, &
      3,   1,   4,   4,   3, &
      1,   5,   2,   2,   1, &
      6,   7,   8 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
           1, &
           2, &
           3, &
           4, &
           5, &
           6, &
           7, &
           8, &
           9, &
          10, &
          30, &
         101, &
         210, &
        1320, &
        1764, &
        2003, &
        2310, &
        2827, &
        8717, &
       12553, &
       30030, &
      510510, &
     9699690 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine owen_values ( n_data, h, a, t )

!*****************************************************************************80
!
!! OWEN_VALUES returns some values of Owen's T function.
!
!  Discussion:
!
!    Owen's T function is useful for computation of the bivariate normal
!    distribution and the distribution of a skewed normal distribution.
!
!    Although it was originally formulated in terms of the bivariate
!    normal function, the function can be defined more directly as
!
!      T(H,A) = 1 / ( 2 * pi ) *
!        integral ( 0 <= X <= A ) e^(-H^2*(1+X^2)/2) / (1+X^2) dX
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mike Patefield, David Tandy,
!    Fast and Accurate Calculation of Owen's T Function,
!    Journal of Statistical Software,
!    Volume 5, Number 5, 2000, pages 1-25.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) H, a parameter.
!
!    Output, real ( kind = 8 ) A, the upper limit of the integral.
!
!    Output, real ( kind = 8 ) T, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.2500000000000000D+00, &
    0.4375000000000000D+00, &
    0.9687500000000000D+00, &
    0.0625000000000000D+00, &
    0.5000000000000000D+00, &
    0.9999975000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.1000000000000000D+02, &
    0.1000000000000000D+03 /)
  real ( kind = 8 ) h
  real ( kind = 8 ), save, dimension ( n_max ) :: h_vec = (/ &
    0.0625000000000000D+00, &
    6.5000000000000000D+00, &
    7.0000000000000000D+00, &
    4.7812500000000000D+00, &
    2.0000000000000000D+00, &
    1.0000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t
  real ( kind = 8 ), save, dimension ( n_max ) :: t_vec = (/ &
    3.8911930234701366D-02, &
    2.0005773048508315D-11, &
    6.3990627193898685D-13, &
    1.0632974804687463D-07, &
    8.6250779855215071D-03, &
    6.6741808978228592D-02, &
    0.4306469112078537D-01, &
    0.6674188216570097D-01, &
    0.7846818699308410D-01, &
    0.7929950474887259D-01, &
    0.6448860284750376D-01, &
    0.1066710629614485D+00, &
    0.1415806036539784D+00, &
    0.1510840430760184D+00, &
    0.7134663382271778D-01, &
    0.1201285306350883D+00, &
    0.1666128410939293D+00, &
    0.1847501847929859D+00, &
    0.7317273327500385D-01, &
    0.1237630544953746D+00, &
    0.1737438887583106D+00, &
    0.1951190307092811D+00, &
    0.7378938035365546D-01, &
    0.1249951430754052D+00, &
    0.1761984774738108D+00, &
    0.1987772386442824D+00, &
    0.2340886964802671D+00, &
    0.2479460829231492D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    h = 0.0D+00
    a = 0.0D+00
    t = 0.0D+00
  else
    h = h_vec(n_data)
    a = a_vec(n_data)
    t = t_vec(n_data)
  end if

  return
end
subroutine partition_count_values ( n_data, n, c )

!*****************************************************************************80
!
!! PARTITION_COUNT_VALUES returns some values of the integer partition count.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  The number of partitions of N is symbolized
!    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
!    following partitions:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!      = 3 + 1 + 1
!      = 2 + 2 + 1
!      = 2 + 1 + 1 + 1
!      = 1 + 1 + 1 + 1 + 1
!
!    In Mathematica, the function can be evaluated by
!
!      PartitionsP[n]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the integer.
!
!    Output, integer ( kind = 4 ) C, the number of partitions of the integer.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
      1, &
      1,   2,   3,   5,   7,  11,  15,  22,  30,  42, &
     56,  77, 101, 135, 176, 231, 297, 385, 490, 627 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  &
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine partition_distinct_count_values ( n_data, n, c )

!*****************************************************************************80
!
!! PARTITION_DISTINCT_COUNT_VALUES returns some values of Q(N).
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  The number of partitions of N is symbolized
!    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
!    following partitions:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!      = 3 + 1 + 1
!      = 2 + 2 + 1
!      = 2 + 1 + 1 + 1
!      = 1 + 1 + 1 + 1 + 1
!
!    However, if we require that each member of the partition
!    be distinct, so that no nonzero summand occurs more than once,
!    we are computing something symbolized by Q(N).
!    The number 5 has Q(N) = 3, because it has the following partitions
!    into distinct parts:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!
!    In Mathematica, the function can be evaluated by
!
!      PartitionsQ[n]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the integer.
!
!    Output, integer ( kind = 4 ) C, the number of partitions of the integer
!    into distinct parts.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
      1, &
      1,   1,   2,   2,   3,   4,   5,   6,   8,  10, &
     12,  15,  18,  22,  27,  32,  38,  46,  54,  64 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  &
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine phi_values ( n_data, n, c )

!*****************************************************************************80
!
!! PHI_VALUES returns some values of the PHI function.
!
!  Discussion:
!
!    PHI(N) is the number of integers between 1 and N which are
!    relatively prime to N.  I and J are relatively prime if they
!    have no common factors.  The function PHI(N) is known as
!    "Euler's totient function".
!
!    By convention, 1 and N are relatively prime.
!
!    In Mathematica, the function can be evaluated by:
!
!      EulerPhi[n]
!
!  First values:
!
!     N  PHI(N)
!
!     1    1
!     2    1
!     3    2
!     4    2
!     5    4
!     6    2
!     7    6
!     8    4
!     9    6
!    10    4
!    11   10
!    12    4
!    13   12
!    14    6
!    15    8
!    16    8
!    17   16
!    18    6
!    19   18
!    20    8
!
!  Formula:
!
!    PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
!
!    PHI(P^K) = P^(K-1) * ( P - 1 ) if P is prime.
!
!    PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
!
!    N = Sum ( D divides N ) PHI(D).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the PHI function.
!
!    Output, integer ( kind = 4 ) C, the value of the PHI function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
      1,   1,   2,   2,   4,   2,   6,   4,   6,   4, &
      8,   8,  16,  20,  16,  40, 148, 200, 200, 648 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     20,  30,  40,  50,  60, 100, 149, 500, 750, 999 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine pi_values ( n_data, n, p )

!*****************************************************************************80
!
!! PI_VALUES returns values of the Pi function.
!
!  Discussion:
!
!    Pi[n] is the number of primes less than or equal to n.
!
!    In Mathematica, the function can be evaluated by:
!
!      PrimePi[n]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument.
!
!    Output, integer ( kind = 4 ) P, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
            10, &
            20, &
            30, &
            40, &
            50, &
            60, &
            70, &
            80, &
            90, &
           100, &
          1000, &
         10000, &
        100000, &
       1000000, &
      10000000, &
     100000000, &
    1000000000 /)
  integer ( kind = 4 ) p
  integer ( kind = 4 ), save, dimension ( n_max ) :: p_vec = (/ &
             4, &
             8, &
            10, &
            12, &
            15, &
            17, &
            19, &
            22, &
            24, &
            25, &
           168, &
          1229, &
          9592, &
         78498, &
        664579, &
       5761455, &
      50847534 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    p = 0
  else
    n = n_vec(n_data)
    p = p_vec(n_data)
  end if

  return
end
subroutine pochhammer_values ( n_data, x, y, fxy )

!*****************************************************************************80
!
!! POCHHAMMER_VALUES returns some values of the Pochhammer function.
!
!  Discussion:
!
!    Pochhammer(X,Y) = Gamma(X+Y) / Gamma(X)
!
!    For integer arguments, Pochhammer(M,N) = ( M + N - 1 )! / ( N - 1 )!
!
!    In Mathematica, the function can be evaluated by:
!
!      Pochhammer[X,Y]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, Y, the arguments of the function.
!
!    Output, real ( kind = 8 ) FXY, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 19

  real ( kind = 8 ), save, dimension ( n_max ) :: f_vec = (/ &
     720.0000000000000D+00, &
       1.875000000000000D+00, &
       1.000000000000000D+00, &
       4.500000000000000D+00, &
      24.75000000000000D+00, &
     110.0000000000000D+00, &
     377.1036305819165D+00, &
       4.362197352456253D+00, &
       1.000000000000000D+00, &
       1.467150493866654D+00, &
       2.180949074356397D+00, &
       3.282686710888467D+00, &
       5.000000000000000D+00, &
       7.702540092799931D+00, &
      11.99521990896018D+00, &
      18.87544858760869D+00, &
      30.00000000000000D+00, &
      48.14087557999957D+00, &
      77.96892940824118D+00 /)
  real ( kind = 8 ) fxy
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     1.00D+00, &
     0.50D+00, &
     4.50D+00, &
     4.50D+00, &
     4.50D+00, &
    10.00D+00, &
    10.00D+00, &
     7.25D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( n_max ) :: y_vec = (/ &
    6.00D+00, &
    3.00D+00, &
    0.00D+00, &
    1.00D+00, &
    2.00D+00, &
    2.00D+00, &
    2.50D+00, &
    0.75D+00, &
    0.00D+00, &
    0.25D+00, &
    0.50D+00, &
    0.75D+00, &
    1.00D+00, &
    1.25D+00, &
    1.50D+00, &
    1.75D+00, &
    2.00D+00, &
    2.25D+00, &
    2.50D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    y = 0.0D+00
    fxy = 0.0D+00
  else
    x = x_vec(n_data)
    y = y_vec(n_data)
    fxy = f_vec(n_data)
  end if

  return
end
subroutine poisson_cdf_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! POISSON_CDF_VALUES returns some values of the Poisson CDF.
!
!  Discussion:
!
!    CDF(X)(A) is the probability of at most X successes in unit time,
!    given that the expected mean number of successes is A.
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = PoissonDistribution [ a ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, integer ( kind = 4 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.02D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    2.00D+00, &
    2.00D+00, &
    2.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.9801986733067553D+00, &
    0.9048374180359596D+00, &
    0.9953211598395555D+00, &
    0.6065306597126334D+00, &
    0.9097959895689501D+00, &
    0.9856123220330293D+00, &
    0.3678794411714423D+00, &
    0.7357588823428846D+00, &
    0.9196986029286058D+00, &
    0.9810118431238462D+00, &
    0.1353352832366127D+00, &
    0.4060058497098381D+00, &
    0.6766764161830635D+00, &
    0.8571234604985470D+00, &
    0.6737946999085467D-02, &
    0.4042768199451280D-01, &
    0.1246520194830811D+00, &
    0.2650259152973617D+00, &
    0.4404932850652124D+00, &
    0.6159606548330631D+00, &
    0.7621834629729387D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) x
  integer ( kind = 4 ), save, dimension ( n_max ) :: x_vec = (/ &
     0, 0, 1, 0, &
     1, 2, 0, 1, &
     2, 3, 0, 1, &
     2, 3, 0, 1, &
     2, 3, 4, 5, &
     6 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine polylogarithm_values ( n_data, n, z, fx )

!*****************************************************************************80
!
!! POLYLOGARITHM_VALUES returns some values of the polylogarithm.
!
!  Discussion:
!
!    The polylogarithm of n and z is defined as
!
!      f[n,z] = Sum ( 1 <= k < +oo ) z^k / k^n
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyLog[n,z]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the exponent of the denominator.
!
!    Output, real ( kind = 8 ) Z, the base of the numerator.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1644934066848226D+01, &
    0.1202056903159594D+01, &
    0.1000994575127818D+01, &
    0.5822405264650125D+00, &
    0.5372131936080402D+00, &
    0.5002463206060068D+00, &
    0.3662132299770635D+00, &
    0.3488278611548401D+00, &
    0.3334424797228716D+00, &
    0.1026177910993911D+00, &
    0.1012886844792230D+00, &
    0.1000097826564961D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     2, 3, 10, 2, 3, 10, 2, 3, 10, 2, 3, 10 /)
  real ( kind = 8 ) z
  real ( kind = 8 ), save, dimension ( n_max ) :: z_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.3333333333333333D+00, &
    0.3333333333333333D+00, &
    0.3333333333333333D+00, &
    0.1000000000000000D+00, &
    0.1000000000000000D+00, &
    0.1000000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    z = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    z = z_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine prandtl_values ( n_data, tc, p, pr )

!*****************************************************************************80
!
!! PRANDTL_VALUES returns some values of the Prandtl number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher, George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) PR, the Prandtl number, dimensionless.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 35

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) pr
  real ( kind = 8 ), save, dimension ( n_max ) :: pr_vec = (/ &
    13.50D+00, &
    13.48D+00, &
    13.46D+00, &
    13.39D+00, &
    13.27D+00, &
    13.15D+00, &
    13.04D+00, &
    12.93D+00, &
    12.83D+00, &
    12.73D+00, &
    12.63D+00, &
    12.53D+00, &
    12.43D+00, &
    12.34D+00, &
    12.25D+00, &
    12.08D+00, &
    11.92D+00, &
    11.77D+00, &
    11.62D+00, &
    11.48D+00, &
    11.36D+00, &
    11.23D+00, &
    11.12D+00, &
    10.91D+00, &
    10.72D+00, &
    10.55D+00, &
     6.137D+00, &
     3.555D+00, &
     2.378D+00, &
     1.000D+00, &
     0.974D+00, &
     0.960D+00, &
     0.924D+00, &
     0.899D+00, &
     0.882D+00 /)
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
       1.0D+00, &
       5.0D+00, &
      10.0D+00, &
      25.0D+00, &
      50.0D+00, &
      75.0D+00, &
     100.0D+00, &
     125.0D+00, &
     150.0D+00, &
     175.0D+00, &
     200.0D+00, &
     225.0D+00, &
     250.0D+00, &
     275.0D+00, &
     300.0D+00, &
     350.0D+00, &
     400.0D+00, &
     450.0D+00, &
     500.0D+00, &
     550.0D+00, &
     600.0D+00, &
     650.0D+00, &
     700.0D+00, &
     800.0D+00, &
     900.0D+00, &
    1000.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
     25.0D+00, &
     50.0D+00, &
     75.0D+00, &
    100.0D+00, &
    150.0D+00, &
    200.0D+00, &
    400.0D+00, &
    600.0D+00, &
    800.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    tc = 0.0D+00
    p = 0.0D+00
    pr = 0.0D+00
  else
    tc = tc_vec(n_data)
    p = p_vec(n_data)
    pr = pr_vec(n_data)
  end if

  return
end
subroutine prime_values ( n_data, n, p )

!*****************************************************************************80
!
!! PRIME_VALUES returns values of the prime function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Prime[n]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the index of the prime.
!
!    Output, integer ( kind = 4 ) P, the value of the prime.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 24

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
          1, &
          2, &
          4, &
          8, &
         16, &
         32, &
         64, &
        128, &
        256, &
        512, &
       1000, &
       2000, &
       4000, &
       8000, &
      16000, &
      32000, &
      64000, &
     128000, &
     256000, &
     512000, &
    1024000, &
    2048000, &
    4096000, &
    8129000 /)
  integer ( kind = 4 ) p
  integer ( kind = 4 ), save, dimension ( n_max ) :: p_vec = (/ &
            2, &
            3, &
            7, &
           19, &
           53, &
          131, &
          311, &
          719, &
         1619, &
         3671, &
         7919, &
        17389, &
        37813, &
        81799, &
       176081, &
       376127, &
       800573, &
      1698077, &
      3588941, &
      7559173, &
     15881419, &
     33283031, &
     69600977, &
    145253029 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    p = 0
  else
    n = n_vec(n_data)
    p = p_vec(n_data)
  end if

  return
end
subroutine psat_values ( n_data, tc, p )

!*****************************************************************************80
!
!! PSAT_VALUES returns some values of the saturation pressure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher, George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the saturation pressure, in bar.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
    0.0061173D+00, &
    0.0065716D+00, &
    0.0087260D+00, &
    0.12344D+00, &
    1.0132D+00, &
    2.3201D+00, &
    4.7572D+00, &
    15.537D+00, &
    39.737D+00, &
    85.838D+00, &
    165.21D+00, &
    220.55D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
    0.100000D-01, &
    0.100000D+01, &
    0.500000D+01, &
    0.500000D+02, &
    0.100000D+03, &
    0.125000D+03, &
    0.150000D+03, &
    0.200000D+03, &
    0.250000D+03, &
    0.300000D+03, &
    0.350000D+03, &
    0.373976D+03 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    tc = 0.0D+00
    p = 0.0D+00
  else
    tc = tc_vec(n_data)
    p = p_vec(n_data)
  end if

  return
end
subroutine psi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! PSI_VALUES returns some values of the Psi or Digamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyGamma[x]
!
!    or
!
!      PolyGamma[0,x]
!
!    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
!
!    PSI(1) = -Euler's constant.
!
!    PSI(X+1) = PSI(X) + 1 / X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.5772156649015329D+00, &
    -0.4237549404110768D+00, &
    -0.2890398965921883D+00, &
    -0.1691908888667997D+00, &
    -0.6138454458511615D-01, &
     0.3648997397857652D-01, &
     0.1260474527734763D+00, &
     0.2085478748734940D+00, &
     0.2849914332938615D+00, &
     0.3561841611640597D+00, &
     0.4227843350984671D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine r8_factorial_values ( n_data, n, fn )

!*****************************************************************************80
!
!! R8_FACTORIAL_VALUES returns values of the real factorial function.
!
!  Discussion:
!
!    0! = 1
!    I! = Product ( 1 <= J <= I ) J
!
!    Although the factorial is an integer valued function, it quickly
!    becomes too large for an integer to hold.  This routine still accepts
!    an integer as the input argument, but returns the function value
!    as a real number.
!
!    In Mathematica, the function can be evaluated by:
!
!      n!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  real ( kind = 8 ) fn
  real ( kind = 8 ), save, dimension ( n_max ) :: fn_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.6000000000000000D+01, &
    0.2400000000000000D+02, &
    0.1200000000000000D+03, &
    0.7200000000000000D+03, &
    0.5040000000000000D+04, &
    0.4032000000000000D+05, &
    0.3628800000000000D+06, &
    0.3628800000000000D+07, &
    0.3991680000000000D+08, &
    0.4790016000000000D+09, &
    0.6227020800000000D+10, &
    0.8717829120000000D+11, &
    0.1307674368000000D+13, &
    0.2092278988800000D+14, &
    0.3556874280960000D+15, &
    0.6402373705728000D+16, &
    0.1216451004088320D+18, &
    0.2432902008176640D+19, &
    0.1551121004333099D+26, &
    0.3041409320171338D+65, &
    0.9332621544394415D+158, &
    0.5713383956445855D+263 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
       0, &
       1, &
       2, &
       3, &
       4, &
       5, &
       6, &
       7, &
       8, &
       9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      25, &
      50, &
     100, &
     150 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    fn = 0.0D+00
  else
    n = n_vec(n_data)
    fn = fn_vec(n_data)
  end if

  return
end
subroutine r8_factorial_log_values ( n_data, n, fn )

!*****************************************************************************80
!
!! R8_FACTORIAL_LOG_VALUES returns values of log(n!).
!
!  Discussion:
!
!    The function log(n!) can be written as
!
!     log(n!) = sum ( 1 <= i <= n ) log ( i )
!
!    In Mathematica, the function can be evaluated by:
!
!      Log[n!]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 27

  real ( kind = 8 ) fn
  real ( kind = 8 ), save, dimension ( n_max ) :: fn_vec = (/ &
    0.0000000000000000D+00, &
    0.0000000000000000D+00, &
    0.6931471805599453D+00, &
    0.1791759469228055D+01, &
    0.3178053830347946D+01, &
    0.4787491742782046D+01, &
    0.6579251212010101D+01, &
    0.8525161361065414D+01, &
    0.1060460290274525D+02, &
    0.1280182748008147D+02, &
    0.1510441257307552D+02, &
    0.1750230784587389D+02, &
    0.1998721449566189D+02, &
    0.2255216385312342D+02, &
    0.2519122118273868D+02, &
    0.2789927138384089D+02, &
    0.3067186010608067D+02, &
    0.3350507345013689D+02, &
    0.3639544520803305D+02, &
    0.3933988418719949D+02, &
    0.4233561646075349D+02, &
    0.5800360522298052D+02, &
    0.1484777669517730D+03, &
    0.3637393755555635D+03, &
    0.6050201058494237D+03, &
    0.2611330458460156D+04, &
    0.5912128178488163D+04 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
       0, &
       1, &
       2, &
       3, &
       4, &
       5, &
       6, &
       7, &
       8, &
       9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      25, &
      50, &
     100, &
     150, &
     500, &
    1000 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    fn = 0.0D+00
  else
    n = n_vec(n_data)
    fn = fn_vec(n_data)
  end if

  return
end
subroutine rayleigh_cdf_values ( n_data, sigma, x, fx )

!*****************************************************************************80
!
!! RAYLEIGH_CDF_VALUES returns some values of the Rayleigh CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = RayleighDistribution [ sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) SIGMA, the shape parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 9

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8646647167633873D+00, &
    0.9996645373720975D+00, &
    0.9999999847700203D+00, &
    0.999999999999987D+00, &
    0.8646647167633873D+00, &
    0.3934693402873666D+00, &
    0.1992625970831920D+00, &
    0.1175030974154046D+00, &
    0.7688365361336422D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine secvir_values ( n_data, tc, vir )

!*****************************************************************************80
!
!! SECVIR_VALUES returns some values of the second virial coefficient.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher, George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) VIR, the second virial coefficient, in
!    m^3/kg.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 19

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
       0.0D+00, &
       5.0D+00, &
      10.0D+00, &
      20.0D+00, &
      30.0D+00, &
      40.0D+00, &
      60.0D+00, &
      90.0D+00, &
     120.0D+00, &
     150.0D+00, &
     180.0D+00, &
     210.0D+00, &
     240.0D+00, &
     300.0D+00, &
     400.0D+00, &
     500.0D+00, &
     700.0D+00, &
    1000.0D+00, &
    2000.0D+00 /)
  real ( kind = 8 ) vir
  real ( kind = 8 ), save, dimension ( n_max ) :: vir_vec = (/ &
    -98.96D+00, &
    -90.08D+00, &
    -82.29D+00, &
    -69.36D+00, &
    -59.19D+00, &
    -51.07D+00, &
    -39.13D+00, &
    -27.81D+00, &
    -20.83D+00, &
    -16.21D+00, &
    -12.98D+00, &
    -10.63D+00, &
     -8.85D+00, &
     -6.39D+00, &
     -4.03D+00, &
     -2.71D+00, &
     -1.32D+00, &
     -0.39D+00, &
      0.53D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    tc = 0.0D+00
    vir = 0.0D+00
  else
    tc = tc_vec(n_data)
    vir = vir_vec(n_data)
  end if

  return
end
subroutine shi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! SHI_VALUES returns some values of the hyperbolic sine integral function.
!
!  Discussion:
!
!    SHI(X) = integral ( 0 <= T <= X ) sinh ( T ) / T dt
!
!    In Mathematica, the function can be evaluated by:
!
!      SinhIntegral[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5069967498196672D+00, &
    0.6121303965633808D+00, &
    0.7193380189288998D+00, &
    0.8289965633789345D+00, &
    0.9414978265114335D+00, &
    1.057250875375729D+00, &
    1.300250361022057D+00, &
    1.561713388361002D+00, &
    1.845814141358504D+00, &
    2.157290343425901D+00, &
    2.501567433354976D+00, &
    3.549340406224435D+00, &
    4.973440475859807D+00, &
    6.966162067504942D+00, &
    9.817326911233034D+00, &
    13.96788504934715D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     1.0D+00, &
     1.2D+00, &
     1.4D+00, &
     1.6D+00, &
     1.8D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00, &
     4.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine si_values ( n_data, x, fx )

!*****************************************************************************80
!
!! SI_VALUES returns some values of the sine integral function.
!
!  Discussion:
!
!    SI(X) = integral ( 0 <= T <= X ) sin ( T ) / T dt
!
!    In Mathematica, the function can be evaluated by:
!
!      SinIntegral[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.4931074180430667D+00, &
    0.5881288096080801D+00, &
    0.6812222391166113D+00, &
    0.7720957854819966D+00, &
    0.8604707107452929D+00, &
    0.9460830703671830D+00, &
    0.1108047199013719D+01, &
    0.1256226732779218D+01, &
    0.1389180485870438D+01, &
    0.1505816780255579D+01, &
    0.1605412976802695D+01, &
    0.1778520173443827D+01, &
    0.1848652527999468D+01, &
    0.1833125398665997D+01, &
    0.1758203138949053D+01, &
    0.1654140414379244D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.5D+00, &
     0.6D+00, &
     0.7D+00, &
     0.8D+00, &
     0.9D+00, &
     1.0D+00, &
     1.2D+00, &
     1.4D+00, &
     1.6D+00, &
     1.8D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00, &
     4.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine sigma_values ( n_data, n, c )

!*****************************************************************************80
!
!! SIGMA_VALUES returns some values of the Sigma function.
!
!  Discussion:
!
!    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
!
!    In Mathematica, the function can be evaluated by:
!
!      DivisorSigma[1,n]
!
!  First values:
!
!     N  SIGMA(N)
!
!     1    1
!     2    3
!     3    4
!     4    7
!     5    6
!     6   12
!     7    8
!     8   15
!     9   13
!    10   18
!    11   12
!    12   28
!    13   14
!    14   24
!    15   24
!    16   31
!    17   18
!    18   39
!    19   20
!    20   42
!
!  Formula:
!
!    SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
!
!    SIGMA(P^K) = ( P^(K+1) - 1 ) / ( P - 1 ) if P is prime.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the Sigma function.
!
!    Output, integer ( kind = 4 ) C, the value of the Sigma function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
     1,    3,    4,    7,    6,   12,    8,   15,   13,   18, &
    72,  128,  255,  176,  576, 1170,  618,  984, 2232, 2340 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,   10, &
     30, 127, 128, 129, 210, 360, 617, 815, 816, 1000 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine sin_values ( n_data, x, fx )

!*****************************************************************************80
!
!! SIN_VALUES returns some values of the sine function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sin[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.00000000000000000000D+00, &
     0.25881904510252076235D+00, &
     0.47942553860420300027D+00, &
     0.50000000000000000000D+00, &
     0.70710678118654752440D+00, &
     0.84147098480789650665D+00, &
     0.86602540378443864676D+00, &
     1.00000000000000000000D+00, &
     0.90929742682568169540D+00, &
     0.14112000805986722210D+00, &
     0.00000000000000000000D+00, &
    -0.75680249530792825137D+00, &
    -0.95892427466313846889D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000000D+00, &
    0.26179938779914943654D+00, &
    0.50000000000000000000D+00, &
    0.52359877559829887308D+00, &
    0.78539816339744830962D+00, &
    1.0000000000000000000D+00, &
    1.0471975511965977462D+00, &
    1.5707963267948966192D+00, &
    2.0000000000000000000D+00, &
    3.0000000000000000000D+00, &
    3.1415926535897932385D+00, &
    4.0000000000000000000D+00, &
    5.0000000000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine sin_degree_values ( n_data, x, fx )

!*********************************************************************72
!
!! SIN_DEGREE_VALUES: the sine function with argument in degrees.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sin[x Degree]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.087155742747658173558D+00, &
     0.000000000000000000000D+00, &
     0.017452406437283512819D+00, &
     0.034899496702500971646D+00, &
     0.052335956242943832722D+00, &
     0.069756473744125300776D+00, &
     0.087155742747658173558D+00, &
     0.17364817766693034885D+00, &
     0.25881904510252076235D+00, &
     0.50000000000000000000D+00, &
     0.70710678118654752440D+00, &
     0.86602540378443864676D+00, &
     0.96592582628906828675D+00, &
     0.99619469809174553230D+00, &
     0.99756405025982424761D+00, &
     0.99862953475457387378D+00, &
     0.99939082701909573001D+00, &
     0.99984769515639123916D+00, &
     1.0000000000000000000D+00, &
     0.99984769515639123916D+00, &
     0.96592582628906828675D+00, &
     0.00000000000000000000D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     -5.0D+00, &
      0.0D+00, &
      1.0D+00, &
      2.0D+00, &
      3.0D+00, &
      4.0D+00, &
      5.0D+00, &
     10.0D+00, &
     15.0D+00, &
     30.0D+00, &
     45.0D+00, &
     60.0D+00, &
     75.0D+00, &
     85.0D+00, &
     86.0D+00, &
     87.0D+00, &
     88.0D+00, &
     89.0D+00, &
     90.0D+00, &
     91.0D+00, &
    105.0D+00, &
    180.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine sin_power_int_values ( n_data, a, b, n, fx )

!*****************************************************************************80
!
!! SIN_POWER_INT_VALUES returns some values of the sine power integral.
!
!  Discussion:
!
!    The function has the form
!
!      SIN_POWER_INT(A,B,N) = integral ( A <= T <= B ) ( sin(T) )^N dt
!
!    In Mathematica, the function can be evaluated by:
!
!      Integrate [ ( Sin[x] )^n, { x, a, b } ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, integer ( kind = 4 ) N, the power.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.10D+02, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.10D+01, &
     0.00D+00, &
     0.00D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.20D+02, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.20D+01, &
     0.20D+01, &
     0.10D+01, &
     0.10D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.10000000000000000000D+02, &
    0.45969769413186028260D+00, &
    0.27267564329357957615D+00, &
    0.17894056254885809051D+00, &
    0.12402556531520681830D+00, &
    0.88974396451575946519D-01, &
    0.90393123848149944133D+00, &
    0.81495684202992349481D+00, &
    0.21887522421729849008D-01, &
    0.17023439374069324596D-01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, &
     1, &
     2, &
     3, &
     4, &
     5, &
     5, &
     5, &
    10, &
    11 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    n = 0
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine sinh_values ( n_data, x, fx )

!*****************************************************************************80
!
!! SINH_VALUES returns some values of the hyperbolic sine function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sinh[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
      -74.203210577788758977D+00, &
       -1.1752011936438014569D+00, &
        0.00000000000000000000D+00, &
        0.10016675001984402582D+00, &
        0.20133600254109398763D+00, &
        0.30452029344714261896D+00, &
        0.41075232580281550854D+00, &
        0.52109530549374736162D+00, &
        0.63665358214824127112D+00, &
        0.75858370183953350346D+00, &
        0.88810598218762300657D+00, &
        1.0265167257081752760D+00, &
        1.1752011936438014569D+00, &
        3.6268604078470187677D+00, &
       10.017874927409901899D+00, &
       27.289917197127752449D+00, &
       74.203210577788758977D+00, &
    11013.232874703393377D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
   -5.0D+00, &
   -1.0D+00, &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
   10.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine six_j_values ( n_data, j1, j2, j3, j4, j5, j6, fx )

!*****************************************************************************80
!
!! SIX_J_VALUES returns some values of the Wigner 6J function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      SixJSymbol[{j1,j2,j3},{j4,j5,j6}]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) J1, J2, J3, J4, J5, J6, the arguments
!    of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.03490905138373300D+00,  &
    -0.03743025039659792D+00,  &
     0.01890866390959560D+00,  &
     0.007342448254928643D+00, &
    -0.02358935185081794D+00,  &
     0.01913476955215437D+00,  &
     0.001288017397724172D+00, &
    -0.01930018366290527D+00,  &
     0.01677305949382889D+00,  &
     0.005501147274850949D+00, &
    -0.02135439790896831D+00,  &
     0.003460364451435387D+00, &
     0.02520950054795585D+00,  &
     0.01483990561221713D+00,  &
     0.002708577680633186D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) j1
  real ( kind = 8 ), save, dimension ( n_max ) :: j1_vec = (/ &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    6.0D+00, &
    7.0D+00, &
    8.0D+00, &
    9.0D+00, &
   10.0D+00, &
   11.0D+00, &
   12.0D+00, &
   13.0D+00, &
   14.0D+00, &
   15.0D+00 /)
  real ( kind = 8 ) j2
  real ( kind = 8 ), save, dimension ( n_max ) :: j2_vec = (/ &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00, &
    8.0D+00 /)
  real ( kind = 8 ) j3
  real ( kind = 8 ), save, dimension ( n_max ) :: j3_vec = (/ &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00 /)
  real ( kind = 8 ) j4
  real ( kind = 8 ), save, dimension ( n_max ) :: j4_vec = (/ &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00, &
    6.5D+00 /)
  real ( kind = 8 ) j5
  real ( kind = 8 ), save, dimension ( n_max ) :: j5_vec = (/ &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00 /)
  real ( kind = 8 ) j6
  real ( kind = 8 ), save, dimension ( n_max ) :: j6_vec = (/ &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00, &
    7.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    j1 = 0.0D+00
    j2 = 0.0D+00
    j3 = 0.0D+00
    j4 = 0.0D+00
    j5 = 0.0D+00
    j6 = 0.0D+00
    fx = 0.0D+00
  else
    j1 = j1_vec(n_data)
    j2 = j2_vec(n_data)
    j3 = j3_vec(n_data)
    j4 = j4_vec(n_data)
    j5 = j5_vec(n_data)
    j6 = j6_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine sound_values ( n_data, tc, p, c )

!*****************************************************************************80
!
!! SOUND_VALUES returns some values of the speed of sound.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher, George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) C, the speed of sound, in m/s.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) c
  real ( kind = 8 ), save, dimension ( n_max ) :: c_vec = (/ &
    1401.0D+00, &
     472.8D+00, &
     533.7D+00, &
     585.7D+00, &
     609.5D+00, &
     632.2D+00, &
     674.6D+00, &
     713.9D+00, &
     802.0D+00, &
     880.1D+00, &
    1017.8D+00, &
    1115.9D+00, &
    1401.7D+00, &
    1402.6D+00, &
    1409.6D+00, &
    1418.1D+00, &
    1443.1D+00, &
    1484.6D+00, &
    1577.1D+00, &
    1913.4D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       5.0D+00, &
      10.0D+00, &
      50.0D+00, &
     100.0D+00, &
     250.0D+00, &
     500.0D+00, &
    1000.0D+00, &
    2500.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
       0.0D+00, &
     100.0D+00, &
     200.0D+00, &
     300.0D+00, &
     350.0D+00, &
     400.0D+00, &
     500.0D+00, &
     600.0D+00, &
     850.0D+00, &
    1100.0D+00, &
    1600.0D+00, &
    2000.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    tc = 0.0D+00
    p = 0.0D+00
    c = 0.0D+00
  else
    tc = tc_vec(n_data)
    p = p_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine sphere_unit_area_values ( n_data, n, area )

!*****************************************************************************80
!
!! SPHERE_UNIT_AREA_VALUES returns some areas of the unit sphere in ND.
!
!  Discussion:
!
!    The formula for the surface area of the unit sphere in N dimensions is:
!
!      Sphere_Unit_Area ( N ) = 2 * pi^(N/2) / Gamma ( N / 2 )
!
!    Some values of the function include:
!
!       N   Area
!
!       2    2        * PI
!       3  ( 4 /    ) * PI
!       4  ( 2 /   1) * PI^2
!       5  ( 8 /   3) * PI^2
!       6  ( 1 /   1) * PI^3
!       7  (16 /  15) * PI^3
!       8  ( 1 /   3) * PI^4
!       9  (32 / 105) * PI^4
!      10  ( 1 /  12) * PI^5
!
!    For the unit sphere, Area(N) = N * Volume(N)
!
!    In Mathematica, the function can be evaluated by:
!
!      2 * Pi^(n/2) / Gamma[n/2]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, real ( kind = 8 ) AREA, the area of the unit sphere
!    in that dimension.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) area
  real ( kind = 8 ), save, dimension ( n_max ) :: area_vec = (/ &
    0.2000000000000000D+01, &
    0.6283185307179586D+01, &
    0.1256637061435917D+02, &
    0.1973920880217872D+02, &
    0.2631894506957162D+02, &
    0.3100627668029982D+02, &
    0.3307336179231981D+02, &
    0.3246969701133415D+02, &
    0.2968658012464836D+02, &
    0.2550164039877345D+02, &
    0.2072514267328890D+02, &
    0.1602315322625507D+02, &
    0.1183817381218268D+02, &
    0.8389703410491089D+01, &
    0.5721649212349567D+01, &
    0.3765290085742291D+01, &
    0.2396678817591364D+01, &
    0.1478625959000308D+01, &
    0.8858104195716824D+00, &
    0.5161378278002812D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     1, &
     2, &
     3, &
     4, &
     5, &
     6, &
     7, &
     8, &
     9, &
    10, &
    11, &
    12, &
    13, &
    14, &
    15, &
    16, &
    17, &
    18, &
    19, &
    20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    area = 0.0D+00
  else
    n = n_vec(n_data)
    area = area_vec(n_data)
  end if

  return
end
subroutine sphere_unit_volume_values ( n_data, n, volume )

!*****************************************************************************80
!
!! SPHERE_UNIT_VOLUME_VALUES returns some volumes of the unit sphere in ND.
!
!  Discussion:
!
!    The formula for the volume of the unit sphere in N dimensions is
!
!      Volume(N) = 2 * pi^(N/2) / ( N * Gamma ( N / 2 ) )
!
!    This function satisfies the relationships:
!
!      Volume(N) = 2 * pi * Volume(N-2) / N
!      Volume(N) = Area(N) / N
!
!    Some values of the function include:
!
!       N  Volume
!
!       1    1
!       2    1        * PI
!       3  ( 4 /   3) * PI
!       4  ( 1 /   2) * PI^2
!       5  ( 8 /  15) * PI^2
!       6  ( 1 /   6) * PI^3
!       7  (16 / 105) * PI^3
!       8  ( 1 /  24) * PI^4
!       9  (32 / 945) * PI^4
!      10  ( 1 / 120) * PI^5
!
!    In Mathematica, the function can be evaluated by:
!
!      2 * Pi^(n/2) / ( n * Gamma[n/2] )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the unit
!    sphere in that dimension.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     1,  2, &
     3,  4, &
     5,  6, &
     7,  8, &
     9, 10, &
    11, 12, &
    13, 14, &
    15, 16, &
    17, 18, &
    19, 20 /)
  real ( kind = 8 ) volume
  real ( kind = 8 ), save, dimension ( n_max ) :: volume_vec = (/ &
    0.2000000000000000D+01, &
    0.3141592653589793D+01, &
    0.4188790204786391D+01, &
    0.4934802200544679D+01, &
    0.5263789013914325D+01, &
    0.5167712780049970D+01, &
    0.4724765970331401D+01, &
    0.4058712126416768D+01, &
    0.3298508902738707D+01, &
    0.2550164039877345D+01, &
    0.1884103879389900D+01, &
    0.1335262768854589D+01, &
    0.9106287547832831D+00, &
    0.5992645293207921D+00, &
    0.3814432808233045D+00, &
    0.2353306303588932D+00, &
    0.1409811069171390D+00, &
    0.8214588661112823D-01, &
    0.4662160103008855D-01, &
    0.2580689139001406D-01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    volume = 0.0D+00
  else
    n = n_vec(n_data)
    volume = volume_vec(n_data)
  end if

  return
end
subroutine spherical_harmonic_values ( n_data, l, m, theta, phi, yr, yi )

!*****************************************************************************80
!
!! SPHERICAL_HARMONIC_VALUES returns values of spherical harmonic functions.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by
!
!      SphericalHarmonicY [ l, m, theta, phi ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) L, integer M, real ( kind = 8 ) THETA, PHI,
!    the arguments of the function.
!
!    Output, real ( kind = 8 ) YR, YI, the real and imaginary parts of
!    the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) l
  integer ( kind = 4 ), save, dimension ( n_max ) :: l_vec = (/ &
     0,  1,  2,  &
     3,  4,  5,  &
     5,  5,  5,  &
     5,  4,  4,  &
     4,  4,  4,  &
     3,  3,  3,  &
     3,  3 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
     0,  0,  1,  &
     2,  3,  5,  &
     4,  3,  2,  &
     1,  2,  2,  &
     2,  2,  2,  &
    -1, -1, -1,  &
    -1, -1 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) phi
  real ( kind = 8 ), save, dimension ( n_max ) :: phi_vec = (/ &
    0.1047197551196598D+01, 0.1047197551196598D+01, 0.1047197551196598D+01, &
    0.1047197551196598D+01, 0.1047197551196598D+01, 0.6283185307179586D+00, &
    0.6283185307179586D+00, 0.6283185307179586D+00, 0.6283185307179586D+00, &
    0.6283185307179586D+00, 0.7853981633974483D+00, 0.7853981633974483D+00, &
    0.7853981633974483D+00, 0.7853981633974483D+00, 0.7853981633974483D+00, &
    0.4487989505128276D+00, 0.8975979010256552D+00, 0.1346396851538483D+01, &
    0.1795195802051310D+01, 0.2243994752564138D+01 /)
  real ( kind = 8 ) theta
  real ( kind = 8 ), save, dimension ( n_max ) :: theta_vec = (/ &
    0.5235987755982989D+00, 0.5235987755982989D+00, 0.5235987755982989D+00, &
    0.5235987755982989D+00, 0.5235987755982989D+00, 0.2617993877991494D+00, &
    0.2617993877991494D+00, 0.2617993877991494D+00, 0.2617993877991494D+00, &
    0.2617993877991494D+00, 0.6283185307179586D+00, 0.1884955592153876D+01, &
    0.3141592653589793D+01, 0.4398229715025711D+01, 0.5654866776461628D+01, &
    0.3926990816987242D+00, 0.3926990816987242D+00, 0.3926990816987242D+00, &
    0.3926990816987242D+00, 0.3926990816987242D+00 /)
  real ( kind = 8 ) yi
  real ( kind = 8 ), save, dimension ( n_max ) :: yi_vec = (/ &
    0.0000000000000000D+00,  0.0000000000000000D+00, -0.2897056515173922D+00, &
    0.1916222768312404D+00,  0.0000000000000000D+00,  0.0000000000000000D+00, &
    0.3739289485283311D-02, -0.4219517552320796D-01,  0.1876264225575173D+00, &
   -0.3029973424491321D+00,  0.4139385503112256D+00, -0.1003229830187463D+00, &
    0.0000000000000000D+00, -0.1003229830187463D+00,  0.4139385503112256D+00, &
   -0.1753512375142586D+00, -0.3159720118970196D+00, -0.3940106541811563D+00, &
   -0.3940106541811563D+00, -0.3159720118970196D+00 /)
  real ( kind = 8 ) yr
  real ( kind = 8 ), save, dimension ( n_max ) :: yr_vec = (/ &
   0.2820947917738781D+00,  0.4231421876608172D+00, -0.1672616358893223D+00, &
  -0.1106331731112457D+00,  0.1354974113737760D+00,  0.5390423109043568D-03, &
  -0.5146690442951909D-02,  0.1371004361349490D-01,  0.6096352022265540D-01, &
  -0.4170400640977983D+00,  0.0000000000000000D+00,  0.0000000000000000D+00, &
   0.0000000000000000D+00,  0.0000000000000000D+00,  0.0000000000000000D+00, &
   0.3641205966137958D+00,  0.2519792711195075D+00,  0.8993036065704300D-01, &
  -0.8993036065704300D-01, -0.2519792711195075D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    l = 0
    m = 0
    theta = 0.0D+00
    phi = 0.0D+00
    yr = 0.0D+00
    yi = 0.0D+00
  else
    l = l_vec(n_data)
    m = m_vec(n_data)
    theta = theta_vec(n_data)
    phi = phi_vec(n_data)
    yr = yr_vec(n_data)
    yi = yi_vec(n_data)
  end if

  return
end
subroutine sqrt_values ( n_data, x, fx )

!*****************************************************************************80
!
!! SQRT_VALUES returns some values of the square root function.
!
!  Discussion:
!
!    SQRT(X) = positive real number Y such that Y * Y = X.
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 14

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.9000000040950000D-04, &
    0.3000000000000000D+00, &
    0.3162277660168379D+00, &
    0.6324555320336759D+00, &
    0.1000000000000000D+01, &
    0.1414213562373095D+01, &
    0.1732050807568877D+01, &
    0.1772453850905516D+01, &
    0.4358898943540674D+01, &
    0.5385164807134504D+01, &
    0.8426149773176359D+01, &
    0.9848857801796105D+01, &
    0.1111111106055556D+05 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
    0.8100000073710001D-08, &
    0.9000000000000000D-01, &
    0.1000000000000000D+00, &
    0.4000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3141592653589793D+01, &
    0.1900000000000000D+02, &
    0.2900000000000000D+02, &
    0.7100000000000000D+02, &
    0.9700000000000000D+02, &
    0.1234567890000000D+09 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine stirling1_values ( n_data, n, m, fx )

!*****************************************************************************80
!
!! STIRLING1_VALUES returns some values of the Stirling numbers, kind 1.
!
!  Discussion:
!
!    The absolute value of the Stirling number S1(N,M) gives the number
!    of permutations on N objects having exactly M cycles, while the
!    sign of the Stirling number records the sign (odd or even) of
!    the permutations.  For example, there are six permutations on 3 objects:
!
!      A B C   3 cycles (A) (B) (C)
!      A C B   2 cycles (A) (BC)
!      B A C   2 cycles (AB) (C)
!      B C A   1 cycle  (ABC)
!      C A B   1 cycle  (ABC)
!      C B A   2 cycles (AC) (B)
!
!    There are
!
!      2 permutations with 1 cycle, and S1(3,1) = 2
!      3 permutations with 2 cycles, and S1(3,2) = -3,
!      1 permutation with 3 cycles, and S1(3,3) = 1.
!
!    Since there are N! permutations of N objects, the sum of the absolute
!    values of the Stirling numbers in a given row,
!
!      sum ( 1 <= I <= N ) abs ( S1(N,I) ) = N!
!
!  First terms:
!
!    N/M:  1     2      3     4     5    6    7    8
!
!    1     1     0      0     0     0    0    0    0
!    2    -1     1      0     0     0    0    0    0
!    3     2    -3      1     0     0    0    0    0
!    4    -6    11     -6     1     0    0    0    0
!    5    24   -50     35   -10     1    0    0    0
!    6  -120   274   -225    85   -15    1    0    0
!    7   720 -1764   1624  -735   175  -21    1    0
!    8 -5040 13068 -13132  6769 -1960  322  -28    1
!
!    In Mathematica, the function can be evaluated by:
!
!      StirlingS1[n,m]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, M, the arguments of the function.
!
!    Output, integer ( kind = 4 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  integer ( kind = 4 ) fx
  integer ( kind = 4 ), save, dimension ( n_max ) :: fx_vec = (/ &
           0, &
           1, &
          -3, &
          11, &
         -50, &
         274, &
       -1764, &
       13068, &
     -109584, &
     1026576, &
      -13132, &
        6769, &
       -1960, &
         322, &
         -28, &
           1 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 8, 8, 8, 8, 8, 8 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    fx = 0
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine stirling2_values ( n_data, n, m, fx )

!*****************************************************************************80
!
!! STIRLING2_VALUES returns some values of the Stirling numbers, kind 2.
!
!  Discussion:
!
!    S2(N,M) represents the number of distinct partitions of N elements
!    into M nonempty sets.  For a fixed N, the sum of the Stirling
!    numbers S2(N,M) is represented by B(N), called "Bell's number",
!    and represents the number of distinct partitions of N elements.
!
!    For example, with 4 objects, there are:
!
!    1 partition into 1 set:
!
!      (A,B,C,D)
!
!    7 partitions into 2 sets:
!
!      (A,B,C) (D)
!      (A,B,D) (C)
!      (A,C,D) (B)
!      (A) (B,C,D)
!      (A,B) (C,D)
!      (A,C) (B,D)
!      (A,D) (B,C)
!
!    6 partitions into 3 sets:
!
!      (A,B) (C) (D)
!      (A) (B,C) (D)
!      (A) (B) (C,D)
!      (A,C) (B) (D)
!      (A,D) (B) (C)
!      (A) (B,D) (C)
!
!    1 partition into 4 sets:
!
!      (A) (B) (C) (D)
!
!    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
!
!
!  First terms:
!
!    N/M: 1    2    3    4    5    6    7    8
!
!    1    1    0    0    0    0    0    0    0
!    2    1    1    0    0    0    0    0    0
!    3    1    3    1    0    0    0    0    0
!    4    1    7    6    1    0    0    0    0
!    5    1   15   25   10    1    0    0    0
!    6    1   31   90   65   15    1    0    0
!    7    1   63  301  350  140   21    1    0
!    8    1  127  966 1701 1050  266   28    1
!
!    In Mathematica, the function can be evaluated by:
!
!      StirlingS2[n,m]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, M, the arguments of the function.
!
!    Output, integer ( kind = 4 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  integer ( kind = 4 ) fx
  integer ( kind = 4 ), save, dimension ( n_max ) :: fx_vec = (/ &
           0, &
           1, &
           3, &
           7, &
          15, &
          31, &
          63, &
         127, &
         255, &
         511, &
         966, &
        1701, &
        1050, &
         266, &
          28, &
           1 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 8, 8, 8, 8, 8, 8 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    fx = 0
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine stromgen_values ( n_data, x, fx )

!*****************************************************************************80
!
!! STROMGEN_VALUES returns some values of the Stromgen function.
!
!  Discussion:
!
!    The function is defined by:
!
!      STROMGEN(X) = integral ( 0 <= t <= X ) t^7 * exp(2*t) / (exp(t)-1)^3 dt
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.21901065985698662316D-15, &
    0.22481399438625244761D-12, &
    0.23245019579558857124D-09, &
    0.24719561475975007037D-06, &
    0.28992610989833245669D-03, &
    0.10698146390809715091D-01, &
    0.89707650964424730705D-01, &
    0.40049605719592888440D+00, &
    0.30504104398079096598D+01, &
    0.11367704858439426431D+02, &
    0.12960679405324786954D+02, &
    0.18548713944748505675D+02, &
    0.27866273821903121400D+02, &
    0.51963334071699323351D+02, &
    0.10861016747891228129D+03, &
    0.15378903316556621624D+03, &
    0.19302665532558721516D+03, &
    0.19636850166006541482D+03, &
    0.19651946766008214217D+03, &
    0.19651956920868316152D+03 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0078125000D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.1250000000D+00, &
      4.5000000000D+00, &
      5.0000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine struve_h0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! STRUVE_H0_VALUES returns some values of the Struve H0 function.
!
!  Discussion:
!
!    The function is defined by:
!
!      HO(x) = 2/pi * integral ( 0 <= t <= pi/2 ) sin ( x * cos ( t ) ) dt
!
!    In Mathematica, the function can be evaluated by:
!
!      StruveH[0,x]
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.12433974658847434366D-02, &
    -0.49735582423748415045D-02, &
     0.39771469054536941564D-01, &
    -0.15805246001653314198D+00, &
     0.56865662704828795099D+00, &
     0.66598399314899916605D+00, &
     0.79085884950809589255D+00, &
    -0.13501457342248639716D+00, &
     0.20086479668164503137D+00, &
    -0.11142097800261991552D+00, &
    -0.17026804865989885869D+00, &
    -0.13544931808186467594D+00, &
     0.94393698081323450897D-01, &
    -0.10182482016001510271D+00, &
     0.96098421554162110012D-01, &
    -0.85337674826118998952D-01, &
    -0.76882290637052720045D-01, &
     0.47663833591418256339D-01, &
    -0.70878751689647343204D-01, &
     0.65752908073352785368D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
       0.0019531250D+00, &
      -0.0078125000D+00, &
       0.0625000000D+00, &
      -0.2500000000D+00, &
       1.0000000000D+00, &
       1.2500000000D+00, &
       2.0000000000D+00, &
      -4.0000000000D+00, &
       7.5000000000D+00, &
      11.0000000000D+00, &
      11.5000000000D+00, &
     -16.0000000000D+00, &
      20.0000000000D+00, &
      25.0000000000D+00, &
     -30.0000000000D+00, &
      50.0000000000D+00, &
      75.0000000000D+00, &
     -80.0000000000D+00, &
     100.0000000000D+00, &
    -125.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine struve_h1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! STRUVE_H1_VALUES returns some values of the Struve H1 function.
!
!  Discussion:
!
!    The function is defined by:
!
!      H1(x) = 2*x/pi * integral ( 0 <= t <= pi/2 )
!        sin ( x * cos ( t ) )^2 * sin ( t ) dt
!
!    In Mathematica, the function can be evaluated by:
!
!      StruveH[1,x]
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.80950369576367526071D-06, &
    0.12952009724113229165D-04, &
    0.82871615165407083021D-03, &
    0.13207748375849572564D-01, &
    0.19845733620194439894D+00, &
    0.29853823231804706294D+00, &
    0.64676372828356211712D+00, &
    0.10697266613089193593D+01, &
    0.38831308000420560970D+00, &
    0.74854243745107710333D+00, &
    0.84664854642567359993D+00, &
    0.58385732464244384564D+00, &
    0.80600584524215772824D+00, &
    0.53880362132692947616D+00, &
    0.72175037834698998506D+00, &
    0.58007844794544189900D+00, &
    0.60151910385440804463D+00, &
    0.70611511147286827018D+00, &
    0.61631110327201338454D+00, &
    0.62778480765443656489D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
       0.0019531250D+00, &
      -0.0078125000D+00, &
       0.0625000000D+00, &
      -0.2500000000D+00, &
       1.0000000000D+00, &
       1.2500000000D+00, &
       2.0000000000D+00, &
      -4.0000000000D+00, &
       7.5000000000D+00, &
      11.0000000000D+00, &
      11.5000000000D+00, &
     -16.0000000000D+00, &
      20.0000000000D+00, &
      25.0000000000D+00, &
     -30.0000000000D+00, &
      50.0000000000D+00, &
      75.0000000000D+00, &
     -80.0000000000D+00, &
     100.0000000000D+00, &
    -125.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine struve_l0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! STRUVE_L0_VALUES returns some values of the Struve L0 function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      StruveL[0,x]
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.12433985199262820188D-02, &
    -0.19896526647882937004D-01, &
     0.79715713253115014945D-01, &
    -0.32724069939418078025D+00, &
     0.71024318593789088874D+00, &
     0.19374337579914456612D+01, &
    -0.11131050203248583431D+02, &
     0.16850062034703267148D+03, &
    -0.28156522493745948555D+04, &
     0.89344618796978400815D+06, &
     0.11382025002851451057D+07, &
    -0.23549701855860190304D+07, &
     0.43558282527641046718D+08, &
     0.49993516476037957165D+09, &
    -0.57745606064408041689D+10, &
     0.78167229782395624524D+12, &
    -0.14894774793419899908D+17, &
     0.29325537838493363267D+21, &
     0.58940770556098011683D+25, &
    -0.12015889579125463605D+30 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
     -0.0312500000D+00, &
      0.1250000000D+00, &
     -0.5000000000D+00, &
      1.0000000000D+00, &
      2.0000000000D+00, &
     -4.0000000000D+00, &
      7.0000000000D+00, &
    -10.0000000000D+00, &
     16.0000000000D+00, &
     16.2500000000D+00, &
    -17.0000000000D+00, &
     20.0000000000D+00, &
     22.5000000000D+00, &
    -25.0000000000D+00, &
     30.0000000000D+00, &
    -40.0000000000D+00, &
     50.0000000000D+00, &
     60.0000000000D+00, &
    -70.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine struve_l1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! STRUVE_L1_VALUES returns some values of the Struve L1 function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      StruveL[1,x]
!
!    The data was reported by McLeod.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.80950410749865126939D-06, &
    0.20724649092571514607D-03, &
    0.33191834066894516744D-02, &
    0.53942182623522663292D-01, &
    0.22676438105580863683D+00, &
    0.11027597873677158176D+01, &
    0.91692778117386847344D+01, &
    0.15541656652426660966D+03, &
    0.26703582852084829694D+04, &
    0.86505880175304633906D+06, &
    0.11026046613094942620D+07, &
    0.22846209494153934787D+07, &
    0.42454972750111979449D+08, &
    0.48869614587997695539D+09, &
    0.56578651292431051863D+10, &
    0.76853203893832108948D+12, &
    0.14707396163259352103D+17, &
    0.29030785901035567967D+21, &
    0.58447515883904682813D+25, &
    0.11929750788892311875D+30 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
     -0.0312500000D+00, &
      0.1250000000D+00, &
     -0.5000000000D+00, &
      1.0000000000D+00, &
      2.0000000000D+00, &
     -4.0000000000D+00, &
      7.0000000000D+00, &
    -10.0000000000D+00, &
     16.0000000000D+00, &
     16.2500000000D+00, &
    -17.0000000000D+00, &
     20.0000000000D+00, &
     22.5000000000D+00, &
    -25.0000000000D+00, &
     30.0000000000D+00, &
    -40.0000000000D+00, &
     50.0000000000D+00, &
     60.0000000000D+00, &
    -70.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine student_cdf_values ( n_data, c, x, fx )

!*****************************************************************************80
!
!! STUDENT_CDF_VALUES returns some values of the Student CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = StudentTDistribution [ c ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) c
  real ( kind = 8 ), save, dimension ( n_max ) :: c_vec = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, &
    5.0D+00, 2.0D+00, 5.0D+00, 2.0D+00, &
    5.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, &
    5.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6000231200328521D+00, &
    0.6001080279134390D+00, &
    0.6001150934648930D+00, &
    0.6000995134721354D+00, &
    0.5999341989834830D+00, &
    0.7498859393137811D+00, &
    0.7500879487671045D+00, &
    0.9500004222186464D+00, &
    0.9499969138365968D+00, &
    0.9900012348724744D+00, &
    0.9900017619355059D+00, &
    0.9900004567580596D+00, &
    0.9900007637471291D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.325D+00, &
    0.289D+00, &
    0.277D+00, &
    0.271D+00, &
    0.267D+00, &
    0.816D+00, &
    0.727D+00, &
    2.920D+00, &
    2.015D+00, &
    6.965D+00, &
    4.541D+00, &
    3.747D+00, &
    3.365D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    c = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    c = c_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine student_noncentral_cdf_values ( n_data, df, lambda, x, fx )

!*****************************************************************************80
!
!! STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NoncentralStudentTDistribution [ df, lambda ]
!      CDF [ dist, x ]
!
!    Mathematica seems to have some difficulty computing this function
!    to the desired number of digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) DF, real ( kind = 8 ) LAMBDA, the parameters
!    of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30

  integer ( kind = 4 ) df
  integer ( kind = 4 ), save, dimension ( n_max ) :: df_vec = (/ &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
    15, 20, 25, &
     1,  2,  3, &
    10, 10, 10, &
    10, 10, 10, &
    10, 10, 10 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8975836176504333D+00, &
    0.9522670169D+00, &
    0.9711655571887813D+00, &
    0.8231218864D+00, &
    0.9049021510D+00, &
    0.9363471834D+00, &
    0.7301025986D+00, &
    0.8335594263D+00, &
    0.8774010255D+00, &
    0.5248571617D+00, &
    0.6293856597D+00, &
    0.6800271741D+00, &
    0.20590131975D+00, &
    0.2112148916D+00, &
    0.2074730718D+00, &
    0.9981130072D+00, &
    0.9994873850D+00, &
    0.9998391562D+00, &
    0.168610566972D+00, &
    0.16967950985D+00, &
    0.1701041003D+00, &
    0.9247683363D+00, &
    0.7483139269D+00, &
    0.4659802096D+00, &
    0.9761872541D+00, &
    0.8979689357D+00, &
    0.7181904627D+00, &
    0.9923658945D+00, &
    0.9610341649D+00, &
    0.8688007350D+00 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    4.0D+00, &
    4.0D+00, &
    4.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
    15.00D+00, &
    15.00D+00, &
    15.00D+00, &
     0.05D+00, &
     0.05D+00, &
     0.05D+00, &
     4.00D+00, &
     4.00D+00, &
     4.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     6.00D+00, &
     6.00D+00, &
     6.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    df = 0
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    df = df_vec(n_data)
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine subfactorial_values ( n_data, n, fn )

!*****************************************************************************80
!
!! SUBFACTORIAL_VALUES returns values of the subfactorial function.
!
!  Discussion:
!
!    The subfactorial function Subfactorial(N) counts the number of
!    permutations of N objects which leave no object unchanged.
!
!    Such a permutation is known as a derangement.
!
!    In Mathematica, the function can be evaluated by:
!
!      << DiscreteMath`CombinatorialFunctions`
!      Subfactorial[n]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, integer ( kind = 4 ) FN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  integer ( kind = 4 ), save, dimension ( n_max ) :: fn_vec = (/ &
            1, &
            0, &
            1, &
            2, &
            9, &
           44, &
          265, &
         1854, &
        14833, &
       133496, &
      1334961, &
     14684570, &
    176214841 /)
  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2,  3, &
     4,  5,  6,  7, &
     8,  9, 10, 11, &
    12 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    fn = 0
  else
    n = n_vec(n_data)
    fn = fn_vec(n_data)
  end if

  return
end
subroutine surten_values ( n_data, tc, sigma )

!*****************************************************************************80
!
!! SURTEN_VALUES returns some values of the surface tension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher, George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) SIGMA, the surface tension,
!    in Pascal * m = Newton / m.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 14

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    74.22D+00, &
    72.74D+00, &
    71.20D+00, &
    69.60D+00, &
    67.95D+00, &
    58.92D+00, &
    48.75D+00, &
    37.68D+00, &
    26.05D+00, &
    14.37D+00, &
     8.78D+00, &
     3.67D+00, &
     0.40D+00, &
     0.00D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
     10.000D+00, &
     20.000D+00, &
     30.000D+00, &
     40.000D+00, &
     50.000D+00, &
    100.000D+00, &
    150.000D+00, &
    200.000D+00, &
    250.000D+00, &
    300.000D+00, &
    325.000D+00, &
    350.000D+00, &
    370.000D+00, &
    373.976D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    tc = 0.0D+00
    sigma = 0.0D+00
  else
    tc = tc_vec(n_data)
    sigma = sigma_vec(n_data)
  end if

  return
end
subroutine synch1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! SYNCH1_VALUES returns some values of the synchrotron radiation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      SYNCH1(x) = x * integral ( x <= t < +oo ) K(5/3)(t) dt
!
!    where K(5/3) is a modified Bessel function of order 5/3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
      0.26514864547487397044D+00, &
      0.62050129979079045645D+00, &
      0.85112572132368011206D+00, &
      0.87081914687546885094D+00, &
      0.65142281535536396975D+00, &
      0.45064040920322354579D+00, &
      0.30163590285073940285D+00, &
      0.19814490804441305867D+00, &
      0.12856571000906381300D+00, &
      0.52827396697866818297D-01, &
      0.42139298471720305542D-01, &
      0.21248129774981984268D-01, &
      0.13400258907505536491D-01, &
      0.84260797314108699935D-02, &
      0.12884516186754671469D-02, &
      0.19223826430086897418D-03, &
      0.28221070834007689394D-04, &
      0.15548757973038189372D-05, &
      0.11968634456097453636D-07, &
      0.89564246772237127742D-10 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     12.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     25.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine synch2_values ( n_data, x, fx )

!*****************************************************************************80
!
!! SYNCH2_VALUES returns some values of the synchrotron radiation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      SYNCH2(x) = x * K(2/3)(x)
!
!    where K(2/3) is a modified Bessel function of order 2/3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.13430727275667378338D+00, &
    0.33485265272424176976D+00, &
    0.50404224110911078651D+00, &
    0.60296523236016785113D+00, &
    0.49447506210420826699D+00, &
    0.36036067860473360389D+00, &
    0.24967785497625662113D+00, &
    0.16813830542905833533D+00, &
    0.11117122348556549832D+00, &
    0.46923205826101330711D-01, &
    0.37624545861980001482D-01, &
    0.19222123172484106436D-01, &
    0.12209535343654701398D-01, &
    0.77249644268525771866D-02, &
    0.12029044213679269639D-02, &
    0.18161187569530204281D-03, &
    0.26884338006629353506D-04, &
    0.14942212731345828759D-05, &
    0.11607696854385161390D-07, &
    0.87362343746221526073D-10 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     12.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     25.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tan_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TAN_VALUES returns some values of the tangent function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Tan[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.00000000000000000000D+00, &
     0.26794919243112270647D+00, &
     0.54630248984379051326D+00, &
     0.57735026918962576451D+00, &
     1.0000000000000000000D+00, &
     1.5574077246549022305D+00, &
     1.7320508075688772935D+00, &
     3.7320508075688772935D+00, &
     7.5957541127251504405D+00, &
    15.257051688265539110D+00, &
    -2.1850398632615189916D+00, &
    -0.14254654307427780530D+00, &
     0.0000000000000000000D+00, &
     1.1578212823495775831D+00, &
    -3.3805150062465856370D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00000000000000000000D+00, &
    0.26179938779914943654D+00, &
    0.50000000000000000000D+00, &
    0.52359877559829887308D+00, &
    0.78539816339744830962D+00, &
    1.0000000000000000000D+00, &
    1.0471975511965977462D+00, &
    1.3089969389957471827D+00, &
    1.4398966328953219010D+00, &
    1.5053464798451092601D+00, &
    2.0000000000000000000D+00, &
    3.0000000000000000000D+00, &
    3.1415926535897932385D+00, &
    4.0000000000000000000D+00, &
    5.0000000000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tanh_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TANH_VALUES returns some values of the hyperbolic tangent function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Tanh[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
   -0.99990920426259513121D+00, &
   -0.76159415595576488812D+00, &
    0.00000000000000000000D+00, &
    0.099667994624955817118D+00, &
    0.19737532022490400074D+00, &
    0.29131261245159090582D+00, &
    0.37994896225522488527D+00, &
    0.46211715726000975850D+00, &
    0.53704956699803528586D+00, &
    0.60436777711716349631D+00, &
    0.66403677026784896368D+00, &
    0.71629787019902442081D+00, &
    0.76159415595576488812D+00, &
    0.96402758007581688395D+00, &
    0.99505475368673045133D+00, &
    0.99932929973906704379D+00, &
    0.99990920426259513121D+00, &
    0.99999999587769276362D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
   -5.0D+00, &
   -1.0D+00, &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
   10.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tau_values ( n_data, n, c )

!*****************************************************************************80
!
!! TAU_VALUES returns some values of the Tau function.
!
!  Discussion:
!
!    TAU(N) is the number of divisors of N, including 1 and N.
!
!    In Mathematica, the function can be evaluated by:
!
!      DivisorSigma[1,n]
!
!  First values:
!
!     N   TAU(N)
!
!     1    1
!     2    2
!     3    2
!     4    3
!     5    2
!     6    4
!     7    2
!     8    4
!     9    3
!    10    4
!    11    2
!    12    6
!    13    2
!    14    4
!    15    4
!    16    5
!    17    2
!    18    6
!    19    2
!    20    6
!
!  Formula:
!
!    If the prime factorization of N is
!
!      N = P1**E1 * P2**E2 * ... * PM**EM,
!
!    then
!
!      TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the Tau function.
!
!    Output, integer ( kind = 4 ) C, the value of the Tau function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
    1,  2,  2,  3,  2,  4,  2,  4,  3,  4, &
    2, 12, 12,  4, 18, 24,  2,  8, 14, 28 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     23,  72, 126, 226, 300, 480, 521, 610, 832, 960 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine thercon_values ( n_data, tc, p, lambda )

!*****************************************************************************80
!
!! THERCON_VALUES returns some values of the thermal conductivity.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher, George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) LAMBDA, the thermal conductivity, in
!    mW/(m degrees Kelvin).
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 35

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
    561.00D+00, &
    561.30D+00, &
    561.50D+00, &
    562.40D+00, &
    563.70D+00, &
    565.10D+00, &
    566.50D+00, &
    567.90D+00, &
    569.30D+00, &
    570.60D+00, &
    572.00D+00, &
    573.40D+00, &
    574.80D+00, &
    576.10D+00, &
    577.50D+00, &
    580.20D+00, &
    582.90D+00, &
    585.50D+00, &
    588.10D+00, &
    590.70D+00, &
    593.30D+00, &
    595.80D+00, &
    598.30D+00, &
    603.10D+00, &
    607.80D+00, &
    612.20D+00, &
    607.20D+00, &
    643.60D+00, &
    666.80D+00, &
     25.08D+00, &
     28.85D+00, &
     33.28D+00, &
     54.76D+00, &
     79.89D+00, &
    107.30D+00 /)
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
       1.0D+00, &
       5.0D+00, &
      10.0D+00, &
      25.0D+00, &
      50.0D+00, &
      75.0D+00, &
     100.0D+00, &
     125.0D+00, &
     150.0D+00, &
     175.0D+00, &
     200.0D+00, &
     225.0D+00, &
     250.0D+00, &
     275.0D+00, &
     300.0D+00, &
     350.0D+00, &
     400.0D+00, &
     450.0D+00, &
     500.0D+00, &
     550.0D+00, &
     600.0D+00, &
     650.0D+00, &
     700.0D+00, &
     800.0D+00, &
     900.0D+00, &
    1000.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
     25.0D+00, &
     50.0D+00, &
     75.0D+00, &
    100.0D+00, &
    150.0D+00, &
    200.0D+00, &
    400.0D+00, &
    600.0D+00, &
    800.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    tc = 0.0D+00
    p = 0.0D+00
    lambda = 0.0D+00
  else
    tc = tc_vec(n_data)
    p = p_vec(n_data)
    lambda = lambda_vec(n_data)
  end if

  return
end
subroutine three_j_values ( n_data, j1, j2, j3, m1, m2, m3, fx )

!*****************************************************************************80
!
!! THREE_J_VALUES returns some values of the Wigner 3J function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ThreeJSymbol[{j1,m1},{j2,m2},{j3,m3}]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) J1, J2, J3, M1, M2, M3, the arguments
!    of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 8

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.2788866755113585D+00, &
    -0.09534625892455923D+00, &
    -0.06741998624632421D+00, &
     0.1533110351679666D+00, &
    -0.1564465546936860D+00, &
     0.1099450412156551D+00, &
    -0.05536235693131719D+00, &
     0.01799835451137786D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) j1
  real ( kind = 8 ), save, dimension ( n_max ) :: j1_vec = (/ &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    6.0D+00, &
    7.0D+00, &
    8.0D+00 /)
  real ( kind = 8 ) j2
  real ( kind = 8 ), save, dimension ( n_max ) :: j2_vec = (/ &
    4.5D+00, &
    4.5D+00, &
    4.5D+00, &
    4.5D+00, &
    4.5D+00, &
    4.5D+00, &
    4.5D+00, &
    4.5D+00 /)
  real ( kind = 8 ) j3
  real ( kind = 8 ), save, dimension ( n_max ) :: j3_vec = (/ &
    3.5D+00, &
    3.5D+00, &
    3.5D+00, &
    3.5D+00, &
    3.5D+00, &
    3.5D+00, &
    3.5D+00, &
    3.5D+00 /)
  real ( kind = 8 ) m1
  real ( kind = 8 ), save, dimension ( n_max ) :: m1_vec = (/ &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00 /)
  real ( kind = 8 ) m2
  real ( kind = 8 ), save, dimension ( n_max ) :: m2_vec = (/ &
    -3.5D+00, &
    -3.5D+00, &
    -3.5D+00, &
    -3.5D+00, &
    -3.5D+00, &
    -3.5D+00, &
    -3.5D+00, &
    -3.5D+00 /)
  real ( kind = 8 ) m3
  real ( kind = 8 ), save, dimension ( n_max ) :: m3_vec = (/ &
    2.5D+00, &
    2.5D+00, &
    2.5D+00, &
    2.5D+00, &
    2.5D+00, &
    2.5D+00, &
    2.5D+00, &
    2.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    j1 = 0.0D+00
    j2 = 0.0D+00
    j3 = 0.0D+00
    m1 = 0.0D+00
    m2 = 0.0D+00
    m3 = 0.0D+00
    fx = 0.0D+00
  else
    j1 = j1_vec(n_data)
    j2 = j2_vec(n_data)
    j3 = j3_vec(n_data)
    m1 = m1_vec(n_data)
    m2 = m2_vec(n_data)
    m3 = m3_vec(n_data)
    fx = fx_vec(n_data)
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

  character ( len = 8 ) ampm
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
subroutine tran02_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRAN02_VALUES returns some values of the order 2 transportation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN02(x) = integral ( 0 <= t <= x ) t^2 exp(t) / ( exp(t) - 1 )^2 dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.19531247930394515480D-02, &
    0.31249152314331109004D-01, &
    0.12494577194783451032D+00, &
    0.49655363615640595865D+00, &
    0.97303256135517012845D+00, &
    0.14121978695932525805D+01, &
    0.18017185674405776809D+01, &
    0.21350385339277043015D+01, &
    0.24110500490169534620D+01, &
    0.28066664045631179931D+01, &
    0.28777421863296234131D+01, &
    0.30391706043438554330D+01, &
    0.31125074928667355940D+01, &
    0.31656687817738577185D+01, &
    0.32623520367816009184D+01, &
    0.32843291144979517358D+01, &
    0.32897895167775788137D+01, &
    0.32898672226665499687D+01, &
    0.32898681336064325400D+01, &
    0.32898681336964528724D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tran03_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRAN03_VALUES returns some values of the order 3 transportation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN03(x) = integral ( 0 <= t <= x ) t^3 * exp(t) / ( exp(t) - 1 )^2 dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.19073483296476379584D-05, &
    0.48826138243180786081D-03, &
    0.78074163848431205820D-02, &
    0.12370868718812031049D+00, &
    0.47984100657241749994D+00, &
    0.10269431622039754738D+01, &
    0.17063547219458658863D+01, &
    0.24539217444475937661D+01, &
    0.32106046629422467723D+01, &
    0.45792174372291563703D+01, &
    0.48722022832940370805D+01, &
    0.56143866138422732286D+01, &
    0.59984455864575470009D+01, &
    0.63033953673480961120D+01, &
    0.69579908688361166266D+01, &
    0.71503227120085929750D+01, &
    0.72110731475871876393D+01, &
    0.72123221966388461839D+01, &
    0.72123414161609465119D+01, &
    0.72123414189575656868D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tran04_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRAN04_VALUES returns some values of the order 4 transportation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN04(x) = integral ( 0 <= t <= x ) t^4 * exp(t) / ( exp(t) - 1 )^2 dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.24835263919461834041D-08, &
    0.10172029353616724881D-04, &
    0.65053332405940765479D-03, &
    0.41150448004155727767D-01, &
    0.31724404523442648241D+00, &
    0.10079442901142373591D+01, &
    0.22010881024333408363D+01, &
    0.38846508619156545210D+01, &
    0.59648223973714765245D+01, &
    0.10731932392998622219D+02, &
    0.11940028876819364777D+02, &
    0.15359784316882182982D+02, &
    0.17372587633093742893D+02, &
    0.19122976016053166969D+02, &
    0.23583979156921941515D+02, &
    0.25273667677030441733D+02, &
    0.25955198214572256372D+02, &
    0.25975350935212241910D+02, &
    0.25975757522084093747D+02, &
    0.25975757609067315288D+02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tran05_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRAN05_VALUES returns some values of the order 5 transportation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN05(x) = integral ( 0 <= t <= x ) t^5 * exp(t) / ( exp(t) - 1 )^2 dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.36379780361036116971D-11, &
    0.23840564453948442379D-06, &
    0.60982205372226969189D-04, &
    0.15410004586376649337D-01, &
    0.23661587923909478926D+00, &
    0.11198756851307629651D+01, &
    0.32292901663684049171D+01, &
    0.70362973105160654056D+01, &
    0.12770557691044159511D+02, &
    0.29488339015245845447D+02, &
    0.34471340540362254586D+02, &
    0.50263092218175187785D+02, &
    0.60819909101127165207D+02, &
    0.70873334429213460498D+02, &
    0.10147781242977788097D+03, &
    0.11638074540242071077D+03, &
    0.12409623901262967878D+03, &
    0.12442270155632550228D+03, &
    0.12443132790838589548D+03, &
    0.12443133061720432435D+03 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tran06_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRAN06_VALUES returns some values of the order 6 transportation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN06(x) = integral ( 0 <= t <= x ) t^6 * exp(t) / ( exp(t) - 1 )^2 dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.56843405953641209574D-14, &
    0.59601180165247401484D-08, &
    0.60978424397580572815D-05, &
    0.61578909866319494394D-02, &
    0.18854360275680840514D+00, &
    0.13319251347921659134D+01, &
    0.50857202271697616755D+01, &
    0.13729222365466557122D+02, &
    0.29579592481641441292D+02, &
    0.88600835706899853768D+02, &
    0.10916037113373004909D+03, &
    0.18224323749575359518D+03, &
    0.23765383125586756031D+03, &
    0.29543246745959381136D+03, &
    0.50681244381280455592D+03, &
    0.63878231134946125623D+03, &
    0.72699203556994876111D+03, &
    0.73230331643146851717D+03, &
    0.73248692015882096369D+03, &
    0.73248700462879996604D+03 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tran07_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRAN07_VALUES returns some values of the order 7 transportation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN07(x) = integral ( 0 <= t <= x ) t^7 * exp(t) / ( exp(t) - 1 )^2 dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.92518563327283409427D-17, &
    0.15521095556949867541D-09, &
    0.63516238373841716290D-06, &
    0.25638801246626135714D-02, &
    0.15665328993811649746D+00, &
    0.16538225039181097423D+01, &
    0.83763085709508211054D+01, &
    0.28078570717830763747D+02, &
    0.72009676046751991365D+02, &
    0.28174905701691911450D+03, &
    0.36660227975327792529D+03, &
    0.70556067982603601123D+03, &
    0.99661927562755629434D+03, &
    0.13288914430417403901D+04, &
    0.27987640273169129925D+04, &
    0.39721376409416504325D+04, &
    0.49913492839319899726D+04, &
    0.50781562639825019000D+04, &
    0.50820777202028708434D+04, &
    0.50820803580047164618D+04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tran08_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRAN08_VALUES returns some values of the order 8 transportation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN08(x) = integral ( 0 <= t <= x ) t^8 * exp(t) / ( exp(t) - 1 )^2 dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.15488598634539359463D-19, &
    0.41574269117845953797D-11, &
    0.68050651245227411689D-07, &
    0.10981703519563009836D-02, &
    0.13396432776187883834D+00, &
    0.21153387806998617182D+01, &
    0.14227877028750735641D+02, &
    0.59312061431647843226D+02, &
    0.18139614577043147745D+03, &
    0.93148001928992220863D+03, &
    0.12817928112604611804D+04, &
    0.28572838386329242218D+04, &
    0.43872971687877730010D+04, &
    0.62993229139406657611D+04, &
    0.16589426277154888511D+05, &
    0.27064780798797398935D+05, &
    0.38974556062543661284D+05, &
    0.40400240716905025786D+05, &
    0.40484316504120655568D+05, &
    0.40484399001892184901D+05 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tran09_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRAN09_VALUES returns some values of the order 9 transportation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN09(x) = integral ( 0 <= t <= x ) t^9 * exp(t) / ( exp(t) - 1 )^2 dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Allan McLeod,
!    Algorithm 757:
!    MISCFUN: A software package to compute uncommon special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.26469772870084897671D-22, &
    0.11367943653594246210D-12, &
    0.74428246255329800255D-08, &
    0.48022728485415366194D-03, &
    0.11700243014358676725D+00, &
    0.27648973910899914391D+01, &
    0.24716631405829192997D+02, &
    0.12827119828849828583D+03, &
    0.46842894800662208986D+03, &
    0.31673967371627895718D+04, &
    0.46140886546630195390D+04, &
    0.11952718545392302185D+05, &
    0.20001612666477027728D+05, &
    0.31011073271851366554D+05, &
    0.10352949905541130133D+06, &
    0.19743173017140591390D+06, &
    0.33826030414658460679D+06, &
    0.36179607036750755227D+06, &
    0.36360622124777561525D+06, &
    0.36360880558827162725D+06 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
      0.0019531250D+00, &
      0.0312500000D+00, &
      0.1250000000D+00, &
      0.5000000000D+00, &
      1.0000000000D+00, &
      1.5000000000D+00, &
      2.0000000000D+00, &
      2.5000000000D+00, &
      3.0000000000D+00, &
      4.0000000000D+00, &
      4.2500000000D+00, &
      5.0000000000D+00, &
      5.5000000000D+00, &
      6.0000000000D+00, &
      8.0000000000D+00, &
     10.0000000000D+00, &
     15.0000000000D+00, &
     20.0000000000D+00, &
     30.0000000000D+00, &
     50.0000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine trigamma_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRIGAMMA_VALUES returns some values of the TriGamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyGamma[1,x]
!
!    TriGamma(X) = d^2 ln ( Gamma ( X ) ) / d X^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1644934066848226D+01, &
     0.1433299150792759D+01, &
     0.1267377205423779D+01, &
     0.1134253434996619D+01, &
     0.1025356590529597D+01, &
     0.9348022005446793D+00, &
     0.8584318931245799D+00, &
     0.7932328301639984D+00, &
     0.7369741375017002D+00, &
     0.6879720582426356D+00, &
     0.6449340668482264D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine tsat_values ( n_data, p, tc )

!*****************************************************************************80
!
!! TSAT_VALUES returns some values of the saturation temperature.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher, George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) TC, the saturation temperature, in
!    degrees Celsius.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
      0.0061173D+00, &
      0.012D+00, &
      0.025D+00, &
      0.055D+00, &
      0.080D+00, &
      0.110D+00, &
      0.160D+00, &
      0.250D+00, &
      0.500D+00, &
      0.750D+00, &
      1.000D+00, &
      1.500D+00, &
      2.000D+00, &
      5.000D+00, &
     10.000D+00, &
     20.000D+00, &
     50.000D+00, &
    100.000D+00, &
    200.000D+00, &
    220.550D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
      0.010D+00, &
      9.655D+00, &
     21.080D+00, &
     34.589D+00, &
     41.518D+00, &
     47.695D+00, &
     55.327D+00, &
     64.980D+00, &
     81.339D+00, &
     91.783D+00, &
     99.632D+00, &
    111.378D+00, &
    120.443D+00, &
    151.866D+00, &
    179.916D+00, &
    212.417D+00, &
    263.977D+00, &
    311.031D+00, &
    365.800D+00, &
    373.976D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    p = 0.0D+00
    tc = 0.0D+00
  else
    p = p_vec(n_data)
    tc = tc_vec(n_data)
  end if

  return
end
subroutine van_der_corput_values ( n_data, base, seed, value )

!*****************************************************************************80
!
!! VAN_DER_CORPUT_VALUES returns some values of the van der Corput sequence.
!
!  Discussion:
!
!    The van der Corput sequence is often used to generate a "subrandom"
!    sequence of points which have a better covering property
!    than pseudorandom points.
!
!    The van der Corput sequence generates a sequence of points in [0,1]
!    which (theoretically) never repeats.  Except for SEED = 0, the
!    elements of the van der Corput sequence are strictly between 0 and 1.
!
!    The van der Corput sequence writes an integer in a given base B,
!    and then its digits are "reflected" about the decimal point.
!    This maps the numbers from 1 to N into a set of numbers in [0,1],
!    which are especially nicely distributed if N is one less
!    than a power of the base.
!
!    Hammersley suggested generating a set of N nicely distributed
!    points in two dimensions by setting the first component of the
!    Ith point to I/N, and the second to the van der Corput
!    value of I in base 2.
!
!    Halton suggested that in many cases, you might not know the number
!    of points you were generating, so Hammersley's formulation was
!    not ideal.  Instead, he suggested that to generate a nicely
!    distributed sequence of points in M dimensions, you simply
!    choose the first M primes, P(1:M), and then for the J-th component of
!    the I-th point in the sequence, you compute the van der Corput
!    value of I in base P(J).
!
!    Thus, to generate a Halton sequence in a 2 dimensional space,
!    it is typical practice to generate a pair of van der Corput sequences,
!    the first with prime base 2, the second with prime base 3.
!    Similarly, by using the first K primes, a suitable sequence
!    in K-dimensional space can be generated.
!
!    The generation is quite simple.  Given an integer SEED, the expansion
!    of SEED in base BASE is generated.  Then, essentially, the result R
!    is generated by writing a decimal point followed by the digits of
!    the expansion of SEED, in reverse order.  This decimal value is actually
!    still in base BASE, so it must be properly interpreted to generate
!    a usable value.
!
!  Example:
!
!    BASE = 2
!
!    SEED     SEED      van der Corput
!    decimal  binary    binary   decimal
!    -------  ------    ------   -------
!        0  =     0  =>  .0     = 0.0D+00
!        1  =     1  =>  .1     = 0.5
!        2  =    10  =>  .01    = 0.25
!        3  =    11  =>  .11    = 0.75
!        4  =   100  =>  .001   = 0.125
!        5  =   101  =>  .101   = 0.625
!        6  =   110  =>  .011   = 0.375
!        7  =   111  =>  .111   = 0.875
!        8  =  1000  =>  .0001  = 0.0625
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, pages 84-90, 1960.
!
!    John Hammersley,
!    Monte Carlo methods for solving multivariable problems,
!    Proceedings of the New York Academy of Science,
!    Volume 86, pages 844-874, 1960.
!
!    Johannes van der Corput,
!    Verteilungsfunktionen,
!    Proc Akad Amsterdam,
!    Volume 38, 1935,
!    Volume 39, 1936.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) BASE, the base of the sequence.
!
!    Output, integer ( kind = 4 ) SEED, the index of the element of
!    the sequence.
!
!    Output, real ( kind = 8 ) VALUE, the value of the SEED-th element of the
!    van der Corput sequence in base BASE.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 75

  integer ( kind = 4 ) base
  integer ( kind = 4 ), save, dimension ( n_max ) :: base_vec = (/ &
     2,   2,   2,   2,   2, &
     2,   2,   2,   2,   3, &
     3,   3,   3,   3,   3, &
     3,   3,   3,   4,   4, &
     4,   4,   4,   4,   4, &
     4,   4,   2,   3,   4, &
     5,   7,  11,  13,   2, &
     3,   4,   5,   7,  11, &
    13,   2,   3,   4,   5, &
     7,  11,  13,   2,   3, &
     4,   5,   7,  11,  13, &
    29,  29,  29,  29,  29, &
    71,  71,  71,  71,  71, &
   173, 173, 173, 173, 173, &
   409, 409, 409, 409, 409 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save, dimension ( n_max ) :: seed_vec = (/ &
        0,     1,     2,     3,     4, &
        5,     6,     7,     8,     0, &
        1,     2,     3,     4,     5, &
        6,     7,     8,     0,     1, &
        2,     3,     4,     5,     6, &
        7,     8,    10,    10,    10, &
       10,    10,    10,    10,   100, &
      100,   100,   100,   100,   100, &
      100,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000, 10000, 10000, &
    10000, 10000, 10000, 10000, 10000, &
     1000,  1001,  1002,  1003,  1004, &
     1000,  1001,  1002,  1003,  1004, &
     1000,  1001,  1002,  1003,  1004, &
     1000,  1001,  1002,  1003,  1004 /)
  real ( kind = 8 ) value
  real ( kind = 8 ), save, dimension ( n_max ) :: value_vec = (/ &
    0.0000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2500000000000000D+00, &
    0.7500000000000000D+00, &
    0.1250000000000000D+00, &
    0.6250000000000000D+00, &
    0.3750000000000000D+00, &
    0.8750000000000000D+00, &
    0.0625000000000000D+00, &
    0.0000000000000000D+00, &
    0.3333333333333333D+00, &
    0.6666666666666666D+00, &
    0.1111111111111111D+00, &
    0.4444444444444444D+00, &
    0.7777777777777777D+00, &
    0.2222222222222222D+00, &
    0.5555555555555556D+00, &
    0.8888888888888888D+00, &
    0.0000000000000000D+00, &
    0.2500000000000000D+00, &
    0.5000000000000000D+00, &
    0.7500000000000000D+00, &
    0.0625000000000000D+00, &
    0.3125000000000000D+00, &
    0.5625000000000000D+00, &
    0.8125000000000000D+00, &
    0.1250000000000000D+00, &
    0.3125000000000000D+00, &
    0.3703703703703703D+00, &
    0.6250000000000000D+00, &
    0.0800000000000000D+00, &
    0.4489795918367347D+00, &
    0.9090909090909092D+00, &
    0.7692307692307693D+00, &
    0.1484375000000000D+00, &
    0.4115226337448559D+00, &
    0.0976562500000000D+00, &
    0.0320000000000000D+00, &
    0.2915451895043731D+00, &
    0.1652892561983471D+00, &
    0.7337278106508875D+00, &
    0.0927734375000000D+00, &
    0.3475080018289895D+00, &
    0.1708984375000000D+00, &
    0.0051200000000000D+00, &
    0.9162848812994586D+00, &
    0.9316303531179565D+00, &
    0.9904415111515704D+00, &
    0.0347290039062500D+00, &
    0.3861200020322105D+00, &
    0.0189208984375000D+00, &
    0.0005120000000000D+00, &
    0.5749985125245433D+00, &
    0.1529950140017758D+00, &
    0.2459297643639929D+00, &
    0.4887449259912255D+00, &
    0.5232276846119153D+00, &
    0.5577104432326049D+00, &
    0.5921932018532945D+00, &
    0.6266759604739842D+00, &
    0.0872842689942472D+00, &
    0.1013687760365007D+00, &
    0.1154532830787542D+00, &
    0.1295377901210077D+00, &
    0.1436222971632613D+00, &
    0.7805138828560928D+00, &
    0.7862942296769020D+00, &
    0.7920745764977113D+00, &
    0.7978549233185205D+00, &
    0.8036352701393298D+00, &
    0.4449997309915651D+00, &
    0.4474447187666262D+00, &
    0.4498897065416874D+00, &
    0.4523346943167484D+00, &
    0.4547796820918096D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    base = 0
    seed = 0
    value = 0.0D+00
  else
    base = base_vec(n_data)
    seed = seed_vec(n_data)
    value = value_vec(n_data)
  end if

  return
end
subroutine viscosity_values ( n_data, tc, p, eta )

!*****************************************************************************80
!
!! VISCOSITY_VALUES returns some values of the viscosity function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher, George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    ISBN: 0-89116-353-0,
!    LC: TJ270.H3.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) ETA, the viscosity, in MegaPascal seconds.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 34

  real ( kind = 8 ) eta
  real ( kind = 8 ), save, dimension ( n_max ) :: eta_vec = (/ &
    1792.0D+00, &
    1791.0D+00, &
    1790.0D+00, &
    1786.0D+00, &
    1780.0D+00, &
    1775.0D+00, &
    1769.0D+00, &
    1764.0D+00, &
    1759.0D+00, &
    1754.0D+00, &
    1749.0D+00, &
    1744.0D+00, &
    1739.0D+00, &
    1735.0D+00, &
    1731.0D+00, &
    1722.0D+00, &
    1714.0D+00, &
    1707.0D+00, &
    1700.0D+00, &
    1694.0D+00, &
    1687.0D+00, &
    1682.0D+00, &
    1676.0D+00, &
    1667.0D+00, &
    1659.0D+00, &
    1653.0D+00, &
     890.8D+00, &
     547.1D+00, &
     378.4D+00, &
     12.28D+00, &
     16.18D+00, &
     24.45D+00, &
     32.61D+00, &
     40.38D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
       1.0D+00, &
       5.0D+00, &
      10.0D+00, &
      25.0D+00, &
      50.0D+00, &
      75.0D+00, &
     100.0D+00, &
     125.0D+00, &
     150.0D+00, &
     175.0D+00, &
     200.0D+00, &
     225.0D+00, &
     250.0D+00, &
     275.0D+00, &
     300.0D+00, &
     350.0D+00, &
     400.0D+00, &
     450.0D+00, &
     500.0D+00, &
     550.0D+00, &
     600.0D+00, &
     650.0D+00, &
     700.0D+00, &
     800.0D+00, &
     900.0D+00, &
    1000.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00, &
       1.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_max ) :: tc_vec = (/ &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
      0.0D+00, &
     25.0D+00, &
     50.0D+00, &
     75.0D+00, &
    100.0D+00, &
    200.0D+00, &
    400.0D+00, &
    600.0D+00, &
    800.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    tc = 0.0D+00
    p = 0.0D+00
    eta = 0.0D+00
  else
    tc = tc_vec(n_data)
    p = p_vec(n_data)
    eta = eta_vec(n_data)
  end if

  return
end
subroutine von_mises_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! VON_MISES_CDF_VALUES returns some values of the von Mises CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kanti Mardia, Peter Jupp,
!    Directional Statistics,
!    Wiley, 2000,
!    LC: QA276.M335
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Output, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 23

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
    -0.2D+01, &
    -0.1D+01, &
     0.0D+01, &
     0.1D+01, &
     0.2D+01, &
     0.3D+01, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.3D+01, &
     0.3D+01, &
     0.3D+01, &
     0.3D+01, &
     0.3D+01, &
     0.3D+01, &
     0.0D+00, &
     0.1D+01, &
     0.2D+01, &
     0.3D+01, &
     0.4D+01, &
     0.5D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.2535089956281180D-01, &
    0.1097539041177346D+00, &
    0.5000000000000000D+00, &
    0.8043381312498558D+00, &
    0.9417460124555197D+00, &
    0.5000000000000000D+00, &
    0.6018204118446155D+00, &
    0.6959356933122230D+00, &
    0.7765935901304593D+00, &
    0.8410725934916615D+00, &
    0.8895777369550366D+00, &
    0.9960322705517925D+00, &
    0.9404336090170247D+00, &
    0.5000000000000000D+00, &
    0.5956639098297530D-01, &
    0.3967729448207649D-02, &
    0.2321953958111930D-03, &
    0.6250000000000000D+00, &
    0.7438406999109122D+00, &
    0.8369224904294019D+00, &
    0.8941711407897124D+00, &
    0.9291058600568743D+00, &
    0.9514289900655436D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.2617993977991494D+01, &
    -0.1570796326794897D+01, &
     0.0000000000000000D+00, &
     0.1047197551196598D+01, &
     0.2094395102393195D+01, &
     0.1000000000000000D+01, &
     0.1200000000000000D+01, &
     0.1400000000000000D+01, &
     0.1600000000000000D+01, &
     0.1800000000000000D+01, &
     0.2000000000000000D+01, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine weekday_values ( n_data, y, m, d, w )

!*****************************************************************************80
!
!! WEEKDAY_VALUES returns the day of the week for various dates.
!
!  Discussion:
!
!    The CE or Common Era calendar is used, under the
!    hybrid Julian/Gregorian Calendar, with a transition from Julian
!    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
!
!    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
!    years BC/BCE are indicated by a negative year value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Reingold, Nachum Dershowitz,
!    Calendrical Calculations: The Millennium Edition,
!    Cambridge University Press, 2001,
!    ISBN: 0 521 77752 6
!    LC: CE12.R45.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) Y, M, D, the Common Era date.
!
!    Output, integer ( kind = 8 ) W, the day of the week.  Sunday = 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 34

  integer ( kind = 4 ) d
  integer ( kind = 4 ), save, dimension ( n_max ) :: d_vec = (/ &
    30, &
     8, &
    26, &
     3, &
     7, &
    18, &
     7, &
    19, &
    14, &
    18, &
    16, &
     3, &
    26, &
    20, &
     4, &
    25, &
    31, &
     9, &
    24, &
    10, &
    30, &
    24, &
    19, &
     2, &
    27, &
    19, &
    25, &
    29, &
    19, &
     7, &
    17, &
    25, &
    10, &
    18 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
     7, &
    12, &
     9, &
    10, &
     1, &
     5, &
    11, &
     4, &
    10, &
     5, &
     3, &
     3, &
     3, &
     4, &
     6, &
     1, &
     3, &
     9, &
     2, &
     6, &
     6, &
     7, &
     6, &
     8, &
     3, &
     4, &
     8, &
     9, &
     4, &
    10, &
     3, &
     2, &
    11, &
     7 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) w
  integer ( kind = 4 ), save, dimension ( n_max ) :: w_vec = (/ &
    1, &
    4, &
    4, &
    1, &
    4, &
    2, &
    7, &
    1, &
    7, &
    1, &
    6, &
    7, &
    6, &
    1, &
    1, &
    4, &
    7, &
    7, &
    7, &
    4, &
    1, &
    6, &
    1, &
    2, &
    4, &
    1, &
    1, &
    2, &
    2, &
    5, &
    3, &
    1, &
    4, &
    1 /)
  integer ( kind = 4 ) y
  integer ( kind = 4 ), save, dimension ( n_max ) :: y_vec = (/ &
    - 587, &
    - 169, &
       70, &
      135, &
      470, &
      576, &
      694, &
     1013, &
     1066, &
     1096, &
     1190, &
     1240, &
     1288, &
     1298, &
     1391, &
     1436, &
     1492, &
     1553, &
     1560, &
     1648, &
     1680, &
     1716, &
     1768, &
     1819, &
     1839, &
     1903, &
     1929, &
     1941, &
     1943, &
     1943, &
     1992, &
     1996, &
     2038, &
     2094 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    y = 0
    m = 0
    d = 0
    w = 0
  else
    y = y_vec(n_data)
    m = m_vec(n_data)
    d = d_vec(n_data)
    w = w_vec(n_data)
  end if

  return
end
subroutine weibull_cdf_values ( n_data, alpha, beta, x, fx )

!*****************************************************************************80
!
!! WEIBULL_CDF_VALUES returns some values of the Weibull CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = WeibullDistribution [ alpha, beta ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) ALPHA, the first parameter of the distribution.
!
!    Output, real ( kind = 8 ) BETA, the second parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) alpha
  real ( kind = 8 ), save, dimension ( n_max ) :: alpha_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  real ( kind = 8 ) beta
  real ( kind = 8 ), save, dimension ( n_max ) :: beta_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8646647167633873D+00, &
    0.9816843611112658D+00, &
    0.9975212478233336D+00, &
    0.9996645373720975D+00, &
    0.6321205588285577D+00, &
    0.4865828809674080D+00, &
    0.3934693402873666D+00, &
    0.3296799539643607D+00, &
    0.8946007754381357D+00, &
    0.9657818816883340D+00, &
    0.9936702845725143D+00, &
    0.9994964109502630D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    alpha = 0.0D+00
    beta = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    alpha = alpha_vec(n_data)
    beta = beta_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine zeta_values ( n_data, n, zeta )

!*****************************************************************************80
!
!! ZETA_VALUES returns some values of the Riemann Zeta function.
!
!  Discussion:
!
!    ZETA(N) = sum ( 1 <= I < +oo ) 1 / I^N
!
!    In Mathematica, the function can be evaluated by:
!
!      Zeta[n]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the Zeta function.
!
!    Output, real ( kind = 8 ) ZETA, the value of the Zeta function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     2, &
     3, &
     4, &
     5, &
     6, &
     7, &
     8, &
     9, &
    10, &
    11, &
    12, &
    16, &
    20, &
    30, &
    40 /)
  real ( kind = 8 ) zeta
  real ( kind = 8 ), save, dimension ( n_max ) :: zeta_vec = (/ &
    0.164493406684822643647D+01, &
    0.120205690315959428540D+01, &
    0.108232323371113819152D+01, &
    0.103692775514336992633D+01, &
    0.101734306198444913971D+01, &
    0.100834927738192282684D+01, &
    0.100407735619794433939D+01, &
    0.100200839292608221442D+01, &
    0.100099457512781808534D+01, &
    0.100049418860411946456D+01, &
    0.100024608655330804830D+01, &
    0.100001528225940865187D+01, &
    0.100000095396203387280D+01, &
    0.100000000093132743242D+01, &
    0.100000000000090949478D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    zeta = 0.0D+00
  else
    n = n_vec(n_data)
    zeta = zeta_vec(n_data)
  end if

  return
end
