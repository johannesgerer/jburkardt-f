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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
function alngam ( x )

!*****************************************************************************80
!
!! ALNGAM computes the log of the absolute value of the Gamma function.
!
!  Discussion:
!
!    The Gamma function is defined as
!
!      GAMMA(Z) = INTEGRAL ( 0 <= T < Infinity) T**(Z-1) EXP ( -T ) DT
!
!    If Z is a positive integer ( kind = 4 ), GAMMA(Z) = (Z-1)!, the factorial.
!
!    There is a special value:
!
!      GAMMA(0.5) = SQRT ( PI ).
!
!  Modified:
!
!    31 May 2000
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the gamma function.
!
!    Output, real ( kind = 8 ) ALNGAM, the logarithm of the absolute
!    value of GAMMA(X).
!
  implicit none

  real ( kind = 8 ) alngam
  real ( kind = 8 ) d9lgmc
  real ( kind = 8 ), save :: dxrel = 0.0D+00
  real ( kind = 8 ) gamma
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sinpiy
  real ( kind = 8 ), parameter :: sq2pil = 0.91893853320467274D+00
  real ( kind = 8 ), parameter :: sqpi2l = 0.22579135264472743D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xmax = 0.0D+00
  real ( kind = 8 ) y

  if ( xmax == 0.0D+00 ) then
    xmax = huge ( xmax ) / log ( huge ( xmax ) )
    dxrel = sqrt ( epsilon ( dxrel ) )
  end if

  y = abs ( x )

  if ( y <= 10.0D+00 ) then
    alngam = log ( abs ( gamma ( x ) ) )
    return
  end if

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ALNGAM - Fatal error!'
    write ( *, '(a)' ) '  |X| is so big that ALNGAM will overflow.'
    stop
  end if

  if ( 0.0D+00 < x ) then
    alngam = sq2pil + ( x - 0.5D+00 ) * log ( x ) - x + d9lgmc ( y )
    return
  end if

  sinpiy = abs ( sin ( pi * y ) )

  if ( sinpiy == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ALNGAM - Fatal error!'
    write ( *, '(a)' ) '  X is a negative integer ( kind = 4 ).'
    stop
  end if

  if ( abs ( ( x - real ( int ( x - 0.5D+00 ), kind = 8 ) ) / x ) < dxrel ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ALNGAM - Warning:'
    write ( *, '(a)' ) '  Answer has less than half usual precision.'
    write ( *, '(a)' ) '  X is very near a negative integer ( kind = 4 ).'
  end if

  alngam = sqpi2l + ( x - 0.5D+00 ) * log ( y ) - x - log ( sinpiy ) &
    - d9lgmc ( y )

  return
end
subroutine asyjy ( funjy, x, fnu, flgjy, in, y, wk, iflw )

!*****************************************************************************80
!
!! ASYJY computes high order Bessel functions J and Y.
!
!  Discussion:
!
!    ASYJY implements the uniform asymptotic expansion of
!    the J and Y Bessel functions for 35 <= FNU and 0.0 < X.
!
!    The forms are identical except for a change
!    in sign of some of the terms.  This change in sign is
!    accomplished by means of the flag FLGJY = 1 or -1.
!
!    On FLGJY = 1 the Airy functions AI(X) and DAI(X) are
!    supplied by the external function JAIRY, and on
!    FLGJY = -1 the Airy functions BI(X) and DBI(X) are
!    supplied by the external funtion YAIRY.
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Donald Amos
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, external FUNJY, is the function JAIRY or YAIRY.
!
!    Input, real ( kind = 8 ) X, the argument, which must be greater than 0.
!
!    Input, real ( kind = 8 ) FNU, the order of the first Bessel function.
!    FNU is generally at least 35.
!
!    Input, real ( kind = 8 ) FLGJY, a selection flag
!     1.0D+00 gives the J function
!    -1.0D+00 gives the Y function
!
!    Input, integer ( kind = 4 ) IN, the number of functions desired, which should be
!    1 or 2.
!
!    Output, real ( kind = 8 ) Y(IN), contains the desired function values.
!
!    Output, integer ( kind = 4 ) IFLW, a flag indicating underflow or overflow
!    return variables for BESJ only.
!
!    Output, real ( kind = 8 ) WK(7), contains the following values:
!
!      wk(1) = 1 - (x/fnu)**2 = w**2
!      wk(2) = sqrt ( abs ( wk(1) ) )
!      wk(3) = abs ( wk(2) - atan ( wk(2) ) )  or
!              abs ( ln((1 + wk(2) )/ ( x / fnu ) ) - wk(2))
!            = abs ( (2/3)*zeta**(3/2))
!      wk(4) = fnu*wk(3)
!      wk(5) = (1.5*wk(3) * fnu)**(1/3) = sqrt ( zeta ) * fnu**(1/3)
!      wk(6) = sign ( 1.0, w**2 ) * wk(5)**2
!            = sign ( 1.0, w**2 ) * zeta * fnu**(2/3)
!      wk(7) = fnu**(1/3)
!
  implicit none

  real ( kind = 8 ) abw2
  real ( kind = 8 ) akm
  real ( kind = 8 ) alfa(26,4)
  real ( kind = 8 ) alfa1(26,2)
  real ( kind = 8 ) alfa2(26,2)
  real ( kind = 8 ) ap
  real ( kind = 8 ), parameter, dimension ( 8 ) :: ar = (/ &
    8.35503472222222D-02, 1.28226574556327D-01, &
    2.91849026464140D-01, 8.81627267443758D-01, 3.32140828186277D+00, &
    1.49957629868626D+01, 7.89230130115865D+01, 4.74451538868264D+02 /)
  real ( kind = 8 ) asum
  real ( kind = 8 ) az
  real ( kind = 8 ) beta(26,5)
  real ( kind = 8 ) beta1(26,2)
  real ( kind = 8 ) beta2(26,2)
  real ( kind = 8 ) beta3(26,1)
  real ( kind = 8 ) br(10)
  real ( kind = 8 ) bsum
  real ( kind = 8 ) c(65)
  real ( kind = 8 ), parameter :: con1 = 6.66666666666667D-01
  real ( kind = 8 ), parameter :: con2 = 3.33333333333333D-01
  real ( kind = 8 ), parameter :: con548 = 1.04166666666667D-01
  real ( kind = 8 ) cr(10)
  real ( kind = 8 ) crz32
  real ( kind = 8 ) d1mach
  real ( kind = 8 ) dfi
  real ( kind = 8 ) elim
  real ( kind = 8 ) dr(10)
  real ( kind = 8 ) fi
  real ( kind = 8 ) flgjy
  real ( kind = 8 ) fn
  real ( kind = 8 ) fnu
  real ( kind = 8 ) fn2
  external funjy
  real ( kind = 8 ) gama(26)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) iflw
  integer ( kind = 4 ) in
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jn
  integer ( kind = 4 ) jr
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kb
  integer ( kind = 4 ) klast
  integer ( kind = 4 ) kmax(5)
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) ks
  integer ( kind = 4 ) ksp1
  integer ( kind = 4 ) kstemp
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrp1
  real ( kind = 8 ) phi
  real ( kind = 8 ) rcz
  real ( kind = 8 ) rden
  real ( kind = 8 ) relb
  real ( kind = 8 ) rfn2
  real ( kind = 8 ) rtz
  real ( kind = 8 ) rzden
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) suma
  real ( kind = 8 ) sumb
  real ( kind = 8 ) s1
  real ( kind = 8 ) ta
  real ( kind = 8 ) tau
  real ( kind = 8 ) tb
  real ( kind = 8 ) tfn
  real ( kind = 8 ) tol
  real ( kind = 8 ), save :: tols = -6.90775527898214D+00
  real ( kind = 8 ) t2
  real ( kind = 8 ) upol(10)
  real ( kind = 8 ) wk(*)
  real ( kind = 8 ) x
  real ( kind = 8 ) xx
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) z
  real ( kind = 8 ) z32

  equivalence (alfa(1,1),alfa1(1,1))
  equivalence (alfa(1,3),alfa2(1,1))
  equivalence (beta(1,1),beta1(1,1))
  equivalence (beta(1,3),beta2(1,1))
  equivalence (beta(1,5),beta3(1,1))

  data  br(1), br(2), br(3), br(4), br(5), br(6), br(7), br(8), br(9), br(10) &
    /-1.45833333333333D-01,-9.87413194444444D-02, &
     -1.43312053915895D-01,-3.17227202678414D-01,-9.42429147957120D-01, &
     -3.51120304082635D+00,-1.57272636203680D+01,-8.22814390971859D+01, &
     -4.92355370523671D+02,-3.31621856854797D+03/

  data c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10), &
     c(11), c(12), c(13), c(14), c(15), c(16), c(17), c(18), &
     c(19), c(20), c(21), c(22), c(23), c(24)/ &
     -2.08333333333333D-01,        1.25000000000000D-01, &
      3.34201388888889D-01,       -4.01041666666667D-01, &
      7.03125000000000D-02,       -1.02581259645062D+00, &
      1.84646267361111D+00,       -8.91210937500000D-01, &
      7.32421875000000D-02,        4.66958442342625D+00, &
     -1.12070026162230D+01,        8.78912353515625D+00, &
     -2.36408691406250D+00,        1.12152099609375D-01, &
     -2.82120725582002D+01,        8.46362176746007D+01, &
     -9.18182415432400D+01,        4.25349987453885D+01, &
     -7.36879435947963D+00,        2.27108001708984D-01, &
      2.12570130039217D+02,       -7.65252468141182D+02, &
      1.05999045252800D+03,       -6.99579627376133D+02/

  data c(25), c(26), c(27), c(28), c(29), c(30), c(31), c(32), &
          c(33), c(34), c(35), c(36), c(37), c(38), c(39), c(40), &
          c(41), c(42), c(43), c(44), c(45), c(46), c(47), c(48)/ &
             2.18190511744212D+02,       -2.64914304869516D+01, &
             5.72501420974731D-01,       -1.91945766231841D+03, &
             8.06172218173731D+03,       -1.35865500064341D+04, &
             1.16553933368645D+04,       -5.30564697861340D+03, &
             1.20090291321635D+03,       -1.08090919788395D+02, &
             1.72772750258446D+00,        2.02042913309661D+04, &
            -9.69805983886375D+04,        1.92547001232532D+05, &
            -2.03400177280416D+05,        1.22200464983017D+05, &
            -4.11926549688976D+04,        7.10951430248936D+03, &
            -4.93915304773088D+02,        6.07404200127348D+00, &
            -2.42919187900551D+05,        1.31176361466298D+06, &
            -2.99801591853811D+06,        3.76327129765640D+06/

  data c(49), c(50), c(51), c(52), c(53), c(54), c(55), c(56), &
          c(57), c(58), c(59), c(60), c(61), c(62), c(63), c(64), &
          c(65)/ &
            -2.81356322658653D+06,        1.26836527332162D+06, &
            -3.31645172484564D+05,        4.52187689813627D+04, &
            -2.49983048181121D+03,        2.43805296995561D+01, &
             3.28446985307204D+06,       -1.97068191184322D+07, &
             5.09526024926646D+07,       -7.41051482115327D+07, &
             6.63445122747290D+07,       -3.75671766607634D+07, &
             1.32887671664218D+07,       -2.78561812808645D+06, &
             3.08186404612662D+05,       -1.38860897537170D+04, &
             1.10017140269247D+02/

  data alfa1(1,1), alfa1(2,1), alfa1(3,1), alfa1(4,1), alfa1(5,1), &
          alfa1(6,1), alfa1(7,1), alfa1(8,1), alfa1(9,1), alfa1(10,1), &
          alfa1(11,1),alfa1(12,1),alfa1(13,1),alfa1(14,1),alfa1(15,1), &
          alfa1(16,1),alfa1(17,1),alfa1(18,1),alfa1(19,1),alfa1(20,1), &
          alfa1(21,1),alfa1(22,1),alfa1(23,1),alfa1(24,1),alfa1(25,1), &
          alfa1(26,1)     /-4.44444444444444D-03,-9.22077922077922D-04, &
     -8.84892884892885D-05, 1.65927687832450D-04, 2.46691372741793D-04, &
      2.65995589346255D-04, 2.61824297061501D-04, 2.48730437344656D-04, &
      2.32721040083232D-04, 2.16362485712365D-04, 2.00738858762752D-04, &
      1.86267636637545D-04, 1.73060775917876D-04, 1.61091705929016D-04, &
      1.50274774160908D-04, 1.40503497391270D-04, 1.31668816545923D-04, &
      1.23667445598253D-04, 1.16405271474738D-04, 1.09798298372713D-04, &
      1.03772410422993D-04, 9.82626078369363D-05, 9.32120517249503D-05, &
      8.85710852478712D-05, 8.42963105715700D-05, 8.03497548407791D-05/

  data alfa1(1,2), alfa1(2,2), alfa1(3,2), alfa1(4,2), alfa1(5,2), &
          alfa1(6,2), alfa1(7,2), alfa1(8,2), alfa1(9,2), alfa1(10,2), &
          alfa1(11,2),alfa1(12,2),alfa1(13,2),alfa1(14,2),alfa1(15,2), &
          alfa1(16,2),alfa1(17,2),alfa1(18,2),alfa1(19,2),alfa1(20,2), &
          alfa1(21,2),alfa1(22,2),alfa1(23,2),alfa1(24,2),alfa1(25,2), &
          alfa1(26,2)     / 6.93735541354589D-04, 2.32241745182922D-04, &
     -1.41986273556691D-05,-1.16444931672049D-04,-1.50803558053049D-04,&
     -1.55121924918096D-04,-1.46809756646466D-04,-1.33815503867491D-04, &
     -1.19744975684254D-04,-1.06184319207974D-04,-9.37699549891194D-05, &
     -8.26923045588193D-05,-7.29374348155221D-05,-6.44042357721016D-05, &
     -5.69611566009369D-05,-5.04731044303562D-05,-4.48134868008883D-05, &
     -3.98688727717599D-05,-3.55400532972042D-05,-3.17414256609022D-05, &
     -2.83996793904175D-05,-2.54522720634871D-05,-2.28459297164725D-05, &
     -2.05352753106481D-05,-1.84816217627666D-05,-1.66519330021394D-05/

  data alfa2(1,1), alfa2(2,1), alfa2(3,1), alfa2(4,1), alfa2(5,1), &
          alfa2(6,1), alfa2(7,1), alfa2(8,1), alfa2(9,1), alfa2(10,1), &
          alfa2(11,1),alfa2(12,1),alfa2(13,1),alfa2(14,1),alfa2(15,1), &
          alfa2(16,1),alfa2(17,1),alfa2(18,1),alfa2(19,1),alfa2(20,1), &
          alfa2(21,1),alfa2(22,1),alfa2(23,1),alfa2(24,1),alfa2(25,1), &
          alfa2(26,1)     /-3.54211971457744D-04,-1.56161263945159D-04, &
      3.04465503594936D-05, 1.30198655773243D-04, 1.67471106699712D-04, &
      1.70222587683593D-04, 1.56501427608595D-04, 1.36339170977445D-04, &
      1.14886692029825D-04, 9.45869093034688D-05, 7.64498419250898D-05, &
      6.07570334965197D-05, 4.74394299290509D-05, 3.62757512005344D-05, &
      2.69939714979225D-05, 1.93210938247939D-05, 1.30056674793963D-05, &
      7.82620866744497D-06, 3.59257485819352D-06, 1.44040049814252D-07, &
     -2.65396769697939D-06,-4.91346867098486D-06,-6.72739296091248D-06, &
     -8.17269379678658D-06,-9.31304715093561D-06,-1.02011418798016D-05/

  data alfa2(1,2), alfa2(2,2), alfa2(3,2), alfa2(4,2), alfa2(5,2), &
          alfa2(6,2), alfa2(7,2), alfa2(8,2), alfa2(9,2), alfa2(10,2), &
          alfa2(11,2),alfa2(12,2),alfa2(13,2),alfa2(14,2),alfa2(15,2), &
          alfa2(16,2),alfa2(17,2),alfa2(18,2),alfa2(19,2),alfa2(20,2), &
          alfa2(21,2),alfa2(22,2),alfa2(23,2),alfa2(24,2),alfa2(25,2), &
          alfa2(26,2)     / 3.78194199201773D-04, 2.02471952761816D-04, &
     -6.37938506318862D-05,-2.38598230603006D-04,-3.10916256027362D-04, &
     -3.13680115247576D-04,-2.78950273791323D-04,-2.28564082619141D-04, &
     -1.75245280340847D-04,-1.25544063060690D-04,-8.22982872820208D-05, &
     -4.62860730588116D-05,-1.72334302366962D-05, 5.60690482304602D-06, &
      2.31395443148287D-05, 3.62642745856794D-05, 4.58006124490189D-05, &
      5.24595294959114D-05, 5.68396208545815D-05, 5.94349820393104D-05, &
      6.06478527578422D-05, 6.08023907788436D-05, 6.01577894539460D-05, &
      5.89199657344698D-05, 5.72515823777593D-05, 5.52804375585853D-05/

  data beta1(1,1), beta1(2,1), beta1(3,1), beta1(4,1), beta1(5,1), &
          beta1(6,1), beta1(7,1), beta1(8,1), beta1(9,1), beta1(10,1), &
          beta1(11,1),beta1(12,1),beta1(13,1),beta1(14,1),beta1(15,1), &
          beta1(16,1),beta1(17,1),beta1(18,1),beta1(19,1),beta1(20,1), &
          beta1(21,1),beta1(22,1),beta1(23,1),beta1(24,1),beta1(25,1), &
          beta1(26,1)     / 1.79988721413553D-02, 5.59964911064388D-03, &
      2.88501402231133D-03, 1.80096606761054D-03, 1.24753110589199D-03, &
      9.22878876572938D-04, 7.14430421727287D-04, 5.71787281789705D-04, &
      4.69431007606482D-04, 3.93232835462917D-04, 3.34818889318298D-04, &
      2.88952148495752D-04, 2.52211615549573D-04, 2.22280580798883D-04, &
      1.97541838033063D-04, 1.76836855019718D-04, 1.59316899661821D-04, &
      1.44347930197334D-04, 1.31448068119965D-04, 1.20245444949303D-04, &
      1.10449144504599D-04, 1.01828770740567D-04, 9.41998224204238D-05, &
      8.74130545753834D-05, 8.13466262162801D-05, 7.59002269646219D-05/

  data beta1(1,2), beta1(2,2), beta1(3,2), beta1(4,2), beta1(5,2), &
          beta1(6,2), beta1(7,2), beta1(8,2), beta1(9,2), beta1(10,2), &
          beta1(11,2),beta1(12,2),beta1(13,2),beta1(14,2),beta1(15,2), &
          beta1(16,2),beta1(17,2),beta1(18,2),beta1(19,2),beta1(20,2), &
          beta1(21,2),beta1(22,2),beta1(23,2),beta1(24,2),beta1(25,2), &
          beta1(26,2)     /-1.49282953213429D-03,-8.78204709546389D-04, &
     -5.02916549572035D-04,-2.94822138512746D-04,-1.75463996970783D-04, &
     -1.04008550460816D-04,-5.96141953046458D-05,-3.12038929076098D-05, &
     -1.26089735980230D-05,-2.42892608575730D-07, 8.05996165414274D-06, &
      1.36507009262147D-05, 1.73964125472926D-05, 1.98672978842134D-05, &
      2.14463263790823D-05, 2.23954659232457D-05, 2.28967783814713D-05, &
      2.30785389811178D-05, 2.30321976080909D-05, 2.28236073720349D-05, &
      2.25005881105292D-05, 2.20981015361991D-05, 2.16418427448104D-05, &
      2.11507649256221D-05, 2.06388749782171D-05, 2.01165241997082D-05/

  data beta2(1,1), beta2(2,1), beta2(3,1), beta2(4,1), beta2(5,1), &
          beta2(6,1), beta2(7,1), beta2(8,1), beta2(9,1), beta2(10,1), &
          beta2(11,1),beta2(12,1),beta2(13,1),beta2(14,1),beta2(15,1), &
          beta2(16,1),beta2(17,1),beta2(18,1),beta2(19,1),beta2(20,1), &
          beta2(21,1),beta2(22,1),beta2(23,1),beta2(24,1),beta2(25,1), &
          beta2(26,1)     / 5.52213076721293D-04, 4.47932581552385D-04, &
      2.79520653992021D-04, 1.52468156198447D-04, 6.93271105657044D-05, &
      1.76258683069991D-05,-1.35744996343269D-05,-3.17972413350427D-05, &
     -4.18861861696693D-05,-4.69004889379141D-05,-4.87665447413787D-05, &
     -4.87010031186735D-05,-4.74755620890087D-05,-4.55813058138628D-05, &
     -4.33309644511266D-05,-4.09230193157750D-05,-3.84822638603221D-05, &
     -3.60857167535411D-05,-3.37793306123367D-05,-3.15888560772110D-05, &
     -2.95269561750807D-05,-2.75978914828336D-05,-2.58006174666884D-05, &
     -2.41308356761280D-05,-2.25823509518346D-05,-2.11479656768913D-05/

  data beta2(1,2), beta2(2,2), beta2(3,2), beta2(4,2), beta2(5,2), &
          beta2(6,2), beta2(7,2), beta2(8,2), beta2(9,2), beta2(10,2), &
          beta2(11,2),beta2(12,2),beta2(13,2),beta2(14,2),beta2(15,2), &
          beta2(16,2),beta2(17,2),beta2(18,2),beta2(19,2),beta2(20,2), &
          beta2(21,2),beta2(22,2),beta2(23,2),beta2(24,2),beta2(25,2), &
          beta2(26,2)     /-4.74617796559960D-04,-4.77864567147321D-04, &
     -3.20390228067038D-04,-1.61105016119962D-04,-4.25778101285435D-05, &
      3.44571294294968D-05, 7.97092684075675D-05, 1.03138236708272D-04, &
      1.12466775262204D-04, 1.13103642108481D-04, 1.08651634848774D-04, &
      1.01437951597662D-04, 9.29298396593364D-05, 8.40293133016090D-05, &
      7.52727991349134D-05, 6.69632521975731D-05, 5.92564547323195D-05, &
      5.22169308826976D-05, 4.58539485165361D-05, 4.01445513891487D-05, &
      3.50481730031328D-05, 3.05157995034347D-05, 2.64956119950516D-05, &
      2.29363633690998D-05, 1.97893056664022D-05, 1.70091984636413D-05/

  data beta3(1,1), beta3(2,1), beta3(3,1), beta3(4,1), beta3(5,1), &
          beta3(6,1), beta3(7,1), beta3(8,1), beta3(9,1), beta3(10,1), &
          beta3(11,1),beta3(12,1),beta3(13,1),beta3(14,1),beta3(15,1), &
          beta3(16,1),beta3(17,1),beta3(18,1),beta3(19,1),beta3(20,1), &
          beta3(21,1),beta3(22,1),beta3(23,1),beta3(24,1),beta3(25,1), &
          beta3(26,1)     / 7.36465810572578D-04, 8.72790805146194D-04, &
      6.22614862573135D-04, 2.85998154194304D-04, 3.84737672879366D-06, &
     -1.87906003636972D-04,-2.97603646594555D-04,-3.45998126832656D-04, &
     -3.53382470916038D-04,-3.35715635775049D-04,-3.04321124789040D-04, &
     -2.66722723047613D-04,-2.27654214122820D-04,-1.89922611854562D-04, &
     -1.55058918599094D-04,-1.23778240761874D-04,-9.62926147717644D-05, &
     -7.25178327714425D-05,-5.22070028895634D-05,-3.50347750511901D-05, &
     -2.06489761035552D-05,-8.70106096849767D-06, 1.13698686675100D-06, &
      9.16426474122779D-06, 1.56477785428873D-05, 2.08223629482467D-05/

  data gama(1),   gama(2),   gama(3),   gama(4),   gama(5), &
          gama(6),   gama(7),   gama(8),   gama(9),   gama(10), &
          gama(11),  gama(12),  gama(13),  gama(14),  gama(15), &
          gama(16),  gama(17),  gama(18),  gama(19),  gama(20), &
          gama(21),  gama(22),  gama(23),  gama(24),  gama(25), &
          gama(26)        / 6.29960524947437D-01, 2.51984209978975D-01, &
      1.54790300415656D-01, 1.10713062416159D-01, 8.57309395527395D-02, &
      6.97161316958684D-02, 5.86085671893714D-02, 5.04698873536311D-02, &
      4.42600580689155D-02, 3.93720661543510D-02, 3.54283195924455D-02, &
      3.21818857502098D-02, 2.94646240791158D-02, 2.71581677112934D-02, &
      2.51768272973862D-02, 2.34570755306079D-02, 2.19508390134907D-02, &
      2.06210828235646D-02, 1.94388240897881D-02, 1.83810633800683D-02, &
      1.74293213231963D-02, 1.65685837786612D-02, 1.57865285987918D-02, &
      1.50729501494096D-02, 1.44193250839955D-02, 1.38184805735342D-02/
!
!  I1MACH(14) replaces I1MACH(11) in a double precision code
!  I1MACH(15) replaces I1MACH(12) in a double precision code
!
  ta = epsilon ( ta )
  tol = max ( ta, 1.0D-15 )
  tb = d1mach(5)
  ju = i1mach(15)

  if ( flgjy /= 1.0D+00 ) then
    jr = i1mach(14)
    elim = 2.303D+00 * tb * ( real ( - ju ) - real ( jr ) )
  else
    elim = 2.303D+00 * ( tb * real ( - ju ) - 3.0D+00 )
  end if

  fn = fnu
  iflw = 0

  do jn = 1, in

    xx = x / fn
    wk(1) = 1.0D+00 - xx * xx
    abw2 = abs ( wk(1) )
    wk(2) = sqrt ( abw2 )
    wk(7) = fn**con2

    if ( 0.27750D+00 < abw2 ) then
      go to 80
    end if
!
!  Asymptotic expansion.
!
!  Cases near x=fn, abs ( 1-(x/fn)**2 ) <= 0.2775
!  coefficients of asymptotic expansion by series
!
!  ZETA and truncation for a(zeta) and b(zeta) series
!
!  KMAX is the truncation index for a(zeta) and b(zeta) series = max ( 2, sa )
!
    if ( abw2 == 0.0D+00 ) then
      sa = 0.0D+00
    else
      sa = tols / log ( abw2 )
    end if

    sb = sa

    do i = 1, 5
      akm = max ( sa, 2.0D+00 )
      kmax(i) = int ( akm )
      sa = sa + sb
    end do

    kb = kmax(5)
    klast = kb - 1
    sa = gama(kb)

    do k = 1, klast
      kb = kb - 1
      sa = sa * wk(1) + gama(kb)
    end do

    z = wk(1) * sa
    az = abs ( z )
    rtz = sqrt ( az )
    wk(3) = con1 * az * rtz
    wk(4) = wk(3) * fn
    wk(5) = rtz * wk(7)
    wk(6) = - wk(5) * wk(5)

    if ( 0.0D+00 < z ) then

      if ( elim < wk(4) ) then
        iflw = 1
        return
      end if

      wk(6) = -wk(6)

    end if

    phi = sqrt ( sqrt ( sa + sa + sa + sa ) )
!
!  b(zeta) for s=0
!
    kb = kmax(5)
    klast = kb - 1
    sb = beta(kb,1)

    do k = 1, klast
      kb = kb - 1
      sb = sb * wk(1) + beta(kb,1)
    end do

    ksp1 = 1
    fn2 = fn * fn
    rfn2 = 1.0D+00 / fn2
    rden = 1.0D+00
    asum = 1.0D+00
    relb = tol * abs ( sb )
    bsum = sb

    do ks = 1, 4

      ksp1 = ksp1 + 1
      rden = rden * rfn2
!
!  a(zeta) and b(zeta) for s=1,2,3,4
!
      kstemp = 5 - ks
      kb = kmax(kstemp)
      klast = kb - 1
      sa = alfa(kb,ks)
      sb = beta(kb,ksp1)

      do k = 1, klast
        kb = kb - 1
        sa = sa * wk(1) + alfa(kb,ks)
        sb = sb * wk(1) + beta(kb,ksp1)
      end do

      ta = sa * rden
      tb = sb * rden
      asum = asum + ta
      bsum = bsum + tb

      if ( abs ( ta ) <= tol .and. abs ( tb ) <= relb ) then
        exit
      end if

    end do

    bsum = bsum / ( fn * wk(7) )

    go to 160

80  continue

    upol(1) = 1.0D+00
    tau = 1.0D+00 / wk(2)
    t2 = 1.0D+00 / wk(1)
!
!  Cases for sqrt ( 1.2775 ) < (x/fn).
!
    if ( wk(1) < 0.0D+00 ) then

      wk(3) = abs ( wk(2) - atan ( wk(2) ) )
      wk(4) = wk(3) * fn
      rcz = -con1 / wk(4)
      z32 = 1.5D+00 * wk(3)
      rtz = z32**con2
      wk(5) = rtz * wk(7)
      wk(6) = -wk(5) * wk(5)
!
!  Cases for (x/fn) < sqrt ( 0.7225 )
!
    else

      wk(3) = abs ( log ( ( 1.0D+00 + wk(2) ) / xx ) - wk(2) )
      wk(4) = wk(3) * fn
      rcz = con1 / wk(4)

      if ( elim < wk(4) ) then
        iflw = 1
        return
      end if

      z32 = 1.5D+00 * wk(3)
      rtz = z32**con2
      wk(7) = fn**con2
      wk(5) = rtz * wk(7)
      wk(6) = wk(5) * wk(5)

    end if

    phi = sqrt ( ( rtz + rtz ) * tau )
    tb = 1.0D+00
    asum = 1.0D+00
    tfn = tau / fn
    upol(2) = ( c(1) * t2 + c(2) ) * tfn
    crz32 = con548 * rcz
    bsum = upol(2) + crz32
    relb = tol * abs ( bsum )
    ap = tfn
    ks = 0
    kp1 = 2
    rzden = rcz
    l = 2

    do lr = 2, 8, 2
!
!  Compute two U polynomials for next a(zeta) and b(zeta)
!
      lrp1 = lr + 1

      do k = lr, lrp1
        ks = ks + 1
        kp1 = kp1 + 1
        l = l + 1
        s1 = c(l)
        do j = 2, kp1
          l = l + 1
          s1 = s1 * t2 + c(l)
        end do
        ap = ap * tfn
        upol(kp1) = ap * s1
        cr(ks) = br(ks) * rzden
        rzden = rzden * rcz
        dr(ks) = ar(ks) * rzden
      end do

      suma = upol(lrp1)
      sumb = upol(lr+2) + upol(lrp1) * crz32
      ju = lrp1

      do jr = 1, lr
        ju = ju - 1
        suma = suma + cr(jr) * upol(ju)
        sumb = sumb + dr(jr) * upol(ju)
      end do

      tb = -tb
      if ( 0.0D+00 < wk(1) ) then
        tb = abs ( tb )
      end if

      asum = asum + suma * tb
      bsum = bsum + sumb * tb

      if ( abs ( suma ) <= tol .and. abs ( sumb ) <= relb ) then
        exit
      end if

    end do

    tb = wk(5)
    if ( 0.0D+00 < wk(1) ) then
      tb = -tb
    end if

    bsum = bsum / tb

160 continue

    call funjy ( wk(6), wk(5), wk(4), fi, dfi )
    y(jn) = flgjy * phi * ( fi * asum + dfi * bsum ) / wk(7)
    fn = fn - flgjy

  end do

  return
end
subroutine bakslv ( nr, n, a, x, b )

!*****************************************************************************80
!
!! BAKSLV solves A'*x=b where A is a lower triangular matrix.
!
!  Discussion:
!
!    BAKSLV solves the linear equations A'*X=B, where A is a
!    lower triangular matrix and A' is the transpose of A.
!
!    This routine is required by the UNCMIN minimization program.
!
!    If B is no longer required by calling routine, then vectors B and
!    X may share the same storage, and the output value of X will
!    overwrite the input value of B.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in A.
!
!    Input, real ( kind = 8 ) A(NR,N), the N by N matrix, containing the lower
!    triangular matrix.  A is not altered by this routine.
!
!    Output, real ( kind = 8 ) X(N), the solution vector.
!
!    Input, real ( kind = 8 ) B(N), the right hand side vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip1
  real ( kind = 8 ) x(n)
!
!  Solve L' * x = b.
!
  i = n
  x(i) = b(i) / a(i,i)

  if ( n == 1 ) then
    return
  end if

  do

    ip1 = i
    i = i - 1

    x(i) = ( b(i) - dot_product ( x(ip1:n), a(ip1:n,i) ) ) / a(i,i)

    if ( i == 1 ) then
      exit
    end if

  end do

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
!      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)**(N-I) * X**I
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
!    B(1,1)(X) =                X
!
!    B(2,0)(X) =     (1-X)**2
!    B(2,1)(X) = 2 * (1-X)    * X
!    B(2,2)(X) =                X**2
!
!    B(3,0)(X) =     (1-X)**3
!    B(3,1)(X) = 3 * (1-X)**2 * X
!    B(3,2)(X) = 3 * (1-X)    * X**2
!    B(3,3)(X) =                X**3
!
!    B(4,0)(X) =     (1-X)**4
!    B(4,1)(X) = 4 * (1-X)**3 * X
!    B(4,2)(X) = 6 * (1-X)**2 * X**2
!    B(4,3)(X) = 4 * (1-X)    * X**3
!    B(4,4)(X) =                X**4
!
!  Special values:
!
!    B(N,I)(X) has a unique maximum value at X = I/N.
!
!    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
!
!    B(N,I)(1/2) = C(N,K) / 2**N
!
!    For a fixed X and N, the polynomials add up to 1:
!
!      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
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
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
function besi0 ( x )

!*****************************************************************************80
!
!! BESI0 computes the hyperbolic Bessel function of the first kind, order zero.
!
!  Modified:
!
!    25 August 2001
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Bessel function.
!
!    Output, real ( kind = 8 ) BESI0, the value of the Bessel function at X.
!
  implicit none

  real ( kind = 8 ) besi0
  real ( kind = 8 ) besi0e
  real ( kind = 8 ), parameter, dimension ( 12 ) :: bi0cs = (/ &
  -0.07660547252839144951D+00,  1.927337953993808270D+00, &
   0.2282644586920301339D+00,   0.01304891466707290428D+00, &
   0.00043442709008164874D+00,  0.00000942265768600193D+00, &
   0.00000014340062895106D+00,  0.00000000161384906966D+00, &
   0.00000000001396650044D+00,  0.00000000000009579451D+00, &
   0.00000000000000053339D+00,  0.00000000000000000245D+00 /)
  real ( kind = 8 ) csevl
  integer ( kind = 4 ) inits
  integer ( kind = 4 ), save :: nti0 = 0
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xmax = 0.0D+00
  real ( kind = 8 ), save :: xsml = 0.0D+00
  real ( kind = 8 ) y

  if ( nti0 == 0 ) then
    nti0 = inits ( bi0cs, 12, 0.1D+00 * epsilon ( bi0cs ) )
    xsml = 2.0D+00 * sqrt ( epsilon ( xsml ) )
    xmax = log ( huge ( xmax ) )
  end if

  y = abs ( x )

  if ( y <= 3.0D+00 ) then

    if ( xsml < y ) then
      besi0 = 2.75D+00 + csevl ( y * y / 4.5D+00 - 1.0D+00, bi0cs, nti0 )
    else
      besi0 = 1.0D+00
    end if

    return

  end if

  if ( xmax < y ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BESI0 - Fatal error!'
    write ( *, '(a)' ) '  |X| is so big that BESI0 will overflow.'
    stop
  end if

  besi0 = exp ( y ) * besi0e ( x )

  return
end
function besi0e ( x )

!*****************************************************************************80
!
!! BESI0E computes the scaled hyperbolic Bessel function I0(X).
!
!  Discussion:
!
!    BESI0E calculates the exponentially scaled modified hyperbolic
!    Bessel function of the first kind of order zero for real argument X.
!
!      besi0e(x) = exp ( - abs ( x ) ) * i0 ( x ).
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Bessel function.
!
!    Output, real ( kind = 8 ) BESI0E, the value of the Bessel function at X.
!
  implicit none

  real ( kind = 8 ) ai02cs(22)
  real ( kind = 8 ) ai0cs(21)
  real ( kind = 8 ) besi0e
  real ( kind = 8 ) bi0cs(12)
  real ( kind = 8 ) csevl
  integer ( kind = 4 ) inits
  integer ( kind = 4 ), save :: ntai0 = 0
  integer ( kind = 4 ), save :: ntai02 = 0
  integer ( kind = 4 ), save :: nti0 = 0
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xsml = 0.0D+00
  real ( kind = 8 ) y

  data bi0cs( 1) /  -0.07660547252839144951D+00  /
  data bi0cs( 2) /   1.927337953993808270D+00  /
  data bi0cs( 3) /   0.2282644586920301339D+00  /
  data bi0cs( 4) /   0.01304891466707290428D+00  /
  data bi0cs( 5) /    0.00043442709008164874D+00  /
  data bi0cs( 6) /    0.00000942265768600193D+00  /
  data bi0cs( 7) /    0.00000014340062895106D+00  /
  data bi0cs( 8) /    0.00000000161384906966D+00  /
  data bi0cs( 9) /    0.00000000001396650044D+00  /
  data bi0cs(10) /    0.00000000000009579451D+00  /
  data bi0cs(11) /    0.00000000000000053339D+00  /
  data bi0cs(12) /    0.00000000000000000245D+00  /
  data ai0cs( 1) /    0.07575994494023796D+00  /
  data ai0cs( 2) /    0.00759138081082334D+00  /
  data ai0cs( 3) /    0.00041531313389237D+00  /
  data ai0cs( 4) /    0.00001070076463439D+00  /
  data ai0cs( 5) /   -0.00000790117997921D+00  /
  data ai0cs( 6) /   -0.00000078261435014D+00  /
  data ai0cs( 7) /    0.00000027838499429D+00  /
  data ai0cs( 8) /    0.00000000825247260D+00  /
  data ai0cs( 9) /   -0.00000001204463945D+00  /
  data ai0cs(10) /    0.00000000155964859D+00  /
  data ai0cs(11) /    0.00000000022925563D+00  /
  data ai0cs(12) /   -0.00000000011916228D+00  /
  data ai0cs(13) /    0.00000000001757854D+00  /
  data ai0cs(14) /    0.00000000000112822D+00  /
  data ai0cs(15) /   -0.00000000000114684D+00  /
  data ai0cs(16) /    0.00000000000027155D+00  /
  data ai0cs(17) /   -0.00000000000002415D+00  /
  data ai0cs(18) /   -0.00000000000000608D+00  /
  data ai0cs(19) /    0.00000000000000314D+00  /
  data ai0cs(20) /   -0.00000000000000071D+00  /
  data ai0cs(21) /    0.00000000000000007D+00  /
  data ai02cs( 1) /    0.05449041101410882D+00  /
  data ai02cs( 2) /    0.00336911647825569D+00  /
  data ai02cs( 3) /    0.00006889758346918D+00  /
  data ai02cs( 4) /    0.00000289137052082D+00  /
  data ai02cs( 5) /    0.00000020489185893D+00  /
  data ai02cs( 6) /    0.00000002266668991D+00  /
  data ai02cs( 7) /    0.00000000339623203D+00  /
  data ai02cs( 8) /    0.00000000049406022D+00  /
  data ai02cs( 9) /    0.00000000001188914D+00  /
  data ai02cs(10) /   -0.00000000003149915D+00  /
  data ai02cs(11) /   -0.00000000001321580D+00  /
  data ai02cs(12) /   -0.00000000000179419D+00  /
  data ai02cs(13) /    0.00000000000071801D+00  /
  data ai02cs(14) /    0.00000000000038529D+00  /
  data ai02cs(15) /    0.00000000000001539D+00  /
  data ai02cs(16) /   -0.00000000000004151D+00  /
  data ai02cs(17) /   -0.00000000000000954D+00  /
  data ai02cs(18) /    0.00000000000000382D+00  /
  data ai02cs(19) /    0.00000000000000176D+00  /
  data ai02cs(20) /   -0.00000000000000034D+00  /
  data ai02cs(21) /   -0.00000000000000027D+00  /
  data ai02cs(22) /    0.00000000000000003D+00  /

  if ( nti0 == 0 ) then
    nti0 = inits ( bi0cs, 12, 0.1D+00 * epsilon ( bi0cs ) )
    ntai0 = inits ( ai0cs, 21, 0.1D+00 * epsilon ( ai0cs ) )
    ntai02 = inits ( ai02cs, 22, 0.1D+00 * epsilon ( ai02cs ) )
    xsml = 2.0D+00 * sqrt ( epsilon ( xsml ) )
  end if

  y = abs ( x )

       if ( y <= xsml ) then

    besi0e = 1.0D+00

  else if ( y <= 3.0D+00 ) then

    besi0e = exp ( -y ) * &
      ( 2.75D+00 + csevl ( y*y/4.5D+00 - 1.0D+00, bi0cs, nti0 ) )

  else if ( y <= 8.0D+00 ) then

    besi0e = ( 0.375D+00 + &
      csevl ( ( 48.0D+00 / y - 11.0D+00 ) / 5.0D+00, ai0cs, ntai0 ) ) &
      / sqrt ( y )

  else if ( 8.0D+00 < y ) then

    besi0e = ( 0.375D+00 + csevl ( 16.0D+00 / y - 1.0D+00, ai02cs, ntai02 ) ) &
      / sqrt ( y )

  end if

  return
end
subroutine besj ( x, alpha, n, y, nz )

!*****************************************************************************80
!
!! BESJ computes a sequence of J Bessel functions of increasing order.
!
!  Discussion:
!
!    BESJ computes an N member sequence of J Bessel functions
!
!      J(ALPHA+K-1) (X)
!
!    for K=1,..,N for non-negative order ALPHA and X.
!
!    A combination of the power series, the asymptotic expansion for X
!    to infinity and the uniform asymptotic expansion for NU to infinity
!    are applied over subdivisions of the (NU,X) plane.  For values of
!    (NU,X) not covered by one of these formulae, the order is
!    incremented or decremented by integer ( kind = 4 ) values into a region where
!    one of the formulas apply.
!
!    Backward recursion is applied to reduce orders by integer ( kind = 4 ) values
!    except where the entire sequence lies in the oscillatory region.
!    In this case forward recursion is stable and values from the
!    asymptotic expansion for X to infinity start the recursion when it
!    is efficient to do so.
!
!    Leading terms of the series and uniform expansion are tested for
!    underflow.  If a sequence is requested and the last member would
!    underflow, the result is set to zero and the next lower order
!    tried, until a member comes on scale or all members are set
!    to zero.
!
!    Overflow cannot occur.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Donald Amos, SL Daniel, MK Weston,
!    CDC 6600 subroutines IBESS and JBESS for Bessel functions
!    I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0
!    ACM Transactions on Mathematical Software,
!    Volume 3, pages 76-92, 1977.
!
!    Frank Olver,
!    Tables of Bessel Functions of Moderate or Large Orders,
!    NPL Mathematical Tables, Volume 6,
!    Her Majesty's Stationery Office, London, 1962.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Bessel function.
!    X must be nonnegative.
!
!    Input, real ( kind = 8 ) ALPHA, the order of the first member of
!    the sequence.  ALPHA must be at least 0.0.
!
!    Input, integer ( kind = 4 ) N, the number of members in the sequence,
!    N must be at least 1.
!
!    Output, real ( kind = 8 ) Y(N), a vector whose first N components contain
!    values for J(ALPHA+K-1)(X), K=1,...,N
!
!    Output, integer ( kind = 4 ) NZ, the number of components of Y set to zero
!    due to underflow.
!
!    NZ=0, normal return, computation completed
!
!    NZ /= 0, Y(N-NZ+1) through Y(N) were set to 0.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ak
  real ( kind = 8 ) akm
  real ( kind = 8 ) alngam
  real ( kind = 8 ) alpha
  real ( kind = 8 ) ans
  real ( kind = 8 ) ap
  real ( kind = 8 ) arg
  real ( kind = 8 ) coef
  real ( kind = 8 ) d1mach
  real ( kind = 8 ) dalpha
  real ( kind = 8 ) dfn
  real ( kind = 8 ) dtm
  real ( kind = 8 ) earg
  real ( kind = 8 ) elim1
  real ( kind = 8 ) etx
  real ( kind = 8 ) fidal
  real ( kind = 8 ) flgjy
  real ( kind = 8 ) fn
  real ( kind = 8 ) fnf
  real ( kind = 8 ) fni
  real ( kind = 8 ) fnp1
  real ( kind = 8 ) fnu
  real ( kind = 8 ), parameter, dimension ( 2 ) :: fnulim = (/ &
    100.0D+00, 60.0D+00 /)
  real ( kind = 8 ) gln
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ialp
  integer ( kind = 4 ) idalp
  integer ( kind = 4 ) iflw
  integer ( kind = 4 ) in
  integer ( kind = 4 ), parameter :: inlim = 150
  integer ( kind = 4 ) is
  external jairy
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km
  integer ( kind = 4 ) kt
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) ns
  integer ( kind = 4 ) nz
  real ( kind = 8 ), parameter :: pdf = 0.785398163397448D+00
  real ( kind = 8 ), parameter :: pidt = 1.57079632679490D+00
  real ( kind = 8 ), parameter, dimension ( 4 ) :: pp = (/ &
    8.72909153935547D+00, 2.65693932265030D-01, &
    1.24578576865586D-01, 7.70133747430388D-04 /)
  real ( kind = 8 ) rden
  real ( kind = 8 ) relb
  real ( kind = 8 ), parameter :: rttp = 7.97884560802865D-01
  real ( kind = 8 ), parameter :: rtwo = 1.34839972492648D+00
  real ( kind = 8 ) rtx
  real ( kind = 8 ) rzden
  real ( kind = 8 ) s
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) sxo2
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) t
  real ( kind = 8 ) ta
  real ( kind = 8 ) tau
  real ( kind = 8 ) tb
  real ( kind = 8 ) temp(3)
  real ( kind = 8 ) tfn
  real ( kind = 8 ) tm
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolln
  real ( kind = 8 ) trx
  real ( kind = 8 ) tx
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) wk(7)
  real ( kind = 8 ) x
  real ( kind = 8 ) xo2
  real ( kind = 8 ) xo2l
  real ( kind = 8 ) y(n)

  nz = 0
  kt = 1
!
!  I1MACH(14) replaces I1MACH(11) in a double precision code
!  I1MACH(15) replaces I1MACH(12) in a double precision code
!
  ta = epsilon ( ta )
  tol = max ( ta, 1.0D-15 )
  i1 = i1mach(14) + 1
  i2 = i1mach(15)
  tb = d1mach(5)
  elim1 = 2.303D+00 * ( real ( -i2, kind = 8 ) * tb - 3.0D+00 )
!
!  TOLLN = -ln(tol)
!
  tolln = 2.303D+00 * tb * real ( i1, kind = 8 )
  tolln = min ( tolln, 34.5388D+00 )

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BESJ - Fatal error!'
    write ( *, '(a)' ) '  N is less than 1.'
    return
  end if

  if ( n == 1 ) then
    kt = 2
  end if

  nn = n

  if ( x < 0.0D+00 ) then
    call xerror ( 'BESJ - X less than zero.', 2, 1 )
    return
  end if

  if ( x == 0.0D+00 ) then

    if ( alpha < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESJ - Fatal error!'
      write ( *, '(a)' ) '  ALPHA less than zero.'
      return
    end if

    if ( alpha == 0.0D+00 ) then

      y(1) = 1.0D+00

      if ( n == 1 ) then
        return
      end if

      i1 = 2

    else

      i1 = 1

    end if

    y(i1:n) = 0.0D+00

    return

  end if

  if ( alpha < 0.0D+00 ) then
    call xerror ( 'BESJ - order, alpha, less than zero.', 2, 1)
    return
  end if

  ialp = int ( alpha )
  fni = real ( ialp + n - 1, kind = 8 )
  fnf = alpha - real ( ialp, kind = 8 )
  dfn = fni + fnf
  fnu = dfn
  xo2 = x * 0.5D+00
  sxo2 = xo2 * xo2
!
!  Decision tree for region where series, asymptotic expansion for x
!  to infinity and asymptotic expansion for nu to infinity are applied.
!
  if ( sxo2 <= ( fnu+1.0D+00 ) ) then
    go to 90
  end if

  ta = max ( 20.0D+00, fnu )

  if ( ta < x ) then
    go to 120
  end if

  if ( 12.0D+00 < x ) then
    go to 110
  end if

  xo2l = log ( xo2 )
  ns = int ( sxo2 - fnu ) + 1
  go to 100

90 continue

  fn = fnu
  fnp1 = fn + 1.0D+00
  xo2l = log ( xo2 )
  is = kt

  if ( x <= 0.50D+00 ) then
    go to 330
  end if

  ns = 0

100 continue

  fni = fni + real ( ns, kind = 8 )
  dfn = fni + fnf
  fn = dfn
  fnp1 = fn + 1.0D+00
  is = kt

  if ( 0 < n - 1 + ns ) then
    is = 3
  end if

  go to 330

110 continue

  ans = max ( 36.0D+00 - fnu, 0.0D+00 )
  ns = int ( ans )
  fni = fni + real ( ns, kind = 8 )
  dfn = fni + fnf
  fn = dfn
  is = kt

  if ( 0 < n - 1 + ns ) then
    is = 3
  end if

  go to 130

120 continue

  rtx = sqrt ( x )
  tau = rtwo * rtx
  ta = tau + fnulim(kt)

  if ( fnu <= ta ) then
    go to 480
  end if

  fn = fnu
  is = kt
!
!  Uniform asymptotic expansion for NU to infinity.
!
130 continue

  i1 = abs ( 3 - is )
  i1 = max ( i1, 1 )
  flgjy = 1.0D+00

  call asyjy ( jairy, x, fn, flgjy, i1, temp(is), wk, iflw )

  if ( iflw /= 0 ) then
    go to 380
  end if

  go to (320, 450, 620), is

310 continue

  temp(1) = temp(3)
  kt = 1

320 continue

  is = 2
  fni = fni - 1.0D+00
  dfn = fni + fnf
  fn = dfn

  if ( i1 == 2 ) then
    go to 450
  end if

  go to 130
!
!  Series for (x/2)**2<=nu+1
!
330 continue

  gln = alngam ( fnp1 )
  arg = fn * xo2l - gln

  if ( arg < (-elim1) ) then
    go to 400
  end if

  earg = exp ( arg )

340 continue

  s = 1.0D+00

  if ( x < tol ) then
    go to 360
  end if

  ak = 3.0D+00
  t2 = 1.0D+00
  t = 1.0D+00
  s1 = fn

  do k = 1, 17
    s2 = t2 + s1
    t = - t * sxo2 / s2
    s = s + t
    if ( abs ( t ) < tol ) then
      exit
    end if
    t2 = t2 + ak
    ak = ak + 2.0D+00
    s1 = s1 + fn
  end do

360 continue

  temp(is) = s * earg
  go to (370, 450, 610), is

370 continue
  earg = earg * fn / xo2
  fni = fni - 1.0D+00
  dfn = fni + fnf
  fn = dfn
  is = 2
  go to 340
!
!  Set underflow value and update parameters
!
380 continue

  y(nn) = 0.0D+00
  nn = nn - 1
  fni = fni - 1.0D+00
  dfn = fni + fnf
  fn = dfn
  if ( nn-1 ) 440, 390, 130

390 continue

  kt = 2
  is = 2
  go to 130

400 continue

  y(nn) = 0.0D+00
  nn = nn - 1
  fnp1 = fn
  fni = fni - 1.0D+00
  dfn = fni + fnf
  fn = dfn
  if ( nn-1 ) 440, 410, 420

410 continue
  kt = 2
  is = 2

420 continue

  if ( sxo2 <= fnp1 ) then
    go to 430
  end if

  go to 130

430 continue

  arg = arg - xo2l + log ( fnp1 )
  if ( arg < (-elim1) ) then
    go to 400
  end if

  go to 330

440 nz = n - nn
  return
!
!  Backward recursion section
!
450 continue
  nz = n - nn
  if ( kt == 2 ) then
    go to 470
  end if
!
!  Backward recur from index ALPHA+NN-1 to ALPHA.
!
  y(nn) = temp(1)
  y(nn-1) = temp(2)
  if ( nn == 2 ) then
    return
  end if

  trx = 2.0D+00 / x
  dtm = fni
  tm = ( dtm + fnf ) * trx
  k = nn + 1

  do i = 3, nn
    k = k - 1
    y(k-2) = tm * y(k-1) - y(k)
    dtm = dtm - 1.0D+00
    tm = ( dtm + fnf ) * trx
  end do

  return

470 continue

  y(1) = temp(2)
  return
!
!  Asymptotic expansion for X to infinity with forward recursion in
!  oscillatory region max ( 20, NU ) < X, provided the last member
!  of the sequence is also in the region.
!
480 continue

  in = int ( alpha - tau + 2.0D+00 )

  if ( in <= 0 ) then
    go to 490
  end if

  idalp = ialp - in - 1
  kt = 1
  go to 500

490 continue

  idalp = ialp
  in = 0

500 continue

  is = kt
  fidal = real ( idalp, kind = 8 )
  dalpha = fidal + fnf
  arg = x - pidt * dalpha - pdf
  sa = sin ( arg )
  sb = cos ( arg )
  coef = rttp / rtx
  etx = 8.0D+00 * x

510 continue

  dtm = fidal + fidal
  dtm = dtm * dtm
  tm = 0.0D+00

  if ( fidal == 0.0D+00 .and. abs ( fnf ) < tol ) then
    go to 520
  end if

  tm = 4.0D+00 * fnf * ( fidal + fidal + fnf )

520 continue

  trx = dtm - 1.0D+00
  t2 = ( trx + tm ) / etx
  s2 = t2
  relb = tol * abs ( t2 )
  t1 = etx
  s1 = 1.0D+00
  fn = 1.0D+00
  ak = 8.0D+00

  do k = 1, 13
    t1 = t1 + etx
    fn = fn + ak
    trx = dtm - fn
    ap = trx + tm
    t2 = -t2 * ap / t1
    s1 = s1 + t2
    t1 = t1 + etx
    ak = ak + 8.0D+00
    fn = fn + ak
    trx = dtm - fn
    ap = trx + tm
    t2 = t2 * ap / t1
    s2 = s2 + t2
    if ( abs ( t2 ) <= relb ) then
      exit
    end if
    ak = ak + 8.0D+00
  end do

540 continue

  temp(is) = coef * ( s1 * sb - s2 * sa )

  if ( is == 2 ) then
    go to 560
  end if

550 continue

  fidal = fidal + 1.0D+00
  dalpha = fidal + fnf
  is = 2
  tb = sa
  sa = -sb
  sb = tb
  go to 510
!
!  Forward recursion section
!
560 continue

  if ( kt == 2 ) then
    go to 470
  end if

  s1 = temp(1)
  s2 = temp(2)
  tx = 2.0D+00 / x
  tm = dalpha * tx

  if ( in == 0 ) then
    go to 580
  end if
!
!   Forward recursion to index alpha
!
  do i = 1, in
    s = s2
    s2 = tm * s2 - s1
    tm = tm + tx
    s1 = s
  end do

  if ( nn == 1 ) then
    go to 600
  end if

  s = s2
  s2 = tm * s2 - s1
  tm = tm + tx
  s1 = s

580 continue
!
!  Forward recursion from index ALPHA to ALPHA+N-1.
!
  y(1) = s1
  y(2) = s2

  do i = 3, nn
    y(i) = tm * y(i-1) - y(i-2)
    tm = tm + tx
  end do

  return

600 continue

  y(1) = s2
  return
!
!  Backward recursion with normalization by
!  asymptotic expansion for nu to infinity or power series.
!
610 continue
!
!  Computation of last order for series normalization
!
  akm = max ( 3.0D+00 - fn, 0.0D+00 )
  km = int ( akm )
  tfn = fn + real ( km, kind = 8 )
  ta = ( gln + tfn - 0.9189385332D+00 - 0.0833333333D+00 / tfn ) &
    / ( tfn + 0.5D+00 )
  ta = xo2l - ta
  tb = - ( 1.0D+00 -1.5D+00 / tfn ) / tfn
  akm = tolln / ( - ta + sqrt ( ta * ta - tolln * tb ) ) + 1.5D+00
  in = km + int ( akm )
  go to 660

  620 continue
!
!  Computation of last order for asymptotic expansion normalization
!
  gln = wk(3) + wk(2)
  if ( 30.0D+00 < wk(6) ) then
    go to 640
  end if

  rden = ( pp(4) * wk(6) + pp(3) ) * wk(6) + 1.0D+00
  rzden = pp(1) + pp(2) * wk(6)
  ta = rzden / rden

  if ( wk(1) < 0.10D+00 ) then
    go to 630
  end if

  tb = gln / wk(5)
  go to 650

630 continue

  tb = ( 1.259921049D+00 + ( 0.1679894730D+00 + 0.0887944358D+00 &
    * wk(1) ) * wk(1) ) / wk(7)
  go to 650

640 continue
  ta = 0.5D+00 * tolln / wk(4)
  ta=( ( 0.0493827160D+00 * ta - 0.1111111111D+00 ) * ta &
    + 0.6666666667D+00 ) * ta * wk(6)

  if ( wk(1) < 0.10D+00 ) then
    go to 630
  end if

  tb = gln / wk(5)

650 continue

  in = int ( ta / tb + 1.5D+00 )

  if ( inlim < in ) then
    go to 310
  end if

660 continue
  dtm = fni + real ( in, kind = 8 )
  trx = 2.0D+00 / x
  tm = ( dtm + fnf ) * trx
  ta = 0.0D+00
  tb = tol
  kk = 1

670 continue
!
!  Backward recur unindexed
!
  do i = 1, in
    s = tb
    tb = tm * tb - ta
    ta = s
    dtm = dtm - 1.0D+00
    tm = ( dtm + fnf ) * trx
  end do
!
!  Normalization.
!
  if ( kk == 1 ) then

    ta = ( ta / tb ) * temp(3)
    tb = temp(3)
    kk = 2
    in = ns

    if ( ns /= 0 ) then
      go to 670
    end if

  end if

  y(nn) = tb
  nz = n - nn

  if ( nn == 1 ) then
    return
  end if

  k = nn - 1
  y(k) = tm * tb - ta

  if ( nn == 2 ) then
    return
  end if

  dtm = dtm - 1.0D+00
  tm = ( dtm + fnf ) * trx
  km = k - 1
!
!  Backward recur indexed
!
  do i = 1, km
    y(k-1) = tm * y(k) - y(k+1)
    dtm = dtm - 1.0D+00
    tm = ( dtm + fnf ) * trx
    k = k - 1
  end do

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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
subroutine bp01 ( n, x, b )

!*****************************************************************************80
!
!! BP01 evaluates the N+1 Bernstein basis functions of degree N on [0,1].
!
!  Definition:
!
!    The I-th Bernstein basis polynomial of degree N is defined as:
!
!      B(N,I,X)= N!/(I!*(N-I)!) * (1-X)**(N-I) * X**I
!
!    although this is not how the values are computed.
!
!  Modified:
!
!    08 February 2003
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, should be 0 or greater.
!
!    Input, real ( kind = 8 ) X, the point where the functions should be
!    evaluated.
!
!    Output, real ( kind = 8 ) B(0:N), the values of the Bernstein polynomials
!    at the point X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  if ( n == 0 ) then

    b(0) = 1.0D+00

  else if ( 0 < n ) then

    do i = 1, n

       if ( i == 1 ) then
         b(1) = x
       else
         b(i) = x * b(i-1)
       end if

       do j = i-1, 1, -1
         b(j) = x * b(j-1) + ( 1.0D+00 - x ) * b(j)
       end do

       if ( i == 1 ) then
         b(0) = 1.0D+00 - x
       else
         b(0) = ( 1.0D+00 - x ) * b(0)
       end if

    end do

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
!    A C8VEC is a vector of complex ( kind = 8 ) values.
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
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last entries to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title
  complex ( kind = 8 ) x(n)

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = max ( 1, i_lo ), min ( n, i_hi )
    write ( *, '(2x,i8,2x,2g14.6)' ) i, x(i)
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
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    theta = 2.0D+00 * pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

  end do

  return
end
subroutine chfdv ( x1, x2, f1, f2, d1, d2, ne, xe, fe, de, next, ierr )

!*****************************************************************************80
!
!! CHFDV evaluates a cubic polynomial and its derivative given in Hermite form.
!
!  Discussion:
!
!    CHFDV evaluates a cubic polynomial and its first derivative.
!    The cubic polynomial is given in Hermite form.  The evaluation
!    is carried out at an array of points.
!
!    This routine was designed for use by PCHFD, but it may also be
!    useful directly as an evaluator for a piecewise cubic Hermite
!    function in applications, such as graphing, where the interval
!    is known in advance.
!
!    If only function values are required, use CHFEV instead.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the interval of
!    definition of  the cubic.  X1 and X2 must be distinct.
!
!    Input, real ( kind = 8 ) F1, F2, the values of the function at X1 and
!    X2, respectively.
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at the ends
!     of the interval.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), the points at which the functions are to
!    be evaluated.  If any of the XE are outside the interval
!    [X1,X2], a warning error is returned in next.
!
!    Output, real ( kind = 8 ) FE(NE), DE(NE), the values of the cubic
!    function and its derivative at the points XE(*).
!
!    Output, integer ( kind = 4 ) NEXT(2), indicates the number of extrapolation points:
!    NEXT(1) = number of evaluation points to left of interval.
!    NEXT(2) = number of evaluation points to right of interval.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, NE < 1.
!    -2, X1 == X2.
!
  implicit none

  integer ( kind = 4 ) ne

  real ( kind = 8 ) c2
  real ( kind = 8 ) c2t2
  real ( kind = 8 ) c3
  real ( kind = 8 ) c3t3
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) de(ne)
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) delta
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fe(ne)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) next(2)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xe(ne)
  real ( kind = 8 ) xma
  real ( kind = 8 ) xmi
!
!  Check arguments.
!
  if ( ne < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The number of evaluation points was less than 1.'
    stop
  end if

  h = x2 - x1

  if ( h == 0.0D+00 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The interval endpoints are equal.'
    return
  end if
!
!  Initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0D+00, h )
  xma = max ( 0.0D+00, h )
!
!  Compute cubic coefficients expanded about X1.
!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h

  c2 = -( del1 + del1 + del2 )
  c2t2 = c2 + c2
  c3 = ( del1 + del2 ) / h
  c3t3 = c3 + c3 + c3
!
!  Evaluation loop.
!
  do i = 1, ne

    x = xe(i) - x1
    fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
    de(i) = d1 + x * ( c2t2 + x * c3t3 )
!
!  Count extrapolation points.
!
    if ( x < xmi ) then
      next(1) = next(1) + 1
    end if

    if ( xma < x ) then
      next(2) = next(2) + 1
    end if

  end do

  return
end
subroutine chfev ( x1, x2, f1, f2, d1, d2, ne, xe, fe, next, ierr )

!*****************************************************************************80
!
!! CHFEV evaluates a cubic polynomial given in Hermite form.
!
!  Discussion:
!
!    This routine evaluates a cubic polynomial given in Hermite form at an
!    array of points.  While designed for use by PCHFE, it may
!    be useful directly as an evaluator for a piecewise cubic
!    Hermite function in applications, such as graphing, where
!    the interval is known in advance.
!
!    The cubic polynomial is determined by function values
!    F1, F2 and derivatives D1, D2 on the interval [X1,X2].
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the interval of
!    definition of the cubic.  X1 and X2 must be distinct.
!
!    Input, real ( kind = 8 ) F1, F2, the values of the function at X1 and
!    X2, respectively.
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at X1 and
!    X2, respectively.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), the points at which the function is to
!    be evaluated.  If any of the XE are outside the interval
!    [X1,X2], a warning error is returned in NEXT.
!
!    Output, real ( kind = 8 ) FE(NE), the value of the cubic function
!    at the points XE.
!
!    Output, integer ( kind = 4 ) NEXT(2), the number of extrapolation points:
!    NEXT(1) = number of evaluation points to the left of interval.
!    NEXT(2) = number of evaluation points to the right of interval.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, NE < 1.
!    -2, X1 == X2.
!
  implicit none

  integer ( kind = 4 ) ne

  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) delta
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fe(ne)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) next(2)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xe(ne)
  real ( kind = 8 ) xma
  real ( kind = 8 ) xmi

  if ( ne < 1 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFEV - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points is less than 1.'
    write ( *, '(a,i6)' ) '  NE = ', ne
    stop
  end if

  h = x2 - x1

  if ( h == 0.0D+00 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFEV - Fatal error!'
    write ( *, '(a)' ) '  The interval [X1,X2] is of zero length.'
    stop
  end if
!
!  Initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0D+00, h )
  xma = max ( 0.0D+00, h )
!
!  Compute cubic coefficients expanded about X1.
!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h
  c2 = -( del1 + del1 + del2 )
  c3 = ( del1 + del2 ) / h
!
!  Evaluation loop.
!
  do i = 1, ne

    x = xe(i) - x1
    fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
!
!  Count the extrapolation points.
!
    if ( x < xmi ) then
      next(1) = next(1) + 1
    end if

    if ( xma < x ) then
      next(2) = next(2) + 1
    end if

  end do

  return
end
function chfiv ( x1, x2, f1, f2, d1, d2, a, b, ierr )

!*****************************************************************************80
!
!! CHFIV evaluates the integral of a cubic polynomial in Hermite form.
!
!  Discussion:
!
!    CHFIV is called by PCHIA to evaluate the integral of a single cubic (in
!    Hermite form) over an arbitrary interval (A,B).
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VALUE, the value of the requested integral.
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the interval of
!    definition of the cubic.  X1 and X2 must be distinct.
!
!    Input, real ( kind = 8 ) F1, F2, the values of the function at X1
!    and X2, respectively.
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at the ends
!    of the interval.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of interval of integration.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, X1 == X2.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) chfiv
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dterm
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fterm
  real ( kind = 8 ) h
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) phia1
  real ( kind = 8 ) phia2
  real ( kind = 8 ) phib1
  real ( kind = 8 ) phib2
  real ( kind = 8 ) psia1
  real ( kind = 8 ) psia2
  real ( kind = 8 ) psib1
  real ( kind = 8 ) psib2
  real ( kind = 8 ) ta1
  real ( kind = 8 ) ta2
  real ( kind = 8 ) tb1
  real ( kind = 8 ) tb2
  real ( kind = 8 ) ua1
  real ( kind = 8 ) ua2
  real ( kind = 8 ) ub1
  real ( kind = 8 ) ub2
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Check input.
!
  if ( x1 == x2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFIV - Fatal error!'
    write ( *, '(a)' ) '  X1 = X2.'
    stop
  end if

  ierr = 0
!
!  Compute integral.
!
  h = x2 - x1
  ta1 = ( a - x1 ) / h
  ta2 = ( x2 - a ) / h
  tb1 = ( b - x1 ) / h
  tb2 = ( x2 - b ) / h

  ua1 = ta1 * ta1 * ta1
  phia1 = ua1 * ( 2.0D+00 - ta1 )
  psia1 = ua1 * ( 3.0D+00 * ta1 - 4.0D+00 )
  ua2 = ta2 * ta2 * ta2
  phia2 =  ua2 * ( 2.0D+00 - ta2)
  psia2 = -ua2 * ( 3.0D+00 * ta2 - 4.0D+00 )

  ub1 = tb1 * tb1 * tb1
  phib1 = ub1 * ( 2.0D+00 - tb1 )
  psib1 = ub1 * ( 3.0D+00 * tb1 - 4.0D+00 )
  ub2 = tb2 * tb2 * tb2
  phib2 =  ub2 * ( 2.0D+00 - tb2 )
  psib2 = -ub2 * ( 3.0D+00 * tb2 - 4.0D+00 )

  fterm =   f1 * ( phia2 - phib2 ) + f2 * ( phib1 - phia1 )
  dterm = ( d1 * ( psia2 - psib2 ) + d2 * ( psib1 - psia1 ) ) * ( h / 6.0D+00 )

  chfiv = 0.5D+00 * h * ( fterm + dterm )

  return
end
function chfmc ( d1, d2, delta )

!*****************************************************************************80
!
!! CHFMC determines the monotonicity properties of a cubic polynomial.
!
!  Discussion:
!
!    CHFMC is called by PCHMC to determine the monotonicity properties
!    of the cubic with boundary derivative values D1, D2 and chord
!    slope DELTA.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at the ends
!    of the interval.
!
!    Input, real ( kind = 8 ) DELTA, the data slope over that interval.
!
!    Output, integer ( kind = 4 ) CHFMC, indicates the monotonicity of the 
!    cubic segment:
!    -1, if function is strictly decreasing;
!     0, if function is constant;
!     1, if function is strictly increasing;
!     2, if function is non-monotonic;
!     3, if unable to determine.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) chfmc
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) delta
  real ( kind = 8 ) eps
  integer ( kind = 4 ) ismon
  integer ( kind = 4 ) itrue
  real ( kind = 8 ) phi

  eps = 10.0D+00 * epsilon ( eps )
!
!  Make the check.
!
  if ( delta == 0.0D+00 ) then

    if ( d1 == 0.0D+00 .and. d2 == 0.0D+00 ) then
      ismon = 0
    else
      ismon = 2
    end if

  else

     itrue = sign ( 1.0D+00, delta)
     a = d1 / delta
     b = d2 / delta

     if ( a < 0.0D+00 .or. b < 0.0D+00 ) then
       ismon = 2
     else if ( a <= 3.0D+00 - eps  .and. b <= 3.0D+00 -eps ) then
!
!  Inside square (0,3)x(0,3) implies OK.
!
       ismon = itrue
     else if ( 4.0D+00 + eps < a .and. 4.0D+00 + eps < b ) then
!
!  Outside square (0,4)x(0,4) implies nonmonotonic.
!
       ismon = 2
     else
!
!  Must check against boundary of ellipse.
!
        a = a - 2.0D+00
        b = b - 2.0D+00
        phi = ( ( a * a + b * b ) + a * b ) - 3.0D+00

        if ( phi < -eps ) then
          ismon = itrue
        else if ( eps < phi ) then
          ismon = 2
        else
!
!  Too close to boundary to tell,
!  in the presence of round-off errors.
!
           ismon = 3
        end if
    end if
  end if

  chfmc = ismon

  return
end
subroutine chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )

!*****************************************************************************80
!
!! CHKDER checks the gradients of M functions of N variables.
!
!  Discussion:
!
!    CHKDER checks the gradients of M nonlinear functions in N variables,
!    evaluated at a point X, for consistency with the functions themselves.
!
!    The user calls CHKDER twice, first with MODE = 1 and then with MODE = 2.
!
!    MODE = 1.
!      On input,
!        X contains the point of evaluation.
!      On output,
!        XP is set to a neighboring point.
!
!    Now the user must evaluate the function and gradients at X, and the
!    function at XP.  Then the subroutine is called again:
!
!    MODE = 2.
!      On input,
!        FVEC contains the function values at X,
!        FJAC contains the function gradients at X.
!        FVECP contains the functions evaluated at XP.
!      On output,
!        ERR contains measures of correctness of the respective gradients.
!
!    The subroutine does not perform reliably if cancellation or
!    rounding errors cause a severe loss of significance in the
!    evaluation of a function.  Therefore, none of the components
!    of X should be unusually small (in particular, zero) or any
!    other value which may cause loss of significance.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, is the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian is
!    to be evaluated.
!
!    Input, real ( kind = 8 ) FVEC(M), is used only when MODE = 2.
!    In that case, it should contain the function values at X.
!
!    Input, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  When MODE = 2,
!    FJAC(I,J) should contain the value of dF(I)/dX(J).
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of the array FJAC.
!    LDFJAC must be at least M.
!
!    Output, real ( kind = 8 ) XP(N), on output with MODE = 1, is a
!    neighboring point of X, at which the function is to be evaluated.
!
!    Input, real ( kind = 8 ) FVECP(M), on input with MODE = 2, is the
!    function value at XP.
!
!    Input, integer ( kind = 4 ) MODE, should be set to 1 on the first call, and
!    2 on the second.
!
!    Output, real ( kind = 8 ) ERR(M).  On output when MODE = 2, ERR
!    contains measures of correctness of the respective gradients.  If
!    there is no severe loss of significance, then if ERR(I):
!      = 1.0D+00, the I-th gradient is correct,
!      = 0.0D+00, the I-th gradient is incorrect.
!      > 0.5D+00, the I-th gradient is probably correct.
!      < 0.5D+00, the I-th gradient is probably incorrect.
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) eps
  real ( kind = 8 ) epsf
  real ( kind = 8 ) epslog
  real ( kind = 8 ) epsmch
  real ( kind = 8 ) err(m)
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(m)
  real ( kind = 8 ) fvecp(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mode
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xp(n)

  epsmch = epsilon ( epsmch )
  eps = sqrt ( epsmch )
!
!  MODE = 1.
!
  if ( mode == 1 ) then

     do j = 1, n
       temp = eps * abs ( x(j) )
       if ( temp == 0.0D+00 ) then
         temp = eps
       end if
       xp(j) = x(j) + temp
     end do
!
!  MODE = 2.
!
  else if ( mode == 2 ) then

     epsf = 100.0D+00 * epsmch
     epslog = log10 ( eps )

     err = 0.0D+00

     do j = 1, n
       temp = abs ( x(j) )
       if ( temp == 0.0D+00 ) then
         temp = 1.0D+00
       end if
       err(1:m) = err(1:m) + temp * fjac(1:m,j)
     end do

     do i = 1, m

       temp = 1.0D+00

       if ( fvec(i) /= 0.0D+00 .and. fvecp(i) /= 0.0D+00 .and. &
         epsf * abs ( fvec(i) ) <= abs ( fvecp(i) - fvec(i) ) ) then
         temp = eps * abs ( ( fvecp(i) - fvec(i) ) / eps - err(i) ) &
           / ( abs ( fvec(i) ) + abs ( fvecp(i) ) )
       end if

       err(i) = 1.0D+00

       if ( epsmch < temp .and. temp < eps ) then
         err(i) = ( log10 ( temp ) - epslog ) / epslog
       end if

       if ( eps <= temp ) then
         err(i) = 0.0D+00
       end if

     end do

  end if

  return
end
subroutine chlhsn ( nr, n, a, epsm, sx, udiag )

!*****************************************************************************80
!
!! CHLHSN finds the L*L' decomposition of the perturbed model hessian matrix.
!
!  Discussion:
!
!    The perturbed model Hessian matrix has the form
!
!      A + MU * I
!
!    (where 0 <= MU and I is the identity matrix) which is safely
!    positive definite.
!
!    If A is safely positive definite upon entry, then MU=0.
!
!    1. If A has any negative diagonal elements, then choose 0 < MU
!    such that the diagonal of A:=A+MU*I is all positive
!    with the ratio of its smallest to largest element on the
!    order of sqrt ( EPSM ).
!
!    2. A undergoes a perturbed Cholesky decomposition which
!    results in an LL+ decomposition of A+D, where D is a
!    non-negative diagonal matrix which is implicitly added to
!    A during the decomposition if A is not positive definite.
!    A is retained and not changed during this process by
!    copying L into the upper triangular part of A and the
!    diagonal into UDIAG.  Then the Cholesky decomposition routine
!    is called.  On return, ADDMAX contains the maximum element of D.
!
!    3. If ADDMAX=0, A was positive definite going into step 2
!    and return is made to calling program.  Otherwise,
!    the minimum number SDD which must be added to the
!    diagonal of A to make it safely strictly diagonally dominant
!    is calculated.  Since A + ADDMAX * I and A + SDD * I are safely
!    positive definite, choose MU = min ( ADDMAX, SDD ) and decompose
!    A + MU * I to obtain L.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input/output, real A(NR,N), contains an N by N matrix.
!    On input, A is the model hessian.  Only the lower triangular part and
!    diagonal are stored.  On output, A contains the factor L of the
!    LL+ decomposition of the perturbed model hessian in the lower triangular
!    part and diagonal, and contains the hessian in the upper triangular part
!    and UDIAG.
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Output, real ( kind = 8 ) UDIAG(N), the diagonal of the hessian.
!
!  Local variables:
!
!    tol              tolerance
!    diagmn           minimum element on diagonal of a
!    diagmx           maximum element on diagonal of a
!    offmax           maximum off-diagonal element of a
!    offrow           sum of off-diagonal elements in a row of a
!    evmin            minimum eigenvalue of a
!    evmax            maximum eigenvalue of a
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) addmax
  real ( kind = 8 ) amu
  real ( kind = 8 ) diagmx
  real ( kind = 8 ) diagmn
  real ( kind = 8 ) epsm
  real ( kind = 8 ) evmax
  real ( kind = 8 ) evmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) offmax
  real ( kind = 8 ) offrow
  real ( kind = 8 ) posmax
  real ( kind = 8 ) sdd
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) tol
  real ( kind = 8 ) udiag(n)
!
!  Scale the hessian.
!
  do j = 1, n
    do i = j, n
      a(i,j) = a(i,j) / ( sx(i) * sx(j) )
    end do
  end do
!
!  Step1
!
  tol = sqrt ( epsm )

  diagmx = a(1,1)
  diagmn = a(1,1)

  do i = 2, n
    if ( a(i,i) < diagmn ) then
      diagmn = a(i,i)
    end if
    if ( diagmx < a(i,i) ) then
      diagmx = a(i,i)
    end if
  end do

  posmax = max ( diagmx, 0.0D+00 )

  if ( diagmn <= posmax * tol ) then

    amu = tol * ( posmax - diagmn ) - diagmn
!
!  Find the largest off-diagonal element of A.
!
    if ( amu == 0.0D+00 ) then

      offmax = 0.0D+00

      do i = 2, n
        do j = 1, i-1
          if ( offmax < abs ( a(i,j) ) ) then
            offmax = abs ( a(i,j) )
          end if
        end do
      end do

      amu = offmax

      if ( amu == 0.0D+00 ) then
        amu = 1.0D+00
      else
        amu = amu * ( 1.0D+00 + tol )
      end if

    end if
!
!  A = A + MU*I
!
    do i = 1, n
      a(i,i) = a(i,i) + amu
    end do

    diagmx = diagmx + amu

  end if
!
!  Step2
!
!  Copy lower triangular part of A to upper triangular part
!  and diagonal of A to udiag
!
  do j = 1, n
    udiag(j) = a(j,j)
    do i = j + 1, n
      a(j,i) = a(i,j)
    end do
  end do

  call choldc ( nr, n, a, diagmx, tol, addmax )
!
!  Step3
!
!  If ADDMAX=0, A was positive definite going into step 2,
!  the ll+ decomposition has been done, and we return.
!
!  Otherwise, 0 < ADDMAX.  perturb A so that it is safely
!  diagonally dominant and find ll+ decomposition
!
  if ( 0.0D+00 < addmax ) then
!
!  Restore original A (lower triangular part and diagonal)
!
    do j = 1, n
      a(j,j) = udiag(j)
      do i = j+1, n
        a(i,j) = a(j,i)
      end do
    end do
!
!  Find SDD such that A+sdd*i is safely positive definite
!  note:  evmin<0 since A is not positive definite;
!
    evmin = 0.0D+00
    evmax = a(1,1)

    do i = 1, n

      offrow = sum ( abs ( a(i,1:i-1) ) ) + sum ( abs ( a(i+1:n,i) ) )
      evmin = min ( evmin, a(i,i)-offrow )
      evmax = max ( evmax, a(i,i)+offrow )

    end do

    sdd = tol * ( evmax - evmin ) - evmin
!
!  Perturb A and decompose again.
!
    amu = min ( sdd, addmax )

    do i = 1, n
      a(i,i) = a(i,i) + amu
      udiag(i) = a(i,i)
    end do
!
!  A is now guaranteed safely positive definite
!
    call choldc ( nr, n, a, 0.0D+00, tol, addmax )

  end if
!
!  Unscale the hessian and Cholesky decomposition matrix.
!
  do j = 1, n

    a(j:n,j) = sx(j:n) * a(j:n,j)

    do i = 1, j - 1
      a(i,j) = sx(i) * sx(j) * a(i,j)
    end do

    udiag(j) = udiag(j) * sx(j) * sx(j)

  end do

  return
end
subroutine choldc ( nr, n, a, diagmx, tol, addmax )

!*****************************************************************************80
!
!! CHOLDC finds the perturbed L*L' decomposition of A+D.
!
!  Discussion:
!
!    D is a non-negative diagonal matrix added to A if
!    necessary to allow the Cholesky decomposition to continue.
!
!    The normal Cholesky decomposition is performed.  However, if at any
!    point the algorithm would attempt to set
!      L(I,I) = sqrt ( TEMP )
!    with
!      TEMP < TOL * DIAGMX,
!    then L(I,I) is set to sqrt ( TOL * DIAGMX )
!    instead.  This is equivalent to adding TOL * DIAGMX-TEMP to A(I,I)
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input/output, real ( kind = 8 ) A(NR,N), the N by N matrix.
!    On input, the matrix for which to find the perturbed
!    Cholesky decomposition.
!    On output, the lower triangular part contains the L factor,
!    and the diagonal of A.
!
!    Input, real ( kind = 8 ) DIAGMX, the maximum diagonal element of A.
!
!    Input, real ( kind = 8 ) TOL, a tolerance.
!
!    Output, real ( kind = 8 ) ADDMAX, the maximum amount implicitly added to
!    the diagonal of A in forming the Cholesky decomposition of A+D.
!
!  Local variables:
!
!    aminl    smallest element allowed on diagonal of L.
!
!    amnlsq   =aminl**2
!
!    offmax   maximum off-diagonal element in column of a
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) addmax
  real ( kind = 8 ) aminl
  real ( kind = 8 ) amnlsq
  real ( kind = 8 ) diagmx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) offmax
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) tol

  addmax = 0.0D+00
  aminl = sqrt ( diagmx * tol )
  amnlsq = aminl**2
!
!  Form column J of L.
!
  do j = 1, n
!
!  Find diagonal elements of L.
!
    sum2 = sum ( a(j,1:j-1)**2 )

    temp = a(j,j) - sum2

    if ( amnlsq <= temp ) then

      a(j,j) = sqrt ( temp )
!
!  Find maximum off-diagonal element in column.
!
    else

      offmax = 0.0D+00

      do i = j+1, n
        if ( offmax < abs ( a(i,j) ) ) then
          offmax = abs ( a(i,j) )
        end if
      end do

      if ( offmax <= amnlsq ) then
        offmax = amnlsq
      end if
!
!  Add to diagonal element to allow Cholesky decomposition to continue
!
      a(j,j) = sqrt ( offmax )
      addmax = max ( addmax, offmax - temp )

    end if
!
!  Find (I,J) element of lower triangular matrix.
!
    do i = j+1, n
      sum2 = 0.0D+00
      do k = 1, j-1
        sum2 = sum2 + a(i,k) * a(j,k)
      end do
      a(i,j) = ( a(i,j) - sum2 ) / a(j,j)
    end do

  end do

  return
end
subroutine cosqb ( n, x, wsave )

!*****************************************************************************80
!
!! COSQB computes the fast cosine transform of quarter wave data.
!
!  Discussion:
!
!    COSQB computes a sequence from its representation in terms of a cosine
!    series with odd wave numbers.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!
!        4 * X_in(K) * cos ( ( 2 * K - 1 ) * ( I - 1 ) * PI / ( 2 * N ) )
!
!    COSQB is the unnormalized inverse of COSQF since a call of COSQB
!    followed by a call of COSQF will multiply the input sequence X by 4*N.
!
!    The array WSAVE must be initialized by calling COSQI.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array X.  The method is
!    more efficient when N is the product of small primes.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, the cosine series coefficients.
!    On output, the corresponding data vector.
!
!    Input, real WSAVE(3*N+15), contains data, depending on N, and
!    required by the algorithm.  The WSAVE array must be initialized by
!    calling COSQI.  A different WSAVE array must be used for each different
!    value of N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: tsqrt2 = 2.82842712474619D+00
  real ( kind = 8 ) wsave(3*n+15)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1

  if ( n < 2 ) then
    x(1) = 4.0D+00 * x(1)
  else if ( n == 2 ) then
    x1 = 4.0D+00 * ( x(1) + x(2) )
    x(2) = tsqrt2 * ( x(1) - x(2) )
    x(1) = x1
  else
    call cosqb1 ( n, x, wsave(1), wsave(n+1) )
  end if

  return
end
subroutine cosqb1 ( n, x, w, xh )

!*****************************************************************************80
!
!! COSQB1 is a lower level routine used by COSQB.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, the cosine series coefficients.
!    On output, the corresponding data vector.
!
!    Input, real ( kind = 8 ) W(N).
!
!    Input, real ( kind = 8 ) XH(2*N+15).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xh(2*n+15)
  real ( kind = 8 ) xim1

  ns2 = ( n + 1 ) / 2

  do i = 3, n, 2
    xim1 = x(i-1) + x(i)
    x(i) = x(i) - x(i-1)
    x(i-1) = xim1
  end do

  x(1) = x(1) + x(1)

  if ( mod ( n, 2 ) == 0 ) then
    x(n) = 2.0D+00 * x(n)
  end if

  call dfftb ( n, x, xh )

  do k = 2, ns2
    kc = n + 2 - k
    xh(k) = w(k-1) * x(kc) + w(kc-1) * x(k)
    xh(kc) = w(k-1) * x(k) - w(kc-1) * x(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    x(ns2+1) = w(ns2) * ( x(ns2+1) + x(ns2+1) )
  end if

  do k = 2, ns2
    kc = n + 2 - k
    x(k) = xh(k) + xh(kc)
    x(kc) = xh(k) - xh(kc)
  end do

  x(1) = 2.0D+00 * x(1)

  return
end
subroutine cosqf ( n, x, wsave )

!*****************************************************************************80
!
!! COSQF computes the fast cosine transform of quarter wave data.
!
!  Discussion:
!
!    COSQF computes the coefficients in a cosine series representation
!    with only odd wave numbers.
!
!    COSQF is the unnormalized inverse of COSQB since a call of COSQF
!    followed by a call of COSQB will multiply the input sequence X
!    by 4*N.
!
!    The array WSAVE must be initialized by calling COSQI.
!
!    The transform is defined by:
!
!      X_out(I) = X_in(1) + sum ( 2 <= K <= N )
!
!        2 * X_in(K) * cos ( ( 2 * I - 1 ) * ( K - 1 ) * PI / ( 2 * N ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array X.  The method is
!    more efficient when N is the product of small primes.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, the data to be transformed.
!    On output, the transformed data.
!
!    Input, real ( kind = 8 ) WSAVE(3*N+15), contains data, depending on N, and
!    required by the algorithm.  The WSAVE array must be initialized by
!    calling COSQI.  A different WSAVE array must be used for each different
!    value of N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: sqrt2 = 1.4142135623731D+00
  real ( kind = 8 ) tsqx
  real ( kind = 8 ) wsave(3*n+15)
  real ( kind = 8 ) x(n)

  if ( n < 2 ) then

  else if ( n == 2 ) then
    tsqx = sqrt2 * x(2)
    x(2) = x(1) - tsqx
    x(1) = x(1) + tsqx
  else
    call cosqf1 ( n, x, wsave(1), wsave(n+1) )
  end if

  return
end
subroutine cosqf1 ( n, x, w, xh )

!*****************************************************************************80
!
!! COSQF1 is a lower level routine used by COSQF.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, the data to be transformed.
!    On output, the transformed data.
!
!    Input, real ( kind = 8 ) W(N).
!
!    Input, real ( kind = 8 ) XH(2*N+15).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xh(2*n+15)
  real ( kind = 8 ) xim1

  ns2 = ( n + 1 ) / 2

  do k = 2, ns2
    kc = n + 2 - k
    xh(k) = x(k) + x(kc)
    xh(kc) = x(k) - x(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    xh(ns2+1) = x(ns2+1) + x(ns2+1)
  end if

  do k = 2, ns2
    kc = n+2-k
    x(k) = w(k-1) * xh(kc) + w(kc-1) * xh(k)
    x(kc) = w(k-1) * xh(k) - w(kc-1) * xh(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    x(ns2+1) = w(ns2) * xh(ns2+1)
  end if

  call dfftf ( n, x, xh )

  do i = 3, n, 2
    xim1 = x(i-1) - x(i)
    x(i) = x(i-1) + x(i)
    x(i-1) = xim1
  end do

  return
end
subroutine cosqi ( n, wsave )

!*****************************************************************************80
!
!! COSQI initializes WSAVE, used in COSQF and COSQB.
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.
!    The method is more efficient when N is the product of small primes.
!
!    Output, real ( kind = 8 ) WSAVE(3*N+15), contains data, depending on N,
!    and required by the COSQB and COSQF algorithms.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) dt
  integer ( kind = 4 ) k
  real ( kind = 8 ) wsave(3*n+15)

  dt = 0.5D+00 * pi / real ( n, kind = 8 )

  do k = 1, n
    wsave(k) = cos ( real ( k, kind = 8 ) * dt )
  end do

  call dffti ( n, wsave(n+1) )

  return
end
subroutine cost ( n, x, wsave )

!*****************************************************************************80
!
!! COST computes the discrete Fourier cosine transform of an even sequence.
!
!  Discussion:
!
!    COST is the unnormalized inverse of itself since a call of COST
!    followed by another call of COST will multiply the input sequence
!    X by 2*(N-1).
!
!    The array WSAVE must be initialized by calling COSTI.
!
!    The transform is defined by:
!
!      X_out(I) = X_in(1) + (-1) **(I-1) * X_in(N) + sum ( 2 <= K <= N-1 )
!
!        2 * X_in(K) * cos ( ( K - 1 ) * ( I - 1 ) * PI / ( N - 1 ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The method is more efficient when N-1 is the product of 
!    small primes.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real ( kind = 8 ) WSAVE(3*N+15).
!    The WSAVE array must be initialized by calling COSTI.  A different
!    array must be used for each different value of N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tx2
  real ( kind = 8 ) wsave(3*n+15)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1h
  real ( kind = 8 ) x1p3
  real ( kind = 8 ) xi
  real ( kind = 8 ) xim2

  ns2 = n / 2

  if ( n <= 1 ) then
    return
  end if

  if ( n == 2 ) then
    x1h = x(1) + x(2)
    x(2) = x(1) - x(2)
    x(1) = x1h
    return
  end if

  if ( n == 3 ) then
    x1p3 = x(1) + x(3)
    tx2 = x(2) + x(2)
    x(2) = x(1) - x(3)
    x(1) = x1p3 + tx2
    x(3) = x1p3 - tx2
    return
  end if

  c1 = x(1) - x(n)
  x(1) = x(1) + x(n)

  do k = 2, ns2
    kc = n + 1 - k
    t1 = x(k) + x(kc)
    t2 = x(k) - x(kc)
    c1 = c1 + wsave(kc) * t2
    t2 = wsave(k) * t2
    x(k) = t1 - t2
    x(kc) = t1 + t2
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(ns2+1) = x(ns2+1) + x(ns2+1)
  end if

  call dfftf ( n-1, x, wsave(n+1) )

  xim2 = x(2)
  x(2) = c1

  do i = 4, n, 2
    xi = x(i)
    x(i) = x(i-2) - x(i-1)
    x(i-1) = xim2
    xim2 = xi
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n) = xim2
  end if

  return
end
subroutine costi ( n, wsave )

!*****************************************************************************80
!
!! COSTI initializes WSAVE, used in COST.
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The method is more efficient when N-1 is the product of 
!    small primes.
!
!    Output, real ( kind = 8 ) WSAVE(3*N+15), contains data, depending on N, 
!    and required by the COST algorithm.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) dt
  integer ( kind = 4 ) k
  real ( kind = 8 ) wsave(3*n+15)

  if ( n <= 3 ) then
    return
  end if

  dt = pi / real ( n - 1, kind = 8 )

  do k = 2, ( n / 2 )
    wsave(k)     = 2.0D+00 * sin ( real ( k - 1, kind = 8 ) * dt )
    wsave(n+1-k) = 2.0D+00 * cos ( real ( k - 1, kind = 8 ) * dt )
  end do

  call dffti ( n-1, wsave(n+1) )

  return
end
function csevl ( x, cs, n )

!*****************************************************************************80
!
!! CSEVL evaluates an N term Chebyshev series.
!
!  Modified:
!
!    15 April 2003
!
!  Reference:
!
!    R Broucke,
!    Algorithm 446,
!    Communications of the ACM,
!    Volume 16, page 254, 1973.
!
!    Leslie Fox, Ian Parker,
!    Chebyshev Polynomials in Numerical Analysis,
!    Oxford Press, page 56.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument at which the series is to be
!    evaluated.  X must satisfy -1.0 <= X <= 1.0.
!
!    Input, real ( kind = 8 ) CS(N), the array of N terms of a Chebyshev series.
!    In evaluating CS, only half the first coefficient is summed.
!
!    Input, integer ( kind = 4 ) N, the number of terms in array CS.
!    N must be at least 1, and no more than 1000.
!
!    Output, real ( kind = 8 ) CSEVL, the value of the Chebyshev series.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) cs(n)
  real ( kind = 8 ) csevl
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms N is less than 1.'
    stop
  end if

  if ( 1000 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  The number of terms is more than 1000.'
    stop
  end if

  if ( x < -1.0D+00 .or. 1.0D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  The input argument X is outside the interval [-1,1].'
    stop
  end if

  b1 = 0.0D+00
  b0 = 0.0D+00

  do i = n, 1, -1
    b2 = b1
    b1 = b0
    b0 = 2.0D+00 * x * b1 - b2 + cs(i)
  end do

  csevl = 0.5D+00 * ( b0 - b2 )

  return
end
subroutine d1fcn ( n, x, g )

!*****************************************************************************80
!
!! D1FCN is a dummy routine for evaluating the gradient vector.
!
!  Discussion:
!
!    We assume that F is a scalar function of N variables X.  The routine
!    is to compute the vector G where G(I) = d F/d X(I).
!
!  Modified:
!
!    16 April 2001
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of X, and order of A.
!
!    Input, real ( kind = 8 ) X(N), the point at which the gradient
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) G(N), the gradient vector..
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'D1FCN - Fatal error!'
  write ( *, '(a)' ) '  This is a dummy routine.'
  write ( *, '(a)' ) '  The user is required to replace it with a'
  write ( *, '(a)' ) '  routine that computes the gradient of F.'

  stop
end
function d1mach ( i )

!*****************************************************************************80
!
!! D1MACH returns double precision real machine constants.
!
!  Discussion:
!
!    Assuming that the internal representation of a double precision real
!    number is in base B, with T the number of base-B digits in the mantissa,
!    and EMIN the smallest possible exponent and EMAX the largest possible
!    exponent, then
!
!      D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!      D1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
!      D1MACH(3) = B**(-T), the smallest relative spacing.
!      D1MACH(4) = B**(1-T), the largest relative spacing.
!      D1MACH(5) = log10(B).
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer
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
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 5.
!
!    Output, real ( kind = 8 ) D1MACH, the value of the chosen parameter.
!
  implicit none

  real ( kind = 8 ) d1mach
  integer ( kind = 4 ) i

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    d1mach = 0.0D+00
    stop
  else if ( i == 1 ) then
    d1mach = 4.450147717014403D-308
  else if ( i == 2 ) then
    d1mach = 8.988465674311579D+307
  else if ( i == 3 ) then
    d1mach = 1.110223024625157D-016
  else if ( i == 4 ) then
    d1mach = 2.220446049250313D-016
  else if ( i == 5 ) then
    d1mach = 0.301029995663981D+000
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    d1mach = 0.0D+00
    stop
  end if

  return
end
subroutine d1mpyq ( m, n, a, lda, v, w )

!*****************************************************************************80
!
!! D1MPYQ computes A*Q, where Q is the product of Householder transformations.
!
!  Discussion:
!
!    Given an M by N matrix A, this subroutine computes A * Q where
!    Q is the product of 2 * (N - 1) transformations
!
!      GV(N-1) * ... * GV(1) * GW(1) * ... * GW(N-1)
!
!    and GV(I), GW(I) are Givens rotations in the (I,N) plane which
!    eliminate elements in the I-th and N-th planes, respectively.
!    Q itself is not given, rather the information to recover the
!    GV, GW rotations is supplied.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the M by N array.
!    On input, the matrix A to be postmultiplied by the orthogonal matrix Q.
!    On output, the value of A*Q.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must not
!    be less than M.
!
!    Input, real ( kind = 8 ) V(N), W(N), contain the information necessary
!    to recover the Givens rotations GV and GW.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) s
  real ( kind = 8 ) temp
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) w(n)
!
!  Apply the first set of Givens rotations to A.
!
  do j = n-1, 1, -1

     if ( 1.0D+00 < abs ( v(j) ) ) then
       c = 1.0D+00 / v(j)
       s = sqrt ( 1.0D+00 - c * c )
     else
       s = v(j)
       c = sqrt ( 1.0D+00 - s * s )
     end if

     do i = 1, m
        temp =   c * a(i,j) - s * a(i,n)
        a(i,n) = s * a(i,j) + c * a(i,n)
        a(i,j) = temp
     end do

  end do
!
!  Apply the second set of Givens rotations to A.
!
  do j = 1, n-1

     if ( 1.0D+00 < abs ( w(j) ) ) then
       c = 1.0D+00 / w(j)
       s = sqrt ( 1.0D+00 - c * c )
     else
       s = w(j)
       c = sqrt ( 1.0D+00 - s * s )
     end if

     do i = 1, m
        temp =     c * a(i,j) + s * a(i,n)
        a(i,n) = - s * a(i,j) + c * a(i,n)
        a(i,j) = temp
     end do

  end do

  return
end
subroutine d2fcn ( nr, n, x, a )

!*****************************************************************************80
!
!! D2FCN is a dummy version of a routine that computes the second derivative.
!
!  Discussion:
!
!    We assume that F is a scalar function of N variables X.  The routine
!    is to compute the matrix H where H(I,J) = d d F / d X(I) d X(J).
!
!  Modified:
!
!    16 April 2001
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the leading dimension of A, which must be 
!    at least N.
!
!    Input, integer ( kind = 4 ) N, the dimension of X, and order of A.
!
!    Input, real ( kind = 8 ) X(N), the point at which the Hessian matrix
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) A(NR,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'D2FCN - Fatal error!'
  write ( *, '(a)' ) '  This is a dummy routine.'
  write ( *, '(a)' ) '  The user is required to replace it with a'
  write ( *, '(a)' ) '  routine that computes the Hessian matrix of F.'

  stop
end
function d9lgmc ( x )

!*****************************************************************************80
!
!! D9LGMC computes the log gamma correction factor.
!
!  Discussion:
!
!    The routine computes the log gamma correction factor for 10 <= X
!    so that
!
!      log ( gamma ( x ) ) =
!        log ( sqrt ( 2 * pi ) ) + ( x - 0.5 ) * log ( x ) - x + d9lgmc ( x )
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the log gamma function.
!    X must be at least 10.
!
!    Output, real ( kind = 8 ) D9LGMC, the correction.
!
  implicit none

  real ( kind = 8 ), save, dimension ( 6 ) :: algmcs = (/ &
     0.166638948045186D+00, -0.0000138494817606D+00,  0.0000000098108256D+00, &
    -0.0000000000180912D+00, 0.0000000000000622D+00, -0.0000000000000003D+00 /)
  real ( kind = 8 ) arg
  real ( kind = 8 ) csevl
  real ( kind = 8 ) d9lgmc
  integer ( kind = 4 ) inits
  integer ( kind = 4 ), save :: nalgm = 0
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xbig = 0.0D+00
  real ( kind = 8 ), save :: xmax = 0.0D+00

  if ( nalgm == 0 ) then
    nalgm = inits ( algmcs, 6, epsilon ( algmcs ) )
    xbig = 1.0D+00 / sqrt ( epsilon ( xbig ) )
    xmax = exp ( min ( log ( huge ( xmax ) / 12.0D+00 ), &
                     -log ( 12.0D+00 * tiny ( xmax ) ) ) )
  end if

  if ( x < 10.0D+00 ) then

    call xerror ( 'D9LGMC - 10 <= x required', 1, 2 )

  else if ( x < xbig ) then

    arg = 2.0D+00 * ( 10.0D+00 / x )**2 - 1.0D+00
    d9lgmc = csevl ( arg, algmcs, nalgm ) / x

  else if ( x < xmax ) then

    d9lgmc = 1.0D+00 / ( 12.0D+00 * x )

  else

    d9lgmc = 0.0D+00
    call xerror ( 'D9LGMC - X so big d9lgmc underflows', 2, 1)

  end if

  return
end
function damax ( n, x, incx )

!*****************************************************************************80
!
!! DAMAX returns the maximum absolute value of the entries in a vector.
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Output, real ( kind = 8 ) DAMAX, the maximum absolute value of
!    an element of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real ( kind = 8 ) damax
  real ( kind = 8 ) x(*)

  if ( n <= 0 ) then

    damax = 0.0D+00

  else if ( n == 1 ) then

    damax = abs ( x(1) )

  else if ( incx == 1 ) then

    damax = abs ( x(1) )

    do i = 2, n
      if ( damax < abs ( x(i) ) ) then
        damax = abs ( x(i) )
      end if
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    damax = abs ( x(ix) )
    ix = ix + incx

    do i = 2, n
      if ( damax < abs ( x(ix) ) ) then
        damax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function dasum ( n, x, incx )

!*****************************************************************************80
!
!! DASUM takes the sum of the absolute values of a vector.
!
!  Modified:
!
!    15 February 2001
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.  INCX must not be negative.
!
!    Output, real ( kind = 8 ) DASUM, the sum of the absolute values of X.
!
  implicit none

  real ( kind = 8 ) dasum
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) n
  real x(*)

  dasum = sum ( abs ( x(1:1+(n-1)*incx:incx) ) )

  return
end
subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    Uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in DX and DY.
!
!    Input, real ( kind = 8 ) DA, the multiplier of DX.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of DX.
!
!    Input/output, real ( kind = 8 ) DY(*), the second vector.
!    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    entries of DY.
!
  implicit none

  real ( kind = 8 ) da
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    return
  end if

  if ( da  == 0.0D+00 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dy(iy) = dy(iy) + da * dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 4 )

    do i = 1, m
      dy(i) = dy(i) + da * dx(i)
    end do

    do i = m+1, n, 4
      dy(i  ) = dy(i  ) + da * dx(i  )
      dy(i+1) = dy(i+1) + da * dx(i+1)
      dy(i+2) = dy(i+2) + da * dx(i+2)
      dy(i+3) = dy(i+3) + da * dx(i+3)
    end do

  end if

  return
end
subroutine ddcor ( dfdy, el, fa, h, impl, ipvt, matdim, miter, ml, mu, n, &
  nde, nq, t, users, y, yh, ywt, evalfa, save1, save2, a, d, jstate )

!*****************************************************************************80
!
!! DDCOR computes corrections to the Y array of DDRIV3.
!
!  Discussion:
!
!    In the case of functional iteration, update Y directly from the
!    result of the last call to F.
!
!    In the case of the chord method, compute the corrector error and
!    solve the linear system with that as right hand side and DFDY as
!    coefficient matrix, using the lu decomposition if miter is 1, 2, 4,
!    or 5.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) matdim
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(matdim,*)
  real ( kind = 8 ) d
  real ( kind = 8 ) dfdy(matdim,*)
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) el(13,12)
  logical evalfa
  external fa
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) impl
  integer ( kind = 4 ) ipvt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstate
  integer ( kind = 4 ) miter
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) mw
  integer ( kind = 4 ) nde
  integer ( kind = 4 ) nq
  real ( kind = 8 ) save1(*)
  real ( kind = 8 ) save2(*)
  real ( kind = 8 ) t
  external users
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) yh(n,*)
  real ( kind = 8 ) ywt(*)

  if ( miter == 0 ) then

    save1(1:n) = ( h * save2(1:n) - yh(1:n,2) - save1(1:n) ) / ywt(1:n)

    d = dnrm2 ( n, save1, 1 ) / sqrt ( real ( n, kind = 8 ) )

    save1(1:n) = h * save2(1:n) - yh(1:n,2)

  else if ( miter == 1 .or. miter == 2 ) then

    if ( impl == 0 ) then

      save2(1:n) = h * save2(1:n) - yh(1:n,2) - save1(1:n)

    else if ( impl == 1 ) then

      if ( evalfa ) then
        call fa ( n, t, y, a, matdim, ml, mu, nde )
        if ( n == 0 ) then
          jstate = 9
          return
        end if
      else
        evalfa = .true.
      end if

      save2(1:n) = h * save2(1:n)

      do j = 1,n
        save2(1:n) = save2(1:n) - a(1:n,j) * ( yh(j,2) + save1(j) )
      end do

    else if ( impl == 2 ) then

      if ( evalfa ) then
        call fa ( n, t, y, a, matdim, ml, mu, nde )
        if ( n == 0 ) then
          jstate = 9
          return
        end if
      else
        evalfa = .true.
      end if

      save2(1:n) = h * save2(1:n) - a(1:n,1) * (yh(i,2) + save1(1:n))

    end if

    call dgesl ( dfdy, matdim, n, ipvt, save2, 0 )

    save1(1:n) = save1(1:n) + save2(1:n)
    save2(1:n) = save2(1:n) / ywt(1:n)

    d = dnrm2 ( n, save2, 1 ) / sqrt ( real ( n, kind = 8 ) )

  else if ( miter == 4 .or. miter == 5 ) then

    if ( impl == 0 ) then

      save2(1:n) = h * save2(1:n) - yh(1:n,2) - save1(1:n)

    else if ( impl == 1 ) then
      if ( evalfa ) then
        call fa ( n, t, y, a(ml+1,1), matdim, ml, mu, nde )
        if ( n == 0 ) then
          jstate = 9
          return
        end if
      else
        evalfa = .true.
      end if

      save2(1:n) = h * save2(1:n)

      mw = ml + 1 + mu

      do j = 1, n
        i1 = max ( ml+1, mw+1-j )
        i2 = min ( mw+n-j, mw+ml )
        do i = i1,i2
          i3 = i + j - mw
          save2(i3) = save2(i3) - a(i,j)*(yh(j,2) + save1(j))
        end do
      end do

    else if ( impl == 2 ) then

      if ( evalfa ) then
        call fa ( n, t, y, a, matdim, ml, mu, nde )
        if ( n == 0 ) then
          jstate = 9
          return
        end if
      else
        evalfa = .true.
      end if

      save2(1:n) = h * save2(1:n) - a(1:n,1)*(yh(1:n,2) + save1(1:n))

    end if

    call dgbsl ( dfdy, matdim, n, ml, mu, ipvt, save2, 0 )

    save1(1:n) = save1(1:n) + save2(1:n)
    save2(1:n) = save2(1:n) / ywt(1:n)

    d = dnrm2 ( n, save2, 1 ) / sqrt ( real ( n, kind = 8 ) )

  else if ( miter == 3 ) then

    iflag = 2

    call users ( y, yh(1,2), ywt, save1, save2, t, h, el(1,nq), impl, &
      n, nde, iflag )

    if ( n == 0 ) then
      jstate = 10
      return
    end if

    save1(1:n) = save1(1:n) + save2(1:n)
    save2(1:n) = save2(1:n) / ywt(1:n)

    d = dnrm2 ( n, save2, 1) / sqrt ( real ( n, kind = 8 ) )

  end if
end
subroutine ddcst ( maxord, mint, iswflg, el, tq )

!*****************************************************************************80
!
!! DDCST sets coefficients used by the core integrator DDSTP.
!
!  Discussion:
!
!    EL and TQ depend upon MINT, and are calculated
!    for orders 1 to maxord(<= 12).  for each order NQ, the coefficients
!    EL are calculated from the generating polynomial:
!      l(t) = el(1,nq) + el(2,nq) * t + ... + el(nq+1,nq) * t**nq.
!    for the implicit adams methods, l(t) is given by
!      dl/dt = (1+t)*(2+t)* ... *(nq-1+t)/k,   l(-1) = 0,
!      where      k = factorial(nq-1).
!    for the gear methods,
!      l(t) = (1+t)*(2+t)* ... *(nq+t)/k,
!      where      k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!    for each order nq, there are three components of tq.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXORD, the maximum order to calculate.
!
!    Input, integer ( kind = 4 ) MINT, 1 for Adams methods, 2 for Gear.
!
!    Input, integer ( kind = 4 ) ISWFLG, indicates whether the constants used
!    in the stiffness test should be calculated.
!
!    Output, real ( kind = 8 ) EL(13,12), used in specifying the basic method.
!
!    Output, real ( kind = 8 ) TQ(3,12), used in adjusting the stepsize in
!    relation to truncation error.
!
  implicit none

  real ( kind = 8 ) el(13,12)
  real ( kind = 8 ) factrl(12)
  real ( kind = 8 ) gamma(14)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iswflg
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxord
  integer ( kind = 4 ) mint
  integer ( kind = 4 ) mxrd
  real ( kind = 8 ) sum2
  real ( kind = 8 ) tq(3,12)

  factrl(1) = 1.0D+00
  do i = 2, maxord
    factrl(i) = real ( i, kind = 8 ) * factrl(i-1)
  end do
!
!  Compute Adams coefficients
!
  if ( mint == 1 ) then

    gamma(1) = 1.0D+00
    do i = 1, maxord + 1
      sum2 = 0.0D+00
      do j = 1, i
        sum2 = sum2 - gamma(j) / real ( i - j + 2, kind = 8 )
      end do
      gamma(i+1) = sum2
    end do

    el(1,1) = 1.0D+00
    el(2,1) = 1.0D+00
    el(2,2) = 1.0D+00
    el(3,2) = 1.0D+00

    do j = 3, maxord
      el(2,j) = factrl(j-1)
      do i = 3,j
        el(i,j) = real ( j - 1, kind = 8 ) * el(i,j-1) + el(i-1,j-1)
      end do
      el(j+1,j) = 1.0D+00
    end do

    do j = 2, maxord
      el(1,j) = el(1,j-1) + gamma(j)
      el(2,j) = 1.0D+00
      do i = 3, j+1
        el(i,j) = el(i,j) / ( real ( i - 1, kind = 8 ) * factrl(j-1) )
      end do
    end do

    do j = 1, maxord
      tq(1,j) = -1.0D+00 / ( factrl(j) * gamma(j) )
      tq(2,j) = -1.0D+00 / gamma(j+1)
      tq(3,j) = -1.0D+00 / gamma(j+2)
    end do
!
!  Compute Gear coefficients.
!
  else if ( mint == 2 ) then

    el(1,1) = 1.0D+00
    el(2,1) = 1.0D+00
    do j = 2, maxord
      el(1,j) = factrl(j)
      do i = 2, j
        el(i,j) = real ( j, kind = 8 ) * el(i,j-1) + el(i-1,j-1)
      end do
      el(j+1,j) = 1.0D+00
    end do

    sum2 = 1.0D+00
    do j = 2, maxord
      sum2 = sum2 + 1.0D+00 / real ( j, kind = 8 )
      do i = 1, j+1
        el(i,j) = el(i,j) / ( factrl(j) * sum2 )
      end do
    end do

    do j = 1, maxord
      if ( 1 < j ) then
        tq(1,j) = 1.0D+00 / factrl(j-1)
      end if
      tq(2,j) = real ( j + 1, kind = 8 ) / el(1,j)
      tq(3,j) = real ( j + 2, kind = 8 ) / el(1,j)
    end do

  end if
!
!  Compute constants used in the stiffness test.
!  these are the ratio of tq(2,nq) for the gear
!  methods to those for the adams methods.
!
  if ( iswflg == 3 ) then
    mxrd = min ( maxord, 5 )
    if ( mint == 2 ) then
      gamma(1) = 1.0D+00
      do i = 1, mxrd
        sum2 = 0.0D+00
        do j = 1, i
          sum2 = sum2 - gamma(j) / real ( i - j + 2, kind = 8 )
        end do
        gamma(i+1) = sum2
      end do
    end if

    sum2 = 1.0D+00
    do i = 2, mxrd
      sum2 = sum2 + 1.0D+00 / real ( i, kind = 8 )
      el(1+i,1) = - real ( i + 1, kind = 8 ) * sum2 * gamma(i+1)
    end do

  end if
end
subroutine ddntl ( eps, f, fa, hmax, hold, impl, jtask, matdim, maxord, &
  mint, miter, ml, mu, n, nde, save1, t, uround, users, y, ywt, h, mntold, &
  mtrold, nfe, rc, yh, a, convrg, el, fac, ier, ipvt, nq, nwait, rh, rmax, &
  save2, tq, trend, iswflg, jstate )

!*****************************************************************************80
!
!! DDNTL sets parameters for DDSTP.
!
!  Discussion:
!
!    DDNTL is called on the first call to DDSTP, on an internal restart, or
!    when the user has altered MINT, miter, and/or h.
!
!    On the first call, the order is set to 1 and the initial derivatives
!    are calculated.  RMAX is the maximum ratio by which h can be
!    increased in one step.  it is initially rminit to compensate
!    for the small initial h, but then is normally equal to rmnorm.
!    if a failure occurs (in corrector convergence or error test), rmax
!    is set at rmfail for the next increase.
!    if the caller has changed mint, or if jtask = 0, DDCST is called
!    to set the coefficients of the method.  if the caller has changed h,
!    yh must be rescaled.  if h or mint has been changed, nwait is
!    reset to nq + 2 to prevent further increases in h for that many
!    steps.  also, rc is reset.  rc is the ratio of new to old values of
!    the coefficient l(0)*h.  if the caller has changed miter, rc is
!    set to 0 to force the partials to be updated, if partials are used.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) matdim
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(matdim,*)
  logical convrg
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) el(13,12)
  real ( kind = 8 ) eps
  external f
  external fa
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hold
  integer ( kind = 4 ) i
  logical ier
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) impl
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(*)
  integer ( kind = 4 ) iswflg
  integer ( kind = 4 ) jstate
  integer ( kind = 4 ) jtask
  integer ( kind = 4 ) maxord
  integer ( kind = 4 ) mint
  integer ( kind = 4 ) miter
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mntold
  integer ( kind = 4 ) mtrold
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nde
  integer ( kind = 4 ) nfe
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nwait
  real ( kind = 8 ) oldl0
  real ( kind = 8 ) rc
  real ( kind = 8 ) rh
  real ( kind = 8 ) rmax
  real ( kind = 8 ), parameter :: rminit = 10000.0D+00
  real ( kind = 8 ) save1(*)
  real ( kind = 8 ) save2(*)
  real ( kind = 8 ) smax
  real ( kind = 8 ) smin
  real ( kind = 8 ) sum0
  real ( kind = 8 ) sum2
  real ( kind = 8 ) t
  real ( kind = 8 ) tq(3,12)
  real ( kind = 8 ) trend
  real ( kind = 8 ) uround
  external users
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) yh(n,*)
  real ( kind = 8 ) ywt(*)

  ier = .false.

  if ( 0 <= jtask ) then

    if ( jtask == 0 ) then
      call ddcst ( maxord, mint, iswflg, el, tq )
      rmax = rminit
    end if

    rc = 0.0D+00
    convrg = .false.
    trend = 1.0D+00
    nq = 1
    nwait = 3

    call f ( n, t, y, save2 )

    if ( n == 0 ) then
      jstate = 6
      return
    end if

    nfe = nfe + 1

    if ( impl /= 0 ) then
      if ( miter == 3 ) then
        iflag = 0
        call users ( y, yh, ywt, save1, save2, t, h, el, impl, n, nde, iflag )
        if ( n == 0 ) then
          jstate = 10
          return
        end if
      else if ( impl == 1 ) then
        if ( miter == 1 .or. miter == 2 ) then
          call fa ( n, t, y, a, matdim, ml, mu, nde )
          if ( n == 0 ) then
            jstate = 9
            return
          end if
          call dgefa ( a, matdim, n, ipvt, info )
          if ( info /= 0 ) then
            ier = .true.
            return
          end if
          call dgesl ( a, matdim, n, ipvt, save2, 0 )
        else if (miter == 4 .or. miter == 5) then
          call fa ( n, t, y, a(ml+1,1), matdim, ml, mu, nde )
          if ( n == 0 ) then
            jstate = 9
            return
          end if
          call dgbfa ( a, matdim, n, ml, mu, ipvt, info )
          if ( info /= 0 ) then
            ier = .true.
            return
          end if
          call dgbsl ( a, matdim, n, ml, mu, ipvt, save2, 0 )
        end if

      else if ( impl == 2 ) then

        call fa ( n, t, y, a, matdim, ml, mu, nde )
        if ( n == 0 ) then
          jstate = 9
          return
        end if

        do i = 1, nde

          if ( a(i,1) == 0.0D+00 ) then
            ier = .true.
            return
          else
            save2(i) = save2(i) / a(i,1)
          end if

        end do

        do i = nde+1,n
          a(i,1) = 0.0D+00
        end do

      end if

    end if

    save1(1:nde) = save2(1:nde) / ywt(1:nde)

    sum2 = dnrm2 ( nde, save1, 1 )
    sum0 = 1.0D+00 / max ( 1.0D+00 , abs ( t ) )
    smax = max ( sum0, sum2 )
    smin = min ( sum0, sum2 )
    sum2 = smax * sqrt ( 1.0D+00 + ( smin / smax )**2 ) &
      / sqrt ( real ( nde, kind = 8 ) )
    h = sign ( min ( 2.0D+00 * eps / sum2, abs ( h ) ), h)
    yh(1:n,2) = h * save2(1:n)

    if ( miter == 2 .or. miter == 5 .or. iswflg == 3 ) then
      do i = 1,n
        fac(i) = sqrt ( uround )
      end do
    end if

  else

    if ( miter /= mtrold ) then
      mtrold = miter
      rc = 0.0D+00
      convrg = .false.
    end if

    if ( mint /= mntold ) then
      mntold = mint
      oldl0 = el(1,nq)
      call ddcst ( maxord, mint, iswflg, el, tq )
      rc = rc * el(1,nq) / oldl0
      nwait = nq + 2
    end if

    if ( h /= hold ) then
      nwait = nq + 2
      rh = h / hold
      call ddscl ( hmax, n, nq, rmax, hold, rc, rh, yh )
    end if

  end if

  return
end
subroutine ddntp ( h, k, n, nq, t, tout, yh, y )

!*****************************************************************************80
!
!! DDNTP interpolates the K-th derivative of the ODE solution Y at TOUT.
!
!  Discussion:
!
!    The routine uses the data in the YH array.  If K has a value greater
!    than NQ, the NQ-th derivative is calculated.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) factor
  real ( kind = 8 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kused
  integer ( kind = 4 ) nq
  real ( kind = 8 ) r
  real ( kind = 8 ) t
  real ( kind = 8 ) tout
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) yh(n,*)

  if ( k == 0 ) then

    y(1:n) = yh(1:n,nq+1)

    r = ( tout - t ) / h

    do jj = 1, nq
      j = nq + 1 - jj
      y(1:n) = yh(1:n,j) + r * y(1:n)
    end do

  else

    kused = min ( k, nq )
    factor = 1.0D+00
    do kk = 1, kused
      factor = factor * real ( nq + 1 - kk, kind = 8 )
    end do

    y(1:n) = factor * yh(1:n,nq+1)

    do jj = kused+1,nq

      j = k + 1 + nq - jj
      factor = 1.0D+00

      do kk = 1, kused
        factor = factor * real ( j - kk, kind = 8 )
      end do

      y(1:n) = factor * yh(1:n,j) + r * y(1:n)

    end do

    y(1:n) = y(1:n) * h**(-kused)

  end if
end
function ddot ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries
!    in X.
!
!    Input, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive entries
!    in Y.
!
!    Output, real DDOT, the sum of the product of the corresponding
!    entries of X and Y.
!
  implicit none

  real ( kind = 8 ) ddot
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  ddot = 0.0D+00
  dtemp = 0.0D+00

  if ( n <= 0 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dtemp = dtemp + dx(ix) * dy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 5 )

    do i = 1, m
      dtemp = dtemp + dx(i) * dy(i)
    end do

    do i = m+1, n, 5

      dtemp = dtemp + dx(i  ) * dy(i  ) &
                    + dx(i+1) * dy(i+1) &
                    + dx(i+2) * dy(i+2) &
                    + dx(i+3) * dy(i+3) &
                    + dx(i+4) * dy(i+4)
    end do

  end if

  ddot = dtemp

  return
end
subroutine ddpsc ( ksgn, n, nq, yh )

!*****************************************************************************80
!
!! DDPSC computes the predicted YH values.
!
!  Discussion:
!
!    The routine effectively multiplies the YH array by the Pascal
!    triangle matrix when KSGN is +1, and performs the inverse function
!    when KSGN is -1.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KSGN, indicates which operation is to be performed.
!
!    Input, integer ( kind = 4 ) N, ?
!
!    Input, integer ( kind = 4 ) NQ, ?
!
!    Input/output, real ( kind = 8 ) YH(N,*), ?
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) ksgn
  integer ( kind = 4 ) nq
  real ( kind = 8 ) yh(n,*)

  if ( 0 < ksgn ) then

    do j1 = 1, nq
      do j2 = j1, nq
        j = nq - j2 + j1
        yh(1:n,j) = yh(1:n,j) + yh(1:n,j+1)
      end do
    end do

  else

    do j1 = 1, nq
      do j2 = j1, nq
        j = nq - j2 + j1
        yh(1:n,j) = yh(1:n,j) - yh(1:n,j+1)
      end do
    end do
  end if

  return
end
subroutine ddpst ( el, f, fa, h, impl, jacobn, matdim, miter, ml, mu, n, nde, &
  nq, save2, t, users, y, yh, ywt, uround, nfe, nje, a, dfdy, fac, ier, ipvt, &
  save1, iswflg, bnd, jstate )

!*****************************************************************************80
!
!! DDPST is called to reevaluate the partial derivatives.
!
!  Discussion:
!
!    If MITER is 1, 2, 4, or 5, the matrix
!    p = i - l(0)*h*jacobian
!    is stored in dfdy and subjected to LU decomposition, with the results
!    also stored in dfdy.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) matdim
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(matdim,*)
  real ( kind = 8 ) bl
  real ( kind = 8 ) bnd
  real ( kind = 8 ) bp
  real ( kind = 8 ) br
  real ( kind = 8 ) bu
  real ( kind = 8 ) dfdy(matdim,*)
  real ( kind = 8 ) dfdymx
  real ( kind = 8 ) diff
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) dy
  real ( kind = 8 ) el(13,12)
  external f
  external fa
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) facmin
  real ( kind = 8 ), parameter :: facmax = 0.5D+00
  real ( kind = 8 ) factor
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  logical ier
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) impl
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(*)
  integer ( kind = 4 ) iswflg
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  external jacobn
  integer ( kind = 4 ) jstate
  integer ( kind = 4 ) k
  integer ( kind = 4 ) miter
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) mw
  integer ( kind = 4 ) nde
  integer ( kind = 4 ) nfe
  integer ( kind = 4 ) nje
  integer ( kind = 4 ) nq
  real ( kind = 8 ) save1(*)
  real ( kind = 8 ) save2(*)
  real ( kind = 8 ) scale
  real ( kind = 8 ) t
  real ( kind = 8 ) uround
  external users
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) yh(n,*)
  real ( kind = 8 ) yj
  real ( kind = 8 ) ys
  real ( kind = 8 ) ywt(*)

  nje = nje + 1
  ier = .false.

  if ( miter == 1 .or. miter == 2 ) then

    if ( miter == 1 ) then

      call jacobn ( n, t, y, dfdy, matdim, ml, mu )

      if ( n == 0 ) then
        jstate = 8
        return
      end if

      if ( iswflg == 3 ) then
        bnd = dnrm2 ( n*n, dfdy, 1 )
      end if

      dfdy(1:n,1:n) = - el(1,nq) * h * dfdy(1:n,1:n)

    else if ( miter == 2 ) then

      br = uround**(0.875D+00 )
      bl = uround**(0.75D+00 )
      bu = uround**(0.25D+00 )
      bp = uround**(-0.15D+00 )
      facmin = uround**(0.78D+00 )

      do j = 1, n

        ys = max ( abs ( ywt(j) ), abs ( y(j) ) )

 120    continue

        dy = fac(j) * ys

        if ( dy == 0.0D+00 ) then
          if ( fac(j) < facmax ) then
            fac(j) = min ( 100.0D+00 * fac(j), facmax )
            go to 120
          else
            dy = ys
          end if
        end if

        if ( nq == 1 ) then
          dy = sign ( dy, save2(j) )
        else
          dy = sign ( dy, yh(j,3) )
        end if

        dy = (y(j) + dy) - y(j)
        yj = y(j)
        y(j) = y(j) + dy
        call f ( n, t, y, save1 )

        if ( n == 0 ) then
          jstate = 6
          return
        end if

        y(j) = yj
        factor = -el(1,nq) * h / dy
        dfdy(1:n,j) = ( save1(1:n) - save2(1:n) ) * factor
        diff = abs ( save2(1) - save1(1) )
        imax = 1

        do i = 2, n
          if ( diff < abs ( save2(i) - save1(i) ) ) then
            imax = i
            diff = abs ( save2(i) - save1(i) )
          end if
        end do
!
!  Step 2
!
        if ( 0.0D+00 < min ( abs ( save2(imax) ), abs ( save1(imax) ) ) ) then
          scale = max ( abs ( save2(imax) ), abs ( save1(imax) ) )
!
!  Step 3
!
          if ( bu * scale < diff ) then
            fac(j) = max ( facmin, fac(j) * 0.1D+00 )
          else if ( br * scale <= diff .and. diff <= bl * scale ) then
            fac(j) = min ( fac(j) * 10.0D+00, facmax )
!
!  Step 4
!
          else if ( diff < br * scale ) then
            fac(j) = min ( bp * fac(j), facmax )
          end if
        end if

      end do

      if ( iswflg == 3 ) bnd = dnrm2 ( n*n, dfdy, 1) / (-el(1,nq)*h)
      nfe = nfe + n
    end if

    if ( impl == 0 ) then

      do i = 1, n
        dfdy(i,i) = dfdy(i,i) + 1.0D+00
      end do

    else if ( impl == 1 ) then

      call fa ( n, t, y, a, matdim, ml, mu, nde )

      if ( n == 0 ) then
        jstate = 9
        return
      end if

      dfdy(1:n,1:n) = dfdy(1:n,1:n) + a(1:n,1:n)

    else if ( impl == 2 ) then

      call fa ( n, t, y, a, matdim, ml, mu, nde )

      if ( n == 0 ) then
        jstate = 9
        return
      end if

      do i = 1, nde
        dfdy(i,i) = dfdy(i,i) + a(i,1)
      end do

    end if

    call dgefa ( dfdy, matdim, n, ipvt, info )

    if ( info /= 0 ) then
      ier = .true.
    end if

  else if ( miter == 4 .or. miter == 5 ) then

    if ( miter == 4 ) then

      call jacobn ( n, t, y, dfdy(ml+1,1), matdim, ml, mu )

      if ( n == 0 ) then
        jstate = 8
        return
      end if

      factor = -el(1,nq) * h
      mw = ml + mu + 1

      do j = 1, n
        i1 = max ( ml+1, mw+1-j )
        i2 = min ( mw+n-j, mw+ml )
        do i = i1, i2
          dfdy(i,j) = factor * dfdy(i,j)
        end do
      end do

    else if ( miter == 5 ) then

      br = uround**(0.875D+00)
      bl = uround**(0.75D+00)
      bu = uround**(0.25D+00)
      bp = uround**(-0.15D+00)
      facmin = uround**(0.78D+00)
      mw = ml + mu + 1
      j2 = min ( mw, n )

      do j = 1, j2

        do k = j, n, mw

          ys = max ( abs (ywt(k) ), abs ( y(k) ) )

 280      continue

          dy = fac(k) * ys

          if ( dy == 0.0D+00 ) then
            if ( fac(k) < facmax ) then
              fac(k) = min ( 100.0D+00 * fac(k), facmax )
              go to 280
            else
              dy = ys
            end if
          end if

          if ( nq == 1 ) then
            dy = sign ( dy, save2(k) )
          else
            dy = sign ( dy, yh(k,3) )
          end if

          dy = (y(k) + dy) - y(k)
          dfdy(mw,k) = y(k)
          y(k) = y(k) + dy

        end do

        call f ( n, t, y, save1 )

        if ( n == 0 ) then
          jstate = 6
          return
        end if

        do k = j, n, mw

          y(k) = dfdy(mw,k)
          ys = max ( abs (ywt(k) ), abs ( y(k) ) )
          dy = fac(k)*ys
          if ( dy == 0.0D+00 ) dy = ys

          if ( nq == 1 ) then
            dy = sign ( dy, save2(k) )
          else
            dy = sign ( dy, yh(k,3) )
          end if

          dy = (y(k) + dy) - y(k)
          factor = -el(1,nq) * h / dy
          i1 = max ( ml+1, mw+1-k )
          i2 = min ( mw+n-k, mw+ml )
          do i = i1,i2
            i3 = k + i - mw
            dfdy(i,k) = factor*(save1(i3) - save2(i3))
          end do

          imax = max ( 1, k - mu )
          diff = abs ( save2(imax) - save1(imax) )
          i1 = imax
          i2 = min ( k + ml, n )

          do i = i1+1,i2
            if ( diff < abs ( save2(i) - save1(i) ) ) then
              imax = i
              diff = abs ( save2(i) - save1(i) )
            end if
          end do

          if ( 0.0D+00 < min ( abs ( save2(imax) ), abs ( save1(imax) ) ) ) then

            scale = max ( abs ( save2(imax) ), abs ( save1(imax) ) )

            if ( bu * scale < diff ) then
              fac(k) = max ( facmin, fac(k) * 0.1D+00 )
            else if ( br * scale <= diff .and. diff <= bl * scale ) then
              fac(k) = min ( fac(k) * 10.0D+00, facmax )
            else if ( diff < br * scale ) then
              fac(k) = min ( bp * fac(k), facmax )
            end if

          end if

        end do

      end do

      nfe = nfe + j2

    end if

    if ( iswflg == 3 ) then

      dfdymx = 0.0D+00

      do j = 1, n
        i1 = max ( ml+1, mw+1-j )
        i2 = min ( mw+n-j, mw+ml )
        do i = i1, i2
          dfdymx = max ( dfdymx, abs ( dfdy(i,j) ) )
        end do
      end do

      bnd = 0.0D+00
      if ( dfdymx /= 0.0D+00 ) then
        do j = 1,n
          i1 = max ( ml+1, mw+1-j )
          i2 = min ( mw+n-j, mw+ml )
          do i = i1,i2
            bnd = bnd + (dfdy(i,j) / dfdymx)**2
          end do
        end do
        bnd = dfdymx * sqrt ( bnd ) / ( -el(1,nq) * h )
      end if

    end if

    if ( impl == 0 ) then

      dfdy(mw,1:n) = dfdy(mw,1:n) + 1.0D+00

    else if ( impl == 1 ) then

      call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)

      if ( n == 0 ) then
        jstate = 9
        return
      end if

      do j = 1, n
        i1 = max ( ml+1, mw+1-j )
        i2 = min ( mw+n-j, mw+ml )
        dfdy(i1:i2,j) = dfdy(i1:i2,j) + a(i1:i2,j)
      end do

    else if ( impl == 2 ) then

      call fa ( n, t, y, a, matdim, ml, mu, nde )

      if ( n == 0 ) then
        jstate = 9
        return
      end if

      dfdy(mw,1:nde) =  dfdy(mw,1:nde) + a(1:nde,1)

    end if

    call dgbfa ( dfdy, matdim, n, ml, mu, ipvt, info )

    if ( info /= 0 ) then
      ier = .true.
    end if

  else if ( miter == 3 ) then

    iflag = 1
    call users ( y, yh(1,2), ywt, save1, save2, t, h, el(1,nq), impl, n, &
      nde, iflag )

    if ( n == 0 ) then
      jstate = 10
      return
    end if

  end if

  return
end
subroutine ddriv1 ( n, t, y, tout, mstate, eps, work, lenw )

!*****************************************************************************80
!
!! DDRIV1 solves a system of ordinary differential equations.
!
!  Discussion:
!
!    The system has the form:
!
!      dy(i)/dt = f(y(i),t),
!
!    given the initial conditions
!
!      y(i) = yi.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Charles Gear,
!    Numerical Initial Value Problems in Ordinary Differential Equations,
!    Prentice-Hall, 1971.
!
!  i.  choosing the correct routine
!
!     sdriv
!     ddriv
!     cdriv
!           these are the generic names for three packages for solving
!           initial value problems for ordinary differential equations.
!           sdriv uses single precision arithmetic.  ddriv uses double
!           precision arithmetic.  cdriv allows complex-valued
!           differential equations, integrated with respect to a single,
!           real, independent variable.
!
!    as an aid in selecting the proper program, the following is a
!    discussion of the important options or restrictions associated with
!    each program:
!
!      a. ddriv1 should be tried first for those routine problems with
!         no more than 200 differential equations.  internally this
!         routine has two important technical defaults:
!           1. numerical approximation of the jacobian matrix of the
!              right hand side is used.
!           2. the stiff solver option is used.
!         most users of ddriv1 should not have to concern themselves
!         with these details.
!
!      b. ddriv2 should be considered for those problems for which
!         ddriv1 is inadequate (ddriv2 has no explicit restriction on
!         the number of differential equations.)  for example, ddriv1
!         may have difficulty with problems having zero initial
!         conditions and zero derivatives.  in this case ddriv2, with an
!         appropriate value of the parameter ewt, should perform more
!         efficiently.  ddriv2 provides three important additional
!         options:
!           1. the nonstiff equation solver (as well as the stiff
!              solver) is available.
!           2. the root-finding option is available.
!           3. the program can dynamically select either the non-stiff
!              or the stiff methods.
!         internally this routine also defaults to the numerical
!         approximation of the jacobian matrix of the right hand side.
!
!      c. ddriv3 is the most flexible, and hence the most complex, of
!         the programs.  its important additional features include:
!           1. the ability to exploit band structure in the jacobian
!              matrix.
!           2. the ability to solve some implicit differential
!              equations, i.e., those having the form:
!                   a(y,t) * dy/dt = f(y,t).
!           3. the option of integrating in the one step mode.
!           4. the option of allowing the user to provide a routine
!              which computes the analytic jacobian matrix of the right
!              hand side.
!           5. the option of allowing the user to provide a routine
!              which does all the matrix algebra associated with
!              corrections to the solution components.
!
!  ii.  abstract
!
!    the function of ddriv1 is to solve n (200 or fewer) ordinary
!    differential equations of the form dy(i)/dt = f(y(i),t), given the
!    initial conditions y(i) = yi.  ddriv1 is to be called once for each
!    output point.
!
!  iii.  parameters
!
!    the user should use parameter names in the call sequence of ddriv1
!    for those quantities whose value may be altered by ddriv1.  the
!    parameters in the call sequence are:
!
!    n      = (input) the number of differential equations, n <= 200
!
!    t      = the independent variable.  on input for the first call, t
!             is the initial point.  on output, t is the point at which
!             the solution is given.
!
!    y      = the vector of dependent variables.  y is used as input on
!             the first call, to set the initial values.  on output, y
!             is the computed solution vector.  this array y is passed
!             in the call sequence of the user-provided routine f.  thus
!             parameters required by f can be stored in this array in
!             components n+1 and above.  (note: changes by the user to
!             the first n components of this array will take effect only
!             after a restart, i.e., after setting mstate to +1(-1).)
!
!    tout   = (input) the point at which the solution is desired.
!
!    mstate = an integer ( kind = 4 ) describing the status of integration.  the user
!             must initialize mstate to +1 or -1.  if mstate is
!             positive, the routine will integrate past tout and
!             interpolate the solution.  this is the most efficient
!             mode.  if mstate is negative, the routine will adjust its
!             internal step to reach tout exactly (useful if a
!             singularity exists beyond tout.)  the meaning of the
!             magnitude of mstate:
!               1  (input) means the first call to the routine.  this
!                  value must be set by the user.  on all subsequent
!                  calls the value of mstate should be tested by the
!                  user.  unless ddriv1 is to be reinitialized, only the
!                  sign of mstate may be changed by the user.  (as a
!                  convenience to the user who may wish to put out the
!                  initial conditions, ddriv1 can be called with
!                  mstate=+1(-1), and tout=t.  in this case the program
!                  will return with mstate unchanged, i.e.,
!                  mstate=+1(-1).)
!               2  (output) means a successful integration.  if a normal
!                  continuation is desired (i.e., a further integration
!                  in the same direction), simply advance tout and call
!                  again.  all other parameters are automatically set.
!               3  (output)(unsuccessful) means the integrator has taken
!                  1000 steps without reaching tout.  the user can
!                  continue the integration by simply calling ddriv1
!                  again.
!               4  (output)(unsuccessful) means too much accuracy has
!                  been requested.  eps has been increased to a value
!                  the program estimates is appropriate.  the user can
!                  continue the integration by simply calling ddriv1
!                  again.
!               5  (output)(unsuccessful) n has been set to zero in
!                subroutine f.  see description of f in section iv.
!
!    eps    = on input, the requested relative accuracy in all solution
!             components.  on output, the adjusted relative accuracy if
!             the input value was too small.  the value of eps should be
!             set as large as is reasonable, because the amount of work
!             done by ddriv1 increases as eps decreases.
!
!    work
!    lenw   = (input)
!             work is an array of lenw real words used
!             internally for temporary storage.  the user must allocate
!             space for this array in the calling program by a statement
!             such as
!                       real work(...)
!             the length of work should be at least n*n + 11*n + 225
!             and lenw should be set to the value used.  the contents of
!             work should not be disturbed between calls to ddriv1.
!
!  long description
!
!  iv.  usage
!
!                   program sample
!                   real alfa, eps, t, tout
!c                                          n is the number of equations
!                   parameter(alfa = 1.0, n = 3,
!                  8          lenw = n*n + 11*n + 225)
!                   real work(lenw), y(n+1)
!c                                                         initial point
!                   t = 0.00001
!c                                                set initial conditions
!                   y(1) = 10.0D+00
!                   y(2) = 0.0D+00
!                   y(3) = 10.0D+00
!c                                                        pass parameter
!                   y(4) = alfa
!                   tout = t
!                   mstate = 1
!                   eps = .001
!              10   call ddriv1 (n, t, y, tout, mstate, eps, work, lenw)
!                   if ( 2 < mstate ) stop
!                   write(*, '(4e12.3)') tout, y(1:3)
!                   tout = 10.0 * tout
!                   if ( tout < 50.0D+00 ) go to 10
!                 end
!
!    the user must write a subroutine called f to evaluate the right
!    hand side of the differential equations.  it is of the form:
!
!                 subroutine f (n, t, y, ydot)
!                   real alfa, t, y(*), ydot(*)
!                   alfa = y(n+1)
!                   ydot(1) = 1.0D+00 + alfa*(y(2) - y(1)) - y(1)*y(3)
!                   ydot(2) = alfa*(y(1) - y(2)) - y(2)*y(3)
!                   ydot(3) = 1.0D+00 - y(3)*(y(1) + y(2))
!                 end
!
!    this computes ydot = f(y,t), the right hand side of the
!    differential equations.  here y is a vector of length at least n.
!    the actual length of y is determined by the user's declaration in
!    the program which calls ddriv1.  thus the dimensioning of y in f,
!    while required by fortran convention, does not actually allocate
!    any storage.  when this subroutine is called, the first n
!    components of y are intermediate approximations to the solution
!    components.  the user should not alter these values.  here ydot is
!    a vector of length n.  the user should only compute ydot(i) for i
!    from 1 to n.  normally a return from f passes control back to
!    ddriv1.  however, if the user would like to abort the calculation,
!    i.e., return control to the program which calls ddriv1, he should
!    set n to zero.  ddriv1 will signal this by returning a value of
!    mstate equal to +5(-5).  altering the value of n in f has no effect
!    on the value of n in the call sequence of ddriv1.
!
!  v.  other communication to the user
!
!    a. the solver communicates to the user through the parameters
!       above.  in addition it writes diagnostic messages through the
!       standard error handling program XERROR.  that program will
!       terminate the user's run if it detects a probable problem setup
!       error, e.g., insufficient storage allocated by the user for the
!       work array.  for further information see section iii-a of the
!       writeup for ddriv3.
!
!    b. the number of evaluations of the right hand side can be found
!       in the work array in the location determined by:
!            lenw - (n + 21) + 4
!
  implicit none

  integer ( kind = 4 ), parameter :: idliw = 21
  integer ( kind = 4 ), parameter :: mxn = 200
  integer ( kind = 4 ) n

  real ( kind = 8 ) eps
  real ( kind = 8 ) ewt(1)
  external f
  real ( kind = 8 ) hmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: ierror = 2
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: impl = 0
  integer ( kind = 4 ) iwork(idliw+mxn)
  integer ( kind = 4 ) leniw
  integer ( kind = 4 ) lenw
  integer ( kind = 4 ) lenwcm
  integer ( kind = 4 ) lnwchk
  integer ( kind = 4 ), parameter :: mint = 2
  integer ( kind = 4 ), parameter :: miter = 2
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mstate
  real ( kind = 8 ) mu
  integer ( kind = 4 ), parameter :: mxord = 5
  integer ( kind = 4 ), parameter :: mxstep = 1000
  integer ( kind = 4 ) nde
  integer ( kind = 4 ), parameter :: nroot = 0
  integer ( kind = 4 ) nstate
  integer ( kind = 4 ) ntask
  real ( kind = 8 ) t
  real ( kind = 8 ) tout
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) y(n)

  ewt(1) = 1.0D+00

  if ( mxn < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV1 - Fatal error!'
    write ( *, '(a)' ) '  The number of equations is too large.'
    write ( *, '(a,i6)' ) '  The number of equations N = ', n
    write ( *, '(a,i6)' ) '  The maximum is MXN = ', mxn
    stop
  end if

  if ( 0 < mstate ) then
    nstate = mstate
    ntask = 1
  else
    nstate = - mstate
    ntask = 3
  end if

  hmax = 2.0D+00 * abs ( tout - t )
  leniw = n + idliw
  lenwcm = lenw - leniw

  if ( lenwcm < (n*n + 10*n + 204) ) then
    lnwchk = n*n + 10*n + 204 + leniw
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV1 - Fatal error!'
    write ( *, '(a)' ) '  Insufficient work storage.'
    write ( *, '(a,i6)' ) '  The given work storage is = ', lenwcm
    write ( *, '(a,i6)' ) '  The required work storage is = ', lnwchk
    stop
  end if

  if ( nstate /= 1 ) then
    do i = 1, leniw
      ii = i + lenwcm
      iwork(i) = int ( work(ii) )
    end do
  end if

  call ddriv3 (n, t, y, f, nstate, tout, ntask, nroot, eps, ewt, &
    ierror, mint, miter, impl, ml, mu, mxord, hmax, work, &
    lenwcm, iwork, leniw, f, f, nde, mxstep, f, f)

  do i = 1, leniw
    ii = lenwcm + i
    work(ii) = real ( iwork(i), kind = 8 )
  end do

  if ( nstate <= 4 ) then
    mstate = sign ( nstate, mstate )
  else if ( nstate == 6 ) then
    mstate = sign ( 5, mstate )
  end if

  return
end
subroutine ddriv2 ( n, t, y, f, tout, mstate, nroot, eps, ewt, mint, work, &
  lenw, iwork, leniw, g )

!*****************************************************************************80
!
!! DDRIV2 solves a system of ordinary differential equations.
!
!  Discussion:
!
!    DDRIV2 solves a system of N ordinary differential equations.
!
!      dy(i)/dt = f(y(i),t),
!
!    given the initial conditions
!
!      y(i) = yi.
!
!    The program has options to allow the solution of both stiff and
!    non-stiff differential equations.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  i.  abstract
!
!    the function of ddriv2 is to solve n ordinary differential
!    equations of the form dy(i)/dt = f(y(i),t), given the initial
!    conditions y(i) = yi.  the program has options to allow the
!    solution of both stiff and non-stiff differential equations.
!    ddriv2 is to be called once for each output point of t.
!
!  ii.  parameters
!
!    the user should use parameter names in the call sequence of ddriv2
!    for those quantities whose value may be altered by ddriv2.  the
!    parameters in the call sequence are:
!
!    n      = (input) the number of differential equations.
!
!    t      = the independent variable.  on input for the first call, t
!             is the initial point.  on output, t is the point at which
!             the solution is given.
!
!    y      = the vector of dependent variables.  y is used as input on
!             the first call, to set the initial values.  on output, y
!             is the computed solution vector.  this array y is passed
!             in the call sequence of the user-provided routines f and
!             g.  thus parameters required by f and g can be stored in
!             this array in components n+1 and above.  (note: changes
!             by the user to the first n components of this array will
!             take effect only after a restart, i.e., after setting
!             mstate to +1(-1).)
!
!    f      = a subroutine supplied by the user.  the name must be
!             declared external in the user's calling program.  this
!           subroutine is of the form:
!                 subroutine f (n, t, y, ydot)
!                   real y(*), ydot(*)
!                     .
!                     .
!                   ydot(1) = ...
!                     .
!                     .
!                   ydot(n) = ...
!                 end (sample)
!             this computes ydot = f(y,t), the right hand side of the
!             differential equations.  here y is a vector of length at
!             least n.  the actual length of y is determined by the
!             user's declaration in the program which calls ddriv2.
!             thus the dimensioning of y in f, while required by fortran
!             convention, does not actually allocate any storage.  when
!             this subroutine is called, the first n components of y are
!             intermediate approximations to the solution components.
!             the user should not alter these values.  here ydot is a
!             vector of length n.  the user should only compute ydot(i)
!             for i from 1 to n.  normally a return from f passes
!             control back to  ddriv2.  however, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls ddriv2, he should set n to zero.
!             ddriv2 will signal this by returning a value of mstate
!             equal to +6(-6).  altering the value of n in f has no
!             effect on the value of n in the call sequence of ddriv2.
!
!    tout   = (input) the point at which the solution is desired.
!
!    mstate = an integer ( kind = 4 ) describing the status of integration.  the user
!             must initialize mstate to +1 or -1.  if mstate is
!             positive, the routine will integrate past tout and
!             interpolate the solution.  this is the most efficient
!             mode.  if mstate is negative, the routine will adjust its
!             internal step to reach tout exactly (useful if a
!             singularity exists beyond tout.)  the meaning of the
!             magnitude of mstate:
!               1  (input) means the first call to the routine.  this
!                  value must be set by the user.  on all subsequent
!                  calls the value of mstate should be tested by the
!                  user.  unless ddriv2 is to be reinitialized, only the
!                  sign of mstate may be changed by the user.  (as a
!                  convenience to the user who may wish to put out the
!                  initial conditions, ddriv2 can be called with
!                  mstate=+1(-1), and tout=t.  in this case the program
!                  will return with mstate unchanged, i.e.,
!                  mstate=+1(-1).)
!               2  (output) means a successful integration.  if a normal
!                  continuation is desired (i.e., a further integration
!                  in the same direction), simply advance tout and call
!                  again.  all other parameters are automatically set.
!               3  (output)(unsuccessful) means the integrator has taken
!                  1000 steps without reaching tout.  the user can
!                  continue the integration by simply calling ddriv2
!                  again.  other than an error in problem setup, the
!                  most likely cause for this condition is trying to
!                  integrate a stiff set of equations with the non-stiff
!                  integrator option. (see description of mint below.)
!               4  (output)(unsuccessful) means too much accuracy has
!                  been requested.  eps has been increased to a value
!                  the program estimates is appropriate.  the user can
!                  continue the integration by simply calling ddriv2
!                  again.
!               5  (output) a root was found at a point less than tout.
!                  the user can continue the integration toward tout by
!                  simply calling ddriv2 again.
!               6  (output)(unsuccessful) n has been set to zero in
!                subroutine f.
!               7  (output)(unsuccessful) n has been set to zero in
!                function g.  see description of g below.
!
!    nroot  = (input) the number of equations whose roots are desired.
!             if nroot is zero, the root search is not active.  this
!             option is useful for obtaining output at points which are
!             not known in advance, but depend upon the solution, e.g.,
!             when some solution component takes on a specified value.
!             the root search is carried out using the user-written
!           function g (see description of g below.)  ddriv2 attempts
!             to find the value of t at which one of the equations
!             changes sign.  ddriv2 can find at most one root per
!             equation per internal integration step, and will then
!             return the solution either at tout or at a root, whichever
!             occurs first in the direction of integration.  the index
!             of the equation whose root is being reported is stored in
!             the sixth element of iwork.
!             note: nroot is never altered by this program.
!
!    eps    = on input, the requested relative accuracy in all solution
!             components.  eps = 0 is allowed.  on output, the adjusted
!             relative accuracy if the input value was too small.  the
!             value of eps should be set as large as is reasonable,
!             because the amount of work done by ddriv2 increases as
!             eps decreases.
!
!    ewt    = (input) problem zero, i.e., the smallest physically
!             meaningful value for the solution.  this is used inter-
!             nally to compute an array ywt(i) = max ( abs ( y(i) ), ewt ).
!             one step error estimates divided by ywt(i) are kept less
!             than eps.  setting ewt to zero provides pure relative
!             error control.  however, setting ewt smaller than
!             necessary can adversely affect the running time.
!
!    mint   = (input) the integration method flag.
!               mint = 1  means the adams methods, and is used for
!                         non-stiff problems.
!               mint = 2  means the stiff methods of gear (i.e., the
!                         backward differentiation formulas), and is
!                         used for stiff problems.
!               mint = 3  means the program dynamically selects the
!                         adams methods when the problem is non-stiff
!                         and the gear methods when the problem is
!                         stiff.
!             mint may not be changed without restarting, i.e., setting
!             the magnitude of mstate to 1.
!
!    work
!    lenw   = (input)
!             work is an array of lenw real words used
!             internally for temporary storage.  the user must allocate
!             space for this array in the calling program by a statement
!             such as
!                       real work(...)
!             the length of work should be at least
!               16*n + 2*nroot + 204         if mint is 1, or
!               n*n + 10*n + 2*nroot + 204   if mint is 2, or
!               n*n + 17*n + 2*nroot + 204   if mint is 3,
!             and lenw should be set to the value used.  the contents of
!             work should not be disturbed between calls to ddriv2.
!
!    iwork
!    leniw  = (input)
!             iwork is an integer ( kind = 4 ) array of length leniw used internally
!             for temporary storage.  the user must allocate space for
!             this array in the calling program by a statement such as
!                       integer ( kind = 4 ) iwork(...)
!             the length of iwork should be at least
!               21      if mint is 1, or
!               n+21    if mint is 2 or 3,
!             and leniw should be set to the value used.  the contents
!             of iwork should not be disturbed between calls to ddriv2.
!
!    g      = a real fortran function supplied by the user
!             if nroot is not 0.  in this case, the name must be
!             declared external in the user's calling program.  g is
!             repeatedly called with different values of iroot to
!             obtain the value of each of the nroot equations for which
!             a root is desired.  g is of the form:
!                   real function g (n, t, y, iroot)
!                   real y(*)
!                   go to (10, ...), iroot
!              10   g = ...
!                     .
!                     .
!                 end (sample)
!             here, y is a vector of length at least n, whose first n
!             components are the solution components at the point t.
!             the user should not alter these values.  the actual length
!             of y is determined by the user's declaration in the
!             program which calls ddriv2.  thus the dimensioning of y in
!             g, while required by fortran convention, does not actually
!             allocate any storage.  normally a return from g passes
!             control back to  ddriv2.  however, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls ddriv2, he should set n to zero.
!             ddriv2 will signal this by returning a value of mstate
!             equal to +7(-7).  in this case, the index of the equation
!             being evaluated is stored in the sixth element of iwork.
!             altering the value of n in g has no effect on the value of
!             n in the call sequence of ddriv2.
!
!  long description
!
!  iii.  other communication to the user
!
!    a. the solver communicates to the user through the parameters
!       above.  in addition it writes diagnostic messages through the
!       standard error handling program XERROR.  that program will
!       terminate the user's run if it detects a probable problem setup
!       error, e.g., insufficient storage allocated by the user for the
!       work array.  messages are written on the standard error message
!       file.  at installations which have this error handling package
!       the user should determine the standard error handling file from
!       the local documentation.  otherwise the short but serviceable
!       routine, XERROR, available with this package, can be used.
!
!    b. the first three elements of work and the first five elements of
!       iwork will contain the following statistical data:
!         avgh     the average step size used.
!         hused    the step size last used (successfully).
!         avgord   the average order used.
!         imxerr   the index of the element of the solution vector that
!                  contributed most to the last error test.
!         nqused   the order last used (successfully).
!         nstep    the number of steps taken since last initialization.
!         nfe      the number of evaluations of the right hand side.
!         nje      the number of evaluations of the jacobian matrix.
!
!  iv.  remarks
!
!    a. on any return from ddriv2 all information necessary to continue
!       the calculation is contained in the call sequence parameters,
!       including the work arrays.  thus it is possible to suspend one
!       problem, integrate another, and then return to the first.
!
!    b. if this package is to be used in an overlay situation, the user
!       must declare in the primary overlay the variables in the call
!       sequence to ddriv2.
!
!    c. when the routine g is not required, difficulties associated with
!       an unsatisfied external can be avoided by using the name of the
!       routine which calculates the right hand side of the differential
!       equations in place of g in the call sequence of ddriv2.
!
!  v.  usage
!
!               program sample
!               external f
!               parameter(mint = 1, nroot = 0, n = ...,
!              8          lenw = 16*n + 2*nroot + 204, leniw = 21)
!                                           n is the number of equations
!               real eps, ewt, t, tout, work(lenw), y(n)
!               integer ( kind = 4 ) iwork(leniw)
!               open(file='tape6', unit=6, status='new')
!               t = 0.                           initial point
!               y(1:n) = ...                     set initial conditions
!               tout = t
!               ewt = ...
!               mstate = 1
!               eps = ...
!          20   call ddriv2 (n, t, y, f, tout, mstate, nroot, eps, ewt,
!              8             mint, work, lenw, iwork, leniw, f)
!                                          last argument is not the same
!                                          as f if rootfinding is used.
!               if ( 2 < mstate ) stop
!               write(6, 100) tout, y(1:n))
!               tout = tout + 1.
!               if ( tout <= 10.0) go to 20
!          100  format(...)
!             end (sample)
!
!  Reference:
!
!    Charles Gear,
!    Numerical Initial Value Problems in Ordinary Differential Equations,
!    Prentice-Hall, 1971.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) eps
  real ( kind = 8 ) ewt
  real ( kind = 8 ) ewtcom(1)
  external f
  real ( kind = 8 ), external :: g
  real ( kind = 8 ) hmax
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), parameter :: impl = 0
  integer ( kind = 4 ) iwork(*)
  integer ( kind = 4 ) leniw
  integer ( kind = 4 ) lenw
  integer ( kind = 4 ) mint
  integer ( kind = 4 ) miter
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mstate
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) mxord
  integer ( kind = 4 ), parameter :: mxstep = 1000
  integer ( kind = 4 ) nde
  integer ( kind = 4 ) nroot
  integer ( kind = 4 ) nstate
  integer ( kind = 4 ) ntask
  real ( kind = 8 ) t
  real ( kind = 8 ) tout
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) y(n)

  if ( mint < 1 .or. 3 < mint ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV2 - Fatal error!'
    write ( *, '(a)' ) '  Improper value for the integration method.'
    write ( *, '(a,i6)' ) '  MINT = ', mint
    write ( *, '(a)' ) '  MINT should be between 1 and 3.'
    stop
  end if

  if ( 0 <= mstate ) then
    nstate = mstate
    ntask = 1
  else
    nstate = - mstate
    ntask = 3
  end if

  ewtcom(1) = ewt

  if ( ewt /= 0.0D+00 ) then
    ierror = 3
  else
    ierror = 2
  end if

  if ( mint == 1 ) then
    miter = 0
    mxord = 12
  else if ( mint == 2 ) then
    miter = 2
    mxord = 5
  else if ( mint == 3 ) then
    miter = 2
    mxord = 12
  end if

  hmax = 2.0D+00 * abs ( tout - t )

  call ddriv3 ( n, t, y, f, nstate, tout, ntask, nroot, eps, ewtcom, &
    ierror, mint, miter, impl, ml, mu, mxord, hmax, work, &
    lenw, iwork, leniw, f, f, nde, mxstep, g, f )

  if ( 0 <= mstate ) then
    mstate = nstate
  else
    mstate = - nstate
  end if

  return
end
subroutine ddriv3 ( n, t, y, f, nstate, tout, ntask, nroot, eps, ewt, ierror, &
  mint, miter, impl, ml, mu, mxord, hmax, work, lenw, iwork, leniw, jacobn, &
  fa, nde, mxstep, g, users )

!*****************************************************************************80
!
!! DDRIV3 solves a system of ordinary differential equations.
!
!  Discussion:
!
!    DDRIV3 solves a system of N ordinary differential equations
!
!      dy(i)/dt = f(y(i),t),
!
!    given the initial conditions
!
!      y(i) = yi.
!
!    The program has options to solve both stiff and non-stiff differential
!    equations.  Other important options are available.
!
!  Reference:
!
!    Charles Gear,
!    Numerical Initial Value Problems in Ordinary Differential Equations,
!    Prentice Hall, 1971.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!
!  i.  abstract
!
!    the primary function of ddriv3 is to solve n ordinary differential
!    equations of the form dy(i)/dt = f(y(i),t), given the initial
!    conditions y(i) = yi.  the program has options to allow the
!    solution of both stiff and non-stiff differential equations.  in
!    addition, ddriv3 may be used to solve:
!      1. the initial value problem, a * dy(i)/dt = f(y(i),t), where a is
!         a non-singular matrix depending on y and t.
!      2. the hybrid differential/algebraic initial value problem,
!         a * dy(i)/dt = f(y(i),t), where a is a vector (whose values may
!         depend upon y and t) some of whose components will be zero
!         corresponding to those equations which are algebraic rather
!         than differential.
!    ddriv3 is to be called once for each output point of t.
!
!  ii.  parameters
!
!    the user should use parameter names in the call sequence of ddriv3
!    for those quantities whose value may be altered by ddriv3.  the
!    parameters in the call sequence are:
!
!    n      = (input) the number of dependent functions whose solution
!             is desired.  n must not be altered during a problem.
!
!    t      = the independent variable.  on input for the first call, t
!             is the initial point.  on output, t is the point at which
!             the solution is given.
!
!    y      = the vector of dependent variables.  y is used as input on
!             the first call, to set the initial values.  on output, y
!             is the computed solution vector.  this array y is passed
!             in the call sequence of the user-provided routines f,
!             jacobn, fa, users, and g.  thus parameters required by
!             those routines can be stored in this array in components
!             n+1 and above.  (note: changes by the user to the first
!             n components of this array will take effect only after a
!             restart, i.e., after setting nstate to 1 .)
!
!    f      = a subroutine supplied by the user.  the name must be
!             declared external in the user's calling program.  this
!           subroutine is of the form:
!                 subroutine f (n, t, y, ydot)
!                   real y(*), ydot(*)
!                     .
!                     .
!                   ydot(1) = ...
!                     .
!                     .
!                   ydot(n) = ...
!                 end (sample)
!             this computes ydot = f(y,t), the right hand side of the
!             differential equations.  here y is a vector of length at
!             least n.  the actual length of y is determined by the
!             user's declaration in the program which calls ddriv3.
!             thus the dimensioning of y in f, while required by fortran
!             convention, does not actually allocate any storage.  when
!             this subroutine is called, the first n components of y are
!             intermediate approximations to the solution components.
!             the user should not alter these values.  here ydot is a
!             vector of length n.  the user should only compute ydot(i)
!             for i from 1 to n.  normally a return from f passes
!             control back to  ddriv3.  however, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls ddriv3, he should set n to zero.
!             ddriv3 will signal this by returning a value of nstate
!             equal to 6 .  altering the value of n in f has no effect
!             on the value of n in the call sequence of ddriv3.
!
!    nstate = an integer ( kind = 4 ) describing the status of integration.  the
!             meaning of nstate is as follows:
!               1  (input) means the first call to the routine.  this
!                  value must be set by the user.  on all subsequent
!                  calls the value of nstate should be tested by the
!                  user, but must not be altered.  (as a convenience to
!                  the user who may wish to put out the initial
!                  conditions, ddriv3 can be called with nstate=1, and
!                  tout=t.  in this case the program will return with
!                  nstate unchanged, i.e., nstate=1.)
!               2  (output) means a successful integration.  if a normal
!                  continuation is desired (i.e., a further integration
!                  in the same direction), simply advance tout and call
!                  again.  all other parameters are automatically set.
!               3  (output)(unsuccessful) means the integrator has taken
!                  mxstep steps without reaching tout.  the user can
!                  continue the integration by simply calling ddriv3
!                  again.
!               4  (output)(unsuccessful) means too much accuracy has
!                  been requested.  eps has been increased to a value
!                  the program estimates is appropriate.  the user can
!                  continue the integration by simply calling ddriv3
!                  again.
!               5  (output) a root was found at a point less than tout.
!                  the user can continue the integration toward tout by
!                  simply calling ddriv3 again.
!               6  (output)(unsuccessful) n has been set to zero in
!                subroutine f.
!               7  (output)(unsuccessful) n has been set to zero in
!                function g.  see description of g below.
!               8  (output)(unsuccessful) n has been set to zero in
!                subroutine jacobn.  see description of jacobn below.
!               9  (output)(unsuccessful) n has been set to zero in
!                subroutine fa.  see description of fa below.
!              10  (output)(unsuccessful) n has been set to zero in
!                subroutine users.  see description of users below.
!
!    Input, real ( kind = 8 ) TOUT, the point at which the solution is desired.  The
!    position of TOUT relative to T on the first call determines the
!    direction of integration.
!
!    ntask  = (input) an index specifying the manner of returning the
!             solution, according to the following:
!               ntask = 1  means ddriv3 will integrate past tout and
!                          interpolate the solution.  this is the most
!                          efficient mode.
!               ntask = 2  means ddriv3 will return the solution after
!                          each internal integration step, or at tout,
!                          whichever comes first.  in the latter case,
!                          the program integrates exactly to tout.
!               ntask = 3  means ddriv3 will adjust its internal step to
!                          reach tout exactly (useful if a singularity
!                          exists beyond tout.)
!
!    nroot  = (input) the number of equations whose roots are desired.
!             if nroot is zero, the root search is not active.  this
!             option is useful for obtaining output at points which are
!             not known in advance, but depend upon the solution, e.g.,
!             when some solution component takes on a specified value.
!             the root search is carried out using the user-written
!             function g (see description of g below.)  ddriv3 attempts
!             to find the value of t at which one of the equations
!             changes sign.  ddriv3 can find at most one root per
!             equation per internal integration step, and will then
!             return the solution either at tout or at a root, whichever
!             occurs first in the direction of integration.  the index
!             of the equation whose root is being reported is stored in
!             the sixth element of iwork.
!             note: nroot is never altered by this program.
!
!    eps    = on input, the requested relative accuracy in all solution
!             components.  eps = 0 is allowed.  on output, the adjusted
!             relative accuracy if the input value was too small.  the
!             value of eps should be set as large as is reasonable,
!             because the amount of work done by ddriv3 increases as eps
!             decreases.
!
!    ewt    = (input) problem zero, i.e., the smallest, nonzero,
!             physically meaningful value for the solution.  (array,
!             possibly of length one.  see following description of
!             ierror.)  setting ewt smaller than necessary can adversely
!             affect the running time.
!
!    ierror = (input) error control indicator.  a value of 3 is
!             suggested for most problems.  other choices and detailed
!             explanations of ewt and ierror are given below for those
!             who may need extra flexibility.
!
!             these last three input quantities eps, ewt and ierror
!             control the accuracy of the computed solution.  ewt and
!             ierror are used internally to compute an array ywt.  one
!             step error estimates divided by ywt(i) are kept less than
!             eps in root mean square norm.
!                 ierror (set by the user) =
!                 1  means ywt(i) = 1. (absolute error control)
!                                   ewt is ignored.
!                 2  means ywt(i) = abs ( y(i) ),  (relative error control)
!                                   ewt is ignored.
!                 3  means ywt(i) = max ( abs ( y(i) ), ewt(1) ).
!                 4  means ywt(i) = max ( abs ( y(i) ), ewt(i) ).
!                    this choice is useful when the solution components
!                    have differing scales.
!                 5  means ywt(i) = ewt(i).
!             if ierror is 3, ewt need only be dimensioned one.
!             if ierror is 4 or 5, the user must dimension ewt at least
!             n, and set its values.
!
!    mint   = (input) the integration method indicator.
!               mint = 1  means the adams methods, and is used for
!                         non-stiff problems.
!               mint = 2  means the stiff methods of gear (i.e., the
!                         backward differentiation formulas), and is
!                         used for stiff problems.
!               mint = 3  means the program dynamically selects the
!                         adams methods when the problem is non-stiff
!                         and the gear methods when the problem is
!                         stiff.  when using the adams methods, the
!                         program uses a value of miter=0; when using
!                         the gear methods, the program uses the value
!                         of miter provided by the user.  only a value
!                         of impl = 0 and a value of miter = 1, 2, 4, or
!                         5 is allowed for this option.  the user may
!                         not alter the value of mint or miter without
!                         restarting, i.e., setting nstate to 1.
!
!    miter  = (input) the iteration method indicator.
!               miter = 0  means functional iteration.  this value is
!                          suggested for non-stiff problems.
!               miter = 1  means chord method with analytic jacobian.
!                          in this case, the user supplies subroutine
!                          jacobn (see description below).
!               miter = 2  means chord method with jacobian calculated
!                          internally by finite differences.
!               miter = 3  means chord method with corrections computed
!                          by the user-written routine users.  See
!                          description of users below.  this option
!                          allows all matrix algebra and storage
!                          decisions to be made by the user.  when using
!                          a value of miter = 3, the subroutine fa is
!                          not required, even if impl is not 0.  for
!                          further information on using this option, see
!                          section iv-e below.
!               miter = 4  means the same as miter = 1 but the a and
!                          jacobian matrices are assumed to be banded.
!               miter = 5  means the same as miter = 2 but the a and
!                          jacobian matrices are assumed to be banded.
!
!    impl   = (input) the implicit method indicator.
!               impl = 0 means solving dy(i)/dt = f(y(i),t).
!               impl = 1 means solving a * dy(i)/dt = f(y(i),t),
!                        non-singular a.  See description of fa below.
!                        only mint = 1 or 2, and miter = 1, 2, 3, 4, or
!                        5 are allowed for this option.
!               impl = 2 means solving certain systems of hybrid
!                        differential/algebraic equations.  See
!                        description of fa below.  only mint = 2 and
!                        miter = 1, 2, 3, 4, or 5, are allowed for this
!                        option.
!               the value of impl must not be changed during a problem.
!
!    ml     = (input) the lower half-bandwidth in the case of a banded
!             a or jacobian matrix.  (i.e., maximum(r-c) for nonzero
!             a(r,c).)
!
!    mu     = (input) the upper half-bandwidth in the case of a banded
!             a or jacobian matrix.  (i.e., maximum(c-r).)
!
!    mxord  = (input) the maximum order desired. this is <= 12 for
!             the adams methods and <= 5 for the gear methods.  normal
!             value is 12 and 5, respectively.  if mint is 3, the
!             maximum order used will be min ( MXORD, 12 ) when using the
!             adams methods, and min ( MXORD, 5 ) when using the gear
!             methods.  mxord must not be altered during a problem.
!
!    hmax   = (input) the maximum magnitude of the step size that will
!             be used for the problem.  this is useful for ensuring that
!             important details are not missed.  if this is not the
!             case, a large value, such as the interval length, is
!             suggested.
!
!    work
!    lenw   = (input)
!             work is an array of lenw real words used
!             internally for temporary storage.  the user must allocate
!             space for this array in the calling program by a statement
!             such as
!                       real work(...)
!             the following table gives the required minimum value for
!             the length of work, depending on the value of impl and
!             miter.  lenw should be set to the value used.  the
!             contents of work should not be disturbed between calls to
!             ddriv3.
!
!      impl =   0                   1                   2
!              ---------------------------------------------------------
! miter =  0   (mxord+4)*n +       not allowed         not allowed
!              2*nroot + 204
!
!         1,2  n*n+(mxord+5)*n     2*n*n+(mxord+5)*n   n*n+(mxord+6)*n
!              + 2*nroot + 204     + 2*nroot + 204     + 2*nroot + 204
!
!          3   (mxord+4)*n +       (mxord+4)*n +       (mxord+4)*n +
!              2*nroot + 204       2*nroot + 204       2*nroot + 204
!
!         4,5  (2*ml+mu)*n +       (4*ml+2*mu)*n +     (2*ml+mu)*n +
!              (mxord+6)*n +       (mxord+7)*n +       (mxord+7)*n +
!              2*nroot + 204       2*nroot + 204       2*nroot + 204
!              ---------------------------------------------------------
!
!    iwork
!    leniw  = (input)
!             iwork is an integer ( kind = 4 ) array of length leniw used internally
!             for temporary storage.  the user must allocate space for
!             this array in the calling program by a statement such as
!                       integer ( kind = 4 ) iwork(...)
!             the length of iwork should be at least
!               21      if miter is 0 or 3, or
!               n+21    if miter is 1, 2, 4, or 5, or mint is 3,
!             and leniw should be set to the value used.  the contents
!             of iwork should not be disturbed between calls to ddriv3.
!
!    jacobn = a subroutine supplied by the user, if miter is 1 or 4.
!             if this is the case, the name must be declared external in
!             the user's calling program.  given a system of n
!             differential equations, it is meaningful to speak about
!             the partial derivative of the i-th right hand side with
!             respect to the j-th dependent variable.  in general there
!             are n*n such quantities.  often however the equations can
!             be ordered so that the i-th differential equation only
!             involves dependent variables with index near i, e.g., i+1,
!             i-2.  such a system is called banded.  if, for all i, the
!             i-th equation depends on at most the variables
!               y(i-ml), y(i-ml+1), ... , y(i), y(i+1), ... , y(i+mu)
!             then we call ml+mu+1 the bandwith of the system.  in a
!             banded system many of the partial derivatives above are
!             automatically zero.  for the cases miter = 1, 2, 4, and 5,
!             some of these partials are needed.  for the cases
!             miter = 2 and 5 the necessary derivatives are
!             approximated numerically by ddriv3, and we only ask the
!             user to tell ddriv3 the value of ml and mu if the system
!             is banded.  for the cases miter = 1 and 4 the user must
!             derive these partials algebraically and encode them in
!           subroutine jacobn.  by computing these derivatives the
!             user can often save 20-30 per cent of the computing time.
!             usually, however, the accuracy is not much affected and
!             most users will probably forego this option.  the optional
!             user-written subroutine jacobn has the form:
!                 subroutine jacobn (n, t, y, dfdy, matdim, ml, mu)
!                   real y(*), dfdy(matdim,*)
!                     .
!                     .
!                     calculate values of dfdy
!                     .
!                     .
!                 end (sample)
!             here y is a vector of length at least n.  the actual
!             length of y is determined by the user's declaration in the
!             program which calls ddriv3.  thus the dimensioning of y in
!             jacobn, while required by fortran convention, does not
!             actually allocate any storage.  when this subroutine is
!             called, the first n components of y are intermediate
!             approximations to the solution components.  the user
!             should not alter these values.  if the system is not
!             banded (miter=1), the partials of the i-th equation with
!             respect to the j-th dependent function are to be stored in
!             dfdy(i,j).  thus partials of the i-th equation are stored
!             in the i-th row of dfdy.  if the system is banded
!             (miter=4), then the partials of the i-th equation with
!             respect to y(j) are to be stored in dfdy(k,j), where
!             k=i-j+mu+1 .  normally a return from jacobn passes control
!             back to ddriv3.  however, if the user would like to abort
!             the calculation, i.e., return control to the program which
!             calls ddriv3, he should set n to zero.  ddriv3 will signal
!             this by returning a value of nstate equal to +8(-8).
!             altering the value of n in jacobn has no effect on the
!             value of n in the call sequence of ddriv3.
!
!    fa     = a subroutine supplied by the user if impl is 1 or 2, and
!             miter is not 3.  if so, the name must be declared external
!             in the user's calling program.  this subroutine computes
!             the array a, where a * dy(i)/dt = f(y(i),t).
!             there are two cases:
!
!               impl=1.
!               subroutine fa is of the form:
!                 subroutine fa (n, t, y, a, matdim, ml, mu, nde)
!                   real y(*), a(matdim,*)
!                     .
!                     .
!                     calculate all values of a
!                     .
!                     .
!                 end (sample)
!               in this case a is assumed to be a nonsingular matrix,
!               with the same structure as dfdy (see jacobn description
!               above).  programming considerations prevent complete
!               generality.  if miter is 1 or 2, a is assumed to be full
!               and the user must compute and store all values of
!               A(1:N,1:N).  If miter is 4 or 5, a is assumed
!               to be banded with lower and upper half bandwidth ml and
!               mu.  the left hand side of the i-th equation is a linear
!               combination of dy(i-ml)/dt, dy(i-ml+1)/dt, ... ,
!               dy(i)/dt, ... , dy(i+mu-1)/dt, dy(i+mu)/dt.  thus in the
!               i-th equation, the coefficient of dy(j)/dt is to be
!               stored in a(k,j), where k=i-j+mu+1.
!               note: the array a will be altered between calls to fa.
!
!               impl=2.
!               subroutine fa is of the form:
!                 subroutine fa (n, t, y, a, matdim, ml, mu, nde)
!                   real y(*), a(*)
!                     .
!                     .
!                     calculate non-zero values of a(1),...,a(nde)
!                     .
!                     .
!                 end (sample)
!               in this case it is assumed that the system is ordered by
!               the user so that the differential equations appear
!               first, and the algebraic equations appear last.  the
!               algebraic equations must be written in the form:
!               0 = f(y(i),t).  when using this option it is up to the
!               user to provide initial values for the y(i) that satisfy
!               the algebraic equations as well as possible.  it is
!               further assumed that a is a vector of length nde.  all
!               of the components of a, which may depend on t, y(i),
!               etc., must be set by the user to non-zero values.
!             here y is a vector of length at least n.  the actual
!             length of y is determined by the user's declaration in the
!             program which calls ddriv3.  thus the dimensioning of y in
!             fa, while required by fortran convention, does not
!             actually allocate any storage.  when this subroutine is
!             called, the first n components of y are intermediate
!             approximations to the solution components.  the user
!             should not alter these values.  fa is always called
!             immediately after calling f, with the same values of t
!             and y.  normally a return from fa passes control back to
!             ddriv3.  however, if the user would like to abort the
!             calculation, i.e., return control to the program which
!             calls ddriv3, he should set n to zero.  ddriv3 will signal
!             this by returning a value of nstate equal to +9(-9).
!             altering the value of n in fa has no effect on the value
!             of n in the call sequence of ddriv3.
!
!    nde    = (input) the number of differential equations.  this is
!             required only for impl = 2, with nde < n.
!
!    mxstep = (input) the maximum number of internal steps allowed on
!             one call to ddriv3.
!
!    g      = a real fortran function supplied by the user
!             if nroot is not 0.  in this case, the name must be
!             declared external in the user's calling program.  g is
!             repeatedly called with different values of iroot to obtain
!             the value of each of the nroot equations for which a root
!             is desired.  g is of the form:
!                   real function g (n, t, y, iroot)
!                   real y(*)
!                   go to (10, ...), iroot
!              10   g = ...
!                     .
!                     .
!                 end (sample)
!             here, y is a vector of length at least n, whose first n
!             components are the solution components at the point t.
!             the user should not alter these values.  the actual length
!             of y is determined by the user's declaration in the
!             program which calls ddriv3.  thus the dimensioning of y in
!             g, while required by fortran convention, does not actually
!             allocate any storage.  normally a return from g passes
!             control back to  ddriv3.  however, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls ddriv3, he should set n to zero.
!             ddriv3 will signal this by returning a value of nstate
!             equal to +7(-7).  in this case, the index of the equation
!             being evaluated is stored in the sixth element of iwork.
!             altering the value of n in g has no effect on the value of
!             n in the call sequence of ddriv3.
!
!    users  = a subroutine supplied by the user, if miter is 3.
!             if this is the case, the name must be declared external in
!             the user's calling program.  the routine users is called
!             by ddriv3 when certain linear systems must be solved.  the
!             user may choose any method to form, store and solve these
!             systems in order to obtain the solution result that is
!             returned to ddriv3.  in particular, this allows sparse
!             matrix methods to be used.  the call sequence for this
!             routine is:
!
!              subroutine users (y, yh, ywt, save1, save2, t, h, el,
!               8                  impl, n, nde, iflag)
!                real y(*), yh(*), ywt(*), save1(*),
!               8     save2(*), t, h, el
!
!             the input variable iflag indicates what action is to be
!             taken.subroutine users should perform the following
!             operations, depending on the value of iflag and impl.
!
!               iflag = 0
!                 impl = 0.  users is not called.
!                 impl = 1 or 2.  solve the system a*x = save2,
!                   returning the result in save2.  the array save1 can
!                   be used as a work array.
!
!               iflag = 1
!                 impl = 0.  compute, decompose and store the matrix
!                   (i - h * el * j), where i is the identity matrix and j
!                   is the jacobian matrix of the right hand side.  the
!                   array save1 can be used as a work array.
!                 impl = 1 or 2. compute, decompose and store the matrix
!                   (a - h * el * j).  the array save1 can be used as a work
!                   array.
!
!               iflag = 2
!                 impl = 0.   solve the system
!                     (i - h * el * j)*x = h * save2 - yh - save1,
!                   returning the result in save2.
!                 impl = 1 or 2.  solve the system
!                   (a - h * el * j)*x = h * save2 - a * (yh + save1)
!                   returning the result in save2.
!                 the array save1 should not be altered.
!             normally a return from users passes control back to
!             ddriv3.  however, if the user would like to abort the
!             calculation, i.e., return control to the program which
!             calls ddriv3, he should set n to zero.  ddriv3 will signal
!             this by returning a value of nstate equal to +10(-10).
!             altering the value of n in users has no effect on the
!             value of n in the call sequence of ddriv3.
!
!  long description
!
!  iii.  other communication to the user
!
!    a. the solver communicates to the user through the parameters
!       above.  in addition it writes diagnostic messages through the
!       standard error handling program XERROR.  that program will
!       terminate the user's run if it detects a probable problem setup
!       error, e.g., insufficient storage allocated by the user for the
!       work array.  messages are written on the standard error message
!       file.  at installations which have this error handling package
!       the user should determine the standard error handling file from
!       the local documentation.  otherwise the short but serviceable
!       routine, XERROR, available with this package, can be used.
!       following is a list of possible errors.  unless otherwise noted,
!       all messages come from ddriv3:
!
!        no.  type         message
!        ---  ----         -------
!         1   fatal        from ddriv2: the integration method flag has
!                          an illegal value.
!         2   warning      the output point is inconsistent with the
!                          value of ntask and t.
!         3   warning      number of steps to reach tout exceeds mxstep.
!         4   recoverable  requested accuracy is too stringent.
!         5   warning      step size is below the roundoff level.
!         6   fatal        eps is less than zero.
!         7   fatal        n is not positive.
!         8   fatal        insufficient work space provided.
!         9   fatal        improper value for nstate, mint, miter and/or
!                          impl.
!        10   fatal        the iwork array is too small.
!        11   fatal        the step size has gone to zero.
!        12   fatal        excessive amount of work.
!        13   fatal        for impl=1 or 2, the matrix a is singular.
!        14   fatal        mxord is not positive.
!        15   fatal        from ddriv1: n is greater than 200.
!        16   fatal        from ddriv1: the work array is too small.
!
!    b. the first three elements of work and the first five elements of
!       iwork will contain the following statistical data:
!         avgh     the average step size used.
!         hused    the step size last used (successfully).
!         avgord   the average order used.
!         imxerr   the index of the element of the solution vector that
!                  contributed most to the last error test.
!         nqused   the order last used (successfully).
!         nstep    the number of steps taken since last initialization.
!         nfe      the number of evaluations of the right hand side.
!         nje      the number of evaluations of the jacobian matrix.
!
!  iv.  remarks
!
!    a. other routines used:
!         ddntp, ddzro, ddstp, ddntl, ddpst, ddcor, ddcst,
!         ddpsc, and ddscl;
!         dgefa, dgesl, dgbfa, dgbsl, and dnrm2 (from linpack)
!         xerror (from the slatec common math library)
!       the last seven routines above, not having been written by the
!       present authors, are not explicitly part of this package.
!
!    b. on any return from ddriv3 all information necessary to continue
!       the calculation is contained in the call sequence parameters,
!       including the work arrays.  thus it is possible to suspend one
!       problem, integrate another, and then return to the first.
!
!    c. if this package is to be used in an overlay situation, the user
!       must declare in the primary overlay the variables in the call
!       sequence to ddriv3.
!
!    d. changing parameters during an integration.
!       the value of nroot, eps, ewt, ierror, mint, miter, or hmax may
!       be altered by the user between calls to ddriv3.  for example, if
!       too much accuracy has been requested (the program returns with
!       nstate = 4 and an increased value of eps) the user may wish to
!       increase eps further.  in general, prudence is necessary when
!       making changes in parameters since such changes are not
!       implemented until the next integration step, which is not
!       necessarily the next call to ddriv3.  this can happen if the
!       program has already integrated to a point which is beyond the
!       new point tout.
!
!    e. as the price for complete control of matrix algebra, the ddriv3
!       users option puts all responsibility for jacobian matrix
!       evaluation on the user.  it is often useful to approximate
!       numerically all or part of the jacobian matrix.  however this
!       must be done carefully.  the fortran sequence below illustrates
!       the method we recommend.  it can be inserted directly into
!       subroutine users to approximate jacobian elements in rows i1
!       to i2 and columns j1 to j2.
!              real dfdy(n,n), epsj, h, r,
!             8     save1(n), save2(n), t, uround, y(n), yj, ywt(n)
!              uround = epsilon ( uround )
!              epsj = sqrt ( uround )
!              do j = j1,j2
!                r = epsj * max ( abs ( ywt(j) ), abs ( y(j) ) )
!                if (r == 0.0) r = ywt(j)
!                yj = y(j)
!                y(j) = y(j) + r
!                call f (n, t, y, save1)
!                if (n == 0) return
!                y(j) = yj
!                do i = i1,i2
!                  dfdy(i,j) = (save1(i) - save2(i))/r
!                end do
!              end do
!
!       many problems give rise to structured sparse jacobians, e.g.,
!       block banded.  it is possible to approximate them with fewer
!       function evaluations than the above procedure uses; see curtis,
!       powell and reid, j. inst. maths applics, (1974), vol. 13,
!       pp. 117-119.
!
!    f. when any of the routines jacobn, fa, g, or users, is not
!       required, difficulties associated with unsatisfied externals can
!       be avoided by using the name of the routine which calculates the
!       right hand side of the differential equations in place of the
!       corresponding name in the call sequence of ddriv3.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ae
  real ( kind = 8 ) big
  logical convrg
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) eps
  real ( kind = 8 ) ewt(*)
  external f
  external fa
  real ( kind = 8 ), external :: g
  real ( kind = 8 ) glast
  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hsign
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ), parameter :: iavgh = 1
  integer ( kind = 4 ), parameter :: iavgrd = 3
  integer ( kind = 4 ) idfdy
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ignow
  integer ( kind = 4 ), parameter :: ihused = 2
  integer ( kind = 4 ), parameter :: iel = 4
  integer ( kind = 4 ), parameter :: ih = 160
  integer ( kind = 4 ), parameter :: ihmax = 161
  integer ( kind = 4 ), parameter :: ihold = 162
  integer ( kind = 4 ), parameter :: ihsign = 163
  integer ( kind = 4 ), parameter :: irc = 164
  integer ( kind = 4 ), parameter :: irmax = 165
  integer ( kind = 4 ), parameter :: it = 166
  integer ( kind = 4 ), parameter :: itout = 167
  integer ( kind = 4 ), parameter :: itq = 168
  integer ( kind = 4 ), parameter :: itrend = 204
  integer ( kind = 4 ), parameter :: iyh = 205
  integer ( kind = 4 ), parameter :: indmxr = 1
  integer ( kind = 4 ), parameter :: inqusd = 2
  integer ( kind = 4 ), parameter :: instep = 3
  integer ( kind = 4 ), parameter :: infe = 4
  integer ( kind = 4 ), parameter :: inje = 5
  integer ( kind = 4 ), parameter :: inroot = 6
  integer ( kind = 4 ), parameter :: icnvrg = 7
  integer ( kind = 4 ), parameter :: ijroot = 8
  integer ( kind = 4 ), parameter :: ijtask = 9
  integer ( kind = 4 ), parameter :: imntld = 10
  integer ( kind = 4 ), parameter :: imtrld = 11
  integer ( kind = 4 ), parameter :: inq = 12
  integer ( kind = 4 ), parameter :: inrtld = 13
  integer ( kind = 4 ), parameter :: indtrt = 14
  integer ( kind = 4 ), parameter :: inwait = 15
  integer ( kind = 4 ), parameter :: imnt = 16
  integer ( kind = 4 ), parameter :: imtrsv = 17
  integer ( kind = 4 ), parameter :: imtr = 18
  integer ( kind = 4 ), parameter :: imxrds = 19
  integer ( kind = 4 ), parameter :: imxord = 20
  integer ( kind = 4 ), parameter :: indprt = 21
  integer ( kind = 4 ), parameter :: indpvt = 22
  integer ( kind = 4 ) impl
  integer ( kind = 4 ) imxerr
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iroot
  integer ( kind = 4 ) isave1
  integer ( kind = 4 ) isave2
  integer ( kind = 4 ) itroot
  integer ( kind = 4 ) iwork(*)
  integer ( kind = 4 ) iywt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja
  external jacobn
  integer ( kind = 4 ) jaml
  integer ( kind = 4 ) jerror
  integer ( kind = 4 ) jgnow
  integer ( kind = 4 ) jhyp
  integer ( kind = 4 ) jroot
  integer ( kind = 4 ) jsave2
  integer ( kind = 4 ) jstate
  integer ( kind = 4 ) jtroot
  integer ( kind = 4 ) jyh
  integer ( kind = 4 ) jywt
  integer ( kind = 4 ) lenchk
  integer ( kind = 4 ) leniw
  integer ( kind = 4 ) lenw
  integer ( kind = 4 ) liwchk
  integer ( kind = 4 ) matdim
  integer ( kind = 4 ) maxord
  integer ( kind = 4 ) mint
  integer ( kind = 4 ) miter
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) mxord
  integer ( kind = 4 ) mxstep
  integer ( kind = 4 ) nde
  integer ( kind = 4 ) ndecom
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nroot
  real ( kind = 8 ), parameter :: nround = 20.0D+00
  integer ( kind = 4 ) nstate
  integer ( kind = 4 ) nstepl
  integer ( kind = 4 ) ntask
  real ( kind = 8 ) re
  real ( kind = 8 ) size
  real ( kind = 8 ) sum2
  real ( kind = 8 ) t
  real ( kind = 8 ) tlast
  real ( kind = 8 ) tout
  real ( kind = 8 ) troot
  real ( kind = 8 ) uround
  external users
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) y(*)

  npar = n
  uround = epsilon ( uround )

  if ( nroot /= 0 ) then
    ae = tiny ( ae )
    re = uround
  end if

  if ( eps < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV3 - Fatal error!'
    write ( *, '(a)' ) '  Improper value of EPS.'
    write ( *, '(a,g14.6)' ) '  EPS = ', eps
    write ( *, '(a)' ) '  EPS should be nonnegative.'
    stop
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV3 - Fatal error!'
    write ( *, '(a)' ) '  Improper value for the number of equations.'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a)' ) '  N should be positive.'
    stop
  end if

  if ( mxord <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV3 - Fatal error!'
    write ( *, '(a)' ) '  Improper value for the maximum order.'
    write ( *, '(a,i6)' ) '  MXORD = ', mxord
    write ( *, '(a)' ) '  MXORD should be positive.'
    stop
  end if

  if ((mint < 1 .or. 3 < mint ) .or. (mint == 3 .and. &
       (miter == 0 .or. miter == 3 .or. impl /= 0)) &
       .or. (miter < 0 .or. 5 < miter ) .or. &
       (impl /= 0 .and. impl /= 1 .and. impl /= 2) .or. &
       ((impl == 1 .or. impl == 2) .and. miter == 0) .or. &
       (impl == 2 .and. mint == 1) .or. &
       (nstate < 1 .or. 10 < nstate )) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV3 - Fatal error!'
    write ( *, '(a)' ) '  Improper value for some input quantity.'
    write ( *, '(a)' ) '  NSTATE/MSTATE/MINT/MITER/IMPL.'
    stop
  end if

  if ( miter == 0 .or. miter == 3 ) then
    liwchk = indpvt - 1
  else if ( miter == 1 .or. miter == 2 .or. miter == 4 .or. miter == 5 ) then
    liwchk = indpvt + n - 1
  end if

  if ( leniw < liwchk ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV3 - Fatal error!'
    write ( *, '(a)' ) '  Insufficient integer ( kind = 4 ) storage.'
    write ( *, '(a,i6)' ) '  LENIW = ', leniw
    write ( *, '(a,i6)' ) '  Required = ', liwchk
    stop
  end if
!
!  Allocate the work array
!  iyh is the index of yh in work.
!
  if ( mint == 1 .or. mint == 3 ) then
    maxord = min ( mxord, 12 )
  else if ( mint == 2) then
    maxord = min ( mxord, 5 )
  end if
  idfdy = iyh + (maxord + 1) * n
!
!  idfdy is the index of dfdy
!
  if (miter == 0 .or. miter == 3 ) then
    iywt = idfdy
  else if (miter == 1 .or. miter == 2 ) then
    iywt = idfdy + n*n
  else if (miter == 4 .or. miter == 5 ) then
    iywt = idfdy + (2*ml + mu + 1) * n
  end if
!
!  iywt is the index of ywt
!
  isave1 = iywt + n
!
! isave1 is the index of save1
!
  isave2 = isave1 + n
!
!  isave2 is the index of save2
!
  ignow = isave2 + n
!
!  ignow is the index of gnow
!
  itroot = ignow + nroot
!
! itroot is the index of troot
!
  ifac = itroot + nroot
!
!  ifac is the index of fac
!
  if (miter == 2 .or. miter == 5 .or. mint == 3) then
    ia = ifac + n
  else
    ia = ifac
  end if
!
!  ia is the index of a
!
  if (impl == 0 .or. miter == 3) then
    lenchk = ia - 1
  else if (impl == 1 .and. (miter == 1 .or. miter == 2)) then
    lenchk = ia - 1 + n*n
  else if (impl == 1 .and. (miter == 4 .or. miter == 5)) then
    lenchk = ia - 1 + (2*ml + mu + 1)*n
  else if (impl == 2 .and. miter /= 3) then
    lenchk = ia - 1 + n
  end if

  if (lenw < lenchk) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV3 - Fatal error!'
    write ( *, '(a)' ) '  Insufficient real storage.'
    write ( *, '(a,i6)' ) '  LENW = ', lenw
    write ( *, '(a,i6)' ) '  Required = ', lenchk
    stop
  end if

  if (miter == 0 .or. miter == 3) then
    matdim = 1
  else if (miter == 1 .or. miter == 2) then
    matdim = n
  else if (miter == 4 .or. miter == 5) then
    matdim = 2*ml + mu + 1
  end if

  if (impl == 0 .or. impl == 1) then
    ndecom = n
  else if (impl == 2) then
    ndecom = nde
  end if

  if (nstate == 1) then
!
!  initialize parameters.
!
    if (mint == 1 .or. mint == 3) then
      iwork(imxord) = min ( mxord, 12 )
    else if (mint == 2) then
      iwork(imxord) = min ( mxord, 5 )
    end if
    iwork(imxrds) = mxord
    if (mint == 1 .or. mint == 2) then
      iwork(imnt) = mint
      iwork(imtr) = miter
      iwork(imntld) = mint
      iwork(imtrld) = miter
    else if (mint == 3) then
      iwork(imnt) = 1
      iwork(imtr) = 0
      iwork(imntld) = iwork(imnt)
      iwork(imtrld) = iwork(imtr)
      iwork(imtrsv) = miter
    end if

    work(ihmax) = hmax
    h = (tout - t)*(1.0D+00 - 4.0D+00 * uround )
    h = sign ( min ( abs ( h ), hmax ), h )
    work(ih) = h
    hsign = sign ( 1.0D+00, h )
    work(ihsign) = hsign
    iwork(ijtask) = 0
    work(iavgh) = 0.0D+00
    work(ihused) = 0.0D+00
    work(iavgrd) = 0.0D+00
    iwork(indmxr) = 0
    iwork(inqusd) = 0
    iwork(instep) = 0
    iwork(infe) = 0
    iwork(inje) = 0
    iwork(inroot) = 0
    work(it) = t
    iwork(icnvrg) = 0
    iwork(indprt) = 0
!
!  Set initial conditions
!
    do i = 1,n
      jyh = i + iyh - 1
      work(jyh) = y(i)
    end do

    if ( t == tout ) then
      return
    end if

    go to 180

  end if
!
!  On a continuation, check that output points have
!  been or will be overtaken.
!
  if (iwork(icnvrg) == 1) then
    convrg = .true.
  else
    convrg = .false.
  end if
  t = work(it)
  h = work(ih)
  hsign = work(ihsign)

  if (iwork(ijtask) == 0) then
    go to 180
  end if
!
!  iwork(ijroot) flags unreported roots, and is set to the value of
!  ntask when a root was last selected.
!  it is set to zero when all roots
!  have been reported.  iwork(inroot)
!  contains the index and work(itout)
!  contains the value of the root last
!  selected to be reported.
!  iwork(inrtld) contains the value of
!  nroot and iwork(indtrt) contains
!  the value of itroot when the array
!  of roots was last calculated.
!
  if ( nroot /= 0 ) then
    jroot = iwork(ijroot)
    if ( 0 < jroot ) then
!
!  tout has just been reported.
!  if troot <= tout, report troot.
!
      if (nstate /= 5) then

       if ( work(itout)*hsign <= tout * hsign ) then
          troot = work(itout)
          call ddntp ( h, 0, n, iwork(inq), t, troot, work(iyh),  y )
          t = troot
          nstate = 5
          go to 580
        end if
!
!  A root has just been reported.
!  Select the next root.
!
      else
        troot = t
        iroot = 0

        do i = 1, iwork(inrtld)

          jtroot = iwork(indtrt) + i - 1
          if ( work(jtroot) * hsign <= troot * hsign ) then
!
!  Check for multiple roots.
!
            if (work(jtroot) == work(itout) .and. iwork(inroot) < i ) then
              iroot = i
              troot = work(jtroot)
              go to 60
            end if

            if ( work(itout) * hsign < work(jtroot) * hsign ) then
              iroot = i
              troot = work(jtroot)
            end if

          end if

        end do

 60     continue

        iwork(inroot) = iroot
        work(itout) = troot
        iwork(ijroot) = ntask

        if (ntask == 1) then
          if (iroot == 0) then
            iwork(ijroot) = 0
          else
            if ( troot * hsign <= tout * hsign ) then
              call ddntp(h, 0, n, iwork(inq), t, troot,work(iyh),y)
              nstate = 5
              t = troot
              go to 580
            end if
          end if
        else if (ntask == 2 .or. ntask == 3) then
!
!  if there are no more roots, or the
!  user has altered tout to be less
!  than a root, set ijroot to zero.
!
          if (iroot == 0 .or. (tout*hsign < troot*hsign)) then
            iwork(ijroot) = 0
          else
            call ddntp(h, 0, n, iwork(inq), t, troot, work(iyh), y)
            nstate = 5
            t = troot
            go to 580
          end if

        end if
      end if
    end if
  end if

  if (ntask == 1) then
    nstate = 2
    if ( tout * hsign <= t * hsign ) then
      call ddntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
      t = tout
      go to 580
    end if
  else if (ntask == 2) then
!
!  Check if TOUT has been reset.
!
    if ( tout * hsign < t * hsign ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DDRIV3 - Warning!'
      write ( *, '(a)' ) '  The input T was beyond TOUT.'
      write ( *, '(a)' ) '  The solution was obtained by interpolation.'
      call ddntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
      t = tout
      nstate = 2
      go to 580
    end if
!
!  Determine if TOUT has been overtaken.
!
    if ( abs ( tout - t ) <= &
      nround * uround * max ( abs ( t ), abs ( tout ) ) ) then
      t = tout
      nstate = 2
      go to 560
    end if
!
!  if there are no more roots to report, report t.
!
    if (nstate == 5) then
      nstate = 2
      go to 560
    end if

    nstate = 2
!
!  see if tout will be overtaken.
!
    if ( tout * hsign < ( t + h ) * hsign ) then
      h = tout - t
      if ( tout * hsign < ( t + h ) * hsign ) then
        h = h * ( 1.0D+00 - 4.0D+00 *uround )
      end if
      work(ih) = h
      if (h == 0.0D+00 ) then
        go to 670
      end if
      iwork(ijtask) = -1
    end if

  else if (ntask == 3) then

    nstate = 2

    if ( tout * hsign < t * hsign ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DDRIV3 - Warning!'
      write ( *, '(a)' ) '  The input T was beyond TOUT.'
      write ( *, '(a)' ) '  The solution was obtained by interpolation.'
      call ddntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
      t = tout
      go to 580
    end if

    if ( abs ( tout - t ) <= &
      nround * uround * max ( abs ( t ), abs ( tout ) ) ) then
      t = tout
      go to 560
    end if
    if ( tout * hsign < (t + h)*hsign ) then
      h = tout - t
      if ( tout * hsign < (t + h)*hsign ) then
        h = h*(1.0D+00 - 4.0D+00 *uround)
      end if
      work(ih) = h
      if (h == 0.0D+00 ) then
        go to 670
      end if
      iwork(ijtask) = -1
    end if
  end if
!
!  implement changes in mint, miter, and/or hmax.
!
  if ((mint /= iwork(imntld) .or. miter /= iwork(imtrld)) .and. &
    mint /= 3 .and. iwork(imntld) /= 3) iwork(ijtask) = -1
  if (hmax /= work(ihmax)) then
    h = sign ( min ( abs ( h ), hmax ), h )
    if (h /= work(ih)) then
      iwork(ijtask) = -1
      work(ih) = h
    end if
    work(ihmax) = hmax
  end if

 180  nstepl = iwork(instep)

  do i = 1,n
    jyh = iyh + i - 1
    y(i) = work(jyh)
  end do

  if ( nroot /= 0 ) then
    do i = 1,nroot
      jgnow = ignow + i - 1
      work(jgnow) = g (npar, t, y, i)
      if ( npar == 0 ) then
        iwork(inroot) = i
        nstate = 7
        return
      end if
    end do
  end if

  if (ierror == 1) then
    do i = 1,n
      jywt = i + iywt - 1
      work(jywt) = 1.0D+00
    end do
    go to 410
  else if (ierror == 5) then
    do i = 1,n
      jywt = i + iywt - 1
      work(jywt) = ewt(i)
    end do
    go to 410
  end if
!
!  reset ywt array.  looping point.
!
 260  continue

  if (ierror == 2) then
    do i = 1,n
      if (y(i) == 0.0D+00 ) then
        go to 290
      end if
      jywt = i + iywt - 1
      work(jywt) = abs ( y(i) )
    end do
    go to 410

 290    continue

  if (iwork(ijtask) == 0) then
      call f (npar, t, y, work(isave2))
      if (npar == 0) then
        nstate = 6
        return
      end if
      iwork(infe) = iwork(infe) + 1
      if (miter == 3 .and. impl /= 0) then
        iflag = 0
        call users(y, work(iyh), work(iywt), work(isave1), &
          work(isave2), t, h, work(iel), impl, npar, ndecom, iflag)
        if (npar == 0) then
          nstate = 10
          return
        end if
      else if (impl == 1) then
        if (miter == 1 .or. miter == 2) then
          call fa (npar, t, y, work(ia), matdim, ml, mu, ndecom)
          if (npar == 0) then
            nstate = 9
            return
          end if
          call dgefa (work(ia), matdim, n, iwork(indpvt), info)
          if (info /= 0) then
            go to 690
          end if
          call dgesl(work(ia),matdim,n,iwork(indpvt),work(isave2),0)
        else if (miter == 4 .or. miter == 5) then
          jaml = ia + ml
          call fa (npar, t, y, work(jaml), matdim, ml, mu, ndecom)
          if (npar == 0) then
            nstate = 9
            return
          end if
          call dgbfa (work(ia),matdim,n,ml,mu,iwork(indpvt),info)
          if (info /= 0) then
            go to 690
          end if
          call dgbsl (work(ia), matdim, n, ml, mu, iwork(indpvt), &
            work(isave2), 0)

        end if

      else if (impl == 2) then

        call fa (npar, t, y, work(ia), matdim, ml, mu, ndecom)
        if (npar == 0) then
          nstate = 9
          return
        end if

        do i = 1,ndecom
          ja = i + ia - 1
          jsave2 = i + isave2 - 1
          if (work(ja) == 0.0D+00 ) then
            go to 690
          end if
          work(jsave2) = work(jsave2)/work(ja)
        end do

      end if

    end if

    do j = i,n

      jywt = j + iywt - 1

      if (y(j) /= 0.0D+00 ) then
        work(jywt) = abs ( y(j) )
      else
        if (iwork(ijtask) == 0) then
          jsave2 = j + isave2 - 1
          work(jywt) = abs ( h * work(jsave2) )
        else
          jhyp = j + iyh + n - 1
          work(jywt) = abs ( work(jhyp) )
        end if
      end if

      if ( work(jywt) == 0.0D+00 ) then
        work(jywt) = uround
      end if

    end do

  else if ( ierror == 3 ) then

    do i = 1,n
      jywt = i + iywt - 1
      work(jywt) = max ( ewt(1), abs ( y(i) ) )
    end do

  else if (ierror == 4) then

    do i = 1,n
      jywt = i + iywt - 1
      work(jywt) = max ( ewt(i), abs ( y(i) ) )
    end do

  end if

 410  continue

  do i = 1,n
    jywt = i + iywt - 1
    jsave2 = i + isave2 - 1
    work(jsave2) = y(i) / work(jywt)
  end do

  sum2 = dnrm2 ( n, work(isave2), 1) / sqrt ( real ( n, kind = 8 ) )

  if ( eps < sum2 * uround ) then
    eps = sum2 * uround * ( 1.0D+00 + 10.0D+00 * uround )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV3 - Warning!'
    write ( *, '(a)' ) '  The requested accuracy EPS was not obtainable.'
    write ( *, '(a,g14.6)' ) '  EPS has been increased to ', eps
    nstate = 4
    go to 560
  end if

  if ( uround * abs ( t ) <= abs ( h ) ) then
    iwork(indprt) = 0
  else if (iwork(indprt) == 0) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DDRIV3 - Warning!'
    write ( *, '(a)' ) '  The stepsize is smaller than roundoff.'
    write ( *, '(a)' ) '  This may occur when there is an abrupt change'
    write ( *, '(a)' ) '  in the right hand side.'
    iwork(indprt) = 1
  end if

  if ( ntask /= 2 ) then
    if ( mxstep < ( iwork(instep) - nstepl ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DDRIV3 - Warning!'
      write ( *, '(a,i6)' ) '  Number of steps taken = ', mxstep
      write ( *, '(a)' ) '  TOUT not reached.'
      nstate = 3
      go to 560
    end if
  end if

  call ddstp (eps, f, fa, work(ihmax), impl, jacobn, matdim, &
    iwork(imxord), iwork(imnt), iwork(imtr), ml, mu, npar, &
    ndecom, work(iywt), uround, users,  work(iavgh), &
    work(iavgrd), work(ih), work(ihused), iwork(ijtask), &
    iwork(imntld), iwork(imtrld), iwork(infe), iwork(inje), &
    iwork(inqusd), iwork(instep), work(it), y, work(iyh), &
    work(ia), convrg, work(idfdy), work(iel), work(ifac), &
    work(ihold), iwork(indpvt), jstate, iwork(inq), &
    iwork(inwait), work(irc), work(irmax), work(isave1), &
    work(isave2), work(itq), work(itrend), mint, &
    iwork(imtrsv), iwork(imxrds))

  t = work(it)
  h = work(ih)
  go to (470, 670, 680, 690, 690, 660, 660, 660, 660, 660), jstate

 470 continue

  iwork(ijtask) = 1
!
!  determine if a root has been overtaken
!
  if (nroot /= 0) then
    iroot = 0

    do i = 1,nroot

      jtroot = itroot + i - 1
      jgnow = ignow + i - 1
      glast = work(jgnow)
      work(jgnow) = g (npar, t, y, i)

      if (npar == 0) then
        iwork(inroot) = i
        nstate = 7
        return
      end if

      if ( 0.0D+00 < glast * work(jgnow) ) then
        work(jtroot) = t + h
      else
        if (work(jgnow) == 0.0D+00 ) then
          work(jtroot) = t
          iroot = i
        else
          if (glast == 0.0D+00 ) then
            work(jtroot) = t + h
          else
            if ( uround * abs ( t ) <= abs ( work(ihused) ) ) then
              tlast = t - work(ihused)
              iroot = i
              troot = t
              call ddzro (ae, g, h, npar, iwork(inq), iroot, re, t, &
                work(iyh), uround,  troot, tlast, work(jgnow), glast,  y)
              do j = 1,n
                y(j) = work(iyh + j -1)
              end do
              if (npar == 0) then
                iwork(inroot) = i
                nstate = 7
                return
              end if
              work(jtroot) = troot
            else
              work(jtroot) = t
              iroot = i
            end if
          end if
        end if
      end if

    end do

    if (iroot == 0) then
      iwork(ijroot) = 0
!
!  select the first root
!
    else

      iwork(ijroot) = ntask
      iwork(inrtld) = nroot
      iwork(indtrt) = itroot
      troot = t + h

      do i = 1,nroot
        jtroot = itroot + i - 1
        if (work(jtroot)*hsign < troot*hsign) then
          troot = work(jtroot)
          iroot = i
        end if
      end do

      iwork(inroot) = iroot
      work(itout) = troot

      if (troot*hsign <= tout*hsign) then
        call ddntp (h, 0, n, iwork(inq), t, troot, work(iyh),  y)
        nstate = 5
        t = troot
        go to 580
      end if

    end if
  end if
!
!  Test for ntask condition to be satisfied.
!
  nstate = 2
  if (ntask == 1) then
    if (t*hsign < tout*hsign) then
      go to 260
    end if
    call ddntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
    t = tout
    go to 580
!
!  TOUT is assumed to have been attained
!  exactly if t is within twenty roundoff
!  units of tout, relative to max ( tout, t ).
!
  else if (ntask == 2) then
    if (abs ( tout - t ) <= nround * uround * &
      max ( abs ( t ), abs ( tout ) ) ) then
      t = tout
    else
      if ( tout * hsign < (t + h) * hsign ) then
        h = tout - t
        if ( tout * hsign < (t + h) * hsign ) then
          h = h*(1.0D+00 - 4.0D+00 *uround)
        end if
        work(ih) = h
        if (h == 0.0D+00 ) then
          go to 670
        end if
        iwork(ijtask) = -1
      end if
    end if

  else if ( ntask == 3 ) then

    if ( abs ( tout - t ) <= &
      nround * uround * max ( abs ( t ), abs ( tout ) ) ) then
      t = tout
    else
      if ( tout * hsign < (t + h) * hsign ) then
        h = tout - t
        if ( tout * hsign < (t + h) * hsign ) then
          h = h * (1.0D+00 - 4.0D+00 * uround)
        end if

        work(ih) = h
        if ( h == 0.0D+00 ) then
          go to 670
        end if
        iwork(ijtask) = -1
      end if
      go to 260
    end if
  end if
!
!  All returns are made through this section.  imxerr is determined.
!
 560  continue

  do i = 1,n
    jyh = i + iyh - 1
    y(i) = work(jyh)
  end do

 580 continue

  if (convrg) then
    iwork(icnvrg) = 1
  else
    iwork(icnvrg) = 0
  end if

  if (iwork(ijtask) == 0) then
    return
  end if

  big = 0.0D+00
  imxerr = 1
  iwork(indmxr) = imxerr
  do i = 1,n
!
!  size = abs ( error(i) / ywt(i) )
!
    jywt = i + iywt - 1
    jerror = i + isave1 - 1
    size = abs ( work(jerror) / work(jywt) )
    if ( big < size ) then
      big = size
      imxerr = i
      iwork(indmxr) = imxerr
    end if
  end do

  return

 660  nstate = jstate
  return

670 continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DDRIV3 - Fatal error!'
  write ( *, '(a)' ) '  The attempted stepsize has been reduced to zero.'
  write ( *, '(a)' ) '  The problem setup may be incorrect.'
  stop

680 continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DDRIV3 - Fatal error!'
  write ( *, '(a)' ) '  The stepsize has been reduced about 50 times.'
  write ( *, '(a)' ) '  The problem setup may be incorrect.'
  stop

690 continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DDRIV3 - Fatal error!'
  write ( *, '(a)' ) '  Matrix A is singular.'
  stop
end
subroutine ddscl ( hmax, n, nq, rmax, h, rc, rh, yh )

!*****************************************************************************80
!
!! DDSCL rescales the YH array whenever the ODE step size is changed.
!
!  Discussion:
!
!    DDSCL is a utility routine for the DDRIV family of ODE solvers.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nq
  real ( kind = 8 ) r1
  real ( kind = 8 ) rc
  real ( kind = 8 ) rh
  real ( kind = 8 ) rmax
  real ( kind = 8 ) yh(n,*)

  if ( h < 1.0D+00 ) then
    rh = min ( abs ( h ) * rh, abs ( h ) * rmax, hmax ) / abs ( h )
  else
    rh = min ( rh, rmax, hmax / abs ( h ) )
  end if

  r1 = 1.0D+00

  do j = 1, nq
    r1 = r1 * rh
    yh(1:n,j+1) = yh(1:n,j+1) * r1
  end do

  h = h * rh
  rc = rc * rh

  return
end
subroutine ddstp ( eps, f, fa, hmax, impl, jacobn, matdim, maxord, mint, &
  miter, ml, mu, n, nde, ywt, uround, users, avgh, avgord, h, hused, jtask, &
  mntold, mtrold, nfe, nje, nqused, nstep, t, y, yh, a, convrg, dfdy, el, &
  fac, hold, ipvt, jstate, nq, nwait, rc, rmax, save1, save2, tq, trend, &
  iswflg, mtrsv, mxrdsv )

!*****************************************************************************80
!
!! DDSTP performs one step of the integration of an ODE system.
!
!  Discussion:
!
!    DDSTP performs one step of the integration of an initial value
!    problem for a system of ordinary differential equations.
!
!  Parameters:
!
!    yh      an n by maxord+1 array containing the dependent variables
!              and their scaled derivatives.  maxord, the maximum order
!              used, is currently 12 for the adams methods and 5 for the
!              gear methods.  yh(i,j+1) contains the j-th derivative of
!              y(i), scaled by h**j/factorial(j).  only y(i),
!              1 <= i <= n, need be set by the calling program on
!              the first entry.  the yh array should not be altered by
!              the calling program.  when referencing yh as a
!              2-dimensional array, use a column length of n, as this is
!              the value used in DDSTP.
!
!    dfdy    a block of locations used for partial derivatives if miter
!              is not 0.  if miter is 1 or 2 its length must be at least
!              n*n.  if miter is 4 or 5 its length must be at least
!              (2*ml+mu+1)*n.
!
!    ywt     an array of n locations used in convergence and error tests
!
!    Workspace, real ( kind = 8 ) SAVE1(N), SAVE2(N).
!
!    ipvt    an integer ( kind = 4 ) array of length n used by the linear system
!              solvers for the storage of row interchange information.
!
!    a       a block of locations used to store the matrix a, when using
!              the implicit method.  if impl is 1, a is a matdim by n
!              array.  if miter is 1 or 2 matdim is n, and if miter is 4
!              or 5 matdim is 2*ml+mu+1.  if impl is 2 its length is n.
!
!    jtask   an integer ( kind = 4 ) used on input.
!              it has the following values and meanings:
!                 == 0  perform the first step.  this value enables
!                         the subroutine to initialize itself.
!                > 0  take a new step continuing from the last.
!                         assumes the last step was successful and
!                         user has not changed any parameters.
!                 < 0  take a new step with a new value of h and/or
!                         mint and/or miter.
!
!    jstate  a completion code with the following meanings:
!                1  the step was successful.
!                2  a solution could not be obtained with h /= 0.
!                3  a solution was not obtained in mxtry attempts.
!                4  for impl /= 0, the matrix a is singular.
!              on a return with 1 < JSTATE, the values of t and
!              the yh array are as of the beginning of the last
!              step, and h is the last step size attempted.
!
  implicit none

  integer ( kind = 4 ) matdim
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(matdim,*)
  real ( kind = 8 ) avgh
  real ( kind = 8 ) avgord
  real ( kind = 8 ), parameter :: bias1 = 1.30D+00
  real ( kind = 8 ), parameter :: bias2 = 1.20D+00
  real ( kind = 8 ), parameter :: bias3 = 1.40D+00
  real ( kind = 8 ) bnd
  logical convrg
  real ( kind = 8 ) ctest
  real ( kind = 8 ) d
  real ( kind = 8 ) denom
  real ( kind = 8 ) dfdy(matdim,*)
  real ( kind = 8 ) d1
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) el(13,12)
  real ( kind = 8 ) eps
  real ( kind = 8 ) erdn
  real ( kind = 8 ) erup
  real ( kind = 8 ) etest
  logical evalfa
  logical evaljc
  external f
  external fa
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hn
  real ( kind = 8 ) hold
  real ( kind = 8 ) hs
  real ( kind = 8 ) hused
  integer ( kind = 4 ) i
  logical, save :: ier = .false.
  integer ( kind = 4 ) impl
  integer ( kind = 4 ) ipvt(*)
  integer ( kind = 4 ) iswflg
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  external jacobn
  integer ( kind = 4 ) jstate
  integer ( kind = 4 ) jtask
  integer ( kind = 4 ) maxord
  integer ( kind = 4 ) mint
  integer ( kind = 4 ) miter
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mntold
  integer ( kind = 4 ) mtrold
  integer ( kind = 4 ) mtrsv
  integer ( kind = 4 ) mu
  integer ( kind = 4 ), parameter :: mxfail = 3
  integer ( kind = 4 ), parameter :: mxiter = 3
  integer ( kind = 4 ) mxrdsv
  integer ( kind = 4 ), parameter :: mxtry = 50
  integer ( kind = 4 ) nde
  integer ( kind = 4 ) nfail
  integer ( kind = 4 ) nfe
  integer ( kind = 4 ) nje
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nqused
  integer ( kind = 4 ) nstep
  integer ( kind = 4 ) nsv
  integer ( kind = 4 ) ntry
  real ( kind = 8 ) numer
  integer ( kind = 4 ) nwait
  real ( kind = 8 ) rc
  real ( kind = 8 ), parameter :: rctest = 0.30D+00
  real ( kind = 8 ) rh
  real ( kind = 8 ) rh1
  real ( kind = 8 ) rh2
  real ( kind = 8 ) rh3
  real ( kind = 8 ) rmax
  real ( kind = 8 ), parameter :: rmfail = 2.0D+00
  real ( kind = 8 ), parameter :: rmnorm = 10.0D+00
  real ( kind = 8 ) save1(n)
  real ( kind = 8 ) save2(n)
  logical switch
  real ( kind = 8 ) t
  real ( kind = 8 ) told
  real ( kind = 8 ) tq(3,12)
  real ( kind = 8 ) trend
  real ( kind = 8 ), parameter :: trshld = 1.0D+00
  real ( kind = 8 ) uround
  external users
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) yh(n,*)
  real ( kind = 8 ) ywt(*)
  real ( kind = 8 ) y0nrm

  nsv = n
  bnd = 0.0D+00
  switch = .false.
  ntry = 0
  told = t
  nfail = 0

  if ( jtask <= 0 ) then

    call ddntl (eps, f, fa, hmax, hold, impl, jtask, matdim, &
      maxord, mint, miter, ml, mu, n, nde, save1, t, &
      uround, users, y, ywt,  h, mntold, mtrold, nfe, rc, &
      yh,  a, convrg, el, fac, ier, ipvt, nq, nwait, rh, &
      rmax, save2, tq, trend, iswflg, jstate)

    if ( n == 0) then
      go to 440
    end if

    if ( h == 0.0D+00 ) then
      go to 400
    end if

    if ( ier ) then
      go to 420
    end if

  end if

 100  continue

  ntry = ntry + 1

  if ( mxtry < ntry ) then
    go to 410
  end if

  t = t + h
  call ddpsc (1, n, nq,  yh)
  evaljc = ( ( rctest < abs ( rc - 1.0D+00 ) ) .and. ( miter /= 0 ) )
  evalfa = .not. evaljc

 110  continue

  iter = 0

  y(1:n) = yh(1:n,1)

  call f ( n, t, y, save2 )

  if ( n == 0 ) then
    jstate = 6
    go to 430
  end if

  nfe = nfe + 1

  if ( evaljc .or. ier ) then
    call ddpst (el, f, fa, h, impl, jacobn, matdim, miter, ml, &
      mu, n, nde, nq, save2, t, users, y, yh, ywt, uround, &
      nfe, nje,  a, dfdy, fac, ier, ipvt, save1, iswflg, &
      bnd, jstate)

    if ( n == 0 ) then
      go to 430
    end if

    if ( ier ) then
      go to 160
    end if

    convrg = .false.
    rc = 1.0D+00
  end if

  save1(1:n) = 0.0D+00
!
!  Up to mxiter corrector iterations are taken.
!  convergence is tested by requiring the r.m.s.
!  norm of changes to be less than eps.  The sum of
!  the corrections is accumulated in the vector
!  save1(i).  it is approximately equal to the l-th
!  derivative of y multiplied by
!  h**l/(factorial(l-1) * el(l,nq)), and is thus
!  proportional to the actual errors to the lowest
!  power of h present (h**l).  the yh array is not
!  altered in the correction loop.  the norm of the
!  iterate difference is stored in d.  if
!  0 < iter, an estimate of the convergence rate
!  constant is stored in trend, and this is used in
!  the convergence test.
!
 130  continue

  call ddcor (dfdy, el, fa, h, impl, ipvt, matdim, miter, ml, &
        mu, n, nde, nq, t, users, y, yh, ywt,  evalfa, save1, &
        save2,  a, d, jstate)

  if ( n == 0 ) then
    go to 430
  end if

  if ( iswflg == 3 .and. mint == 1 ) then

    if (iter == 0) then

      numer = dnrm2 ( n, save1, 1)

      dfdy(1,1:n) = save1(1:n)

      y0nrm = dnrm2 ( n, yh, 1)

    else

      denom = numer

      dfdy(1,1:n) = save1(1:n) - dfdy(1,1:n)

      numer = dnrm2 ( n, dfdy, matdim)

      if ( el(1,nq) * numer <= 100.0D+00 * uround * y0nrm ) then
        if ( rmax == rmfail ) then
          switch = .true.
          go to 170
        end if
      end if

      dfdy(1,1:n) = save1(1:n)

      if ( denom /= 0.0D+00 ) then
        bnd = max ( bnd, numer / ( denom * abs ( h ) * el(1,nq) ) )
      end if

    end if
  end if

  if ( 0 < iter ) then
    trend = max ( 0.9D+00 * trend, d / d1 )
  end if

  d1 = d
  ctest = min ( 2.0D+00 * trend, 1.0D+00 ) * d

  if ( ctest <= eps ) then
    go to 170
  end if

  iter = iter + 1

  if ( iter < mxiter ) then
    do i = 1,n
      y(i) = yh(i,1) + el(1,nq) * save1(i)
    end do
    call f ( n, t, y, save2 )
    if ( n == 0 ) then
      jstate = 6
      go to 430
    end if
    nfe = nfe + 1
    go to 130
  end if
!
!  The corrector iteration failed to converge in
!  mxiter tries.  if partials are involved but are
!  not up to date, they are reevaluated for the next
!  try.  otherwise the yh array is retracted to its
!  values before prediction, and h is reduced, if
!  possible.  if not, a no-convergence exit is taken.
!
  if ( convrg ) then
    evaljc = .true.
    evalfa = .false.
    go to 110
  end if

 160 continue

  t = told
  call ddpsc (-1, n, nq,  yh)
  nwait = nq + 2

  if ( jtask /= 0 .and. jtask /= 2 ) then
    rmax = rmfail
  end if

  if ( iter == 0 ) then
    rh = 0.3D+00
  else
    rh = 0.9D+00 * ( eps / ctest )**(0.2D+00 )
  end if

  if ( rh * h == 0.0D+00 ) then
    go to 400
  end if

  call ddscl ( hmax, n, nq, rmax,  h, rc, rh, yh )
  go to 100
!
!  The corrector has converged.  convrg is set
!  to .true. if partial derivatives were used,
!  to indicate that they may need updating on
!  subsequent steps.  the error test is made.
!
 170  continue

  convrg = ( miter /= 0 )

  save2(1:nde) = save1(1:nde) / ywt(1:nde)

  etest = dnrm2 ( nde, save2, 1 ) &
    / ( tq(2,nq) * sqrt ( real ( nde, kind = 8 ) ) )
!
!  The error test failed.  nfail keeps track of
!  multiple failures.  restore t and the yh
!  array to their previous values, and prepare
!  to try the step again.  compute the optimum
!  step size for this or one lower order.
!
  if ( eps < etest ) then

    t = told
    call ddpsc (-1, n, nq,  yh)
    nfail = nfail + 1

    if ( nfail < mxfail ) then

      if ( jtask /= 0 .and. jtask /= 2 ) then
        rmax = rmfail
      end if

      rh2 = 1.0D+00 / ( bias2 * ( etest / eps ) &
        **(1.0D+00 / real ( nq + 1, kind = 8 ) ) )

      if ( 1 < nq ) then

        save2(1:nde) = yh(1:nde,nq+1) / ywt(1:nde)

        erdn = dnrm2 ( nde, save2, 1 ) &
          / ( tq(1,nq) * sqrt ( real ( nde, kind = 8 ) ) )
        rh1 = 1.0D+00 / max ( 1.0D+00, &
          bias1 * ( erdn / eps )**( 1.0D+00 / real ( nq, kind = 8 ) ) )

        if ( rh2 < rh1 ) then
          nq = nq - 1
          rc = rc * el(1,nq) / el(1,nq+1)
          rh = rh1
        else
          rh = rh2
        end if

      else

        rh = rh2

      end if
      nwait = nq + 2
      if ( rh * h == 0.0D+00 ) then
        go to 400
      end if

      call ddscl (hmax, n, nq, rmax,  h, rc, rh, yh)
      go to 100
    end if
!
!  Control reaches this section if the error test has
!  failed mxfail or more times.  It is assumed that the
!  derivatives that have accumulated in the yh array have
!  errors of the wrong order.  Hence the first derivative
!  is recomputed, the order is set to 1, and the step is retried.
!
    nfail = 0
    jtask = 2
    y(1:n) = yh(1:n,1)

    call ddntl (eps, f, fa, hmax, hold, impl, jtask, matdim, &
      maxord, mint, miter, ml, mu, n, nde, save1, t, &
      uround, users, y, ywt,  h, mntold, mtrold, nfe, rc, &
      yh,  a, convrg, el, fac, ier, ipvt, nq, nwait, rh, &
      rmax, save2, tq, trend, iswflg, jstate )

    rmax = rmnorm

    if ( n == 0) then
      go to 440
    end if

    if ( h == 0.0D+00 ) then
      go to 400
    end if

    if ( ier ) then
      go to 420
    end if

    go to 100

  end if
!
!  After a successful step, update the yh array.
!
  nstep = nstep + 1
  hused = h
  nqused = nq
  avgh = ( real ( nstep - 1, kind = 8 ) * avgh + h ) &
    / real ( nstep, kind = 8 )
  avgord = ( real ( nstep - 1, kind = 8 ) * avgord &
           + real ( nq, kind = 8 ) ) / real ( nstep, kind = 8 )

  do j = 1, nq+1
    do i = 1, n
      yh(i,j) = yh(i,j) + el(j,nq) * save1(i)
    end do
  end do

  y(1:n) = yh(1:n,1)
!
!  If iswflg is 3, consider changing integration methods.
!
  if ( iswflg == 3 ) then
    if ( bnd /= 0.0D+00 ) then
      if ( mint == 1 .and. nq <= 5 ) then
        hn = abs ( h ) / max ( uround, &
          ( etest / eps )**( 1.0D+00 / real ( nq + 1, kind = 8 ) ) )
        hn = min ( hn, 1.0D+00 / ( 2.0D+00 * el(1,nq) * bnd ) )
        hs = abs ( h ) / max ( uround, &
          ( etest / ( eps * el(nq+1,1) ) )&
          **( 1.0D+00 / real ( nq+1, kind = 8 ) ) )

        if ( 1.2D+00 * hn < hs ) then
          mint = 2
          mntold = mint
          miter = mtrsv
          mtrold = miter
          maxord = min ( mxrdsv, 5 )
          rc = 0.0D+00
          rmax = rmnorm
          trend = 1.0D+00
          call ddcst ( maxord, mint, iswflg, el, tq )
          nwait = nq + 2
        end if
      else if (mint == 2) then
        hs = abs ( h ) / max ( uround, ( etest / eps )&
          **( 1.0D+00 / real ( nq+1, kind = 8) ) )
        hn = abs ( h ) / max ( uround, &
          ( etest * el(nq+1,1) / eps)**(1.0D+00 / real ( nq+1, kind = 8 ) ) )
        hn = min ( hn, 1.0D+00 / ( 2.0D+00 * el(1,nq) * bnd ) )

        if ( hs <= hn ) then
          mint = 1
          mntold = mint
          miter = 0
          mtrold = miter
          maxord = min ( mxrdsv, 12 )
          rmax = rmnorm
          trend = 1.0D+00
          convrg = .false.
          call ddcst (maxord, mint, iswflg, el, tq)
          nwait = nq + 2
        end if
      end if
    end if
  end if

  if ( switch ) then
    mint = 2
    mntold = mint
    miter = mtrsv
    mtrold = miter
    maxord = min ( mxrdsv, 5 )
    nq = min ( nq, maxord )
    rc = 0.0D+00
    rmax = rmnorm
    trend = 1.0D+00
    call ddcst ( maxord, mint, iswflg, el, tq )
    nwait = nq + 2
  end if
!
!  Consider changing H if nwait = 1.  Otherwise
!  decrease nwait by 1.  If nwait is then 1 and
!  nq < maxord, then save1 is saved for use in
!  a possible order increase on the next step.
!
  if (jtask == 0 .or. jtask == 2) then

    rh = 1.0D+00 / max ( uround, bias2 * ( etest / eps )&
      **(1.0D+00 / real ( nq+1, kind = 8 ) ) )

    if ( trshld < rh ) then
      call ddscl ( hmax, n, nq, rmax, h, rc, rh, yh )
    end if

  else if ( 1 < nwait ) then

    nwait = nwait - 1

    if ( nwait == 1 .and. nq < maxord ) then
      do i = 1, nde
        yh(i,maxord+1) = save1(i)
      end do
    end if
!
!  If a change in H is considered, an increase or decrease in
!  order by one is considered also.  A change in H is made
!  only if it is by a factor of at least trshld.  Factors
!  rh1, rh2, and rh3 are computed, by which H could be
!  multiplied at order nq - 1, order nq, or order nq + 1,
!  respectively.  The largest of these is determined and the
!  new order chosen accordingly.  If the order is to be
!  increased, we compute one additional scaled derivative.
!  If there is a change of order, reset nq and the
!  coefficients.  In any case, H is reset according to rh and
!  the yh array is rescaled.
!
  else
    if ( nq == 1 ) then
      rh1 = 0.0D+00
    else
      do i = 1,nde
        save2(i) = yh(i,nq+1) / ywt(i)
      end do
      erdn = dnrm2 ( nde, save2, 1)/(tq(1,nq) * sqrt ( real ( nde, kind = 8 ) ) )
      rh1 = 1.0D+00 / max ( uround, &
        bias1 * ( erdn / eps )**( 1.0D+00 / real ( nq ) ) )
    end if

    rh2 = 1.0D+00 / max ( uround, bias2 * ( etest / eps )&
      **( 1.0D+00 / real ( nq + 1, kind = 8 ) ) )

    if ( nq == maxord ) then
      rh3 = 0.0D+00
    else
      do i = 1, nde
        save2(i) = ( save1(i) - yh(i,maxord+1) ) / ywt(i)
      end do
      erup = dnrm2 ( nde, save2, 1) &
        / (tq(3,nq) * sqrt ( real ( nde, kind = 8 ) ) )
      rh3 = 1.0D+00 / max ( uround, bias3 * ( erup / eps )&
        **( 1.0D+00 / real ( nq + 2, kind = 8 ) ) )
    end if

    if ( rh2 < rh1 .and. rh3 <= rh1 ) then
      rh = rh1
      if ( rh <= trshld ) then
        go to 380
      end if
      nq = nq - 1
      rc = rc * el(1,nq) / el(1,nq+1)
    else if ( rh1 <= rh2 .and. rh3 <= rh2 ) then
      rh = rh2
      if ( rh <= trshld ) then
        go to 380
      end if
    else
      rh = rh3
      if ( rh <= trshld ) then
        go to 380
      end if
      do i = 1,n
        yh(i,nq+2) = save1(i) * el(nq+1,nq) / real ( nq + 1, kind = 8 )
      end do
      nq = nq + 1
      rc = rc * el(1,nq) / el(1,nq-1)
    end if

    if ( iswflg == 3 .and. mint == 1 ) then
      if ( bnd /= 0.0D+00 ) then
        rh = min ( rh, 1.0D+00 / ( 2.0D+00 * el(1,nq) * bnd * abs ( h ) ) )
      end if
    end if

    call ddscl (hmax, n, nq, rmax,  h, rc, rh, yh)
    rmax = rmnorm
 380    nwait = nq + 2
  end if
!
!  All returns are made through this section.  H is saved
!  in HOLD to allow the caller to change H on the next step.
!
  jstate = 1
  hold = h
  return

 400  jstate = 2
  hold = h
  do i = 1,n
    y(i) = yh(i,1)
  end do
  return

 410  jstate = 3
  hold = h
  return

 420  jstate = 4
  hold = h
  return

 430  t = told
  call ddpsc (-1, nsv, nq,  yh)
  do i = 1,nsv
    y(i) = yh(i,1)
  end do

 440 continue

  hold = h

  return
end
subroutine ddzro ( ae, f, h, n, nq, iroot, re, t, yh, uround, b, c, fb, fc, y )

!*****************************************************************************80
!
!! DDZRO searches for a zero of a function in a given interval.
!
!  Discussion:
!
!    The routine searches for a zero of a function f(n, t, y, iroot)
!    between the given values B and C until the width of the
!    interval (B, C) has collapsed to within a tolerance specified
!    by the stopping criterion, abs ( b - c) <= 2 * ( rw * abs ( b ) + ae ).
!
!  Reference:
!
!    Lawrence Shampine, Herman Watts,
!    ZEROIN, a Root-Solving Routine,
!    Technical Report SC-TM-70-631,
!    Sandia National Laboratories, September 1970.
!
!    TJ Dekker,
!    Finding a Zero by Means of Successive Linear Interpolation,
!    "Constructive Aspects of the Fundamental Theorem of Algebra",
!    Edited by B Dejon and P Henrici, 1969.
!
!  Parameters:
!
!    Input, external F, the name of the routine which evaluates the function.
!    It must have the form
!
!      function f ( x )
!      real f
!      real x
!
!    Input/output, real ( kind = 8 ) B, C, the ends of the interval (B, C).
!    On output, both B and C have been reduced, to give a tighter estimate
!    of the root.  B will not necessarily be less than C.  On output, B
!    is the better estimate of the root, in the sense that the function
!    value is smaller there.
!
!    Input, real ( kind = 8 ) RE, the relative error used for RW in the
!    stopping criterion.  If the requested RE is less than machine precision,
!    then RW is set to approximately machine precision.
!
!    Input, real ( kind = 8 ) AE, the absolute error used in the stopping
!    criterion. If the given interval (B, C) contains the origin, then a
!    nonzero value should be chosen for AE.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) acbs
  real ( kind = 8 ) acmb
  real ( kind = 8 ) ae
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cmb
  real ( kind = 8 ) er
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) h
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) iroot
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) nq
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) re
  real ( kind = 8 ) rw
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) uround
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) yh(n,*)

  er = 4.0D+00 * uround
  rw = max ( re, er )
  ic = 0
  acbs = abs ( b - c )
  a = c
  fa = fc
  kount = 0
!
!  Perform interchange
!
 10 continue

  if ( abs ( fc ) < abs ( fb ) ) then
    a = b
    fa = fb
    b = c
    fb = fc
    c = a
    fc = fa
  end if

  cmb = 0.5D+00 * ( c - b )
  acmb = abs ( cmb )
  tol = rw * abs ( b ) + ae
!
!  Test stopping criterion
!
  if ( acmb <= tol ) then
    return
  end if

  if ( 50 < kount ) then
    return
  end if
!
!  Calculate new iterate implicitly as b + p/q, where we arrange 0 <= P.
!  The implicit form is used to prevent overflow.
!
  p = ( b - a ) * fb
  q = fa - fb

  if ( p < 0.0D+00 ) then
    p = -p
    q = -q
  end if
!
!  Update A and check for satisfactory reduction
!  in the size of our bounding interval.
!
  a = b
  fa = fb
  ic = ic + 1
  if ( 4 <= ic ) then
    if ( acbs <= 8.0D+00 * acmb ) then
!
!  bisect
!
      b = 0.5D+00 * ( c + b )
      go to 20
    end if
    ic = 0
  end if
  acbs = acmb
!
!  Test for too small a change.
!
  if ( p <= abs ( q ) * tol ) then
!
!  Increment by tolerance.
!
    b = b + sign ( tol, cmb )
!
!  Root ought to be between b and (c + b)/2.
!
  else if ( p < cmb * q ) then
!
!  Interpolate
!
    b = b + p / q
  else
!
!  Bisect
!
    b = 0.5D+00 * ( c + b )
  end if
!
!  Have completed computation for new iterate B.
!
 20   continue

  call ddntp ( h, 0, n, nq, t, b, yh,  y )
  fb = f ( n, b, y, iroot )

  if ( n == 0 ) then
    return
  end if

  if ( fb == 0.0D+00 ) then
    return
  end if

  kount = kount + 1
!
!  Decide whether next step is interpolation or extrapolation
!
  if ( sign ( 1.0D+00, fb ) == sign ( 1.0D+00, fc ) ) then
    c = a
    fc = fa
  end if

  go to 10
end
subroutine dfault ( n, x, typsiz, fscale, method, iexp, msg, ndigit, itnlim, &
  iagflg, iahflg, ipr, dlt, gradtl, stepmx, steptl )

!*****************************************************************************80
!
!! DFAULT sets default values for the optimization algorithm.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), an initial guess for the solution,
!    used to compute a maximum stepsize.
!
!    Output, real ( kind = 8 ) TYPSIZ(N), a typical size for each component
!    of X.
!
!    Output, real ( kind = 8 ) FSCALE, an estimate of the scale of the
!    minimization function.
!
!    Output, integer ( kind = 4 ) METHOD, specifies the algorithm to use to
!    solve the minimization problem.
!
!    Output, integer ( kind = 4 ) IEXP, set to 0 if minimization function
!    not expensive to evaluate.
!
!    Output, integer ( kind = 4 ) MSG, a message to inhibit certain automatic
!    checks and output.
!
!    Output, integer ( kind = 4 ) NDIGIT, the number of good digits in
!    minimization function.
!
!    Output, integer ( kind = 4 ) ITNLIM, the maximum number of allowable
!    iterations.
!
!    Output, integer ( kind = 4 ) IAGFLG, set to 0, meaning the analytic
!    gradient is not supplied.
!
!    Output, integer ( kind = 4 ) IAHFLG, set to 0, meaning the analytic hessian is
!    not supplied.
!
!    Output, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Output, real ( kind = 8 ) GRADTL, a tolerance at which the gradient is
!    considered close enough to zero to terminate algorithm.
!
!    Output, real ( kind = 8 ) STEPMX, the maximum stepsize, set to 0.0 to trip
!    the default maximum in OPTCHK.
!
!    Output, real ( kind = 8 ) STEPTL, a tolerance at which successive
!    iterates are considered close enough to terminate the algorithm.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dlt
  real ( kind = 8 ) epsm
  real ( kind = 8 ) fscale
  real ( kind = 8 ) gradtl
  integer ( kind = 4 ) iagflg
  integer ( kind = 4 ) iahflg
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) itnlim
  integer ( kind = 4 ) method
  integer ( kind = 4 ) msg
  integer ( kind = 4 ) ndigit
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) steptl
  real ( kind = 8 ) typsiz(n)
  real ( kind = 8 ) x(n)
!
!  Typical size of X and minimization function.
!
  typsiz(1:n) = 1.0D+00
  fscale = 1.0D+00
!
!  Tolerances.
!
  dlt = -1.0D+00
  epsm = epsilon ( epsm )
  gradtl = epsm**( 1.0D+00 / 3.0D+00 )
  stepmx = 0.0D+00
  steptl = sqrt ( epsm )
!
!  Flags.
!
  method = 1
  iexp = 1
  msg = 9
  ndigit = -1
  itnlim = 150
  iagflg = 0
  iahflg = 0
  ipr = 6

  return
end
subroutine dfftb ( n, r, wsave )

!*****************************************************************************80
!
!! DFFTB computes a real periodic sequence from its Fourier coefficients.
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    The transform is unnormalized.  A call to DFFTF followed by a call to
!    DFFTB will multiply the input sequence by N.
!
!    If N is even, the transform is defined by:
!
!      R_out(I) = R_in(1) + (-1)**(I-1) * R_in(N) + sum ( 2 <= K <= N/2 )
!
!        + 2 * R_in(2*K-2) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!        - 2 * R_in(2*K-1) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!    If N is odd, the transform is defined by:
!
!      R_out(I) = R_in(1) + sum ( 2 <= K <= (N+1)/2 )
!
!        + 2 * R_in(2*K-2) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!        - 2 * R_in(2*K-1) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real ( kind = 8 ) R(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real ( kind = 8 ) WSAVE(2*N+15), a work array.  The WSAVE array
!    must be initialized by calling DFFTI.  A different WSAVE array must be used
!    for each different value of N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) wsave(2*n+15)

  if ( n <= 1 ) then
    return
  end if

  call dfftb1 ( n, r, wsave(1), wsave(n+1), wsave(2*n+1) )

  return
end
subroutine dfftb1 ( n, c, ch, wa, ifac )

!*****************************************************************************80
!
!! DFFTB1 is a lower level routine used by DFFTB.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.
!
!    Input/output, real ( kind = 8 ) C(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real ( kind = 8 ) CH(N).
!
!    Input, real ( kind = 8 ) WA(N).
!
!    Input, integer ( kind = 4 ) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) ch(n)
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ifac(15)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(n)

  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = ido * l1

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call radb4 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call radb4 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call radb2 ( ido, l1, c, ch, wa(iw) )
      else
        call radb2 ( ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call radb3 ( ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call radb3 ( ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call radb5 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call radb5 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call radbg ( ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call radbg ( ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( ido == 1 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine dfftf ( n, r, wsave )

!*****************************************************************************80
!
!! DFFTF computes the Fourier coefficients of a real periodic sequence.
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    The transform is unnormalized.  A call to DFFTF followed by a call
!    to DFFTB will multiply the input sequence by N.
!
!    The transform is defined by:
!
!      R_out(1) = sum ( 1 <= I <= N ) R_in(I)
!
!    Letting L = (N+1)/2, then for K = 2,...,L
!
!      R_out(2*K-2) = sum ( 1 <= I <= N )
!
!        R_in(I) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!      R_out(2*K-1) = sum ( 1 <= I <= N )
!
!        -R_in(I) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!    And, if N is even, then:
!
!      R_out(N) = sum ( 1 <= I <= N ) (-1)**(I-1) * R_in(I)
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real ( kind = 8 ) R(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real ( kind = 8 ) WSAVE(2*N+15), a work array.  The WSAVE array
!    must be initialized by calling DFFTI.  A different WSAVE array must be
!    used for each different value of N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) wsave(2*n+15)

  if ( n <= 1 ) then
    return
  end if

  call dfftf1 ( n, r, wsave(1), wsave(n+1), wsave(2*n+1) )

  return
end
subroutine dfftf1 ( n, c, ch, wa, ifac )

!*****************************************************************************80
!
!! DFFTF1 is a lower level routine used by DFFTF and SINT.
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.
!
!    Input/output, real ( kind = 8 ) C(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real ( kind = 8 ) CH(N).
!
!    Input, real ( kind = 8 ) WA(N).
!
!    Input, integer ( kind = 4 ) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) ch(n)
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ifac(15)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kh
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(n)

  nf = ifac(2)
  na = 1
  l2 = n
  iw = n

  do k1 = 1, nf

    kh = nf - k1
    ip = ifac(kh+3)
    l1 = l2 / ip
    ido = n / l2
    idl1 = ido * l1
    iw = iw - ( ip - 1 ) * ido
    na = 1 - na

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call radf4 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call radf4 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call radf2 ( ido, l1, c, ch, wa(iw) )
      else
        call radf2 ( ido, l1, ch, c, wa(iw) )
      end if

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call radf3 ( ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call radf3 ( ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call radf5 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call radf5 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

    else

      if ( ido == 1 ) then
        na = 1 - na
      end if

      if ( na == 0 ) then
        call radfg ( ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
        na = 1
      else
        call radfg ( ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
        na = 0
      end if

    end if

    l2 = l1

  end do

  if ( na /= 1 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine dffti ( n, wsave )

!*****************************************************************************80
!
!! DFFTI initializes WSAVE, used in DFFTF and DFFTB.
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be transformed.
!
!    Output, real ( kind = 8 ) WSAVE(2*N+15), contains data, dependent on the
!    value of N, which is necessary for the DFFTF and DFFTB routines.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) wsave(2*n+15)

  if ( n <= 1 ) then
    return
  end if

  call dffti1 ( n, wsave(n+1), wsave(2*n+1) )

  return
end
subroutine dffti1 ( n, wa, ifac )

!*****************************************************************************80
!
!! DFFTI1 is a lower level routine used by DFFTI.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be transformed.
!
!    Input, real ( kind = 8 ) WA(N).
!
!    Input, integer ( kind = 4 ) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) argh
  real ( kind = 8 ) argld
  real ( kind = 8 ) fi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ifac(15)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) ld
  integer ( kind = 4 ) nf
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) wa(n)

  call i8_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0D+00 * pi / real ( n, kind = 8 )
  is = 0
  l1 = 1

  do k1 = 1, nf-1

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      ld = ld + l1
      i = is
      argld = real ( ld, kind = 8 ) * argh
      fi = 0.0D+00

      do ii = 3, ido, 2
        i = i + 2
        fi = fi + 1.0D+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      is = is + ido

    end do

    l1 = l2

  end do

  return
end
subroutine dgbfa ( abd, lda, n, ml, mu, ipvt, info )

!*****************************************************************************80
!
!! DGBFA factors a real band matrix by elimination.
!
!  Discussion:
!
!    DGBFA is usually called by DGBCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) ABD(LDA,N).  On input, the matrix in band
!    storage.  The columns of the matrix are stored in the columns of ABD
!    and the diagonals of the matrix are stored in rows ML+1 through
!    2*ML+MU+1 of ABD.  On output, an upper triangular matrix in band storage
!    and the multipliers which were used to obtain it.  The factorization
!    can be written A = L*U where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!    2*ML + MU + 1 <= LDA.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above the
!    main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal value.
!    K, if U(K,K) == 0.0D+00.  This is not an error condition for this
!      subroutine, but it does indicate that DGBSL will divide by zero if
!      called.  Use RCOND in DGBCO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mu
  real ( kind = 8 ) t

  m = ml + mu + 1
  info = 0
!
!  Zero initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    do i = i0, ml
      abd(i,jz) = 0.0D+00
    end do
  end do

  jz = j1
  ju = 0
!
!  Gaussian elimination with partial pivoting.
!
  do k = 1, n-1
!
!  Zero next fill-in column.
!
     jz = jz + 1
     if ( jz <= n ) then
       abd(1:ml,jz) = 0.0D+00
     end if
!
!  Find L = pivot index.
!
     lm = min ( ml, n-k )
     l = idamax ( lm+1, abd(m,k), 1 ) + m - 1
     ipvt(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
     if ( abd(l,k) == 0.0D+00 ) then

       info = k
!
!  Interchange if necessary.
!
     else

        if ( l /= m ) then
           t = abd(l,k)
           abd(l,k) = abd(m,k)
           abd(m,k) = t
        end if
!
!  Compute multipliers.
!
        t = -1.0D+00 / abd(m,k)
        call dscal ( lm, t, abd(m+1,k), 1 )
!
!  Row elimination with column indexing.
!
        ju = min ( max ( ju, mu+ipvt(k) ), n )
        mm = m

        do j = k+1, ju
           l = l - 1
           mm = mm - 1
           t = abd(l,j)
           if ( l /= mm ) then
              abd(l,j) = abd(mm,j)
              abd(mm,j) = t
           end if
           call daxpy ( lm, t, abd(m+1,k), 1, abd(mm+1,j), 1 )
        end do

     end if

  end do

  ipvt(n) = n

  if ( abd(m,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine dgbsl ( abd, lda, n, ml, mu, ipvt, b, job )

!*****************************************************************************80
!
!! DGBSL solves a real banded system factored by DGBCO or DGBFA.
!
!  Discussion:
!
!    DGBSL can solve either A * X = B  or  A' * X = B.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if DGBCO has set 0.0 < RCOND
!    or DGBFA has set INFO == 0.
!
!    To compute inverse(A) * C  where C is a matrix with P columns:
!
!      call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
!
!      if ( rcond is too small ) then
!        exit
!      end if
!
!      do j = 1, p
!        call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
!      end do
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ABD(LDA,N), the output from DGBCO or DGBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above the
!    main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGBCO or DGBFA.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Input, integer ( kind = 4 ) JOB, job choice.
!    0, solve A*X=B.
!    nonzero, solve A'*X=B.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t

  m = mu + ml + 1
!
!  JOB = 0, Solve  a * x = b.
!
!  First solve l*y = b.
!
  if ( job == 0 ) then

     if ( 0 < ml ) then

        do k = 1, n-1
           lm = min ( ml, n-k )
           l = ipvt(k)
           t = b(l)
           if ( l /= k ) then
              b(l) = b(k)
              b(k) = t
           end if
           call daxpy ( lm, t, abd(m+1,k), 1, b(k+1), 1 )
        end do

     end if
!
!  Now solve U * x = y.
!
     do k = n, 1, -1
        b(k) = b(k) / abd(m,k)
        lm = min ( k, m ) - 1
        la = m - lm
        lb = k - lm
        t = -b(k)
        call daxpy ( lm, t, abd(la,k), 1, b(lb), 1 )
     end do
!
!  JOB nonzero, solve A' * X = B.
!
!  First solve U'*Y = B.
!
  else

     do k = 1, n
        lm = min ( k, m ) - 1
        la = m - lm
        lb = k - lm
        t = ddot ( lm, abd(la,k), 1, b(lb), 1 )
        b(k) = ( b(k) - t ) / abd(m,k)
     end do
!
!  Now solve L'*X = Y.
!
     if ( 0 < ml ) then

        do k = n-1, 1, -1
           lm = min ( ml, n-k )
           b(k) = b(k) + ddot ( lm, abd(m+1,k), 1, b(k+1), 1 )
           l = ipvt(k)
           if ( l /= k ) then
              t = b(l)
              b(l) = b(k)
              b(k) = t
           end if
        end do

      end if

  end if

  return
end
subroutine dgeco ( a, lda, n, ipvt, rcond, z )

!*****************************************************************************80
!
!! DGECO factors a real matrix and estimates its condition number.
!
!  Discussion:
!
!    If RCOND is not needed, DGEFA is slightly faster.
!
!    To solve A * X = B, follow DGECO by DGESL.
!
!    To compute inverse ( A ) * C, follow DGECO by DGESL.
!
!    To compute determinant ( A ), follow DGECO by DGEDI.
!
!    To compute inverse ( A ), follow DGECO by DGEDI.
!
!    For the system A * X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    EPSILON/RCOND.
!
!    If RCOND is so small that the logical expression
!      1.0D+00 + RCOND == 1.0D+00
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate
!    underflows.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Author:
!
!    Cleve Moler,
!    University of New Mexico / Argonne National Lab.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, a matrix to be
!    factored.  On output, the LU factorization of the matrix.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal condition
!    number of A.
!
!    Output, real ( kind = 8 ) Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) anorm
  real ( kind = 8 ) ek
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) sm
  real ( kind = 8 ) t
  real ( kind = 8 ) wk
  real ( kind = 8 ) wkm
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Compute the L1 norm of A.
!
  anorm = 0.0D+00
  do j = 1, n
    anorm = max ( anorm, sum ( abs ( a(1:n,j) ) ) )
  end do
!
!  Compute the LU factorization.
!
  call dgefa ( a, lda, n, ipvt, info )
!
!  RCOND = 1 / ( norm(A) * (estimate of norm(inverse(A))) )
!
!  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
!
!  where
!    A * Z = Y
!  and
!    A' * Y = E
!
!  The components of E are chosen to cause maximum local growth in the
!  elements of W, where U'*W = E.  The vectors are frequently rescaled
!  to avoid overflow.
!
!  Solve U' * W = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00

  do k = 1, n

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, -z(k) )
    end if

    if ( abs ( a(k,k) ) < abs ( ek - z(k) ) ) then
      s = abs ( a(k,k) ) / abs ( ek - z(k) )
      z(1:n) = s * z(1:n)
      ek = s * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = abs ( wk )
    sm = abs ( wkm )

    if ( a(k,k) /= 0.0D+00 ) then
      wk = wk / a(k,k)
      wkm = wkm / a(k,k)
    else
      wk = 1.0D+00
      wkm = 1.0D+00
    end if

    if ( k+1 <= n ) then

      do j = k+1, n
        sm = sm + abs ( z(j) + wkm * a(k,j) )
        z(j) = z(j) + wk * a(k,j)
        s = s + abs ( z(j) )
      end do

      if ( s < sm ) then
        t = wkm - wk
        wk = wkm
        z(k+1:n) = z(k+1:n) + t * a(k,k+1:n)
      end if

    end if

    z(k) = wk

  end do

  z(1:n) = z(1:n) / sum ( abs ( z(1:n) ) )
!
!  Solve L' * Y = W
!
  do k = n, 1, -1

    z(k) = z(k) + dot_product ( a(k+1:n,k), z(k+1:n) )

    if ( 1.0D+00 < abs ( z(k) ) ) then
      z(1:n) = z(1:n) / abs ( z(k) )
    end if

    l = ipvt(k)

    t = z(l)
    z(l) = z(k)
    z(k) = t

  end do

  z(1:n) = z(1:n) / sum ( abs ( z(1:n) ) )

  ynorm = 1.0D+00
!
!  Solve L * V = Y.
!
  do k = 1, n

    l = ipvt(k)

    t = z(l)
    z(l) = z(k)
    z(k) = t

    z(k+1:n) = z(k+1:n) + t * a(k+1:n,k)

    if ( 1.0D+00 < abs ( z(k) ) ) then
      ynorm = ynorm / abs ( z(k) )
      z(1:n) = z(1:n) / abs ( z(k) )
    end if

  end do

  s = sum ( abs ( z(1:n) ) )
  z(1:n) = z(1:n) / s
  ynorm = ynorm / s
!
!  Solve U * Z = V.
!
  do k = n, 1, -1

    if ( abs ( a(k,k) ) < abs ( z(k) ) ) then
      s = abs ( a(k,k) ) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    if ( a(k,k) /= 0.0D+00 ) then
      z(k) = z(k) / a(k,k)
    else
      z(k) = 1.0D+00
    end if

    z(1:k-1) = z(1:k-1) - z(k) * a(1:k-1,k)

  end do
!
!  Normalize Z in the L1 norm.
!
  s = 1.0D+00 / sum ( abs ( z(1:n) ) )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! DGEFA factors a real matrix.
!
!  Modified:
!
!    07 March 2001
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On intput, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers used to obtain
!    it.  The factorization can be written A=L*U, where L is a product of
!    permutation and unit lower triangular matrices, and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity indicator.
!    0, normal value.
!    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
!    but it does indicate that DGESL or DGEDI will divide by zero if called.
!    Use RCOND in DGECO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find L = pivot index.
!
    l = idamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      cycle
    end if
!
!  Interchange if necessary.
!
    if ( l /= k ) then
      t = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Compute multipliers.
!
    t = -1.0D+00 / a(k,k)
    call dscal ( n-k, t, a(k+1,k), 1 )
!
!  Row elimination with column indexing.
!
    do j = k+1, n
      t = a(l,j)
      if ( l /= k ) then
        a(l,j) = a(k,j)
        a(k,j) = t
      end if
      call daxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
    end do

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine dgefs ( a, lda, n, v, itask, ind, work, iwork, rcond )

!*****************************************************************************80
!
!! DGEFS solves a general N by N system of single precision linear equations.
!
!  Discussion:
!
!    DGEFS uses the LINPACK subroutines DGECO and DGESL.  That is, if A is
!    an N by N real matrix and if X and B are real N vectors, then DGEFS
!    solves the equation
!
!      A * X = B.
!
!    The matrix A is first factored into upper and lower triangular
!    matrices U and L using partial pivoting.  These factors and the
!    pivoting information are used to find the solution vector X.
!    An approximate condition number is calculated to provide a rough
!    estimate of the number of digits of accuracy in the computed solution.
!
!    If the equation A*X=B is to be solved for more than one vector
!    B, the factoring of A does not need to be performed again and
!    the option to only solve (ITASK == 2) will be faster for
!    the succeeding solutions.  In this case, the contents of A,
!    LDA, N and IWORK must not have been altered by the user following
!    factorization (ITASK=1).  IND will not be changed by DGEFS
!    in this case.  Other settings of ITASK are used to solve linear
!    systems involving the transpose of A.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On input, the coefficient matrix.
!    On output, an upper triangular matrix U and the multipliers necessary to
!    construct a matrix L so that A=L*U.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.  The first N elements of
!    the array A are the elements of the first column of the matrix A.
!    N must be greater than or equal to 1.  (Terminal error message IND=-2)
!
!    Input/output, real ( kind = 8 ) V(N).
!    On entry, the right hand side B of a system of simultaneous linear
!    equations A*X=B.
!    On output, V contains the solution vector, X.
!
!    Input, integer ( kind = 4 ) ITASK, indicates the task to carry out.
!    1, the matrix A is factored and the linear equation is solved.
!    2, the equation is solved using the existing factored matrix A and IWORK.
!    3, the matrix is factored and A'*X=b is solved
!    4, the transposed equation is solved using the existing factored matrix
!       A and IWORK.
!
!    Output, integer ( kind = 4 ) IND, accuracy estimate and error flag.
!             gt. 0  ind is a rough estimate of the number of digits
!                     of accuracy in the solution, x.
!             lt. 0  see error message corresponding to ind below.
!    ind=-1  fatal   n is greater than lda.
!    ind=-2  fatal   n is less than 1.
!    ind=-3  fatal   itask is less than 1 or greater than 4.
!    ind=-4  fatal   the matrix a is computationally singular.
!                         a solution has not been computed.
!    ind=-10 warning    the solution has no apparent significance.
!                         the solution may be inaccurate or the matrix
!                         a may be poorly scaled.
!
!    Workspace, real ( kind = 8 ) WORK(N).
!
!    Workspace, integer ( kind = 4 ) IWORK(N).
!
!    Output, real ( kind = 8 ) RCOND, estimate of 1 / condition_number(A).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) itask
  integer ( kind = 4 ) iwork(*)
  integer ( kind = 4 ) job
  real ( kind = 8 ) rcond
  real ( kind = 8 ) v(*)
  real ( kind = 8 ) work(*)

  if ( lda < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGEFS - Fatal error!'
    write ( *, '(a)' ) '  LDA < N.'
    return
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGEFS - Fatal error!'
    write ( *, '(a)' ) '  N <= 0.'
    return
  end if

  if ( itask < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGEFS - Fatal error!'
    write ( *, '(a)' ) '  ITASK < 1.'
    return
  end if

  if ( 4 < itask ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGEFS - Fatal error!'
    write ( *, '(a)' ) '  4 < ITASK.'
    return
  end if
!
!  Factor the matrix.
!
  if ( itask == 1 .or. itask == 3 ) then

    call dgeco ( a, lda, n, iwork, rcond, work )
!
!  Check for computational singularity.
!
    if ( rcond == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGEFS - Error!'
      write ( *, '(a)' ) '  The matrix is compuationally singular.'
      return
    end if
!
!  Estimate the number of significant digits.
!
    ind = - int ( log10 ( epsilon ( rcond ) / rcond ) )

    if ( ind <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGEFS - Warning!'
      write ( *, '(a)' ) '  Solution may have no significant digits.'
    end if

  end if

  if ( itask <= 2 ) then
    job = 0
  else
    job = 1
  end if

  call dgesl ( a, lda, n, iwork, v, job )

  return
end
subroutine dgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! DGESL solves a real general linear system A * X = B.
!
!  Discussion:
!
!    DGESL can solve either of the systems A * X = B or A' * X = B.
!
!    The system matrix must have been factored by DGECO or DGEFA.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if DGECO has set 0.0D+00 < RCOND
!    or DGEFA has set INFO == 0.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Modified:
!
!    07 March 2001
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(LDA,N), the output from DGECO or DGEFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGECO or DGEFA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, solve A * X = B;
!    nonzero, solve A' * X = B.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
!
!  Solve A * X = B.
!
  if ( job == 0 ) then

    do k = 1, n-1

      l = ipvt(k)
      t = b(l)

      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if

      call daxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      call daxpy ( k-1, t, a(1,k), 1, b(1), 1 )
    end do

  else
!
!  Solve A' * X = B.
!
    do k = 1, n
      t = ddot ( k-1, a(1,k), 1, b(1), 1 )
      b(k) = ( b(k) - t ) / a(k,k)
    end do

    do k = n-1, 1, -1

      b(k) = b(k) + ddot ( n-k, a(k+1,k), 1, b(k+1), 1 )
      l = ipvt(k)

      if ( l /= k ) then
        t = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
function dnor ( )

!*****************************************************************************80
!
!! DNOR generates normal random numbers.
!
!  Discussion:
!
!    DNOR generates normal random numbers with zero mean and
!    unit standard deviation, often denoted n(0,1).
!
!    Before the first call to DNOR, you should call DSTART, passing it
!    a nonzero value of ISEED.  This will initialize the routine.
!
!  Modified:
!
!    21 April 2007
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    George Marsaglia, Wai Wan Tsang,
!    A fast, easily implemented method for sampling from decreasing or
!    symmetric unimodal density functions,
!    SIAM Journal of Scientific and Statistical Computing,
!    Volume 5, 1983, pages 349-359.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) DNOR, a normal random number.
!
  implicit none

  real ( kind = 8 ), parameter :: aa = 12.37586D+00
  real ( kind = 8 ), parameter :: b = 0.4878992D+00
  real ( kind = 8 ), parameter :: c = 12.67706D+00
  real ( kind = 8 ), save :: c1 = 0.9689279D+00
  real ( kind = 8 ), save :: c2 = 1.301198D+00
  real ( kind = 8 ) dnor
  real ( kind = 8 ) dstart
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) id
  integer ( kind = 4 ), save :: ii = 17
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) iseed
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: jj = 5
  integer ( kind = 4 ) jjj
  real ( kind = 8 ), save :: pc = 0.01958303D+00
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ), save, dimension ( 17 ) :: u = (/ &
     0.8668672834288D+00,  0.3697986366357D+00,  0.8008968294805D+00, &
     0.4173889774680D+00,  0.8254561579836D+00,  0.9640965269077D+00, &
     0.4508667414265D+00,  0.6451309529668D+00,  0.1645456024730D+00, &
     0.2787901807898D+00,  0.06761531340295D+00, 0.9663226330820D+00, &
     0.01963343943798D+00, 0.02947398211399D+00, 0.1636231515294D+00, &
     0.3976343250467D+00,  0.2631008574685D+00 /)
  real ( kind = 8 ) un
  real ( kind = 8 ), dimension ( 65 ) :: v = (/ &
    0.3409450D+00, 0.4573146D+00, 0.5397793D+00, 0.6062427D+00, 0.6631691D+00, &
    0.7136975D+00, 0.7596125D+00, 0.8020356D+00, 0.8417227D+00, 0.8792102D+00, &
    0.9148948D+00, 0.9490791D+00, 0.9820005D+00, 1.0138492D+00, 1.0447810D+00, &
    1.0749254D+00, 1.1043917D+00, 1.1332738D+00, 1.1616530D+00, 1.1896010D+00, &
    1.2171815D+00, 1.2444516D+00, 1.2714635D+00, 1.2982650D+00, 1.3249008D+00, &
    1.3514125D+00, 1.3778399D+00, 1.4042211D+00, 1.4305929D+00, 1.4569915D+00, &
    1.4834526D+00, 1.5100121D+00, 1.5367061D+00, 1.5635712D+00, 1.5906454D+00, &
    1.6179680D+00, 1.6455802D+00, 1.6735255D+00, 1.7018503D+00, 1.7306045D+00, &
    1.7598422D+00, 1.7896223D+00, 1.8200099D+00, 1.8510770D+00, 1.8829044D+00, &
    1.9155830D+00, 1.9492166D+00, 1.9839239D+00, 2.0198430D+00, 2.0571356D+00, &
    2.0959930D+00, 2.1366450D+00, 2.1793713D+00, 2.2245175D+00, 2.2725185D+00, &
    2.3239338D+00, 2.3795007D+00, 2.4402218D+00, 2.5075117D+00, 2.5834658D+00, &
    2.6713916D+00, 2.7769943D+00, 2.7769943D+00, 2.7769943D+00, 2.7769943D+00 /)
  real ( kind = 8 ) vni
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xn = 2.776994D+00
  real ( kind = 8 ) y
!
!  fast part...
!
!  Basic generator is Fibonacci.
!
  un = u(ii) - u(jj)
  if ( un < 0.0D+00 ) then
    un = un + 1.0D+00
  end if

  u(ii) = un
!
!  u(ii) and un are uniform on [0,1)
!  vni is uniform on [-1,1)
!
  vni = un + un - 1.0D+00
  ii = ii-1
  if ( ii == 0 ) ii = 17
  jj = jj-1
  if ( jj == 0 ) jj = 17
!
!  int ( un(ii) * 128 ) in range [0,127],  j is in range [1,64]
!
  j = mod ( int ( u(ii) * 128 ), 64 ) + 1
!
!  Pick sign as VNI is positive or negative.
!
  dnor = vni * v(j+1)

  if ( abs ( dnor ) <= v(j) ) then
    return
  end if
!
!  slow part; aa is a * f(0)
!
  x = ( abs ( dnor ) - v(j) ) / ( v(j+1) - v(j) )
!
!  Y is uniform on [0,1)
!
  y = u(ii) - u(jj)
  if ( y < 0.0D+00 ) then
    y = y + 1.0D+00
  end if

  u(ii) = y

  ii = ii-1
  if ( ii == 0 ) then
    ii = 17
  end if

  jj = jj-1
  if ( jj == 0 ) then
    jj = 17
  end if

  s = x + y

  if ( c2 < s ) then
    dnor = sign ( b - b * x, dnor )
    return
  end if

  if ( s <= c1 ) then
    return
  end if

  if ( c - aa * exp ( -0.5D+00 * ( b - b * x )**2 ) < y ) then
    dnor = sign ( b - b * x, dnor )
    return
  end if

  if ( exp ( -0.5D+00 * v(j+1)**2 ) + y * pc / v(j+1) <= &
    exp ( -0.5D+00 * dnor**2 ) ) then
    return
  end if
!
!  tail part; 0.3601016 is 1.0/xn
!  y is uniform on [0,1)
!
  do

    y = u(ii) - u(jj)
    if ( y <= 0.0D+00 ) then
      y = y + 1.0D+00
    end if

    u(ii) = y

    ii = ii - 1
    if ( ii == 0 ) then
      ii = 17
    end if

    jj = jj - 1
    if ( jj == 0 ) then
      jj = 17
    end if

    x = 0.3601016D+00 * log ( y )
!
!  Y is uniform on [0,1).
!
    y = u(ii) - u(jj)
    if ( y <= 0.0D+00 ) then
      y = y + 1.0D+00
    end if

    u(ii) = y

    ii = ii - 1
    if ( ii == 0 ) then
      ii = 17
    end if

    jj = jj - 1
    if ( jj == 0 ) then
      jj = 17
    end if

    if ( x * x < -2.0D+00 * log ( y ) ) then
      dnor = sign ( xn - x, dnor )
      return
    end if

  end do
!
!  fill
!
entry dstart ( iseed )

!*****************************************************************************80
!
!! DSTART is an entry point used to initialize DNOR.
!
  if ( iseed /= 0 ) then
!
!  generate random bit pattern in array based on given seed
!
    ii = 17
    jj = 5
    ia = mod ( abs ( iseed ), 32707 )
    ib = 1111
    ic = 1947
    do iii = 1, 17
      s = 0.0D+00
      t = 0.50
!
!  do for each of the bits of mantissa of word
!  loop over 64 bits, enough for all known machines in single precision
!
      do jjj = 1,64

        id = ic - ia
        if ( id < 0 ) then
          id = id + 32707
          s = s + t
        end if

        ia = ib
        ib = ic
        ic = id
        t = 0.5D+00 * t

      end do

      u(iii) = s

    end do

  end if
!
!  return floating echo of iseed.
!
  dstart = iseed

  return
end
function dnrm2 ( n, x, incx )

!*****************************************************************************80
!
!! DNRM2 returns the euclidean norm of a vector.
!
!  Discussion:
!
!     DNRM2 ( X ) = sqrt ( X' * X )
!
!  Author:
!
!    Sven Hammarling
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries of X.
!
!    Output, real ( kind = 8 ) DNRM2, the Euclidean norm of X.
!
  implicit none

  real ( kind = 8 ) absxi
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm
  real ( kind = 8 ) scale
  real ( kind = 8 ) ssq
  real ( kind = 8 ) x(*)

  if ( n < 1 .or. incx < 1 ) then

    norm  = 0.0D+00

  else if ( n == 1 ) then

    norm  = abs ( x(1) )

  else

    scale = 0.0D+00
    ssq = 1.0D+00

    do ix = 1, 1 + ( n - 1 )*incx, incx
      if ( x(ix) /= 0.0D+00 ) then
        absxi = abs ( x(ix) )
        if ( scale < absxi ) then
          ssq = 1.0D+00 + ssq * ( scale / absxi )**2
          scale = absxi
        else
          ssq = ssq + ( absxi / scale )**2
        end if
      end if
    end do
    norm  = scale * sqrt( ssq )
  end if

  dnrm2 = norm

  return
end
subroutine dnsq ( fcn, jac, iopt, n, x, fvec, fjac, ldfjac, xtol, maxfev, ml, &
  mu, epsfcn, diag, mode, factor, nprint, info, nfev, njev, r, lr, qtf, wa1, &
  wa2, wa3, wa4 )

!*****************************************************************************80
!
!! DNSQ finds a zero of a system of N nonlinear functions in N variables.
!
!  Discussion:
!
!    DNSQ uses a modification of the Powell hybrid method.  This code is the
!    combination of the MINPACK codes (argonne) hybrd and hybrdj.
!
!    The purpose of DNSQ is to find a zero of a system of N non-
!    linear functions in N variables by a modification of the powell
!    hybrid method.  The user must provide a subroutine which calcu-
!    lates the functions.  The user has the option of either to
!    provide a subroutine which calculates the jacobian or to let the
!    code calculate it by a forward-difference approximation.
!    This code is the combination of the minpack codes (argonne)
!    hybrd and hybrdj.
!
!  Reference:
!
!    MJD Powell,
!    A Hybrid Method for Nonlinear Equations,
!    Numerical Methods for Nonlinear Algebraic Equations,
!    P. Rabinowitz, editor.
!    Gordon and Breach, 1970.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  FCN must be declared in an external statement
!    in the user calling program, and should be written as follows.
!
!      subroutine fcn(n,x,fvec,iflag)
!      integer ( kind = 4 ) n,iflag
!      real ( kind = 8 ) x(n),fvec(n)
!
!      calculate the functions at x and return this vector in fvec.
!
!      return
!      end
!
!    The value of iflag should not be changed by fcn unless the
!    user wants to terminate execution of DNSQ.  in this case, set
!    iflag to a negative integer ( kind = 4 ).
!
!    Input, external JAC, the name of the user-supplied subroutine which
!    calculates the jacobian.  if iopt=1, then jac must be declared in an
!    external statement in the user calling program, and should be
!    written as follows.
!
!      subroutine jac(n,x,fvec,fjac,ldfjac,iflag)
!      integer ( kind = 4 ) n,ldfjac,iflag
!      real ( kind = 8 ) x(n),fvec(n),fjac(ldfjac,n)
!
!      calculate the jacobian at x and return this
!      matrix in fjac.  fvec contains the function
!      values at x and should not be altered.
!
!      return
!      end
!
!    the value of iflag should not be changed by jac unless the
!    user wants to terminate execution of DNSQ.  in this case, set
!    iflag to a negative integer ( kind = 4 ).
!    if iopt=2, jac can be ignored (treat it as a dummy argument).
!
!    Input, integer ( kind = 4 ) IOPT, specifies how the jacobian will be calculated.
!    1, the user must supply the jacobian through the subroutine jac.
!    2, the code will approximate the jacobian by forward-differencing.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input/output, real X(N).  On input, x must contain an initial
!    estimate of the solution vector.  on output, x contains the
!    final estimate of the solution vector.
!
!       fvec is an output array of length n which contains the functions
!         evaluated at the output x.
!
!       fjac is an output n by n array which contains the orthogonal
!         matrix q produced by the qr factorization of the final approx-
!         imate jacobian.
!
!       ldfjac is a positive integer ( kind = 4 ) input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       xtol is a non-negative input variable.  termination occurs when
!         the relative error between two consecutive iterates is at most
!         xtol.  therefore, xtol measures the relative error desired in
!         the approximate solution.  section 4 contains more details
!         about xtol.
!
!       maxfev is a positive integer ( kind = 4 ) input variable.  termination occurs
!         when the number of calls to fcn is at least maxfev by the end
!         of an iteration.
!
!       ml is a non-negative integer ( kind = 4 ) input variable which specifies the
!         number of subdiagonals within the band of the jacobian matrix.
!         if the jacobian is not banded or iopt=1, set ml to at
!         least n - 1.
!
!       mu is a non-negative integer ( kind = 4 ) input variable which specifies the
!         number of superdiagonals within the band of the jacobian
!         matrix.  if the jacobian is not banded or iopt=1, set mu to at
!         least n - 1.
!
!       epsfcn is an input variable used in determining a suitable step
!         for the forward-difference approximation.  this approximation
!         assumes that the relative errors in the functions are of the
!         order of epsfcn.  if epsfcn is less than the machine preci-
!         sion, it is assumed that the relative errors in the functions
!         are of the order of the machine precision.  if iopt=1, then
!         epsfcn can be ignored (treat it as a dummy argument).
!
!       diag is an array of length n.  if mode = 1 (see below), diag is
!         internally set.  if mode = 2, diag must contain positive
!         entries that serve as implicit (multiplicative) scale factors
!         for the variables.
!
!       mode is an integer ( kind = 4 ) input variable.  if mode = 1, the variables
!         will be scaled internally.  if mode = 2, the scaling is speci-
!         fied by the input diag.  other values of mode are equivalent
!         to mode = 1.
!
!       factor is a positive input variable used in determining the ini-
!         tial step bound.  this bound is set to the product of factor
!         and the euclidean norm of diag*x if nonzero, or else to factor
!         itself.  in most cases factor should lie in the interval
!         (.1,100.).  100. is a generally recommended value.
!
!       nprint is an integer ( kind = 4 ) input variable that enables controlled
!         printing of iterates if it is positive.  in this case, fcn is
!         called with iflag = 0 at the beginning of the first iteration
!         and every nprint iteration thereafter and immediately prior
!         to return, with x and fvec available for printing. appropriate
!         print statements must be added to fcn(see example).  if nprint
!         is not positive, no special calls of fcn with iflag = 0 are
!         made.
!
!       info is an integer ( kind = 4 ) output variable.  if the user has terminated
!         execution, info is set to the (negative) value of iflag.  see
!         description of fcn and jac. otherwise, info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  relative error between two consecutive iterates is
!                   at most xtol.
!
!         info = 2  number of calls to fcn has reached or exceeded
!                   maxfev.
!
!         info = 3  xtol is too small.  no further improvement in the
!                   approximate solution x is possible.
!
!         info = 4  iteration is not making good progress, as measured
!                   by the improvement from the last five jacobian eval-
!                   uations.
!
!         info = 5  iteration is not making good progress, as measured
!                   by the improvement from the last ten iterations.
!
!         sections 4 and 5 contain more details about info.
!
!       nfev is an integer ( kind = 4 ) output variable set to the number of calls to
!         fcn.
!
!       njev is an integer ( kind = 4 ) output variable set to the number of calls to
!         jac. (if iopt=2, then njev is set to zero.)
!
!       r is an output array of length lr which contains the upper
!         triangular matrix produced by the qr factorization of the
!         final approximate jacobian, stored rowwise.
!
!       lr is a positive integer ( kind = 4 ) input variable not less than
!         (n*(n+1))/2.
!
!       qtf is an output array of length n which contains the vector
!         (q transpose) * fvec.
!
!       wa1, wa2, wa3, and wa4 are work arrays of length n.
!
!
! 4. successful completion.
!
!       The accuracy of DNSQ is controlled by the convergence parameter
!       xtol.  this parameter is used in a test which makes a comparison
!       between the approximation x and a solution xsol.  DNSQ termi-
!       nates when the test is satisfied.  if the convergence parameter
!       is less than the machine precision (as defined by the function
!       epsilon), then DNSQ only attempts to satisfy the test
!       defined by the machine precision.  further progress is not
!       usually possible.
!
!       The test assumes that the functions are reasonably well behaved,
!       and, if the jacobian is supplied by the user, that the functions
!       and the jacobian are coded consistently.  if these conditions
!       are not satisfied, then DNSQ may incorrectly indicate conver-
!       gence.  the coding of the jacobian can be checked by the
!       subroutine chkder. if the jacobian is coded correctly or iopt=2,
!       then the validity of the answer can be checked, for example, by
!       rerunning DNSQ with a tighter tolerance.
!
!       Convergence test.  If dnrm2 ( z) denotes the euclidean norm of a
!         vector z and d is the diagonal matrix whose entries are
!         defined by the array diag, then this test attempts to guaran-
!         tee that
!
!               dnrm2 ( d*(x-xsol)) <= xtol * dnrm2(d*xsol).
!
!         if this condition is satisfied with xtol = 10**(-k), then the
!         larger components of d*x have k significant decimal digits and
!         info is set to 1.  there is a danger that the smaller compo-
!         nents of d*x may have large relative errors, but the fast rate
!         of convergence of DNSQ usually avoids this possibility.
!         unless high precision solutions are required, the recommended
!         value for xtol is the square root of the machine precision.
!
!
! 5. unsuccessful completion.
!
!       unsuccessful termination of DNSQ can be due to improper input
!       parameters, arithmetic interrupts, an excessive number of func-
!       tion evaluations, or lack of good progress.
!
!       improper input parameters.  info is set to 0 if iopt < 1,
!         or 2 < iopt, or n <= 0, or ldfjac < n, or
!         xtol < 0.0, or maxfev <= 0, or ml < 0, or mu < 0,
!         or factor <= 0.0, or lr < (n*(n+1))/2.
!
!       arithmetic interrupts.  if these interrupts occur in the fcn
!         subroutine during an early stage of the computation, they may
!         be caused by an unacceptable choice of x by DNSQ.  in this
!         case, it may be possible to remedy the situation by rerunning
!         DNSQ with a smaller value of factor.
!
!       excessive number of function evaluations.  a reasonable value
!         for maxfev is 100*(n+1) for iopt=1 and 200*(n+1) for iopt=2.
!         if the number of calls to fcn reaches maxfev, then this
!         indicates that the routine is converging very slowly as
!         measured by the progress of fvec, and info is set to 2.  this
!         situation should be unusual because, as indicated below, lack
!         of good progress is usually diagnosed earlier by DNSQ,
!         causing termination with info = 4 or info = 5.
!
!       lack of good progress.  DNSQ searches for a zero of the system
!         by minimizing the sum of the squares of the functions.  in so
!         doing, it can become trapped in a region where the minimum
!         does not correspond to a zero of the system and, in this situ-
!         ation, the iteration eventually fails to make good progress.
!         in particular, this will happen if the system does not have a
!         zero.  if the system has a zero, rerunning DNSQ from a dif-
!         ferent starting point may be helpful.
!
!
! 6. characteristics of the algorithm.
!
!       DNSQ is a modification of the Powell hybrid method.  Two of its
!       main characteristics involve the choice of the correction as a
!       convex combination of the Newton and scaled gradient directions,
!       and the updating of the jacobian by the rank-1 method of Broyden.
!       The choice of the correction guarantees (under reasonable
!       conditions) global convergence for starting points far from the
!       solution and a fast rate of convergence.  The jacobian is
!       calculated at the starting point by either the user-supplied
!       subroutine or a forward-difference approximation, but it is not
!       recalculated until the rank-1 method fails to produce satisfactory
!       progress.
!
!       timing.  The time required by DNSQ to solve a given problem
!         depends on N, the behavior of the functions, the accuracy
!         requested, and the starting point.  the number of arithmetic
!         operations needed by DNSQ is about 11.5*(n**2) to process
!         each evaluation of the functions (call to fcn) and 1.3*(n**3)
!         to process each evaluation of the jacobian (call to jac,
!         if iopt = 1).  unless fcn and jac can be evaluated quickly,
!         the timing of DNSQ will be strongly influenced by the time
!         spent in FCN and JAC.
!
!       storage.  DNSQ requires (3*n**2 + 17*n)/2 single precision
!         storage locations, in addition to the storage required by the
!         program.  there are no internally declared storage arrays.
!
!
! 7. example.
!
!       the problem is to determine the values of x(1), x(2), ..., x(9),
!       which solve the system of tridiagonal equations
!
!       (3-2*x(1))*x(1)           -2*x(2)                   = -1
!               -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
!                                   -x(8) + (3-2*x(9))*x(9) = -1
! c     **********
!
!       program test(input,output,tape6=output)
! c
! c     driver for DNSQ example.
! c
!       integer ( kind = 4 ) j,iopt,n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
!       integer ( kind = 4 ) nwrite
!       real ( kind = 8 ) xtol,epsfcn,factor,fnorm
!       real ( kind = 8 ) x(9),fvec(9),diag(9),fjac(9,9),r(45),qtf(9)
!       real ( kind = 8 ) wa1(9),wa2(9),wa3(9),wa4(9)
!       real ( kind = 8 ) dnrm2
!       external fcn
!       data nwrite /6/
!
!       iopt = 2
!       n = 9
! c
! c     the following starting values provide a rough solution.
! c
!       do j = 1, 9
!          x(j) = -1.0D+00
!       end do
!
!       ldfjac = 9
!       lr = 45
! c
! c  set xtol to the square root of the machine precision.
! c  unless high precision solutions are required,
! c  this is the recommended setting.
! c
!       xtol = sqrt ( epsilon ( xtol ) )
!
!       maxfev = 2000
!       ml = 1
!       mu = 1
!       epsfcn = 0.0D+00
!       mode = 2
!       do j = 1, 9
!          diag(j) = 1.0D+00
!       end do
!       factor = 1.e2
!       nprint = 0
!
!       call dnsq (fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,ml,mu,
!      *           epsfcn,diag,mode,factor,nprint,info,nfev,njev,
!      *           r,lr,qtf,wa1,wa2,wa3,wa4)
!       fnorm = dnrm2 ( n,fvec,1)
!       write (nwrite,1000) fnorm,nfev,info, x(1:n)
!       stop
!  1000 format (5x,' final l2 norm of the residuals',e15.7 //
!      *        5x,' number of function evaluations',i10 //
!      *        5x,' exit parameter',16x,i10 //
!      *        5x,' final approximate solution' // (5x,3e15.7))
!       end
!       subroutine fcn(n,x,fvec,iflag)
!       integer ( kind = 4 ) n,iflag
!       real ( kind = 8 ) x(n),fvec(n)
!       integer ( kind = 4 ) k
!       real temp,temp1,temp2
! c
!       if (iflag /= 0) go to 5
! c
! c     insert print statements here when nprint is positive.
! c
!       return
!     5 continue
!       do k = 1, n
!          temp = ( 3.0D+00 - 2.0D+00 * x(k) ) * x(k)
!          temp1 = 0.0D+00
!          if (k /= 1) temp1 = x(k-1)
!          temp2 = 0.0D+00
!          if (k /= n) temp2 = x(k+1)
!          fvec(k) = temp - temp1 - 2.0D+00 * temp2 + 1.0D+00
!       end do
!       return
!       end
!
!       results obtained with different compilers or machines
!       may be slightly different.
!
!       final l2 norm of the residuals  0.1192636e-07
!
!       number of function evaluations        14
!
!       exit parameter                         1
!
!       final approximate solution
!
!       -0.5706545e+00 -0.6816283e+00 -0.7017325e+00
!       -0.7042129e+00 -0.7013690e+00 -0.6918656e+00
!       -0.6657920e+00 -0.5960342e+00 -0.4164121e+00
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) n

  real ( kind = 8 ) actred
  real ( kind = 8 ) delta
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) epsfcn
  real ( kind = 8 ) epsmch
  real ( kind = 8 ) factor
  external fcn
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fnorm
  real ( kind = 8 ) fnorm1
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) iwa(1)
  integer ( kind = 4 ) j
  external jac
  logical jeval
  integer ( kind = 4 ) jm1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) ncfail
  integer ( kind = 4 ) ncsuc
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) njev
  integer ( kind = 4 ) nprint
  integer ( kind = 4 ) nslow1
  integer ( kind = 4 ) nslow2
  real ( kind = 8 ), parameter :: p001 = 0.001D+00
  real ( kind = 8 ), parameter :: p0001 = 0.0001D+00
  real ( kind = 8 ), parameter :: p1 = 0.1D+00
  real ( kind = 8 ), parameter :: p5 = 0.5D+00
  real ( kind = 8 ) pnorm
  real ( kind = 8 ) prered
  real ( kind = 8 ) qtf(n)
  real ( kind = 8 ) r(lr)
  real ( kind = 8 ) ratio
  logical sing
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) xnorm
  real ( kind = 8 ) wa1(n)
  real ( kind = 8 ) wa2(n)
  real ( kind = 8 ) wa3(n)
  real ( kind = 8 ) wa4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtol

  epsmch = epsilon ( epsmch )
  info = 0
  iflag = 0
  nfev = 0
  njev = 0

  if ( iopt < 1 .or. 2 < iopt ) then
    go to 300
  else if ( n <= 0 .or. xtol < 0.0D+00 .or. maxfev <= 0 ) then
    go to 300
  else if ( ml < 0 .or. mu < 0 .or. factor <= 0.0D+00 ) then
    go to 300
  else if ( ldfjac < n .or. lr < ( n * ( n + 1 ) ) / 2 ) then
    go to 300
  end if

  if ( mode /= 2 ) then
    go to 20
  end if

  do j = 1, n
    if ( diag(j) <= 0.0D+00 ) then
      go to 300
    end if
  end do

   20 continue
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn ( n, x, fvec, iflag )
  nfev = 1
  if ( iflag < 0 ) then
    go to 300
  end if

  fnorm = dnrm2 ( n, fvec, 1 )
!
!  Initialize iteration counter and monitors.
!
  iter = 1
  ncsuc = 0
  ncfail = 0
  nslow1 = 0
  nslow2 = 0
!
!  Beginning of the outer loop.
!
   30 continue

     jeval = .true.
!
!  Calculate the jacobian matrix.
!
     if ( iopt /= 2 ) then
!
!  User supplies jacobian
!
        call jac ( n, x, fvec, fjac, ldfjac, iflag )
        njev = njev + 1
!
!  Code approximates the jacobian
!
      else

        iflag = 2

        call fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )

        nfev = nfev + min ( ml+mu+1, n )

     end if

     if ( iflag < 0 ) then
       go to 300
     end if
!
!  Compute the qr factorization of the jacobian.
!
     call qrfac ( n, n, fjac, ldfjac, .false., iwa, 1, wa1, wa2 )
!
!  On the first iteration and if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
     if ( iter /= 1 ) then
       go to 70
     end if

     if ( mode == 2 ) then
       go to 50
     end if

     diag(1:n) = wa2(1:n)
     do j = 1, n
       if ( wa2(j) == 0.0D+00 ) then
         diag(j) = 1.0D+00
       end if
     end do

   50    continue
!
!  On the first iteration, calculate the norm of the scaled x
!  and initialize the step bound delta.
!
     wa3(1:n) = diag(1:n) * x(1:n)
     xnorm = dnrm2 ( n, wa3, 1 )

     delta = factor * xnorm
     if ( delta == 0.0D+00 ) then
       delta = factor
     end if

70   continue
!
!  Form Q' * FVEC and store in QTF.
!
     qtf(1:n) = fvec(1:n)

     do j = 1, n

       if ( fjac(j,j) /= 0.0D+00 ) then

        sum2 = 0.0D+00
        do i = j, n
          sum2 = sum2 + fjac(i,j) * qtf(i)
        end do

        temp = -sum2 / fjac(j,j)
        do i = j, n
          qtf(i) = qtf(i) + fjac(i,j) * temp
        end do

      end if

    end do
!
!  Copy the triangular factor of the qr factorization into r.
!
     sing = .false.

     do j = 1, n

        l = j
        jm1 = j - 1
        do i = 1, jm1
          r(l) = fjac(i,j)
          l = l + n - i
        end do
        r(l) = wa1(j)
        if ( wa1(j) == 0.0D+00 ) then
           sing = .true.
        end if

     end do
!
!  Accumulate the orthogonal factor in FJAC.
!
     call qform ( n, n, fjac, ldfjac )
!
!  Rescale if necessary.
!
     if ( mode /= 2 ) then

       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do

     end if
!
!  beginning of the inner loop.
!
  180    continue
!
!  If requested, call fcn to enable printing of iterates.
!
        if ( 0 < nprint ) then
          iflag = 0
          if ( mod ( iter-1, nprint ) == 0 ) then
            call fcn ( n, x, fvec, iflag )
          end if
          if ( iflag < 0 ) then
            go to 300
          end if

        end if
!
!  Determine the direction P.
!
        call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
!
!  Store the direction p and x + p. calculate the norm of p.
!
        wa1(1:n) = -wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)
        pnorm = dnrm2 ( n,wa3,1)
!
!  on the first iteration, adjust the initial step bound.
!
        if ( iter == 1 ) then
          delta = min (delta,pnorm)
        end if
!
!  Evaluate the function at x + p and calculate its norm.
!
        iflag = 1
        call fcn ( n, wa2, wa4, iflag )
        nfev = nfev + 1

        if ( iflag < 0 ) then
          go to 300
        end if

        fnorm1 = dnrm2 ( n, wa4, 1 )
!
!  Compute the scaled actual reduction.
!
        if ( fnorm1 < fnorm ) then
          actred = 1.0D+00 - ( fnorm1 / fnorm )**2
        else
          actred = -1.0D+00
        end if
!
!  Compute the scaled predicted reduction.
!
        l = 1
        do i = 1, n
           sum2 = 0.0D+00
           do j = i, n
             sum2 = sum2 + r(l) * wa1(j)
             l = l + 1
           end do
           wa3(i) = qtf(i) + sum2
        end do

        temp = dnrm2 (n, wa3, 1 )

        if ( temp < fnorm ) then
          prered = 1.0D+00 - ( temp / fnorm )**2
        else
          prered = 0.0D+00
        end if
!
!  Compute the ratio of the actual to the predicted reduction.
!
        if ( 0.0D+00 < prered ) then
          ratio = actred / prered
        else
          ratio = 0.0D+00
        end if
!
!  Update the step bound.
!
        if ( ratio < p1 ) then

           ncsuc = 0
           ncfail = ncfail + 1
           delta = p5 * delta

        else

           ncfail = 0
           ncsuc = ncsuc + 1

           if ( p5 <= ratio .or. 1 < ncsuc ) then
             delta = max ( delta, pnorm / p5 )
           end if

           if ( abs ( ratio - 1.0D+00 ) <= p1 ) then
             delta = pnorm / p5
           end if

        end if
!
!  Successful iteration.  Update x, fvec, and their norms.
!
        if ( p0001 <= ratio ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:n) = wa4(1:n)
          xnorm = dnrm2 ( n, wa2, 1 )
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  Determine the progress of the iteration.
!
        if ( p001 <= actred ) then
          nslow1 = 0
        else
          nslow1 = nslow1 + 1
        end if

        if ( jeval ) then
          nslow2 = nslow2 + 1
        end if

        if ( p1 <= actred ) then
          nslow2 = 0
        end if
!
!  Test for convergence.
!
        if ( delta <= xtol * xnorm .or. fnorm == 0.0D+00 ) then
          info = 1
        end if

        if ( info /= 0 ) then
          go to 300
        end if
!
!  Tests for termination and stringent tolerances.
!
        if ( maxfev <= nfev ) then
          info = 2
        end if

        if ( p1 * max ( p1 * delta, pnorm ) <= epsmch * xnorm ) then
          info = 3
        end if

        if ( nslow2 == 5 ) info = 4
        if ( nslow1 == 10 ) info = 5

        if ( info /= 0 ) then
          go to 300
        end if
!
!  Criterion for recalculating jacobian
!
        if ( ncfail == 2 ) then
          go to 290
        end if
!
!  Calculate the rank one modification to the jacobian
!  and update QTF if necessary.
!
        do j = 1, n
           sum2 = 0.0D+00
           do i = 1, n
             sum2 = sum2 + fjac(i,j) * wa4(i)
           end do
           wa2(j) = ( sum2 - wa3(j) ) / pnorm
           wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )

           if ( p0001 <= ratio ) then
             qtf(j) = sum2
           end if

        end do
!
!  Compute the QR factorization of the updated jacobian.
!
        call r1updt ( n,n,r,lr,wa1,wa2,wa3,sing )
        call d1mpyq ( n,n,fjac,ldfjac,wa2,wa3 )
        call d1mpyq ( 1, n, qtf, 1, wa2, wa3 )
!
!  end of the inner loop.
!
        jeval = .false.
        go to 180
  290    continue
!
!  end of the outer loop.
!
     go to 30

  300 continue
!
!  termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0
  if ( 0 < nprint ) then
    call fcn ( n, x, fvec, iflag )
  end if

  if ( info < 0 ) call xerror( &
    'DNSQ   -- execution terminated because user set iflag negative.',1,1)

  if (info == 0) call xerror( 'DNSQ   -- invalid input parameter.',2,1)
  if (info == 2) call xerror( 'DNSQ   -- too many function evaluations.',9,1)
  if (info == 3) then
    call xerror( 'DNSQ   -- xtol too small. no further improvement possible.', &
      3,1)
  end if

  if ( 4 < info ) then
    call xerror( 'DNSQ   -- iteration not making good progress.',1,1)
  end if

  return
end
subroutine dnsqe ( fcn, jac, iopt, n, x, fvec, tol, nprint, info, wa, lwa )

!*****************************************************************************80
!
!! DNSQE is the easy-to-use version of DNSQ.
!
!  Discussion:
!
!    DNSQE finds a zero of a system of N non-linear functions in N variables
!    by a modification of the Powell hybrid method.  This is done by using
!    the more general nonlinear equation solver DNSQ.  The user must provide
!    a subroutine which calculates the functions.  The user has the option
!    of either to provide a subroutine which calculates the jacobian or
!    to let the code calculate it by a forward-difference approximation.
!
!    This code is a combination of the MINPACK codes HYBRD1 and HYBRJ1.
!
!  Reference:
!
!    MJD Powell,
!    A Hybrid Method for Nonlinear Equations,
!    in Numerical Methods for Nonlinear Algebraic Equations,
!    edited by P. Rabinowitz,
!    Gordon and Breach, 1970.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which calculates
!         the functions.  fcn must be declared in an external statement
!         in the user calling program, and should be written as follows.
!
!       subroutine fcn(n,x,fvec,iflag)
!         integer ( kind = 4 ) n,iflag
!         real ( kind = 8 ) x(n),fvec(n)
!
!         calculate the functions at x and
!         return this vector in fvec.
!
!         return
!       end
!
!         the value of iflag should not be changed by fcn unless the
!         user wants to terminate execution of DNSQE.  in this case, set
!         iflag to a negative integer ( kind = 4 ).
!
!       jac is the name of the user-supplied subroutine which calculates
!         the jacobian.  if iopt=1, then jac must be declared in an
!         external statement in the user calling program, and should be
!         written as follows.
!
!       subroutine jac(n,x,fvec,fjac,ldfjac,iflag)
!         integer ( kind = 4 ) n,ldfjac,iflag
!         real ( kind = 8 ) x(n),fvec(n),fjac(ldfjac,n)
!
!         calculate the jacobian at x and return this
!         matrix in fjac.  fvec contains the function
!         values at x and should not be altered.
!
!         return
!       end
!
!         the value of iflag should not be changed by jac unless the
!         user wants to terminate execution of DNSQE.  in this case, set
!         iflag to a negative integer ( kind = 4 ).
!
!         if iopt=2, jac can be ignored (treat it as a dummy argument).
!
!    Input, integer ( kind = 4 ) IOPT, specifies how the jacobian will be calculated.
!    1, the user supplies the jacobian through the subroutine JAC.
!    2, the code will approximate the jacobian by forward-differencing.
!
!       n is a positive integer ( kind = 4 ) input variable set to the number of
!       functions and variables.
!
!       x is an array of length n.  on input, x must contain an initial
!         estimate of the solution vector.  on output, x contains the
!         final estimate of the solution vector.
!
!       fvec is an output array of length n which contains the functions
!         evaluated at the output x.
!
!       tol is a non-negative input variable.  termination occurs when
!         the algorithm estimates that the relative error between x and
!         the solution is at most tol.  section 4 contains more details
!         about tol.
!
!       nprint is an integer ( kind = 4 ) input variable that enables controlled
!         printing of iterates if it is positive.  in this case, fcn is
!         called with iflag = 0 at the beginning of the first iteration
!         and every nprint iteration thereafter and immediately prior
!         to return, with x and fvec available for printing. appropriate
!         print statements must be added to fcn (see example). if nprint
!         is not positive, no special calls of fcn with iflag = 0 are
!         made.
!
!       info is an integer ( kind = 4 ) output variable.  if the user has terminated
!         execution, info is set to the (negative) value of iflag.  see
!         description of fcn and jac. otherwise, info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error between
!                   x and the solution is at most tol.
!
!         info = 2  number of calls to fcn has reached or exceeded
!                   100*(n+1) for iopt=1 or 200*(n+1) for iopt=2.
!
!         info = 3  tol is too small.  no further improvement in the
!                   approximate solution x is possible.
!
!         info = 4  iteration is not making good progress.
!
!         sections 4 and 5 contain more details about info.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer ( kind = 4 ) input variable not less than
!         (3*n**2+13*n))/2.
!
!
!  successful completion.
!
!       the accuracy of DNSQE is controlled by the convergence parame-
!       ter tol.  this parameter is used in a test which makes a compar-
!       ison between the approximation x and a solution xsol.  DNSQE
!       terminates when the test is satisfied.  if tol is less than the
!       machine precision (as defined by the function epsilon), then
!       DNSQE attemps only to satisfy the test defined by the machine
!       precision.  further progress is not usually possible.  unless
!       high precision solutions are required, the recommended value
!       for tol is the square root of the machine precision.
!
!       the test assumes that the functions are reasonably well behaved,
!       and, if the jacobian is supplied by the user, that the functions
!       and the jacobian  coded consistently.  if these conditions
!       are not satisfied, DNSQE may incorrectly indicate convergence.
!       the coding of the jacobian can be checked by the subroutine
!       chkder.  if the jacobian is coded correctly or iopt=2, then
!       the validity of the answer can be checked, for example, by
!       rerunning DNSQE with a tighter tolerance.
!
!       convergence test.  if dnrm2 ( z) denotes the euclidean norm of a
!         vector z, then this test attempts to guarantee that
!
!               dnrm2 ( x-xsol) <=  tol * dnrm2(xsol).
!
!         if this condition is satisfied with tol = 10**(-k), then the
!         larger components of x have k significant decimal digits and
!         info is set to 1.  there is a danger that the smaller compo-
!         nents of x may have large relative errors, but the fast rate
!         of convergence of DNSQE usually avoids this possibility.
!
!
!  unsuccessful completion.
!
!       unsuccessful termination of DNSQE can be due to improper input
!       parameters, arithmetic interrupts, an excessive number of func-
!       tion evaluations, errors in the functions, or lack of good prog-
!       ress.
!
!       improper input parameters.  info is set to 0 if iopt < 1, or
!         2 < iopt, or n <= 0, or tol < 0.0, or
!         lwa < (3*n**2+13*n)/2.
!
!       arithmetic interrupts.  if these interrupts occur in the fcn
!       subroutine during an early stage of the computation, they may
!         be caused by an unacceptable choice of x by DNSQE.  in this
!         case, it may be possible to remedy the situation by not evalu-
!         ating the functions here, but instead setting the components
!         of fvec to numbers that exceed those in the initial fvec.
!
!       excessive number of function evaluations.  if the number of
!         calls to fcn reaches 100*(n+1) for iopt=1 or 200*(n+1) for
!         iopt=2, then this indicates that the routine is converging
!         very slowly as measured by the progress of fvec, and info is
!         set to 2.  this situation should be unusual because, as
!         indicated below, lack of good progress is usually diagnosed
!         earlier by DNSQE, causing termination with info = 4.
!
!       errors in the functions.  when iopt=2, the choice of step length
!         in the forward-difference approximation to the jacobian
!         assumes that the relative errors in the functions are of the
!         order of the machine precision.  if this is not the case,
!         DNSQE may fail (usually with info = 4).  the user should
!         then either use DNSQ and set the step length or use iopt=1
!         and supply the jacobian.
!
!       lack of good progress.  DNSQE searches for a zero of the system
!         by minimizing the sum of the squares of the functions.  in so
!         doing, it can become trapped in a region where the minimum
!         does not correspond to a zero of the system and, in this situ-
!         ation, the iteration eventually fails to make good progress.
!         in particular, this will happen if the system does not have a
!         zero.  if the system has a zero, rerunning DNSQE from a dif-
!         ferent starting point may be helpful.
!
!
!  characteristics of the algorithm.
!
!       DNSQE is a modification of the powell hybrid method.  two of
!       its main characteristics involve the choice of the correction as
!       a convex combination of the Newton and scaled gradient direc-
!       tions, and the updating of the jacobian by the rank-1 method of
!       broyden.  the choice of the correction guarantees (under reason-
!       able conditions) global convergence for starting points far from
!       the solution and a fast rate of convergence.  the jacobian is
!       calculated at the starting point by either the user-supplied
!       subroutine or a forward-difference approximation, but it is not
!       recalculated until the rank-1 method fails to produce satis-
!       factory progress.
!
!       timing.  the time required by DNSQE to solve a given problem
!         depends on n, the behavior of the functions, the accuracy
!         requested, and the starting point.  the number of arithmetic
!         operations needed by DNSQE is about 11.5*(n**2) to process
!         each evaluation of the functions (call to fcn) and 1.3*(n**3)
!         to process each evaluation of the jacobian (call to jac,
!         if iopt = 1).  unless fcn and jac can be evaluated quickly,
!         the timing of DNSQE will be strongly influenced by the time
!         spent in fcn and jac.
!
!       storage.  DNSQE requires (3*n**2 + 17*n)/2 single precision
!         storage locations, in addition to the storage required by the
!         program.  there are no internally declared storage arrays.
!
!
!  example.
!
!       the problem is to determine the values of x(1), x(2), ..., x(9),
!       which solve the system of tridiagonal equations
!
!       (3-2*x(1))*x(1)           -2*x(2)                   = -1
!               -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
!                                   -x(8) + (3-2*x(9))*x(9) = -1
!
!       program test(input,output,tape6=output)
! c
! c     driver for DNSQE example.
! c
!       integer ( kind = 4 ) j,n,iopt,nprint,info,lwa,nwrite
!       real ( kind = 8 ) tol,fnorm
!       real ( kind = 8 ) x(9),fvec(9),wa(180)
!       real ( kind = 8 ) dnrm2
!       external fcn
!       data nwrite /6/
! c
!       iopt = 2
!       n = 9
! c
! c     the following starting values provide a rough solution.
! c
!       x(1:9) = -1.0D+00
!
!       lwa = 180
!       nprint = 0
! c
! c     set tol to the square root of the machine precision.
! c     unless high precision solutions are required,
! c     this is the recommended setting.
! c
!       tol = sqrt ( epsilon ( tol ) )
!
!       call dnsqe (fcn,jac,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
!       fnorm = dnrm2 ( n,fvec)
!       write (nwrite,1000) fnorm,info, x(1:n)
!       stop
!  1000 format (5x,' final l2 norm of the residuals',e15.7 //
!      *        5x,' exit parameter',16x,i10 //
!      *        5x,' final approximate solution' // (5x,3e15.7))
!     end
!     subroutine fcn(n,x,fvec,iflag)
!       integer ( kind = 4 ) n,iflag
!       real ( kind = 8 ) x(n),fvec(n)
!       integer ( kind = 4 ) k
!       real ( kind = 8 ) temp,temp1,temp2
!
!       do k = 1, n
!          temp = (3.0D+00 - 2.0D+00 * x(k) ) * x(k)
!          temp1 = 0.0D+00
!          if (k /= 1) temp1 = x(k-1)
!          temp2 = 0.0D+00
!          if (k /= n) temp2 = x(k+1)
!          fvec(k) = temp - temp1 - 2.0D+00 * temp2 + 1.0D+00
!       end do
!
!       return
!     end
!
!       results obtained with different compilers or machines
!       may be slightly different.
!
!       final l2 norm of the residuals  0.1192636e-07
!
!       exit parameter                         1
!
!       final approximate solution
!
!       -0.5706545e+00 -0.6816283e+00 -0.7017325e+00
!       -0.7042129e+00 -0.7013690e+00 -0.6918656e+00
!       -0.6657920e+00 -0.5960342e+00 -0.4164121e+00
!
  implicit none

  integer ( kind = 4 ) lwa
  integer ( kind = 4 ) n

  real ( kind = 8 ) epsfcn
  real ( kind = 8 ), parameter :: factor = 100.0D+00
  external fcn
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) index
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iopt
  external jac
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) njev
  integer ( kind = 4 ) nprint
  real ( kind = 8 ) tol
  real ( kind = 8 ) wa(lwa)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtol

  info = 0
!
!  Check the input parameters for errors.
!
  if ( iopt < 1 ) then
    call xerror ( 'DNSQE  -- invalid input parameter.', 2, 1 )
    return
  end if

  if ( 2 < iopt ) then
    call xerror ( 'DNSQE  -- invalid input parameter.', 2, 1 )
    return
  end if

  if ( n <= 0 .or. tol < 0.0D+00 .or. &
    lwa < (3*n**2 +13*n)/2 ) then
    if ( info == 0 ) then
      call xerror ( 'DNSQE  -- invalid input parameter.', 2, 1 )
    end if
    return
  end if

  maxfev = 100 * ( n + 1 )
  if ( iopt == 2 ) then
    maxfev = 2 * maxfev
  end if

  xtol = tol
  ml = n - 1
  mu = n - 1
  epsfcn = 0.0D+00
  mode = 2
  wa(1:n) = 1.0D+00
  lr = ( n * ( n + 1 ) ) / 2
  index = 6 * n + lr

  call dnsq ( fcn, jac, iopt, n, x, fvec, wa(index+1), n, xtol, maxfev, ml, &
    mu, epsfcn, wa(1), mode, factor, nprint, info, nfev, njev, &
    wa(6*n+1), lr, wa(n+1), wa(2*n+1), wa(3*n+1), wa(4*n+1), wa(5*n+1) )

  if ( info == 5 ) then
    info = 4
  end if

  if ( info == 0 ) then
    call xerror ( 'DNSQE  -- invalid input parameter.', 2, 1 )
  end if

  return
end
subroutine dogdrv ( nr, n, x, f, g, a, p, xpls, fpls, fcn, sx, stepmx, &
  steptl, dlt, iretcd, mxtake, sc, wrk1, wrk2, wrk3, ipr )

!*****************************************************************************80
!
!! DOGDRV finds the next Newton iterate by the double dogleg method.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, "X[K-1]".
!
!    Input, real ( kind = 8 ) F, the function value at the old iterate, "F(X)".
!
!    Input, real ( kind = 8 ) G(N), the gradient at the old iterate.
!
!    Input, real ( kind = 8 ) A(N,N), the Cholesky decomposition of the
!    Hessian matrix in lower triangular part and diagonal.
!
!    Input, real ( kind = 8 ) P(N), the Newton step.
!
!    Output, real ( kind = 8 ) XPLS(N), the new iterate "X[K]".
!
!    Output, real ( kind = 8 ) FPLS, the function value at the new iterate,
!    F(XPLS).
!
!    Input, external FCN, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real x(n)
!      real f
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which
!    successive iterates are considered close enough to terminate algorithm.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!    [retain value between successive calls].
!
!    Output, integer ( kind = 4 ) IRETCD, the return code.
!    0, satisfactory XPLS found
!    1, failed to find satisfactory XPLS sufficiently distinct from X.
!
!    Output, logical MXTAKE, TRUE if a step of maximum length was used.
!
!    Workspace, real ( kind = 8 ) SC(N), holds the current step.
!
!    Workspace, real ( kind = 8 ) WRK1(N).
!
!    Workspace, real ( kind = 8 ) WRK2(N).
!
!    Workspace, real ( kind = 8 ) WRK3(N).
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) cln
  real ( kind = 8 ) dlt
  real ( kind = 8 ) eta
  real ( kind = 8 ) f
  external fcn
  real ( kind = 8 ) fpls
  real ( kind = 8 ) fplsp
  logical fstdog
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) iretcd
  logical mxtake
  logical nwtake
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) rnwtln
  real ( kind = 8 ) sc(n)
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) steptl
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) wrk1(n)
  real ( kind = 8 ) wrk2(n)
  real ( kind = 8 ) wrk3(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xpls(n)

  iretcd = 4
  fstdog = .true.

  rnwtln = sqrt ( sum ( sx(1:n)**2 * p(1:n)**2 ) )

  do
!
!  Find new step by double dogleg algorithm.
!
    call dogstp ( nr, n, g, a, p, sx, rnwtln, dlt, nwtake, fstdog, wrk1, &
      wrk2, cln, eta, sc, ipr, stepmx )
!
!  Check new point and update trust region.
!
    call tregup ( nr, n, x, f, g, a, fcn, sc, sx, nwtake, stepmx, steptl, dlt, &
      iretcd, wrk3, fplsp, xpls, fpls, mxtake, ipr, 2, wrk1 )

    if ( iretcd <= 1 ) then
      exit
    end if

  end do

  return
end
subroutine dogleg ( n, r, lr, diag, qtb, delta, x )

!*****************************************************************************80
!
!! DOGLEG finds the minimizing combination of Gauss-Newton and gradient steps.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA, the
!    problem is to determine the convex combination X of the
!    Gauss-Newton and scaled gradient directions that minimizes
!    (A*X - B) in the least squares sense, subject to the
!    restriction that the euclidean norm of D*X be at most DELTA.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization of A.  That is, if A = Q*R, where Q has
!    orthogonal columns and R is an upper triangular matrix,
!    then DOGLEG expects the full upper triangle of R and
!    the first N components of Q'*B.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix R.
!
!    Input, real ( kind = 8 ) R(LR), the upper triangular matrix R stored
!    by rows.
!
!    Input, integer ( kind = 4 ) LR, the size of the R array, which must be no less
!    than (N*(N+1))/2.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'* B.
!
!    Input, real ( kind = 8 ) DELTA, is a positive upper bound on the
!    euclidean norm of D*X(1:N).
!
!    Output, real ( kind = 8 ) X(N), the desired convex combination of the
!    Gauss-Newton direction and the scaled gradient direction.
!
  implicit none

  integer ( kind = 4 ) lr
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) bnorm
  real ( kind = 8 ) delta
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) enorm
  real ( kind = 8 ) epsmch
  real ( kind = 8 ) gnorm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) qnorm
  real ( kind = 8 ) qtb(n)
  real ( kind = 8 ) r(lr)
  real ( kind = 8 ) sgnorm
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa1(n)
  real ( kind = 8 ) wa2(n)
  real ( kind = 8 ) x(n)

  epsmch = epsilon ( epsmch )
!
!  Calculate the Gauss-Newton direction.
!
  jj = ( n * ( n + 1 ) ) / 2 + 1

  do k = 1, n

     j = n - k + 1
     jj = jj - k
     l = jj + 1

     sum2 = 0.0D+00
     do i = j+1, n
       sum2 = sum2 + r(l) * x(i)
       l = l + 1
     end do

     temp = r(jj)

     if ( temp == 0.0D+00 ) then

       l = j
       do i = 1, j
         temp = max ( temp, abs ( r(l)) )
         l = l + n - i
       end do

       if ( temp == 0.0D+00 ) then
         temp = epsmch
       else
         temp = epsmch * temp
       end if

     end if

     x(j) = ( qtb(j) - sum2 ) / temp

  end do
!
!  Test whether the Gauss-Newton direction is acceptable.
!
  wa1(1:n) = 0.0D+00
  wa2(1:n) = diag(1:n) * x(1:n)
  qnorm = enorm ( n, wa2 )

  if ( qnorm <= delta ) then
    return
  end if
!
!  The Gauss-Newton direction is not acceptable.
!  Calculate the scaled gradient direction.
!
  l = 1
  do j = 1, n
    temp = qtb(j)
    do i = j, n
      wa1(i) = wa1(i) + r(l) * temp
      l = l + 1
    end do
    wa1(j) = wa1(j) / diag(j)
  end do
!
!  Calculate the norm of the scaled gradient.
!  Test for the special case in which the scaled gradient is zero.
!
  gnorm = enorm ( n, wa1 )
  sgnorm = 0.0D+00
  alpha = delta / qnorm

  if ( gnorm /= 0.0D+00 ) then
!
!  Calculate the point along the scaled gradient which minimizes the quadratic.
!
    wa1(1:n) = ( wa1(1:n) / gnorm ) / diag(1:n)

    l = 1
    do j = 1, n
      sum2 = 0.0D+00
      do i = j, n
        sum2 = sum2 + r(l) * wa1(i)
        l = l + 1
      end do
      wa2(j) = sum2
    end do

    temp = enorm ( n, wa2 )
    sgnorm = ( gnorm / temp ) / temp
!
!  Test whether the scaled gradient direction is acceptable.
!
    alpha = 0.0D+00
!
!  The scaled gradient direction is not acceptable.
!  Calculate the point along the dogleg at which the quadratic is minimized.
!
    if ( sgnorm < delta ) then

      bnorm = enorm ( n, qtb )
      temp = ( bnorm / gnorm ) * ( bnorm / qnorm ) * ( sgnorm / delta )
      temp = temp - ( delta / qnorm ) * ( sgnorm / delta)**2 &
        + sqrt ( ( temp - ( delta / qnorm ) )**2 &
        + ( 1.0D+00 - ( delta / qnorm )**2 ) &
        * ( 1.0D+00 - ( sgnorm / delta )**2 ) )

      alpha = ( ( delta / qnorm ) * ( 1.0D+00 - ( sgnorm / delta )**2 ) ) / temp

    end if

  end if
!
!  Form appropriate convex combination of the Gauss-Newton
!  direction and the scaled gradient direction.
!
  temp = ( 1.0D+00 - alpha ) * min ( sgnorm, delta )

  x(1:n) = temp * wa1(1:n) + alpha * x(1:n)

  return
end
subroutine dogstp ( nr, n, g, a, p, sx, rnwtln, dlt, nwtake, fstdog, ssd, v, &
  cln, eta, sc, ipr, stepmx )

!*****************************************************************************80
!
!! DOGSTP finds a new step by the double dogleg algorithm.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) G(N), the gradient at the current iterate.
!
!    Input, real ( kind = 8 ) A(NR,N), the Cholesky decomposition of the
!    hessian in the lower triangle and diagonal.
!
!    Input, real ( kind = 8 ) P(N), the Newton step.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) RNWTLN, the Newton step length.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Input/output, logical NWTAKE, TRUE if a Newton step was taken.
!
!    Input/output, logical FSTDOG, TRUE if on first leg of dogleg.
!
!    Input/output, real ( kind = 8 ) SSD(N), workspace [cauchy step to
!    the minimum of the quadratic model in the scaled steepest descent
!    direction] [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) V(N), workspace  [retain value
!    between successive calls]
!
!    Workspace, real ( kind = 8 ) CLN, the cauchy length.
!    [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) ETA, [retain value between successive calls]
!
!    Output, real ( kind = 8 ) SC(N), the current step.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!  Local variables:
!
!    CLN, the length of cauchy step
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) alam
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) cln
  real ( kind = 8 ) dlt
  real ( kind = 8 ) dot1
  real ( kind = 8 ) dot2
  real ( kind = 8 ) eta
  logical fstdog
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) j
  logical nwtake
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) rnwtln
  real ( kind = 8 ) sc(n)
  real ( kind = 8 ) ssd(n)
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) tmp
  real ( kind = 8 ) v(n)
!
!  Can we take a Newton step?
!
  if ( rnwtln <= dlt ) then

    nwtake = .true.
    sc(1:n) = p(1:n)
    dlt = rnwtln

  else
!
!  The Newton step is too long.
!  The Cauchy step is on double dogleg curve.
!
    nwtake = .false.

    if ( fstdog ) then
!
!  Calculate double dogleg curve, SSD.
!
      fstdog = .false.
      alpha = sum ( ( g(1:n) / sx(1:n) )**2 )
      beta = 0.0D+00
      do i = 1, n
        tmp = 0.0D+00
        do j = i, n
          tmp = tmp + ( a(j,i) * g(j) ) / ( sx(j) * sx(j) )
        end do
        beta = beta + tmp * tmp
      end do

      ssd(1:n) = - ( alpha / beta ) * g(1:n) / sx(1:n)

      cln = alpha * sqrt ( alpha ) / beta

      eta = 0.2D+00 + ( 0.8D+00 * alpha * alpha ) / &
        ( - beta * dot_product ( g, p ) )

      v(1:n) = eta * sx(1:n) * p(1:n) - ssd(1:n)

      if ( dlt == - 1.0D+00 ) then
        dlt = min ( cln, stepmx )
      end if

    end if
!
!  Take a partial step in the Newton direction.
!
      if ( eta * rnwtln <= dlt ) then

        sc(1:n) = ( dlt / rnwtln ) * p(1:n)
!
!  Take a step in steepest descent direction.
!
      else if ( dlt <= cln ) then

        sc(1:n) = ( dlt / cln ) * ssd(1:n) / sx(1:n)
!
!  Convex combination of SSD and eta*p which has scaled length DLT.
!
      else

        dot1 = dot_product ( v, ssd )
        dot2 = dot_product ( v, v )
        alam = ( -dot1 + sqrt ( ( dot1 * dot1 ) &
          - dot2 * ( cln * cln - dlt * dlt ) ) ) / dot2
        sc(1:n) = ( ssd(1:n) + alam * v(1:n) ) / sx(1:n)

      end if

  end if

  return
end
subroutine dqrank ( a, lda, m, n, tol, kr, jpvt, qraux, work )

!*****************************************************************************80
!
!! DQRANK computes the QR factorization of a rectangular matrix.
!
!  Discussion:
!
!    This routine is used in conjunction with sqrlss to solve
!    overdetermined, underdetermined and singular linear systems
!    in a least squares sense.
!
!    DQRANK uses the LINPACK subroutine DQRDC to compute the QR
!    factorization, with column pivoting, of an M by N matrix A.
!    The numerical rank is determined using the tolerance TOL.
!
!    Note that on output, ABS ( A(1,1) ) / ABS ( A(KR,KR) ) is an estimate
!    of the condition number of the matrix of independent columns,
!    and of R.  This estimate will be <= 1/TOL.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the matrix whose
!    decomposition is to be computed.  On output, the information from DQRDC.
!    The triangular matrix R of the QR factorization is contained in the
!    upper triangle and information needed to recover the orthogonal
!    matrix Q is stored below the diagonal in A and in the vector QRAUX.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be at least M.
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, real ( kind = 8 ) TOL, a relative tolerance used to determine the
!    numerical rank.  The problem should be scaled so that all the elements
!    of A have roughly the same absolute accuracy, EPS.  Then a reasonable
!    value for TOL is roughly EPS divided by the magnitude of the largest
!    element.
!
!    Output, integer ( kind = 4 ) KR, the numerical rank.
!
!    Output, integer ( kind = 4 ) JPVT(N), the pivot information from DQRDC.
!    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
!    independent to within the tolerance TOL and the remaining columns
!    are linearly dependent.
!
!    Output, real ( kind = 8 ) QRAUX(N), will contain extra information defining
!    the QR factorization.
!
!    Workspace, real ( kind = 8 ) WORK(N).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jpvt(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kr
  integer ( kind = 4 ) m
  real ( kind = 8 ) qraux(n)
  real ( kind = 8 ) tol
  real ( kind = 8 ) work(n)

  jpvt(1:n) = 0

  call dqrdc ( a, lda, m, n, qraux, jpvt, work, 1 )

  kr = 0
  k = min ( m, n )

  do j = 1, k
    if ( abs ( a(j,j) ) <= tol * abs ( a(1,1) ) ) then
      return
    end if
    kr = j
  end do

  return
end
subroutine dqrdc ( a, lda, n, p, qraux, jpvt, work, job )

!*****************************************************************************80
!
!! DQRDC computes the QR factorization of a real rectangular matrix.
!
!  Discussion:
!
!    DQRDC uses Householder transformations.
!
!    Column pivoting based on the 2-norms of the reduced columns may be
!    performed at the user's option.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,P).  On input, the N by P matrix
!    whose decomposition is to be computed.  On output, A contains in
!    its upper triangle the upper triangular matrix R of the QR
!    factorization.  Below its diagonal A contains information from
!    which the orthogonal part of the decomposition can be recovered.
!    Note that if pivoting has been requested, the decomposition is not that
!    of the original matrix A but that of A with its columns permuted
!    as described by JPVT.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix A.
!
!    Input, integer ( kind = 4 ) P, the number of columns of the matrix A.
!
!    Output, real ( kind = 8 ) QRAUX(P), contains further information required
!    to recover the orthogonal part of the decomposition.
!
!    Input/output, integer ( kind = 4 ) JPVT(P).  On input, JPVT contains
!    integers that control the selection of the pivot columns.  The K-th
!    column A(*,K) of A is placed in one of three classes according to the
!    value of JPVT(K).
!      > 0, then A(K) is an initial column.
!      = 0, then A(K) is a free column.
!      < 0, then A(K) is a final column.
!    Before the decomposition is computed, initial columns are moved to
!    the beginning of the array A and final columns to the end.  Both
!    initial and final columns are frozen in place during the computation
!    and only free columns are moved.  At the K-th stage of the
!    reduction, if A(*,K) is occupied by a free column it is interchanged
!    with the free column of largest reduced norm.  JPVT is not referenced
!    if JOB == 0.  On output, JPVT(K) contains the index of the column of the
!    original matrix that has been interchanged into the K-th column, if
!    pivoting was requested.
!
!    Workspace, real ( kind = 8 ) WORK(P).  WORK is not referenced if JOB == 0.
!
!    Input, integer ( kind = 4 ) JOB, initiates column pivoting.
!    0, no pivoting is done.
!    nonzero, pivoting is done.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  real ( kind = 8 ) a(lda,p)
  integer ( kind = 4 ) jpvt(p)
  real ( kind = 8 ) qraux(p)
  real ( kind = 8 ) work(p)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lup
  integer ( kind = 4 ) maxj
  real ( kind = 8 ) maxnrm
  real ( kind = 8 ) nrmxl
  integer ( kind = 4 ) pl
  integer ( kind = 4 ) pu
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dnrm2
  logical swapj
  real ( kind = 8 ) t
  real ( kind = 8 ) tt

  pl = 1
  pu = 0
!
!  If pivoting is requested, rearrange the columns.
!
  if ( job /= 0 ) then

    do j = 1, p

      swapj = 0 < jpvt(j)

      if ( jpvt(j) < 0 ) then
        jpvt(j) = - j
      else
        jpvt(j) = j
      end if

      if ( swapj ) then

        if ( j /= pl ) then
          call dswap ( n, a(1,pl), 1, a(1,j), 1 )
        end if

        jpvt(j) = jpvt(pl)
        jpvt(pl) = j
        pl = pl + 1

      end if

    end do

    pu = p

    do j = p, 1, -1

      if ( jpvt(j) < 0 ) then

        jpvt(j) = - jpvt(j)

        if ( j /= pu ) then
          call dswap ( n, a(1,pu), 1, a(1,j), 1 )
          jp = jpvt(pu)
          jpvt(pu) = jpvt(j)
          jpvt(j) = jp
        end if

        pu = pu - 1

      end if

    end do

  end if
!
!  Compute the norms of the free columns.
!
  do j = pl, pu
    qraux(j) = dnrm2 ( n, a(1,j), 1 )
  end do

  work(pl:pu) = qraux(pl:pu)
!
!  Perform the Householder reduction of A.
!
  lup = min ( n, p )

  do l = 1, lup
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pl <= l .and. l < pu ) then

      maxnrm = 0.0D+00
      maxj = l
      do j = l, pu
        if ( maxnrm < qraux(j) ) then
          maxnrm = qraux(j)
          maxj = j
        end if
      end do

      if ( maxj /= l ) then
        call dswap ( n, a(1,l), 1, a(1,maxj), 1 )
        qraux(maxj) = qraux(l)
        work(maxj) = work(l)
        jp = jpvt(maxj)
        jpvt(maxj) = jpvt(l)
        jpvt(l) = jp
      end if

    end if
!
!  Compute the Householder transformation for column L.
!
    qraux(l) = 0.0D+00

    if ( l /= n ) then

      nrmxl = dnrm2 ( n-l+1, a(l,l), 1 )

      if ( nrmxl /= 0.0D+00 ) then

        if ( a(l,l) /= 0.0D+00 ) then
          nrmxl = sign ( nrmxl, a(l,l) )
        end if

        call dscal ( n-l+1, 1.0D+00 / nrmxl, a(l,l), 1 )
        a(l,l) = 1.0D+00 + a(l,l)
!
!  Apply the transformation to the remaining columns, updating the norms.
!
        do j = l + 1, p

          t = - ddot ( n-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
          call daxpy ( n-l+1, t, a(l,l), 1, a(l,j), 1 )

          if ( pl <= j .and. j <= pu ) then

            if ( qraux(j) /= 0.0D+00 ) then

              tt = 1.0D+00 - ( abs ( a(l,j) ) / qraux(j) )**2
              tt = max ( tt, 0.0D+00 )
              t = tt
              tt = 1.0D+00 + 0.05D+00 * tt * ( qraux(j) / work(j) )**2

              if ( tt /= 1.0D+00 ) then
                qraux(j) = qraux(j) * sqrt ( t )
              else
                qraux(j) = dnrm2 ( n-l, a(l+1,j), 1 )
                work(j) = qraux(j)
              end if

            end if

          end if

        end do
!
!  Save the transformation.
!
        qraux(l) = a(l,l)
        a(l,l) = - nrmxl

      end if

    end if

  end do

  return
end
subroutine dqrls ( a, lda, m, n, tol, kr, b, x, rsd, jpvt, qraux, work, &
  itask, ind )

!*****************************************************************************80
!
!! DQRLS factors and solves a linear system in the least squares sense.
!
!  Discussion:
!
!    The linear system may be overdetermined, underdetermined or singular.
!    The solution is obtained using a QR factorization of the
!    coefficient matrix.
!
!    DQRLS can be efficiently used to solve several least squares
!    problems with the same matrix A.  The first system is solved
!    with ITASK = 1.  The subsequent systems are solved with
!    ITASK = 2, to avoid the recomputation of the matrix factors.
!    The parameters KR, JPVT, and QRAUX must not be modified
!    between calls to DQRLS.
!
!    DQRLS is used to solve in a least squares sense
!    overdetermined, underdetermined and singular linear systems.
!    The system is A*X approximates B where A is M by N.
!    B is a given M-vector, and X is the N-vector to be computed.
!    A solution X is found which minimimzes the sum of squares (2-norm)
!    of the residual,  A*X - B.
!
!    The numerical rank of A is determined using the tolerance TOL.
!
!    DQRLS uses the LINPACK subroutine DQRDC to compute the QR
!    factorization, with column pivoting, of an M by N matrix A.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N), an M by N matrix.
!    On input, the matrix whose decomposition is to be computed.
!    In a least squares data fitting problem, A(I,J) is the
!    value of the J-th basis (model) function at the I-th data point.
!    On output, A contains the output from DQRDC.  The triangular matrix R
!    of the QR factorization is contained in the upper triangle and
!    information needed to recover the orthogonal matrix Q is stored
!    below the diagonal in A and in the vector QRAUX.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    M <= LDA.
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, real ( kind = 8 ) TOL, a relative tolerance used to determine the
!    numerical rank.  The problem should be scaled so that all the elements
!    of A have roughly the same absolute accuracy EPS.  Then a reasonable
!    value for TOL is roughly EPS divided by the magnitude of the largest
!    element.
!
!    Output, integer ( kind = 4 ) KR, the numerical rank.
!
!    Input, real ( kind = 8 ) B(M), the right hand side of the linear system.
!    In a least squares data fitting problem B(I) contains the
!    value of the I-th observation.
!
!    Output, real ( kind = 8 ) X(N), a least squares solution to the linear
!    system.
!
!    Output, real ( kind = 8 ) RSD(M), the residual, B - A*X.  RSD may
!    overwrite B.
!
!    Workspace, integer ( kind = 4 ) JPVT(N), required if ITASK = 1.
!    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
!    independent to within the tolerance TOL and the remaining columns
!    are linearly dependent.  ABS ( A(1,1) ) / ABS ( A(KR,KR) ) is an estimate
!    of the condition number of the matrix of independent columns,
!    and of R.  This estimate will be <= 1/TOL.
!
!    Workspace, real ( kind = 8 ) QRAUX(N), required if ITASK = 1.
!
!    Workspace, real ( kind = 8 ) WORK(N), required if ITASK = 1.
!
!    Input, integer ( kind = 4 ) ITASK.
!    1, DQRLS factors the matrix A and solves the least squares problem.
!    2, DQRLS assumes that the matrix A was factored with an earlier
!       call to DQRLS, and only solves the least squares problem.
!
!    Output, integer ( kind = 4 ) IND, error code.
!    0:  no error
!    -1: LDA < N   (fatal error)
!    -2: N < 1     (fatal error)
!    -3: ITASK < 1 (fatal error)
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) itask
  integer ( kind = 4 ) jpvt(n)
  integer ( kind = 4 ) kr
  real ( kind = 8 ) qraux(n)
  real ( kind = 8 ) rsd(m)
  real ( kind = 8 ) tol
  real ( kind = 8 ) work(n)
  real ( kind = 8 ) x(n)

  if ( lda < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DQRLS - Fatal error!'
    write ( *, '(a)' ) '  LDA < M.'
    stop
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DQRLS - Fatal error!'
    write ( *, '(a)' ) '  N <= 0.'
    stop
  end if

  if ( itask < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DQRLS - Fatal error!'
    write ( *, '(a)' ) '  ITASK < 1.'
    stop
  end if

  ind = 0
!
!  Factor the matrix.
!
  if ( itask == 1 ) then
    call dqrank ( a, lda, m, n, tol, kr, jpvt, qraux, work )
  end if
!
!  Solve the least-squares problem.
!
  call dqrlss ( a, lda, m, n, kr, b, x, rsd, jpvt, qraux )

  return
end
subroutine dqrlss ( a, lda, m, n, kr, b, x, rsd, jpvt, qraux )

!*****************************************************************************80
!
!! DQRLSS solves a linear system in a least squares sense.
!
!  Discussion:
!
!    DQRLSS must be preceeded by a call to DQRANK.
!
!    The system is to be solved is
!      A * X = B
!    where
!      A is an M by N matrix with rank KR, as determined by DQRANK,
!      B is a given M-vector,
!      X is the N-vector to be computed.
!
!    A solution X, with at most KR nonzero components, is found which
!    minimizes the 2-norm of the residual (A*X-B).
!
!    Once the matrix A has been formed, DQRANK should be
!    called once to decompose it.  Then, for each right hand
!    side B, DQRLSS should be called once to obtain the
!    solution and residual.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(LDA,N), the QR factorization information
!    from DQRANK.  The triangular matrix R of the QR factorization is
!    contained in the upper triangle and information needed to recover
!    the orthogonal matrix Q is stored below the diagonal in A and in
!    the vector QRAUX.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be at least M.
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, integer ( kind = 4 ) KR, the rank of the matrix, as estimated
!    by DQRANK.
!
!    Input, real ( kind = 8 ) B(M), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), a least squares solution to the
!    linear system.
!
!    Output, real ( kind = 8 ) RSD(M), the residual, B - A*X.  RSD may
!    overwite B.
!
!    Input, integer ( kind = 4 ) JPVT(N), the pivot information from DQRANK.
!    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
!    independent to within the tolerance TOL and the remaining columns
!    are linearly dependent.
!
!    Input, real ( kind = 8 ) QRAUX(N), auxiliary information from DQRANK
!    defining the QR factorization.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jpvt(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kr
  real ( kind = 8 ) qraux(n)
  real ( kind = 8 ) rsd(m)
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  if ( kr /= 0 ) then
    call dqrsl ( a, lda, m, kr, qraux, b, rsd, rsd, x, rsd, rsd, 110, info )
  end if

  jpvt(1:n) = - jpvt(1:n)

  x(kr+1:n) = 0.0D+00

  do j = 1, n

    if ( jpvt(j) <= 0 ) then

      k = -jpvt(j)
      jpvt(j) = k

      do while ( k /= j )
        t = x(j)
        x(j) = x(k)
        x(k) = t
        jpvt(k) = -jpvt(k)
        k = jpvt(k)
      end do

    end if

  end do

  return
end
subroutine dqrsl ( a, lda, n, k, qraux, y, qy, qty, b, rsd, ab, job, info )

!*****************************************************************************80
!
!! DQRSL computes transformations, projections, and least squares solutions.
!
!  Discussion:
!
!    DQRSL requires the output of DQRDC.
!
!    For K <= min(N,P), let AK be the matrix
!
!      AK = ( A(JPVT(1)), A(JPVT(2)), ..., A(JPVT(K)) )
!
!    formed from columns JPVT(1), ..., JPVT(K) of the original
!    N by P matrix A that was input to DQRDC.  If no pivoting was
!    done, AK consists of the first K columns of A in their
!    original order.  DQRDC produces a factored orthogonal matrix Q
!    and an upper triangular matrix R such that
!
!      AK = Q * (R)
!               (0)
!
!    This information is contained in coded form in the arrays
!    A and QRAUX.
!
!    The parameters QY, QTY, B, RSD, and AB are not referenced
!    if their computation is not requested and in this case
!    can be replaced by dummy variables in the calling program.
!    To save storage, the user may in some cases use the same
!    array for different parameters in the calling sequence.  A
!    frequently occuring example is when one wishes to compute
!    any of B, RSD, or AB and does not need Y or QTY.  In this
!    case one may identify Y, QTY, and one of B, RSD, or AB, while
!    providing separate arrays for anything else that is to be
!    computed.
!
!    Thus the calling sequence
!
!      call dqrsl ( a, lda, n, k, qraux, y, dum, y, b, y, dum, 110, info )
!
!    will result in the computation of B and RSD, with RSD
!    overwriting Y.  More generally, each item in the following
!    list contains groups of permissible identifications for
!    a single calling sequence.
!
!      1. (Y,QTY,B) (RSD) (AB) (QY)
!
!      2. (Y,QTY,RSD) (B) (AB) (QY)
!
!      3. (Y,QTY,AB) (B) (RSD) (QY)
!
!      4. (Y,QY) (QTY,B) (RSD) (AB)
!
!      5. (Y,QY) (QTY,RSD) (B) (AB)
!
!      6. (Y,QY) (QTY,AB) (B) (RSD)
!
!    In any group the value returned in the array allocated to
!    the group corresponds to the last member of the group.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(LDA,P), contains the output of DQRDC.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix AK.  It must
!    have the same value as N in DQRDC.
!
!    Input, integer ( kind = 4 ) K, the number of columns of the matrix AK.  K
!    must not be greater than min(N,P), where P is the same as in the
!    calling sequence to DQRDC.
!
!    Input, real ( kind = 8 ) QRAUX(P), the auxiliary output from DQRDC.
!
!    Input, real ( kind = 8 ) Y(N), a vector to be manipulated by DQRSL.
!
!    Output, real ( kind = 8 ) QY(N), contains Q * Y, if requested.
!
!    Output, real ( kind = 8 ) QTY(N), contains Q' * Y, if requested.
!
!    Output, real ( kind = 8 ) B(K), the solution of the least squares problem
!      minimize norm2 ( Y - AK * B),
!    if its computation has been requested.  Note that if pivoting was
!    requested in DQRDC, the J-th component of B will be associated with
!    column JPVT(J) of the original matrix A that was input into DQRDC.
!
!    Output, real ( kind = 8 ) RSD(N), the least squares residual Y - AK * B,
!    if its computation has been requested.  RSD is also the orthogonal
!    projection of Y onto the orthogonal complement of the column space
!    of AK.
!
!    Output, real ( kind = 8 ) AB(N), the least squares approximation Ak * B,
!    if its computation has been requested.  AB is also the orthogonal
!    projection of Y onto the column space of A.
!
!    Input, integer ( kind = 4 ) JOB, specifies what is to be computed.  JOB has
!    the decimal expansion ABCDE, with the following meaning:
!
!      if A /= 0, compute QY.
!      if B /= 0, compute QTY.
!      if C /= 0, compute QTY and B.
!      if D /= 0, compute QTY and RSD.
!      if E /= 0, compute QTY and AB.
!
!    Note that a request to compute B, RSD, or AB automatically triggers
!    the computation of QTY, for which an array must be provided in the
!    calling sequence.
!
!    Output, integer ( kind = 4 ) INFO, is zero unless the computation of B has
!    been requested and R is exactly singular.  In this case, INFO is the
!    index of the first zero diagonal element of R, and B is left unaltered.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) ab(n)
  real ( kind = 8 ) b(k)
  logical cab
  logical cb
  logical cqty
  logical cqy
  logical cr
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) kp1
  real ( kind = 8 ) qraux(*)
  real ( kind = 8 ) qty(n)
  real ( kind = 8 ) qy(n)
  real ( kind = 8 ) rsd(n)
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) y(n)
!
!  set info flag.
!
  info = 0
!
!  Determine what is to be computed.
!
  cqy =        job / 10000         /= 0
  cqty = mod ( job,  10000 )       /= 0
  cb =   mod ( job,   1000 ) / 100 /= 0
  cr =   mod ( job,    100 ) /  10 /= 0
  cab =  mod ( job,     10 )       /= 0

  ju = min ( k, n-1 )
!
!  Special action when N = 1.
!
  if ( ju == 0 ) then

    if ( cqy ) then
      qy(1) = y(1)
    end if

    if ( cqty ) then
      qty(1) = y(1)
    end if

    if ( cab ) then
      ab(1) = y(1)
     end if

    if ( cb ) then

      if ( a(1,1) == 0.0D+00 ) then
        info = 1
      else
        b(1) = y(1) / a(1,1)
      end if

    end if

    if ( cr ) then
      rsd(1) = 0.0D+00
    end if

    return

  end if
!
!  Set up to compute QY or QTY.
!
  if ( cqy ) then
    qy(1:n) = y(1:n)
  end if

  if ( cqty ) then
    qty(1:n) = y(1:n)
  end if
!
!  Compute QY.
!
  if ( cqy ) then

    do jj = 1, ju

      j = ju - jj + 1

      if ( qraux(j) /= 0.0D+00 ) then
        temp = a(j,j)
        a(j,j) = qraux(j)
        t = - ddot ( n-j+1, a(j,j), 1, qy(j), 1 ) / a(j,j)
        call daxpy ( n-j+1, t, a(j,j), 1, qy(j), 1 )
        a(j,j) = temp
      end if

    end do

  end if
!
!  Compute Q'*Y.
!
     if ( cqty ) then

        do j = 1, ju
           if ( qraux(j) /= 0.0D+00 ) then
              temp = a(j,j)
              a(j,j) = qraux(j)
              t = - ddot ( n-j+1, a(j,j), 1, qty(j), 1 ) / a(j,j)
              call daxpy ( n-j+1, t, a(j,j), 1, qty(j), 1 )
              a(j,j) = temp
           end if
        end do

     end if
!
!  Set up to compute B, RSD, or AB.
!
     if ( cb ) then
       b(1:k) = qty(1:k)
     end if

     kp1 = k + 1

     if ( cab ) then
       ab(1:k) = qty(1:k)
     end if

     if ( cr .and. k < n ) then
       rsd(k+1:n) = qty(k+1:n)
     end if

     if ( cab .and. k+1 <= n ) then
        ab(k+1:n) = 0.0D+00
     end if

     if ( cr ) then
        rsd(1:k) = 0.0D+00
     end if
!
!  Compute B.
!
     if ( cb ) then

        do jj = 1, k

           j = k - jj + 1

           if ( a(j,j) == 0.0D+00 ) then
              info = j
              exit
           end if

           b(j) = b(j)/a(j,j)

           if ( j /= 1 ) then
              t = -b(j)
              call daxpy ( j-1, t, a(1,j), 1, b, 1 )
           end if

        end do

     end if

     if ( cr .or. cab ) then
!
!  Compute RSD or AB as required.
!
        do jj = 1, ju

           j = ju - jj + 1

           if ( qraux(j) /= 0.0D+00 ) then

              temp = a(j,j)
              a(j,j) = qraux(j)

              if ( cr ) then
                 t = - ddot ( n-j+1, a(j,j), 1, rsd(j), 1 ) / a(j,j)
                 call daxpy ( n-j+1, t, a(j,j), 1, rsd(j), 1 )
              end if

              if ( cab ) then
                 t = - ddot ( n-j+1, a(j,j), 1, ab(j), 1 ) / a(j,j)
                 call daxpy ( n-j+1, t, a(j,j), 1, ab(j), 1 )
              end if

              a(j,j) = temp

           end if

        end do

  end if

  return
end
subroutine drot ( n, x, incx, y, incy, c, s )

!*****************************************************************************80
!
!! DROT applies a plane rotation.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 8 ) X(*), one of the vectors to be rotated.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries 
!    of X.
!
!    Input/output, real ( kind = 8 ) Y(*), one of the vectors to be rotated.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive elements 
!    of Y.
!
!    Input, real ( kind = 8 ) C, S, parameters (presumably the cosine and
!    sine of some angle) that define a plane rotation.
!
  implicit none

  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n
  real ( kind = 8 ) s
  real ( kind = 8 ) stemp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    do i = 1, n
      stemp = c * x(i) + s * y(i)
      y(i) = c * y(i) - s * x(i)
      x(i) = stemp
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      stemp = c * x(ix) + s * y(iy)
      y(iy) = c * y(iy) - s * x(ix)
      x(ix) = stemp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine drotg ( sa, sb, c, s )

!*****************************************************************************80
!
!! DROTG constructs a Givens plane rotation.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) SA, SB, ...
!
!    Output, real ( kind = 8 ) C, S, ...
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) r
  real ( kind = 8 ) roe
  real ( kind = 8 ) s
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) scale
  real ( kind = 8 ) z

  if ( abs ( sb ) < abs ( sa ) ) then
    roe = sa
  else
    roe = sb
  end if

  scale = abs ( sa ) + abs ( sb )

  if ( scale == 0.0D+00 ) then
    c = 1.0D+00
    s = 0.0D+00
    r = 0.0D+00
  else
    r = scale * sqrt ( ( sa / scale )**2 + ( sb / scale )**2 )
    r = sign ( 1.0D+00, roe ) * r
    c = sa / r
    s = sb / r
  end if

  if ( 0.0D+00 < abs ( c ) .and. abs ( c ) <= s ) then
    z = 1.0D+00 / c
  else
    z = s
  end if

  sa = r
  sb = z

  return
end
subroutine dscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! DSCAL scales a vector by a constant.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input/output, real ( kind = 8 ) X(*), the vector to be scaled.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) x(*)

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    x(1:m) = sa * x(1:m)

    do i = m+1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
subroutine dsftb ( n, r, azero, a, b )

!*****************************************************************************80
!
!! DSFTB computes a "slow" backward Fourier transform of real data.
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Output, real ( kind = 8 ) R(N), the reconstructed data sequence.
!
!    Input, real ( kind = 8 ) AZERO, the constant Fourier coefficient.
!
!    Input, real ( kind = 8 ) A(N/2), B(N/2), the Fourier coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n/2)
  real ( kind = 8 ) azero
  real ( kind = 8 ) b(n/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) theta

  r(1:n) = azero
  do i = 1, n
    do k = 1, n/2
      theta = real ( k * ( i - 1 ) * 2 ) * pi / real ( n, kind = 8 )
      r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
    end do
  end do

  return
end
subroutine dsftf ( n, r, azero, a, b )

!*****************************************************************************80
!
!! DSFTF computes a "slow" forward Fourier transform of real data.
!
!  Modified:
!
!    13 March 2001
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) R(N), the data to be transformed.
!
!    Output, real ( kind = 8 ) AZERO, = sum ( 1 <= I <= N ) R(I) / N.
!
!    Output, real ( kind = 8 ) A(N/2), B(N/2), the Fourier coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(1:n/2)
  real ( kind = 8 ) azero
  real ( kind = 8 ) b(1:n/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) theta

  azero = sum ( r(1:n) ) / real ( n, kind = 8 )

  do i = 1, n / 2

    a(i) = 0.0D+00
    b(i) = 0.0D+00

    do j = 1, n
      theta = real ( 2 * i * ( j - 1 ) ) * pi / real ( n, kind = 8 )
      a(i) = a(i) + r(j) * cos ( theta )
      b(i) = b(i) + r(j) * sin ( theta )
    end do

    a(i) = a(i) / real ( n, kind = 8 )
    b(i) = b(i) / real ( n, kind = 8 )

    if ( i /= ( n / 2 ) ) then
      a(i) = 2.0D+00 * a(i)
      b(i) = 2.0D+00 * b(i)
    end if

  end do

  return
end
subroutine dsvdc ( x, ldx, n, p, s, e, u, ldu, v, ldv, work, job, info )

!*****************************************************************************80
!
!! DSVDC computes the singular value decomposition of a real rectangular matrix.
!
!  Discussion:
!
!    DSVDC reduces a real ( kind = 8 ) N by P matrix X to diagonal form by
!    orthogonal transformations U and V.  The diagonal elements S(I) are
!    the singular values of X.  The columns of U are the corresponding
!    left singular vectors, and the columns of V the right singular vectors.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(LDX,P).  On input, the matrix whose
!    singular value decomposition is to be computed.  On output, the matrix
!    has been destroyed.
!
!    Input, integer ( kind = 4 ) LDX, the leading dimension of the array X.
!    LDX must be at least N.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) P, the number of columns of the matrix X.
!
!    Output, real ( kind = 8 ) S(MM), where MM = min(N+1,P).  The first min(N,P)
!    entries of S contain the singular values of X arranged in descending
!    order of magnitude.
!
!    Output, real ( kind = 8 ) E(P), ordinarily contains zeros.  However see the
!    discussion of INFO for exceptions.
!
!    Output, real ( kind = 8 ) U(LDU,K).  If JOBA = 1 then K = N; if
!    2 <= JOBA, then K = min(N,P).  U contains the matrix of left singular
!    vectors.  U is not referenced if JOBA = 0.  If N <= P or if JOBA = 2, then
!    U may be identified with X in the subroutine call.
!
!    Input, integer ( kind = 4 ) LDU, the leading dimension of the array U.
!    LDU must be at least N.
!
!    Output, real ( kind = 8 ) V(LDV,P), the matrix of right singular vectors.
!    V is not referenced if JOB is 0.  If P <= N, then V may be identified
!    with X in the subroutine call.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of the array V.
!    LDV must be at least P.
!
!    Workspace, real ( kind = 8 ) WORK(N).
!
!    Input, integer ( kind = 4 ) JOB, controls the computation of the singular
!    vectors.  It has the decimal expansion AB with the following meaning:
!
!      A =  0, do not compute the left singular vectors.
!      A =  1, return the N left singular vectors in U.
!      A >= 2, return the first min(N,P) singular vectors in U.
!      B =  0, do not compute the right singular vectors.
!      B =  1, return the right singular vectors in V.
!
!    Output, integer ( kind = 4 ) INFO, status indicator.
!    The singular values (and their corresponding singular vectors)
!    S(INFO+1), S(INFO+2),...,S(M) are correct (here M=min(N,P)).
!    Thus if INFO is 0, all the singular values and their vectors are
!    correct.  In any event, the matrix B = U'*X*V is the bidiagonal matrix
!    with the elements of S on its diagonal and the elements of E on
!    its superdiagonal (U' is the transpose of U).  Thus the singular
!    values of X and B are the same.
!
  implicit none

  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) ldx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cs
  real ( kind = 8 ) e(p)
  real ( kind = 8 ) el
  real ( kind = 8 ) emm1
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jobu
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kase
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) lls
  integer ( kind = 4 ) ls
  integer ( kind = 4 ) lu
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: maxit = 30
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mm1
  integer ( kind = 4 ) mp1
  integer ( kind = 4 ) nct
  integer ( kind = 4 ) nctp1
  integer ( kind = 4 ) ncu
  integer ( kind = 4 ) nrt
  integer ( kind = 4 ) nrtp1
  real ( kind = 8 ) s(*)
  real ( kind = 8 ) scale
  real ( kind = 8 ) ddot
  real ( kind = 8 ) shift
  real ( kind = 8 ) sl
  real ( kind = 8 ) sm
  real ( kind = 8 ) smm1
  real ( kind = 8 ) sn
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) test
  real ( kind = 8 ) u(ldu,*)
  real ( kind = 8 ) v(ldv,p)
  logical wantu
  logical wantv
  real ( kind = 8 ) work(n)
  real ( kind = 8 ) x(ldx,p)
  real ( kind = 8 ) ztest
!
!  Determine what is to be computed.
!
  wantu = .false.
  wantv = .false.
  jobu = mod ( job, 100 ) / 10

  if ( 1 < jobu ) then
    ncu = min ( n, p )
  else
    ncu = n
  end if

  if ( jobu /= 0 ) then
    wantu = .true.
  end if

  if ( mod ( job, 10 ) /= 0 ) then
    wantv = .true.
  end if
!
!  Reduce X to bidiagonal form, storing the diagonal elements
!  in S and the super-diagonal elements in E.
!
  info = 0
  nct = min ( n-1, p )
  nrt = max ( 0, min ( p-2, n ) )
  lu = max ( nct, nrt )

  do l = 1, lu
!
!  Compute the transformation for the L-th column and
!  place the L-th diagonal in S(L).
!
     if ( l <= nct ) then

        s(l) = dnrm2 ( n-l+1, x(l,l), 1 )

        if ( s(l) /= 0.0D+00 ) then
          if (x(l,l) /= 0.0D+00 ) then
            s(l) = sign ( s(l), x(l,l) )
          end if
          call dscal ( n-l+1, 1.0D+00 / s(l), x(l,l), 1 )
          x(l,l) = 1.0D+00 + x(l,l)
        end if

        s(l) = -s(l)

     end if

     do j = l+1, p
!
!  Apply the transformation.
!
        if ( l <= nct .and. s(l) /= 0.0D+00 ) then
           t = - ddot ( n-l+1, x(l,l), 1, x(l,j), 1 ) / x(l,l)
           call daxpy ( n-l+1, t, x(l,l), 1, x(l,j), 1 )
        end if
!
!  Place the L-th row of X into E for the
!  subsequent calculation of the row transformation.
!
        e(j) = x(l,j)
     end do
!
!  Place the transformation in U for subsequent back multiplication.
!
     if ( wantu .and. l <= nct ) then
        u(l:n,l) = x(l:n,l)
     end if

     if ( l <= nrt ) then
!
!  Compute the L-th row transformation and place the
!  L-th superdiagonal in E(L).
!
        e(l) = dnrm2 ( p-l, e(l+1), 1 )

        if ( e(l) /= 0.0D+00 ) then
          if ( e(l+1) /= 0.0D+00 ) then
            e(l) = sign ( e(l), e(l+1) )
          end if
          call dscal ( p-l, 1.0D+00 / e(l), e(l+1), 1 )
          e(l+1) = 1.0D+00 + e(l+1)
        end if

        e(l) = -e(l)
!
!  Apply the transformation.
!
        if ( l+1 <= n .and. e(l) /= 0.0D+00 ) then

           work(l+1:n) = 0.0D+00

           do j = l+1, p
              call daxpy ( n-l, e(j), x(l+1,j), 1, work(l+1), 1 )
           end do

           do j = l+1, p
              call daxpy ( n-l, -e(j)/e(l+1), work(l+1), 1, x(l+1,j), 1 )
           end do

        end if
!
!  Place the transformation in V for subsequent back multiplication.
!
        if ( wantv ) then
          v(l+1:p,l) = e(l+1:p)
        end if

     end if

  end do
!
!  Set up the final bidiagonal matrix of order M.
!
  m = min ( p, n+1 )
  nctp1 = nct + 1
  nrtp1 = nrt + 1

  if ( nct < p ) then
    s(nctp1) = x(nctp1,nctp1)
  end if

  if ( n < m ) then
    s(m) = 0.0D+00
  end if

  if ( nrtp1 < m ) then
    e(nrtp1) = x(nrtp1,m)
  end if

  e(m) = 0.0D+00
!
!  If required, generate U.
!
  if ( wantu ) then

    u(1:n,nctp1:ncu) = 0.0D+00

    do j = nctp1, ncu
      u(j,j) = 1.0D+00
    end do

    do ll = 1, nct

      l = nct - ll + 1

      if ( s(l) /= 0.0D+00 ) then

        do j = l+1, ncu
          t = - ddot ( n-l+1, u(l,l), 1, u(l,j), 1 ) / u(l,l)
          call daxpy ( n-l+1, t, u(l,l), 1, u(l,j), 1 )
        end do

        call dscal ( n-l+1, -1.0D+00, u(l,l), 1 )
        u(l,l) = 1.0D+00 + u(l,l)
        u(1:l-1,l) = 0.0D+00

      else

        u(1:n,l) = 0.0D+00
        u(l,l) = 1.0D+00

      end if

    end do

  end if
!
!  If it is required, generate V.
!
  if ( wantv ) then

     do ll = 1, p

        l = p - ll + 1

        if ( l <= nrt .and. e(l) /= 0.0D+00 ) then

           do j = l+1, p
              t = - ddot ( p-l, v(l+1,l), 1, v(l+1,j), 1 ) / v(l+1,l)
              call daxpy ( p-l, t, v(l+1,l), 1, v(l+1,j), 1 )
           end do

        end if

        v(1:p,l) = 0.0D+00
        v(l,l) = 1.0D+00

     end do

  end if
!
!  Main iteration loop for the singular values.
!
  mm = m
  iter = 0

  do while ( 0 < m )
!
!  If too many iterations have been performed, set flag and return.
!
     if ( maxit <= iter ) then
        info = m
        return
     end if
!
!  This section of the program inspects for
!  negligible elements in the S and E arrays.
!
!  On completion the variables KASE and L are set as follows:
!
!  KASE = 1     if S(M) and E(L-1) are negligible and L<M
!  KASE = 2     if S(L) is negligible and L<M
!  KASE = 3     if E(L-1) is negligible, L<M, and
!               S(L), ..., s(M) are not negligible (QR step).
!  KASE = 4     if E(M-1) is negligible (convergence).
!
     do ll = 1, m

        l = m - ll

        if ( l == 0 ) then
          exit
        end if

        test = abs ( s(l) ) + abs ( s(l+1) )
        ztest = test + abs ( e(l) )

        if ( ztest == test ) then
           e(l) = 0.0D+00
           exit
        end if

     end do

     if ( l == m - 1 ) then

        kase = 4

     else

        mp1 = m + 1

        do lls = l+1, m+1

           ls = m - lls + l+1

           if ( ls == l ) then
             exit
           end if

           test = 0.0D+00
           if ( ls /= m ) then
             test = test + abs ( e(ls) )
           end if

           if ( ls /= l + 1 ) then
             test = test + abs ( e(ls-1) )
           end if

           ztest = test + abs ( s(ls) )

           if ( ztest == test ) then
              s(ls) = 0.0D+00
              exit
           end if

        end do

        if ( ls == l ) then
          kase = 3
        else if ( ls == m ) then
          kase = 1
        else
          kase = 2
          l = ls
        end if

     end if

     l = l + 1
!
!  Deflate negligible S(M).
!
     if ( kase == 1 ) then

        mm1 = m - 1
        f = e(m-1)
        e(m-1) = 0.0D+00

        do kk = l, mm1

           k = mm1 - kk + l
           t1 = s(k)
           call drotg ( t1, f, cs, sn )
           s(k) = t1

           if ( k /= l ) then
             f = - sn * e(k-1)
             e(k-1) = cs * e(k-1)
           end if

           if ( wantv ) then
             call drot ( p, v(1,k), 1, v(1,m), 1, cs, sn )
           end if

        end do
!
!  Split at negligible S(L).
!
     else if ( kase == 2 ) then

        f = e(l-1)
        e(l-1) = 0.0D+00

        do k = l, m

           t1 = s(k)
           call drotg ( t1, f, cs, sn )
           s(k) = t1
           f = - sn * e(k)
           e(k) = cs * e(k)
           if ( wantu ) then
             call drot ( n, u(1,k), 1, u(1,l-1), 1, cs, sn )
           end if

        end do
!
!  Perform one QR step.
!
     else if ( kase == 3 ) then
!
!  Calculate the shift.
!
        scale = max ( abs ( s(m) ), abs ( s(m-1) ), abs ( e(m-1) ), &
                      abs ( s(l) ), abs ( e(l) ) )

        sm = s(m) / scale
        smm1 = s(m-1) / scale
        emm1 = e(m-1) / scale
        sl = s(l) / scale
        el = e(l) / scale
        b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1**2 ) / 2.0D+00
        c = ( sm * emm1 )**2
        shift = 0.0D+00

        if ( b /= 0.0D+00 .or. c /= 0.0D+00 ) then
           shift = sqrt ( b**2 + c )
           if ( b < 0.0D+00 ) then
             shift = - shift
           end if
           shift = c / ( b + shift )
        end if

        f = ( sl + sm ) * ( sl - sm ) - shift
        g = sl * el
!
!  Chase zeros.
!
        mm1 = m - 1

        do k = l, mm1

           call drotg ( f, g, cs, sn )

           if (k /= l) then
             e(k-1) = f
           end if

           f = cs * s(k) + sn * e(k)
           e(k) = cs * e(k) - sn * s(k)
           g = sn * s(k+1)
           s(k+1) = cs * s(k+1)

           if ( wantv ) then
             call drot ( p, v(1,k), 1, v(1,k+1), 1, cs, sn )
           end if

           call drotg ( f, g, cs, sn )
           s(k) = f
           f = cs * e(k) + sn * s(k+1)
           s(k+1) = -sn * e(k) + cs * s(k+1)
           g = sn * e(k+1)
           e(k+1) = cs * e(k+1)

           if ( wantu .and. k < n ) then
             call drot ( n, u(1,k), 1, u(1,k+1), 1, cs, sn )
           end if

        end do

        e(m-1) = f
        iter = iter + 1
!
!  Convergence.
!
     else if ( kase == 4 ) then
!
!  Make the singular value nonnegative.
!
        if ( s(l) < 0.0D+00 ) then
           s(l) = -s(l)
           if ( wantv ) then
             call dscal ( p, -1.0D+00, v(1,l), 1 )
           end if
        end if
!
!  Order the singular value.
!
  590   continue

        if ( l /= mm ) then

          if ( s(l) < s(l+1) ) then

            t = s(l)
            s(l) = s(l+1)
            s(l+1) = t

            if ( wantv .and. l < p ) then
              call dswap ( p, v(1,l), 1, v(1,l+1), 1 )
            end if

            if ( wantu .and. l < n ) then
              call dswap ( n, u(1,l), 1, u(1,l+1), 1 )
            end if

            l = l + 1

            go to 590

          end if

        end if

        iter = 0
        m = m - 1

    end if

  end do

  return
end
subroutine dswap ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! DSWAP interchanges two vectors.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 8 ) X(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    elements of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 3 )

    do i = 1, m
      temp = x(i)
      x(i) = y(i)
      y(i) = temp
    end do

    do i = m+1, n, 3

      temp = x(i)
      x(i) = y(i)
      y(i) = temp

      temp = x(i+1)
      x(i+1) = y(i+1)
      y(i+1) = temp

      temp = x(i+2)
      x(i+2) = y(i+2)
      y(i+2) = temp

    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      temp = x(ix)
      x(ix) = y(iy)
      y(iy) = temp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine ea ( newflg, svalue, limexp, result, abserr, epstab, ierr )

!*****************************************************************************80
!
!! EA performs extrapolation to accelerate the convergence of a sequence.
!
!  Discussion:
!
!    Given a slowly convergent sequence, this routine attempts
!    to extrapolate nonlinearly to a better estimate of the
!    sequence's limiting value, thus improving the rate of
!    convergence.
!
!    The routine is based on the epsilon algorithm of P. Wynn.
!
!    An estimate of the absolute error is also given.
!
!    The routine can be called repeatedly, using the results of
!    previous calls to efficiently compute new values.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, logical NEWFLG, is TRUE if this is the first call for this
!    data.
!
!    Input, real ( kind = 8 ) SVALUE, ?
!
!    Input, integer ( kind = 4 ) LIMEXP, the size of the epsilon table that can be
!    generated.  LIMEXP must be at least 3.
!
!    Output, real ( kind = 8 ) RESULT, the estimate of the sequence's
!    limiting value.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of the absolute error.
!
!    Input/output, real ( kind = 8 ) EPSTAB(), ?
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no error occurred.
!    Nonzero, an error occurred.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) delta1
  real ( kind = 8 ) delta2
  real ( kind = 8 ) delta3
  real ( kind = 8 ) eprn
  real ( kind = 8 ) epstab(*)
  real ( kind = 8 ) error
  real ( kind = 8 ) err1
  real ( kind = 8 ) err2
  real ( kind = 8 ) err3
  real ( kind = 8 ) e0
  real ( kind = 8 ) e1
  real ( kind = 8 ) e2
  real ( kind = 8 ) e3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ib2
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) in
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) limexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) newelm
  integer ( kind = 4 ) num
  integer ( kind = 4 ) nres
  logical newflg
  real ( kind = 8 ) relpr
  real ( kind = 8 ) res
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) ss
  real ( kind = 8 ) svalue
  real ( kind = 8 ) tol1
  real ( kind = 8 ) tol2
  real ( kind = 8 ) tol3

  if ( limexp < 3 ) then
    ierr = 1
    call xerror ( 'LIMEXP is less than 3', 1, 1 )
    return
  end if

  ierr = 0
  res3la(1) = epstab(limexp+5)
  res3la(2) = epstab(limexp+6)
  res3la(3) = epstab(limexp+7)
  result = svalue

  if ( newflg ) then
    n = 1
    nres = 0
    newflg = .false.
    epstab(n) = svalue
    abserr = abs ( result )
    go to 100
  else
    n = int ( epstab(limexp+3) )
    nres = int ( epstab(limexp+4) )
    if ( n == 2 ) then
      epstab(n) = svalue
      abserr = 6.0D+00 * abs ( result - epstab(1) )
      go to 100
    end if
  end if

  epstab(n) = svalue
  relpr = epsilon ( relpr )
  eprn = 10.0D+00 * relpr
  epstab(n+2) = epstab(n)
  newelm = ( n - 1 ) / 2
  num = n
  k1 = n

  do i = 1, newelm

    k2 = k1 - 1
    k3 = k1 - 2
    res = epstab(k1+2)
    e0 = epstab(k3)
    e1 = epstab(k2)
    e2 = res
    delta2 = e2 - e1
    err2 = abs ( delta2 )
    tol2 = max ( abs ( e2 ), abs ( e1 ) ) * relpr
    delta3 = e1 - e0
    err3 = abs ( delta3 )
    tol3 = max ( abs ( e1 ), abs ( e0 ) ) * relpr
!
!  If e0, e1 and e2 are equal to within machine accuracy,
!  convergence is assumed.
!
    if ( err2 <= tol2 .and. err3 <= tol3 ) then
      result = res
      abserr = err2 + err3
      go to 50
    end if

    if ( i /= 1 ) then

      e3 = epstab(k1)
      epstab(k1) = e1
      delta1 = e1 - e3
      err1 = abs ( delta1 )
      tol1 = max ( abs ( e1 ), abs ( e3 ) ) * relpr
!
!  If two elements are very close to each other, omit
!  a part of the table by adjusting the value of N.
!
      if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) then
        go to 20
      end if

      ss = 1.0D+00 / delta1 + 1.0D+00 / delta2 - 1.0D+00 / delta3

    else

      epstab(k1) = e1

      if ( err2 <= tol2 .or. err3 <= tol3 ) then
        go to 20
      end if

      ss = 1.0D+00 / delta2 - 1.0D+00 / delta3

    end if
!
!  Test to detect irregular behavior in the table, and
!  eventually omit a part of the table adjusting the value of N.
!
    if ( 0.1D-03 < abs ( ss * e1 ) ) then
      go to 30
    end if

   20  continue

    n = i + i - 1

    if ( nres == 0 ) then
      abserr = err2 + err3
      result = res
    else if ( nres == 1 ) then
      result = res3la(1)
    else if ( nres == 2 ) then
      result = res3la(2)
    else
      result = res3la(3)
    end if

    go to 50
!
!  Compute a new element and eventually adjust the value of RESULT.
!
   30  continue

    res = e1 + 1.0D+00 / ss
    epstab(k1) = res
    k1 = k1-2

    if ( nres == 0 ) then
      abserr = err2 + abs ( res - e2 ) + err3
      result = res
      go to 40
    else if ( nres == 1 ) then
      error = 6.0D+00 * ( abs ( res - res3la(1) ) )
    else if ( nres == 2 ) then
      error = 2.0D+00 * ( abs ( res - res3la(2) ) + abs ( res - res3la(1) ) )
    else
      error = abs ( res - res3la(3) ) + abs ( res - res3la(2) ) &
        + abs ( res - res3la(1) )
    end if

    if ( error <= 10.0D+00 * abserr ) then
      abserr = error
      result = res
    end if

   40 continue

  end do
!
!  Compute error estimate.
!
    if ( nres == 1 ) then
      abserr = 6.0D+00 * ( abs ( result - res3la(1) ) )
    else if ( nres == 2 ) then
      abserr = 2.0D+00 * abs ( result - res3la(2) ) + abs ( result - res3la(1) )
    else if ( 2 < nres ) then
      abserr = abs ( result - res3la(3) ) + abs ( result - res3la(2) ) &
        + abs ( result - res3la(1) )
    end if
!
!  Shift the table.
!
   50 continue

  if ( n == limexp ) then
    n = 2 * ( limexp / 2 ) - 1
  end if

  ib = 1
  if ( ( num / 2 ) * 2 == num ) then
    ib = 2
  end if

  ie = newelm+1

  do i = 1, ie
    ib2 = ib+2
    epstab(ib) = epstab(ib2)
    ib = ib2
  end do

  if ( num /= n ) then

    in = num - n + 1

    do i = 1, n
      epstab(i) = epstab(in)
      in = in + 1
    end do

  end if
!
!  Update RES3LA.
!
  if ( nres == 0 ) then
    res3la(1) = result
  else if ( nres == 1 ) then
    res3la(2) = result
  else if ( nres == 2 ) then
    res3la(3) = result
  else
    res3la(1) = res3la(2)
    res3la(2) = res3la(3)
    res3la(3) = result
  end if

90 continue

  abserr = max ( abserr, eprn * abs ( result ) )
  nres = nres + 1

100 continue

  n = n + 1
  epstab(limexp+3) = real ( n, kind = 8 )
  epstab(limexp+4) = real ( nres, kind = 8 )
  epstab(limexp+5) = res3la(1)
  epstab(limexp+6) = res3la(2)
  epstab(limexp+7) = res3la(3)

110 continue

  return
end
function enorm ( n, x )

!*****************************************************************************80
!
!! ENORM computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    The Euclidean norm is computed by accumulating the sum of
!    squares in three different sums.  The sums of squares for the
!    small and large components are scaled so that no overflows
!    occur.  Non-destructive underflows are permitted.  Underflows
!    and overflows do not occur in the computation of the unscaled
!    sum of squares for the intermediate components.
!
!    The definitions of small, intermediate and large components
!    depend on two constants, RDWARF and RGIANT.  The main
!    restrictions on these constants are that RDWARF**2 not
!    underflow and RGIANT**2 not overflow.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the length of the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
!
!    Output, real ( kind = 8 ) ENORM, the Euclidean norm of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) agiant
  real ( kind = 8 ) enorm
  integer ( kind = 4 ) i
  real ( kind = 8 ) rdwarf
  real ( kind = 8 ) rgiant
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) s3
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xabs
  real ( kind = 8 ) x1max
  real ( kind = 8 ) x3max

  rdwarf = sqrt ( tiny ( rdwarf ) )
  rgiant = sqrt ( huge ( rgiant ) )

  s1 = 0.0D+00
  s2 = 0.0D+00
  s3 = 0.0D+00
  x1max = 0.0D+00
  x3max = 0.0D+00
  agiant = rgiant / real ( n, kind = 8 )

  do i = 1, n

    xabs = abs ( x(i) )

    if ( xabs <= rdwarf ) then

      if ( x3max < xabs ) then
        s3 = 1.0D+00 + s3 * ( x3max / xabs )**2
        x3max = xabs
      else if ( xabs /= 0.0D+00 ) then
        s3 = s3 + ( xabs / x3max )**2
      end if

    else if ( agiant <= xabs ) then

      if ( x1max < xabs ) then
        s1 = 1.0D+00 + s1 * ( x1max / xabs )**2
        x1max = xabs
      else
        s1 = s1 + ( xabs / x1max )**2
      end if

    else

      s2 = s2 + xabs**2

    end if

  end do
!
!  Calculation of norm.
!
  if ( s1 /= 0.0D+00 ) then

    enorm = x1max * sqrt ( s1 + ( s2 / x1max ) / x1max )

  else if ( s2 /= 0.0D+00 ) then

    if ( x3max <= s2 ) then
      enorm = sqrt ( s2 * ( 1.0D+00 + ( x3max / s2 ) * ( x3max * s3 ) ) )
    else
      enorm = sqrt ( x3max * ( ( s2 / x3max ) + ( x3max * s3 ) ) )
    end if

  else

    enorm = x3max * sqrt ( s3 )

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
!      ERROR_F(X) = ( 2 / sqrt ( PI )
!        * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
!
!    In Mathematica, the function can be evaluated by:
!
!      Erf[x]
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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
function error_f ( x )

!*****************************************************************************80
!
!! ERROR_F computes the error function.
!
!  Definition:
!
!    ERROR_F(X) = ( 2 / SQRT ( PI ) )
!      * Integral ( 0 <= T <= X ) EXP ( -T**2 ) dT
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    Wayne Fullerton
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) ERROR_F, the value of the error function at X.
!
  implicit none

  real ( kind = 8 ) csevl
  real ( kind = 8 ) error_f
  real ( kind = 8 ) error_fc
  real ( kind = 8 ), parameter, dimension ( 13 ) :: erfcs = (/ &
    -0.049046121234691808D+00, -0.14226120510371364D+00, &
     0.010035582187599796D+00, -0.000576876469976748D+00, &
     0.000027419931252196D+00, -0.000001104317550734D+00, &
     0.000000038488755420D+00, -0.000000001180858253D+00, &
     0.000000000032334215D+00, -0.000000000000799101D+00, &
     0.000000000000017990D+00, -0.000000000000000371D+00, &
     0.000000000000000007D+00 /)
  integer ( kind = 4 ) inits
  integer ( kind = 4 ), save :: nterf = 0
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ), parameter :: sqrtpi = 1.7724538509055160D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xbig = 0.0D+00
  real ( kind = 8 ) y
!
!  Initialize the Chebyshev series.
!
  if ( nterf == 0 ) then
    nterf = inits ( erfcs, 13, 0.1D+00 * epsilon ( erfcs ) )
    xbig = sqrt ( - log ( sqrtpi * epsilon ( xbig ) ) )
    sqeps = sqrt ( 2.0D+00 * epsilon ( sqeps ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    error_f = 2.0D+00 * x / sqrtpi
  else if ( y <= 1.0D+00 ) then
    error_f = x * ( 1.0D+00 + csevl ( 2.0D+00 * x**2 - 1.0D+00, erfcs, nterf ) )
  else if ( y <= xbig ) then
    error_f = sign ( 1.0D+00 - error_fc ( y ), x )
  else
    error_f = sign ( 1.0D+00, x )
  end if

  return
end
function error_fc ( x )

!*****************************************************************************80
!
!! ERROR_FC computes the complementary error function.
!
!  Modified:
!
!    26 August 2001
!
!  Author:
!
!    Wayne Fullerton
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) ERROR_FC, the value of the complementary
!    error function at X.
!
  implicit none

  real ( kind = 8 ) csevl
  real ( kind = 8 ) error_fc
  real ( kind = 8 ), parameter, dimension ( 13 ) :: erfcs = (/ &
   -0.049046121234691808D+00, -0.14226120510371364D+00, &
    0.010035582187599796D+00, -0.000576876469976748D+00, &
    0.000027419931252196D+00, -0.000001104317550734D+00, &
    0.000000038488755420D+00, -0.000000001180858253D+00, &
    0.000000000032334215D+00, -0.000000000000799101D+00, &
    0.000000000000017990D+00, -0.000000000000000371D+00, &
    0.000000000000000007D+00 /)
  real ( kind = 8 ), parameter, dimension ( 24 ) :: erfccs = (/ &
    0.0715179310202925D+00, &
   -0.026532434337606719D+00, &
    0.001711153977920853D+00, &
   -0.000163751663458512D+00, &
    0.000019871293500549D+00, &
   -0.000002843712412769D+00, &
    0.000000460616130901D+00, &
   -0.000000082277530261D+00, &
    0.000000015921418724D+00, &
   -0.000000003295071356D+00, &
    0.000000000722343973D+00, &
   -0.000000000166485584D+00, &
    0.000000000040103931D+00, &
   -0.000000000010048164D+00, &
    0.000000000002608272D+00, &
   -0.000000000000699105D+00, &
    0.000000000000192946D+00, &
   -0.000000000000054704D+00, &
    0.000000000000015901D+00, &
   -0.000000000000004729D+00, &
    0.000000000000001432D+00, &
   -0.000000000000000439D+00, &
    0.000000000000000138D+00, &
   -0.000000000000000048D+00 /)
  real ( kind = 8 ), parameter, dimension ( 23 ) :: erc2cs = (/ &
   -0.069601346602309501D+00, &
   -0.041101339362620893D+00, &
    0.003914495866689626D+00, &
   -0.000490639565054897D+00, &
    0.000071574790013770D+00, &
   -0.000011530716341312D+00, &
    0.000001994670590201D+00, &
   -0.000000364266647159D+00, &
    0.000000069443726100D+00, &
   -0.000000013712209021D+00, &
    0.000000002788389661D+00, &
   -0.000000000581416472D+00, &
    0.000000000123892049D+00, &
   -0.000000000026906391D+00, &
    0.000000000005942614D+00, &
   -0.000000000001332386D+00, &
    0.000000000000302804D+00, &
   -0.000000000000069666D+00, &
    0.000000000000016208D+00, &
   -0.000000000000003809D+00, &
    0.000000000000000904D+00, &
   -0.000000000000000216D+00, &
    0.000000000000000052D+00 /)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) inits
  integer ( kind = 4 ), save :: nterc2 = 0
  integer ( kind = 4 ), save :: nterf = 0
  integer ( kind = 4 ), save :: nterfc = 0
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ), parameter :: sqrtpi = 1.7724538509055160D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xmax = 0.0D+00
  real ( kind = 8 ), save :: xsml = 0.0D+00
  real ( kind = 8 ) y

  if ( nterf == 0 ) then

    eta = 0.1D+00 * epsilon ( eta )
    nterf = inits ( erfcs, 13, eta )
    nterfc = inits ( erfccs, 24, eta )
    nterc2 = inits ( erc2cs, 23, eta )

    xsml = -sqrt ( - log ( sqrtpi * epsilon ( xsml ) ) )
    xmax = sqrt ( - log ( sqrtpi * tiny ( xmax ) ) )
    xmax = xmax - 0.5D+00 * log ( xmax ) / xmax - 0.01D+00
    sqeps = sqrt ( 2.0D+00 * epsilon ( sqeps ) )

  end if

  if ( x <= xsml ) then
    error_fc = 2.0D+00
    return
  end if
!
!  X so big that ERROR_FC will underflow.
!
  if ( xmax < x ) then
    error_fc = 0.0D+00
    return
  end if

  y = abs ( x )
!
!  error_fc(x) = 1.0D+00 - error_f(x) for -1 <= x <= 1.
!
  if ( y <= 1.0D+00 ) then

    if ( y < sqeps ) then
      error_fc = 1.0D+00 - 2.0D+00 * x / sqrtpi
    else if ( sqeps <= y ) then
      error_fc = 1.0D+00 - x * ( 1.0D+00 + &
        csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf ) )
    end if

    return

  end if
!
!  For 1 < |x| <= xmax, error_fc(x) = 1.0D+00 - error_f(x)
!
  y = y * y

  if ( y <= 4.0D+00 ) then
    error_fc = exp ( -y ) / abs ( x ) * ( 0.5D+00 &
      + csevl ( ( 8.0D+00 / y - 5.0D+00 ) / 3.0D+00, erc2cs, nterc2 ) )
  else
    error_fc = exp ( -y ) / abs ( x ) * ( 0.5D+00 &
      + csevl ( 8.0D+00 / y - 1.0D+00, erfccs, nterfc ) )
  end if

  if ( x < 0.0D+00 ) then
    error_fc = 2.0D+00 - error_fc
  end if

  return
end
subroutine ezfftb ( n, r, azero, a, b, wsave )

!*****************************************************************************80
!
!! EZFFTB computes a real periodic sequence from its Fourier coefficients.
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    EZFFTB is a simplified but slower version of DFFTB.
!
!    The transform is defined by:
!
!      R(I) = AZERO + sum ( 1 <= K <= N/2 )
!
!          A(K) * cos ( K * ( I - 1 ) * 2 * PI / N )
!        + B(K) * sin ( K * ( I - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the output array.  The
!    method is more efficient when N is the product of small primes.
!
!    Output, real ( kind = 8 ) R(N), the reconstructed data sequence.
!
!    Input, real ( kind = 8 ) AZERO, the constant Fourier coefficient.
!
!    Input, real ( kind = 8 ) A(N/2), B(N/2), the Fourier coefficients.
!
!    Input, real ( kind = 8 ) WSAVE(3*N+15), a work array.  The WSAVE array
!    must be initialized by calling EZFFFTI.  A different WSAVE array must
!    be used for each different value of N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n/2)
  real ( kind = 8 ) azero
  real ( kind = 8 ) b(n/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) wsave(3*n+15)

  if ( n < 2 ) then

    r(1) = azero

  else if ( n == 2 ) then

    r(1) = azero + a(1)
    r(2) = azero - a(1)

  else

    ns2 = ( n - 1 ) / 2

    do i = 1, ns2
      r(2*i) = 0.5D+00 * a(i)
      r(2*i+1) = -0.5D+00 * b(i)
    end do

    r(1) = azero

    if ( mod ( n, 2 ) == 0 ) then
      r(n) = a(ns2+1)
    end if

    call dfftb ( n, r, wsave(n+1) )

  end if

  return
end
subroutine ezfftf ( n, r, azero, a, b, wsave )

!*****************************************************************************80
!
!! EZFFTF computes the Fourier coefficients of a real periodic sequence.
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    EZFFTF is a simplified but slower version of DFFTF.
!
!    The transform is defined by:
!
!      AZERO = sum ( 1 <= I <= N ) R(I) / N,
!
!    and, for K = 1 to (N-1)/2,
!
!      A(K) = sum ( 1 <= I <= N )
!        ( 2 / N ) * R(I) * cos ( K * ( I - 1 ) * 2 * PI / N )
!
!    and, if N is even, then
!
!      A(N/2) = sum ( 1 <= I <= N ) (-1) **(I-1) * R(I) / N
!
!    For K = 1 to (N-1)/2,
!
!      B(K) = sum ( 1 <= I <= N )
!        ( 2 / N ) * R(I) * sin ( K * ( I - 1 ) * 2 * PI / N )
!
!    and, if N is even, then
!
!      B(N/2) = 0.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input, real ( kind = 8 ) R(N), the sequence to be transformed.
!
!    Input, real ( kind = 8 ) WSAVE(3*N+15), a work array.  The WSAVE array
!    must be initialized by calling EZFFTI.  A different WSAVE array must
!    be used for each different value of N.
!
!    Output, real ( kind = 8 ) AZERO, the constant Fourier coefficient.
!
!    Output, real ( kind = 8 ) A(N/2), B(N/2), the Fourier coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n/2)
  real ( kind = 8 ) azero
  real ( kind = 8 ) b(n/2)
  real ( kind = 8 ) cf
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) wsave(3*n+15)

  if ( n < 2 ) then

    azero = r(1)

  else if ( n == 2 ) then

    azero = 0.5D+00 * ( r(1) + r(2) )
    a(1) = 0.5D+00 * ( r(1) - r(2) )

  else

    wsave(1:n) = r(1:n)

    call dfftf ( n, wsave(1), wsave(n+1) )

    cf = 2.0D+00 / real ( n, kind = 8 )
    azero = 0.5D+00 * cf * wsave(1)
    ns2 = ( n + 1 ) / 2

    do i = 1, ns2-1
      a(i) = cf * wsave(2*i)
      b(i) = -cf * wsave(2*i+1)
    end do

    if ( mod ( n, 2 ) /= 1 ) then
      a(ns2) = 0.5D+00 * cf * wsave(n)
      b(ns2) = 0.0D+00
    end if

  end if

  return
end
subroutine ezffti ( n, wsave )

!*****************************************************************************80
!
!! EZFFTI initializes WSAVE, used in EZFFTF and EZFFTB.
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Output, real ( kind = 8 ) WSAVE(3*N+15), contains data, dependent on
!    the value of N, which is necessary for the EZFFTF or EZFFTB routines.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) wsave(3*n+15)

  if ( n <= 1 ) then
    return
  end if

  call ezffti1 ( n, wsave(2*n+1), wsave(3*n+1) )

  return
end
subroutine ezffti1 ( n, wa, ifac )

!*****************************************************************************80
!
!! EZFFTI1 is a lower level routine used by EZFFTI.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.
!
!    Output, real ( kind = 8 ) WA(N).
!
!    Input, integer ( kind = 4 ) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg1
  real ( kind = 8 ) argh
  real ( kind = 8 ) ch1
  real ( kind = 8 ) ch1h
  real ( kind = 8 ) dch1
  real ( kind = 8 ) dsh1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ifac(15)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) nf
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sh1
  real ( kind = 8 ) wa(n)

  call i8_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0D+00 * pi / real ( n, kind = 8 )
  is = 0
  l1 = 1

  do k1 = 1, nf-1

    ip = ifac(k1+2)
    l2 = l1 * ip
    ido = n / l2
    arg1 = real ( l1, kind = 8 ) * argh
    ch1 = 1.0D+00
    sh1 = 0.0D+00
    dch1 = cos ( arg1 )
    dsh1 = sin ( arg1 )

    do j = 1, ip-1

      ch1h = dch1 * ch1 - dsh1 * sh1
      sh1  = dch1 * sh1 + dsh1 * ch1
      ch1 = ch1h
      i = is + 2
      wa(i-1) = ch1
      wa(i) = sh1

      do ii = 5, ido, 2
        i = i + 2
        wa(i-1) = ch1 * wa(i-3) - sh1 * wa(i-2)
        wa(i)   = ch1 * wa(i-2) + sh1 * wa(i-3)
      end do

      is = is + ido

    end do

    l1 = l2

  end do

  return
end
subroutine fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )

!*****************************************************************************80
!
!! FDJAC1 estimates an N by N jacobian matrix using forward differences.
!
!  Discussion:
!
!    This subroutine computes a forward-difference approximation
!    to the N by N jacobian matrix associated with a specified
!    problem of N functions in N variables. If the jacobian has
!    a banded form, then function evaluations are saved by only
!    approximating the nonzero terms.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!
!      subroutine fcn ( n, x, fvec, iflag )
!
!      integer ( kind = 4 ) n
!
!      real fvec(n)
!      integer ( kind = 4 ) iflag
!      real x(n)
!
!    The value of IFLAG should not be changed by FCN unless
!    the user wants to terminate execution of the routine.
!    In this case set IFLAG to a negative integer ( kind = 4 ).
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the jacobian is evaluated.
!
!    Input, real ( kind = 8 ) FVEC(N), the functions evaluated at X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), the N by N approximate
!    jacobian matrix.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC, which must
!    not be less than N.
!
!    Output, integer ( kind = 4 ) IFLAG, is an error flag returned by FCN.  If FCN
!    returns a nonzero value of IFLAG, then this routine returns immediately
!    to the calling program, with the value of IFLAG.
!
!    Input, integer ( kind = 4 ) ML, MU, specify the number of subdiagonals and
!    superdiagonals within the band of the jacobian matrix.  If the
!    jacobian is not banded, set ML and MU to N-1.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable step
!    length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed
!    that the relative errors in the functions are of the order of the machine
!    precision.
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) n

  real ( kind = 8 ) eps
  real ( kind = 8 ) epsfcn
  real ( kind = 8 ) epsmch
  external fcn
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(n)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) msum
  integer ( kind = 4 ) mu
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa1(n)
  real ( kind = 8 ) wa2(n)
  real ( kind = 8 ) x(n)

  epsmch = epsilon ( epsmch )

  eps = sqrt ( max ( epsfcn, epsmch ) )
  msum = ml + mu + 1
!
!  Computation of dense approximate jacobian.
!
  if ( n <= msum ) then

     do j = 1, n

       temp = x(j)
       h = eps * abs ( temp )
       if ( h == 0.0D+00 ) then
         h = eps
       end if
       x(j) = temp + h
       call fcn ( n, x, wa1, iflag )

       if ( iflag < 0 ) then
         exit
       end if

       x(j) = temp
       fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h

     end do

  else
!
!  Computation of banded approximate jacobian.
!
     do k = 1, msum

        do j = k, n, msum
          wa2(j) = x(j)
          h = eps * abs ( wa2(j) )
          if ( h == 0.0D+00 ) then
            h = eps
          end if
          x(j) = wa2(j) + h
        end do

        call fcn ( n, x, wa1, iflag )

        if ( iflag < 0 ) then
          exit
        end if

        do j = k, n, msum

          x(j) = wa2(j)
          h = eps * abs ( wa2(j) )
          if ( h == 0.0D+00 ) then
            h = eps
          end if

          fjac(1:n,j) = 0.0D+00

          do i = 1, n
            if ( j - mu <= i .and. i <= j + ml ) then
              fjac(i,j) = ( wa1(i) - fvec(i) ) / h
            end if
          end do

        end do

     end do

  end if

  return
end
function fmin ( a, b, f, tol )

!*****************************************************************************80
!
!! FMIN seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    FMIN seeks an approximation to the point where F attains a minimum
!    value (strictly) inside the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324....
!
!    The function F is never evaluated at two points closer together
!    than EPS * ABS ( XSTAR ) + (TOL/3), where EPS is approximately the
!    square root of the relative machine precision.  If F is a unimodal
!    function and the computed values of F are always unimodal when
!    separated by at least EPS * ABS ( XSTAR ) + (TOL/3), then FMIN
!    approximates the abcissa of the global minimum of F on the
!    interval [A, B] with an error less than 3 * EPS * ABS ( FMIN ) + TOL.
!    If F is not unimodal, then FMIN may approximate a local, but
!    perhaps non-global, minimum to the same accuracy.
!
!    It is worth stating explicitly that this routine will NOT be
!    able to detect a minimizer that occurs at either initial endpoint
!    A or B.  If this is a concern to the user, then the user must
!    either ensure that the initial interval is larger, or to check
!    the function value at the returned minimizer against the values
!    at either endpoint.
!
!  Modified:
!
!    16 October 2004
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters
!
!    Input/output, real ( kind = 8 ) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper
!    bounds for the minimizer.  It is required that A < B.
!
!    Input, external F, a real function of the form
!      function f ( x )
!      real f
!      real x
!    which evaluates F(X) for any X in the interval [A,B].
!
!    Input, real ( kind = 8 ) TOL, the desired length of the interval of
!    uncertainty of the final result.  TOL must not be negative.
!
!    Output, real ( kind = 8 ) FMIN, the approximate minimizer of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fmin
  real ( kind = 8 ) fu
  real ( kind = 8 ) fv
  real ( kind = 8 ) fw
  real ( kind = 8 ) fx
  real ( kind = 8 ) midpoint
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) tol
  real ( kind = 8 ) tol1
  real ( kind = 8 ) tol2
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x

  if ( b <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FMIN - Fatal error!'
    write ( *, '(a)' ) '  A < B is required, but'
    write ( *, '(a,g14.6)' ) '  A = ', a
    write ( *, '(a,g14.6)' ) '  B = ', b
    stop
  end if

  c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )
!
!  C is the squared inverse of the golden ratio.
!
!  EPS is the square root of the relative machine precision.
!
  eps = sqrt ( epsilon ( eps ) )
!
!  Initialization.
!
  v = a + c * ( b - a )
  w = v
  x = v
  e = 0.0D+00
  fx = f(x)
  fv = fx
  fw = fx
!
!  The main loop starts here.
!
  do

    midpoint = 0.5D+00 * ( a + b )
    tol1 = eps * abs ( x ) + tol / 3.0D+00
    tol2 = 2.0D+00 * tol1
!
!  Check the stopping criterion.
!
    if ( abs ( x - midpoint ) <= ( tol2 - 0.5D+00 * ( b - a ) ) ) then
      exit
    end if
!
!  Is golden-section necessary?
!
    if ( abs ( e ) <= tol1 ) then
      if ( midpoint <= x ) then
        e = a - x
      else
        e = b - x
      end if

      d = c * e
!
!  Consider fitting a parabola.
!
    else

      r = ( x - w ) * ( fx - fv )
      q = ( x - v ) * ( fx - fw )
      p = ( x - v ) * q - ( x - w ) * r
      q = 2.0D+00 * ( q - r )
      if ( 0.0D+00 < q ) then
        p = -p
      end if
      q = abs ( q )
      r = e
      e = d
!
!  Choose a golden-section step if the parabola is not advised.
!
      if ( &
        ( abs ( 0.5D+00 * q * r ) <= abs ( p ) ) .or. &
        ( p <= q * ( a - x ) ) .or. &
        ( q * ( b - x ) <= p ) ) then

        if ( midpoint <= x ) then
          e = a - x
        else
          e = b - x
        end if

        d = c * e
!
!  Choose a parabolic interpolation step.
!
      else

        d = p / q
        u = x + d

        if ( ( u - a ) < tol2 ) then
          d = sign ( tol1, midpoint - x )
        end if

        if ( ( b - u ) < tol2 ) then
          d = sign ( tol1, midpoint - x )
        end if

      end if

    end if
!
!  F must not be evaluated too close to X.
!
    if ( tol1 <= abs ( d ) ) then
      u = x + d
    end if

    if ( abs ( d ) < tol1 ) then
      u = x + sign ( tol1, d )
    end if

    fu = f(u)
!
!  Update the data.
!
    if ( fu <= fx ) then

      if ( x <= u ) then
        a = x
      else
        b = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      cycle

    end if

    if ( u < x ) then
      a = u
    else
      b = u
    end if

    if ( fu <= fw .or. w == x ) then
      v = w
      fv = fw
      w = u
      fw = fu
    else if ( fu <= fv .or. v == x .or. v == w ) then
      v = u
      fv = fu
    end if

  end do

  fmin = x

  return
end
subroutine fmin_rc ( a, b, arg, status, value )

!*****************************************************************************80
!
!! FMIN_RC seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    FMIN_RC seeks an approximation to the point where F attains a minimum on
!    the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324....
!
!    The routine is a revised version of the Brent FMIN algorithm,
!    which now uses reverse communication.
!
!    It is worth stating explicitly that this routine will NOT be
!    able to detect a minimizer that occurs at either initial endpoint
!    A or B.  If this is a concern to the user, then the user must
!    either ensure that the initial interval is larger, or to check
!    the function value at the returned minimizer against the values
!    at either endpoint.
!
!  Modified:
!
!    22 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters
!
!    Input/output, real ( kind = 8 ) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper
!    bounds for an interval containing the minimizer.  It is required
!    that A < B.
!
!    Output, real ( kind = 8 ) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
!    estimate for the function minimizer.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between the user
!    and the routine.  The user only sets STATUS to zero on the first call,
!    to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated minimizer.
!
!    Input, real ( kind = 8 ) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
!  Local parameters:
!
!    C is the squared inverse of the golden ratio.
!
!    EPS is the square root of the relative machine precision.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ), save :: c
  real ( kind = 8 ), save :: d
  real ( kind = 8 ), save :: e
  real ( kind = 8 ), save :: eps
  real ( kind = 8 ), save :: fu
  real ( kind = 8 ), save :: fv
  real ( kind = 8 ), save :: fw
  real ( kind = 8 ), save :: fx
  real ( kind = 8 ), save :: midpoint
  real ( kind = 8 ), save :: p
  real ( kind = 8 ), save :: q
  real ( kind = 8 ), save :: r
  integer ( kind = 4 ) status
  real ( kind = 8 ), save :: tol
  real ( kind = 8 ), save :: tol1
  real ( kind = 8 ), save :: tol2
  real ( kind = 8 ), save :: u
  real ( kind = 8 ), save :: v
  real ( kind = 8 ) value
  real ( kind = 8 ), save :: w
  real ( kind = 8 ), save :: x
!
!  STATUS (INPUT) = 0, startup.
!
  if ( status == 0 ) then

    if ( b <= a ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FMIN_RC - Fatal error!'
      write ( *, '(a)' ) '  A < B is required, but'
      write ( *, '(a,g14.6)' ) '  A = ', a
      write ( *, '(a,g14.6)' ) '  B = ', b
      status = -1
      stop
    end if

    c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

    eps = sqrt ( epsilon ( eps ) )
    tol = epsilon ( tol )

    v = a + c * ( b - a )
    w = v
    x = v
    e = 0.0D+00

    status = 1
    arg = x

    return
!
!  STATUS (INPUT) = 1, return with initial function value of FX.
!
  else if ( status == 1 ) then

    fx = value
    fv = fx
    fw = fx
!
!  STATUS (INPUT) = 2 or more, update the data.
!
  else if ( 2 <= status ) then

    fu = value

    if ( fu <= fx ) then

      if ( x <= u ) then
        a = x
      else
        b = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu

    else

      if ( u < x ) then
        a = u
      else
        b = u
      end if

      if ( fu <= fw .or. w == x ) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if ( fu <= fv .or. v == x .or. v == w ) then
        v = u
        fv = fu
      end if

    end if

  end if
!
!  Take the next step.
!
  midpoint = 0.5D+00 * ( a + b )
  tol1 = eps * abs ( x ) + tol / 3.0D+00
  tol2 = 2.0D+00 * tol1
!
!  If the stopping criterion is satisfied, we can exit.
!
  if ( abs ( x - midpoint ) <= ( tol2 - 0.5D+00 * ( b - a ) ) ) then
    status = 0
    return
  end if
!
!  Is golden-section necessary?
!
  if ( abs ( e ) <= tol1 ) then
    if ( midpoint <= x ) then
      e = a - x
    else
      e = b - x
    end if

    d = c * e
!
!  Consider fitting a parabola.
!
  else

    r = ( x - w ) * ( fx - fv )
    q = ( x - v ) * ( fx - fw )
    p = ( x - v ) * q - ( x - w ) * r
    q = 2.0D+00 * ( q - r )
    if ( 0.0D+00 < q ) then
      p = -p
    end if
    q = abs ( q )
    r = e
    e = d
!
!  Choose a golden-section step if the parabola is not advised.
!
    if ( &
      ( abs ( 0.5D+00 * q * r ) <= abs ( p ) ) .or. &
      ( p <= q * ( a - x ) ) .or. &
      ( q * ( b - x ) <= p ) ) then

      if ( midpoint <= x ) then
        e = a - x
      else
        e = b - x
      end if

      d = c * e
!
!  Choose a parabolic interpolation step.
!
    else

      d = p / q
      u = x + d

      if ( ( u - a ) < tol2 ) then
        d = sign ( tol1, midpoint - x )
      end if

      if ( ( b - u ) < tol2 ) then
        d = sign ( tol1, midpoint - x )
      end if

    end if

  end if
!
!  F must not be evaluated too close to X.
!
  if ( tol1 <= abs ( d ) ) then
    u = x + d
  end if

  if ( abs ( d ) < tol1 ) then
    u = x + sign ( tol1, d )
  end if
!
!  Request value of F(U).
!
  arg = u
  status = status + 1

  return
end
subroutine forslv ( nr, n, a, x, b )

!*****************************************************************************80
!
!! FORSLV solves A*x=b where A is lower triangular matrix.
!
!  Discussion:
!
!    If B is no longer required by calling routine,
!    then vectors B and X may share the same storage.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) A(NR,N), the N by N lower triangular matrix.
!
!    Output, real ( kind = 8 ) X(N), the solution.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  x(1) = b(1) / a(1,1)

  do i = 2, n
    x(i) = ( b(i) - dot_product ( a(i,1:i-1), x(1:i-1) ) ) / a(i,i)
  end do

  return
end
subroutine fstocd ( n, x, fcn, sx, rnoise, g )

!*****************************************************************************80
!
!! FSTOCD approximates the gradient of a function using central differences.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the point at which the gradient is to
!    be approximated.
!
!    Input, external FCN, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real x(n)
!      real f
!
!    Input, real ( kind = 8 ) SX(N), the scaling factors for X.
!
!    Input, real ( kind = 8 ) RNOISE, the relative noise in FCN [F(X)].
!
!    Output, real ( kind = 8 ) G(N), a central difference approximation
!    to the gradient.
!
  implicit none

  integer ( kind = 4 ) n

  external fcn
  real ( kind = 8 ) fminus
  real ( kind = 8 ) fplus
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) rnoise
  real ( kind = 8 ) stepi
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) third
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtempi
!
!  Find I-th stepsize, evaluate two neighbors in direction of I-th
!  unit vector, and evaluate I-th component of gradient.
!
  third = 1.0D+00 / 3.0D+00

  do i = 1, n
    stepi = rnoise**third * max ( abs ( x(i) ), 1.0D+00 / sx(i) )
    xtempi = x(i)
    x(i) = xtempi + stepi
    call fcn ( n, x, fplus )
    x(i) = xtempi - stepi
    call fcn ( n, x, fminus )
    x(i) = xtempi
    g(i) = ( fplus - fminus ) / ( 2.0D+00 * stepi )
  end do

  return
end
subroutine fstofd ( nr, m, n, xpls, fcn, fpls, a, sx, rnoise, fhat, icase )

!*****************************************************************************80
!
!! FSTOFD approximates a derivative by a first order approximation.
!
!  Discussion:
!
!    The routine finds the first order forward finite difference
!    approximation A to the first derivative of the function defined
!    by the subprogram "fname" evaluated at the new iterate "xpls".
!
!    For optimization, use this routine to estimate:
!
!    * the first derivative (gradient) of the optimization function "fcn"
!      if no analytic user routine has been supplied;
!
!    * the second derivative (hessian) of the optimization function
!      if no analytic user routine has been supplied for the hessian but
!      one has been supplied for the gradient ("fcn") and if the
!      optimization function is inexpensive to evaluate.
!
!    m=1 (optimization) algorithm estimates the gradient of the function
!    (fcn).   fcn(x) # f: r(n)-->r(1)
!
!    m=n (systems) algorithm estimates the jacobian of the function
!    fcn(x) # f: r(n)-->r(n).
!
!    m=n (optimization) algorithm estimates the hessian of the optimization
!    function, where the hessian is the first derivative of "fcn"
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix A.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A, and the dimension
!    of the problem.
!
!    Input, real ( kind = 8 ) XPLS(N), the point at which the derivative is
!    to be estimated.
!
!    Input, external FCN, the name of the subroutine to evaluate
!    the function, of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real x(n)
!      real f
!
!    Input, real ( kind = 8 ) FPLS(M).
!    If M is 1, (optimization), then this is the function value at
!    the new iterate.
!    If M = N for optimization, then this is the value of the first
!    derivative of the function.
!    If M = N for nonlinear systems, then this is the value
!    of the associated minimization function.
!
!    Output, real ( kind = 8 ) A(NR,N), the N by N finite difference
!    approximation.  Only the lower triangular matrix and diagonal are
!    computed.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) RNOISE, the relative noise or inaccuracy in the
!    function value FCN.
!
!    Workspace, real ( kind = 8 ) FHAT(M).
!
!    Input, integer ( kind = 4 ) ICASE, problem specifier:
!    1, optimization (gradient)
!    2, systems
!    3, optimization (hessian)
!
!  Local variables:
!
!    real STEPSZ, the stepsize in the J-th variable direction
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  external fcn
  real ( kind = 8 ) fhat(m)
  real ( kind = 8 ) fpls(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icase
  integer ( kind = 4 ) j
  real ( kind = 8 ) rnoise
  real ( kind = 8 ) stepsz
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) xpls(n)
  real ( kind = 8 ) xtmpj
!
!  Find the J-th column of A.
!  Each column is the derivative of f(fcn) with respect to xpls(j).
!
  do j = 1, n
    stepsz = sqrt ( rnoise ) * max ( abs ( xpls(j) ), 1.0D+00 / sx(j) )
    xtmpj = xpls(j)
    xpls(j) = xtmpj + stepsz
    call fcn ( n, xpls, fhat )
    xpls(j) = xtmpj
    a(1:m,j) = ( fhat(1:m) - fpls(1:m) ) / stepsz
  end do

  if ( icase /= 3 ) then
    return
  end if
!
!  If computing the hessian, A must be symmetric.
!
  do j = 1, n-1
    do i = j+1, m
      a(i,j) = ( a(i,j) + a(j,i) ) / 2.0D+00
    end do
  end do

  return
end
subroutine fzero ( f, b, c, r, re, ae, iflag )

!*****************************************************************************80
!
!! FZERO searches for a zero of a function F(X) in a given interval.
!
!  Discussion:
!
!    FZERO searches for a zero of a function F(X) between
!    the given values B and C until the width of the interval
!    (B,C) has collapsed to within a tolerance specified by
!    the stopping criterion, abs ( B - C ) <= 2 * ( RW * abs ( B ) + AE ).
!    The method used is an efficient combination of bisection
!    and the secant rule.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Lawrence Shampine, H A Watts,
!    FZERO, a Root-Solving Code,
!    sc-tm-70-631, September 1970.
!
!    TJ Dekker,
!    Finding a Zero by Means of Successive Linear Interpolation,
!    'Constructive Aspects of the Fundamental Theorem of Algebra',
!    edited by B. Dejon, P. Henrici, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external F, the name of the function.
!    This name must be in an external statement in the calling
!    program.   F must be a function of one real argument.
!
!    Input/output, real ( kind = 8 ) B, one end of the interval (B,C).  The
!    value returned for B usually is the better approximation to
!    a zero of F.
!
!    Input/output, real ( kind = 8 ) C, the other end of the interval (B,C).
!
!    Input, real ( kind = 8 ) R, a (better) guess of a zero of F which could
!    help in speeding up convergence.  If F(B) and F(R) have opposite signs,
!    a root will be found in the interval (B,R); if not, but F(R) and F(C)
!    have opposite signs, a root will be found in the interval (R,C);
!    otherwise, the interval (B,C) will be searched for a possible root.
!    When no better guess is known, it is recommended that R be set to B or C;
!    because if R is not interior to the interval (B,C), it will be ignored.
!
!    Input, real ( kind = 8 ) RE, the relative error used for RW in the
!    stopping criterion.  If the input RE is less than machine precision,
!    then RW is set to approximately machine precision.
!
!    Input, real ( kind = 8 ) AE, the absolute error used in the stopping
!    criterion.  If the given interval (B,C) contains the origin, then a
!    nonzero value should be chosen for AE.
!
!    Output, integer ( kind = 4 ) IFLAG, a status code.  The user must check IFLAG
!    after each call.  Control returns to the user in all cases.
!
!    1, B is within the requested tolerance of a zero.
!      the interval (b,c) collapsed to the requested
!      tolerance, the function changes sign in (b,c), and
!      f(x) decreased in magnitude as (b,c) collapsed.
!
!    2, F(B) = 0.  however, the interval (b,c) may not have
!      collapsed to the requested tolerance.
!
!    3, B may be near a singular point of f(x).
!      the interval (b,c) collapsed to the requested tolerance and
!      the function changes sign in (b,c), but
!      f(x) increased in magnitude as (b,c) collapsed,i.e.
!      max ( ABS ( f(b in) ), ABS ( f(c in) ) ) < ABS ( f(b out) ).
!
!    4, no change in sign of f(x) was found although the
!      interval (b,c) collapsed to the requested tolerance.
!      the user must examine this case and decide whether
!      b is near a local minimum of f(x), or B is near a
!      zero of even multiplicity, or neither of these.
!
!    5, too many (more than 500) function evaluations used.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) acbs
  real ( kind = 8 ) acmb
  real ( kind = 8 ) ae
  real ( kind = 8 ) aw
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cmb
  real ( kind = 8 ) er
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) fx
  real ( kind = 8 ) fz
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) kount
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) re
  real ( kind = 8 ) rw
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) z

  er = 2.0D+00 * epsilon ( er )
!
!  Initialize.
!
  z = r
  if ( r <= min ( b, c ) .or. max ( b, c ) <= r ) then
    z = c
  end if

  rw = max ( re, er )
  aw = max ( ae, 0.0D+00 )
  ic = 0
  t = z
  fz = f(t)
  fc = fz
  t = b
  fb = f(t)
  kount = 2

  if ( sign ( 1.0D+00, fz ) /= sign ( 1.0D+00, fb ) ) then

    c = z

  else if ( z /= c ) then

    t = c
    fc = f(t)
    kount = 3

    if ( sign ( 1.0D+00, fz ) /= sign ( 1.0D+00, fc ) ) then
      b = z
      fb = fz
    end if

  end if

  a = c
  fa = fc
  acbs = abs ( b - c )
  fx = max ( abs ( fb ), abs ( fc ) )

  do
!
!  Interchange
!
    if ( abs ( fc ) < abs ( fb ) ) then
      a = b
      fa = fb
      b = c
      fb = fc
      c = a
      fc = fa
    end if

    cmb = 0.5D+00 * ( c - b )
    acmb = abs ( cmb )
    tol = rw * abs ( b ) + aw
!
!  Test stopping criterion and function count.
!
    if ( acmb <= tol ) then
      exit
    end if

    if ( fb == 0.0D+00 ) then
      iflag = 2
      return
    end if

    if ( 500 <= kount ) then
      iflag = 5
      return
    end if
!
!  Calculate new iterate implicitly as b+p/q
!  where we arrange 0 <= p.
!  The implicit form is used to prevent overflow.
!
    p = ( b - a ) * fb
    q = fa - fb

    if ( p < 0.0D+00 ) then
      p = -p
      q = -q
    end if
!
!  Update A and check for satisfactory reduction
!  in the size of the bracketing interval.
!  If not, perform bisection.
!
5   continue

    a = b
    fa = fb
    ic = ic + 1

    if ( ic < 4 ) then
      go to 6
    end if

    if ( acbs <= 8.0D+00 * acmb ) then
      b = b + cmb
      go to 9
    end if

    ic = 0
    acbs = acmb
!
!  Test for too small a change
!
6   continue

    if ( abs ( q ) * tol < p ) then
      go to 7
    end if
!
!  Increment by tolerance
!
    b = b + sign ( tol, cmb )
    go to 9
!
!  Root ought to be between b and (c+b)/2.
!
7   continue
!
!  Use the secant rule or bisection.
!
    if ( p < cmb * q ) then
      b = b + p / q
    else
      b = b + cmb
    end if
!
!  Have completed computation for new iterate B.
!
9   continue

    t = b
    fb = f(t)
    kount = kount + 1
!
!  Decide whether next step is interpolation or extrapolation.
!
    if ( sign ( 1.0D+00, fb ) == sign ( 1.0D+00, fc ) ) then
      c = a
      fc = fa
    end if

  end do
!
!  Finished.  Process results for proper setting of IFlAG.
!
  if ( sign ( 1.0D+00, fb ) == sign ( 1.0D+00, fc ) ) then
    iflag = 4
    return
  end if

  if ( fx < abs ( fb ) ) then
    iflag = 3
    return
  end if

  iflag = 1

  return
end
subroutine gamlim ( xmin, xmax )

!*****************************************************************************80
!
!! GAMLIM computes the minimum and maximum bounds for X in GAMMA(X).
!
!  Discussion:
!
!    GAMLIM calculates the minimum and maximum legal bounds for X in GAMMA(X).
!
!  Modified:
!
!    11 August 2001
!
!  Parameters:
!
!    Output, real ( kind = 8 ) XMIN, the minimum legal value of X in GAMMA(X).
!    Any smaller value might result in underflow.
!
!    Output, real ( kind = 8 ) XMAX, the maximum legal value of X in GAMMA(X).
!    Any larger value will cause overflow.
!
  implicit none

  real ( kind = 8 ) alnbig
  real ( kind = 8 ) alnsml
  logical converged
  integer ( kind = 4 ) i
  real ( kind = 8 ) xln
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xold

  alnsml = log ( tiny ( alnsml ) )
  xmin = -alnsml

  converged = .false.

  do i = 1, 10

    xold = xmin
    xln = log ( xmin )

    xmin = xmin - xmin * &
      ( ( xmin + 0.5D+00 ) * xln - xmin - 0.2258D+00 + alnsml ) &
      / ( xmin * xln + 0.5D+00 )

    if ( abs ( xmin - xold ) < 0.005D+00 ) then
      converged = .true.
      exit
    end if

  end do

  if ( .not. converged ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMLIM - Fatal error!'
    write ( *, '(a)' ) '  Unable to determine XMIN.'
    stop
  end if

  xmin = -xmin + 0.01D+00

  alnbig = log ( huge ( alnbig ) )
  xmax = alnbig

  converged = .false.

  do i = 1, 10

    xold = xmax
    xln = log ( xmax )
    xmax = xmax - xmax * &
      ( ( xmax - 0.5D+00 ) * xln - xmax + 0.9189D+00 - alnbig ) &
      / ( xmax * xln - 0.5D+00 )

    if ( abs ( xmax - xold ) < 0.005D+00 ) then
      converged = .true.
      exit
    end if

  end do

  if ( .not. converged ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMLIM - Fatal error!'
    write ( *, '(a)' ) '  Unable to determine XMAX.'
    stop
  end if

  xmax = xmax - 0.01D+00
  xmin = max ( xmin, -xmax + 1.0D+00 )

  return
end
function gamma ( x )

!*****************************************************************************80
!
!! GAMMA computes the gamma function.
!
!  Parameters:
!
!    Input, real X, the argument of the gamma function, which must not
!    be 0, -1, or any other negative integral value.
!
!    Output, real GAMMA, the value of the gamma function of X.
!
  implicit none

  real ( kind = 8 ) csevl
  real ( kind = 8 ) d9lgmc
  real ( kind = 8 ), save :: dxrel = 0.0D+00
  real ( kind = 8 ) gamma
  real ( kind = 8 ), parameter, dimension ( 23 ) :: gcs = (/ &
   0.008571195590989331D+00, &
   0.004415381324841007D+00, &
   0.05685043681599363D+00, &
  -0.004219835396418561D+00, &
   0.001326808181212460D+00, &
  -0.0001893024529798880D+00, &
   0.0000360692532744124D+00, &
  -0.0000060567619044608D+00, &
   0.0000010558295463022D+00, &
  -0.0000001811967365542D+00, &
   0.0000000311772496471D+00, &
  -0.0000000053542196390D+00, &
   0.0000000009193275519D+00, &
  -0.0000000001577941280D+00, &
   0.0000000000270798062D+00, &
  -0.0000000000046468186D+00, &
   0.0000000000007973350D+00, &
  -0.0000000000001368078D+00, &
   0.0000000000000234731D+00, &
  -0.0000000000000040274D+00, &
   0.0000000000000006910D+00, &
  -0.0000000000000001185D+00, &
   0.0000000000000000203D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inits
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: ngcs = 0
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sinpiy
  real ( kind = 8 ), parameter :: sq2pil = 0.91893853320467274D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xmax = 0.0D+00
  real ( kind = 8 ), save :: xmin = 0.0D+00
  real ( kind = 8 ) y
!
!  Initialize.  Find legal bounds for X, and determine the number of
!  terms in the series required to attain an accuracy ten times better
!  than machine precision.
!
  if ( ngcs == 0 ) then
    ngcs = inits ( gcs, 23, 0.1D+00 * epsilon ( gcs ) )
    call gamlim ( xmin, xmax )
    dxrel = sqrt ( epsilon ( dxrel ) )
  end if

  y = abs ( x )
  if ( 10.0D+00 < y ) then
    go to 50
  end if
!
!  Compute gamma(x) for abs ( x ) <= 10.0.  Reduce interval and
!  find gamma(1+y) for 0 <= y < 1.
!
  n = int ( x )

  if ( x < 0.0D+00 ) then
    n = n - 1
  end if

  y = x - real ( n, kind = 8 )
  n = n - 1
  gamma = 0.9375D+00 + csevl ( 2.0D+00 * y - 1.0D+00, gcs, ngcs )

  if ( n == 0 ) then

    return

  else if ( n < 0 ) then

    n = -n
    if ( x == 0.0D+00 ) then
      call xerror ( 'GAMMA - x is 0', 4, 2)
    end if

    if ( x < 0.0D+00 .and. x + real ( n - 2, kind = 8 ) == 0.0D+00 ) then
       call xerror (  'GAMMA - x is a negative integer ( kind = 4 )', 4, 2 )
    end if

    if ( x < -0.5D+00 .and. &
      abs ( ( x - aint ( x - 0.5D+00) ) / x ) < dxrel ) then
      call xerror ( &
        'GAMMA - answer < half precision because x too near negative integer', &
        1, 1 )
    end if

    do i = 1, n
      gamma = gamma / ( x + real ( i - 1, kind = 8 ) )
    end do

    return

  else

    do i = 1, n
      gamma = ( y + real ( i, kind = 8 ) ) * gamma
    end do

    return

  end if
!
!  Compute gamma(x) for 10 < abs ( x ).  Recall y = abs ( x ).
!
 50   continue

  if ( xmax < x ) then
    call xerror ( 'GAMMA - x so big gamma overflows', 3, 2)
  end if

  gamma = 0.0D+00

  if ( x < xmin ) then
    call xerror ( 'GAMMA - x so small gamma underflows', 2, 1)
    return
  end if

  gamma = exp ( ( y - 0.5D+00 ) * log ( y ) - y + sq2pil + d9lgmc ( y ) )

  if ( 0.0D+00 < x ) then
    return
  end if

  if ( abs ( ( x - aint ( x - 0.5D+00 ) ) / x ) < dxrel ) then
    call xerror ( &
      'GAMMA - answer lt half precision, x too near negative integer', &
      1, 1)
  end if

  sinpiy = sin ( pi * y )

  if ( sinpiy == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMA - Fatal error!'
    write ( *, '(a)' ) '  The argument X is a negative integer ( kind = 4 ).'
    stop
  end if

  gamma = -pi / ( y * sinpiy * gamma )

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
!      Gamma(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) exp(-T) dT
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
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
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
    -0.20D+00, &
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
subroutine gl15t ( f, a, b, xl, xr, r, ae, ra, rasc, fmin, fmax )

!*****************************************************************************80
!
!! GL15T estimates the integral of a function over a finite interval.
!
!  Discussion:
!
!    GL15T is a utility routine for Q1DAX, and is not called directly
!    by the user.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, external F, the name of the routine that evaluates the function,
!    of the form
!
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!
!    The function G(X) is defined to be
!
!      G(X) = F(PHI(X)) * PHIP(X)
!
!    where PHI(X) is the cubic given by the arithmetic statement function
!    below and PHIP(X) is its derivative.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of integration.
!
!    Input, real ( kind = 8 ) XL, XR, the lower and upper limits of the
!    parent interval of which [A,B] is a part.
!
!    Output, real ( kind = 8 ) R, approximation to the integral I.  R is
!    computed by applying the 15-point Kronrod rule RESK obtained by optimal
!    addition of abscissae to the 7-point Gauss rule RESG.
!
!    Output, real ( kind = 8 ) AE, estimate of the modulus of the absolute
!    error, which should not exceed abs ( I - R ).
!
!    Output, real ( kind = 8 ) RA, approximation to the integral J.
!
!    Output, real ( kind = 8 ) RASC, approximation to the integral of
!    ABS ( G - I / (B-A) ) over (A,B).
!
!    Output, real ( kind = 8 ) FMAX, FMIN, the maximum and minimum values of the
!    function F(X) encountered during the computation.
!
!  Local variables:
!
!    centr  - mid point of the interval
!    hlgth  - half-length of the interval
!    absc   - abscissa
!    fval*  - function value
!    resg   - r of the 7-point Gauss formula
!    resk   - r of the 15-point Kronrod formula
!    reskh  - approximation to the mean value of f over (a,b), i.e. to i/(b-a)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) ae
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  real ( kind = 8 ), save :: epmach = 0.0D+00
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fc
  real ( kind = 8 ) fmax
  real ( kind = 8 ) fmin
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fv1(7)
  real ( kind = 8 ) fv2(7)
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) hlgth
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtw
  integer ( kind = 4 ) jtwm1
  real ( kind = 8 ) phi
  real ( kind = 8 ) phip
  real ( kind = 8 ) phiu
  real ( kind = 8 ) r
  real ( kind = 8 ) ra
  real ( kind = 8 ) rasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) sl
  real ( kind = 8 ) sr
  real ( kind = 8 ) u
  real ( kind = 8 ), save :: uflow = 0.0D+00
  real ( kind = 8 ), parameter, dimension ( 4 ) :: wg = (/ &
    0.129484966168869693270611432679082D+00, &
    0.279705391489276667901467771423780D+00, &
    0.381830050505118944950369775488975D+00, &
    0.417959183673469387755102040816327D+00 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: wgk = (/ &
    0.022935322010529224963732008058970D+00, &
    0.063092092629978553290700663189204D+00, &
    0.104790010322250183839876322541518D+00, &
    0.140653259715525918745189590510238D+00, &
    0.169004726639267902826583426598550D+00, &
    0.190350578064785409913256402421014D+00, &
    0.204432940075298892414161999234649D+00, &
    0.209482141084727828012999174891714D+00 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: xgk = (/ &
    0.991455371120812639206854697526329D+00, &
    0.949107912342758524526189684047851D+00, &
    0.864864423359769072789712788640926D+00, &
    0.741531185599394439863864773280788D+00, &
    0.586087235467691130294144838258730D+00, &
    0.405845151377397166906606412076961D+00, &
    0.207784955007898467600689403773245D+00, &
    0.000000000000000000000000000000000D+00 /)
  double precision xl
  double precision xr
!
!  Statement functions (ugh)
!
  phi ( u ) = xr - ( xr - xl ) * u * u * ( 2.0D+00 * u + 3.0D+00 )
  phip ( u ) = -6.0D+00 * u * ( u + 1.0D+00 )

  if ( epmach == 0.0D+00 ) then
    epmach = epsilon ( 1.0D+00 )
    uflow = tiny ( uflow )
  end if

  if ( xl < xr ) then
    sl = xl
    sr = xr
  else
    sl = xr
    sr = xl
  end if

  hlgth = 0.5D+00 * ( b - a )
  centr = a + hlgth
  dhlgth = abs ( hlgth )
!
!  Compute the 15-point Kronrod approximation to
!  the integral, and estimate the absolute error.
!
  u = ( centr - xr ) / ( xr - xl )
  phiu = phi(u)
  if ( phiu <= sl .or. sr <= phiu ) then
    phiu = centr
  end if

  fmin = f(phiu)
  fmax = fmin
  fc = fmin * phip(u)
  resg = fc * wg(4)
  resk = fc * wgk(8)
  ra = abs ( resk )

  do j = 1, 3

    jtw = j * 2
    absc = hlgth * xgk(jtw)
    u = ( centr - absc - xr ) / ( xr - xl )
    phiu = phi(u)
    if ( phiu <= sl .or. sr <= phiu ) then
      phiu = centr
    end if
    fval1 = f(phiu)
    fmax = max ( fmax, fval1 )
    fmin = min ( fmin, fval1 )
    fval1 = fval1 * phip(u)
    u = ( centr + absc - xr ) / ( xr - xl )
    phiu = phi(u)
    if ( phiu <= sl .or. sr <= phiu ) then
      phiu = centr
    end if
    fval2 = f(phiu)
    fmax = max ( fmax, fval2 )
    fmin = min ( fmin, fval2 )
    fval2 = fval2 * phip(u)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1 + fval2
    resg = resg + wg(j) * fsum
    resk = resk + wgk(jtw) * fsum
    ra = ra + wgk(jtw) * ( abs ( fval1 ) + abs ( fval2 ) )

  end do

  do j = 1, 4

    jtwm1 = j * 2 - 1
    absc = hlgth * xgk(jtwm1)
    u = ( centr - absc - xr ) / ( xr - xl )
    phiu = phi(u)
    if ( phiu <= sl .or. sr <= phiu ) then
      phiu = centr
    end if
    fval1 = f(phiu)
    fmax = max ( fmax, fval1 )
    fmin = min ( fmin, fval1 )
    fval1 = fval1 * phip(u)
    u = ( centr + absc - xr ) / ( xr - xl )

    phiu = phi(u)
    if ( phiu <= sl .or. sr <= phiu ) then
      phiu = centr
    end if

    fval2 = f(phiu)
    fmax = max ( fmax, fval2 )
    fmin = min ( fmin, fval2 )
    fval2 = fval2 * phip(u)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1 + fval2
    resk = resk + wgk(jtwm1) * fsum
    ra = ra + wgk(jtwm1) * ( abs ( fval1 ) + abs ( fval2 ) )

  end do

  reskh = resk * 0.5D+00
  rasc = wgk(8) * abs ( fc - reskh )

  do j = 1, 7
    rasc = rasc + wgk(j) * ( abs ( fv1(j) - reskh ) + abs ( fv2(j) - reskh ) )
  end do

  r = resk * hlgth
  ra = ra * dhlgth
  rasc = rasc * dhlgth
  ae = abs ( ( resk - resg ) * hlgth )

  if ( rasc /= 0.0D+00 .and. ae /= 0.0D+00 ) then
    ae = rasc * min ( 1.0D+00, ( 0.2D+03 * ae / rasc )**1.5D+00 )
  end if

  if ( uflow / ( 0.5D+02 * epmach ) < ra ) then
    ae = max ( ( epmach * 0.5D+02 ) * ra, ae )
  end if

  return
end
subroutine grdchk ( n, x, fcn, f, g, typsiz, sx, fscale, rnf, analtl, &
  wrk1, msg, ipr )

!*****************************************************************************80
!
!! GRDCHK checks an analytic gradient against an estimated gradient.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), a point at which the gradient is
!    to be checked.
!
!    Input, external FCN, the name of the subroutine that evaluates
!    the optimization function, of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real f
!      real x(n)
!
!    Input, real ( kind = 8 ) F, the function value at X.
!
!    Input, real ( kind = 8 ) G(N), the gradient value at X.
!
!    Input, real ( kind = 8 ) TYPSIZ(N), a typical size for each component of X.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling values:
!    SX(1:N)=1.0/TYPSIZ(1:N)
!
!    Input, real ( kind = 8 ) FSCALE, an estimate of the scale of the
!    objective function FCN.
!
!    Input, real ( kind = 8 ) RNF, the relative noise in the optimization
!    function FCN.
!
!    Input, real ( kind = 8 ) ANALTL, a tolerance for comparison of estimated
!    and analytical gradients
!
!    Workspace, real ( kind = 8 ) WRK1(N).
!
!    Output, integer ( kind = 4 ) MSG, message or error code.
!      0: no error detected.
!    -21: probable coding error of gradient
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) analtl
  real ( kind = 8 ) f
  external fcn
  real ( kind = 8 ) fscale
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) gs
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) ker
  integer ( kind = 4 ) msg
  real ( kind = 8 ) rnf
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) typsiz(n)
  real ( kind = 8 ) value(1)
  real ( kind = 8 ) wrk(1)
  real ( kind = 8 ) wrk1(n)
  real ( kind = 8 ) x(n)

  msg = 0
!
!  Compute the first order finite difference gradient;
!  compare it to the analytic gradient.
!
  value(1) = f
  call fstofd ( 1, 1, n, x, fcn, value, wrk1, sx, rnf, wrk, 1 )

  ker = 0

  do i = 1, n

    gs = max ( abs ( f ), fscale ) / max ( abs ( x(i) ), typsiz(i) )

    if ( max ( abs ( g(i) ), gs ) * analtl < abs ( g(i) - wrk1(i) ) ) then
      ker = 1
    end if

  end do

  if ( ker /= 0 ) then
    write ( ipr, * ) ' '
    write ( ipr, * ) 'GRDCHK - probable error in analytic gradient.'
    write ( ipr, * ) ' '
    write ( ipr, * ) ' grdchk     comp            analytic            est'
    write ( ipr, 902 ) ( i, g(i), wrk1(i), i = 1, n )
    msg = -21
  end if

  return

  902 format(' grdchk    ',i5,3x,e20.13,3x,e20.13)
end
subroutine heschk ( nr, n, x, fcn, d1fcn, d2fcn, f, g, a, typsiz, sx, rnf, &
  analtl, iagflg, udiag, wrk1, wrk2, msg, ipr )

!*****************************************************************************80
!
!! HESCHK checks an analytic hessian against a computed estimate.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), a point at which the Hessian is to
!    be checked.
!
!    Input, external FCN, the name of the subroutine that evaluates
!    the optimization function, of the form:
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!    Input, external D1FCN, the name of the subroutine to evaluate the
!    gradient of the function, of the form:
!
!      subroutine d1fcn ( n, x, g )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) g(n)
!      real ( kind = 8 ) x(n)
!
!    Input, external D2FCN, the name of the subroutine to evaluate the
!    Hessian of the function, of the form:
!
!      subroutine d2fcn ( nr, n, x, h )
!      integer ( kind = 4 ) nr
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) h(nr,n)
!      real ( kind = 8 ) x(n)
!
!    Input, real ( kind = 8 ) F, the function value at X.
!
!    Output, real ( kind = 8 ) G(N), the value of the gradient at X.
!
!    Output, real ( kind = 8 ) A(NR,N), the analytic Hessian matrix will
!    be stored in the lower triangle and diagonal.
!
!    Input, real ( kind = 8 ) TYPSIZ(N), a typical size for each component of X.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix.
!
!    Input, real ( kind = 8 ) RNF, the relative noise in the optimization
!    function FCN.
!
!    Input, real ( kind = 8 ) ANALTL, a tolerance for comparison of the
!    estimated and analytic gradients.
!
!    Input, integer ( kind = 4 ) IAGFLG, is 1 if the analytic gradient is supplied.
!
!    Workspace, real ( kind = 8 ) UDIAG(N).
!
!    Workspace, real ( kind = 8 ) WRK1(N).
!
!    Workspace, real ( kind = 8 ) WRK2(N).
!
!    Input/output, integer ( kind = 4 ) MSG, message or error code
!    on input : if =1xx do not compare analytic + estimated hessian.
!    on output: =-22, probable coding error of hessian.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) analtl
  external d1fcn
  external d2fcn
  real ( kind = 8 ) f
  external fcn
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) hs
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iagflg
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ker
  integer ( kind = 4 ) msg
  real ( kind = 8 ) rnf
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) typsiz(n)
  real ( kind = 8 ) udiag(n)
  real ( kind = 8 ) wrk1(n)
  real ( kind = 8 ) wrk2(n)
  real ( kind = 8 ) x(n)
!
!  Compute the finite difference approximation A to the hessian.
!
  if ( iagflg == 1 ) then
    call fstofd ( nr, n, n, x, d1fcn, g, a, sx, rnf, wrk1, 3 )
  else
    call sndofd ( nr, n, x, fcn, f, a, sx, rnf, wrk1, wrk2 )
  end if

  ker = 0
!
!  Copy lower triangular part of A to upper triangular part
!  and diagonal of A to UDIAG.
!
  do j = 1, n
    udiag(j) = a(j,j)
    do i = j+1, n
      a(j,i) = a(i,j)
    end do
  end do
!
!  Compute analytic hessian and compare to finite difference approximation.
!
  call d2fcn ( nr, n, x, a )

  do j = 1, n

    hs = max ( abs ( g(j) ), 1.0D+00 ) / max ( abs ( x(j) ), typsiz(j) )

    if ( max ( abs ( udiag(j) ), hs ) * analtl &
       < abs ( a(j,j) - udiag(j) ) ) then
      ker = 1
    end if

    do i = j+1, n
      if ( max ( abs ( a(i,j) ), hs ) * analtl &
        < abs ( a(i,j) - a(j,i) ) ) then
        ker = 1
      end if
    end do

  end do

  if ( ker /= 0 ) then

    write ( ipr, '(a)' ) ' '
    write ( ipr, '(a)' ) 'HESCHK:'
    write ( ipr, '(a)' ) '  Probable error in coding of analytic hessian.'
    write ( ipr, '(a)' ) '            row  col              analytic' &
      // '              (estimate)'
    write ( ipr, '(a)' ) ' '

    do i = 1, n
      do j = 1, i-1
        write(ipr,902) i, j, a(i,j), a(j,i)
      end do
      write(ipr,902) i, i, a(i,i), udiag(i)
    end do

    msg = -22

  end if

  return
  902 format('heschk    ',2i5,2x,e20.13,2x,'(',e20.13,')')
end
subroutine hookdr ( nr, n, x, f, g, a, udiag, p, xpls, fpls, fcn, sx, stepmx, &
  steptl, dlt, iretcd, mxtake, amu, dltp, phi, phip0, sc, xplsp, wrk0, epsm, &
  itncnt, ipr )

!*****************************************************************************80
!
!! HOOKDR finds the next Newton iterate by the More-Hebdon method.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, "X[K-1]".
!
!    Input, real ( kind = 8 ) F, the function value at the old iterate.
!
!    Input, real ( kind = 8 ) G(N), the gradient or an approximation, at
!    the old iterate.
!
!    Input, real ( kind = 8 ) A(NR,N), the Cholesky decomposition of hessian
!    in lower triangular part and diagonal.  Hessian in upper triangular
!    part and UDIAG.
!
!    Input, real ( kind = 8 ) UDIAG(N), the diagonal of the hessian matrix.
!
!    Input, real ( kind = 8 ) P(N), the Newton step.
!
!    Output, real ( kind = 8 ) XPLS(N), the new iterate X[K].
!
!    Output, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Input, external FCN, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) x(n)
!      real ( kind = 8 ) f
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which successive
!    iterates are considered close enough to terminate algorithm.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Output, integer ( kind = 4 ) IRETCD, return code
!    0, satisfactory xpls found
!    1, failed to find satisfactory xpls sufficiently distinct from x.
!
!    Output, logical MXTAKE, is TRUE if a step of maximum length was used.
!
!    Workspace, real ( kind = 8 ) AMU, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) DLTP, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) PHI, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) PHIP0, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) SC(N).
!
!    Workspace, real ( kind = 8 ) XPLSP(N).
!
!    Workspace, real ( kind = 8 ) WRK0(N).
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon
!
!    Input, integer ( kind = 4 ) ITNCNT, the iteration count.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) amu
  real ( kind = 8 ) beta
  real ( kind = 8 ) dlt
  real ( kind = 8 ) dltp
  real ( kind = 8 ) epsm
  real ( kind = 8 ) f
  external fcn
  real ( kind = 8 ) fpls
  real ( kind = 8 ) fplsp
  logical fstime
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) iretcd
  integer ( kind = 4 ) itncnt
  integer ( kind = 4 ) j
  logical mxtake
  logical nwtake
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) phi
  real ( kind = 8 ) phip0
  real ( kind = 8 ) rnwtln
  real ( kind = 8 ) sc(n)
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) steptl
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) tmp
  real ( kind = 8 ) udiag(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xpls(n)
  real ( kind = 8 ) xplsp(n)
  real ( kind = 8 ) wrk0(n)

  iretcd = 4
  fstime = .true.

  rnwtln = sqrt ( sum ( sx(1:n)**2 * p(1:n)**2 ) )
!
!  If first iteration and trust region not provided by user,
!  compute initial trust region.
!
  if ( itncnt <= 1 ) then

    amu = 0.0D+00

    if ( dlt == -1.0D+00 ) then

      alpha = sum ( ( g(1:n) / sx(1:n) )**2 )

      beta = 0.0D+00
      do i = 1, n
        tmp = 0.0D+00
        do j = i, n
          tmp = tmp + ( a(j,i) * g(j) ) / ( sx(j) * sx(j) )
        end do
        beta = beta + tmp * tmp
      end do

      dlt = alpha * sqrt ( alpha ) / beta
      dlt = min ( dlt, stepmx )

    end if

  end if
!
!  Find the new step by More-Hebdon algorithm.
!
  do

    call hookst ( nr, n, g, a, udiag, p, sx, rnwtln, dlt, amu, dltp, phi, &
      phip0, fstime, sc, nwtake, wrk0, epsm, ipr )

    dltp = dlt
!
!  Check the new point and update trust region.
!
    call tregup ( nr, n, x, f, g, a, fcn, sc, sx, nwtake, stepmx, steptl, &
      dlt, iretcd, xplsp, fplsp, xpls, fpls, mxtake, ipr, 3, udiag )

    if ( iretcd <= 1 ) then
      exit
    end if

  end do

  return
end
subroutine hookst ( nr, n, g, a, udiag, p, sx, rnwtln, dlt, amu, &
  dltp, phi, phip0, fstime, sc, nwtake, wrk0, epsm, ipr )

!*****************************************************************************80
!
!! HOOKST finds the new step by the More-Hebdon algorithm.
!
!  Modified:
!
!    15 May 2005
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) G(N), the gradient at the current iterate.
!
!    Input, real ( kind = 8 ) A(NR,N), an N by N array.  It contains the
!    Cholesky decomposition of the hessian in the lower triangular
!    part and diagonal; the hessian or approximation in the upper
!    triangular part.
!
!    Input, real ( kind = 8 ) UDIAG(N), the diagonal of the hessian. whose lower
!    triangular part is stored in A.
!
!    Input, real ( kind = 8 ) P(N), the Newton step.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) RNWTLN, the Newton step length.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Workspace, real ( kind = 8 ) AMU, [retain value between successive calls]
!
!    Input, real ( kind = 8 ) DLTP, the trust region radius at last exit
!    from this routine.
!
!    Workspace, real ( kind = 8 ) PHI, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) PHIP0, [retain value between successive calls]
!
!    Input/output, logical FSTIME, TRUE if first entry to this routine
!    during k-th iteration.
!
!    Output, real ( kind = 8 ) SC(N), the current step.
!
!    Output, logical NWTAKE, is TRUE if a Newton step taken.
!
!    Workspace, real ( kind = 8 ) WRK0(N).
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) addmax
  real ( kind = 8 ), parameter :: alo = 0.75D+00
  real ( kind = 8 ) amu
  real ( kind = 8 ) amulo
  real ( kind = 8 ) amuup
  real ( kind = 8 ) dlt
  real ( kind = 8 ) dltp
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) epsm
  logical fstime
  real ( kind = 8 ) g(n)
  real ( kind = 8 ), parameter :: hi = 1.50D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) j
  logical nwtake
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) phi
  real ( kind = 8 ) phip
  real ( kind = 8 ) phip0
  real ( kind = 8 ) rnwtln
  real ( kind = 8 ) sc(n)
  real ( kind = 8 ) stepln
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) udiag(n)
  real ( kind = 8 ) wrk0(n)
!
!  Take a Newton step?
!
  if ( rnwtln <= hi * dlt ) then
    nwtake = .true.
    sc(1:n) = p(1:n)
    dlt = min ( dlt, rnwtln )
    amu = 0.0D+00
    return
  end if
!
!  Newton step not taken.
!
  nwtake = .false.

  if ( 0.0D+00 < amu ) then
    amu = amu - ( phi + dltp ) * ( ( dltp - dlt ) + phi ) / ( dlt * phip )
  end if

  phi = rnwtln - dlt

  if ( fstime ) then

    wrk0(1:n) = sx(1:n) * sx(1:n) * p(1:n)
!
!  Solve L * Y = (SX**2)*P
!
    call forslv ( nr, n, a, wrk0, wrk0 )

    phip0 = -dnrm2 ( n, wrk0, 1 )**2 / rnwtln
    fstime = .false.

  end if

  phip = phip0
  amulo = -phi / phip
  amuup = 0.0D+00
  do i = 1, n
    amuup = amuup + ( g(i) * g(i) ) / ( sx(i) * sx(i) )
  end do
  amuup = sqrt ( amuup ) / dlt
!
!  Test the value of amu; generate next amu if necessary.
!
  do

    if ( amu < amulo .or. amuup < amu ) then
      amu = max ( sqrt ( amulo * amuup ), amuup * 1.0D-03 )
    end if
!
!  Copy (h,udiag) to L
!  where h <-- h + amu*(sx**2) [do not actually change (h,udiag)]
!
    do j = 1, n
      a(j,j) = udiag(j) + amu * sx(j) * sx(j)
      a(j+1:n,j) = a(j,j+1:n)
    end do
!
!  Factor h=l(l+)
!
    call choldc ( nr, n, a, 0.0D+00, sqrt ( epsm ), addmax )
!
!  Solve h*p = l(l+) * sc = -g
!
    wrk0(1:n) = -g(1:n)

    call lltslv ( nr, n, a, sc, wrk0 )
!
!  Reset H.  Note since UDIAG has not been destroyed, we need do
!  nothing here.  H is in the upper part and in UDIAG, still intact
!
    stepln = sqrt ( dot_product ( sx(1:n)**2, sc(1:n)**2 ) )

    phi = stepln - dlt

    wrk0(1:n) = sx(1:n)**2 * sc(1:n)

    call forslv ( nr, n, a, wrk0, wrk0 )

    phip = -dnrm2 ( n, wrk0, 1 )**2 / stepln
!
!  If SC not acceptable hookstep, then select new AMU.
!
    if ( ( stepln < alo * dlt .or. hi * dlt < stepln ) .and. &
      ( 0.0D+00 < amuup - amulo ) ) then

      amulo = max ( amulo, amu - ( phi / phip ) )

      if ( phi < 0.0D+00 ) then
        amuup = min ( amuup, amu )
      end if

      amu = amu - ( stepln * phi ) / ( dlt * phip )
!
!  SC is acceptable hookstep.
!
    else

      exit

    end if

  end do

  return
end
subroutine hsnint ( nr, n, a, sx, method )

!*****************************************************************************80
!
!! HSNINT provides initial hessian when using secant updates.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Output, real ( kind = 8 ) A(NR,N), the initial N by N Hessian.  Only the
!    lower triangle of the matrix is assigned values.
!
!    Input, real ( kind = 8 ) SX(N), the scaling factors for X.
!
!    Input, integer ( kind = 4 ) METHOD, specifies the algorithm to use to solve
!    the minimization problem.
!    1 or 2: factored secant method used
!    3:  unfactored secant method used
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) method
  real ( kind = 8 ) sx(n)

  do j = 1, n

    if ( method == 3 ) then
      a(j,j) = sx(j)**2
    else
      a(j,j) = sx(j)
    end if

    a(j+1:n,j) = 0.0D+00

  end do

  return
end
function i1mach ( i )

!*****************************************************************************80
!
!! I1MACH returns integer ( kind = 4 ) machine constants.
!
!  I/O unit numbers.
!
!    I1MACH(1) = the standard input unit.
!    I1MACH(2) = the standard output unit.
!    I1MACH(3) = the standard punch unit.
!    I1MACH(4) = the standard error message unit.
!
!  Words.
!
!    I1MACH(5) = the number of bits per integer ( kind = 4 ) storage unit.
!    I1MACH(6) = the number of characters per integer ( kind = 4 ) storage unit.
!
!  Integers.
!
!  Assume integer ( kind = 4 )s are represented in the S digit base A form:
!
!  Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
!  where 0<=X(I)<A for I=0 to S-1.
!
!    I1MACH(7) = A, the base.
!    I1MACH(8) = S, the number of base A digits.
!    I1MACH(9) = A**S-1, the largest integer ( kind = 4 ).
!
!  Floating point numbers
!
!  Assume floating point numbers are represented in the T digit base B form:
!
!    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
!
!  where 0<=X(I)<B for I = 1 to T, 0<X(1) and EMIN<=E<=EMAX
!
!    I1MACH(10) = B, the base.
!
!  Single precision
!
!    I1MACH(11) = T, the number of base B digits.
!    I1MACH(12) = EMIN, the smallest exponent E.
!    I1MACH(13) = EMAX, the largest exponent E.
!
!  Double precision
!
!    I1MACH(14) = T, the number of base B digits.
!    I1MACH(15) = EMIN, the smallest exponent E.
!    I1MACH(16) = EMAX, the largest exponent E.
!
!  To alter this function for a particular environment, the desired set of DATA
!  statements should be activated by removing the C from column 1.  On rare
!  machines, a STATIC statement may need to be added, but probably more systems
!  prohibit than require it.
!
!  Also, the values of I1MACH(1) through I1MACH(4) should be checked for
!  consistency with the local operating system.  For FORTRAN 77, you may wish
!  to adjust the data statement so imach(6) is set to 1, and then to comment
!  out the executable test on I.EQ.6 below.
!
!  For IEEE-arithmetic machines (binary standard), the first set of constants
!  below should be appropriate, except perhaps for IMACH(1) - IMACH(4).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) imach(16)
  integer ( kind = 4 ) output

  equivalence (imach(4),output)
!
!  IEEE arithmetic machines, such as the ATT 3B series, Motorola
!  68000 based machines such as the SUN 3 and ATT PC 7300, and
!  8087 based micros such asthe IBM PC and ATT 6300.
!
   data imach( 1) /    5 /
   data imach( 2) /    6 /
   data imach( 3) /    7 /
   data imach( 4) /    6 /
   data imach( 5) /   32 /
   data imach( 6) /    4 /
   data imach( 7) /    2 /
   data imach( 8) /   31 /
   data imach( 9) / 2147483647 /
   data imach(10) /    2 /
   data imach(11) /   24 /
   data imach(12) / -125 /
   data imach(13) /  128 /
   data imach(14) /   53 /
   data imach(15) / -1021 /
   data imach(16) /  1024 /
!
!  ALLIANT FX/8 UNIX FORTRAN compiler.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) /2147483647/
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1024 /
!
!  AMDAHL machines.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  63 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  63 /
!
!  BURROUGHS 1700 system.
!
!      data imach( 1) /    7 /
!      data imach( 2) /    2 /
!      data imach( 3) /    2 /
!      data imach( 4) /    2 /
!      data imach( 5) /   36 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   33 /
!      data imach( 9) / Z1FFFFFFFF /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -256 /
!      data imach(13) /  255 /
!      data imach(14) /   60 /
!      data imach(15) / -256 /
!      data imach(16) /  255 /
!
!  BURROUGHS 5700 system.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  48 /
!      data imach( 6) /   6 /
!      data imach( 7) /   2 /
!      data imach( 8) /  39 /
!      data imach( 9) / O0007777777777777 /
!      data imach(10) /   8 /
!      data imach(11) /  13 /
!      data imach(12) / -50 /
!      data imach(13) /  76 /
!      data imach(14) /  26 /
!      data imach(15) / -50 /
!      data imach(16) /  76 /
!
!  BURROUGHS 6700/7700 systems.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  48 /
!      data imach( 6) /   6 /
!      data imach( 7) /   2 /
!      data imach( 8) /  39 /
!      data imach( 9) / O0007777777777777 /
!      data imach(10) /   8 /
!      data imach(11) /  13 /
!      data imach(12) / -50 /
!      data imach(13) /  76 /
!      data imach(14) /  26 /
!      data imach(15) / -32754 /
!      data imach(16) /  32780 /
!
!  CDC CYBER 170/180 series using NOS
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / O"00007777777777777777" /
!      data imach(10) /    2 /
!      data imach(11) /   48 /
!      data imach(12) / -974 /
!      data imach(13) / 1070 /
!      data imach(14) /   96 /
!      data imach(15) / -927 /
!      data imach(16) / 1070 /
!
!  CDC CYBER 170/180 series using NOS/VE
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     7 /
!      data imach( 4) /     6 /
!      data imach( 5) /    64 /
!      data imach( 6) /     8 /
!      data imach( 7) /     2 /
!      data imach( 8) /    63 /
!      data imach( 9) / 9223372036854775807 /
!      data imach(10) /     2 /
!      data imach(11) /    47 /
!      data imach(12) / -4095 /
!      data imach(13) /  4094 /
!      data imach(14) /    94 /
!      data imach(15) / -4095 /
!      data imach(16) /  4094 /
!
!  CDC CYBER 200 series
!
!      data imach( 1) /      5 /
!      data imach( 2) /      6 /
!      data imach( 3) /      7 /
!      data imach( 4) /      6 /
!      data imach( 5) /     64 /
!      data imach( 6) /      8 /
!      data imach( 7) /      2 /
!      data imach( 8) /     47 /
!      data imach( 9) / X'00007FFFFFFFFFFF' /
!      data imach(10) /      2 /
!      data imach(11) /     47 /
!      data imach(12) / -28625 /
!      data imach(13) /  28718 /
!      data imach(14) /     94 /
!      data imach(15) / -28625 /
!      data imach(16) /  28718 /
!
!  CDC 6000/7000 series using FTN4.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / 00007777777777777777B /
!      data imach(10) /    2 /
!      data imach(11) /   47 /
!      data imach(12) / -929 /
!      data imach(13) / 1070 /
!      data imach(14) /   94 /
!      data imach(15) / -929 /
!      data imach(16) / 1069 /
!
!  CDC 6000/7000 series using FTN5.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / O"00007777777777777777" /
!      data imach(10) /    2 /
!      data imach(11) /   47 /
!      data imach(12) / -929 /
!      data imach(13) / 1070 /
!      data imach(14) /   94 /
!      data imach(15) / -929 /
!      data imach(16) / 1069 /
!
!  CONVEX C-1.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   53 /
!      data imach(15) /-1024 /
!      data imach(16) / 1023 /
!
!  CONVEX C-120 (native mode) without -R8 option
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    0 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   53 /
!      data imach(15) / -1023 /
!      data imach(16) /  1023 /
!
!  CONVEX C-120 (native mode) with -R8 option
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     0 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    53 /
!      data imach(12) / -1023 /
!      data imach(13) /  1023 /
!      data imach(14) /    53 /
!      data imach(15) / -1023 /
!      data imach(16) /  1023 /
!
!  CONVEX C-120 (IEEE mode) without -R8 option
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    0 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -125 /
!      data imach(13) /  128 /
!      data imach(14) /   53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  CONVEX C-120 (IEEE mode) with -R8 option
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     0 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    53 /
!      data imach(12) / -1021 /
!      data imach(13) /  1024 /
!      data imach(14) /    53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  CRAY 1, 2, XMP and YMP.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /   102 /
!      data imach( 4) /     6 /
!      data imach( 5) /    64 /
!      data imach( 6) /     8 /
!      data imach( 7) /     2 /
!      data imach( 8) /    63 /
!      data imach( 9) /  777777777777777777777B /
!      data imach(10) /     2 /
!      data imach(11) /    47 /
!      data imach(12) / -8189 /
!      data imach(13) /  8190 /
!      data imach(14) /    94 /
!      data imach(15) / -8099 /
!      data imach(16) /  8190 /
!
!  DATA GENERAL ECLIPSE S/200.
!
!      data imach( 1) /   11 /
!      data imach( 2) /   12 /
!      data imach( 3) /    8 /
!      data imach( 4) /   10 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) /32767 /
!      data imach(10) /   16 /
!      data imach(11) /    6 /
!      data imach(12) /  -64 /
!      data imach(13) /   63 /
!      data imach(14) /   14 /
!      data imach(15) /  -64 /
!      data imach(16) /   63 /
!
!  ELXSI 6400
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  HARRIS 220
!
!      data imach( 1) /       5 /
!      data imach( 2) /       6 /
!      data imach( 3) /       0 /
!      data imach( 4) /       6 /
!      data imach( 5) /      24 /
!      data imach( 6) /       3 /
!      data imach( 7) /       2 /
!      data imach( 8) /      23 /
!      data imach( 9) / 8388607 /
!      data imach(10) /       2 /
!      data imach(11) /      23 /
!      data imach(12) /    -127 /
!      data imach(13) /     127 /
!      data imach(14) /      38 /
!      data imach(15) /    -127 /
!      data imach(16) /     127 /
!
!  HARRIS SLASH 6 and SLASH 7.
!
!      data imach( 1) /       5 /
!      data imach( 2) /       6 /
!      data imach( 3) /       0 /
!      data imach( 4) /       6 /
!      data imach( 5) /      24 /
!      data imach( 6) /       3 /
!      data imach( 7) /       2 /
!      data imach( 8) /      23 /
!      data imach( 9) / 8388607 /
!      data imach(10) /       2 /
!      data imach(11) /      23 /
!      data imach(12) /    -127 /
!      data imach(13) /     127 /
!      data imach(14) /      38 /
!      data imach(15) /    -127 /
!      data imach(16) /     127 /
!
!  HONEYWELL DPS 8/70 and 600/6000 series.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /   43 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / O377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   63 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  HP 2100, 3 word double precision option with FTN4
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    4 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   23 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   39 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  HP 2100, 4 word double precision option with FTN4
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    4 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   23 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   55 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  HP 9000
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     7 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1015 /
!      data imach(16) /  1017 /
!
!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!  and PERKIN ELMER (INTERDATA) 3230.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / Z7FFFFFFF /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  63 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  63 /
!
!  IBM PC - Microsoft FORTRAN
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!
!      data imach( 1) /     4 /
!      data imach( 2) /     7 /
!      data imach( 3) /     7 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!  For the INTERDATA FORTRAN VII compiler, replace the Z's specifying hex
!  constants with Y's.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   6 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / Z'7FFFFFFF' /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  62 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  62 /
!
!  PDP-10 (KA processor).
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    5 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / "377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   54 /
!      data imach(15) / -101 /
!      data imach(16) /  127 /
!
!  PDP-10 (KI processor).
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    5 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / "377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   62 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  PDP-11 FORTRANS supporting 32-bit integer ( kind = 4 ) arithmetic.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  PDP-11 FORTRANS supporting 16-bit integer ( kind = 4 ) arithmetic.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  PRIME 50 series systems with 32-bit integer ( kind = 4 )s and 64V MODE instructions,
!  supplied by Igor Bray.
!
!      data imach( 1) /            1 /
!      data imach( 2) /            1 /
!      data imach( 3) /            2 /
!      data imach( 4) /            1 /
!      data imach( 5) /           32 /
!      data imach( 6) /            4 /
!      data imach( 7) /            2 /
!      data imach( 8) /           31 /
!      data imach( 9) / :17777777777 /
!      data imach(10) /            2 /
!      data imach(11) /           23 /
!      data imach(12) /         -127 /
!      data imach(13) /         +127 /
!      data imach(14) /           47 /
!      data imach(15) /       -32895 /
!      data imach(16) /       +32637 /
!
!  SEQUENT BALANCE 8000.
!
!      data imach( 1) /     0 /
!      data imach( 2) /     0 /
!      data imach( 3) /     7 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     1 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) /  2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -125 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  SUN Microsystems UNIX F77 compiler.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) /2147483647/
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1024 /
!
!  SUN 3 (68881 or FPA)
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    6 /
!      data imach( 4) /    0 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -125 /
!      data imach(13) /  128 /
!      data imach(14) /   53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  UNIVAC 1100 series.
!  Note that the punch unit, I1MACH(3), has been set to 7, which is appropriate
!  for the UNIVAC-FOR system.  If you have the UNIVAC-FTN system, set it to 1
!  instead.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    6 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / O377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   60 /
!      data imach(15) /-1024 /
!      data imach(16) / 1023 /
!
!  VAX.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  Z80 microprocessor.
!
!      data imach( 1) /    1 /
!      data imach( 2) /    1 /
!      data imach( 3) /    0 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
  if ( i < 1 .or. 16 < i )then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a,i8)' ) '  I is out of bounds:', i
    i1mach = 0
    stop
  else
    i1mach = imach(i)
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integer ( kind = 4 ) values.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i8_factor ( n, ifac )

!*****************************************************************************80
!
!! I8_FACTOR factors an integer ( kind = 4 ).
!
!  Modified:
!
!    14 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number to be factored.
!
!    Output, integer ( kind = 4 ) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifac(15)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry

  ifac(1) = n

  nf = 0
  nl = n

  if ( n == 0 ) then
    nf = 1
    ifac(2) = nf
    ifac(2+nf) = 0
    return
  end if

  if ( n < 1 ) then
    nf = nf + 1
    ifac(2+nf) = -1
    nl = - n
  end if

  if ( nl == 1 ) then
    nf = nf + 1
    ifac(2) = nf
    ifac(2+nf) = 1
    return
  end if

  j = 0

  do while ( 1 < nl )

    j = j + 1
!
!  Choose a trial divisor, NTRY.
!
    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if
!
!  Divide by the divisor as many times as possible.
!
    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nl = nq
      nf = nf + 1
!
!  Make sure factors of 2 appear in the front of the list.
!
      if ( ntry /= 2 ) then

        ifac(2+nf) = ntry

      else

        do i = nf, 2, -1
          ifac(i+2) = ifac(i+1)
        end do
        ifac(3) = 2

      end if

    end do

  end do

  ifac(2) = nf

  return
end
function idamax ( n, dx, incx )

!*****************************************************************************80
!
!! IDAMAX finds the index of the vector element of maximum absolute value.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries of SX.
!
!    Output, integer ( kind = 4 ) IDAMAX, the index of the element of SX of maximum
!    absolute value.
!
  implicit none


  real ( kind = 8 ) dmax
  real ( kind = 8 ) dx(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n

  idamax = 0

  if ( n < 1 .or. incx <= 0 ) then
    return
  end if

  idamax = 1

  if ( n == 1 ) then
    return
  end if

  if ( incx == 1 ) then

    dmax = abs ( dx(1) )

    do i = 2, n
      if ( dmax < abs ( dx(i) ) ) then
        idamax = i
        dmax = abs ( dx(i) )
      end if
    end do

  else

    ix = 1
    dmax = abs ( dx(1) )
    ix = ix + incx

    do i = 2, n
      if ( dmax < abs ( dx(ix) ) ) then
        idamax = i
        dmax = abs ( dx(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function inbin ( x, nbins, xmin, xmax, width )

!*****************************************************************************80
!
!! INBIN takes a real value X and finds the correct bin for it.
!
!  Discussion:
!
!    Values below XMIN come back in 1.  Values above XMAX come back
!    in NBINS.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, a value to be binned.
!
!    Input, integer ( kind = 4 ) NBINS, the number of bins.
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the minimum and maximum bin limits.
!
!    Input, real ( kind = 8 ) WIDTH, the width of each bin.
!
!    Output, integer ( kind = 4 ) INBIN, the index of the bin containing X.
!
  implicit none

  integer ( kind = 4 ) inbin
  integer ( kind = 4 ) nbins
  real ( kind = 8 ) width
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  if ( x < xmin ) then
    inbin = 1
  else if ( xmax <= x ) then
    inbin = nbins
  else
    inbin = 2 + int ( ( x - xmin ) / width )
  end if

  return
end
function inits ( os, nos, eta )

!*****************************************************************************80
!
!! INITS estimates the order of an orthogonal series for a given accuracy.
!
!  Discussion:
!
!    Because this is a function, it is not possible to print out
!    warning error messages.  Therefore, if an error condition is
!    detected, a bogus value of INITS is returned.
!
!  Modified:
!
!    15 April 2003
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) OS(NOS), the coefficients in the series.
!
!    Input, integer ( kind = 4 ) NOS, the number of coefficients.  NOS must be
!    at least 1, or an error condition arises.
!
!    Input, real ( kind = 8 ) ETA, the requested accuracy of the series.
!    Ordinarily, ETA will be chosen to be one-tenth machine precision.
!
!    Output, integer ( kind = 4 ) INITS, the order of the series guaranteeing the
!    given accuracy.  However, on error, INITS will be returned
!    as a negative huge number.
!
  implicit none

  integer ( kind = 4 ) nos

  real ( kind = 8 ) err
  real ( kind = 8 ) eta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inits
  real ( kind = 8 ) os(nos)
!
!  Fatal error.  Number of coefficients less than 1.
!
  if ( nos < 1 ) then
    inits = - huge ( i )
    return
  end if

  err = 0.0D+00

  i = 0

  do ii = 1, nos

    i = nos + 1 - ii
    err = err + abs ( os(i) )

    if ( eta < err ) then
      i = nos + 1 - ii
      exit
    end if

  end do

  if ( i == 0 ) then
    i = nos
    call xerror ( 'INITS - ETA may be too small', 1, 2)
  end if

  inits = i

  return
end
function j4save ( iwhich, ivalue, iset )

!*****************************************************************************80
!
!! J4SAVE saves variables needed by the library error handling routines.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IWHICH, the index of the item desired.
!    1, the current error number.
!    2, the current error control flag.
!    3, the current unit number to which error messages are sent.
!       (0 means use standard.)
!    4, the maximum times any message is printed (as set by xermax).
!    5, the number of units to which each error message is written.
!    6, the 2nd unit for error messages.
!    7, the 3rd unit for error messages.
!    8, the 4th unit for error messages.
!    9, the 5th unit for error messages.
!
!    Input, integer ( kind = 4 ) IVALUE, the value to be set for the IWHICH-th 
!    parameter, if ISET is TRUE.
!
!    Input, logical ISET.
!    TRUE: the IWHICH-th parameter will be given the value, IVALUE.
!
!    Output, integer ( kind = 4 ) J4SAVE, the old value of the IWHICH-th 
!    parameter.
!
  implicit none

  integer ( kind = 4 ), save, dimension ( 9 ) :: iparam = (/ &
    0, 2, 0, 10, 1, 0, 0, 0, 0 /)
  logical iset
  integer ( kind = 4 ) ivalue
  integer ( kind = 4 ) iwhich
  integer ( kind = 4 ) j4save

  j4save = iparam(iwhich)

  if ( iset ) then
    iparam(iwhich) = ivalue
  end if

  return
end
subroutine jairy ( x, rx, c, ai, dai )

!*****************************************************************************80
!
!! JAIRY computes the Airy function and its derivative.
!
!  Author:
!
!    D E Amos,
!    S L Daniel,
!    M K Weston.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, real ( kind = 8 ) RX, sqrt ( abs ( X ) ).
!
!    Input, real ( kind = 8 ) C, = 2 * ( ABS ( X )**1.5 ) / 3, computed
!    by ASYJY.
!
!    Output, real ( kind = 8 ) AI, the value of the Airy function.
!
!    Output, real ( kind = 8 ) DAI, the derivative of the Airy function.
!
  implicit none

  real ( kind = 8 ) a(15)
  real ( kind = 8 ) ai
  real ( kind = 8 ) ajn(19)
  real ( kind = 8 ) ajp(19)
  real ( kind = 8 ), parameter, dimension ( 14 ) :: ak1 = (/ &
    2.20423090987793D-01,-1.25290242787700D-01, 1.03881163359194D-02, &
    8.22844152006343D-04,-2.34614345891226D-04, 1.63824280172116D-05, &
    3.06902589573189D-07,-1.29621999359332D-07, 8.22908158823668D-09, &
    1.53963968623298D-11, -3.39165465615682D-11, 2.03253257423626D-12, &
   -1.10679546097884D-14, -5.16169497785080D-15 /)
  real ( kind = 8 ), parameter, dimension ( 23 ) :: ak2 = (/ &
      2.74366150869598D-01, 5.39790969736903D-03, &
     -1.57339220621190D-03, 4.27427528248750D-04,-1.12124917399925D-04, &
      2.88763171318904D-05,-7.36804225370554D-06, 1.87290209741024D-06, &
     -4.75892793962291D-07, 1.21130416955909D-07,-3.09245374270614D-08, &
      7.92454705282654D-09,-2.03902447167914D-09, 5.26863056595742D-10, &
     -1.36704767639569D-10, 3.56141039013708D-11,-9.31388296548430D-12, &
      2.44464450473635D-12,-6.43840261990955D-13, 1.70106030559349D-13, &
     -4.50760104503281D-14, 1.19774799164811D-14,-3.19077040865066D-15 /)
  real ( kind = 8 ) ak3(14)
  real ( kind = 8 ) b(15)
  real ( kind = 8 ) c
  real ( kind = 8 ) ccv
  real ( kind = 8 ), parameter :: con2 = 5.03154716196777D+00
  real ( kind = 8 ), parameter :: con3 = 3.80004589867293D-01
  real ( kind = 8 ), parameter :: con4 = 8.33333333333333D-01
  real ( kind = 8 ), parameter :: con5 = 8.66025403784439D-01
  real ( kind = 8 ) cv
  real ( kind = 8 ) da(15)
  real ( kind = 8 ) dai
  real ( kind = 8 ) dajn(19)
  real ( kind = 8 ) dajp(19)
  real ( kind = 8 ) dak1(14)
  real ( kind = 8 ) dak2(24)
  real ( kind = 8 ) dak3(14)
  real ( kind = 8 ) db(15)
  real ( kind = 8 ) ec
  real ( kind = 8 ) e1
  real ( kind = 8 ) e2
  real ( kind = 8 ), parameter :: fpi12 = 1.30899693899575D+00
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: m1 = 12
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ), parameter :: m2 = 21
  integer ( kind = 4 ) m2d
  integer ( kind = 4 ), parameter :: m3 = 17
  integer ( kind = 4 ) m3d
  integer ( kind = 4 ), parameter :: m4 = 13
  integer ( kind = 4 ) m4d
  integer ( kind = 4 ), parameter :: n1 = 14
  integer ( kind = 4 ) n1d
  integer ( kind = 4 ), parameter :: n2 = 23
  integer ( kind = 4 ) n2d
  integer ( kind = 4 ), parameter :: n3 = 19
  integer ( kind = 4 ) n3d
  integer ( kind = 4 ), parameter :: n4 = 15
  integer ( kind = 4 ) n4d
  real ( kind = 8 ) rtrx
  real ( kind = 8 ) rx
  real ( kind = 8 ) scv
  real ( kind = 8 ) t
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  real ( kind = 8 ) tt
  real ( kind = 8 ) x

  data ak3(1), ak3(2), ak3(3), ak3(4), ak3(5), ak3(6), ak3(7), &
          ak3(8), ak3(9), ak3(10),ak3(11),ak3(12),ak3(13), &
          ak3(14)         / 2.80271447340791D-01,-1.78127042844379D-03, &
      4.03422579628999D-05,-1.63249965269003D-06, 9.21181482476768D-08, &
     -6.52294330229155D-09, 5.47138404576546D-10,-5.24408251800260D-11, &
      5.60477904117209D-12,-6.56375244639313D-13, 8.31285761966247D-14, &
     -1.12705134691063D-14, 1.62267976598129D-15,-2.46480324312426D-16/

  data ajp(1), ajp(2), ajp(3), ajp(4), ajp(5), ajp(6), ajp(7), &
          ajp(8), ajp(9), ajp(10),ajp(11),ajp(12),ajp(13),ajp(14), &
          ajp(15),ajp(16),ajp(17),ajp(18), &
          ajp(19)         / 7.78952966437581D-02,-1.84356363456801D-01, &
      3.01412605216174D-02, 3.05342724277608D-02,-4.95424702513079D-03, &
     -1.72749552563952D-03, 2.43137637839190D-04, 5.04564777517082D-05, &
     -6.16316582695208D-06,-9.03986745510768D-07, 9.70243778355884D-08, &
      1.09639453305205D-08,-1.04716330588766D-09,-9.60359441344646D-11, &
      8.25358789454134D-12, 6.36123439018768D-13,-4.96629614116015D-14, &
     -3.29810288929615D-15, 2.35798252031104D-16/

  data ajn(1), ajn(2), ajn(3), ajn(4), ajn(5), ajn(6), ajn(7), &
          ajn(8), ajn(9), ajn(10),ajn(11),ajn(12),ajn(13),ajn(14), &
          ajn(15),ajn(16),ajn(17),ajn(18), &
          ajn(19)         / 3.80497887617242D-02,-2.45319541845546D-01, &
      1.65820623702696D-01, 7.49330045818789D-02,-2.63476288106641D-02, &
     -5.92535597304981D-03, 1.44744409589804D-03, 2.18311831322215D-04, &
     -4.10662077680304D-05,-4.66874994171766D-06, 7.15218807277160D-07, &
      6.52964770854633D-08,-8.44284027565946D-09,-6.44186158976978D-10, &
      7.20802286505285D-11, 4.72465431717846D-12,-4.66022632547045D-13, &
     -2.67762710389189D-14, 2.36161316570019D-15/

  data a(1),   a(2),   a(3),   a(4),   a(5),   a(6),   a(7), &
          a(8),   a(9),   a(10),  a(11),  a(12),  a(13),  a(14), &
          a(15)           / 4.90275424742791D-01, 1.57647277946204D-03, &
     -9.66195963140306D-05, 1.35916080268815D-07, 2.98157342654859D-07, &
     -1.86824767559979D-08,-1.03685737667141D-09, 3.28660818434328D-10, &
     -2.57091410632780D-11,-2.32357655300677D-12, 9.57523279048255D-13, &
     -1.20340828049719D-13,-2.90907716770715D-15, 4.55656454580149D-15, &
     -9.99003874810259D-16/

  data b(1),   b(2),   b(3),   b(4),   b(5),   b(6),   b(7), &
          b(8),   b(9),   b(10),  b(11),  b(12),  b(13),  b(14), &
          b(15)           / 2.78593552803079D-01,-3.52915691882584D-03, &
     -2.31149677384994D-05, 4.71317842263560D-06,-1.12415907931333D-07, &
     -2.00100301184339D-08, 2.60948075302193D-09,-3.55098136101216D-11, &
     -3.50849978423875D-11, 5.83007187954202D-12,-2.04644828753326D-13, &
     -1.10529179476742D-13, 2.87724778038775D-14,-2.88205111009939D-15, &
     -3.32656311696166D-16/

  data n1d,n2d,n3d,n4d/14,24,19,15/
  data m1d,m2d,m3d,m4d/12,22,17,13/

  data dak1(1), dak1(2), dak1(3), dak1(4), dak1(5), dak1(6), &
          dak1(7), dak1(8), dak1(9), dak1(10),dak1(11),dak1(12), &
         dak1(13),dak1(14)/ 2.04567842307887D-01,-6.61322739905664D-02, &
     -8.49845800989287D-03, 3.12183491556289D-03,-2.70016489829432D-04, &
     -6.35636298679387D-06, 3.02397712409509D-06,-2.18311195330088D-07, &
     -5.36194289332826D-10, 1.13098035622310D-09,-7.43023834629073D-11, &
      4.28804170826891D-13, 2.23810925754539D-13,-1.39140135641182D-14/

  data dak2(1), dak2(2), dak2(3), dak2(4), dak2(5), dak2(6), &
          dak2(7), dak2(8), dak2(9), dak2(10),dak2(11),dak2(12), &
          dak2(13),dak2(14),dak2(15),dak2(16),dak2(17),dak2(18), &
          dak2(19),dak2(20),dak2(21),dak2(22),dak2(23), &
          dak2(24)        / 2.93332343883230D-01,-8.06196784743112D-03, &
      2.42540172333140D-03,-6.82297548850235D-04, 1.85786427751181D-04, &
     -4.97457447684059D-05, 1.32090681239497D-05,-3.49528240444943D-06, &
      9.24362451078835D-07,-2.44732671521867D-07, 6.49307837648910D-08, &
     -1.72717621501538D-08, 4.60725763604656D-09,-1.23249055291550D-09, &
      3.30620409488102D-10,-8.89252099772401D-11, 2.39773319878298D-11, &
     -6.48013921153450D-12, 1.75510132023731D-12,-4.76303829833637D-13, &
      1.29498241100810D-13,-3.52679622210430D-14, 9.62005151585923D-15, &
     -2.62786914342292D-15/

  data dak3(1), dak3(2), dak3(3), dak3(4), dak3(5), dak3(6), &
          dak3(7), dak3(8), dak3(9), dak3(10),dak3(11),dak3(12), &
         dak3(13),dak3(14)/ 2.84675828811349D-01, 2.53073072619080D-03, &
     -4.83481130337976D-05, 1.84907283946343D-06,-1.01418491178576D-07, &
      7.05925634457153D-09,-5.85325291400382D-10, 5.56357688831339D-11, &
     -5.90889094779500D-12, 6.88574353784436D-13,-8.68588256452194D-14, &
      1.17374762617213D-14,-1.68523146510923D-15, 2.55374773097056D-16/

  data dajp(1), dajp(2), dajp(3), dajp(4), dajp(5), dajp(6), &
          dajp(7), dajp(8), dajp(9), dajp(10),dajp(11),dajp(12), &
          dajp(13),dajp(14),dajp(15),dajp(16),dajp(17),dajp(18), &
          dajp(19)        / 6.53219131311457D-02,-1.20262933688823D-01, &
      9.78010236263823D-03, 1.67948429230505D-02,-1.97146140182132D-03, &
     -8.45560295098867D-04, 9.42889620701976D-05, 2.25827860945475D-05, &
     -2.29067870915987D-06,-3.76343991136919D-07, 3.45663933559565D-08, &
      4.29611332003007D-09,-3.58673691214989D-10,-3.57245881361895D-11, &
      2.72696091066336D-12, 2.26120653095771D-13,-1.58763205238303D-14, &
     -1.12604374485125D-15, 7.31327529515367D-17/

  data dajn(1), dajn(2), dajn(3), dajn(4), dajn(5), dajn(6), &
          dajn(7), dajn(8), dajn(9), dajn(10),dajn(11),dajn(12), &
          dajn(13),dajn(14),dajn(15),dajn(16),dajn(17),dajn(18), &
          dajn(19)        / 1.08594539632967D-02, 8.53313194857091D-02, &
     -3.15277068113058D-01,-8.78420725294257D-02, 5.53251906976048D-02, &
      9.41674060503241D-03,-3.32187026018996D-03,-4.11157343156826D-04, &
      1.01297326891346D-04, 9.87633682208396D-06,-1.87312969812393D-06, &
     -1.50798500131468D-07, 2.32687669525394D-08, 1.59599917419225D-09, &
     -2.07665922668385D-10,-1.24103350500302D-11, 1.39631765331043D-12, &
      7.39400971155740D-14,-7.32887475627500D-15/

  data da(1),  da(2),  da(3),  da(4),  da(5),  da(6),  da(7), &
          da(8),  da(9),  da(10), da(11), da(12), da(13), da(14), &
          da(15)          / 4.91627321104601D-01, 3.11164930427489D-03, &
      8.23140762854081D-05,-4.61769776172142D-06,-6.13158880534626D-08, &
      2.87295804656520D-08,-1.81959715372117D-09,-1.44752826642035D-10, &
      4.53724043420422D-11,-3.99655065847223D-12,-3.24089119830323D-13, &
      1.62098952568741D-13,-2.40765247974057D-14, 1.69384811284491D-16, &
      8.17900786477396D-16/

  data db(1),  db(2),  db(3),  db(4),  db(5),  db(6),  db(7), &
          db(8),  db(9),  db(10), db(11), db(12), db(13), db(14), &
          db(15)          /-2.77571356944231D-01, 4.44212833419920D-03, &
     -8.42328522190089D-05,-2.58040318418710D-06, 3.42389720217621D-07, &
     -6.24286894709776D-09,-2.36377836844577D-09, 3.16991042656673D-10, &
     -4.40995691658191D-12,-5.18674221093575D-12, 9.64874015137022D-13, &
     -4.90190576608710D-14,-1.77253430678112D-14, 5.55950610442662D-15, &
     -7.11793337579530D-16/

  if ( x < 0.0D+00 ) then
    go to 90
  end if

  if ( 5.0D+00 < c ) then
    go to 60
  end if

  if ( 1.20D+00 < x ) then
    go to 30
  end if

  t = ( x + x - 1.2D+00 ) * con4
  tt = t + t
  j = n1
  f1 = ak1(j)
  f2 = 0.0D+00

  do i = 1, m1
    j = j - 1
    temp1 = f1
    f1 = tt * f1 - f2 + ak1(j)
    f2 = temp1
  end do

  ai = t * f1 - f2 + ak1(1)

  j = n1d
  f1 = dak1(j)
  f2 = 0.0D+00
  do i = 1, m1d
    j = j - 1
    temp1 = f1
    f1 = tt * f1 - f2 + dak1(j)
    f2 = temp1
  end do

  dai = -( t * f1 - f2 + dak1(1) )

  return

30 continue

  t = ( x + x - con2 ) * con3
  tt = t + t
  j = n2
  f1 = ak2(j)
  f2 = 0.0D+00

  do i = 1, m2
    j = j - 1
    temp1 = f1
    f1 = tt * f1 - f2 + ak2(j)
    f2 = temp1
  end do

  rtrx = sqrt ( rx )
  ec = exp ( -c )
  ai = ec * ( t * f1 - f2 + ak2(1) ) / rtrx
  j = n2d
  f1 = dak2(j)
  f2 = 0.0D+00

  do i = 1, m2d
    j = j - 1
    temp1 = f1
    f1 = tt * f1 - f2 + dak2(j)
    f2 = temp1
  end do

  dai = -ec * ( t * f1 - f2 + dak2(1) ) * rtrx
  return

60 continue

  t = 10.0D+00 / c - 1.0D+00
  tt = t + t
  j = n1
  f1 = ak3(j)
  f2 = 0.0D+00

  do i = 1, m1
    j = j - 1
    temp1 = f1
    f1 = tt * f1 - f2 + ak3(j)
    f2 = temp1
  end do

  rtrx = sqrt ( rx )
  ec = exp ( -c )
  ai = ec * ( t * f1 - f2 + ak3(1) ) / rtrx
  j = n1d
  f1 = dak3(j)
  f2 = 0.0D+00

  do i = 1, m1d
    j = j - 1
    temp1 = f1
    f1 = tt * f1 - f2 + dak3(j)
    f2 = temp1
  end do

  dai = -rtrx * ec * ( t * f1 - f2 + dak3(1) )
  return

90 continue

  if ( 5.0D+00 < c ) then
    go to 120
  end if

  t = 0.4D+00 * c - 1.0D+00
  tt = t + t
  j = n3
  f1 = ajp(j)
  e1 = ajn(j)
  f2 = 0.0D+00
  e2 = 0.0D+00

  do i = 1, m3
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt * f1 - f2 + ajp(j)
    e1 = tt * e1 - e2 + ajn(j)
    f2 = temp1
    e2 = temp2
  end do

  ai = ( t * e1 - e2 + ajn(1) ) - x * ( t * f1 - f2 + ajp(1) )
  j = n3d
  f1 = dajp(j)
  e1 = dajn(j)
  f2 = 0.0D+00
  e2 = 0.0D+00

  do i = 1, m3d
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt * f1 - f2 + dajp(j)
    e1 = tt * e1 - e2 + dajn(j)
    f2 = temp1
    e2 = temp2
  end do

  dai = x * x * ( t * f1 - f2 + dajp(1) ) + ( t * e1 - e2 + dajn(1) )

  return

120 continue

  t = 10.0D+00 / c - 1.0D+00
  tt = t + t
  j = n4
  f1 = a(j)
  e1 = b(j)
  f2 = 0.0D+00
  e2 = 0.0D+00

  do i = 1, m4
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt * f1 - f2 + a(j)
    e1 = tt * e1 - e2 + b(j)
    f2 = temp1
    e2 = temp2
  end do

  temp1 = t * f1 - f2 + a(1)
  temp2 = t * e1 - e2 + b(1)
  rtrx = sqrt ( rx )
  cv = c - fpi12
  ccv = cos ( cv )
  scv = sin ( cv )
  ai = ( temp1 * ccv - temp2 * scv ) / rtrx
  j = n4d
  f1 = da(j)
  e1 = db(j)
  f2 = 0.0D+00
  e2 = 0.0D+00

  do i = 1, m4d
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt * f1 - f2 + da(j)
    e1 = tt * e1 - e2 + db(j)
    f2 = temp1
    e2 = temp2
  end do

  temp1 = t * f1 - f2 + da(1)
  temp2 = t * e1 - e2 + db(1)
  e1 = ccv * con5 + 0.5D+00 * scv
  e2 = scv * con5 - 0.5D+00 * ccv
  dai = ( temp1 * e1 - temp2 * e2 ) * rtrx

  return
end
subroutine lltslv ( nr, n, a, x, b )

!*****************************************************************************80
!
!! LLTSLV solves A*x=b where A = L * L'.
!
!  Discussion:
!
!    L is a lower triangular matrix.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) A(NR,N), contains the lower triangular matrix L.
!
!    Output, real X(N), the solution vector.
!
!    Input, real B(N), the right hand side vector.  If B is not required by
!    the calling program, then B and X may share the same storage.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)
!
!  Forward solve, result in X.
!
  call forslv ( nr, n, a, x, b )
!
!  Back solve, result in X.
!
  call bakslv ( nr, n, a, x, x )

  return
end
subroutine lnsrch ( n, x, f, g, p, xpls, fpls, fcn, mxtake, iretcd, stepmx, &
  steptl, sx, ipr )

!*****************************************************************************80
!
!! LNSRCH finds a next Newton iterate by line search.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, sometimes called X[K-1].
!
!    Input, real ( kind = 8 ) F, the function value at the old iterate.
!
!    Input, real ( kind = 8 ) G(N), the gradient at the old iterate, or an
!    approximation to that value.
!
!    Input, real ( kind = 8 ) P(N), the (non-zero) Newton step.
!
!    Output, real ( kind = 8 ) XPLS(N), the new iterate.
!
!    Output, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Input, external FCN, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) x(n)
!      real ( kind = 8 ) f
!
!    Output, integer ( kind = 4 ) IRETCD, the return code.
!
!    Output, logical MXTAKE, is TRUE if a step of maximum size was used.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which successive
!    iterates are considered close enough to terminate algorithm.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
!  Local variables:
!
!    sln, the Newton length.
!
!    rln, the relative length of Newton step
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) almbda
  real ( kind = 8 ) b
  real ( kind = 8 ) disc
  real ( kind = 8 ) f
  external fcn
  real ( kind = 8 ) fpls
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) iretcd
  logical mxtake
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) pfpls
  real ( kind = 8 ) plmbda
  real ( kind = 8 ) rln
  real ( kind = 8 ) rmnlmb
  real ( kind = 8 ) scl
  real ( kind = 8 ) sln
  real ( kind = 8 ) slp
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) steptl
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) tlmbda
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xpls(n)

  mxtake = .false.
  iretcd = 2

  sln = sqrt ( sum ( ( sx(1:n) * p(1:n) )**2 ) )
!
!  Newton step longer than maximum allowed.
!
  if ( stepmx < sln ) then
    scl = stepmx / sln
    p(1:n) = p(1:n) * stepmx / sln
    sln = stepmx
  end if

  slp = dot_product ( g, p )

  rln = 0.0D+00
  do i = 1, n
    rln = max ( rln, abs ( p(i) ) / max ( abs ( x(i) ), 1.0D+00 / sx(i) ) )
  end do

  rmnlmb = steptl / rln
  almbda = 1.0D+00
!
!  Check if new iterate satisfactory.  Generate new lambda if necessary.
!
  do

    if ( iretcd < 2 ) then
      exit
    end if

    xpls(1:n) = x(1:n) + almbda * p(1:n)

    call fcn ( n, xpls, fpls )

    if ( f + slp * 0.0001D+00 * almbda < fpls ) then
      go to 130
    end if
!
!  Solution found.
!
    iretcd = 0

    if ( almbda == 1.0D+00 .and. 0.99D+00 * stepmx < sln ) then
      mxtake = .true.
    end if

    cycle
!
!  Solution not (yet) found.
!
130  continue
!
!  No satisfactory XPLS found sufficiently distinct from X.
!
    if ( almbda < rmnlmb ) then
      iretcd = 1
      cycle
    end if
!
!  Calculate new lambda.
!
!  First backtrack: quadratic fit.
!
    if ( almbda == 1.0D+00 ) then
      tlmbda = -slp / ( 2.0D+00 * ( fpls - f - slp ) )
      go to 170
    end if
!
!  All subsequent backtracks: cubic fit.
!
150 continue

    t1 = fpls - f - almbda * slp
    t2 = pfpls - f - plmbda * slp
    t3 = 1.0D+00 / ( almbda - plmbda )
    a = t3 * ( t1 / ( almbda * almbda ) - t2 / ( plmbda * plmbda ) )
    b = t3 * ( t2 *  almbda / ( plmbda * plmbda ) &
      - t1 * plmbda / ( almbda * almbda ) )
    disc = b * b - 3.0D+00 * a * slp

    if ( disc <= b * b ) then
      go to 160
    end if
!
!  Only one positive critical point, must be minimum.
!
    tlmbda = ( - b + sign ( 1.0D+00, a ) * sqrt ( disc ) ) / ( 3.0D+00 * a )
    go to 165
!
!  Both critical points positive, first is minimum.
!
160 continue

    tlmbda = ( -b - sign ( 1.0D+00, a ) * sqrt ( disc ) ) / ( 3.0D+00 * a )

165 continue

    if ( 0.5D+00 * almbda < tlmbda ) then
      tlmbda = 0.5D+00 * almbda
    end if

170 continue

    plmbda = almbda
    pfpls = fpls

    if ( almbda * 0.1D+00 <= tlmbda ) then
      almbda = tlmbda
    else
      almbda = almbda * 0.1D+00
    end if

  end do

  return
end
subroutine mvmltl ( nr, n, a, x, y )

!*****************************************************************************80
!
!! MVMLTL computes y = L * x where L is a lower triangular matrix stored in A.
!
!  Discussion:
!
!    Note that X and Y cannot share storage.
!
!  Modified:
!
!    29 May 2001
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) A(NR,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) Y(N), the result.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n
    y(i) = dot_product ( a(i,1:i), x(1:i) )
  end do

  return
end
subroutine mvmlts ( nr, n, a, x, y )

!*****************************************************************************80
!
!! MVMLTS computes y = A * x where A is a symmetric matrix.
!
!  Discussion:
!
!    A is a symmetric N by N matrix stored in its lower triangular part
!    and X and Y are N vectors.
!
!    X and Y cannot share storage.
!
!  Modified:
!
!    25 August 2001
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) A(NR,N), the symmetric N by N matrix.  The entries
!    of A are stored in the lower half of the array.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) Y(N), the result.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n

    y(i) = dot_product ( a(i,1:i), x(1:i) ) &
         + dot_product ( a(i+1:n,i), x(i+1:n) )
  end do

  return
end
subroutine mvmltu ( nr, n, a, x, y )

!*****************************************************************************80
!
!! MVMLTU computes y = L' * x where L is a lower triangular matrix.
!
!  Discussion:
!
!    Note that X and Y cannot share storage.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) A(NR,N), the N by N lower triangular matrix,
!
!    Input, real ( kind = 8 ) X(N), the matrix to be multiplied.
!
!    Output, real ( kind = 8 ) Y(N), the result vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n
    y(i) = dot_product ( x(i:n), a(i:n,i) )
  end do

  return
end
function numxer ( nerr )

!*****************************************************************************80
!
!! NUMXER returns the most recent error number.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NERR, the most recent error number.
!
!    Output, integer ( kind = 4 ) NUMXER, the most recent error number.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) numxer

  nerr = j4save ( 1, 0, .false. )
  numxer = nerr

  return
end
subroutine optchk ( n, x, typsiz, sx, fscale, gradtl, itnlim, ndigit, epsm, &
  dlt, method, iexp, iagflg, iahflg, stepmx, msg, ipr )

!*****************************************************************************80
!
!! OPTCHK checks the input to the optimization routine.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), an approximate solution of the problem.
!
!    Input/output, real ( kind = 8 ) TYPSIZ(N), a typical size for each
!    component of X.  If TYPSIZ(I) is zero, it is reset to 1.
!
!    Input, real ( kind = 8 ) SX(N), the  diagonal scaling matrix for X.
!
!    Input/output, real ( kind = 8 ) FSCALE, an estimate of the scale of
!    the objective function FCN.
!
!    Input, real ( kind = 8 ) GRADTL, the tolerance at which the gradient
!    is considered close enough to zero to terminate the algorithm.
!
!    Input/output, integer ( kind = 4 ) ITNLIM, the maximum number of allowable
!    iterations.
!
!    Input/output, integer ( kind = 4 ) NDIGIT, the number of good digits in
!    optimization function FCN.
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Input/output, integer ( kind = 4 ) METHOD, the algorithm indicator.
!
!    Input/output, integer ( kind = 4 ) IEXP, the expense flag.
!
!    Input/output, integer ( kind = 4 ) IAGFLG, = 1 if analytic gradient supplied.
!
!    Input/output, integer ( kind = 4 ) IAHFLG, = 1 if analytic hessian supplied.
!
!    Input/output, real ( kind = 8 ) STEPMX, the maximum step size.
!
!    Input/output, integer ( kind = 4 ) MSG, the message and error code.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dlt
  real ( kind = 8 ) epsm
  real ( kind = 8 ) fscale
  real ( kind = 8 ) gradtl
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iagflg
  integer ( kind = 4 ) iahflg
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) itnlim
  integer ( kind = 4 ) method
  integer ( kind = 4 ) msg
  integer ( kind = 4 ) ndigit
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) stpsiz
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) typsiz(n)
  real ( kind = 8 ) x(n)
!
!  Check that parameters only take on acceptable values.
!  if not, set them to default values.
!
  if ( method < 1 .or. 3 < method ) then
    method = 1
  end if

  if ( iagflg /= 1 ) then
    iagflg = 0
  end if

  if ( iahflg /= 1 ) then
    iahflg = 0
  end if

  if ( iexp /= 0 ) then
    iexp = 1
  end if

  if ( mod ( msg/2, 2 ) == 1 .and. iagflg == 0 ) then
    write ( ipr, 906 ) msg, iagflg
    msg = -6
    return
  end if

  if ( mod ( msg/4, 2 ) == 1 .and. iahflg == 0 ) then
    write ( ipr, 907 ) msg, iahflg
    msg = -7
    return
  end if
!
!  Check N.
!
  if ( n <= 0 ) then
    write ( ipr, * ) ' '
    write ( ipr, * ) 'OPTCHK - Fatal error!'
    write ( ipr, * ) '  Illegal nonpositive value of N = ', n
    msg = -1
    return
  end if

  if ( n == 1 .and. mod ( msg, 2 ) == 0 ) then
    write ( ipr, 902 )
    msg = -2
    return
  end if
!
!  Compute the scale matrix.
!
  do i = 1, n
    if ( typsiz(i) == 0.0D+00 ) then
      typsiz(i) = 1.0D+00
    end if
  end do

  typsiz(1:n) = abs ( typsiz(1:n) )
  sx(1:n) = 1.0D+00 / typsiz(1:n)
!
!  Check maximum step size.
!
  if ( stepmx <= 0.0D+00 ) then

    stpsiz = sqrt ( sum ( x(1:n)**2 * sx(1:n)**2 ) )

    stepmx = max ( 1.0D+03 * stpsiz, 1.0D+03 )

  end if
!
!  Check the function scale.
!
  if ( fscale == 0.0D+00 ) then
    fscale = 1.0D+00
  end if

  if ( fscale < 0.0D+00 ) then
    fscale = -fscale
  end if
!
!  Check gradient tolerance
!
  if ( gradtl < 0.0D+00 ) then
    write ( ipr, 903 ) gradtl
    msg = -3
    return
  end if
!
!  Check iteration limit
!
  if ( itnlim <= 0 ) then
    write ( ipr, 904 ) itnlim
    msg = -4
    return
  end if
!
!  Check number of digits of accuracy in function FCN.
!
  if ( ndigit == 0 ) then
    write ( ipr, 905 ) ndigit
    msg = -5
    return
  end if

  if ( ndigit < 0 ) then
    ndigit = -log10 ( epsm )
  end if
!
!  Check trust region radius.
!
  if ( dlt <= 0.0D+00 ) then
    dlt = -1.0D+00
  end if

  if ( stepmx < dlt ) then
    dlt = stepmx
  end if

  902 format(' optchk    +++ warning +++  this package is inefficient', &
    'for problems of size n=1.'/ &
    ' optchk    check installation libraries for more appropriate routines.'/ &
    ' optchk    if none, set msg and resubmit.')
  903 format(' optchk    illegal tolerance.  gradtl=',e20.13)
  904 format(' optchk    illegal iteration limit.  itnlim=',i5)
  905 format(' optchk    minimization function has no good digits.', &
     'ndigit=',i5)
  906 format(' optchk    user requests that analytic gradient be', &
     ' accepted as properly coded (msg=',i5, '),'/ &
     ' optchk    but analytic gradient not supplied (iagflg=',i5, ').')
  907 format(' optchk    user requests that analytic hessian be', &
     ' accepted as properly coded (msg=',i5, '),'/ &
     ' optchk    but analytic hessian not supplied(iahflg=',i5, ').')

  return
end
subroutine optdrv ( nr, n, x, fcn, d1fcn, d2fcn, typsiz, fscale, method, &
  iexp, msg, ndigit, itnlim, iagflg, iahflg, ipr, dlt, gradtl, stepmx, &
  steptl, xpls, fpls, gpls, itrmcd, a, udiag, g, p, sx, wrk0, wrk1, wrk2, &
  wrk3 )

!*****************************************************************************80
!
!! OPTDRV is a driver for the nonlinear optimization package.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, an rough solution of
!    the problem.  On output, the computed solution.
!
!    Input, external FCN, the name of the subroutine that evaluates
!    the optimization function, of the form:
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!    Input, external D1FCN, the name of the subroutine to evaluate gradient
!    of FCN, of the form:
!      subroutine d1fcn ( n, x, g )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) g(n)
!      real ( kind = 8 ) x(n)
!
!    Input, external D2FCN, the name of the subroutine to evaluate the
!    Hessian of the function, of the form:
!      subroutine d2fcn ( nr, n, x, h )
!      integer ( kind = 4 ) nr
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) h(nr,n)
!      real ( kind = 8 ) x(n)
!
!    Input, real ( kind = 8 ) TYPSIZ(N), a typical size for each component of X.
!
!    Input, real ( kind = 8 ) FSCALE, an estimate of the scale of the objective
!    function.
!
!    Input, integer ( kind = 4 ) METHOD, the algorithm to use to solve
!    minimization problem:
!    1, line search
!    2, double dogleg
!    3, More-Hebdon
!
!    Input, integer ( kind = 4 ) IEXP, function expense flag.
!    Set IEXP to 1 if optimization function fcn is expensive to
!    evaluate,  and 0 otherwise.  If set then hessian will
!    be evaluated by secant update instead of
!    analytically or by finite differences.
!
!    Input/output, integer ( kind = 4 ) MSG.
!    On input, set it positive to inhibit certain automatic checks
!    On output. < 0 indicates an error occurred.
!
!    Input, integer ( kind = 4 ) NDIGIT, the number of good digits in optimization
!    function fcn.
!
!    Input, integer ( kind = 4 ) ITNLIM, the maximum number of allowable iterations.
!
!    Input, integer ( kind = 4 ) IAGFLG, is 1 if analytic gradient supplied.
!
!    Input, integer ( kind = 4 ) IAHFLG, is 1 if analytic hessian supplied.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Input, real ( kind = 8 ) DLT, the trust region radius.
!
!    Input, real ( kind = 8 ) GRADTL, the tolerance at which gradient
!    considered close enough to zero to terminate algorithm.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which
!    successive iterates are considered close enough to terminate algorithm
!
!    Input/output, real ( kind = 8 ) XPLS(N); on exit, XPLS is the local
!    minimizer.
!
!    Input/output, real ( kind = 8 ) FPLS; on exit, the function value at XPLS.
!
!    Input/output, real ( kind = 8 ) GPLS(N); on exit, the gradient at XPLS.
!
!    Output, integer ( kind = 4 ) ITRMCD, the termination code.
!
!    Workspace, real ( kind = 8 ) A(NR,N), workspace for hessian (or estimate)
!    and its Cholesky decomposition.
!
!    Workspace, real ( kind = 8 ) UDIAG(N), workspace for diagonal of hessian.
!
!    Workspace, real ( kind = 8 ) G(N), workspace for gradient at current
!    iterate.
!
!    Workspace, real ( kind = 8 ) P(N), workspace for the step.
!
!    Workspace, real ( kind = 8 ) SX(N), workspace for diagonal scaling matrix.
!
!    Workspace, real ( kind = 8 ) WRK0(N), WRK1(N), WRK2(N), WRK3(N).
!
!  Local variables:
!
!    analtl, tolerance for gradient and hessian checking.
!
!    epsm, machine epsilon.
!
!    f, function value: fcn(x).
!
!    itncnt, current iteration, k
!
!    rnf, relative noise in optimization function fcn.
!
!    noise=10.**(-ndigit)
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) amu
  real ( kind = 8 ) amusav
  real ( kind = 8 ) analtl
  external d1fcn
  external d2fcn
  real ( kind = 8 ) dlpsav
  real ( kind = 8 ) dlt
  real ( kind = 8 ) dltp
  real ( kind = 8 ) dltsav
  real ( kind = 8 ) epsm
  real ( kind = 8 ) f
  external fcn
  real ( kind = 8 ) fpls
  real ( kind = 8 ) fscale
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) gpls(n)
  real ( kind = 8 ) gradtl
  integer ( kind = 4 ) iagflg
  integer ( kind = 4 ) iahflg
  integer ( kind = 4 ) icscmx
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) iretcd
  integer ( kind = 4 ) itncnt
  integer ( kind = 4 ) itnlim
  integer ( kind = 4 ) itrmcd
  integer ( kind = 4 ) method
  integer ( kind = 4 ) msg
  logical mxtake
  integer ( kind = 4 ) ndigit
  logical noupdt
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) phi
  real ( kind = 8 ) phip0
  real ( kind = 8 ) phisav
  real ( kind = 8 ) phpsav
  real ( kind = 8 ) rnf
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) steptl
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) typsiz(n)
  real ( kind = 8 ) udiag(n)
  real ( kind = 8 ) value(1)
  real ( kind = 8 ) wrk(1)
  real ( kind = 8 ) wrk0(n)
  real ( kind = 8 ) wrk1(n)
  real ( kind = 8 ) wrk2(n)
  real ( kind = 8 ) wrk3(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xpls(n)
!
!  Initialization.
!
  p(1:n) = 0.0D+00
  itncnt = 0
  iretcd = -1
  epsm = epsilon ( epsm )

  call optchk ( n, x, typsiz, sx, fscale, gradtl, itnlim, ndigit, epsm, &
    dlt, method, iexp, iagflg, iahflg, stepmx, msg, ipr )

  if ( msg < 0 ) then
    return
  end if

  rnf = max ( 10.0D+00**(-ndigit), epsm )

  analtl = max ( 1.0D-02, sqrt ( rnf ) )

  if ( mod ( msg / 8, 2 ) == 0 ) then
    write ( ipr, 901 )
    write ( ipr, 900 ) typsiz(1:n)
    write ( ipr, 902 )
    write ( ipr, 900 ) sx(1:n)
    write ( ipr, 903 ) fscale
    write ( ipr, 904 ) ndigit,iagflg,iahflg,iexp,method,itnlim,epsm
    write ( ipr, 905 ) stepmx,steptl,gradtl,dlt,rnf,analtl
  end if
!
!  Evaluate fcn(x)
!
  call fcn ( n, x, f )
!
!  Evaluate analytic or finite difference gradient and check analytic
!  gradient, if requested.
!
  if ( iagflg /= 1 ) then

    value(1) = f
    call fstofd ( 1, 1, n, x, fcn, value, g, sx, rnf, wrk, 1 )

  else

    call d1fcn ( n, x, g )

    if ( mod ( msg/2, 2 ) /= 1 ) then

      call grdchk ( n, x, fcn, f, g, typsiz, sx, fscale, rnf, analtl, wrk1, &
        msg, ipr )

      if ( msg < 0 ) then
        return
      end if

    end if

  end if

  call optstp ( n, x, f, g, wrk1, itncnt, icscmx, itrmcd, gradtl, steptl, &
    sx, fscale, itnlim, iretcd, mxtake, ipr, msg )

  if ( itrmcd /= 0 ) then
    go to 700
  end if

  if ( iexp /= 1 ) then
    go to 80
  end if
!
!  If optimization function expensive to evaluate (iexp=1), then
!  hessian will be obtained by secant updates.  Get initial hessian.
!
  call hsnint ( nr, n, a, sx, method )
  go to 90

80 continue
!
!  Evaluate analytic or finite difference hessian and check analytic
!  hessian if requested (only if user-supplied analytic hessian
!  routine d2fcn fills only lower triangular part and diagonal of a).
!
  if ( iahflg == 1 ) then
    go to 82
  end if

  if ( iagflg == 1 ) then
    call fstofd ( nr, n, n, x, d1fcn, g, a, sx, rnf, wrk1, 3 )
  else
    call sndofd ( nr, n, x, fcn, f, a, sx, rnf, wrk1, wrk2 )
  end if

  go to 88

82 continue

  if ( mod ( msg / 4, 2 ) == 0 ) then
    go to 85
  end if

  call d2fcn ( nr, n, x, a )
  go to 88

85 continue

  call heschk ( nr, n, x, fcn, d1fcn, d2fcn, f, g, a, typsiz, &
    sx, rnf, analtl, iagflg, udiag, wrk1, wrk2, msg, ipr )
!
!  HESCHK evaluates d2fcn and checks it against the finite
!  difference hessian which it calculates by calling fstofd
!  (if iagflg == 1) or sndofd (otherwise).
!
  if ( msg < 0 ) then
    return
  end if

88 continue

90 continue

  if ( mod ( msg / 8, 2 ) == 0 ) then
    call result ( nr, n, x, f, g, a, p, itncnt, 1, ipr )
  end if
!
!  iteration
!
  100 continue

  itncnt = itncnt + 1
!
!  Find perturbed local model hessian and its ll+ decomposition
!  (skip this step if line search or dogstep techniques being used with
!  secant updates.  Cholesky decomposition l already obtained from
!  secfac.)
!
  if ( iexp == 1 .and. method /= 3 ) then
    go to 105
  end if

  103   continue

  call chlhsn ( nr, n, a, epsm, sx, udiag )
  105 continue
!
!  Solve for Newton step:  ap = -g
!
  wrk1(1:n) = - g(1:n)

  call lltslv ( nr, n, a, p, wrk1 )
!
!  Decide whether to accept Newton step  xpls = x + p
!  or to choose xpls by a global strategy.
!
  if ( iagflg == 0 .and. method /= 1 ) then

    dltsav = dlt

    if ( method /= 2 ) then
      amusav = amu
      dlpsav = dltp
      phisav = phi
      phpsav = phip0
    end if

  end if

  if ( method == 1 ) then

    call lnsrch ( n, x, f, g, p, xpls, fpls, fcn, mxtake, iretcd, &
      stepmx, steptl, sx, ipr )

  else if ( method == 2 ) then

    call dogdrv ( nr, n, x, f, g, a, p, xpls, fpls, fcn, sx, stepmx, &
      steptl, dlt, iretcd, mxtake, wrk0, wrk1, wrk2, wrk3, ipr )

  else if ( method == 3 ) then

    call hookdr ( nr, n, x, f, g, a, udiag, p, xpls, fpls, fcn, sx, stepmx, &
      steptl, dlt, iretcd, mxtake, amu, dltp, phi, phip0, wrk0, &
      wrk1, wrk2, epsm, itncnt, ipr )

  end if
!
!  If could not find satisfactory step and forward difference
!  gradient was used, retry using central difference gradient.
!
  if ( iretcd /= 1 .or. iagflg /= 0 ) then
    go to 112
  end if
!
!  Set iagflg for central differences.
!
     iagflg = -1
     write(ipr,906) itncnt
     call fstocd ( n, x, fcn, sx, rnf, g )

     if ( method == 1 ) then
       go to 105
     end if

     dlt = dltsav

     if ( method == 2 ) then
       go to 105
     end if

     amu = amusav
     dltp = dlpsav
     phi = phisav
     phip0 = phpsav
     go to 103
!
!  Calculate step for output
!
  112 continue

  p(1:n) = xpls(1:n) - x(1:n)
!
!  Calculate the gradient at XPLS.
!
  if ( iagflg == -1 ) then
    call fstocd ( n, xpls, fcn, sx, rnf, gpls )
  else if ( iagflg == 0 ) then
    value(1) = fpls
    call fstofd ( 1, 1, n, xpls, fcn, value, gpls, sx, rnf, wrk, 1 )
  else
    call d1fcn ( n, xpls, gpls )
  end if
!
!  Check whether stopping criteria satisfied.
!
  call optstp ( n, xpls, fpls, gpls, x, itncnt, icscmx, itrmcd, gradtl, &
    steptl, sx, fscale, itnlim, iretcd, mxtake, ipr, msg )

  if ( itrmcd /= 0 ) then
    go to 690
  end if
!
!  Evaluate hessian at xpls
!
  if ( iexp == 0 ) then
    go to 130
  end if

  if ( method == 3 ) then
     call secunf ( nr, n, x, g, a, udiag, xpls, gpls, epsm, itncnt, rnf, &
       iagflg, noupdt, wrk1, wrk2, wrk3 )
  else
    call secfac ( nr, n, x, g, a, xpls, gpls, epsm, itncnt, rnf, iagflg, &
      noupdt, wrk0, wrk1, wrk2, wrk3 )
  end if

  go to 150

  130 continue

  if ( iahflg == 1 ) then
    go to 140
  end if

  if ( iagflg == 1 ) then
    call fstofd ( nr, n, n, xpls, d1fcn, gpls, a, sx, rnf, wrk1, 3 )
  else
    call sndofd ( nr, n, xpls, fcn, fpls, a, sx, rnf, wrk1, wrk2 )
  end if

  go to 150

  140 continue

  call d2fcn ( nr, n, xpls, a )

  150 continue

  if ( mod ( msg / 16, 2 ) == 1 ) then
    call result ( nr, n, xpls, fpls, gpls, a, p, itncnt, 1, ipr )
  end if
!
!  x <-- xpls
!  g <-- gpls
!  f <-- fpls
!
  f = fpls
  x(1:n) = xpls(1:n)
  g(1:n) = gpls(1:n)

  go to 100
!
!  Termination.
!
!  Reset XPLS, FPLS, GPLS, if previous iterate solution
!
  690 if ( itrmcd /= 3 ) then
    go to 710
  end if

  700 continue

  fpls = f
  xpls(1:n) = x(1:n)
  gpls(1:n) = g(1:n)
!
!  Print results
!
  710 continue

  if ( mod ( msg / 8, 2 ) == 0) then
    call result ( nr, n, xpls, fpls, gpls, a, p, itncnt, 0, ipr )
  end if

  msg = 0

  900 format(' optdrv       ',5(e20.13,3x))
  901 format('0optdrv    typical x')
  902 format(' optdrv    diagonal scaling matrix for x')
  903 format(' optdrv    typical f =',e20.13)
  904 format('0optdrv    number of good digits in fcn=',i5/ &
             ' optdrv    gradient flag  =',i5,'   (=1 if analytic', &
             ' gradient supplied)'/ &
             ' optdrv    hessian flag   =',i5,'   (=1 if analytic', &
             ' hessian supplied)'/ &
             ' optdrv    expense flag   =',i5, '   (=1 if', &
             ' minimization function expensive to evaluate)'/ &
             ' optdrv    method to use  =',i5,'   (=1,2,3 for line', &
             ' search, double dogleg, more-hebdon respectively)'/ &
             ' optdrv    iteration limit=',i5/ &
             ' optdrv    machine epsilon=',e20.13)

  905 format('0optdrv    maximum step size =',e20.13/ &
             ' optdrv    step tolerance    =',e20.13/ &
             ' optdrv    gradient tolerance=',e20.13/ &
             ' optdrv    trust reg radius  =',e20.13/ &
             ' optdrv    rel noise in fcn  =',e20.13/ &
             ' optdrv    anal-fd tolerance =',e20.13)

  906 format(' optdrv    shift from forward to central differences', &
     ' in iteration ', i5)

  return
end
subroutine optif0 ( n, x, fcn, xpls, fpls, gpls, itrmcd )

!*****************************************************************************80
!
!! OPTIF0 provides a simple interface to the minimization package.
!
!  Modified:
!
!    29 May 2001
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, an rough solution of
!    the problem.  On output, the computed solution.
!
!    Input, external FCN, the name of the subroutine that evaluates
!    the optimization function, of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!    Output, real ( kind = 8 ) XPLS(N), estimated local minimizer of
!    the function.
!
!    Output, real ( kind = 8 ) FPLS, the function value at XPLS.
!
!    Output, real ( kind = 8 ) GPLS(N), the gradient at XPLS.
!
!    Output, integer ( kind = 4 ) ITRMCD, the termination code.
!    1, relative gradient close to zero.
!       The current iterate is probably solution.
!    2, successive iterates within tolerance.
!       The current iterate is probably solution.
!    3, the last global step failed to locate a point lower than X.
!       Either x is an approximate local minimum of the function,
!       the function is too non-linear for this algorithm,
!       or STEPTL is too large.
!    4, iteration limit exceeded.  The algorithm failed.
!    5, maximum step size exceeded 5 consecutive times.
!       Either the function is unbounded below, becomes asymptotic to a
!       finite value from above in some direction, or STEPMX is too small.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  external d1fcn
  external d2fcn
  real ( kind = 8 ) dlt
  external fcn
  real ( kind = 8 ) fscale
  real ( kind = 8 ) fpls
  real ( kind = 8 ) gpls(n)
  real ( kind = 8 ) gradtl
  integer ( kind = 4 ) iagflg
  integer ( kind = 4 ) iahflg
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) itnlim
  integer ( kind = 4 ) itrmcd
  integer ( kind = 4 ) method
  integer ( kind = 4 ) msg
  integer ( kind = 4 ) ndigit
  integer ( kind = 4 ) nr
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) steptl
  real ( kind = 8 ) wrk(n,9)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xpls(n)
!
! equivalence wrk(n,1) = udiag(n)
!             wrk(n,2) = g(n)
!             wrk(n,3) = p(n)
!             wrk(n,4) = typsiz(n)
!             wrk(n,5) = sx(n)
!             wrk(n,6) = wrk0(n)
!             wrk(n,7) = wrk1(n)
!             wrk(n,8) = wrk2(n)
!             wrk(n,9) = wrk3(n)
!
  nr = n

  call dfault ( n, x, wrk(1,4), fscale, method, iexp, msg, ndigit, &
    itnlim, iagflg, iahflg, ipr, dlt, gradtl, stepmx, steptl )

  call optdrv ( nr, n, x, fcn, d1fcn, d2fcn, wrk(1,4), fscale, &
    method, iexp, msg, ndigit, itnlim, iagflg, iahflg, ipr, &
    dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, &
    a, wrk(1,1), wrk(1,2), wrk(1,3), wrk(1,5), wrk(1,6), &
    wrk(1,7), wrk(1,8), wrk(1,9) )

  return
end
subroutine optstp ( n, xpls, fpls, gpls, x, itncnt, icscmx, &
  itrmcd, gradtl, steptl, sx, fscale, itnlim, iretcd, mxtake, ipr, msg )

!*****************************************************************************80
!
!! OPTSTP: unconstrained minimization stopping criteria
!
!  Discussion:
!
!    OPSTP determines whether the optimization algorithm should terminate,
!    due to any of the following:
!    1) the problem has been solved to the user's tolerance;
!    2) convergence within user tolerance;
!    3) iteration limit reached;
!    4) divergence or too restrictive maximum step (stepmx) suspected;
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) XPLS(N), the new iterate, X[K].
!
!    Input, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Input, real ( kind = 8 ) GPLS(N), the gradient at the new iterate, or an
!    approximation of that value.
!
!    Input, real ( kind = 8 ) X(N), the old iterate X[K-1].
!
!    Input, integer ( kind = 4 ) ITNCNT, the current iteration K.
!
!    Input/output, integer ( kind = 4 ) ICSCMX, the number of consecutive steps
!    greater than or equal to STEPMX.
!    [retain value between successive calls].
!
!    Output, integer ( kind = 4 ) ITRMD, the termination code.
!
!    Input, real ( kind = 8 ) GRADTL, the tolerance at which relative gradient
!    considered close enough to zero to terminate algorithm.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which
!    successive iterates are considered close enough to terminate algorithm.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) FSCALE, an estimate of the scale of objective
!    function.
!
!    Input, integer ( kind = 4 ) ITNLIM, the maximum number of allowable iterations.
!
!    Output, integer ( kind = 4 ) IRETCD, the return code.
!
!    Input, logical MXTAKE, TRUE if a step of maximum length was used.
!
!    Output, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Input, integer ( kind = 4 ) MSG, if includes a term 8, suppress output.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  real ( kind = 8 ) fpls
  real ( kind = 8 ) fscale
  real ( kind = 8 ) gpls(n)
  real ( kind = 8 ) gradtl
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icscmx
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) iretcd
  integer ( kind = 4 ) itncnt
  integer ( kind = 4 ) itnlim
  integer ( kind = 4 ) itrmcd
  integer ( kind = 4 ) jtrmcd
  integer ( kind = 4 ) msg
  logical mxtake
  real ( kind = 8 ) relgrd
  real ( kind = 8 ) relstp
  real ( kind = 8 ) rgx
  real ( kind = 8 ) rsx
  real ( kind = 8 ) steptl
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xpls(n)

  itrmcd = 0
!
!  Last global step failed to locate a point lower than X.
!
  if ( iretcd == 1 ) then
    jtrmcd = 3
    go to 600
  end if
!
!  Find direction in which relative gradient maximum.
!  Check whether within tolerance
!
  d = max ( abs ( fpls ), fscale )

  rgx = 0.0D+00
  do i = 1, n
    relgrd = abs ( gpls(i) ) * max ( abs ( xpls(i) ), 1.0D+00 / sx(i) ) / d
    rgx = max ( rgx, relgrd )
  end do

  jtrmcd = 1
  if ( rgx <= gradtl ) then
    go to 600
  end if

  if ( itncnt == 0 ) then
    return
  end if
!
!  Find direction in which relative stepsize is maximum.
!  Check whether within tolerance.
!
  rsx = 0.0D+00
  do i = 1, n
    relstp = abs ( xpls(i) - x(i) ) / max ( abs ( xpls(i) ), 1.0D+00 / sx(i) )
    rsx = max ( rsx, relstp )
  end do

  jtrmcd = 2
  if ( rsx <= steptl ) then
    go to 600
  end if
!
!  Check iteration limit.
!
  jtrmcd = 4
  if ( itnlim <= itncnt ) then
    go to 600
  end if
!
!  Check number of consecutive steps \ stepmx
!
  if ( .not. mxtake ) then
    icscmx = 0
    return
  else
    if ( mod ( msg / 8, 2 ) == 0 ) then
      write ( ipr, 900 )
    end if
    icscmx = icscmx + 1
    if ( icscmx < 5 ) then
      return
    end if
    jtrmcd = 5
  end if
!
!  Print termination code
!
  600 continue

  itrmcd = jtrmcd

  if ( itrmcd == 1 ) then
    write ( ipr, 901 )
  else if ( itrmcd == 2 ) then
    write(ipr,902)
  else if ( itrmcd == 3 ) then
    write(ipr,903)
  else if ( itrmcd == 4 ) then
    write(ipr,904)
  else if ( itrmcd == 5 ) then
    write(ipr,905)
  end if

  900 format('0optstp    step of maximum length (stepmx) taken')
  901 format('0optstp    relative gradient close to zero.'/ &
             ' optstp    current iterate is probably solution.')
  902 format('0optstp    successive iterates within tolerance.'/ &
             ' optstp    current iterate is probably solution')
  903 format('0optstp    last global step failed to locate a point', &
             ' lower than x.'/ &
             ' optstp    either x is an approximate local minimum', &
             ' of the function',/ &
             ' optstp    the function is too non-linear for this algorithm,'/ &
             ' optstp    or steptl is too large.')
  904 format('optstp    iteration limit exceeded.'/'optstp    algorithm failed.')
  905 format('0optstp    maximum step size exceeded 5 consecutive times.'/ &
             ' optstp    either the function is unbounded below',/ &
             ' optstp    becomes asymptotic to a finite value', &
             ' from above in some direction',/ &
             ' optstp    or stepmx is too small')

  return
end
subroutine passb ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )

!*****************************************************************************80
!
!! PASSB is a lower level routine used by CFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) c1(ido,l1,ip)
  real ( kind = 8 ) c2(idl1,ip)
  real ( kind = 8 ) cc(ido,ip,l1)
  real ( kind = 8 ) ch(ido,l1,ip)
  real ( kind = 8 ) ch2(idl1,ip)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idj
  integer ( kind = 4 ) idl
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) idp
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nac
  integer ( kind = 4 ) nt
  real ( kind = 8 ) wa(*)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

  nt = ip * idl1
  ipph = ( ip + 1 ) / 2
  idp = ip * ido

  if ( l1 <= ido ) then

    do j = 2, ipph
      jc = ip + 2 - j
      do k = 1, l1
        ch(1:ido,k,j)  = cc(1:ido,j,k) + cc(1:ido,jc,k)
        ch(1:ido,k,jc) = cc(1:ido,j,k) - cc(1:ido,jc,k)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      do i = 1, ido
        ch(i,1:l1,j)  = cc(i,j,1:l1) + cc(i,jc,1:l1)
        ch(i,1:l1,jc) = cc(i,j,1:l1) - cc(i,jc,1:l1)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l) = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =            wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j
      idlj = idlj + inc
      if ( idp < idlj ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) + wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  nac = 1

  if ( ido == 2 ) then
    return
  end if

  nac = 0
  c2(1:idl1,1) = ch2(1:idl1,1)
  c1(1:2,1:l1,2:ip) = ch(1:2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) - wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   + wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) - wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   + wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passb2 ( ido, l1, cc, ch, wa1 )

!*****************************************************************************80
!
!! PASSB2 is a lower level routine used by CFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,2,l1)
  real ( kind = 8 ) ch(ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa1(ido)

  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)
        ch(i,k,1)   = cc(i,1,k)   + cc(i,2,k)
        ti2         = cc(i,1,k)   - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 + wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 - wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passb3 ( ido, l1, cc, ch, wa1, wa2 )

!*****************************************************************************80
!
!! PASSB3 is a lower level routine used by CFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,3,l1)
  real ( kind = 8 ) ch(ido,l1,3)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: taui = 0.866025403784439D+00
  real ( kind = 8 ), parameter :: taur = -0.5D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k) - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passb4 ( ido, l1, cc, ch, wa1, wa2, wa3 )

!*****************************************************************************80
!
!! PASSB4 is a lower level routine used by CFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,4,l1)
  real ( kind = 8 ) ch(ido,l1,4)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,4,k) - cc(2,2,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,2,k) - cc(1,4,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k) - cc(i,3,k)
        ti2 = cc(i,1,k) + cc(i,3,k)
        ti3 = cc(i,2,k) + cc(i,4,k)
        tr4 = cc(i,4,k) - cc(i,2,k)

        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,2,k) - cc(i-1,4,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i,k,1) = ti2 + ti3
        ci3 = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 - wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 + wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 - wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 + wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 - wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 + wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passb5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )

!*****************************************************************************80
!
!! PASSB5 is a lower level routine used by CFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,5,l1)
  real ( kind = 8 ) ch(ido,l1,5)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: ti11 = 0.951056516295154D+00
  real ( kind = 8 ), parameter :: ti12 = 0.587785252292473D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: tr11 = 0.309016994374947D+00
  real ( kind = 8 ), parameter :: tr12 = -0.809016994374947D+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 - wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 + wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 - wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 + wa4(i) * dr5

      end do
    end do

  end if

  return
end
subroutine passf ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )

!*****************************************************************************80
!
!! PASSF is a lower level routine used by CFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) c1(ido,l1,ip)
  real ( kind = 8 ) c2(idl1,ip)
  real ( kind = 8 ) cc(ido,ip,l1)
  real ( kind = 8 ) ch(ido,l1,ip)
  real ( kind = 8 ) ch2(idl1,ip)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idj
  integer ( kind = 4 ) idl
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) idp
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nac
  integer ( kind = 4 ) nt
  real ( kind = 8 ) wa(*)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

  nt = ip * idl1
  ipph = ( ip + 1 ) / 2
  idp = ip * ido

  if ( l1 <= ido ) then

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l)  = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =           - wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j

      idlj = idlj + inc
      if ( idp < idlj ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) - wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  if ( ido == 2 ) then
    nac = 1
    return
  end if

  nac = 0

  c2(1:idl1,1)    = ch2(1:idl1,1)
  c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)
  c1(2,1:l1,2:ip) = ch(2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) + wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   - wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) + wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   - wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passf2 ( ido, l1, cc, ch, wa1 )

!*****************************************************************************80
!
!! PASSF2 is a lower level routine used by CFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,2,l1)
  real ( kind = 8 ) ch(ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa1(ido)

  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)

        ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
        ti2       = cc(i,1,k) - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 - wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 + wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passf3 ( ido, l1, cc, ch, wa1, wa2 )

!*****************************************************************************80
!
!! PASSF3 is a lower level routine used by CFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,3,l1)
  real ( kind = 8 ) ch(ido,l1,3)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: taui = -0.866025403784439D+00
  real ( kind = 8 ), parameter :: taur = -0.5D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k)   - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )

!*****************************************************************************80
!
!! PASSF4 is a lower level routine used by CFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,4,l1)
  real ( kind = 8 ) ch(ido,l1,4)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,2,k) - cc(2,4,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,4,k) - cc(1,2,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k)   - cc(i,3,k)
        ti2 = cc(i,1,k)   + cc(i,3,k)
        ti3 = cc(i,2,k)   + cc(i,4,k)
        tr4 = cc(i,2,k)   - cc(i,4,k)
        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,4,k) - cc(i-1,2,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3         = tr2 - tr3
        ch(i,k,1)   = ti2 + ti3
        ci3         = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 + wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 - wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 + wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 - wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 + wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 - wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )

!*****************************************************************************80
!
!! PASSF5 is a lower level routine used by CFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,5,l1)
  real ( kind = 8 ) ch(ido,l1,5)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: ti11 = -0.951056516295154D+00
  real ( kind = 8 ), parameter :: ti12 = -0.587785252292473D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ), parameter :: tr11 =  0.309016994374947D+00
  real ( kind = 8 ), parameter :: tr12 = -0.809016994374947D+00
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 + wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 - wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 + wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 - wa4(i) * dr5

      end do
    end do

  end if

  return
end
subroutine pchce ( ic, vc, n, x, h, slope, d, incfd, ierr )

!*****************************************************************************80
!
!! PCHCE is called by PCHIC to set end derivatives as requested by the user.
!
!  Discussion:
!
!    PCHCE must be called after interior derivative values have been set.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D array.
!
!    One could reduce the number of arguments and amount of local
!    storage, at the expense of reduced code clarity, by passing in
!    the array WK, rather than splitting it into H and SLOPE, and
!    increasing its length enough to incorporate STEMP and XTEMP.
!
!    The two monotonicity checks only use the sufficient conditions.
!    thus, it is possible (but unlikely) for a boundary condition to
!    be changed, even though the original interpolant was monotonic.
!    At least the result is a continuous function of the data.
!
!  Author:
!
!    FORTRAN77 original version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IC(2), specifies the desired boundary conditions:
!    IC(1) = IBEG, desired condition at beginning of data.
!    IC(2) = IEND, desired condition at end of data.
!    See the prologue to PCHIC for details.
!
!    Input, real ( kind = 8 ) VC(2), specifies desired boundary values, as
!    indicated above.  VC(1) need be set only if IC(1) = 2 or 3.
!    VC(2) need be set only if IC(2) = 2 or 3.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) H(N), interval lengths.  H(I) = X(I+1)-X(I),
!    for I = 1 to N-1.
!
!    Input, real ( kind = 8 ) SLOPE(N), the data slopes.
!    SLOPE(I) = ( Y(I+1) - Y(I) ) / H(I), for I = 1 to N-1.
!
!    Input/output, real ( kind = 8 ) D(INCFD,N), the derivative values at the
!    data points.  The value corresponding to X(I) must be stored in
!    D(1+(I-1)*INCFD).  On output, the value of D at X(1) and/or X(N) is
!    changed, if necessary, to produce the requested boundary conditions.
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in D.
!    This argument is provided primarily for 2-d applications.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    1, if IBEG < 0 and D(1) had to be adjusted for monotonicity.
!    2, if IEND < 0 and D(1+(N-1)*INCFD) had to be adjusted for monotonicity.
!    3, if both of the above are true.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) h(n)
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ic(2)
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ierf
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) pchdf
  real ( kind = 8 ) pchst
  real ( kind = 8 ) slope(n)
  real ( kind = 8 ) stemp(3)
  real ( kind = 8 ) vc(2)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtemp(4)

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0
!
!  Set to default boundary conditions if N is too small.
!
  if ( n < abs ( ibeg ) ) then
    ibeg = 0
  end if

  if ( n < abs ( iend ) ) then
    iend = 0
  end if
!
!  Treat beginning boundary condition.
!
  if ( ibeg == 0 ) then
    go to 2000
  end if

  k = abs ( ibeg )

  if ( k == 1 ) then
!
!  Boundary value provided.
!
     d(1,1) = vc(1)
  else if ( k == 2 ) then
!
!  Boundary second derivative provided.
!
     d(1,1) = 0.5D+00 *( (3.0D+00 * slope(1) - d(1,2)) - 0.5D+00 * vc(1) * h(1) )
  else if ( k < 5 ) then
!
!  Use K-point derivative formula.
!  Pick up first K points, in reverse order.
!
     do j = 1, k
       index = k-j+1
       xtemp(j) = x(index)
       if ( j < k ) then
         stemp(j) = slope(index-1)
       end if
     end do

     d(1,1) = pchdf ( k, xtemp, stemp, ierf )

     if ( ierf /= 0 ) then
       ierr = -1
       call xerror ('PCHCE -- error return from pchdf', ierr, 1)
       return
     end if

  else
!
!  Use 'not a knot' condition.
!
     d(1,1) = ( 3.0D+00 * ( h(1) * slope(2) + h(2) * slope(1) ) &
       - 2.0D+00 * ( h(1) + h(2) ) * d(1,2) - h(1) * d(1,3) ) / h(2)
  end if
!
!  Check d(1,1) for compatibility with monotonicity.
!
  if ( ibeg <= 0 ) then

    if ( slope(1) == 0.0D+00 ) then
      if ( d(1,1) /= 0.0D+00 ) then
        d(1,1) = 0.0D+00
        ierr = ierr + 1
      end if
    else if ( pchst ( d(1,1), slope(1) ) < 0.0D+00 ) then
      d(1,1) = 0.0D+00
      ierr = ierr + 1
    else if ( 3.0D+00 * abs ( slope(1) ) < abs ( d(1,1) ) ) then
      d(1,1) = 3.0D+00 * slope(1)
      ierr = ierr + 1
    end if

  end if

2000 continue
!
!  Treat end boundary condition.
!
  if ( iend == 0 ) then
    return
  end if

  k = abs ( iend )
  if ( k == 1 ) then
!
!  Boundary value provided.
!
     d(1,n) = vc(2)

  else if ( k == 2 ) then
!
!  Boundary second derivative provided.
!
     d(1,n) = 0.5D+00 * ( (3.0D+00 * slope(n-1) - d(1,n-1)) &
       + 0.5D+00 * vc(2) * h(n-1) )

  else if ( k < 5 ) then
!
!  Use K-point derivative formula.  Pick up last K points.
!
     do j = 1, k
       index = n - k + j
       xtemp(j) = x(index)
       if ( j < k ) then
         stemp(j) = slope(index)
       end if
     end do

     d(1,n) = pchdf ( k, xtemp, stemp, ierf )

     if ( ierf /= 0 ) then
       ierr = -1
       call xerror ('pchce -- error return from pchdf', ierr, 1)
       return
     end if

  else
!
!  Use 'not a knot' condition.
!
     d(1,n) = ( 3.0D+00 * ( h(n-1) * slope(n-2) + h(n-2) * slope(n-1) ) &
       - 2.0D+00 * ( h(n-1) + h(n-2)) * d(1,n-1) - h(n-1) * d(1,n-2) ) / h(n-2)
  end if

  if ( 0 < iend ) then
    return
  end if
!
!  Check D(1,n) for compatibility with monotonicity.
!
  if ( slope(n-1) == 0.0D+00 ) then
    if ( d(1,n) /= 0.0D+00 ) then
       d(1,n) = 0.0D+00
       ierr = ierr + 2
    end if
  else if ( pchst ( d(1,n), slope(n-1) ) < 0.0D+00 ) then
    d(1,n) = 0.0D+00
    ierr = ierr + 2
  else if ( 3.0D+00 * abs ( slope(n-1) ) < abs ( d(1,n) ) ) then
    d(1,n) = 3.0D+00 * slope(n-1)
    ierr = ierr + 2
  end if

  return
end
subroutine pchci ( n, h, slope, d, incfd )

!*****************************************************************************80
!
!! PCHCI sets derivatives for a monotone piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  If the data are only piecewise monotonic, the
!    interpolant will have an extremum at each point where monotonicity
!    switches direction.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D array.
!
!    The resulting piecewise cubic Hermite function should be identical
!    (within roundoff error) to that produced by PCHIM.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) H(N), interval lengths.  H(I) = X(I+1)-X(I),
!    for I = 1 to N-1.
!
!    Input, real ( kind = 8 ) SLOPE(N), the data slopes.
!    SLOPE(I) = ( Y(I+1) - Y(I) ) / H(I), for I = 1 to N-1.
!
!    Output, real ( kind = 8 ) D(INCFD,N), the derivative values at the data
!    points.  If the data are monotonic, these values will determine a monotone
!    cubic Hermite function.  The value corresponding to X(I) is stored in
!    D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in D.
!    This argument is provided primarily for 2D applications.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) dmax
  real ( kind = 8 ) dmin
  real ( kind = 8 ) drat1
  real ( kind = 8 ) drat2
  real ( kind = 8 ) h(n)
  real ( kind = 8 ) hsum
  real ( kind = 8 ) hsumt3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nless1
  real ( kind = 8 ) pchst
  real ( kind = 8 ) slope(n)
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2

  nless1 = n - 1
  del1 = slope(1)
!
!  Special case N=2 -- use linear interpolation.
!
  if ( nless1 <= 1 ) then
    d(1,1) = del1
    d(1,n) = del1
    return
  end if
!
!  Normal case, 3 <= N.
!
  del2 = slope(2)
!
!  Set D(1) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  hsum = h(1) + h(2)
  w1 = ( h(1) + hsum ) / hsum
  w2 = -h(1) / hsum
  d(1,1) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,1), del1 ) <= 0.0D+00 ) then
    d(1,1) = 0.0D+00
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
!
!  Need do this check only if monotonicity switches.
!
    dmax = 3.0D+00 * del1

    if ( abs ( dmax ) < abs ( d(1,1) ) ) then
      d(1,1) = dmax
    end if

  end if
!
!  Loop through interior points.
!
  do i = 2, nless1

    if ( i /= 2 ) then
      hsum = h(i-1) + h(i)
      del1 = del2
      del2 = slope(i)
    end if
!
!  Set D(I)=0 unless data are strictly monotonic.
!
    d(1,i) = 0.0D+00
!
!  Use Brodlie modification of Butland formula.
!
    if ( 0.0D+00 < pchst ( del1, del2 ) ) then

      hsumt3 = hsum + hsum + hsum
      w1 = ( hsum + h(i-1)) / hsumt3
      w2 = ( hsum + h(i)  ) / hsumt3
      dmax = max ( abs ( del1 ), abs ( del2 ) )
      dmin = min ( abs ( del1 ), abs ( del2 ) )
      drat1 = del1 / dmax
      drat2 = del2 / dmax
      d(1,i) = dmin / ( w1 * drat1 + w2 * drat2 )

    end if

  end do
!
!  Set D(N) via non-centered three point formula, adjusted to
!  be shape preserving.
!
  w1 = -h(n-1) / hsum
  w2 = ( h(n-1) + hsum ) / hsum
  d(1,n) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,n), del2 ) <= 0.0D+00 ) then
    d(1,n) = 0.0D+00
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
    dmax = 3.0D+00 * del2
    if ( abs ( dmax ) < abs ( d(1,n) ) ) then
      d(1,n) = dmax
    end if
  end if

  return
end
subroutine pchcs ( switch, n, h, slope, d, incfd, ierr )

!*****************************************************************************80
!
!! PCHCS adjusts the curve produced by PCHIM so it is more "visually pleasing".
!
!  Discussion:
!
!    PCHCS is called by PCHIC to adjust the values of D in the vicinity of a
!    switch in direction of monotonicity, to produce a more "visually
!    pleasing" curve than that given by PCHIM.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) SWITCH, indicates the amount of control desired
!    over local excursions from data.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 3.
!
!    Input, real ( kind = 8 ) H(N), interval lengths.  H(I) = X(I+1)-X(I),
!    for I = 1 to N-1.
!
!    Input, real ( kind = 8 ) SLOPE(N), the data slopes.
!    SLOPE(I) = ( Y(I+1) - Y(I) ) / H(I), for I = 1 to N-1.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, the derivative values at
!    the data points, as determined by PCHIC.  On output, derivatives in the
!    vicinity of switches in direction of monotonicity may be adjusted to
!    produce a more "visually pleasing" curve.  The value corresponding to
!    X(I) is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in D.
!    This argument is provided primarily for 2D applications.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    negative, trouble in PCHSW.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) del(3)
  real ( kind = 8 ) dext
  real ( kind = 8 ) dfloc
  real ( kind = 8 ) dfmx
  real ( kind = 8 ) fact
  real ( kind = 8 ) h(n)
  real ( kind = 8 ), parameter :: fudge = 4.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nless1
  real ( kind = 8 ) pchst
  real ( kind = 8 ) slmax
  real ( kind = 8 ) slope(n)
  real ( kind = 8 ) switch
  real ( kind = 8 ) wtave(2)
!
!  Define inline function for weighted average of slopes.
!
  real ( kind = 8 ) pchsd, s1, s2, h1, h2
  pchsd ( s1, s2, h1, h2 ) = ( h2 / ( h1 + h2 ) ) * s1 + ( h1 / ( h1 + h2 ) ) * s2
!
!  Initialize.
!
  ierr = 0
  nless1 = n - 1
!
!  Loop over segments.
!
  do i = 2, nless1

     if ( pchst ( slope(i-1), slope(i) ) )  100, 300, 900

  100    continue
!
!  Slope switches monotonicity at i-th point
!
!  Do not change D if 'up-down-up'.
!
        if ( 2 < i ) then
          if ( 0.0D+00 < pchst ( slope(i-2), slope(i) ) ) then
            cycle
          end if
        end if

        if ( i < nless1 ) then
          if ( 0.0D+00 < pchst ( slope(i+1), slope(i-1) ) ) then
            cycle
          end if
        end if
!
!  Compute provisional value for D(1,i).
!
        dext = pchsd ( slope(i-1), slope(i), h(i-1), h(i) )
!
!  Determine which interval contains the extremum.
!
        if ( pchst ( dext, slope(i-1) ) )  200, 900, 250

  200       continue
!
!  DEXT and slope(i-1) have opposite signs.
!  extremum is in (x(i-1),x(i)).
!
           k = i - 1
!
!  Set up to compute new values for D(1,i-1) and D(1,i).
!
           wtave(2) = dext
           if ( 1 < k ) then
             wtave(1) = pchsd (slope(k-1), slope(k), h(k-1), h(k))
           end if
           go to 400

  250       continue
!
!  DEXT and SLOPE(I) have opposite signs.
!  The extremum is in (x(i),x(i+1)).
!
           k = i
!
!  Set up to compute new values for D(1,i) and D(1,i+1).
!
           wtave(1) = dext
           if ( k < nless1 ) then
             wtave(2) = pchsd ( slope(k), slope(k+1), h(k), h(k+1) )
           end if
           go to 400

  300    continue
!
!  At least one of SLOPE(I-1) and slope(i) is zero.
!  Check for flat-topped peak
!
        if ( i == nless1 ) then
          cycle
        end if

        if ( 0.0D+00 <= pchst ( slope(i-1), slope(i+1) ) ) then
          cycle
        end if
!
!  We have flat-topped peak on (x(i),x(i+1)).
!
        k = i
!
!  Set up to compute new values for d(1,i) and d(1,i+1).
!
        wtave(1) = pchsd ( slope(k-1), slope(k), h(k-1), h(k) )
        wtave(2) = pchsd ( slope(k), slope(k+1), h(k), h(k+1) )

  400    continue
!
!  At this point we have determined that there will be an extremum
!  on (x(k),x(k+1)), where k=i or i-1, and have set array WTAVE.
!  wtave(1) is a weighted average of slope(k-1) and slope(k), if 1 < K
!  wtave(2) is a weighted average of slope(k) and slope(k+1), if k<n-1
!
     slmax = abs ( slope(k) )
     if ( 1 < k ) then
       slmax = max ( slmax, abs ( slope(k-1) ) )
     end if

     if ( k < nless1 ) then
       slmax = max ( slmax, abs ( slope(k+1) ) )
     end if

     if ( 1 < k ) then
       del(1) = slope(k-1) / slmax
     end if

     del(2) = slope(k) / slmax

     if ( k < nless1 ) then
       del(3) = slope(k+1) / slmax
     end if

     if ( 1 < k .and. k < nless1 ) then
!
!  Normal case -- extremum is not in a boundary interval.
!
        fact = fudge * abs ( del(3) * ( del(1) - del(2) ) &
          * ( wtave(2) / slmax ) )

        d(1,k) = d(1,k) + min ( fact, 1.0D+00 ) * ( wtave(1) - d(1,k) )
        fact = fudge * abs ( del(1) * ( del(3) - del(2) ) &
          * ( wtave(1) / slmax ) )
        d(1,k+1) = d(1,k+1) + min ( fact, 1.0D+00 ) * ( wtave(2) - d(1,k+1) )
     else
!
!  Special case K=1 (which can occur only if I=2) or
!  k=nless1 (which can occur only if i=nless1).
!
        fact = fudge * abs ( del(2) )
        d(1,i) = min ( fact, 1.0D+00 ) * wtave(i-k+1)
!
!  Note that i-k+1 = 1 if k=i  (=nless1),
!            i-k+1 = 2 if k=i-1(=1).
!
     end if
!
!  Adjust if necessary to limit excursions from data.
!
     if ( switch <= 0.0D+00 ) then
       cycle
     end if

     dfloc = h(k) * abs ( slope(k) )

     if ( 1 < k ) then
       dfloc = max ( dfloc, h(k-1) * abs ( slope(k-1) ) )
     end if

     if ( k < nless1 ) then
       dfloc = max ( dfloc, h(k+1) * abs ( slope(k+1) ) )
     end if

     dfmx = switch * dfloc
     indx = i-k+1
!
!  INDX = 1 if K = I,
!  INDX = 2 if K = I-1.
!
     call pchsw ( dfmx, indx, d(1,k), d(1,k+1), h(k), slope(k), ierr )

     if ( ierr /= 0 ) then
       return
     end if

  900 continue

  end do

  return
end
function pchdf ( k, x, s, ierr )

!*****************************************************************************80
!
!! PCHDF approximates a derivative using divided differences of data.
!
!  Discussion:
!
!    The routine uses a divided difference formulation to compute a K-point
!    approximation to the derivative at X(K) based on the data in X and S.
!
!    It is called by PCHCE and PCHSP to compute 3 and 4 point boundary
!    derivative approximations.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Pages 10-16,
!    Springer-Verlag, 1978.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, is the order of the desired derivative approximation.
!    K must be at least 3.
!
!    Input, real ( kind = 8 ) X(K), contains the K values of the independent
!    variable.  X need not be ordered, but the values must be distinct.
!
!    Input/output, real ( kind = 8 ) S(K-1).  On input, the associated slope
!    values:
!      S(I) = ( F(I+1)-F(I))/(X(I+1)-X(I))
!    On output, S is overwritten.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no error.
!    -1, if K < 2.
!
!    Output, real ( kind = 8 ) PCHDF, the desired derivative approximation if
!    IERR=0 or to zero if IERR=-1.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  real ( kind = 8 ) pchdf
  real ( kind = 8 ) s(k-1)
  real ( kind = 8 ) value
  real ( kind = 8 ) x(k)
!
!  Check for legal value of K.
!
  if ( k < 3 ) then
    ierr = -1
    call xerror ( 'pchdf -- k less than three', ierr, 1 )
    pchdf = 0.0D+00
    return
  end if
!
!  Compute coefficients of interpolating polynomial.
!
  do j = 2, k-1
    do i = 1, k-j
      s(i) = ( s(i+1) - s(i) ) / ( x(i+j) - x(i) )
    end do
  end do
!
!  Evaluate the derivative at X(K).
!
  value = s(1)

  do i = 2, k-1
    value = s(i) + value * ( x(k) - x(i) )
  end do

  ierr = 0
  pchdf = value

  return
end
subroutine pchev ( n, x, f, d, nval, xval, fval, dval, ierr )

!*****************************************************************************80
!
!! PCHEV evaluates a piecewise cubic Hermite or spline function.
!
!  Discussion:
!
!    PCHEV evaluates the function and first derivative of a piecewise
!    cubic Hermite or spline function at an array of points XVAL.
!
!    The evaluation will be most efficient if the elements of XVAL are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XVAL(J)
!    implies
!      X(I) <= XVAL(K).
!
!    If any of the XVAL are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(N), the derivative values.  D(i) is the value
!    corresponding to X(I).
!
!    Input, integer ( kind = 4 ) NVAL, the number of points at which the functions are
!    to be evaluated.
!
!    Input, real ( kind = 8 ) XVAL(NVAL), the points at which the functions
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) FVAL(NVAL), the values of the cubic Hermite
!    function at XVAL.
!
!    Output, real ( kind = 8 ) DVAL(NVAL), the derivatives of the cubic
!    Hermite function at XVAL.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -3, if the X array is not strictly increasing.
!    -4, if NVAL < 1.
!    -5, if an error has occurred in CHFDV.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nval

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dval(nval)
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fval(nval)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), save :: incfd = 1
  logical, save :: skip = .true.
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval(nval)

  call pchfd ( n, x, f, d, incfd, skip, nval, xval, fval, dval, ierr )

  return
end
subroutine pchez ( n, x, f, d, spline, wk, lwk, ierr )

!*****************************************************************************80
!
!! PCHEZ carries out easy to use spline or cubic Hermite interpolation.
!
!  Discussion:
!
!    This routine sets derivatives for spline (two continuous derivatives)
!    or Hermite cubic (one continuous derivative) interpolation.
!    Spline interpolation is smoother, but may not "look" right if the
!    data contains both "steep" and "flat" sections.  Hermite cubics
!    can produce a "visually pleasing" and monotone interpolant to
!    monotone data.
!
!    This routine is an easy to use driver for the PCHIP routines.
!    Various boundary conditions are set to default values by PCHEZ.
!    Many other choices are available in the subroutines PCHIC,
!    PCHIM and PCHSP.
!
!    Use PCHEV to evaluate the resulting function and its derivative.
!
!    If SPLINE is TRUE, the interpolating spline satisfies the default
!    "not-a-knot" boundary condition, with a continuous third derivative
!    at X(2) and X(N-1).
!
!    If SPLINE is FALSE, the interpolating Hermite cubic will be monotone
!    if the input data is monotone.  Boundary conditions are computed from
!    the derivative of a local quadratic unless this alters monotonicity.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!    Carl deBoor,
!    A Practical Guide to Splines, Chapter IV,
!    Springer-Verlag,
!    New York, 1978.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Output, real ( kind = 8 ) D(N), the derivative values at the data points.
!
!    Input, logical SPLINE, specifies if the interpolant is to be a spline
!    with two continuous derivatives (SPLINE is TRUE), or a Hermite cubic
!    interpolant with one continuous derivative (SPLINE is FALSE).
!
!    Workspace, real ( kind = 8 ) WK(LWK), required only if SPLINE is TRUE.
!
!    Input, integer ( kind = 4 ) LWK, the length of the work array WK, which must
!    be at least 2*N.  However, WK is not needed if SPLINE is FALSE,
!    and in this case LWK is arbitrary.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, can only occur when SPLINE is FALSE,  means that there were
!      IERR switches in the direction of monotonicity.  When SPLINE is
!      FALSE, PCHEZ guarantees that if the input data is monotone, the
!      interpolant will be too.  This warning is to alert you to the fact
!      that the input data was not monotone.
!    -1, if N < 2.
!    -3, if the X array is not strictly increasing.
!    -7, if LWK is less than 2*N and SPLINE is TRUE.
!
  implicit none

  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ), save, dimension ( 2 ) :: ic = (/ 0, 0 /)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), parameter :: incfd = 1
  logical spline
  real ( kind = 8 ) vc(2)
  real ( kind = 8 ) wk(lwk)
  real ( kind = 8 ) x(n)

  if ( spline ) then
    call pchsp ( ic, vc, n, x, f, d, incfd, wk, lwk, ierr )
  else
    call pchim ( n, x, f, d, incfd, ierr )
  end if

  return
end
subroutine pchfd ( n, x, f, d, incfd, skip, ne, xe, fe, de, ierr )

!*****************************************************************************80
!
!! PCHFD evaluates a piecewise cubic Hermite function.
!
!  Discsussion:
!
!    PCHFD evaluates a piecewise cubic Hermite function and its first
!    derivative at an array of points.  PCHFD may be used by itself
!    for Hermite interpolation, or as an evaluator for PCHIM
!    or PCHIC.
!
!    PCHFD evaluates the cubic Hermite function and its first derivative
!    at the points XE.
!
!    If only function values are required, use PCHFE instead.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!    Most of the coding between the call to CHFDV and the end of
!    the IR loop could be eliminated if it were permissible to
!    assume that XE is ordered relative to X.
!
!    CHFDV does not assume that X1 is less than X2.  Thus, it would
!    be possible to write a version of PCHFD that assumes a strictly
!    decreasing X array by simply running the IR loop backwards
!    and reversing the order of appropriate tests.
!
!    The present code has a minor bug, which I have decided is not
!    worth the effort that would be required to fix it.
!    If XE contains points in [X(N-1),X(N)], followed by points less than
!    X(N-1), followed by points greater than X(N), the extrapolation points
!    will be counted (at least) twice in the total returned in IERR.
!
!    The evaluation will be most efficient if the elements of XE are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XE(J)
!    implies
!      X(I) <= XE(K).
!
!    If any of the XE are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Modified:
!
!    13 August 2005
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with 0 <= IERR, SKIP will be set to TRUE.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), points at which the function is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) FE(NE), the values of the cubic Hermite
!    function at XE.
!
!    Output, real ( kind = 8 ) DE(NE), the derivative of the cubic
!    Hermite function at XE.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if NE < 1.
!    -5, if an error has occurred in the lower-level routine CHFDV.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) de(ne)
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) fe(ne)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_first
  integer ( kind = 4 ) j_new
  integer ( kind = 4 ) j_save
  integer ( kind = 4 ) next(2)
  integer ( kind = 4 ) nj
  logical skip
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xe(ne)
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFD - Fatal error!'
      write ( *, '(a)' ) '  Number of data points less than 2.'
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFD - Fatal error!'
      write ( *, '(a)' ) '  Increment less than 1.'
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFD - Fatal error!'
        write ( *, '(a)' ) '  X array not strictly increasing.'
        return
      end if
    end do

  end if

  if ( ne < 1 ) then
    ierr = -4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHFD - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points less than 1.'
    return
  end if

  ierr = 0
  skip = .true.
!
!  Loop over intervals.
!  The interval index is IL = IR - 1.
!  The interval is X(IL) <= X < X(IR).
!
  j_first = 1
  ir = 2

  do
!
!  Skip out of loop if have processed all evaluation points.
!
    if ( ne < j_first ) then
      exit
    end if
!
!  Locate all points in interval.
!
    j_save = ne + 1

    do j = j_first, ne
      if ( x(ir) <= xe(j) ) then
        j_save = j
        if ( ir == n ) then
          j_save = ne + 1
        end if
        exit
      end if
    end do
!
!  Have located first point beyond interval.
!
    j = j_save

    nj = j - j_first
!
!  Skip evaluation if no points in interval.
!
    if ( nj /= 0 ) then
!
!  Evaluate cubic at XE(J_FIRST:J-1).
!
      call chfdv ( x(ir-1), x(ir), f(1,ir-1), f(1,ir), d(1,ir-1), d(1,ir), &
        nj, xe(j_first), fe(j_first), de(j_first), next, ierc )

      if ( ierc < 0 ) then
        ierr = -5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFD - Fatal error!'
        write ( *, '(a)' ) '  Error return from CHFDV.'
        return
      end if
!
!  In the current set of XE points, there are NEXT(2) to the right of X(IR).
!
      if ( next(2) /= 0 ) then

        if ( ir < n ) then
          ierr = -5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PCHFD - Fatal error!'
          write ( *, '(a)' ) '  IR < N.'
          return
        end if
!
!  These are actually extrapolation points.
!
        ierr = ierr + next(2)

      end if
!
!  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
!
      if ( next(1) /= 0 ) then
!
!  These are actually extrapolation points.
!
        if ( ir <= 2 ) then
          ierr = ierr + next(1)
!
!  XE is not ordered relative to X, so must adjust evaluation interval.
!
!  First, locate first point to left of X(IR-1).
!
        else

          j_new = -1

          do i = j_first, j-1
            if ( xe(i) < x(ir-1) ) then
              j_new = i
              exit
            end if
          end do

          if ( j_new == -1 ) then
            ierr = -5
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PCHFD - Fatal error!'
            write ( *, '(a)' ) '  Could not bracket the data point.'
            return
          end if
!
!  Reset J.  This will be the new J_FIRST.
!
          j = j_new
!
!  Now find out how far to back up in the X array.
!
          do i = 1, ir-1
            if ( xe(j) < x(i) ) then
              exit
            end if
          end do
!
!  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
!
!  Reset IR, recognizing that it will be incremented before cycling.
!
          ir = max ( 1, i-1 )

        end if

      end if

      j_first = j

    end if

    ir = ir + 1

    if ( n < ir ) then
      exit
    end if

  end do

  return
end
subroutine pchfe ( n, x, f, d, incfd, skip, ne, xe, fe, ierr )

!*****************************************************************************80
!
!! PCHFE evaluates a piecewise cubic Hermite function at an array of points.
!
!  Description:
!
!    PCHFE may be used by itself for Hermite interpolation, or as an
!    evaluator for PCHIM or PCHIC.
!
!    PCHFE evaluates the cubic Hermite function at the points XE.
!
!    To provide compatibility with PCHIM and PCHIC, the routine includes an
!    increment between successive values of the F and D arrays.
!
!    Most of the coding between the call to CHFEV and the end of
!    the IR loop could be eliminated if it were permissible to
!    assume that XE is ordered relative to X.
!
!    CHFEV does not assume that X1 is less than X2.  Thus, it would
!    be possible to write a version of pchfe that assumes a strictly
!    decreasing X array by simply running the IR loop backwards
!    and reversing the order of appropriate tests.
!
!    The present code has a minor bug, which I have decided is not
!    worth the effort that would be required to fix it.
!    If XE contains points in [X(N-1),X(N)], followed by points less than
!    X(N-1), followed by points greater than X(N), the extrapolation points
!    will be counted (at least) twice in the total returned in IERR.
!
!    The evaluation will be most efficient if the elements of XE are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XE(J)
!    implies
!      X(I) <= XE(K).
!
!    If any of the XE are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Modified:
!
!    12 August 2005
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with 0 <= IERR, SKIP will be set to TRUE.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), points at which the function is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) FE(NE), the values of the cubic Hermite
!    function at XE.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if NE < 1.
!    -5, error in CHFEV.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) fe(ne)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_first
  integer ( kind = 4 ) j_new
  integer ( kind = 4 ) j_save
  integer ( kind = 4 ) next(2)
  integer ( kind = 4 ) nj
  logical skip
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xe(ne)
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFE - Fatal error!'
      write ( *, '(a)' ) '  Number of data points less than 2.'
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFE - Fatal error!'
      write ( *, '(a)' ) '  Increment less than 1.'
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFE - Fatal error!'
        write ( *, '(a)' ) '  X array not strictly increasing.'
        return
      end if
    end do

  end if

  if ( ne < 1 ) then
    ierr = -4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHFE - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points less than 1.'
    return
  end if

  ierr = 0
  skip = .true.
!
!  Loop over intervals.
!  The interval index is IL = IR-1.
!  The interval is X(IL) <= X < X(IR).
!
  j_first = 1
  ir = 2

  do
!
!  Skip out of the loop if have processed all evaluation points.
!
    if ( ne < j_first ) then
      exit
    end if
!
!  Locate all points in the interval.
!
    j_save = ne + 1

    do j = j_first, ne
      if ( x(ir) <= xe(j) ) then
        j_save = j
        if ( ir == n ) then
          j_save = ne + 1
        end if
        exit
      end if
    end do
!
!  Have located first point beyond interval.
!
    j = j_save

    nj = j - j_first
!
!  Skip evaluation if no points in interval.
!
    if ( nj /= 0 ) then
!
!  Evaluate cubic at XE(J_FIRST:J-1).
!
      call chfev ( x(ir-1), x(ir), f(1,ir-1), f(1,ir), d(1,ir-1), d(1,ir), &
        nj, xe(j_first), fe(j_first), next, ierc )

      if ( ierc < 0 ) then
        ierr = -5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFE - Fatal error!'
        write ( *, '(a)' ) '  Error return from CHFEV.'
        return
      end if
!
!  In the current set of XE points, there are NEXT(2) to the right of X(IR).
!
      if ( next(2) /= 0 ) then

        if ( ir < n ) then
          ierr = -5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PCHFE - Fatal error!'
          write ( *, '(a)' ) '  IR < N.'
          return
        end if
!
!  These are actually extrapolation points.
!
        ierr = ierr + next(2)

      end if
!
!  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
!
      if ( next(1) /= 0 ) then
!
!  These are actually extrapolation points.
!
        if ( ir <= 2 ) then
          ierr = ierr + next(1)
        else

          j_new = -1

          do i = j_first, j-1
            if ( xe(i) < x(ir-1) ) then
              j_new = i
              exit
            end if
          end do

          if ( j_new == -1 ) then
            ierr = -5
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PCHFE - Fatal error!'
            write ( *, '(a)' ) '  Could not bracket the data point.'
            return
          end if
!
!  Reset J.  This will be the new J_FIRST.
!
          j = j_new
!
!  Now find out how far to back up in the X array.
!
          do i = 1, ir-1
            if ( xe(j) < x(i) ) then
              exit
            end if
          end do
!
!  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
!
!  Reset IR, recognizing that it will be incremented before cycling.
!
          ir = max ( 1, i-1 )

        end if

      end if

      j_first = j

    end if

    ir = ir + 1

    if ( n < ir ) then
      exit
    end if

  end do

  return
end
function pchia ( n, x, f, d, incfd, skip, a, b, ierr )

!*****************************************************************************80
!
!! PCHIA evaluates the integral of a piecewise cubic Hermite function.
!
!  Description:
!
!    PCHIA evaluates the definite integral of a cubic Hermite function
!    over the interval [A, B].
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VALUE, the value of the requested integral.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with 0 <= IERR, SKIP will be set to TRUE.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.  The
!    integration interval is normally contained within [X(1),X(N)], but
!    this is not required.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    1, if A is outside the interval [X(1),X(N)].
!    2, if B is outside the interval [X(1),X(N)].
!    3, if both of the above are true.  This means that either [A,B] contains
!       the data interval or the intervals do not intersect at all.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) chfiv
  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ierd
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ierv
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ir
  real ( kind = 8 ) pchia
  real ( kind = 8 ) pchid
  logical skip
  real ( kind = 8 ) value
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      call xerror ('pchia -- number of data points less than two', ierr, 1)
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      call xerror ('pchia -- increment less than one', ierr, 1)
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        call xerror ('pchia -- x-array not strictly increasing', ierr, 1)
        return
      end if
    end do

    skip = .true.

  end if

  ierr = 0

  if ( a < x(1) .or. x(n) < a ) then
    ierr = ierr + 1
  end if

  if ( b < x(1) .or. x(n) < b ) then
    ierr = ierr + 2
  end if
!
!  Compute integral value.
!
  if ( a == b ) then

    value = 0.0D+00

  else

    xa = min (a, b)
    xb = max (a, b)
!
!  Interval is to left of X(2), so use first cubic.
!
    if ( xb <= x(2) ) then

      value = chfiv ( x(1), x(2), f(1,1), f(1,2), &
        d(1,1), d(1,2), a, b, ierv )

      if ( ierv < 0 ) then
        ierr = -4
        call xerror ('pchia -- trouble in chfiv', ierr, 1)
        return
      end if
!
!  Interval is to right of x(n-1), so use last cubic.
!
    else if ( x(n-1) <= xa ) then

      value = chfiv ( x(n-1), x(n), f(1,n-1), f(1,n), &
        d(1,n-1), d(1,n), a, b, ierv )

      if ( ierv < 0 ) then
        ierr = -4
        call xerror ('pchia -- trouble in chfiv', ierr, 1)
        return
      end if
!
!  Normal case -- xa<xb, xa<x(n-1), x(2) < xb.
!
!  Locate ia and ib such that
!  x(ia-1) < xa <= x(ia) <= x(ib) <= xb <= x(ib+1)
!
    else

      ia = 1
      do i = 1, n-1
        if ( x(i) < xa ) then
          ia = i + 1
        end if
      end do
!
!  IA = 1 implies xa<x(1) .  Otherwise,
!  ia is largest index such that x(ia-1)<xa,.
!
      ib = n
      do i = n, ia, -1
        if ( xb < x(i) ) then
          ib = i - 1
        end if
      end do
!
!  IB = N implies X(N) < XB.  Otherwise,
!  ib is smallest index such that xb<x(ib+1) .
!
!  Compute the integral.
!
      ierv = 0
      if ( ib < ia ) then
!
!  This means IB = IA-1 and (A,B) is a subset of (x(ib),x(ia)).
!
           value = chfiv ( x(ib), x(ia), f(1,ib), f(1,ia), &
             d(1,ib), d(1,ia), a, b, ierv )

           if ( ierv < 0 ) then
             ierr = -4
             call xerror ('pchia -- trouble in chfiv', ierr, 1)
             return
           end if

        else
!
!  First compute integral over (x(ia),x(ib)).
!
           if ( ib == ia ) then
              value = 0.0D+00
           else

              value = pchid ( n, x, f, d, incfd, skip, ia, ib, ierd )

              if ( ierd < 0 ) then
                ierr = -5
                call xerror ('pchia -- trouble in pchid', ierr, 1)
                return
              end if

           end if
!
!  Then add on integral over ( XA, X(IA) ).
!
           if ( xa < x(ia) ) then
              il = max ( 1, ia-1 )
              ir = il + 1

              value = value + chfiv ( x(il), x(ir), f(1,il), f(1,ir), &
                d(1,il), d(1,ir), xa, x(ia), ierv )

              if ( ierv < 0 ) then
                ierr = -4
                call xerror ('pchia -- trouble in chfiv', ierr, 1)
                return
              end if

           end if
!
!  Then add on integral over ( X(IB), XB ).
!
           if ( x(ib) < xb ) then
              ir = min ( ib+1, n )
              il = ir - 1

              value = value + chfiv ( x(il), x(ir), f(1,il), f(1,ir), &
                d(1,il), d(1,ir), x(ib), xb, ierv )

              if ( ierv < 0 ) then
                ierr = -4
                call xerror ('pchia -- trouble in chfiv', ierr, 1)
                return
              end if

           end if
!
!  Adjust sign if necessary.
!
           if ( b < a ) then
             value = -value
           end if

        end if
     end if
  end if

  pchia = value

  return
end
subroutine pchic ( ic, vc, switch, n, x, f, d, incfd, wk, nwk, ierr )

!*****************************************************************************80
!
!! PCHIC sets derivatives for a piecewise monotone cubic Hermite interpolant.
!
!  Description:
!
!    PCHIC sets derivatives needed to determine a piecewise monotone
!    piecewise cubic Hermite interpolant to the data given in X and F
!    satisfying the boundary conditions specified by IC and VC.
!
!    The treatment of points where monotonicity switches direction is
!    controlled by argument switch.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!    User control is available over boundary conditions and
!    treatment of points where monotonicity switches direction.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IC(2), specifies desired boundary conditions:
!    IC(1) = IBEG, desired condition at beginning of data.
!    IC(2) = IEND, desired condition at end of data.
!
!    IBEG = 0  for the default boundary condition (the same as used by PCHIM).
!    If IBEG/=0, then its sign indicates whether the boundary derivative is
!    to be adjusted, if necessary, to be compatible with monotonicity:
!    0 < IBEG, if no adjustment is to be performed.
!    IBEG < 0, if the derivative is to be adjusted for monotonicity.
!
!    Allowable values for the magnitude of IBEG are:
!    1, if first derivative at x(1) is given in VC(1).
!    2, if second derivative at x(1) is given in VC(1).
!    3, to use the 3-point difference formula for D(1).
!       This reverts to the default boundary condition if N < 3.
!    4, to use the 4-point difference formula for D(1).
!       This reverts to the default boundary condition if N < 4.
!    5, to set D(1) so that the second derivative is continuous at X(2).
!       This reverts to the default boundary condition if N < 4.
!       This option is somewhat analogous to the "not a knot"
!       boundary condition provided by PCHSP.
!
!    An error return is taken if 5 < |IBEG|.
!    Only in case IBEG <= 0 is it guaranteed that the interpolant will be
!    monotonic in the first interval.  If the returned value of D(1) lies
!    between zero and 3 * SLOPE(1), the interpolant will be monotonic.  This
!    is not checked if 0 < IBEG.
!    If IBEG < 0 and D(1) had to be changed to achieve monotonicity, a
!    warning error is returned.
!
!    IEND may take on the same values as IBEG, but applied to the derivative
!    at X(N).  In case IEND = 1 or 2, the value is given in VC(2).
!
!    An error return is taken if 5 < |IEND|.
!    Only in case IEND <= 0  is it guaranteed that the interpolant will be
!    monotonic in the last interval.  If the returned value of
!    D(1+(N-1)*INCFD) lies between zero and 3 * SLOPE(N-1), the interpolant
!    will be monotonic.  This is not checked if 0 < IEND.
!    If IEND < 0 and D(1+(N-1)*INCFD) had to be changed to achieve
!    monotonicity, a warning error is returned.
!
!    Input, real ( kind = 8 ) VC(2), specifies desired boundary values,
!    as indicated above.
!    VC(1) need be set only if IC(1) = 1 or 2.
!    VC(2) need be set only if IC(2) = 1 or 2.
!
!    Input, integer ( kind = 4 ) SWITCH, indicates the desired treatment of points where
!    the direction of monotonicity switches:
!    * set SWITCH to zero if the interpolant is required to be monotonic in
!      each interval, regardless of monotonicity of data.  This will cause D
!      to be set to zero at all switch points, thus forcing extrema there.
!      The result of using this option with the default boundary conditions
!      will be identical to using PCHIM, but will generally cost more
!      compute time.  This option is provided only to facilitate comparison
!      of different switch and/or boundary conditions.
!    * set SWITCH nonzero to use a formula based on the 3-point difference
!      formula in the vicinity of switch points.  If SWITCH is positive, the
!      interpolant on each interval containing an extremum is controlled to
!      not deviate from the data by more than SWITCH * DFLOC, where DFLOC is the
!      maximum of the change of F on this interval and its two immediate
!      neighbors.  If SWITCH is negative, no such control is to be imposed.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the dependent values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Output, real ( kind = 8 ) D(INCFD,N), the derivative values at the data
!    points.  These values will determine a monotone cubic Hermite function
!    on each subinterval on which the data are monotonic, except possibly
!    adjacent to switches in monotonicity.  The value corresponding to X(I)
!    is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in F and D.
!
!    Workspace, real ( kind = 8 ) WK(NWK).  The user may wish to know that
!    the returned values, for I = 1 to N-1, are:
!      WK(I)     = H(I)     = X(I+1) - X(I)
!      WK(N-1+I) = SLOPE(I) = ( F(1,I+1) - F(1,I)) / H(I)
!
!    Input, integer ( kind = 4 ) NWK, the length of the work array, which must be at
!    least 2 * ( N - 1 ).
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    1, if IBEG < 0 and D(1) had to be adjusted for monotonicity.
!    2, if IEND < 0 and D(1+(N-1)*INCFD) had to be adjusted for monotonicity.
!    3, if both of the above are true.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if 5 < abs ( IBEG ).
!    -5, if 5 < abs ( IEND ).
!    -6, if both of the above are true.
!    -7, if NWK < 2 * ( N - 1 ).
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nwk

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ic(2)
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) nless1
  real ( kind = 8 ) switch
  real ( kind = 8 ) vc(2)
  real ( kind = 8 ) wk(nwk)
  real ( kind = 8 ) x(n)

  if ( n < 2 ) then
    ierr = -1
    call xerror ('pchic -- number of data points less than two', ierr, 1)
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    call xerror ('pchic -- increment less than one', ierr, 1)
    return
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      call xerror ('pchic -- x-array not strictly increasing', ierr, 1)
      return
    end if
  end do

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0

  if ( 5 < abs ( ibeg ) ) then
    ierr = ierr - 1
  end if

  if ( 5 < abs ( iend ) ) then
    ierr = ierr - 2
  end if

  if ( ierr < 0 ) then
    ierr = ierr - 3
    call xerror ('pchic -- ic out of range', ierr, 1)
    return
  end if
!
!  Function definition is ok -- go on.
!
  nless1 = n - 1
  if ( nwk < 2 * nless1 ) then
    ierr = -7
    call xerror ('pchic -- work array too small', ierr, 1)
    return
  end if
!
!  Set up H and slope arrays.
!
  do i = 1, nless1
    wk(i) = x(i+1) - x(i)
    wk(nless1+i) = (f(1,i+1) - f(1,i)) / wk(i)
  end do
!
!  Special case n=2 -- use linear interpolation.
!
  if ( 1 < nless1 ) then
    go to 1000
  end if

  d(1,1) = wk(2)
  d(1,n) = wk(2)
  go to 3000
!
!  Normal case  (3 <= N) .
!
 1000 continue
!
!  Set interior derivatives and default end conditions.
!
  call pchci ( n, wk(1), wk(n), d, incfd )
!
!  Set derivatives at points where monotonicity switches direction.
!
  if ( switch /= 0.0D+00 ) then

    call pchcs (switch, n, wk(1), wk(n), d, incfd, ierr)

    if ( ierr /= 0 ) then
      ierr = -8
      call xerror ('pchic -- error return from pchcs', ierr, 1)
      return
    end if

  end if
!
!  Set end conditions.
!
 3000 continue

  if ( ibeg == 0 .and. iend == 0 ) then
    return
  end if

  call pchce ( ic, vc, n, x, wk(1), wk(n), d, incfd, ierr )

  if ( ierr < 0 ) then
    ierr = -9
    call xerror ('pchic -- error return from pchce', ierr, 1)
    return
  end if

  return
end
function pchid ( n, x, f, d, incfd, skip, ia, ib, ierr )

!*****************************************************************************80
!
!! PCHID evaluates the definite integral of a piecewise cubic Hermite function.
!
!  Description:
!
!    PCHID evaluates the definite integral of a cubic Hermite function
!    over the interval [X(IA), X(IB)].  The endpoints of the integration
!    interval must be data points.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) F(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in F and D.
!
!    Input/output, logical SKIP, should be set to TRUE if the user wishes to
!    skip checks for validity of preceding parameters, or to FALSE otherwise.
!    This will save time in case these checks have already been performed
!    say, in PCHIM or PCHIC.  SKIP will be set to TRUE on return with
!    IERR = 0 or -4.
!
!    Input, integer ( kind = 4 ) IA, IB, the indices in the X array for the limits of
!    integration.  Both must be in the range [1,N].
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if IA or IB is out of range.
!
!    Output, real ( kind = 8 ) PCHID, the value of the requested integral.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 )  n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) low
  real ( kind = 8 ) pchid
  logical skip
  real ( kind = 8 ) sum2
  real ( kind = 8 ) value
  real ( kind = 8 ) x(n)

  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      call xerror ('pchid -- number of data points less than two', ierr, 1)
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      call xerror ('pchid -- increment less than one', ierr, 1)
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        call xerror ('pchid -- x-array not strictly increasing', ierr, 1)
        return
      end if
    end do

  end if

  skip = .true.

  if ( ia < 1 .or. n < ia ) then
    go to 5004
  end if

  if ( ib < 1 .or. n < ib ) then
    go to 5004
  end if

  ierr = 0
!
!  Compute integral value.
!
  if ( ia == ib ) then

    value = 0.0D+00

  else

    low = min ( ia, ib )
    iup = max ( ia, ib ) - 1
    sum2 = 0.0D+00

    do i = low, iup
      h = x(i+1) - x(i)
      sum2 = sum2 + h * &
        ( ( f(1,i) + f(1,i+1) ) + ( d(1,i) - d(1,i+1) ) * ( h / 6.0D+00 ) )
    end do

    value = 0.5D+00 * sum2

    if ( ib < ia ) then
      value = -value
    end if

  end if

  pchid = value

  return
!
!  error returns.
!
 5004 continue
!
!  ia or ib out of range return.
!
  ierr = -4
  call xerror ('pchid -- ia or ib out of range', ierr, 1)
  return
end
subroutine pchim ( n, x, f, d, incfd, ierr )

!*****************************************************************************80
!
!! PCHIM sets derivatives for a piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    The routine set derivatives needed to determine a monotone piecewise
!    cubic Hermite interpolant to given data.
!
!    The interpolant will have an extremum at each point where
!    monotonicity switches direction.  See PCHIC if user control is desired
!    over boundary or switch conditions.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), dependent variable values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!    PCHIM is designed for monotonic data, but it will work for any F-array.
!    It will force extrema at points where monotonicity switches direction.
!    If some other treatment of switch points is desired, PCHIC should be
!    used instead.
!
!    Output, real ( kind = 8 ) D(INCFD,N), the derivative values at the
!    data points.  If the data are monotonic, these values will determine
!    a monotone cubic Hermite function.  The value corresponding to X(I)
!    is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, the increment between successive
!    values in F and D.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means IERR switches in the direction of monotonicity
!    were detected.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) dmax
  real ( kind = 8 ) dmin
  real ( kind = 8 ) drat1
  real ( kind = 8 ) drat2
  real ( kind = 8 ) dsave
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) hsum
  real ( kind = 8 ) hsumt3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) nless1
  real ( kind = 8 ) pchst
  real ( kind = 8 ) temp
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x(n)
!
!  Check the arguments.
!
  if ( n < 2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHIM - Fatal error!'
    write ( *, '(a)' ) '  Number of data points less than 2.'
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHIM - Fatal error!'
    write ( *, '(a)' ) '  Increment less than 1.'
    return
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHIM - Fatal error!'
      write ( *, '(a)' ) '  X array not strictly increasing.'
      return
    end if
  end do

  ierr = 0
  nless1 = n - 1
  h1 = x(2) - x(1)
  del1 = ( f(1,2) - f(1,1) ) / h1
  dsave = del1
!
!  Special case N=2, use linear interpolation.
!
  if ( n == 2 ) then
    d(1,1) = del1
    d(1,n) = del1
    return
  end if
!
!  Normal case, 3 <= N.
!
  h2 = x(3) - x(2)
  del2 = ( f(1,3) - f(1,2) ) / h2
!
!  Set D(1) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  hsum = h1 + h2
  w1 = ( h1 + hsum ) / hsum
  w2 = -h1 / hsum
  d(1,1) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,1), del1 ) <= 0.0D+00 ) then

    d(1,1) = 0.0D+00
!
!  Need do this check only if monotonicity switches.
!
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then

     dmax = 3.0D+00 * del1

     if ( abs ( dmax ) < abs ( d(1,1) ) ) then
       d(1,1) = dmax
     end if

  end if
!
!  Loop through interior points.
!
  do i = 2, nless1

    if ( 2 < i ) then
      h1 = h2
      h2 = x(i+1) - x(i)
      hsum = h1 + h2
      del1 = del2
      del2 = ( f(1,i+1) - f(1,i) ) / h2
    end if
!
!  Set D(I)=0 unless data are strictly monotonic.
!
    d(1,i) = 0.0D+00

    temp = pchst ( del1, del2 )

    if ( temp < 0.0D+00 ) then

      ierr = ierr + 1
      dsave = del2
!
!  Count number of changes in direction of monotonicity.
!
    else if ( temp == 0.0D+00 ) then

      if ( del2 /= 0.0D+00 ) then
        if ( pchst ( dsave, del2 ) < 0.0D+00 ) then
          ierr = ierr + 1
        end if
        dsave = del2
      end if
!
!  Use Brodlie modification of Butland formula.
!
    else

      hsumt3 = 3.0D+00 * hsum
      w1 = ( hsum + h1 ) / hsumt3
      w2 = ( hsum + h2 ) / hsumt3
      dmax = max ( abs ( del1 ), abs ( del2 ) )
      dmin = min ( abs ( del1 ), abs ( del2 ) )
      drat1 = del1 / dmax
      drat2 = del2 / dmax
      d(1,i) = dmin / ( w1 * drat1 + w2 * drat2 )

    end if

  end do
!
!  Set D(N) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  w1 = -h2 / hsum
  w2 = ( h2 + hsum ) / hsum
  d(1,n) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,n), del2 ) <= 0.0D+00 ) then
    d(1,n) = 0.0D+00
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
!
!  Need do this check only if monotonicity switches.
!
    dmax = 3.0D+00 * del2

    if ( abs ( dmax ) < abs ( d(1,n) ) ) then
      d(1,n) = dmax
    end if

  end if

  return
end
subroutine pchmc ( n, x, f, d, incfd, skip, ismon, ierr )

!*****************************************************************************80
!
!! PCHMC: piecewise cubic Hermite monotonicity checker.
!
!  Discussion:
!
!    PCHMC checks a cubic Hermite function for monotonicity.
!
!    To provide compatibility with PCHIM and PCHIC, the routine includes an
!    increment between successive values of the F and D arrays.
!
!  Modified:
!
!    31 August 2002
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at
!    least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in F and D.
!
!    Input/output, logical SKIP.  On input, should be set to TRUE if the
!    user wishes to skip checks for validity of preceding parameters, or
!    to FALSE otherwise.  This will save time in case these checks have
!    already been performed.  SKIP will be set to TRUE on normal return.
!
!    Output, integer ( kind = 4 ) ISMON(N), indicates the intervals on which the
!    piecewise cubic Hermite function is monotonic.
!    For data interval [X(I),X(I+1)], and 1 <= I <= N-1, ISMON(I) is
!    -1, if function is strictly decreasing;
!     0, if function is constant;
!     1, if function is strictly increasing;
!     2, if function is non-monotonic;
!     3, if unable to determine.  This means that the D values are near the
!       boundary of the monotonicity region.  A small increase produces
!       non-monotonicity; decrease, strict monotonicity.
!    ISMON(N) indicates whether the entire function is monotonic on [X(1),X(N)].
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 )  n

  integer ( kind = 4 ) chfmc
  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) delta
  real ( kind = 8 ) f(incfd,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ismon(n)
  integer ( kind = 4 ) nseg
  logical skip
  real ( kind = 8 ) x(n)

  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      call xerror ('pchmc -- number of data points less than two', ierr, 1)
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      call xerror ('pchmc -- increment less than one', ierr, 1)
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        call xerror ('pchmc -- x-array not strictly increasing', ierr, 1)
        return
      end if
    end do

    skip = .true.

  end if

  nseg = n - 1

  do i = 1, nseg

     delta = ( f(1,i+1) - f(1,i) ) / ( x(i+1) - x(i) )

     ismon(i) = chfmc ( d(1,i), d(1,i+1), delta )

     if ( i == 1 ) then
        ismon(n) = ismon(1)
     else
!
!  Need to figure out cumulative monotonicity from 'multiplication table'--
!
!                    *      i s m o n (i)
!                     *  -1   0   1   2   3
!               i      *--------------------*
!               s   -1 i -1  -1   2   2   3 i
!               m    0 i -1   0   1   2   3 i
!               o    1 i  2   1   1   2   3 i
!               n    2 i  2   2   2   2   2 i
!              (n)   3 i  3   3   3   2   3 i
!                      *--------------------*
!
!  If equal or already declared nonmonotonic, no change needed.
!
        if ( ismon(i) /= ismon(n) .and. ismon(n) /= 2 ) then
           if ( 1 < max ( ismon(i), ismon(n) ) ) then
!
!  At least one is either 'no' or 'maybe'.
!
              if ( ismon(i) == 2 ) then
                 ismon(n) = 2
              else
                 ismon(n) = 3
              end if
           else if ( ismon(i) * ismon(n) < 0 ) then
!
!  Both monotonic, but in opposite senses.
!
              ismon(n) = 2
           else
!
!  At this point, one is zero, the other is +-1.
!
              ismon(n) = ismon(n) + ismon(i)
           end if
        end if
     end if

  end do

  ierr = 0

  return
end
function pchqa ( n, x, f, d, a, b, ierr )

!*****************************************************************************80
!
!! PCHQA: easy to use cubic Hermite or spline integration.
!
!  Discussion:
!
!    PCHQA evaluates the definite integral of a cubic Hermite or spline
!    function over the interval [A, B].  This is an easy to use driver
!    for the routine PCHIA.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(N), the derivative values.  D(I) is the value
!    corresponding to X(I).
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.  The
!    interval [A,B] is normally contained in [X(1),X(N)], but this is
!    not required.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors).
!    1, if A is outside the interval [X(1),X(N)].
!    2, if B is outside the interval [X(1),X(N)].
!    3, if both of the above are true.  This means that either [A,B] contains
!       the data interval or the intervals do not intersect at all.
!    -1, if N < 2 .
!    -3, if the X array is not strictly increasing.
!
!    Output, real ( kind = 8 ) PCHQA, the value of the requested integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), save :: incfd = 1
  real ( kind = 8 ) pchia
  real ( kind = 8 ) pchqa
  logical, save :: skip = .true.
  real ( kind = 8 ) x(n)

  pchqa  =  pchia ( n, x, f, d, incfd, skip, a, b, ierr )

  return
end
subroutine pchsp ( ic, vc, n, x, f, d, incfd, wk, nwk, ierr )

!*****************************************************************************80
!
!! PCHSP sets derivatives needed for Hermite cubic spline interpolant.
!
!  Description:
!
!    PCHSP sets derivatives needed to determine the Hermite representation
!    of the cubic spline interpolant to given data, with specified boundary
!    conditions.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!    This is a modified version of Carl de Boor's cubic spline routine CUBSPL.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer-Verlag (new york, 1978).
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IC(2), specifies desired boundary conditions:
!    IC(1) = IBEG, desired condition at beginning of data.
!    0, to set D(1) so that the third derivative is continuous at X(2).
!      This is the "not a knot" condition provided by de Boor's cubic spline
!      routine CUBSPL, and is the default boundary condition here.
!    1, if first derivative at X(1) is given in VC(1).
!    2, if second derivative at X(1) is given in VC(1).
!    3, to use the 3-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 3.
!    4, to use the 4-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 4.
!    For the "natural" boundary condition, use ibeg=2 and vc(1)=0.
!    IC(2) = IEND, desired condition at end of data.
!    IEND may take on the same values as IBEG, but applied to derivative at
!    X(N).  In case IEND = 1 or 2, the value is given in VC(2).
!
!    Input, real ( kind = 8 ) VC(2), specifies desired boundary values,
!    as indicated above.  VC(1) need be set only if IC(1) = 1 or 2.
!    VC(2) need be set only if IC(2) = 1 or 2.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the dependent values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Output, real ( kind = 8 ) D(INCFD,N), the derivative values at the
!    data points.  These values will determine the cubic spline interpolant
!    with the requested boundary conditions.  The value corresponding to
!    X(I) is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in F and D.
!
!    Workspace, real ( kind = 8 ) WK(NWK).
!
!    Input, integer ( kind = 4 ) NWK, the size of WK, which must be at least 2 * N.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if IBEG < 0 or 4 < IBEG.
!    -5, if IEND < 0 or 4 < IEND.
!    -6, if both of the above are true.
!    -7, if NWK is too small.
!    -8, in case of trouble solving the linear system
!        for the interior derivative values.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) g
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ic(2)
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nwk
  real ( kind = 8 ) pchdf
  real ( kind = 8 ) stemp(3)
  real ( kind = 8 ) vc(2)
  real ( kind = 8 ) wk(2,n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtemp(4)

  if ( n < 2 ) then
    ierr = -1
    call xerror ('pchsp -- number of data points less than two', ierr, 1)
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    call xerror ('pchsp -- increment less than one', ierr, 1)
    return
  end if

  do j = 2, n
    if ( x(j) <= x(j-1) ) then
      ierr = -3
      call xerror ('pchsp -- x-array not strictly increasing', ierr, 1)
      return
    end if
  end do

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0

  if ( ibeg < 0 .or. 4 < ibeg ) then
    ierr = ierr - 1
  end if

  if ( iend < 0 .or. 4 < iend ) then
    ierr = ierr - 2
  end if

  if ( ierr < 0 ) then
    go to 5004
  end if
!
!  Function definition is ok -- go on.
!
  if ( nwk < 2 * n ) then
    go to 5007
  end if
!
!  Compute first differences of X sequence and store in wk(1,.). also,
!  compute first divided difference of data and store in wk(2,.).
!
  do j = 2, n
    wk(1,j) = x(j) - x(j-1)
    wk(2,j) = ( f(1,j) - f(1,j-1) ) / wk(1,j)
  end do
!
!  Set to default boundary conditions if N is too small.
!
  if ( n < ibeg ) then
    ibeg = 0
  end if

  if ( n < iend ) then
    iend = 0
  end if
!
!  Set up for boundary conditions.
!
  if ( ibeg == 1 .or. ibeg == 2 ) then
     d(1,1) = vc(1)
  else if ( 2 < ibeg ) then
!
!  Pick up first IBEG points, in reverse order.
!
     do j = 1, ibeg
       index = ibeg - j + 1
       xtemp(j) = x(index)
       if ( j < ibeg ) then
         stemp(j) = wk(2,index)
       end if
     end do

     d(1,1) = pchdf ( ibeg, xtemp, stemp, ierr )
     if ( ierr /= 0 ) then
       go to 5009
     end if

     ibeg = 1
  end if

  if ( iend == 1 .or. iend == 2 ) then
     d(1,n) = vc(2)
  else if ( 2 < iend ) then
!
!  Pick up last IEND points.
!
     do j = 1, iend
       index = n - iend + j
       xtemp(j) = x(index)
       if ( j < iend ) then
         stemp(j) = wk(2,index+1)
       end if
     end do

     d(1,n) = pchdf ( iend, xtemp, stemp, ierr )

     if ( ierr /= 0 ) then
       go to 5009
     end if

     iend = 1

  end if
!
!  Begin coding from cubspl
!
!  A tridiagonal linear system for the unknown slopes S(1:N) of
!  F at X(1:N) is generated and then solved by Gauss elimination,
!  with s(j) ending up in d(1,j), all j.
!  wk(1,.) and wk(2,.) are used for temporary storage.
!
!  Construct first equation from first boundary condition, of the form
!    wk(2,1) * s(1) + wk(1,1) * s(2) = D(1,1)
!
  if ( ibeg == 0 ) then
!
!  No condition at left end and N = 2.
!
     if ( n == 2 ) then
        wk(2,1) = 1.0D+00
        wk(1,1) = 1.0D+00
        d(1,1) = 2.0D+00 * wk(2,2)
!
!  Not-a-knot condition at left end and 2 < N.
!
     else
        wk(2,1) = wk(1,3)
        wk(1,1) = wk(1,2) + wk(1,3)
        d(1,1) =(( wk(1,2) + 2.0D+00 * wk(1,1) ) * wk(2,2) * wk(1,3) &
                             + wk(1,2)**2 * wk(2,3) ) / wk(1,1)
     end if
  else if ( ibeg == 1 ) then
!
!  Slope prescribed at left end.
!
     wk(2,1) = 1.0D+00
     wk(1,1) = 0.0D+00
  else
!
!  Second derivative prescribed at left end.
!
     wk(2,1) = 2.0D+00
     wk(1,1) = 1.0D+00
     d(1,1) = 3.0D+00 * wk(2,2) - 0.5D+00 * wk(1,2) * d(1,1)
  end if
!
!  If there are interior knots, generate the corresponding equations and
!  carry out the forward pass of Gauss elimination, after which the J-th
!  equation reads
!
!    wk(2,j) * s(j) + wk(1,j) * s(j+1) = d(1,j).
!
  if ( 1 < n-1 ) then
    do j = 2, n-1
        if ( wk(2,j-1) == 0.0D+00 ) then
          go to 5008
        end if
        g = -wk(1,j+1) / wk(2,j-1)
        d(1,j) = g * d(1,j-1) + 3.0D+00 &
          * ( wk(1,j) * wk(2,j+1) + wk(1,j+1) * wk(2,j) )
        wk(2,j) = g * wk(1,j-1) + 2.0D+00 * ( wk(1,j) + wk(1,j+1) )
    end do
  end if
!
!  Construct last equation from second boundary condition, of the form
!
!    (-g * wk(2,n-1)) * s(n-1) + wk(2,n) * s(n) = d(1,n)
!
!  If slope is prescribed at right end, one can go directly to back-
!  substitution, since arrays happen to be set up just right for it
!  at this point.
!
  if ( iend == 1 ) then
    go to 30
  end if

  if ( iend == 0 ) then
     if ( n == 2 .and. ibeg == 0 ) then
!
!  Not-a-knot at right endpoint and at left endpoint and N = 2.
!
        d(1,2) = wk(2,2)
        go to 30
     else if ( n == 2 .or. ( n == 3 .and. ibeg == 0 ) ) then
!
!  Either ( N = 3 and not-a-knot also at left) or (N=2 and *not*
!  not-a-knot at left end point).
!
        d(1,n) = 2.0D+00 * wk(2,n)
        wk(2,n) = 1.0D+00
        if ( wk(2,n-1) == 0.0D+00 ) then
          go to 5008
        end if
        g = -1.0D+00 / wk(2,n-1)
     else
!
!  Not-a-knot and 3 <= N, and either 3 < N or also not-a-
!  knot at left end point.
!
        g = wk(1,n-1) + wk(1,n)
!
!  Do not need to check following denominators (x-differences).
!
        d(1,n) = ( ( wk(1,n) + 2.0D+00 * g ) * wk(2,n) * wk(1,n-1) &
          + wk(1,n)**2 * ( f(1,n-1) - f(1,n-2) ) / wk(1,n-1) ) / g
        if ( wk(2,n-1) == 0.0D+00 ) then
          go to 5008
        end if
        g = -g / wk(2,n-1)
        wk(2,n) = wk(1,n-1)
     end if
  else
!
!  Second derivative prescribed at right endpoint.
!
     d(1,n) = 3.0D+00 *wk(2,n) + 0.5D+00 * wk(1,n) * d(1,n)
     wk(2,n) = 2.0D+00
     if ( wk(2,n-1) == 0.0D+00 ) then
       go to 5008
     end if
     g = -1.0D+00 / wk(2,n-1)
  end if
!
!  Complete forward pass of Gauss elimination.
!
  wk(2,n) = g * wk(1,n-1) + wk(2,n)

  if ( wk(2,n) == 0.0D+00 ) then
    go to 5008
  end if

  d(1,n) = ( g * d(1,n-1) + d(1,n) ) / wk(2,n)
!
!  Carry out back substitution.
!
   30 continue

  do j = n-1, 1, -1
    if ( wk(2,j) == 0.0D+00 ) then
      go to 5008
    end if
    d(1,j) = ( d(1,j) - wk(1,j) * d(1,j+1) ) / wk(2,j)
  end do

  return
!
!  error returns.
!
 5004 continue
!
!  ic out of range return.
!
  ierr = ierr - 3
  call xerror ('pchsp -- ic out of range', ierr, 1)
  return

 5007 continue
!
!  nwk too small return.
!
  ierr = -7
  call xerror ('pchsp -- work array too small', ierr, 1)
  return

 5008 continue
!
!  singular system.
!  theoretically, this can only occur if successive x-values
!  are equal, which should already have been caught (ierr=-3).
!
  ierr = -8
  call xerror ('pchsp -- singular linear system', ierr, 1)
  return
!
 5009 continue
!
!  error return from pchdf.
!  this case should never occur.
!
  ierr = -9
  call xerror ('pchsp -- error return from pchdf', ierr, 1)

  return
end
function pchst ( arg1, arg2 )

!*****************************************************************************80
!
!! PCHST: PCHIP sign-testing routine.
!
!  Discussion:
!
!    This routine essentially computes the sign of ARG1 * ARG2.
!
!    The object is to do this without multiplying ARG1 * ARG2, to avoid
!    possible over/underflow problems.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG1, ARG2, two values to check.
!
!    Output, real ( kind = 8 ) PCHST,
!    -1.0, if ARG1 and ARG2 are of opposite sign.
!     0.0, if either argument is zero.
!    +1.0, if ARG1 and ARG2 are of the same sign.
!
  implicit none

  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) pchst

  pchst = sign ( 1.0D+00, arg1 ) * sign ( 1.0D+00, arg2 )

  if ( arg1 == 0.0D+00 .or. arg2 == 0.0D+00 ) then
    pchst = 0.0D+00
  end if

  return
end
subroutine pchsw ( dfmax, iextrm, d1, d2, h, slope, ierr )

!*****************************************************************************80
!
!! PCHSW: the PCHCS switch excursion limiter.
!
!  Discussion:
!
!    This routine is called by PCHCS to adjust D1 and D2 if necessary to
!    insure that the extremum on this interval is not further than DFMAX
!    from the extreme data value.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DFMAX, the maximum allowed difference between
!    F(IEXTRM) and the cubic determined by the derivative values D1 and D2.
!    DFMAX should be nonnegative.
!
!    Input, integer ( kind = 4 ) IEXTRM, the index of the extreme data value,
!    which should be 1 or 2.
!
!    Input/output, real ( kind = 8 ) D1, D2, the derivative values at the
!    ends of the interval.  It is assumed that D1 * D2 <= 0.  On output,
!    the values may be modified if necessary to meet the restriction
!    imposed by DFMAX.
!
!    Input, real ( kind = 8 ) H, interval length.  H should be positive.
!
!    Input, real ( kind = 8 ) SLOPE, the data slope on the interval.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, assumption on D1 and D2 is not satisfied.
!    -2, quadratic equation locating extremum has negative descriminant
!        (should never occur).
!
!  Local variables:
!
!    RHO is the ratio of the data slope to the derivative being tested.
!
!    LAMBDA is the ratio of D2 to D1.
!
!    THAT = T-hat(rho) is the normalized location of the extremum.
!
!    PHI is the normalized value of P(X)-f1 at X = xhat = x-hat(rho),
!      where  that = (xhat - x1)/h .
!      that is, p(xhat)-f1 = D * H * PHI,  where d=d1 or d2.
!      similarly,  p(xhat)-f2 = d*h*(phi-rho) .
!
!    SMALL should be a few orders of magnitude greater than macheps.
!
  implicit none

  real ( kind = 8 ) cp
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dfmax
  real ( kind = 8 ) dmax
  real ( kind = 8 ), parameter :: fact = 100.0D+00
  real ( kind = 8 ) h
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iextrm
  real ( kind = 8 ) lambda
  real ( kind = 8 ) nu
  real ( kind = 8 ) phi
  real ( kind = 8 ) radcal
  real ( kind = 8 ) rho
  real ( kind = 8 ) sigma
  real ( kind = 8 ) slope
  real ( kind = 8 ) small
  real ( kind = 8 ) that
  real ( kind = 8 ), parameter :: third = 0.33333D+00

  small = fact * epsilon ( 1.0D+00 )

  if ( d1 == 0.0D+00 ) then
!
!  Special case -- D1 == 0.0D+00 .
!
!  If D2 is also zero, this routine should not have been called.
!
     if ( d2 == 0.0D+00 ) then
       ierr = -1
       call xerror ('pchsw -- d1 and/or d2 invalid', ierr, 1)
       return
     end if

     rho = slope / d2
!
!  Extremum is outside interval when 1/3 <= RHO.
!
     if ( third <= rho ) then
       ierr = 0
       return
     end if

     that = ( 2.0D+00 * ( 3.0D+00 * rho - 1.0D+00 ) ) &
       / ( 3.0D+00 * ( 2.0D+00 * rho - 1.0D+00 ) )

     phi = that**2 * ( ( 3.0D+00 * rho - 1.0D+00 ) / 3.0D+00 )
!
!  Convert to distance from F2 if IEXTRM /= 1.
!
     if ( iextrm /= 1 ) then
       phi = phi - rho
     end if
!
!  Test for exceeding limit, and adjust accordingly.
!
     dmax = dfmax / ( h * abs ( phi ) )
     if ( dmax < abs ( d2 ) ) then
       d2 = sign ( dmax, d2 )
     end if

  else

     rho = slope / d1
     lambda = -d2 / d1
     if ( d2 == 0.0D+00 ) then
!
!  Special case -- D2 == 0.0D+00 .
!
!  Extremum is outside interval when 1/3 <= RHO.
!
        if ( third <= rho ) then
          ierr = 0
          return
        end if

        cp = 2.0D+00 - 3.0D+00 * rho
        nu = 1.0D+00 - 2.0D+00 * rho
        that = 1.0D+00 / ( 3.0D+00 * nu )

     else

        if ( lambda <= 0.0D+00 ) then
          ierr = -1
          call xerror ('pchsw -- d1 and/or d2 invalid', ierr, 1)
          return
        end if
!
!  Normal case, D1 and D2 both nonzero, opposite signs.
!
        nu = 1.0D+00 - lambda - 2.0D+00 * rho
        sigma = 1.0D+00 - rho
        cp = nu + sigma

        if ( small < abs ( nu ) ) then

          radcal = ( nu - ( 2.0D+00 * rho + 1.0D+00 ) ) * nu + sigma**2

          if ( radcal < 0.0D+00 ) then
            ierr = -2
            call xerror ( 'pchsw -- negative radical', ierr, 1)
            return
          end if

          that = ( cp - sqrt ( radcal ) ) / ( 3.0D+00 * nu )

        else

          that = 1.0D+00 / ( 2.0D+00 * sigma )

        end if

     end if

     phi = that * ( ( nu * that - cp ) * that + 1.0D+00 )
!
!  Convert to distance from F2 if IEXTRM /= 1.
!
     if ( iextrm /= 1 ) then
       phi = phi - rho
     end if
!
!  Test for exceeding limit, and adjust accordingly.
!
     dmax = dfmax / ( h * abs ( phi ) )

     if ( dmax < abs ( d1 ) ) then
        d1 = sign ( dmax, d1 )
        d2 = -lambda * d1
     end if

  end if

  ierr = 0

  return
end
subroutine q1da ( f, a, b, eps, r, e, kf, iflag )

!*****************************************************************************80
!
!! Q1DA approximates the definite integral of a function of one variable.
!
!  Discussion:
!
!    A small amount of randomization is built into this program.
!    Calling Q1DA a few times in succession will give different
!    but hopefully consistent results.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, external F, the name of the function which evaluates the integrand.
!    F must have the form:
!
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the integration interval.
!
!    Input, real ( kind = 8 ) EPS, the desired accuracy.
!
!    Output, real ( kind = 8 ) R, the estimate of the value of the integral.
!
!    Output, real ( kind = 8 ) E, an estimate of the error, |integral-R|.
!
!    Output, integer ( kind = 4 ) KF, the cost of the integration, measured in
!    number of evaluations of your integrand.  KF will always be at least 30.
!
!    Output, integer ( kind = 4 ) IFLAG, termination flag.
!    0, normal completion, E < EPS and E < EPS * |R|.
!    1, normal completion, E < EPS, but EPS * |R| < E.
!    2, normal completion, E < EPS * |R|, but EPS < E.
!    3, normal completion, but EPS was too small to satisfy absolute
!       or relative error request.
!    4, aborted calculation because of serious rounding error.  Probably E
!       and R are consistent.
!    5, aborted calculation because of insufficient storage.  R and E
!       are consistent.
!    6, aborted calculation because of serious difficulties meeting your
!       error request.
!    7, aborted calculation because EPS was set <= 0.0
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 50

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) e
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fmax
  real ( kind = 8 ) fmin
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) kf
  integer ( kind = 4 ) nint
  real ( kind = 8 ) r
  logical rst
  real ( kind = 8 ) w(nmax,6)

  nint = 1
  rst = .false.
  call q1dax ( f, a, b, eps, r, e, nint, rst, w, nmax, fmin, fmax, kf, iflag )

  return
end
subroutine q1dax ( f, a, b, eps, r, e, nint, rst, w, nmax, fmin, fmax, kf, &
  iflag )

!*****************************************************************************80
!
!! Q1DAX approximates the integral of a function of one variable.
!
!  Discussion:
!
!    For an easier to use routine see Q1DA.
!
!    Capabilities of Q1DAX, beyond those of Q1DA, include:
!    * the ability to restart a calculation to greater accuracy without
!      penalty...
!    * the ability to specify an initial partition of the integration
!      interval...
!    * the ability to increase the work space to handle more difficult
!      problems...
!    * output of largest/smallest integrand value for applications such
!      as scaling graphs...
!
!    A small amount of randomization is built into this program.
!    Calling Q1DAX a few times in succession will give different
!    but hopefully consistent results.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function, of the form
!
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the integration interval.
!
!    Input, real ( kind = 8 ) EPS, the desired accuracy.
!
!    Output, real ( kind = 8 ) R, the estimate of the value of the integral.
!
!    Output, real ( kind = 8 ) E, an estimate of the error, |integral-R|.
!
!    Input/output, integer ( kind = 4 ) NINT.  On input, NINT must be set to the
!    number of subintervals in the initial partition of [a,b].
!    For most problems this is just 1, for the interval [a,b] itself.
!    NINT must be less than NMAX.  NINT is useful if you would like to
!    help Q1DAX locate a difficult spot on [a,b].  In this regard NINT is
!    used along with the array w (see below).  if you set
!    NINT=1 it is not necessary to be concerned
!    with W, except that it must be dimensioned...
!    as an example of more general applications,
!    if [a,b]=[0,1] but the integrand jumps at 0.3,
!    it would be wise to set NINT=2 and then set
!      w(1,1)=0.0D+00  (left endpoint)
!      w(2,1)=0.3  (singular point)
!      w(3,1)=1.0D+00  (right endpoint)
!    If you set NINT greater than 1, be sure to set
!      W(1,1)      = A
!      W(NINT+1,1) = B
!
!    As an output quantity, NINT gives the number of subintervals in the
!    final partition of [A,B].
!
!    Input, logical RST, indicates first call or restart.
!    FALSE, for the initial call to Q1DAX.
!    TRUE for a subsequent call, for the same problem, for which
!    more accuracy is desired (a smaller EPS).  A restart only
!    makes sense if the preceding call returned with a value of IFLAG
!    less than 3.  On a restart you may not change the values of any
!    other arguments in the call sequence, except EPS.
!
!    Input, real ( kind = 8 ) W(NMAX,6).  an adequate value of
!    nmax is 50.  if you set NINT greater than 1, you must also
!    initialize w, see nint above.
!
!    Input, integer ( kind = 4 ) NMAX, the leading dimension of W.  This is
!    also equal to the maximum number of subintervals permitted in the internal
!    partition of [A,B].  A value of 50 is ample for most problems.
!
!    Output, real ( kind = 8 ) FMIN, FMAX, the smallest and largest values of
!    the integrand which occurred during the calculation.  The actual
!    integrand range on [A,B] may, of course, be greater but probably not
!    by more than 10%.
!
!    Output, integer ( kind = 4 ) KF, the cost of the integration, measured in
!    number of evaluations of your integrand.  KF will always be at least 30.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, normal completion, e < eps  and  e < eps * abs ( r )
!    1, normal completion, e < eps, but eps * abs ( r ) < E;
!    2. normal completion, e < eps * abs ( r ), but EPS < E.
!    3. normal completion but eps was too small to satisfy absolute
!      or relative error request.
!    4, aborted calculation because of serious rounding
!      error.  probably e and r are consistent.
!    5, aborted calculation because of insufficient storage.
!      r and e are consistent.  perhaps increasing nmax
!      will produce better results.
!    6, aborted calculation because of serious difficulties
!      meeting your error request.
!    7, aborted calculation because either eps, nint or nmax
!      has been set to an illegal value.
!    8, aborted calculation because you set NINT greater than 1, but forgot
!      to set w(1,1)=a  and  w(nint+1,1)=b
!
  implicit none

  integer ( kind = 4 ) nmax

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) e
  real ( kind = 8 ) eb
  real ( kind = 8 ) epmach
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fmax
  real ( kind = 8 ) fmaxl
  real ( kind = 8 ) fmaxr
  real ( kind = 8 ) fmin
  real ( kind = 8 ) fminl
  real ( kind = 8 ) fminr
  real ( kind = 8 ) fmn
  real ( kind = 8 ) fmx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iroff
  integer ( kind = 4 ) kf
  integer ( kind = 4 ) loc
  integer ( kind = 4 ) mxtry
  integer ( kind = 4 ) nint
  real ( kind = 8 ) r
  real ( kind = 8 ) rab
  real ( kind = 8 ) rabs
  real ( kind = 8 ) rav
  logical rst
  real ( kind = 8 ) t
  real ( kind = 8 ) te
  real ( kind = 8 ) te1
  real ( kind = 8 ) te2
  real ( kind = 8 ) tr
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) uflow
  real ( kind = 8 ) uni
  real ( kind = 8 ) w(nmax,6)
  real ( kind = 8 ) xm

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  mxtry = nmax / 2
!
!  In case there is no more room, we can toss out easy intervals,
!  at most MXTRY times.
!
  if ( a == b ) then
    r = 0.0D+00
    e = 0.0D+00
    nint = 0
    iflag = 0
    kf = 1
    fmin = f(a)
    fmax = fmin
    go to 20
  end if

  if ( rst ) then
    if ( iflag < 3 ) then
      eb = max ( 100.0D+00 * uflow, max ( eps, 50.0D+00 * epmach ) * abs ( r ) )
      do i = 1, nint
        if ( ( eb * ( w(i,2) - w(i,1) ) / ( b - a ) ) < abs ( w(i,3) ) ) then
          w(i,3) = abs ( w(i,3) )
        else
          w(i,3) = -abs ( w(i,3) )
        end if
      end do
      go to 15
    else
      go to 20
    end if
  end if

  kf = 0

  if ( eps <= 0.0D+00 .or. nint <= 0 .or. nmax <= nint ) then
    iflag = 7
    go to 20
  end if

  if ( nint == 1 ) then

      w(1,1) = a
      w(2,2) = b
      w(1,5) = a
      w(1,6) = b
      w(2,5) = a
      w(2,6) = b
!
!  Select the first subdivision randomly.
!
      w(1,2) = a + ( b - a ) / 2.0D+00 * ( 2.0D+00 * uni() + 7.0D+00 ) / 8.0D+00
      w(2,1) = w(1,2)
      nint = 2

    else

      if ( w(1,1) /= a .or. w(nint+1,1) /= b ) then
        iflag = 8
        go to 20
      end if

      w(1,5) = a

      do i = 1, nint
        w(i,2) = w(i+1,1)
        w(i,5) = w(i,1)
        w(i,6) = w(i,2)
      end do

  end if

  iflag = 0
  iroff = 0
  rabs = 0.0D+00

  do i = 1, nint

    call gl15t ( f, w(i,1), w(i,2), w(i,5), w(i,6), &
      w(i,4), w(i,3), rab, rav, fmn, fmx )

    kf = kf+15

    if ( i == 1 ) then
      r = w(i,4)
      e = w(i,3)
      rabs = rabs + rab
      fmin = fmn
      fmax = fmx
    else
      r = r + w(i,4)
      e = e + w(i,3)
      rabs = rabs + rab
      fmax = max ( fmax, fmx )
      fmin = min ( fmin, fmn )
    end if

  end do

  w(nint+1:nmax,3) = 0.0D+00

   15 continue
!
!  main subprogram loop
!
  if ( abs ( r ) <= 100.0D+00 * epmach * rabs .and. e < eps ) then
    go to 20
  end if

  eb = max ( 100.0D+00 * uflow, max ( eps, 50.0D+00 * epmach ) * abs ( r ) )

  if ( e <= eb ) then
    go to 20
  end if

  if ( nint < nmax ) then

    nint = nint + 1
    c = nint

  else

    c = 0

16  continue

    if ( c == nmax .or. mxtry <= 0 ) then
      iflag = 5
      go to 20
    end if

    c = c + 1

    if ( 0.0D+00 < w(c,3) ) then
      go to 16
    end if
!
!  Found an interval to throw out.
!
    mxtry = mxtry - 1

  end if

  loc = idamax ( nint, w(1,3), 1 )
  xm = w(loc,1) + ( w(loc,2) - w(loc,1) ) / 2.0D+00

  if ( ( max ( abs ( w(loc,1) ), abs ( w(loc,2) ) ) ) > &
    ( ( 1.0D+00 + 100.0D+00 * epmach ) &
    * ( abs ( xm ) + 0.1D+04 * uflow ) ) ) then

        call gl15t ( f, w(loc,1), xm, w(loc,5), w(loc,6), &
          tr1, te1, rab, rav, fminl, fmaxl )

        kf = kf + 15

        if ( te1 < ( eb * ( xm - w(loc,1) ) / ( b - a ) ) ) then
          te1 = -te1
        end if

        call gl15t ( f, xm, w(loc,2), w(loc,5), w(loc,6), &
          tr2, te2, rab, rav, fminr, fmaxr )

        kf = kf + 15
        fmin = min ( fmin, fminl, fminr )
        fmax = max ( fmax, fmaxl, fmaxr )

        if ( te2 < ( eb * ( w(loc,2) - xm ) / ( b - a ) ) ) then
          te2 = -te2
        end if

        te = abs ( w(loc,3) )
        tr = w(loc,4)
        w(c,3) = te2
        w(c,4) = tr2
        w(c,1) = xm
        w(c,2) = w(loc,2)
        w(c,5) = w(loc,5)
        w(c,6) = w(loc,6)
        w(loc,3) = te1
        w(loc,4) = tr1
        w(loc,2) = xm
        e = e - te + ( abs ( te1 ) + abs ( te2 ) )
        r = r - tr + ( tr1 + tr2 )
        if ( abs ( abs ( te1 ) + abs ( te2 ) - te ) < 0.001D+00 * te ) then
            iroff = iroff + 1
            if ( 10 <= iroff ) then
              iflag = 4
              go to 20
            end if
        end if
      else

        if ( w(loc,3) < eb ) then
          w(loc,3) = 0.0D+00
        else
          iflag = 6
          go to 20
        end if

  end if

  go to 15

   20 continue

  if ( 4 <= iflag ) then
    return
  end if

  t = eps * abs ( r )

  if ( eps < e .and. t < e ) then
    iflag = 3
    return
  end if

  if ( eps < e .and. e < t ) then
    iflag = 2
    return
  end if

  if ( e < eps .and. t < e ) then
    iflag = 1
    return
  end if

  iflag = 0

  return
end
subroutine qagi ( f, bound, inf, epsabs, epsrel, result, abserr, neval, &
  ier, limit, lenw, last, iwork, work )

!*****************************************************************************80
!
!! QAGI approximates an integral over an infinite or semi-infinite interval.
!
!  Discussion:
!
!    QAGI calculates an approximation RESULT to a given integral:
!
!      I = integral of F(X) over (bound,+infinity), or
!      I = integral of F(X) over (-infinity,bound), or
!      I = integral of F(X) over (-infinity,+infinity)
!
!    hopefully satisfying following claim for accuracy:
!
!      abs ( i - result ) <= max ( epsabs, epsrel * abs ( i ) ).
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external F, the name of the routine that evaluates the function,
!    of the form
!
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!
!    Input, real ( kind = 8 ) BOUND, the value of the finite endpoint of
!    the integration range, if any, that is, if INF is 1 or -1.
!
!    Input, integer ( kind = 4 ) INF, indicates the type of integration range.
!    1:  (  BOUND,    +Infinity),
!    -1: ( -Infinity,  BOUND),
!    2:  ( -Infinity, +Infinity).
!
!    Input, real ( kind = 8 ) EPSABS, the absolute accuracy requested.
!
!    Input, real ( kind = 8 ) EPSREL, the relative accuracy requested
!    If EPSABS <= 0 and EPSREL < max ( 50 * epsilon, 0.5D-14),
!    the routine will end with IER = 6.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, the estimate of the modulus of the
!    absolute error, which should equal or exceed | I - RESULT |.
!
!    Output, integer ( kind = 4 ) NEVAL, the number of integrand evaluations.
!
!         on return
!
!            ier    - integer ( kind = 4 )
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                   - 0 < ier, abnormal termination of the routine. the
!                             estimates for result and error are less
!                             reliable. it is assumed that the requested
!                             accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. if
!                             the position of a local difficulty can be
!                             determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table.
!                             it is assumed that the requested tolerance
!                             cannot be achieved, and that the returned
!                             result is the best which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             (epsabs <= 0 and
!                              epsrel < max ( 50 * EPSILON, 0.5d-28 ) )
!                              or limit<1 or leniw<limit*4.
!                             result, abserr, neval, last are set to
!                             zero. exept when limit or leniw is
!                             invalid, iwork(1), work(limit*2+1) and
!                             work(limit*3+1) are set to zero, work(1)
!                             is set to a and work(limit+1) to b.
!
!         dimensioning parameters
!            limit - integer ( kind = 4 )
!                    dimensioning parameter for iwork
!                    limit determines the maximum number of subintervals
!                    in the partition of the given integration interval
!                    (a,b), 1 <= LIMIT. in many cases limit = 100 is ok.
!                    if limit<1, the routine will end with ier = 6.
!
!            lenw  - integer ( kind = 4 )
!                    dimensioning parameter for work
!                    lenw must be at least limit*4.
!                    if lenw<limit*4, the routine will end
!                    with ier = 6.
!
!            last  - integer ( kind = 4 )
!                    on return, last equals the number of subintervals
!                    produced in the subdivision process, which
!                    determines the number of significant elements
!                    actually in the work arrays.
!
!         work arrays
!            iwork - integer ( kind = 4 )
!                    vector of dimension at least limit, the first
!                    k elements of which contain pointers
!                    to the error estimates over the subintervals,
!                    such that work(limit*3+iwork(1)),... ,
!                    work(limit*3+iwork(k)) form a decreasing
!                    sequence, with k = last if last<=(limit/2+2), and
!                    k = limit+1-last otherwise
!
!            work  - real
!                    vector of dimension at least lenw
!                    on return
!                    work(1), ..., work(last) contain the left
!                   end points of the subintervals in the
!                     partition of (a,b),
!                    work(limit+1), ..., work(limit+last) contain
!                     the right end points,
!                    work(limit*2+1), ...,work(limit*2+last) contain the
!                     integral approximations over the subintervals,
!                    work(limit*3+1), ..., work(limit*3)
!                     contain the error estimates.
!
  implicit none

  integer ( kind = 4 ) lenw
  integer ( kind = 4 ) limit

  real ( kind = 8 ) abserr
  real ( kind = 8 ) bound
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inf
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lvl
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) l3
  integer ( kind = 4 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) work(lenw)

  ier = 6
  lvl = 0
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00

  if ( limit < 1 .or. lenw < 4 * limit ) then
    lvl = 1
    call xerror ( 'QAGI - abnormal return',ier, lvl )
    return
  end if

  l1 = limit + 1
  l2 = limit + l1
  l3 = limit + l2

  call qagie ( f, bound, inf, epsabs, epsrel, limit, result, abserr, &
    neval, ier, work(1), work(l1), work(l2), work(l3), iwork, last )

  if ( ier /= 0 ) then
    call xerror ( 'QAGI - abnormal return from QAGIE', ier, lvl )
  end if

  return
end
subroutine qagie ( f, bound, inf, epsabs, epsrel, limit, result, abserr, &
  neval, ier, alist, blist, rlist, elist, iord, last )

!*****************************************************************************80
!
!! QAGIE is called by QAGI to compute integrals on an infinite interval.
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, external F, the name of the user-supplied function, of the form:
!
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!
!    Input, real ( kind = 8 ) BOUND, the value of the finite endpoint of
!    the integration range, if any, that is, if INF is 1 or -1.
!
!    Input, integer ( kind = 4 ) INF, indicates the type of integration range.
!    1:  (  BOUND,    +Infinity),
!    -1: ( -Infinity,  BOUND),
!    2:  ( -Infinity, +Infinity).
!
!    Input, real ( kind = 8 ) EPSABS, the absolute accuracy requested.
!
!    Input, real ( kind = 8 ) EPSREL, the relative accuracy requested
!    If EPSABS <= 0 and EPSREL < max ( 50 * epsilon, 0.5D-14),
!    the routine will end with IER = 6.
!
!    Input, integer ( kind = 4 ) LIMIT, the maximum number of subintervals.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, the estimate of the modulus of the
!    absolute error, which should equal or exceed | I - RESULT |.
!
!    Output, integer ( kind = 4 ) NEVAL, the number of integrand evaluations.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, no error occurred.
!
!    ?, real ( kind = 8 ) ALIST(LIMIT), ?
!
!    ?, real ( kind = 8 ) BLIST(LIMIT), ?
!
!    ?, real ( kind = 8 ) RLIST(LIMIT), ?
!
!    ?, real ( kind = 8 ) ELIST(LIMIT), ?
!
!    ?, integer ( kind = 4 ) IORD(LIMIT), ?
!
!    ?, integer ( kind = 4 ) LAST, ?
!
!  Local parameters:
!
!    LIMEXP is the size of the epsilon table that can be generated in EA.
!
  implicit none

  integer ( kind = 4 ), parameter :: limexp = 50
  integer ( kind = 4 ) limit

  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) boun
  real ( kind = 8 ) bound
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) correc
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epmach
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extrap
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ierro
  integer ( kind = 4 ) inf
  integer ( kind = 4 ) iord(limit)
  integer ( kind = 4 ) iroff1
  integer ( kind = 4 ) iroff2
  integer ( kind = 4 ) iroff3
  integer ( kind = 4 ) jupbnd
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ksgn
  integer ( kind = 4 ) ktmin
  integer ( kind = 4 ) last
  logical lerr
  integer ( kind = 4 ) maxerr
  integer ( kind = 4 ) neval
  logical newflg
  logical noext
  integer ( kind = 4 ) nrmax
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(limexp+7)
  real ( kind = 8 ) small
  real ( kind = 8 ) uflow

  epmach = epsilon ( epmach )
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  alist(1) = 0.0D+00
  blist(1) = 1.0D+00
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00
  iord(1) = 0
  newflg = .true.

  if ( epsabs <= 0.0D+00 .and. epsrel < max ( 0.5D+02 * epmach, 0.5D-14 ) ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral
!
!  Determine the interval to be mapped onto (0,1).
!  if inf = 2 the integral is computed as i = i1+i2, where
!  i1 = integral of f over (-infinity,0),
!  i2 = integral of f over (0,+infinity).
!
  if ( inf == 2 ) then
    boun = 0.0D+00
  else
    boun = bound
  end if

  call qk15i ( f, boun, inf, 0.0D+00, 1.0D+00, result, abserr, defabs, resabs )
!
!  Test on accuracy
!
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )

  if ( abserr <= 1.0D+02 * epmach * defabs .and. errbnd < abserr ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 ) then
    go to 130
  end if

  if ( ( abserr <= errbnd .and. abserr /= resabs ) .or. abserr == 0.0D+00 ) then
    go to 130
  end if
!
!  Initialization.
!
  uflow = 2.0D+00 * tiny ( uflow )
  lerr = .false.
  call ea ( newflg, result, limexp, reseps, abseps, rlist2, ierr )
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  ktmin = 0
  extrap = .false.
  noext = .false.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( ( 1.0D+00 - 0.5D+02 * epmach ) * defabs <= dres ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5D+00 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk15i ( f, boun, inf, a1, b1, area1, error1, resabs, defab1 )
    call qk15i ( f, boun, inf, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral
!  and error and test for accuracy.
!
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)

    if ( defab1 == error1 .or. defab2 == error2 ) then
      go to 15
    end if

    if ( 0.1D-04 * abs ( area12 ) < abs ( rlist(maxerr) - area12 ) .or. &
      erro12 < 0.99D+00 * errmax ) then
      go to 10
    end if

    if ( extrap ) then
      iroff2 = iroff2+1
    else
      iroff1 = iroff1+1
    end if

10  continue

    if ( 10 < last .and. errmax < erro12 ) then
      iroff3 = iroff3+1
    end if

15  continue

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel * abs ( area ) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( 10 <= iroff1 + iroff2 .or. 20 <= iroff3 ) then
      ier = 2
    end if

    if ( 5 <= iroff2 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at some points of the integration range.
!
    if ( max ( abs ( a1 ), abs ( b2 ) ) <= ( 1.0D+00 + 0.1D+03 * epmach )* &
      ( abs ( a2 ) + 0.1D+04 * uflow ) ) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QPSRT to maintain the descending ordering
!  in the list of error estimates and select the
!  subinterval with nrmax-th largest error estimate (to be
!  bisected next).
!
    call qpsrt ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) then
      go to 115
    end if

    if ( ier /= 0 ) then
      go to 100
    end if

    if ( last == 2 ) then
      go to 80
    end if

    if ( noext ) then
      go to 90
    end if

    erlarg = erlarg - erlast

    if ( small < abs ( b1 - a1 ) ) then
      erlarg = erlarg + erro12
    end if

    if ( extrap ) then
      go to 40
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
      go to 90
    end if

    extrap = .true.
    nrmax = 2

40  continue

    if ( ierro == 3 .or. erlarg <= ertest ) then
      go to 60
    end if
!
!  the smallest interval has the largest error.
!  before bisecting decrease the sum of the errors
!  over the larger intervals (erlarg) and perform
!  extrapolation.
!
    id = nrmax
    jupbnd = last

    if ( (2+limit/2) < last ) then
      jupbnd = limit+3-last
    end if

    do k = id, jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
        go to 90
      end if
      nrmax = nrmax + 1
    end do
!
!  Perform extrapolation.
!
   60   continue

    call ea ( newflg, area, limexp, reseps, abseps, rlist2, ierr )
    ktmin = ktmin+1

    if ( 5 < ktmin .and. ( abserr < 0.1D-02 * errsum ) .and. ( lerr ) ) then
      ier = 5
    end if

    if ( abserr <= abseps .and. lerr ) then
      go to 70
    end if

    ktmin = 0
    abserr = abseps
    lerr = .true.
    result = reseps
    correc = erlarg
    ertest = max ( epsabs, epsrel * abs ( reseps ) )

    if ( abserr <= ertest .and. lerr ) then
      go to 100
    end if
!
!  Prepare bisection of the smallest interval.
!
   70   continue

    if ( rlist2(limexp+3) == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      go to 100
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 0.5D+00
    erlarg = errsum
    go to 90

   80   continue

    small = 0.375D+00
    erlarg = errsum
    ertest = errbnd
    call ea ( newflg, area, limexp, reseps, abseps, rlist2, ierr )

   90 continue

  end do
!
!  Set final result and error estimate.
!
  100 continue

  if ( .not. lerr ) then
    go to 115
  end if

  if ( ier + ierro == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0D+00 .and. area /= 0.0D+00 ) then
    go to 105
  end if

  if ( errsum < abserr ) then
    go to 115
  end if

  if ( area == 0.0D+00 ) then
    go to 130
  end if

  go to 110

  105 continue

  if ( errsum / abs ( area ) < abserr / abs ( result ) ) then
    go to 115
  end if
!
!  Test on divergence
!
  110 continue

    if ( ksgn == (-1) .and. &
      max ( abs ( result ), abs ( area ) ) <= defabs * 0.1D-01 ) then
      go to 130
    end if

  if ( ( result / area ) < 0.1D-01 .or. &
    0.1D+03 < ( result / area ) .or. &
    abs ( area ) < errsum ) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
  115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum
  130 continue

  neval = 30 * last - 15

  if ( inf == 2 ) then
    neval = 2 * neval
  end if

  if ( 2 < ier ) then
    ier = ier-1
  end if

  return
end
subroutine qform ( m, n, q, ldq )

!*****************************************************************************80
!
!! QFORM produces the explicit QR factorization of a matrix.
!
!  Discussion:
!
!    The QR factorization of a matrix is usually accumulated in implicit
!    form, that is, as a series of orthogonal transformations of the
!    original matrix.  This routine carries out those transformations,
!    to explicitly exhibit the factorization construced by QRFAC.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A and the order of Q.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) Q(LDQ,M).  Q is an M by M array.
!    On input the full lower trapezoid in the first min ( M, N ) columns of Q
!    contains the factored form.
!    On output, Q has been accumulated into a square matrix.
!
!    Input, integer ( kind = 4 ) LDQ, not less
!    than M which specifies the leading dimension of the array Q.
!
  implicit none

  integer ( kind = 4 ) ldq
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) minmn
  real ( kind = 8 ) q(ldq,m)
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa(m)

  minmn = min ( m, n )

  do j = 2, minmn
    q(1:j-1,j) = 0.0D+00
  end do
!
!  Initialize remaining columns to those of the identity matrix.
!
  do j = n+1, m
    q(1:m,j) = 0.0D+00
    q(j,j) = 1.0D+00
  end do
!
!  Accumulate Q from its factored form.
!
  do l = 1, minmn

    k = minmn - l + 1

    wa(k:m) = q(k:m,k)

    q(k:m,k) = 0.0D+00
    q(k,k) = 1.0D+00

    if ( wa(k) /= 0.0D+00 ) then

      do j = k, m
        temp = dot_product ( wa(k:m), q(k:m,j) ) / wa(k)
        q(k:m,j) = q(k:m,j) - temp * wa(k:m)
      end do

    end if

  end do

  return
end
subroutine qk15 ( f, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external F, the name of the user-supplied function, of the form:
!
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimate of the integral I.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of the modulus of the
!    absolute error, which should not exceed |I-RESULT|.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of
!    the absolute value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral
!    | F-I/(B-A) | over [A,B].
!
!  Local parameters:
!
!    centr  - mid point of the interval
!    hlgth  - half-length of the interval
!    absc   - abscissa
!    fval*  - function value
!    resg   - result of the 7-point Gauss formula
!    resk   - result of the 15-point Kronrod formula
!    reskh  - approximation to the mean value of f over (a,b), i.e. to i/(b-a)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  real ( kind = 8 ) epmach
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(7)
  real ( kind = 8 ) fv2(7)
  real ( kind = 8 ) hlgth
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtw
  integer ( kind = 4 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) uflow
  real ( kind = 8 ), parameter, dimension ( 4 ) :: wg = (/ &
    0.1294849661688697D+00, &
    0.2797053914892767D+00, &
    0.3818300505051189D+00, &
    0.4179591836734694D+00 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: wgk = (/ &
     0.2293532201052922D-01,   0.6309209262997855D-01, &
     0.1047900103222502D+00,   0.1406532597155259D+00, &
     0.1690047266392679D+00,   0.1903505780647854D+00, &
     0.2044329400752989D+00,   0.2094821410847278D+00 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: xgk = (/ &
     0.9914553711208126D+00,   0.9491079123427585D+00, &
     0.8648644233597691D+00,   0.7415311855993944D+00, &
     0.5860872354676911D+00,   0.4058451513773972D+00, &
     0.2077849550078985D+00,   0.0D+00              /)

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00 * ( a + b )
  hlgth = 0.5D+00 * ( b - a )
  dhlgth = abs ( hlgth )
!
!  Compute the 15-point Kronrod approximation to
!  the integral, and estimate the absolute error.
!
  fc = f(centr)
  resg = fc * wg(4)
  resk = fc * wgk(8)
  resabs = abs ( resk )

  do j = 1, 3
    jtw = j * 2
    absc = hlgth * xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1 + fval2
    resg = resg + wg(j) * fsum
    resk = resk + wgk(jtw) * fsum
    resabs = resabs + wgk(jtw) * ( abs ( fval1 ) + abs ( fval2 ) )
  end do

  do j = 1, 4
    jtwm1 = j * 2 - 1
    absc = hlgth * xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1 + fval2
    resk = resk + wgk(jtwm1) * fsum
    resabs = resabs + wgk(jtwm1) * ( abs ( fval1 ) + abs ( fval2 ) )
  end do

  reskh = resk * 0.5D+00
  resasc = wgk(8) * abs ( fc - reskh )

  do j = 1, 7
    resasc = resasc + wgk(j) * &
      ( abs ( fv1(j) - reskh ) + abs ( fv2(j) - reskh ) )
  end do

  result = resk * hlgth
  resabs = resabs * dhlgth
  resasc = resasc * dhlgth
  abserr = abs ( ( resk - resg ) * hlgth )

  if ( resasc /= 0.0D+00 .and. abserr /= 0.0D+00 ) then
    abserr = resasc * min ( 1.0D+00, ( 0.2D+03 * abserr / resasc )**1.5D+00 )
  end if

  if ( uflow / ( 0.5D+02 * epmach ) < resabs ) then
    abserr = max ( ( epmach * 0.5D+02 ) * resabs, abserr )
  end if

  return
end
subroutine qk15i ( f, boun, inf, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK15I applies a 15 point Gauss Kronrod quadrature rule.
!
!  Discussion:
!
!    The integration interval is assumed to be infinite or semi-infinite.
!    A mapping is applied to the interval to bring it into [0,1].
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external F, the name of the user-supplied function, of the form:
!
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!
!    Input, real ( kind = 8 ) BOUN, the value of the finite endpoint of the
!    integration range, if any, that is, if INF is 1 or -1.
!
!    Input, integer ( kind = 4 ) INF, indicates the type of integration range.
!    1:  (  BOUND,    +Infinity),
!    -1: ( -Infinity,  BOUND),
!    2:  ( -Infinity, +Infinity).
!
!    Input, real ( kind = 8 ) A, B, ?
!
!    Output, real ( kind = 8 ) RESULT, the estimate of the integral I.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of the modulus of the
!    absolute error, which should not exceed |I-RESULT|.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral J.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral
!    |F-I/(B-A)| over [A,B].
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) absc1
  real ( kind = 8 ) absc2
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) boun
  real ( kind = 8 ) centr
  real ( kind = 8 ) dinf
  real ( kind = 8 ) epmach
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(7)
  real ( kind = 8 ) fv2(7)
  real ( kind = 8 ) hlgth
  integer ( kind = 4 ) inf
  integer ( kind = 4 ) j
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) tabsc1
  real ( kind = 8 ) tabsc2
  real ( kind = 8 ) uflow
  real ( kind = 8 ), dimension(8) :: wg = (/ &
    0.0000000000000000D+00,     0.1294849661688697D+00, &
    0.0000000000000000D+00,     0.2797053914892767D+00, &
    0.0000000000000000D+00,     0.3818300505051189D+00, &
    0.0000000000000000D+00,     0.4179591836734694D+00 /)
  real ( kind = 8 ), dimension(8) :: wgk = (/ &
    0.2293532201052922D-01,     0.6309209262997855D-01, &
    0.1047900103222502D+00,     0.1406532597155259D+00, &
    0.1690047266392679D+00,     0.1903505780647854D+00, &
    0.2044329400752989D+00,     0.2094821410847278D+00 /)
  real ( kind = 8 ), dimension(8) :: xgk = (/ &
    0.9914553711208126D+00,     0.9491079123427585D+00, &
    0.8648644233597691D+00,     0.7415311855993944D+00, &
    0.5860872354676911D+00,     0.4058451513773972D+00, &
    0.2077849550078985D+00,     0.0000000000000000D+00 /)

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  dinf = min ( 1, inf )
  centr = 0.5D+00 * ( a + b )
  hlgth = 0.5D+00 * ( b - a )
  tabsc1 = boun + dinf * ( 1.0D+00 - centr ) / centr
  fval1 = f(tabsc1)

  if ( inf == 2 ) then
    fval1 = fval1 + f(-tabsc1)
  end if

  fc = ( fval1 / centr ) / centr
!
!  Compute the 15-point Kronrod approximation to
!  the integral, and estimate the error.
!
  resg = wg(8) * fc
  resk = wgk(8) * fc
  resabs = abs ( resk )

  do j = 1, 7
    absc = hlgth * xgk(j)
    absc1 = centr - absc
    absc2 = centr + absc
    tabsc1 = boun + dinf * ( 1.0D+00 - absc1 ) / absc1
    tabsc2 = boun + dinf * ( 1.0D+00 -absc2 ) / absc2
    fval1 = f(tabsc1)
    fval2 = f(tabsc2)
    if ( inf == 2 ) then
      fval1 = fval1 + f(-tabsc1)
      fval2 = fval2 + f(-tabsc2)
    end if
    fval1 = ( fval1 / absc1 ) / absc1
    fval2 = ( fval2 / absc2 ) / absc2
    fv1(j) = fval1
    fv2(j) = fval2
    fsum = fval1 + fval2
    resg = resg + wg(j) * fsum
    resk = resk + wgk(j) * fsum
    resabs = resabs + wgk(j) * ( abs ( fval1 ) + abs ( fval2 ) )
  end do

  reskh = resk * 0.5D+00
  resasc = wgk(8) * abs ( fc - reskh )
  do j = 1, 7
    resasc = resasc + wgk(j) * &
      ( abs ( fv1(j) - reskh ) + abs ( fv2(j) - reskh ) )
  end do

  result = resk * hlgth
  resasc = resasc * hlgth
  resabs = resabs * hlgth
  abserr = abs ( ( resk - resg ) * hlgth )

  if ( resasc /= 0.0D+00 .and. abserr /= 0.0D+00 ) then
    abserr = resasc * min ( 1.0D+00, ( 0.2D+03 * abserr / resasc )**1.5D+00 )
  end if

  if ( uflow < ( 50.0D+00 * epmach ) * resabs ) then
    abserr = max ( abserr, ( 50.0D+00 * epmach ) * resabs )
  end if

  return
end
subroutine qpsrt ( limit, last, maxerr, ermax, elist, iord, nrmax )

!*****************************************************************************80
!
!! QPSRT maintains the descending ordering in a list of integral error estimates.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LIMIT, the maximum number of error estimates 
!    the list can contain.
!
!    Input, integer ( kind = 4 ) LAST, the number of error estimates.
!
!    Input/output, integer ( kind = 4 ) MAXERR, the index in the list of the 
!    NRMAX-th largest error.
!
!    Output, real ( kind = 8 ) ERMAX, the NRMAX-th largest error,
!    = ELIST(MAXERR).
!
!    Input, real ( kind = 8 ) ELIST(LIMIT), contains the error estimates.
!
!    Input/output, integer ( kind = 4 ) IORD(LAST).  The first K elements 
!    contain pointers to the error estimates such that ELIST(IORD(1)) through
!    ELIST(IORD(K)) form a decreasing sequence, with
!      K = LAST
!    if
!      LAST <= (LIMIT/2+2),
!    and otherwise
!      K = LIMIT+1-LAST.
!
!    Input/output, integer ( kind = 4 ) NRMAX.
!
  implicit none

  integer ( kind = 4 ) last

  real ( kind = 8 ) elist(last)
  real ( kind = 8 ) ermax
  real ( kind = 8 ) errmax
  real ( kind = 8 ) errmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) iord(last)
  integer ( kind = 4 ) isucc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jbnd
  integer ( kind = 4 ) jupbn
  integer ( kind = 4 ) k
  integer ( kind = 4 ) limit
  integer ( kind = 4 ) maxerr
  integer ( kind = 4 ) nrmax
!
!  Check whether the list contains more than two error estimates.
!
  if ( last <= 2 ) then
    iord(1) = 1
    iord(2) = 2
    go to 90
  end if
!
!  This part of the routine is only executed if, due to a difficult
!  integrand, subdivision increased the error estimate.  In the normal
!  case the insert procedure should start after the nrmax-th largest
!  error estimate.
!
  errmax = elist(maxerr)

  ido = nrmax - 1

  do i = 1, ido
    isucc = iord(nrmax-1)
    if ( errmax <= elist(isucc) ) then
      exit
    end if
    iord(nrmax) = isucc
    nrmax = nrmax - 1
  end do
!
!  Compute the number of elements in the list to
!  be maintained in descending order. this number
!  depends on the number of subdivisions still allowed.
!
  jupbn = last
  if ( ( limit / 2 + 2 ) < last ) then
    jupbn = limit + 3 - last
  end if

  errmin = elist(last)
!
!  Insert ERRMAX by traversing the list top-down,
!  starting comparison from the element elist(iord(nrmax+1)).
!
  jbnd = jupbn - 1
  ibeg = nrmax + 1

  do i = ibeg, jbnd
    isucc = iord(i)
    if ( elist(isucc) <= errmax ) then
      go to 60
    end if
    iord(i-1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
!
!  Insert ERRMIN by traversing the list bottom-up.
!
60 continue

  iord(i-1) = maxerr
  k = jbnd
  do j = i, jbnd
    isucc = iord(k)
    if ( errmin < elist(isucc) ) then
      go to 80
    end if
    iord(k+1) = isucc
    k = k-1
  end do

  iord(i) = last
  go to 90

80 continue

  iord(k+1) = last
!
!  Set maxerr and ermax.
!
90 continue

  maxerr = iord(nrmax)
  ermax = elist(maxerr)

  return
end
subroutine qraux1 ( nr, n, r, i )

!*****************************************************************************80
!
!! QRAUX1 interchanges two rows of an upper Hessenberg matrix.
!
!  Discussion:
!
!    QRAUX1 interchanges rows I and I+1 of the upper Hessenberg matrix
!    R, columns I to N.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the leading dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input/output, real ( kind = 8 ) R(NR,N), the N by N upper Hessenberg
!    matrix.
!
!    Input, integer ( kind = 4 ) I, the index of the first row to interchange.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(nr,n)

  do j = i, n
    call r8_swap ( r(i,j), r(i+1,j) )
  end do

  return
end
subroutine qraux2 ( nr, n, r, i, a, b )

!*****************************************************************************80
!
!! QRAUX2 pre-multiplies an upper Hessenberg matrix by a Jacobi rotation.
!
!  Discussion:
!
!    QRAUX2 pre-multiplies an upper Hessenberg matrix by a Jacobi rotation
!    J(I,I+1,A,B)
!
!  Modified:
!
!    15 December 2001
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) R(NR,N), the N by N upper Hessenberg
!    matrix.
!
!    Input, integer ( kind = 4 ) I, the index of the row.
!
!    Input, real ( kind = 8 ) A, B, scalars that define the rotation.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) den
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(nr,n)
  real ( kind = 8 ) s
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  den = sqrt ( a * a + b * b )
  c = a / den
  s = b / den

  do j = i, n
    y = r(i,j)
    z = r(i+1,j)
    r(i,j) = c * y - s * z
    r(i+1,j) = s * y + c * z
  end do

  return
end
subroutine qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm )

!*****************************************************************************80
!
!! QRFAC computes a QR factorization using Householder transformations.
!
!  Discussion:
!
!    This subroutine uses Householder transformations with column
!    pivoting (optional) to compute a QR factorization of the
!    M by N matrix A.  That is, QRFAC determines an orthogonal
!    matrix Q, a permutation matrix P, and an upper trapezoidal
!    matrix R with diagonal elements of nonincreasing magnitude,
!    such that A*P = Q*R.  The Householder transformation for
!    column K, K = 1,2,...,min ( M, N ), is of the form
!
!      I - ( 1 / U(K) ) * U * U'
!
!    where U has zeros in the first K-1 positions.  The form of
!    this transformation and the method of pivoting first
!    appeared in the corresponding LINPACK subroutine.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the M by N array.
!    On input, A contains the matrix for which the QR factorization is to
!    be computed.  On output, the strict upper trapezoidal part of A contains
!    the strict upper trapezoidal part of R, and the lower trapezoidal
!    part of A contains a factored form of Q (the non-trivial elements of
!    the U vectors described above).
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be no less than M.
!
!    Input, logical PIVOT, is TRUE if column pivoting is to be carried out.
!
!    Output, integer ( kind = 4 ) IPVT(LIPVT), defines the permutation matrix P 
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the identity 
!    matrix.  If PIVOT is false, IPVT is not referenced.
!
!    Input, integer ( kind = 4 ) LIPVT, the dimension of IPVT, which should be 
!    N if pivoting is used.
!
!    Output, real ( kind = 8 ) RDIAG(N), contains the diagonal elements of R.
!
!    Output, real ( kind = 8 ) ACNORM(N), the norms of the corresponding
!    columns of the input matrix A.  If this information is not needed,
!    then ACNORM can coincide with RDIAG.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) lipvt
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) acnorm(n)
  real ( kind = 8 ) ajnorm
  real ( kind = 8 ) enorm
  real ( kind = 8 ) epsmch
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(lipvt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) m
  integer ( kind = 4 ) minmn
  logical pivot
  real ( kind = 8 ) rdiag(n)
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa(n)

  epsmch = epsilon ( epsmch )
!
!  Compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = enorm ( m, a(1,j) )
  end do

  rdiag(1:n) = acnorm(1:n)
  wa(1:n) = acnorm(1:n)

  if ( pivot ) then
    do j = 1, n
      ipvt(j) = j
    end do
  end if
!
!  Reduce A to R with Householder transformations.
!
  minmn = min ( m, n )

  do j = 1, minmn
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pivot ) then

      kmax = j

      do k = j, n
        if ( rdiag(kmax) < rdiag(k) ) then
          kmax = k
        end if
      end do

      if ( kmax /= j ) then

        do i = 1, m
          call r8_swap ( a(i,j), a(i,kmax) )
        end do

        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)

        call i4_swap ( ipvt(j), ipvt(kmax) )

      end if

    end if
!
!  Compute the Householder transformation to reduce the
!  J-th column of A to a multiple of the J-th unit vector.
!
    ajnorm = enorm ( m-j+1, a(j,j) )

    if ( ajnorm /= 0.0D+00 ) then

      if ( a(j,j) < 0.0D+00 ) then
        ajnorm = -ajnorm
      end if

      a(j:m,j) = a(j:m,j) / ajnorm
      a(j,j) = a(j,j) + 1.0D+00
!
!  Apply the transformation to the remaining columns and update the norms.
!
      do k = j+1, n

        temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)

        a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

        if ( pivot .and. rdiag(k) /= 0.0D+00 ) then

          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * sqrt ( max ( 0.0D+00, 1.0D+00-temp**2 ) )

          if ( 0.05D+00 * ( rdiag(k) / wa(k) )**2 <= epsmch ) then
            rdiag(k) = enorm ( m-j, a(j+1,k) )
            wa(k) = rdiag(k)
          end if

        end if

      end do

    end if

    rdiag(j) = -ajnorm

  end do

  return
end
subroutine qrupdt ( nr, n, a, u, v )

!*****************************************************************************80
!
!! QRUPDT updates a QR factorization.
!
!  Discussion:
!
!    The routine finds an orthogonal N by N matrix Q* and an upper triangular
!    N by N matrix R* such that (Q*)(R*) = R + U*V'
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(NR,N), on input, contains the original QR
!    factorization.  On output, contains the revised factorization.
!
!    Input, real ( kind = 8 ) U(N), V(N), vectors that describe the rank
!    one update applied to the original matrix A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
!
!  Determine the last non-zero in U.
!
  k = n

  do while ( u(k) == 0.0D+00 .and. 1 < k )
    k = k - 1
  end do
!
!  (k-1) Jacobi rotations transform
!    r + u(v+) --> (r*) + ( u(1) * e1 ) (v+)
!  which is upper Hessenberg
!
  if ( 1 < k ) then

    do i = k-1, 1, -1

      if ( u(i) == 0.0D+00 ) then
        call qraux1 ( nr, n, a, i )
        u(i) = u(i+1)
      else
        call qraux2 ( nr, n, a, i, u(i), -u(i+1) )
        u(i) = sqrt ( u(i) * u(i) + u(i+1) * u(i+1) )
      end if

    end do

  end if
!
!  R <-- R + ( u(1) * e1 ) (v+)
!
  a(1,1:n) = a(1,1:n) + u(1) * v(1:n)
!
!  (k-1) Jacobi rotations transform upper Hessenberg R
!  to upper triangular (R*)
!
    do i = 1, k-1

      if ( a(i,i) == 0.0D+00 ) then
        call qraux1 ( nr, n, a, i )
      else
        t1 = a(i,i)
        t2 = -a(i+1,i)
        call qraux2 ( nr, n, a, i, t1, t2 )
      end if

    end do

  return
end
function r1mach ( i )

!*****************************************************************************80
!
!! R1MACH returns single precision real machine constants.
!
!  Discussion:
!
!    Assume that single precision real numbers are stored with a mantissa
!    of T digits in base B, with an exponent whose value must lie
!    between EMIN and EMAX.  Then for values of I between 1 and 5,
!    R1MACH will return the following values:
!
!      R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!      R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
!      R1MACH(3) = B**(-T), the smallest relative spacing.
!      R1MACH(4) = B**(1-T), the largest relative spacing.
!      R1MACH(5) = log10(B)
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer
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
!    Output, real R1MACH, the value of the chosen parameter.
!
  implicit none

  integer ( kind = 4 ) i
  real r1mach

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r1mach = 0.0E+00
    stop
  else if ( i == 1 ) then
    r1mach = 1.1754944E-38
  else if ( i == 2 ) then
    r1mach = 3.4028235E+38
  else if ( i == 3 ) then
    r1mach = 5.9604645E-08
  else if ( i == 4 ) then
    r1mach = 1.1920929E-07
  else if ( i == 5 ) then
    r1mach = 0.3010300E+00
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r1mach = 0.0E+00
    stop
  end if

  return
end
subroutine r1updt ( m, n, s, ls, u, v, w, sing )

!*****************************************************************************80
!
!! R1UPDT retriangularizes a matrix after a rank one update.
!
!  Discussion:
!
!    Given an M by N lower trapezoidal matrix S, an M-vector U, and an
!    N-vector V, the problem is to determine an orthogonal matrix Q such that
!
!      (S + U * V' ) * Q
!
!    is again lower trapezoidal.
!
!    This subroutine determines Q as the product of 2 * (N - 1)
!    transformations
!
!      GV(N-1) * ... * GV(1) * GW(1) * ... * GW(N-1)
!
!    where GV(I), GW(I) are Givens rotations in the (I,N) plane
!    which eliminate elements in the I-th and N-th planes,
!    respectively.  Q itself is not accumulated, rather the
!    information to recover the GV and GW rotations is returned.
!
!  Reference:
!
!    Jorge More, Burton Garbow and Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of S.
!
!    Input, integer ( kind = 4 ) N, the number of columns of S.  
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) S(LS).  On input, the lower trapezoidal
!    matrix S stored by columns.  On output S contains the lower trapezoidal
!    matrix produced as described above.
!
!    Input, integer ( kind = 4 ) LS, the length of the S array.  LS must be at least
!    (N*(2*M-N+1))/2.
!
!    Input, real ( kind = 8 ) U(M), the U vector.
!
!    Input/output, real ( kind = 8 ) V(N).  On input, V must contain the
!    vector V.  On output V contains the information necessary to recover
!    the Givens rotations GV described above.
!
!    Output, real ( kind = 8 ) W(M), contains information necessary to
!    recover the Givens rotations GW described above.
!
!    Output, logical SING, is set to TRUE if any of the diagonal elements
!    of the output S are zero.  Otherwise SING is set FALSE.
!
  implicit none

  integer ( kind = 4 ) ls
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) cos
  real ( kind = 8 ) cotan
  real ( kind = 8 ) giant
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) l
  real ( kind = 8 ) s(ls)
  real ( kind = 8 ) sin
  logical sing
  real ( kind = 8 ) tan
  real ( kind = 8 ) tau
  real ( kind = 8 ) temp
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) w(m)
!
!  GIANT is the largest magnitude.
!
  giant = huge ( 1.0D+00 )
!
!  Initialize the diagonal element pointer.
!
  jj = ( n * ( 2 * m - n + 1 ) ) / 2 - ( m - n )
!
!  Move the nontrivial part of the last column of S into W.
!
  l = jj
  do i = n, m
    w(i) = s(l)
    l = l + 1
  end do
!
!  Rotate the vector V into a multiple of the N-th unit vector
!  in such a way that a spike is introduced into W.
!
  do j = n-1, 1, -1

    jj = jj - ( m - j + 1 )
    w(j) = 0.0D+00

    if ( v(j) /= 0.0D+00 ) then
!
!  Determine a Givens rotation which eliminates the
!  J-th element of V.
!
      if ( abs ( v(n) ) < abs ( v(j) ) ) then
        cotan = v(n) / v(j)
        sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan**2 )
        cos = sin * cotan
        tau = 1.0D+00
        if ( 1.0D+00 < abs ( cos ) * giant ) then
          tau = 1.0D+00 / cos
        end if
      else
        tan = v(j) / v(n)
        cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan**2 )
        sin = cos * tan
        tau = sin
      end if
!
!  Apply the transformation to V and store the information
!  necessary to recover the Givens rotation.
!
      v(n) = sin * v(j) + cos * v(n)
      v(j) = tau
!
!  Apply the transformation to S and extend the spike in W.
!
      l = jj
      do i = j, m
        temp = cos * s(l) - sin * w(i)
        w(i) = sin * s(l) + cos * w(i)
        s(l) = temp
        l = l + 1
      end do

    end if

  end do
!
!  Add the spike from the rank 1 update to W.
!
   w(1:m) = w(1:m) + v(n) * u(1:m)
!
!  Eliminate the spike.
!
  sing = .false.

  do j = 1, n-1

    if ( w(j) /= 0.0D+00 ) then
!
!  Determine a Givens rotation which eliminates the
!  J-th element of the spike.
!
      if ( abs ( s(jj) ) < abs ( w (j) ) ) then
        cotan = s(jj) / w(j)
        sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan**2 )
        cos = sin * cotan
        tau = 1.0D+00
        if ( 1.0D+00 < abs ( cos ) * giant ) then
          tau = 1.0D+00 / cos
        end if
      else
        tan = w(j) / s(jj)
        cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan**2 )
        sin = cos * tan
        tau = sin
      end if
!
!  Apply the transformation to S and reduce the spike in W.
!
      l = jj
      do i = j, m
        temp = cos * s(l) + sin * w(i)
        w(i) = - sin * s(l) + cos * w(i)
        s(l) = temp
        l = l + 1
      end do
!
!  Store the information necessary to recover the Givens rotation.
!
      w(j) = tau

    end if
!
!  Test for zero diagonal elements in the output S.
!
    if ( s(jj) == 0.0D+00 ) then
      sing = .true.
    end if

    jj = jj + ( m - j + 1 )

  end do
!
!  Move W back into the last column of the output S.
!
  l = jj
  do i = n, m
    s(l) = w(i)
    l = l + 1
  end do

  if ( s(jj) == 0.0D+00 ) then
    sing = .true.
  end if

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Modified:
!
!    22 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices to print.
!    The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,g14.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_reverse ( n, a )

!*****************************************************************************80
!
!! R8VEC_REVERSE reverses the elements of an R8VEC.
!
!  Example:
!
!    Input:
!
!      N = 5, A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
!
!    Output:
!
!      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call r8_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine radb2 ( ido, l1, cc, ch, wa1 )

!*****************************************************************************80
!
!! RADB2 is a lower level routine used by RFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,2,l1)
  real ( kind = 8 ) ch(ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa1(ido)

  ch(1,1:l1,1) = cc(1,1,1:l1) + cc(ido,2,1:l1)
  ch(1,1:l1,2) = cc(1,1,1:l1) - cc(ido,2,1:l1)

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        ch(i-1,k,1) = cc(i-1,1,k) + cc(ic-1,2,k)
        tr2         = cc(i-1,1,k) - cc(ic-1,2,k)
        ch(i,k,1)   = cc(i,1,k)   - cc(ic,2,k)
        ti2         = cc(i,1,k)   + cc(ic,2,k)

        ch(i-1,k,2) = wa1(i-2) * tr2 - wa1(i-1) * ti2
        ch(i,k,2)   = wa1(i-2) * ti2 + wa1(i-1) * tr2

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  ch(ido,1:l1,1) =    cc(ido,1,1:l1) + cc(ido,1,1:l1)
  ch(ido,1:l1,2) = -( cc(1,2,1:l1)   + cc(1,2,1:l1) )

  return
end
subroutine radb3 ( ido, l1, cc, ch, wa1, wa2 )

!*****************************************************************************80
!
!! RADB3 is a lower level routine used by RFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,3,l1)
  real ( kind = 8 ) ch(ido,l1,3)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: taui =  0.866025403784439D+00
  real ( kind = 8 ), parameter :: taur = -0.5D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

  do k = 1, l1

    tr2 = cc(ido,2,k) + cc(ido,2,k)
    cr2 = cc(1,1,k) + taur * tr2
    ch(1,k,1) = cc(1,1,k) + tr2
    ci3 = taui * ( cc(1,3,k) + cc(1,3,k) )

    ch(1,k,2) = cr2 - ci3
    ch(1,k,3) = cr2 + ci3

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
      cr2 = cc(i-1,1,k) + taur * tr2
      ch(i-1,k,1) = cc(i-1,1,k) + tr2

      ti2 = cc(i,3,k) - cc(ic,2,k)
      ci2 = cc(i,1,k) + taur * ti2
      ch(i,k,1) = cc(i,1,k) + ti2

      cr3 = taui * ( cc(i-1,3,k) - cc(ic-1,2,k) )
      ci3 = taui * ( cc(i,3,k)   + cc(ic,2,k) )

      dr2 = cr2 - ci3
      dr3 = cr2 + ci3
      di2 = ci2 + cr3
      di3 = ci2 - cr3

      ch(i-1,k,2) = wa1(i-2) * dr2 - wa1(i-1) * di2
      ch(i,k,2)   = wa1(i-2) * di2 + wa1(i-1) * dr2
      ch(i-1,k,3) = wa2(i-2) * dr3 - wa2(i-1) * di3
      ch(i,k,3)   = wa2(i-2) * di3 + wa2(i-1) * dr3

    end do
  end do

  return
end
subroutine radb4 ( ido, l1, cc, ch, wa1, wa2, wa3 )

!*****************************************************************************80
!
!! RADB4 is a lower level routine used by RFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,4,l1)
  real ( kind = 8 ) ch(ido,l1,4)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: sqrt2 = 1.414213562373095D+00
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  do k = 1, l1

    tr1 = cc(1,1,k) - cc(ido,4,k)
    tr2 = cc(1,1,k) + cc(ido,4,k)
    tr3 = cc(ido,2,k) + cc(ido,2,k)
    tr4 = cc(1,3,k) + cc(1,3,k)

    ch(1,k,1) = tr2 + tr3
    ch(1,k,2) = tr1 - tr4
    ch(1,k,3) = tr2 - tr3
    ch(1,k,4) = tr1 + tr4

  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        ti1 = cc(i,1,k) + cc(ic,4,k)
        ti2 = cc(i,1,k) - cc(ic,4,k)
        ti3 = cc(i,3,k) - cc(ic,2,k)
        tr4 = cc(i,3,k) + cc(ic,2,k)

        tr1 = cc(i-1,1,k) - cc(ic-1,4,k)
        tr2 = cc(i-1,1,k) + cc(ic-1,4,k)
        ti4 = cc(i-1,3,k) - cc(ic-1,2,k)
        tr3 = cc(i-1,3,k) + cc(ic-1,2,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3         = tr2 - tr3
        ch(i,k,1)   = ti2 + ti3
        ci3         = ti2 - ti3

        cr2 = tr1 - tr4
        cr4 = tr1 + tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-2) * cr2 - wa1(i-1) * ci2
        ch(i,k,2)   = wa1(i-2) * ci2 + wa1(i-1) * cr2
        ch(i-1,k,3) = wa2(i-2) * cr3 - wa2(i-1) * ci3
        ch(i,k,3)   = wa2(i-2) * ci3 + wa2(i-1) * cr3
        ch(i-1,k,4) = wa3(i-2) * cr4 - wa3(i-1) * ci4
        ch(i,k,4)   = wa3(i-2) * ci4 + wa3(i-1) * cr4

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1

    ti1 = cc(1,2,k)   + cc(1,4,k)
    ti2 = cc(1,4,k)   - cc(1,2,k)
    tr1 = cc(ido,1,k) - cc(ido,3,k)
    tr2 = cc(ido,1,k) + cc(ido,3,k)

    ch(ido,k,1) = tr2 + tr2
    ch(ido,k,2) = sqrt2 * ( tr1 - ti1 )
    ch(ido,k,3) = ti2 + ti2
    ch(ido,k,4) = -sqrt2 * ( tr1 + ti1 )

  end do

  return
end
subroutine radb5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )

!*****************************************************************************80
!
!! RADB5 is a lower level routine used by RFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,5,l1)
  real ( kind = 8 ) ch(ido,l1,5)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: ti11 =  0.951056516295154D+00
  real ( kind = 8 ), parameter :: ti12 =  0.587785252292473D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: tr11 =  0.309016994374947D+00
  real ( kind = 8 ), parameter :: tr12 = -0.809016994374947D+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

  do k = 1, l1

    ti5 = cc(1,3,k) + cc(1,3,k)
    ti4 = cc(1,5,k) + cc(1,5,k)
    tr2 = cc(ido,2,k) + cc(ido,2,k)
    tr3 = cc(ido,4,k) + cc(ido,4,k)

    ch(1,k,1) = cc(1,1,k) + tr2 + tr3
    cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
    cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
    ci5 = ti11 * ti5 + ti12 * ti4
    ci4 = ti12 * ti5 - ti11 * ti4

    ch(1,k,2) = cr2 - ci5
    ch(1,k,3) = cr3 - ci4
    ch(1,k,4) = cr3 + ci4
    ch(1,k,5) = cr2 + ci5

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      ti5 = cc(i,3,k) + cc(ic,2,k)
      ti2 = cc(i,3,k) - cc(ic,2,k)
      ti4 = cc(i,5,k) + cc(ic,4,k)
      ti3 = cc(i,5,k) - cc(ic,4,k)
      tr5 = cc(i-1,3,k) - cc(ic-1,2,k)
      tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
      tr4 = cc(i-1,5,k) - cc(ic-1,4,k)
      tr3 = cc(i-1,5,k) + cc(ic-1,4,k)

      ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
      ch(i,k,1)   = cc(i,1,k) + ti2 + ti3

      cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
      cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      dr3 = cr3 - ci4
      dr4 = cr3 + ci4
      di3 = ci3 + cr4
      di4 = ci3 - cr4
      dr5 = cr2 + ci5
      dr2 = cr2 - ci5
      di5 = ci2 - cr5
      di2 = ci2 + cr5

      ch(i-1,k,2) = wa1(i-2) * dr2 - wa1(i-1) * di2
      ch(i,k,2)   = wa1(i-2) * di2 + wa1(i-1) * dr2
      ch(i-1,k,3) = wa2(i-2) * dr3 - wa2(i-1) * di3
      ch(i,k,3)   = wa2(i-2) * di3 + wa2(i-1) * dr3
      ch(i-1,k,4) = wa3(i-2) * dr4 - wa3(i-1) * di4
      ch(i,k,4)   = wa3(i-2) * di4 + wa3(i-1) * dr4
      ch(i-1,k,5) = wa4(i-2) * dr5 - wa4(i-1) * di5
      ch(i,k,5)   = wa4(i-2) * di5 + wa4(i-1) * dr5

    end do
  end do

  return
end
subroutine radbg ( ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )

!*****************************************************************************80
!
!! RADBG is a lower level routine used by RFFTB1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ai1
  real ( kind = 8 ) ai2
  real ( kind = 8 ) ar1
  real ( kind = 8 ) ar1h
  real ( kind = 8 ) ar2
  real ( kind = 8 ) ar2h
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1(ido,l1,ip)
  real ( kind = 8 ) c2(idl1,ip)
  real ( kind = 8 ) cc(ido,ip,l1)
  real ( kind = 8 ) ch(ido,l1,ip)
  real ( kind = 8 ) ch2(idl1,ip)
  real ( kind = 8 ) dc2
  real ( kind = 8 ) dcp
  real ( kind = 8 ) ds2
  real ( kind = 8 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) wa(*)

  arg = 2.0D+00 * pi / real ( ip, kind = 8 )
  dcp = cos ( arg )
  dsp = sin ( arg )
  nbd = ( ido - 1 ) / 2
  ipph = ( ip + 1 ) / 2
  ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  do j = 2, ipph
    jc = ip + 2 - j
    j2 = j + j
    ch(1,1:l1,j) =  cc(ido,j2-2,1:l1) + cc(ido,j2-2,1:l1)
    ch(1,1:l1,jc) = cc(1,j2-1,1:l1)   + cc(1,j2-1,1:l1)
  end do

  if ( ido /= 1 ) then

    if ( l1 <= nbd ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            ic = ido + 2 - i
            ch(i-1,k,j)  = cc(i-1,2*j-1,k) + cc(ic-1,2*j-2,k)
            ch(i-1,k,jc) = cc(i-1,2*j-1,k) - cc(ic-1,2*j-2,k)
            ch(i,k,j)    = cc(i,2*j-1,k)   - cc(ic,2*j-2,k)
            ch(i,k,jc)   = cc(i,2*j-1,k)   + cc(ic,2*j-2,k)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          ic = ido + 2 - i
          ch(i-1,1:l1,j)  = cc(i-1,2*j-1,1:l1) + cc(ic-1,2*j-2,1:l1)
          ch(i-1,1:l1,jc) = cc(i-1,2*j-1,1:l1) - cc(ic-1,2*j-2,1:l1)
          ch(i,1:l1,j)    = cc(i,2*j-1,1:l1)   - cc(ic,2*j-2,1:l1)
          ch(i,1:l1,jc)   = cc(i,2*j-1,1:l1)   + cc(ic,2*j-2,1:l1)
        end do
      end do

    end if

  end if

  ar1 = 1.0D+00
  ai1 = 0.0D+00

  do l = 2, ipph

    lc = ip + 2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      c2(ik,l)  = ch2(ik,1) + ar1 * ch2(ik,2)
      c2(ik,lc) =             ai1 * ch2(ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ip + 2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2  = dc2 * ai2 + ds2 * ar2
      ar2 = ar2h

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + ar2 * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) + ai2 * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    ch(1,1:l1,j)  = c1(1,1:l1,j) - c1(1,1:l1,jc)
    ch(1,1:l1,jc) = c1(1,1:l1,j) + c1(1,1:l1,jc)
  end do

  if ( ido /= 1 ) then

    if ( l1 <= nbd ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            ch(i-1,k,j)  = c1(i-1,k,j) - c1(i,k,jc)
            ch(i-1,k,jc) = c1(i-1,k,j) + c1(i,k,jc)
            ch(i,k,j)    = c1(i,k,j)   + c1(i-1,k,jc)
            ch(i,k,jc)   = c1(i,k,j)   - c1(i-1,k,jc)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          ch(i-1,1:l1,j)  = c1(i-1,1:l1,j) - c1(i,1:l1,jc)
          ch(i-1,1:l1,jc) = c1(i-1,1:l1,j) + c1(i,1:l1,jc)
          ch(i,1:l1,j)    = c1(i,1:l1,j)   + c1(i-1,1:l1,jc)
          ch(i,1:l1,jc)   = c1(i,1:l1,j)   - c1(i-1,1:l1,jc)
        end do
      end do

    end if

  end if

  if ( ido == 1 ) then
    return
  end if

  c2(1:idl1,1) = ch2(1:idl1,1)
  c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)

  if ( nbd <= l1 ) then

    is = -ido

    do j = 2, ip
      is = is + ido
      idij = is
      do i = 3, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) - wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   + wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    is = -ido
    do j = 2, ip
      is = is + ido
      do k = 1, l1
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) - wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   + wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine radf2 ( ido, l1, cc, ch, wa1 )

!*****************************************************************************80
!
!! RADF2 is a lower level routine used by RFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,l1,2)
  real ( kind = 8 ) ch(ido,2,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa1(ido)

  ch(1,1,1:l1)   = cc(1,1:l1,1) + cc(1,1:l1,2)
  ch(ido,2,1:l1) = cc(1,1:l1,1) - cc(1,1:l1,2)

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        tr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ti2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)

        ch(i,1,k) = cc(i,k,1) + ti2
        ch(ic,2,k) = ti2 - cc(i,k,1)
        ch(i-1,1,k) = cc(i-1,k,1) + tr2
        ch(ic-1,2,k) = cc(i-1,k,1) - tr2

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  ch(1,2,1:l1) = -cc(ido,1:l1,2)
  ch(ido,1,1:l1) = cc(ido,1:l1,1)

  return
end
subroutine radf3 ( ido, l1, cc, ch, wa1, wa2 )

!*****************************************************************************80
!
!! RADF3 is a lower level routine used by RFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,l1,3)
  real ( kind = 8 ) ch(ido,3,l1)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) cr2
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: taui = 0.866025403784439D+00
  real ( kind = 8 ), parameter :: taur = -0.5D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

  do k = 1, l1
    cr2 = cc(1,k,2) + cc(1,k,3)
    ch(1,1,k) = cc(1,k,1) + cr2
    ch(1,3,k) = taui * ( cc(1,k,3) - cc(1,k,2) )
    ch(ido,2,k) = cc(1,k,1) + taur * cr2
  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)

      cr2 = dr2 + dr3
      ci2 = di2 + di3

      ch(i-1,1,k) = cc(i-1,k,1) + cr2
      ch(i,1,k)   = cc(i,k,1) + ci2

      tr2 = cc(i-1,k,1) + taur * cr2
      ti2 = cc(i,k,1) + taur * ci2
      tr3 = taui * ( di2 - di3 )
      ti3 = taui * ( dr3 - dr2 )

      ch(i-1,3,k) = tr2 + tr3
      ch(ic-1,2,k) = tr2 - tr3
      ch(i,3,k) = ti2 + ti3
      ch(ic,2,k) = ti3 - ti2

    end do
  end do

  return
end
subroutine radf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )

!*****************************************************************************80
!
!! RADF4 is a lower level routine used by RFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,l1,4)
  real ( kind = 8 ) ch(ido,4,l1)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ), parameter :: hsqt2 = 0.7071067811865475D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  do k = 1, l1
    tr1 = cc(1,k,2) + cc(1,k,4)
    tr2 = cc(1,k,1) + cc(1,k,3)
    ch(1,1,k) = tr1 + tr2
    ch(ido,4,k) = tr2 - tr1
    ch(ido,2,k) = cc(1,k,1) - cc(1,k,3)
    ch(1,3,k) = cc(1,k,4) - cc(1,k,2)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        cr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ci2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
        cr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
        ci3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
        cr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
        ci4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)

        tr1 = cr2+cr4
        tr4 = cr4-cr2
        ti1 = ci2+ci4
        ti4 = ci2-ci4
        ti2 = cc(i,k,1) + ci3
        ti3 = cc(i,k,1) - ci3
        tr2 = cc(i-1,k,1) + cr3
        tr3 = cc(i-1,k,1) - cr3

        ch(i-1,1,k)  = tr1 + tr2
        ch(ic-1,4,k) = tr2 - tr1
        ch(i,1,k)    = ti1 + ti2
        ch(ic,4,k)   = ti1 - ti2
        ch(i-1,3,k)  = ti4 + tr3
        ch(ic-1,2,k) = tr3 - ti4
        ch(i,3,k)    = tr4 + ti3
        ch(ic,2,k)   = tr4 - ti3

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1

    ti1 = -hsqt2 * ( cc(ido,k,2) + cc(ido,k,4) )
    tr1 =  hsqt2 * ( cc(ido,k,2) - cc(ido,k,4) )

    ch(ido,1,k) = tr1 + cc(ido,k,1)
    ch(ido,3,k) = cc(ido,k,1) - tr1

    ch(1,2,k) = ti1 - cc(ido,k,3)
    ch(1,4,k) = ti1 + cc(ido,k,3)

  end do

  return
end
subroutine radf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )

!*****************************************************************************80
!
!! RADF5 is a lower level routine used by RFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(ido,l1,5)
  real ( kind = 8 ) ch(ido,5,l1)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: ti11 =  0.951056516295154D+00
  real ( kind = 8 ), parameter :: ti12 =  0.587785252292473D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: tr11 =  0.309016994374947D+00
  real ( kind = 8 ), parameter :: tr12 = -0.809016994374947D+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

  do k = 1, l1

    cr2 = cc(1,k,5) + cc(1,k,2)
    ci5 = cc(1,k,5) - cc(1,k,2)
    cr3 = cc(1,k,4) + cc(1,k,3)
    ci4 = cc(1,k,4) - cc(1,k,3)

    ch(1,1,k)   = cc(1,k,1) + cr2 + cr3
    ch(ido,2,k) = cc(1,k,1) + tr11 * cr2 + tr12 * cr3
    ch(1,3,k)   = ti11 * ci5 + ti12 * ci4
    ch(ido,4,k) = cc(1,k,1) + tr12 * cr2 + tr11 * cr3
    ch(1,5,k)   = ti12 * ci5 - ti11 * ci4

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
      dr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
      di4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)
      dr5 = wa4(i-2) * cc(i-1,k,5) + wa4(i-1) * cc(i,k,5)
      di5 = wa4(i-2) * cc(i,k,5)   - wa4(i-1) * cc(i-1,k,5)

      cr2 = dr2 + dr5
      ci5 = dr5 - dr2
      cr5 = di2 - di5
      ci2 = di2 + di5
      cr3 = dr3 + dr4
      ci4 = dr4 - dr3
      cr4 = di3 - di4
      ci3 = di3 + di4

      ch(i-1,1,k) = cc(i-1,k,1) + cr2 + cr3
      ch(i,1,k)   = cc(i,k,1)   + ci2 + ci3

      tr2 = cc(i-1,k,1) + tr11 * cr2 + tr12 * cr3
      ti2 = cc(i,k,1)   + tr11 * ci2 + tr12 * ci3
      tr3 = cc(i-1,k,1) + tr12 * cr2 + tr11 * cr3
      ti3 = cc(i,k,1)   + tr12 * ci2 + tr11 * ci3

      tr5 = ti11 * cr5 + ti12 * cr4
      ti5 = ti11 * ci5 + ti12 * ci4
      tr4 = ti12 * cr5 - ti11 * cr4
      ti4 = ti12 * ci5 - ti11 * ci4

      ch(i-1,3,k)  = tr2 + tr5
      ch(ic-1,2,k) = tr2 - tr5
      ch(i,3,k)    = ti2 + ti5
      ch(ic,2,k)   = ti5 - ti2
      ch(i-1,5,k)  = tr3 + tr4
      ch(ic-1,4,k) = tr3 - tr4
      ch(i,5,k)    = ti3 + ti4
      ch(ic,4,k)   = ti4 - ti3

    end do
  end do

  return
end
subroutine radfg ( ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )

!*****************************************************************************80
!
!! RADFG is a lower level routine used by RFFTF1.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ai1
  real ( kind = 8 ) ai2
  real ( kind = 8 ) ar1
  real ( kind = 8 ) ar1h
  real ( kind = 8 ) ar2
  real ( kind = 8 ) ar2h
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1(ido,l1,ip)
  real ( kind = 8 ) c2(idl1,ip)
  real ( kind = 8 ) cc(ido,ip,l1)
  real ( kind = 8 ) ch(ido,l1,ip)
  real ( kind = 8 ) ch2(idl1,ip)
  real ( kind = 8 ) dc2
  real ( kind = 8 ) dcp
  real ( kind = 8 ) ds2
  real ( kind = 8 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) wa(*)

  arg = 2.0D+00 * pi / real ( ip, kind = 8 )
  dcp = cos ( arg )
  dsp = sin ( arg )
  ipph = ( ip + 1 ) / 2
  nbd = ( ido - 1 ) / 2

  if ( ido == 1 ) then

    c2(1:idl1,1) = ch2(1:idl1,1)

  else

    ch2(1:idl1,1) = c2(1:idl1,1)
    ch(1,1:l1,2:ip) = c1(1,1:l1,2:ip)

    if ( nbd <= l1 ) then

      is = -ido
      do j = 2, ip
        is = is + ido
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          do k = 1, l1
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    else

      is = -ido

      do j = 2, ip
        is = is + ido
        do k = 1, l1
          idij = is
          do i = 3, ido, 2
            idij = idij + 2
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    end if

    if ( l1 <= nbd ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            c1(i-1,k,j)  = ch(i-1,k,j)  + ch(i-1,k,jc)
            c1(i-1,k,jc) = ch(i,k,j)    - ch(i,k,jc)
            c1(i,k,j)    = ch(i,k,j)    + ch(i,k,jc)
            c1(i,k,jc)   = ch(i-1,k,jc) - ch(i-1,k,j)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          c1(i-1,1:l1,j)  = ch(i-1,1:l1,j)  + ch(i-1,1:l1,jc)
          c1(i-1,1:l1,jc) = ch(i,1:l1,j)    - ch(i,1:l1,jc)
          c1(i,1:l1,j)    = ch(i,1:l1,j)    + ch(i,1:l1,jc)
          c1(i,1:l1,jc)   = ch(i-1,1:l1,jc) - ch(i-1,1:l1,j)
        end do
      end do

    end if

  end if

  do j = 2, ipph
    jc = ip + 2 - j
    c1(1,1:l1,j)  = ch(1,1:l1,j)  + ch(1,1:l1,jc)
    c1(1,1:l1,jc) = ch(1,1:l1,jc) - ch(1,1:l1,j)
  end do

  ar1 = 1.0D+00
  ai1 = 0.0D+00

  do l = 2, ipph

    lc = ip + 2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      ch2(ik,l) = c2(ik,1) + ar1 * c2(ik,2)
      ch2(ik,lc) =           ai1 * c2(ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ip + 2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2 =  dc2 * ai2 + ds2 * ar2
      ar2 = ar2h

      do ik = 1, idl1
        ch2(ik,l) =  ch2(ik,l)  + ar2 * c2(ik,j)
        ch2(ik,lc) = ch2(ik,lc) + ai2 * c2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + c2(1:idl1,j)
  end do

  cc(1:ido,1,1:l1) = ch(1:ido,1:l1,1)

  do j = 2, ipph
    jc = ip + 2 - j
    j2 = j + j
    cc(ido,j2-2,1:l1) = ch(1,1:l1,j)
    cc(1,j2-1,1:l1)   = ch(1,1:l1,jc)
  end do

  if ( ido == 1 ) then
    return
  end if

  if ( l1 <= nbd ) then

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do k = 1, l1
        do i = 3, ido, 2
          ic = ido + 2 - i
          cc(i-1,j2-1,k)  = ch(i-1,k,j) + ch(i-1,k,jc)
          cc(ic-1,j2-2,k) = ch(i-1,k,j) - ch(i-1,k,jc)
          cc(i,j2-1,k)    = ch(i,k,j)   + ch(i,k,jc)
          cc(ic,j2-2,k)   = ch(i,k,jc)  - ch(i,k,j)
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do i = 3, ido, 2
        ic = ido + 2 - i
        cc(i-1,j2-1,1:l1)  = ch(i-1,1:l1,j) + ch(i-1,1:l1,jc)
        cc(ic-1,j2-2,1:l1) = ch(i-1,1:l1,j) - ch(i-1,1:l1,jc)
        cc(i,j2-1,1:l1)    = ch(i,1:l1,j)   + ch(i,1:l1,jc)
        cc(ic,j2-2,1:l1)   = ch(i,1:l1,jc)  - ch(i,1:l1,j)
      end do
    end do

  end if

  return
end
subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a seed value.
!
  implicit none

  integer ( kind = 4 ) date_time(8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 8 ) t
!
!  Initialize the random number seed.
!
  call random_seed ( )
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )
!
!  Get the current date and time.
!
  call date_and_time ( values = date_time )
!
!  Construct a slightly random value.
!
  seed = 0
  do i = 1, 8
    seed = ieor ( seed, date_time(i) )
  end do
!
!  Make slightly random assignments to SEED_VECTOR.
!
  do i = 1, seed_size
    seed_vector(i) = ieor ( seed, i )
  end do
!
!  Set the random number seed value.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )

  do i = 1, 100
    call random_number ( harvest = t )
  end do

  return
end
subroutine result ( nr, n, x, f, g, a, p, itncnt, iflg, ipr )

!*****************************************************************************80
!
!! RESULT prints information about the optimization process.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the current iterate.
!
!    Input, real ( kind = 8 ) F, the function value at X.
!
!    Input, real ( kind = 8 ) G(N), the gradient at X.
!
!    Input, real ( kind = 8 ) A(NR,N), the N by N Hessian matrix at X.
!
!    Input, real ( kind = 8 ) P(N), the step taken.
!
!    Input, integer ( kind = 4 ) ITNCNT, the iteration number.
!
!    Input, integer ( kind = 4 ) IFLG, the flag controlling the amount of printout.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflg
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) itncnt
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) x(n)

  write ( ipr, 903 ) itncnt

  if ( iflg /= 0 ) then
    write ( ipr, * ) ' result       step'
    write ( ipr,905) p(1:n)
  end if

  write ( ipr, * ) ' result       x(k)'
  write ( ipr, 905) x(1:n)
  write ( ipr, * ) ' result     function at x(k)'
  write ( ipr, 905) f
  write ( ipr, * ) ' result       gradient at x(k)'
  write ( ipr, 905) g(1:n)

  if ( iflg /= 0 ) then

    write ( ipr, * ) ' result       Hessian at x(k)'
    do i = 1, n
      write ( ipr, 900) i
      write ( ipr, 902) a(i,1:i)
    end do

  end if

  return

  900 format(' result     row',i5)
  902 format(' result       ',5(2x,e20.13))
  903 format(/'0result    iterate k=',i5)
  905 format(' result               ',5(2x,e20.13) )
end
function runge ( x )

!*****************************************************************************80
!
!! RUNGE evaluates Runge's function.
!
!  Discussion:
!
!    Runge's function is a common test case for interpolation
!    and approximation, over the interval [-1,1].
!
!    RUNGE(X) = 1 / ( 1 + 25 * X**2 )
!
!  Modified:
!
!    11 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of Runge's function.
!
!    Output, real ( kind = 8 ) RUNGE, the value of Runge's function.
!
  implicit none

  real ( kind = 8 ) runge
  real ( kind = 8 ) x

  runge = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * x * x )

  return
end
function rungep ( x )

!*****************************************************************************80
!
!! RUNGEP evaluates the derivative of Runge's function.
!
!  Discussion:
!
!    Runge's function is a common test case for interpolation
!    and approximation, over the interval [-1,1].
!
!    RUNGE(X) = 1 / ( 1 + 25 * X**2 )
!
!  Modified:
!
!    11 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) RUNGEP, the derivative of Runge's function.
!
  implicit none

  real ( kind = 8 ) rungep
  real ( kind = 8 ) x

  rungep = ( -50.0D+00 * x ) / ( 1.0D+00 + 25.0D+00 * x * x )**2

  return
end
subroutine secfac ( nr, n, x, g, a, xpls, gpls, epsm, itncnt, rnf, &
  iagflg, noupdt, s, y, u, w )

!*****************************************************************************80
!
!! SECFAC updates the hessian by the BFGS factored method.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, X[K-1].
!
!    Input, real ( kind = 8 ) G(N), the gradient or an approximation,
!    at the old iterate.
!
!    Input/output, real ( kind = 8 ) A(NR,N).
!    On input, the Cholesky decomposition of hessian in lower part and diagonal.
!    On output, the updated Cholesky decomposition of hessian
!    in lower triangular part and diagonal
!
!    Input, real ( kind = 8 ) XPLS(N), the new iterate, X[K].
!
!    Input, real ( kind = 8 ) GPLSN(N), gradient, or an approximation,
!    at the new iterate.
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input, integer ( kind = 4 ) ITNCNT, the iteration count.
!
!    Input, real ( kind = 8 ) RNF, the relative noise in the optimization
!    function FCN.
!
!    Input, integer ( kind = 4 ) IAGFLG, 1 if analytic gradient supplied.
!
!    Input/output, logical NOUPDT, is TRUE if there has been no update
!    yet.  The user should retain the output value between successive
!    calls.
!
!    Workspace, real ( kind = 8 ) S(N).
!
!    Workspace, real ( kind = 8 ) Y(N).
!
!    Workspace, real ( kind = 8 ) U(N).
!
!    Workspace, real ( kind = 8 ) W(N).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) alp
  real ( kind = 8 ) den1
  real ( kind = 8 ) den2
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) epsm
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) gpls(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iagflg
  integer ( kind = 4 ) itncnt
  integer ( kind = 4 ) j
  logical noupdt
  real ( kind = 8 ) reltol
  real ( kind = 8 ) rnf
  real ( kind = 8 ) s(n)
  logical skpupd
  real ( kind = 8 ) snorm2
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xpls(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ynrm2

  if ( itncnt == 1 ) then
    noupdt = .true.
  end if

  s(1:n) = xpls(1:n) - x(1:n)
  y(1:n) = gpls(1:n) - g(1:n)

  den1 = dot_product ( s, y )

  snorm2 = dnrm2 ( n, s, 1)

  ynrm2 = dnrm2 ( n, y, 1)

  if ( den1 < sqrt ( epsm ) * snorm2 * ynrm2 ) then
    return
  end if

  call mvmltu ( nr, n, a, s, u )

  den2 = dot_product ( u, u )
!
!  L <-- sqrt ( den1 / den2 ) * L
!
  alp = sqrt ( den1 / den2 )

  if ( noupdt ) then

    u(1:n) = alp * u(1:n)

    do j = 1, n
      do i = j, n
        a(i,j) = alp * a(i,j)
      end do
    end do

    noupdt = .false.
    den2 = den1
    alp = 1.0D+00

  end if

  skpupd = .true.
!
!  W = l(l+)s = hs
!
  call mvmltl ( nr, n, a, u, w )
  i = 1

  if ( iagflg == 0 ) then
    reltol = sqrt ( rnf )
  else
    reltol = rnf
  end if

60  continue

  if ( i <= n .and. skpupd ) then

    if ( abs ( y(i) - w(i) ) < reltol * &
      max ( abs ( g(i) ), abs ( gpls(i) ) ) ) then
      i = i + 1
    else
      skpupd = .false.
    end if
    go to 60
  end if

  if ( skpupd ) then
    return
  end if
!
!  W = y-alp*l(l+)s
!
  w(1:n) = y(1:n) - alp * w(1:n)
!
!  ALP = 1 / sqrt ( den1 * den2 )
!
  alp = alp / den1
!
!  U = (l+) / sqrt ( den1 * den2 ) = (l+)s/ sqrt ( ( y+ ) s * (s+) l (l+) s )
!
  u(1:n) = alp * u(1:n)
!
!  Copy L into upper triangular part.  Zero L.
!
  do i = 2, n
    do j = 1, i-1
      a(j,i) = a(i,j)
      a(i,j) = 0.0D+00
    end do
  end do
!
!  Find Q, (l+) such that  q(l+) = (l+) + u(w+)
!
  call qrupdt ( nr, n, a, u, w )
!
!  Upper triangular part and diagonal of a now contain updated
!  Cholesky decomposition of hessian.  Copy back to lower triangular part.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = a(j,i)
    end do
  end do

  return
end
subroutine secunf ( nr, n, x, g, a, udiag, xpls, gpls, epsm, itncnt, &
  rnf, iagflg, noupdt, s, y, t )

!*****************************************************************************80
!
!! SECUNF updates a Hessian matrix by the BFGS unfactored method.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, X[K-1].
!
!    Input, real ( kind = 8 ) G(N), the gradient, or an approximate value,
!    at the  old iterate.
!
!    Input/output, real ( kind = 8 ) A(NR,N).
!    on entry: approximate hessian at old iterate
!    in upper triangular part (and udiag)
!    on exit:  updated approx hessian at new iterate
!    in lower triangular part and diagonal
!    [lower triangular part of symmetric matrix]
!
!    Input, real ( kind = 8 ) UDIAG(N), the diagonal entries of the hessian.
!
!    Input, real ( kind = 8 ) XPLS(N), the new iterate, X[K].
!
!    Input, real ( kind = 8 ) GPLS(N), the gradient or an approximate value, at
!    the new iterate
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input, integer ( kind = 4 ) ITNCNT, the iteration count.
!
!    Input, real ( kind = 8 ) RNF, the relative noise in the optimization
!    function.
!
!    Input, integer ( kind = 4 ) IAGFLG, =1 if analytic gradient supplied, =0 otherwise
!
!    Input/output, logical NOUPDT, TRUE if no update yet.
!    [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) S(N).
!
!    Workspace, real ( kind = 8 ) Y(N).
!
!    Workspace, real ( kind = 8 ) T(N).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) den1
  real ( kind = 8 ) den2
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) epsm
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) gam
  real ( kind = 8 ) gpls(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iagflg
  integer ( kind = 4 ) itncnt
  integer ( kind = 4 ) j
  logical noupdt
  real ( kind = 8 ) rnf
  real ( kind = 8 ) s(n)
  logical skpupd
  real ( kind = 8 ) snorm2
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) tol
  real ( kind = 8 ) udiag(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xpls(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ynrm2
!
!  Copy hessian in upper triangular part and UDIAG to
!  lower triangular part and diagonal.
!
  do j = 1, n
    a(j,j) = udiag(j)
    do i = j+1, n
      a(i,j) = a(j,i)
    end do
  end do

  if ( itncnt == 1 ) then
    noupdt = .true.
  end if

  s(1:n) = xpls(1:n) - x(1:n)
  y(1:n) = gpls(1:n) - g(1:n)

  den1 = dot_product ( s, y )

  snorm2 = dnrm2 ( n, s, 1 )

  ynrm2 = dnrm2 ( n, y, 1 )

  if ( den1 < sqrt ( epsm ) * snorm2 * ynrm2 ) then
    return
  end if

  call mvmlts ( nr, n, a, s, t )

  den2 = dot_product ( s, t )

  if ( noupdt ) then
!
!  H <-- [(s+)y/(s+)hs]h
!
    gam = den1 / den2
    den2 = gam * den2
    do j = 1, n
      t(j) = gam * t(j)
      do i = j, n
        a(i,j) = gam * a(i,j)
      end do
    end do
    noupdt = .false.

  end if

  skpupd = .true.
!
!  Check update condition on row I.
!
  do i = 1, n

    tol = rnf * max ( abs ( g(i) ), abs ( gpls(i) ) )
    if ( iagflg == 0 ) then
      tol = tol / sqrt ( rnf )
    end if

    if ( tol <= abs ( y(i) - t(i) ) ) then
      skpupd = .false.
      exit
    end if

  end do

  if ( skpupd ) then
    return
  end if
!
!  BFGS update
!
  do j = 1, n
    do i = j, n
      a(i,j) = a(i,j) + y(i) * y(j) / den1 - t(i) * t(j) / den2
    end do
  end do

  return
end
subroutine sinqb ( n, x, wsave )

!*****************************************************************************80
!
!! SINQB computes the fast sine transform of quarter wave data.
!
!  Discussion:
!
!    SINQB computes a sequence from its representation in terms of a sine
!    series with odd wave numbers.
!
!    SINQF is the unnormalized inverse of SINQB since a call of SINQB
!    followed by a call of SINQF will multiply the input sequence X by 4*N.
!
!    The array WSAVE must be initialized by calling SINQI.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!
!        4 * X_in(K) * sin ( ( 2 * K - 1 ) * I * PI / ( 2 * N ) )
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real ( kind = 8 ) WSAVE(3*N+15), a work array.  The WSAVE array
!    must be initialized by calling SINQI.  A different WSAVE array must be
!    used for each different value of N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) wsave(3*n+15)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    x(1) = 4.0D+00 * x(1)
    return
  end if

  x(2:n:2) = -x(2:n:2)

  call cosqb ( n, x, wsave )
!
!  Reverse the X vector.
!
  call r8vec_reverse ( n, x )

  return
end
subroutine sinqf ( n, x, wsave )

!*****************************************************************************80
!
!! SINQF computes the fast sine transform of quarter wave data.
!
!  Discussion:
!
!    SINQF computes the coefficients in a sine series representation with
!    only odd wave numbers.
!
!    SINQB is the unnormalized inverse of SINQF since a call of SINQF
!    followed by a call of SINQB will multiply the input sequence X by 4*N.
!
!    The array WSAVE, which is used by SINQF, must be initialized by
!    calling SINQI.
!
!    The transform is defined by:
!
!      X_out(I) = (-1)**(I-1) * X_in(N) + sum ( 1 <= K <= N-1 )
!        2 * X_in(K) * sin ( ( 2 * I - 1 ) * K * PI / ( 2 * N ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real ( kind = 8 ) WSAVE(3*N+15), a work array.  The WSAVE array
!    must be initialized by calling SINQI.  A different WSAVE array must be
!    used for each different value of N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) wsave(3*n+15)
  real ( kind = 8 ) x(n)

  if ( n <= 1 ) then
    return
  end if
!
!  Reverse the X vector.
!
  call r8vec_reverse ( n, x )

  call cosqf ( n, x, wsave )

  x(2:n:2) = -x(2:n:2)

  return
end
subroutine sinqi ( n, wsave )

!*****************************************************************************80
!
!! SINQI initializes WSAVE, used in SINQF and SINQB.
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the array to be transformed.
!
!    Output, real ( kind = 8 ) WSAVE(3*N+15), contains data, dependent on
!    the value of N, which is necessary for the SINQF or SINQB routines.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) wsave(3*n+15)

  call cosqi ( n, wsave )

  return
end
subroutine sint ( n, x, wsave )

!*****************************************************************************80
!
!! SINT computes the discrete Fourier sine transform of an odd sequence.
!
!  Discussion:
!
!    SINT is the unnormalized inverse of itself since a call of SINT
!    followed by another call of SINT will multiply the input sequence
!    X by 2*(N+1).
!
!    The array WSAVE must be initialized by calling SINTI.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!        2 * X_in(K) * sin ( K * I * PI / ( N + 1 ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be transformed.
!    The method is most efficient when N+1 is the product of small primes.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real ( kind = 8 ) WSAVE((5*N+30)/2), a work array.  The WSAVE
!    array must be initialized by calling SINTI.  A different WSAVE array
!    must be used for each different value of N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) iw2
  integer ( kind = 4 ) iw3
  real ( kind = 8 ) wsave((5*n+30)/2)
  real ( kind = 8 ) x(n)

  iw1 = ( n / 2 ) + 1
  iw2 = iw1 + n + 1
  iw3 = iw2 + n + 1

  call sint1 ( n, x, wsave(1), wsave(iw1), wsave(iw2), wsave(iw3) )

  return
end
subroutine sint1 ( n, war, was, xh, x, ifac )

!*****************************************************************************80
!
!! SINT1 is a lower level routine used by SINT.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be transformed.
!
!    Input/output, real ( kind = 8 ) WAR(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real ( kind = 8 ) WAS(N/2).
!
!    Input, real ( kind = 8 ) XH(N).
!
!    Input, real ( kind = 8 ) X(N+1).
!
!    Input, integer ( kind = 4 ) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifac(15)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ns2
  real ( kind = 8 ), parameter :: sqrt3 = 1.73205080756888D+00
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) war(n)
  real ( kind = 8 ) was(n/2)
  real ( kind = 8 ) x(n+1)
  real ( kind = 8 ) xh(n)
  real ( kind = 8 ) xhold

  xh(1:n) = war(1:n)
  war(1:n) = x(1:n)

  if ( n <= 1 ) then
    xh(1) = 2.0D+00 * xh(1)
    return
  end if

  if ( n == 2 ) then
    xhold = sqrt3 * ( xh(1) + xh(2) )
    xh(2) = sqrt3 * ( xh(1) - xh(2) )
    xh(1) = xhold
    return
  end if

  ns2 = n / 2
  x(1) = 0.0D+00

  do k = 1, n/2
    t1 = xh(k) - xh(n+1-k)
    t2 = was(k) * ( xh(k) + xh(n+1-k) )
    x(k+1) = t1 + t2
    x(n+2-k) = t2 - t1
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n/2+2) = 4.0D+00 * xh(n/2+1)
  end if

  call dfftf1 ( n+1, x, xh, war, ifac )

  xh(1) = 0.5D+00 * x(1)
  do i = 3, n, 2
    xh(i-1) = -x(i)
    xh(i) = xh(i-2) + x(i-1)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    xh(n) = -x(n+1)
  end if

  x(1:n) = war(1:n)
  war(1:n) = xh(1:n)

  return
end
subroutine sinti ( n, wsave )

!*****************************************************************************80
!
!! SINTI initializes WSAVE, used in SINT.
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be transformed.
!    The method is most efficient when N+1 is a product of small primes.
!
!    Output, real ( kind = 8 ) WSAVE((5*N+30)/2), contains data, dependent
!    on the value of N, which is necessary for the SINT routine.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) dt
  integer ( kind = 4 ) k
  real ( kind = 8 ) wsave((5*n+30)/2)

  if ( n <= 1 ) then
    return
  end if

  dt = pi / real ( n + 1, kind = 8 )

  do k = 1, n/2
    wsave(k) = 2.0D+00 * sin ( real ( k, kind = 8 ) * dt )
  end do

  call dffti ( n+1, wsave((n/2)+1) )

  return
end
subroutine sndofd ( nr, n, xpls, fcn, fpls, a, sx, rnoise, stepsz, anbr )

!*****************************************************************************80
!
!! SNDOFD approximates a Hessian with a second order finite difference.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) XPLS(N), the new iterate, X[K].
!
!    Input, external FCN, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real x(n)
!      real f
!
!    Input, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Output, real ( kind = 8 ) A(NR,N), the N by N  finite difference
!    approximation to the hessian matrix.  Only the lower triangular matrix and
!    diagonal are returned.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling for X.
!
!    Input, real ( kind = 8 ) RNOISE, the relative noise in the function.
!
!    Workspace, real ( kind = 8 ) STEPSZ(N), used for the stepsize.
!
!    Workspace, real ( kind = 8 ) ANBR(N), holds neighbors.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) anbr(n)
  external fcn
  real ( kind = 8 ) fhat
  real ( kind = 8 ) fpls
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) ov3
  real ( kind = 8 ) rnoise
  real ( kind = 8 ) stepsz(n)
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) xpls(n)
  real ( kind = 8 ) xtmpi
  real ( kind = 8 ) xtmpj
!
!  Find I-th stepsize and evaluate neighbor in direction of I-th unit vector.
!
  ov3 = 1.0D+00 / 3.0D+00

  do i = 1, n
    stepsz(i) = rnoise**ov3 * max ( abs ( xpls(i) ), 1.0D+00 / sx(i) )
    xtmpi = xpls(i)
    xpls(i) = xtmpi + stepsz(i)
    call fcn ( n, xpls, anbr(i) )
    xpls(i) = xtmpi
  end do
!
!  Calculate column I of A.
!
  do i = 1, n

    xtmpi = xpls(i)
    xpls(i) = xtmpi + 2.0D+00 * stepsz(i)
    call fcn ( n, xpls, fhat )
    a(i,i) = ( ( fpls - anbr(i) ) &
      + ( fhat - anbr(i) ) ) / ( stepsz(i) * stepsz(i) )
!
!  Calculate sub-diagonal elements of column.
!
    if ( i /= n ) then

      xpls(i) = xtmpi + stepsz(i)

      do j = i + 1, n
        xtmpj = xpls(j)
        xpls(j) = xtmpj + stepsz(j)
        call fcn ( n, xpls, fhat )
        a(j,i) = ( ( fpls - anbr(i) ) + ( fhat - anbr(j) ) ) &
          / ( stepsz(i) * stepsz(j) )
        xpls(j) = xtmpj
      end do

    end if

    xpls(i) = xtmpi

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tregup ( nr, n, x, f, g, a, fcn, sc, sx, nwtake, stepmx, steptl, &
  dlt, iretcd, xplsp, fplsp, xpls, fpls, mxtake, ipr, method, udiag )

!*****************************************************************************80
!
!! TREGUP decides whether to accept the next optimization iterate.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate X[K-1].
!
!    Input, real ( kind = 8 ) F, the function value at the old iterate.
!
!    Input, real ( kind = 8 ) G(N), the gradient at the old iterate, or
!    an approximate value of it.
!
!    Input, real ( kind = 8 ) A(NR,N), the Cholesky decomposition of hessian in
!    lower triangular part and diagonal.  Hessian or approximation in
!    upper triangular part.
!
!    Input, external FCN, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) x(n)
!      real ( kind = 8 ) f
!
!    Input, real ( kind = 8 ) SC(N), the current step.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, logical NWTAKE, is TRUE if a Newton step is taken.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which successive
!    iterates are considered close enough to terminate algorithm.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Input/output, integer ( kind = 4 ) IRETCD, the status code.
!    0, xpls accepted as next iterate;  dlt trust region for next iteration.
!    1, xpls unsatisfactory but accepted as next iterate because
!      xpls-x < smallest allowable step length.
!    2, f(xpls) too large.  continue current iteration with new reduced dlt.
!    3, f(xpls) sufficiently small, but quadratic model predicts f(xpls)
!      sufficiently well to continue current iteration with new doubled dlt.
!
!    Workspace, real ( kind = 8 ) XPLSP(N), [value needs to be retained between
!    succesive calls of k-th global step].
!
!    Worskpace, real ( kind = 8 ) FPLSP, [retain value between successive
!    calls].
!
!    Output, real ( kind = 8 ) XPLS(N), the new iterate x[k].
!
!    Output, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Output, logical MXTAKE, is TRUE if a step of maximum length was taken.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Input, integer ( kind = 4 ) METHOD, the algorithm to use.
!    1, line search,
!    2, double dogleg,
!    3, More-Hebdon.
!
!    Input, real ( kind = 8 ) UDIAG(N), the diagonal of hessian in a(.,.)
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a(nr,n)
  real ( kind = 8 ) dlt
  real ( kind = 8 ) dltf
  real ( kind = 8 ) dltfp
  real ( kind = 8 ) dltmp
  real ( kind = 8 ) f
  external fcn
  real ( kind = 8 ) fpls
  real ( kind = 8 ) fplsp
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) iretcd
  integer ( kind = 4 ) j
  integer ( kind = 4 ) method
  logical mxtake
  logical nwtake
  real ( kind = 8 ) rln
  real ( kind = 8 ) sc(n)
  real ( kind = 8 ) slp
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) steptl
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) temp
  real ( kind = 8 ) udiag(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xpls(n)
  real ( kind = 8 ) xplsp(n)

  mxtake = .false.
  xpls(1:n) = x(1:n) + sc(1:n)
  call fcn ( n, xpls, fpls )
  dltf = fpls - f
  slp = dot_product ( g(1:n), sc(1:n) )
!
!  Next statement added for case of compilers which do not optimize
!  evaluation of next "if" statement (in which case fplsp could be
!  undefined).
!
  if ( iretcd == 4 ) then
    fplsp = 0.0D+00
  end if
!
!  Reset XPLS to XPLSP and terminate global step.
!
  if ( iretcd == 3 .and. ( fplsp <= fpls .or. 1.0D-04 * slp < dltf ) ) then
    iretcd = 0
    xpls(1:n) = xplsp(1:n)
    fpls = fplsp
    dlt = 0.5D+00 * dlt
    return
  end if

  if ( dltf <= 1.0D-04 * slp ) then
    go to 170
  end if

  rln = 0.0D+00

  do i = 1, n

    rln = max (                                        &
                rln,                                   &
                abs ( sc(i) ) / max (                  &
                                      abs ( xpls(i) ), &
                                      1.0D+00 / sx(i)  &
                                    )                  &
              )
  end do
!
!  Cannot find satisfactory xpls sufficiently distinct from x
!
  if ( rln < steptl ) then
    iretcd = 1
    return
  end if
!
!  Reduce trust region and continue global step
!
        iretcd = 2
        dltmp = -slp * dlt / ( 2.0D+00 * ( dltf - slp ) )

        if ( 0.1D+00 * dlt <= dltmp ) then
          go to 155
        end if

          dlt = 0.1D+00 * dlt
          go to 160

  155   continue
        dlt = dltmp

  160   continue
        return
!
!  FPLS sufficiently small
!
  170     continue

      dltfp = 0.0D+00

      if ( method == 2 ) then

        do i = 1, n
          temp = dot_product ( sc(i:n), a(i:n,i) )
          dltfp = dltfp + temp**2
        end do

      else

        do i = 1, n
          dltfp = dltfp + udiag(i) * sc(i) * sc(i)

          temp = 0.0D+00
          do j = i+1, n
            temp = temp + a(i,j) * sc(i) * sc(j)
          end do
          dltfp = dltfp + 2.0D+00 * temp
        end do

      end if

      dltfp = slp + dltfp / 2.0D+00

      if ( iretcd == 2 .or. &
         0.1D+00 * abs ( dltf ) < abs ( dltfp - dltf ) .or. &
         nwtake .or. &
         0.99D+00 * stepmx < dlt ) then
        go to 210
      end if
!
!  Double trust region and continue global step
!
        iretcd = 3
        xplsp(1:n) = xpls(1:n)
        fplsp = fpls
        dlt = min ( 2.0D+00 * dlt, stepmx )
        return
!
!  Accept XPLS as the next iterate.  Choose the new trust region.
!
  210       continue

  iretcd = 0

  if ( 0.99D+00 * stepmx < dlt ) then
    mxtake = .true.
  end if

  if ( dltf < 0.1D+00 * dltfp ) then
    if ( dltf <= 0.75D+00 * dltfp ) then
      dlt = min ( 2.0D+00 * dlt, stepmx )
    end if
  else
    dlt = 0.5D+00 * dlt
  end if

  return
end
subroutine uncmin ( n, x0, fcn, x, f, info, w, lw )

!*****************************************************************************80
!
!! UNCMIN minimizes a smooth nonlinear function of N variables.
!
!  Discussion:
!
!    A subroutine that computes the function value at any point
!    must be supplied.  Derivative values are not required.
!    This subroutine provides the simplest interface to the uncmin
!    minimization package.  The user has no control over options.
!
!    This routine uses a quasi-Newton algorithm with line search
!    to minimize the function represented by the subroutine fcn.
!    At each iteration, the nonlinear function is approximated
!    by a quadratic function derived from a taylor series.
!    The quadratic function is minimized to obtain a search direction,
!    and an approximate minimum of the nonlinear function along
!    the search direction is found using a line search.  The
!    algorithm computes an approximation to the second derivative
!    matrix of the nonlinear function using quasi-Newton techniques.
!
!    The uncmin package is quite general, and provides many options
!    for the user.  However, this subroutine is designed to be
!    easy to use, with few choices allowed.  For example:
!
!    1.  only function values need be computed.  first derivative
!    values are obtained by finite differencing.  this can be
!    very costly when the number of variables is large.
!
!    2.  it is assumed that the function values can be obtained
!    accurately (to an accuracy comparable to the precision of
!    the computer arithmetic).
!
!    3.  at most 150 iterations are allowed.
!
!    4.  it is assumed that the function values are well-scaled,
!    that is, that the optimal function value is not pathologically
!    large or small.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    Robert Schnabel, John Koontz, B E Weiss,
!    A modular system of algorithms for unconstrained minimization,
!    Report cu-cs-240-82,
!    Computer Science Department,
!    University of Colorado at Boulder, 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X0(N), an initial estimate of the minimizer.
!
!    Input, external FCN, the name of the routine to evaluate the minimization
!    function, of the form
!
!      subroutine fcn ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) x(n)
!      real ( kind = 8 ) f
!
!    Output, real ( kind = 8 ) X(N), the local minimizer.
!
!    Output, real ( kind = 8 ) F, the function value at X.
!
!    Output, integer ( kind = 4 ) INFO, termination code.
!    0:  optimal solution found.
!    1:  terminated with gradient small, X is probably optimal.
!    2:  terminated with stepsize small, X is probably optimal.
!    3:  lower point cannot be found, X is probably optimal.
!    4:  iteration limit (150) exceeded.
!    5:  too many large steps, function may be unbounded.
!    -1:  insufficient workspace.
!
!    Workspace, real ( kind = 8 ) W(LW).
!
!    Input, integer ( kind = 4 ) LW, the size of the workspace vector W, which
!    must be at least N * ( N + 10 ).
!
  implicit none

  integer ( kind = 4 ) lw
  integer ( kind = 4 ) n

  external d1fcn
  external d2fcn
  real ( kind = 8 ) dlt
  real ( kind = 8 ) epsm
  real ( kind = 8 ) f
  external fcn
  real ( kind = 8 ) fscale
  real ( kind = 8 ) gradtl
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) iagflg
  integer ( kind = 4 ) iahflg
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ig
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipr
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itnlim
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) iw2
  integer ( kind = 4 ) iw3
  integer ( kind = 4 ) iw4
  integer ( kind = 4 ) iw5
  integer ( kind = 4 ) iw6
  integer ( kind = 4 ) iw7
  integer ( kind = 4 ) iw8
  integer ( kind = 4 ) lwmin
  integer ( kind = 4 ) method
  integer ( kind = 4 ) msg
  integer ( kind = 4 ) ndigit
  integer ( kind = 4 ) nr
  real ( kind = 8 ) stepmx
  real ( kind = 8 ) steptl
  real ( kind = 8 ) w(lw)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0(n)
!
!  Subdivide workspace
!
  ig  = 1
  it  = ig  + n
  iw1 = it  + n
  iw2 = iw1 + n
  iw3 = iw2 + n
  iw4 = iw3 + n
  iw5 = iw4 + n
  iw6 = iw5 + n
  iw7 = iw6 + n
  iw8 = iw7 + n
  ia  = iw8 + n
  lwmin = ia + n*n-1

  if ( lw < lwmin ) then
    info = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNCMIN - Fatal error!'
    write ( *, '(a)' ) '  Insufficient workspace.'
    write ( *, '(a)' ) '  LW < LWMIN.'
    write ( *, '(a,i6)' ) '  LW = ', lw
    write ( *, '(a,i6)' ) '  LWMIN = ', lwmin
    stop
  end if
!
!  Set up parameters for OPTDRV.
!
!  parameters that should not be changed when using condensed code
!
! nr     = parameter used to divide workspace
! method = 1 (line search) -- do not change
! msg    = 9 => no printing, n=1 allowed
! iahflg = 1 => analytic hessian  supplied (0 otherwise)
! ipr    = device for output (irrelevant in current version)
! dlt    = (irrelevant parameter for method = 1)
! epsm   = machine epsilon
!
  nr = n
  method = 1
  msg = 9
  iahflg = 0
  ipr = 6
  dlt = -1.0D+00
  epsm = epsilon ( epsm )
!
! parameters that may be changed:
!
! iexp   = 1 means function expensive to evaluate (iexp = 0 => cheap)
! iagflg = 1 means analytic gradient supplied (0 otherwise)
! ndigit = -1 means optdrv assumes f is fully accurate
! itnlim = 150 = maximum number of iterations allowed
! gradtl = zero tolerance for gradient, for convergence tests
! stepmx = maximum allowable step size
! steptl = zero tolerance for step, for convergence tests
! fscale = typical order-of-magnitude size of function
! typsiz = typical order-of-magnitude size of x (stored in w(lt))
!
  iexp = 1
  iagflg = 0
  ndigit = -1
  itnlim = 150
  gradtl = epsm**(1.0D+00 / 3.0D+00 )
  stepmx = 0.0D+00
  steptl = sqrt ( epsm )
  fscale = 1.0D+00
  w(it:it+n-1) = 1.0D+00
!
!  Minimize function
!
  call optdrv ( nr, n, x0, fcn, d1fcn, d2fcn, w(it), fscale, method, iexp, &
    msg, ndigit, itnlim, iagflg, iahflg,  ipr, dlt, gradtl, stepmx, steptl, &
    x, f, w(ig), info, w(ia), w(iw1), w(iw2), w(iw3), w(iw4), &
    w(iw5), w(iw6), w(iw7), w(iw8))

  if ( info == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNCMIN - Note!'
    write ( *, '(a)' ) '  INFO = 1.'
    write ( *, '(a)' ) '  The iteration probably converged.'
    write ( *, '(a)' ) '  The gradient is very small.'
    return
  end if

  if ( info == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNCMIN - Note!'
    write ( *, '(a)' ) '  INFO = 2.'
    write ( *, '(a)' ) '  The iteration probably converged.'
    write ( *, '(a)' ) '  The stepsize is very small.'
    return
  end if

  if ( info == 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNCMIN - Warning!'
    write ( *, '(a)' ) '  INFO = 3.'
    write ( *, '(a)' ) '  Cannot find a point with lower value.'
    write ( *, '(a)' ) '  (But not completely happy with the current value.)'
    return
  end if

  if ( info == 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNCMIN - Warning!'
    write ( *, '(a)' ) '  INFO = 4.'
    write ( *, '(a)' ) '  Too many iterations.'
    return
  end if

  if ( info == 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNCMIN - Warning!'
    write ( *, '(a)' ) '  INFO = 5.'
    write ( *, '(a)' ) '  Too many large steps.'
    write ( *, '(a)' ) '  The function may be unbounded.'
    return
  end if

  return
end
function uni ( )

!*****************************************************************************80
!
!! UNI generates real uniform random numbers on [0,1).
!
!  Discussion:
!
!    Before calling UNI, initialize the generator by calling USTART with
!    a seed value.
!
!    Users can choose to run UNI in its default mode (requiring no user action)
!    which will generate the same sequence of numbers on any computer supporting
!    floating point numbers with at least 24 bit mantissas, or in a mode that
!    will generate numbers with a longer period on computers with
!    larger mantissas.
!
!    To exercise this option, before invoking USTART, insert the instruction
!
!      UBITS = UNIB ( K )
!
!    with
!
!      24 <= K
!
!    where K is the number of bits in the mantissa of your floating
!    point word (K=48 for cray, cyber 205).  UNIB returns the
!    floating point value of K that it actually used.
!      K input as <= 24, then UBITS=24.
!      K input as > 24, then UBITS=float(K)
!
!    If 24 < K, the sequence of numbers generated by UNI may differ
!    from one computer to another.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    George Marsaglia,
!    "Comments on the perfect uniform random number generator",
!    Unpublished notes,
!    Washington State University.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) UNI, a random number in the interval [0,1].
!
  implicit none

  real ( kind = 8 ), save :: c = 362436.0D+00 / 16777216.0D+00
  real ( kind = 8 ), parameter :: cd = 7654321.0D+00 / 16777216.0D+00
  real ( kind = 8 ), parameter :: cm = 16777213.0D+00 / 16777216.0D+00
  real ( kind = 8 ), parameter :: csave = 362436.0D+00 / 16777216.0D+00
  integer ( kind = 4 ), save :: i = 17
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iseed
  integer ( kind = 4 ), save :: j = 5
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) jj
  integer ( kind = 4 ), save :: k = 24
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m1
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ), save, dimension ( 17 ) :: u = (/ &
     0.8668672834288D+00,  0.3697986366357D+00,  0.8008968294805D+00, &
     0.4173889774680D+00,  0.8254561579836D+00,  0.9640965269077D+00, &
     0.4508667414265D+00,  0.6451309529668D+00,  0.1645456024730D+00, &
     0.2787901807898D+00,  0.06761531340295D+00, 0.9663226330820D+00, &
     0.01963343943798D+00, 0.02947398211399D+00, 0.1636231515294D+00, &
     0.3976343250467D+00,  0.2631008574685D+00 /)
  real ( kind = 8 ) uni
  real ( kind = 8 ) unib
  real ( kind = 8 ) ustart
!
!  The basic generator is Fibonacci.
!
  uni = u(i) - u(j)
  if ( uni < 0.0D+00 ) then
    uni = uni + 1.0D+00
  end if

  u(i) = uni

  i = i - 1
  if ( i == 0 ) then
    i = 17
  end if

  j = j - 1
  if ( j == 0 ) then
    j = 17
  end if
!
!  The second generator is congruential.
!
  c = c - cd
  if ( c < 0.0D+00 ) then
    c = c + cm
  end if
!
!  Combination generator.
!
  uni = uni - c
  if ( uni < 0.0D+00 ) then
    uni = uni + 1.0D+00
  end if

  return

entry ustart ( iseed )

!*****************************************************************************80
!
!! USTART is an entry into UNI that allows the user to specify the seed.
!
!          set up ...
!          convert iseed to four smallish positive integer ( kind = 4 )s.
!
    i1 = mod ( abs ( iseed ), 177 ) + 1
    j1 = mod ( abs ( iseed ), 167 ) + 1
    k1 = mod ( abs ( iseed ), 157 ) + 1
    l1 = mod ( abs ( iseed ), 147 ) + 1
!
!  Generate random bit pattern in array based on given seed.
!
    do ii = 1, 17

      s = 0.0D+00
      t = 0.5D+00
!
!  Do for each of the bits of mantissa of word
!  loop  over k bits, where k is defaulted to 24 but can
!  be changed by user call to unib(k)
!
      do jj = 1, k
        m1 = mod ( mod ( i1 * j1, 179 ) * k1, 179 )
        i1 = j1
        j1 = k1
        k1 = m1
        l1 = mod ( 53 * l1 + 1, 169 )
        if ( 32 <= mod ( l1 * m1, 64 ) ) then
          s = s + t
        end if
        t = 0.5D+00 * t
       end do

       u(ii) = s

    end do

    ustart = real ( iseed, kind = 8 )
    return

entry unib ( kk )

!*****************************************************************************80
!
!! UNIB ?
!
  if ( kk <= 24 ) then
    k = 24
  else
    k = kk
  end if

  unib = real ( k, kind = 8 )

  return
end
subroutine xerabt ( messg, nmessg )

!*****************************************************************************80
!
!! XERABT aborts program execution and prints an error message.
!
!  Discussion:
!
!    XERABT aborts the execution of the program.  The error message causing
!    the abort is given in the calling sequence, in case one needs it for
!    printing on a dayfile, for example.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed,
!    containing no more than 72 characters.
!
!    Input, integer ( kind = 4 ) NMESSG, the actual number of characters 
!    in MESSG.  If NMESSG is 0, no message is being supplied.
!
  implicit none

  character ( len = * ) messg
  integer ( kind = 4 ) nmessg

  if ( 0 < nmessg ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XERABT - Termination after fatal error!'
    write ( *, '(a)' ) trim ( messg )
  end if

  stop
end
subroutine xerclr

!*****************************************************************************80
!
!! XERCLR resets the current error number to zero.
!
!  Discussion:
!
!    This routine simply resets the current error number to zero.
!    This may be necessary to do in order to determine that
!    a certain error has occurred again since the last time
!    NUMXER was referenced.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk

  junk = j4save ( 1, 0, .true. )

  return
end
subroutine xerctl ( messg1, nmessg, nerr, level, kontrl )

!*****************************************************************************80
!
!! XERCTL allows user control over handling of individual errors.
!
!  Discussion:
!
!    Allows user control over handling of individual errors.
!    Just after each message is recorded, but before it is
!    processed any further (i.e., before it is printed or
!    a decision to abort is made), a call is made to XERCTL.
!    If the user has provided his own version of XERCTL, he
!    can then override the value of KONTRL used in processing
!    this message by redefining its value.
!
!    KONTRL may be set to any value from -2 to 2.
!    The meanings for KONTRL are the same as in XSETF, except
!    that the value of KONTRL changes only for this message.
!    If KONTRL is set to a value outside the range from -2 to 2,
!    it will be moved back into that range.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG1, the first word (only) of the error
!    message.
!
!    Input, integer ( kind = 4 ) NMESSG, same as in the call to XERROR 
!    or XERRWV.
!
!    Input, integer ( kind = 4 ) NERR, same as in the call to XERROR or XERRWV.
!
!    Input, integer ( kind = 4 ) LEVEL, same as in the call to XERROR or XERRWV.
!
!    Input/output, integer ( kind = 4 ) KONTRL.  On input, the current value of 
!    the control flag as set by a call to XSETF.  On output, the new value of
!    kontrl.  If KONTRL is not defined, it will remain at its original value.
!    This changed value affects only the current occurrence of the current
!    message.
!
  implicit none

  integer ( kind = 4 ) kontrl
  integer ( kind = 4 ) level
  character ( len = * ) messg1
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) nmessg

  return
end
subroutine xerdmp

!*****************************************************************************80
!
!! XERDMP prints the error tables and then clears them.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) kount

  call xersav ( ' ', 0, 0, 0, kount )

  return
end
subroutine xermax ( maxnum )

!*****************************************************************************80
!
!! XERMAX sets the maximum number of times any error message is to be printed.
!
!  Discussion:
!
!    XERMAX sets the maximum number of times any message is to be printed.
!    That is, non-fatal messages are not to be printed after they have
!    occured MAXNUM times.  Such non-fatal messages may be printed less than
!    MAXNUM times even if they occur MAXNUM times, if error suppression mode
!    (KONTRL=0) is ever in effect.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXNUM, the maximum number of times any one 
!    message is to be printed.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) maxnum

  junk = j4save ( 4, maxnum, .true. )

  return
end
subroutine xerprt ( messg, nmessg )

!*****************************************************************************80
!
!! XERPRT prints a message on each file indicated by xgetua.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be printed.
!
!    Input, integer ( kind = 4 ) NMESSG, the actual number of characters 
!    in MESSG.
!
  implicit none

  integer ( kind = 4 ) ichar
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lenmes
  integer ( kind = 4 ) lun(5)
  character ( len = * ) messg
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nunit
!
!  Obtain unit numbers and write line to each unit
!
  call xgetua ( lun, nunit )

  lenmes = len ( messg )

  do kunit = 1, nunit

     iunit = lun(kunit)

     do ichar = 1, lenmes, 72
        last = min ( ichar+71 , lenmes )
        if ( iunit == 0 ) then
          write (*,'(1x,a)') messg(ichar:last)
        else
          write (iunit,'(1x,a)') messg(ichar:last)
        end if
    end do

  end do

  return
end
subroutine xerror ( messg, nerr, level )

!*****************************************************************************80
!
!! XERROR processes a diagnostic error message.
!
!  Discussion:
!
!    XERROR processes a diagnostic message, in a manner
!    determined by the value of level and the current value
!    of the library error control flag, kontrl.
!    See subroutine xsetf for details.
!
!  Example:
!
!    call xerror('smooth -- num was zero.',1,2)
!
!    call xerror('integ  -- less than full accuracy achieved.',2,1)
!
!    call xerror( &
!      'rooter -- actual zero of f found before interval fully collapsed.',3,0)
!
!    call xerror('exp    -- underflows being set to zero.',1,-1)
!
!  Modified:
!
!    13 August 2005
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed,
!    containing no more than 72 characters.
!
!    Input, integer ( kind = 4 ) NERR, the error number associated with this 
!    message.  NERR must not be zero.
!
!    Input, integer ( kind = 4 ) LEVEL, the error category.
!    2 means this is an unconditionally fatal error.
!    1 means this is a recoverable error.  (i.e., it is
!      non-fatal if XSETF has been appropriately called.)
!    0 means this is a warning message only.
!    -1 means this is a warning message which is to be printed at most once,
!      regardless of how many times this call is executed.
!
  implicit none

  integer ( kind = 4 ) level
  character ( len = * ) messg
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) nmessg

  nmessg = len ( messg )

  call xerrwv ( messg, nmessg, nerr, level, 0, 0, 0, 0, 0.0D+00, 0.0D+00 )

  return
end
subroutine xerrwv ( messg, nmessg, nerr, level, ni, i1, i2, nr, r1, r2 )

!*****************************************************************************80
!
!! XERRWV processes an error message that includes numeric information.
!
!  Discussion:
!
!    XERRWV processes a diagnostic message, in a manner
!    determined by the value of level and the current value
!    of the library error control flag, kontrl.
!    (see subroutine xsetf for details.)
!    in addition, up to two integer ( kind = 4 ) values and two real
!    values may be printed along with the message.
!
!  Example:
!
!    call xerrwv ( 'smooth -- num (=i1) was zero.',29,1,2,1,num,0,0,0.,0.)
!
!    call xerrwv ( &
!      'quadxy -- requested error (r1) less than minimum(r2).', &
!      54,77,1,0,0,0,2,errreq,errmin)
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed.
!
!    Input, integer ( kind = 4 ) NMESSG, the number of characters in MESSG.
!
!    Input, integer ( kind = 4 ) NERR, the error number associated with
!    this message.  NERR must not be zero.
!
!    Input, integer ( kind = 4 ) LEVEL, the error category.
!    2 means this is an unconditionally fatal error.
!    1 means this is a recoverable error.  (i.e., it is
!      non-fatal if xsetf has been appropriately called.)
!    0 means this is a warning message only.
!    -1 means this is a warning message which is to be printed at most
!      once, regardless of how many times this call is executed.
!
!    Input, integer ( kind = 4 ) NI, the number of integer values to be
!    printed. (0 to 2)
!
!    Input, integer ( kind = 4 ) I1, I2, the first and second integer values.
!
!    Input, integer ( kind = 4 ) NR, the number of real values to be
!    printed. (0 to 2)
!
!    Input, real ( kind = 8 ) R1, R2, the first and second real values.
!
  implicit none

  character ( len = 37 ) form
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ifatal
  integer ( kind = 4 ) isizei
  integer ( kind = 4 ) isizef
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) kdummy
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) lerr
  integer ( kind = 4 ) level
  character ( len = 20 ) lfirst
  integer ( kind = 4 ) lkntrl
  integer ( kind = 4 ) llevel
  integer ( kind = 4 ) lmessg
  integer ( kind = 4 ) lun(5)
  integer ( kind = 4 ) maxmes
  character ( len = * ) messg
  integer ( kind = 4 ) mkntrl
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nunit
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
!
!  Get flags
!
  lkntrl = j4save ( 2, 0, .false. )
  maxmes = j4save ( 4, 0, .false. )
!
!  Check for valid input
!
  if ( 0 < nmessg .and. nerr /= 0 .and. -1 <= level .and. level <= 2 ) then
    go to 10
  end if

    if ( 0 < lkntrl ) then
      call xerprt('fatal error in...',17)
    end if

    call xerprt( 'XERROR -- invalid input', 23 )

    if ( 0 < lkntrl ) then
      call xerprt('job abort due to fatal error.',29)
    end if

    if ( 0 < lkntrl ) then
      call xersav ( ' ', 0, 0, 0, kdummy )
    end if

    call xerabt('XERROR -- invalid input',23)
    return

   10 continue
!
!  Record the message.
!
  junk = j4save(1,nerr,.true.)
  call xersav ( messg, nmessg, nerr, level, kount )
!
!  Let user override
!
  lfirst = messg
  lmessg = nmessg
  lerr = nerr
  llevel = level
  call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
!
!  Reset to original values.
!
  lmessg = nmessg
  lerr = nerr
  llevel = level
  lkntrl = max ( -2, min ( 2, lkntrl ) )
  mkntrl = abs ( lkntrl )
!
!  Decide whether to print message
!
  if ( llevel < 2 .and. lkntrl == 0 ) then
    go to 100
  end if

  if (((llevel == (-1)) .and. ( min ( 1, maxmes ) < kount ) ) &
    .or.((llevel == 0) .and. ( maxmes < kount )) &
    .or.((llevel == 1) .and. ( maxmes < kount ).and.(mkntrl==1) ) &
    .or.((llevel == 2) .and. ( max ( 1, maxmes ) < kount ) ) ) then
    go to 100
  end if

  if ( 0 < lkntrl ) then

    call xerprt(' ',1)

    if ( llevel == -1 ) then

      call xerprt &
      ( 'warning message...this message will only be printed once.',57)

    end if

    if ( llevel == 0 ) then
      call xerprt ( 'warning in...', 13 )
    else if ( llevel == 1 ) then
      call xerprt ( 'recoverable error in...', 23 )
    else if ( llevel == 2 ) then
      call xerprt ( 'fatal error in...', 17 )
    end if

  end if
!
!  Message
!
     call xerprt(messg,lmessg)
     call xgetua(lun,nunit)
     isizei = log10 ( real ( i1mach(9), kind = 8 ) ) + 1.0D+00
     isizef = log10 ( real ( i1mach(10), kind = 8 )**i1mach(14) ) + 1.0D+00

     do kunit = 1, nunit

        iunit = lun(kunit)

        do i = 1, min ( ni, 2 )
           write (form,21) i,isizei
   21          format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
           if ( iunit == 0 ) then
             if (i == 1) write (*,form) i1
             if (i == 2) write (*,form) i2
           else
             if (i == 1) write (iunit,form) i1
             if (i == 2) write (iunit,form) i2
           end if
        end do

        do i = 1, min ( nr, 2 )
           write (form,23) i,isizef+10,isizef
   23          format ('(11x,21hin above message, r',i1,'=,e',i2,'.',i2,')')
           if ( iunit == 0 ) then
             if ( i == 1 ) write (*,form) r1
             if ( i == 2 ) write (*,form) r2
           else
             if (i == 1) write (iunit,form) r1
             if (i == 2) write (iunit,form) r2
           end if
        end do

        if ( lkntrl <= 0 ) then
          go to 40
        end if
!
!  error number
!
           if ( iunit == 0 ) then
             write(*,30) lerr
           else
             write (iunit,30) lerr
           end if
   30          format (15h error number =,i10)
   40       continue

     end do
!
!  Traceback
!
  100 continue
  ifatal = 0
  if ((llevel == 2).or.((llevel==1) .and. (mkntrl==2))) then
    ifatal = 1
  end if
!
!  quit here if message is not fatal
!
  if ( ifatal <= 0 ) then
    return
  end if

  if ( lkntrl <= 0 .or. max ( 1, maxmes ) < kount ) then
    go to 120
  end if
!
!  Print reason for abort
!
     if ( llevel == 1 ) then
       call xerprt ('job abort due to unrecovered error.',35)
     end if

     if ( llevel == 2 ) then
       call xerprt('job abort due to fatal error.',29)
     end if
!
!  Print error summary
!
     call xersav ( ' ', -1, 0, 0, kdummy )

  120 continue
!
!  Abort
!
  if ( llevel == 2 .and. max ( 1, maxmes ) < kount ) then
    lmessg = 0
  end if

  call xerabt ( messg, lmessg )

  return
end
subroutine xersav ( messg, nmessg, nerr, level, icount )

!*****************************************************************************80
!
!! XERSAV records that an error occurred.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, as in XERROR.
!
!    Input, integer ( kind = 4 ) NMESSG, as in XERROR, except that, when 
!    NMESSG = 0, the tables will be dumped and cleared; and when NMESSG < 0,
!    the tables will be dumped, but not cleared.
!
!    Input, integer ( kind = 4 ) NERR, as in XERROR.
!
!    Input, integer ( kind = 4 ) LEVEL, as in XERROR.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of times this message has
!    been seen, or zero if the table has overflowed and
!    does not contain this message specifically.
!    when nmessg=0, icount will not be altered.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ), save, dimension ( 10 ) :: kount = (/ &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer ( kind = 4 ), save :: kountx = 0
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) level
  integer ( kind = 4 ), save, dimension ( 10 ) :: levtab
  integer ( kind = 4 ) lun(5)
  character ( len = 20 ) mes
  character ( len = * ) messg
  character ( len = 20 ), save, dimension ( 10 ) :: mestab
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ), save, dimension ( 10 ) :: nertab
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nunit
!
!  Dump the table
!
  if ( nmessg <= 0 ) then

     if ( kount(1) == 0 ) then
       return
     end if
!
!  Print to each unit
!
     call xgetua ( lun, nunit )

     do kunit = 1, nunit

       iunit = lun(kunit)

       if ( iunit == 0 ) then
         iunit = i1mach(4)
       end if
!
!  Print table header
!
       write ( iunit, 10 )
   10  format ('0          error message summary'/ &
       ' message start             nerr     level     count')
!
!  print body of table
!
        do i = 1, 10
          if ( kount(i) == 0 ) then
            exit
          end if
          write (iunit,15) mestab(i),nertab(i),levtab(i),kount(i)
   15     format (1x,a20,3i10)
        end do
!
!  Print number of other errors
!
        if ( kountx /= 0 ) then
          write (iunit,40) kountx
        end if

   40       format (41h0other errors not individually tabulated=,i10)
        write ( iunit, '(a)' ) ' '
     end do
!
!  Clear the error tables
!
    if ( nmessg == 0 ) then
      kount(1:10) = 0
      kountx = 0
    end if

    return

  end if
!
!  process a message...
!  search for this message, or else an empty slot for this messg,
!  or else determine that the error table is full.
!
  mes = messg

  do i = 1, 10

    ii = i

    if ( kount(i) == 0 ) then
      mestab(ii) = mes
      nertab(ii) = nerr
      levtab(ii) = level
      kount(ii)  = 1
      icount = 1
      return
    end if

    if ( mes /= mestab(i) ) then
      go to 90
    end if

    if (nerr /= nertab(i) ) then
      go to 90
    end if

    if (level /= levtab(i) ) then
      go to 90
    end if

    go to 100

90  continue

  end do
!
!  The table is full.
!
  kountx = kountx + 1
  icount = 1
  return
!
!  Message found in table
!
  100    continue

     kount(ii) = kount(ii) + 1
     icount = kount(ii)

  return
end
subroutine xgetf ( kontrl )

!*****************************************************************************80
!
!! XGETF returns current value of error control flag.
!
!  Discussion:
!
!    See subroutine XSETF for flag value meanings.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KONTRL, the current value of the error
!    control flag.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) kontrl

  kontrl = j4save ( 2, 0, .false. )

  return
end
subroutine xgetua ( iunita, n )

!*****************************************************************************80
!
!! XGETUA returns the unit number(s) to which error messages are being sent.
!
!  Discussion:
!
!    XGETUA may be called to determine the unit number or numbers to which
!    error messages are being sent.  These unit numbers may have been set
!    by a call to XSETUN, or a call to XSETUA, or may be a default value.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNITA(N),  an array unit numbers,
!    A value of zero refers to the default unit, as defined by the
!    I1MACH machine constant routine.  Only IUNITA(1),..., IUNITA(N) are
!    defined by XGETUA.  The values of IUNITA(N+1),..., IUNITA(5) are
!    not defined (for N < 5) or altered in any way by XGETUA.
!
!    Output, integer ( kind = 4 ) N, the number of units to which copies of the
!    error messages are being sent.  N will be in the range from 1 to 5.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) iunita(5)
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) n

  n = j4save ( 5, 0, .false. )

  do i = 1, n

    index = i+4
    if ( i == 1 ) then
      index = 3
    end if

    iunita(i) = j4save ( index, 0, .false. )

  end do

  return
end
subroutine xgetun ( iunit )

!*****************************************************************************80
!
!! XGETUN returns the (first) output file to which messages are being sent.
!
!  Discussion:
!
!    XGETUN gets the (first) output file to which error messages
!    are being sent.  To find out if more than one file is being
!    used, one must use the XGETUA routine.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the logical unit number of the
!    (first) unit to which error messages are being sent.  A value of zero
!    means that the default file, as defined by the I1MACH routine, is
!    being used.
!
  implicit none

  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j4save

  iunit = j4save ( 3, 0, .false. )

  return
end
subroutine xsetf ( kontrl )

!*****************************************************************************80
!
!! XSETF sets the error control flag.
!
!  Discussion:
!
!    XSETF sets the error control flag value to KONTRL.
!
!    The following table shows how each message is treated,
!    depending on the values of KONTRL and LEVEL.  See XERROR
!    for description of LEVEL.
!
!    if kontrl is zero or negative, no information other than the
!    message itself (including numeric values, if any) will be
!    printed.  if kontrl is positive, introductory messages,
!    tracebacks, etc., will be printed in addition to the message.
!
!              abs ( kontrl )
!        level        0              1              2
!        value
!          2        fatal          fatal          fatal
!
!          1     not printed      printed         fatal
!
!          0     not printed      printed        printed
!
!         -1     not printed      printed        printed
!                                  only           only
!                                  once           once
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KONTRL, the desired value of the error
!    control flag.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) kontrl

  if ( kontrl < -2 .or. 2 < kontrl ) then

    call xerrwv ( 'xsetf  -- invalid value of kontrl (i1).', 33, 1, 2, &
      1, kontrl, 0, 0, 0.0D+00, 0.0D+00 )

    return
  end if

  junk = j4save ( 2, kontrl, .true. )

  return
end
subroutine xsetua ( iunita, n )

!*****************************************************************************80
!
!! XSETUA sets up to 5 unit numbers to which messages are to be sent.
!
!  Discussion:
!
!    XSETUA may be called to declare a list of up to five
!    logical units, each of which is to receive a copy of
!    each error message processed by this package.
!    The purpose of XSETUA is to allow simultaneous printing
!    of each error message on, say, a main output file,
!    an interactive terminal, and other files such as graphics
!    communication files.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT(N), unit numbers, which should normally
!    be distinct.
!
!    Input, integer ( kind = 4 ) N, the number of unit numbers provided in IUNIT.
!    N must be between 1 and 5.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) iunita(n)
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk

  if ( n < 1 .or. 5 < n ) then
    call xerrwv ( 'xsetua -- invalid value of n (i1).', 34, 1, 2, 1, n, 0, &
      0, 0.0D+00, 0.0D+00 )
    return
  end if

  do i = 1, n
    index = i + 4
    if ( i  == 1 ) then
      index = 3
    end if
    junk = j4save ( index, iunita(i), .true. )
  end do

  junk = j4save ( 5, n, .true. )

  return
end
subroutine xsetun ( iunit )

!*****************************************************************************80
!
!! XSETUN sets the output file to which error messages are to be sent.
!
!  Discussion:
!
!    XSETUN sets the output file to which error messages are to be sent.
!    Only one file will be used.  See XSETUA for how to declare more than
!    one file.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the logical unit number to which error
!    messages are to be sent.
!
  implicit none

  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk

  junk = j4save ( 3, iunit, .true. )
  junk = j4save ( 5, 1, .true. )

  return
end
subroutine zfftb ( n, c, wsave )

!*****************************************************************************80
!
!! ZFFTB computes the backward complex discrete Fourier transform.
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    ZFFTB computes a complex periodic sequence from its Fourier coefficients.
!
!    A call of ZFFTF followed by a call of ZFFTB will multiply the
!    sequence by N.  In other words, the transforms are not normalized.
!
!    The array WSAVE must be initialized by ZFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N )
!        C_in(K) * exp ( sqrt ( - 1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be transformed.
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex ( kind = 8 ) C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, real ( kind = 8 ) WSAVE(4*N+15).  The array must be initialized
!    by calling ZFFTI.  A different WSAVE array must be used for each different
!    value of N.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  real ( kind = 8 ) wsave(4*n+15)

  if ( n <= 1 ) then
    return
  end if

  call zfftb1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine zfftb1 ( n, c, ch, wa, ifac )

!*****************************************************************************80
!
!! ZFFTB1 is a lower-level routine used by ZFFTB.
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.
!
!    Input/output, complex ( kind = 8 ) C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, complex ( kind = 8 ) CH(N).
!
!    Input, real ( kind = 8 ) WA(2*N).
!
!    Input, integer ( kind = 4 ) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  complex ( kind = 8 ) ch(n)
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ifac(15)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nac
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(2*n)

  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido

      if ( na == 0 ) then
        call passb4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passb4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passb2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passb2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido

      if ( na == 0 ) then
        call passb3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passb3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passb5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passb5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passb ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passb ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine zfftb_2d ( ldf, n, f, wsave )

!*****************************************************************************80
!
!! ZFFTB_2D computes a backward two dimensional complex fast Fourier transform.
!
!  Discussion:
!
!    The routine computes the backward two dimensional fast Fourier transform,
!    of a complex N by N matrix of data.
!
!    The output is unscaled, that is, a call to ZFFTB_2D followed by a call
!    to ZFFTF_2D will return the original data multiplied by N*N.
!
!    For some applications it is desirable to have the transform scaled so
!    the center of the N by N frequency square corresponds to zero
!    frequency.  The user can do this replacing the original input data
!    F(I,J) by F(I,J) * (-1.0)**(I+J),  I,J =0,...,N-1.
!
!    Before calling ZFFTF_2D or ZFFTB_2D, it is necessary to initialize
!    the array WSAVE by calling ZFFTI.
!
!  Modified:
!
!    12 March 2001
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDF, the leading dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrix.
!
!    Input/output, complex ( kind = 8 ) F(LDF,N),
!    On input, an N by N array of complex values to be transformed.
!    On output, the transformed values.
!
!    Input, real ( kind = 8 ) WSAVE(4*N+15), a work array whose values
!    depend on N, and which must be initialized by calling ZFFTI.
!
  implicit none

  integer ( kind = 4 ) ldf
  integer ( kind = 4 ) n

  complex ( kind = 8 ) f(ldf,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) wsave(4*n+15)
!
!  Row transforms:
!
  f(1:n,1:n) = transpose ( f(1:n,1:n) )

  do i = 1, n
    call zfftb ( n, f(1,i), wsave )
  end do

  f(1:n,1:n) = transpose ( f(1:n,1:n) )
!
!  Column transforms:
!
  do i = 1, n
    call zfftb ( n, f(1,i), wsave )
  end do

  return
end
subroutine zfftf ( n, c, wsave )

!*****************************************************************************80
!
!! ZFFTF computes the forward complex discrete Fourier transform.
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    ZFFTF computes the Fourier coefficients of a complex periodic sequence.
!
!    The transform is not normalized.  To obtain a normalized transform,
!    the output must be divided by N.  Otherwise a call of ZFFTF
!    followed by a call of ZFFTB will multiply the sequence by N.
!
!    The array WSAVE must be initialized by calling ZFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N )
!        C_in(K) * exp ( - sqrt ( -1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be transformed.
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex ( kind = 8 ) C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, real ( kind = 8 ) WSAVE(4*N+15).  The array must be initialized
!    by calling ZFFTI.  A different WSAVE array must be used for each different
!    value of N.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  real ( kind = 8 ) wsave(4*n+15)

  if ( n <= 1 ) then
    return
  end if

  call zfftf1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine zfftf1 ( n, c, ch, wa, ifac )

!*****************************************************************************80
!
!! ZFFTF1 is a lower level routine used by ZFFTF.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be transformed.
!
!    Input/output, complex ( kind = 8 ) C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, complex ( kind = 8 ) CH(N).
!
!    Input, real ( kind = 8 ) WA(2*N).
!
!    Input, integer ( kind = 4 ) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  complex ( kind = 8 ) ch(n)
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ifac(15)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nac
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(2*n)

  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido

      if ( na == 0 ) then
        call passf4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passf4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passf2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passf2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido

      if ( na == 0 ) then
        call passf3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passf3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passf5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passf5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passf ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passf ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine zfftf_2d ( ldf, n, f, wsave )

!*****************************************************************************80
!
!! ZFFTF_2D computes a two dimensional complex fast Fourier transform.
!
!  Discussion:
!
!    The routine computes the forward two dimensional fast Fourier transform,
!    of a complex N by N matrix of data.
!
!    The output is unscaled, that is, a call to ZFFTF_2D,
!    followed by a call to ZFFTB_2D will return the original data
!    multiplied by N*N.
!
!    For some applications it is desirable to have the transform scaled so
!    the center of the N by N frequency square corresponds to zero
!    frequency.  The user can do this replacing the original input data
!    F(I,J) by F(I,J) *(-1.0)**(I+J),  I,J =0,...,N-1.
!
!    Before calling ZFFTF_2D or ZFFTB_2D, it is necessary to initialize
!    the array WSAVE by calling ZFFTI.
!
!  Modified:
!
!    12 March 2001
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDF, the leading dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in 
!    the matrix.
!
!    Input/output, complex ( kind = 8 ) F(LDF,N),
!    On input, an N by N array of complex values to be transformed.
!    On output, the transformed values.
!
!    Input, real ( kind = 8 ) WSAVE(4*N+15), a work array whose values depend
!    on N, and which must be initialized by calling CFFTI.
!
  implicit none

  integer ( kind = 4 ) ldf
  integer ( kind = 4 ) n

  complex ( kind = 8 ) f(ldf,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) wsave(4*n+15)
!
!  Row transforms:
!
  f(1:n,1:n) = transpose ( f(1:n,1:n) )

  do i = 1, n
    call zfftf ( n, f(1,i), wsave )
  end do

  f(1:n,1:n) = transpose ( f(1:n,1:n) )
!
!  Column transforms:
!
  do i = 1, n
    call zfftf ( n, f(1,i), wsave )
  end do

  return
end
subroutine zffti ( n, wsave )

!*****************************************************************************80
!
!! ZFFTI initializes WSAVE, used in ZFFTF and ZFFTB.
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Paul Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    Bill Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to 
!    be transformed.
!
!    Output, real ( kind = 8 ) WSAVE(4*N+15), contains data, dependent on
!    the value of N, which is necessary for the CFFTF or CFFTB routines.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) wsave(4*n+15)

  if ( n <= 1 ) then
    return
  end if

  call zffti1 ( n, wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine zffti1 ( n, wa, ifac )

!*****************************************************************************80
!
!! ZFFTI1 is a lower level routine used by ZFFTI.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be transformed.
!
!    Input, real ( kind = 8 ) WA(2*N).
!
!    Input, integer ( kind = 4 ) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) argh
  real ( kind = 8 ) argld
  real ( kind = 8 ) fi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ifac(15)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) ld
  integer ( kind = 4 ) nf
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) wa(2*n)

  call i8_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0D+00 * pi / real ( n, kind = 8 )
  i = 2
  l1 = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      i1 = i
      wa(i-1) = 1.0D+00
      wa(i) = 0.0D+00
      ld = ld + l1
      fi = 0.0D+00
      argld = real ( ld, kind = 8 ) * argh

      do ii = 4, 2*ido+2, 2
        i = i + 2
        fi = fi + 1.0D+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      if ( 5 < ip ) then
        wa(i1-1) = wa(i-1)
        wa(i1) = wa(i)
      end if

    end do

    l1 = l2

  end do

  return
end
