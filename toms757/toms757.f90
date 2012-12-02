function abram0 ( xvalue )

!*****************************************************************************80
!
!! ABRAM0 evaluates the Abramowitz function of order 0.
!
!  Discussion:
!
!    The function is defined by:
!
!      ABRAM0(x) = Integral ( 0 <= t < infinity ) exp ( -t^2 - x / t ) dt
!
!    The code uses Chebyshev expansions with the coefficients
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) ABRAM0, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension(0:8) :: ab0f = (/ &
    -0.68121927093549469816d0, &
    -0.78867919816149252495d0, &
     0.5121581776818819543d-1, &
    -0.71092352894541296d-3, &
     0.368681808504287d-5, &
    -0.917832337237d-8, &
     0.1270202563d-10, &
    -0.1076888d-13, &
     0.599d-17 /)
  real ( kind = 8 ), dimension(0:8) :: ab0g = (/ &
     -0.60506039430868273190d0, &
     -0.41950398163201779803d0, &
      0.1703265125190370333d-1, &
     -0.16938917842491397d-3, &
      0.67638089519710d-6, &
     -0.135723636255d-8, &
      0.156297065d-11, &
     -0.112887d-14, &
      0.55d-18 /)
  real ( kind = 8 ), dimension(0:8) :: ab0h = (/ &
      1.38202655230574989705d0, &
     -0.30097929073974904355d0, &
      0.794288809364887241d-2, &
     -0.6431910276847563d-4, &
      0.22549830684374d-6, &
     -0.41220966195d-9, &
      0.44185282d-12, &
     -0.30123d-15, &
      0.14d-18 /)
  real ( kind = 8 ), dimension(0:27) :: ab0as = (/ &
     1.97755499723693067407d+0, &
    -0.1046024792004819485d-1, &
     0.69680790253625366d-3, &
    -0.5898298299996599d-4, &
     0.577164455305320d-5, &
    -0.61523013365756d-6, &
     0.6785396884767d-7, &
    -0.723062537907d-8, &
     0.63306627365d-9, &
    -0.989453793d-11, &
    -0.1681980530d-10, &
     0.673799551d-11, &
    -0.200997939d-11, &
     0.54055903d-12, &
    -0.13816679d-12, &
     0.3422205d-13, &
    -0.826686d-14, &
     0.194566d-14, &
    -0.44268d-15, &
     0.9562d-16, &
    -0.1883d-16, &
     0.301d-17, &
    -0.19d-18, &
    -0.14d-18, &
     0.11d-18, &
    -0.4d-19, &
     0.2d-19, &
    -0.1d-19 /)
  real ( kind = 8 ) abram0
  real ( kind = 8 ) asln
  real ( kind = 8 ) asval
  real ( kind = 8 ) cheval
  real ( kind = 8 ) fval
  real ( kind = 8 ) gval
  real ( kind = 8 ), parameter :: gval0 = 0.13417650264770070909D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  real ( kind = 8 ) hval
  real ( kind = 8 ), parameter :: lnxmin = -708.3964D+00
  integer, parameter :: nterma = 22
  integer, parameter :: ntermf = 8
  integer, parameter :: ntermg = 8
  integer, parameter :: ntermh = 8
  real ( kind = 8 ), parameter :: onerpi = 0.56418958354775628695D+00
  real ( kind = 8 ), parameter :: rt3bpi = 0.97720502380583984317D+00
  real ( kind = 8 ), parameter :: rtpib2 = 0.88622692545275801365D+00
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) v
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xlow1 = 1.490116D-08
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ABRAM0 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    abram0 = zero
  else if ( x == zero ) then
    abram0 = rtpib2
  else if ( x < xlow1 ) then
    abram0 = rtpib2 + x * ( log ( x ) - gval0 )
  else if ( x <= two ) then
    t =  ( x * x / two - half ) - half
    fval = cheval ( ntermf, ab0f, t )
    gval = cheval ( ntermg, ab0g, t )
    hval = cheval ( ntermh, ab0h, t )
    abram0 = fval / onerpi + x * ( log ( x ) * hval - gval )
  else
    v = three * ( ( x / two ) ** ( two / three ) )
    t = ( six / v - half ) - half
    asval = cheval ( nterma, ab0as, t )
    asln = log ( asval / rt3bpi ) - v
    if ( asln < lnxmin ) then
      abram0 = zero
    else
      abram0 = exp ( asln )
    end if
  end if

  return
end
subroutine abram0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ABRAM0_VALUES returns some values of the Abramowitz0 function.
!
!  Discussion:
!
!    The function is defined by:
!
!      ABRAM0(x) = Integral ( 0 <= t < infinity ) exp ( -t^2 - x / t ) dt
!
!    The data was reported by McLeod.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function abram1 ( xvalue )

!*****************************************************************************80
!
!! ABRAM1 evaluates the Abramowitz function of order 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      ABRAM1(x) = Integral ( 0 <= t < infinity ) t * exp ( -t^2 - x / t ) dt
!
!    The code uses Chebyshev expansions with the coefficients
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) ABRAM1, the value of the function.
!
  implicit none

  real ( kind = 8 ) ab1as(0:27)
  real ( kind = 8 ) ab1f(0:9)
  real ( kind = 8 ) ab1g(0:8)
  real ( kind = 8 ) ab1h(0:8)
  real ( kind = 8 ) abram1
  real ( kind = 8 ) asln
  real ( kind = 8 ) asval
  real ( kind = 8 ) cheval
  real ( kind = 8 ) fval
  real ( kind = 8 ) gval
  real ( kind = 8 ), parameter :: half = 0.5D+00
  real ( kind = 8 ) hval
  real ( kind = 8 ) lnxmin
  integer, parameter :: nterma = 23
  integer, parameter :: ntermf = 9
  integer, parameter :: ntermg = 8
  integer, parameter :: ntermh = 8
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) onerpi
  real ( kind = 8 ), parameter :: rt3bpi = 0.97720502380583984317D+00
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) v
  real ( kind = 8 ) x
  real ( kind = 8 ) xlow
  real ( kind = 8 ) xlow1
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  data ab1f/1.47285192577978807369d0, &
            0.10903497570168956257d0, &
           -0.12430675360056569753d0, &
            0.306197946853493315d-2, &
           -0.2218410323076511d-4, &
            0.6989978834451d-7, &
           -0.11597076444d-9, &
            0.11389776d-12, &
           -0.7173d-16, &
            0.3d-19/
  data ab1g/0.39791277949054503528d0, &
           -0.29045285226454720849d0, &
            0.1048784695465363504d-1, &
           -0.10249869522691336d-3, &
            0.41150279399110d-6, &
           -0.83652638940d-9, &
            0.97862595d-12, &
           -0.71868d-15, &
            0.35d-18/
  data ab1h/0.84150292152274947030d0, &
           -0.7790050698774143395d-1, &
            0.133992455878390993d-2, &
           -0.808503907152788d-5, &
            0.2261858281728d-7, &
           -0.3441395838d-10, &
            0.3159858d-13, &
           -0.1884d-16, &
            0.1d-19/
  data ab1as(0)/  2.13013643429065549448d0/
  data ab1as(1)/  0.6371526795218539933d-1/
  data ab1as(2)/ -0.129334917477510647d-2/
  data ab1as(3)/  0.5678328753228265d-4/
  data ab1as(4)/ -0.279434939177646d-5/
  data ab1as(5)/  0.5600214736787d-7/
  data ab1as(6)/  0.2392009242798d-7/
  data ab1as(7)/ -0.750984865009d-8/
  data ab1as(8)/  0.173015330776d-8/
  data ab1as(9)/ -0.36648877955d-9/
  data ab1as(10)/ 0.7520758307d-10/
  data ab1as(11)/-0.1517990208d-10/
  data ab1as(12)/ 0.301713710d-11/
  data ab1as(13)/-0.58596718d-12/
  data ab1as(14)/ 0.10914455d-12/
  data ab1as(15)/-0.1870536d-13/
  data ab1as(16)/ 0.262542d-14/
  data ab1as(17)/-0.14627d-15/
  data ab1as(18)/-0.9500d-16/
  data ab1as(19)/ 0.5873d-16/
  data ab1as(20)/-0.2420d-16/
  data ab1as(21)/ 0.868d-17/
  data ab1as(22)/-0.290d-17/
  data ab1as(23)/ 0.93d-18/
  data ab1as(24)/-0.29d-18/
  data ab1as(25)/ 0.9d-19/
  data ab1as(26)/-0.3d-19/
  data ab1as(27)/ 0.1d-19/
  data onerpi/ 0.56418958354775628695d0/
!
!  Machine-dependent constants (suitable for IEEE machines)
!
  data xlow,xlow1,lnxmin/1.11023d-16,1.490116d-8,-708.3964d0/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ABRAM1 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    abram1 = zero
  else if ( x == zero ) then
    abram1 = half
  else if ( x < xlow ) then
    abram1 = half
  else if ( x < xlow1 ) then
    abram1 = ( one - x / onerpi - x * x * log ( x ) ) * half
  else if ( x <= two ) then
    t = ( x * x / two - half ) - half
    fval = cheval ( ntermf, ab1f, t )
    gval = cheval ( ntermg, ab1g, t )
    hval = cheval ( ntermh, ab1h, t )
    abram1 = fval - x * ( gval / onerpi + x * log ( x ) * hval )
  else
    v = three *  ( ( x / two ) ** ( two / three ) )
    t =  ( six / v - half ) - half
    asval = cheval ( nterma, ab1as, t )
    asln = log ( asval * sqrt ( v / three ) / rt3bpi ) - v
    if ( asln < lnxmin ) then
      abram1 = zero
    else
      abram1 = exp ( asln )
    end if
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
!      ABRAM1(x) = Integral ( 0 <= t < infinity ) t * exp ( -t^2 - x / t ) dt
!
!    The data was reported by McLeod.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function abram2 ( xvalue )

!*****************************************************************************80
!
!! ABRAM2 evaluates the Abramowitz function of order 2.
!
!  Discussion:
!
!    The function is defined by:
!
!      ABRAM2(x) = Integral ( 0 <= t < infinity ) t^2 * exp ( -t^2 - x / t ) dt
!
!    The code uses Chebyshev expansions with the coefficients
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) ABRAM2, the value of the function.
!
  implicit none

  real ( kind = 8 ) abram2
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer, parameter :: nterma = 23
  integer, parameter :: ntermf = 9
  integer, parameter :: ntermg = 8
  integer, parameter :: ntermh = 7
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) ab2f(0:9),ab2g(0:8),ab2h(0:7),ab2as(0:26), &
       asln,asval,fval,gval,hval,lnxmin, &
       onerpi,rtpib4,rt3bpi,t, &
       v,xlow,xlow1
  data ab2f/1.03612162804243713846d0, &
            0.19371246626794570012d0, &
           -0.7258758839233007378d-1, &
            0.174790590864327399d-2, &
           -0.1281223233756549d-4, &
            0.4115018153651d-7, &
           -0.6971047256d-10, &
            0.6990183d-13, &
           -0.4492d-16, &
            0.2d-19/
  data ab2g/1.46290157198630741150d0, &
            0.20189466883154014317d0, &
           -0.2908292087997129022d-1, &
            0.47061049035270050d-3, &
           -0.257922080359333d-5, &
            0.656133712946d-8, &
           -0.914110203d-11, &
            0.774276d-14, &
           -0.429d-17/
  data ab2h/0.30117225010910488881d0, &
           -0.1588667818317623783d-1, &
            0.19295936935584526d-3, &
           -0.90199587849300d-6, &
            0.206105041837d-8, &
           -0.265111806d-11, &
            0.210864d-14, &
           -0.111d-17/
  data ab2as(0)/  2.46492325304334856893d0/
  data ab2as(1)/  0.23142797422248905432d0/
  data ab2as(2)/ -0.94068173010085773d-3/
  data ab2as(3)/  0.8290270038089733d-4/
  data ab2as(4)/ -0.883894704245866d-5/
  data ab2as(5)/  0.106638543567985d-5/
  data ab2as(6)/ -0.13991128538529d-6/
  data ab2as(7)/  0.1939793208445d-7/
  data ab2as(8)/ -0.277049938375d-8/
  data ab2as(9)/  0.39590687186d-9/
  data ab2as(10)/-0.5408354342d-10/
  data ab2as(11)/ 0.635546076d-11/
  data ab2as(12)/-0.38461613d-12/
  data ab2as(13)/-0.11696067d-12/
  data ab2as(14)/ 0.6896671d-13/
  data ab2as(15)/-0.2503113d-13/
  data ab2as(16)/ 0.785586d-14/
  data ab2as(17)/-0.230334d-14/
  data ab2as(18)/ 0.64914d-15/
  data ab2as(19)/-0.17797d-15/
  data ab2as(20)/ 0.4766d-16/
  data ab2as(21)/-0.1246d-16/
  data ab2as(22)/ 0.316d-17/
  data ab2as(23)/-0.77d-18/
  data ab2as(24)/ 0.18d-18/
  data ab2as(25)/-0.4d-19/
  data ab2as(26)/ 0.1d-19/
  data rt3bpi/ 0.97720502380583984317d0/
  data rtpib4/ 0.44311346272637900682d0/
  data onerpi/ 0.56418958354775628695d0/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow,xlow1,lnxmin/2.22045d-16,1.490116d-8,-708.3964d0/

  x = xvalue

  if ( x < zero ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ABRAM2 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    abram2 = zero

  else if ( x == zero ) then

    abram2 = rtpib4

  else if ( x < xlow ) then

    abram2 = rtpib4

  else if ( x < xlow1 ) then

    abram2 = rtpib4 - half * x + x * x * x * log ( x ) / six

  else if ( x <= 2.0D+00 ) then

    t =  ( x * x / two - half ) - half
    fval = cheval ( ntermf, ab2f, t )
    gval = cheval ( ntermg, ab2g, t )
    hval = cheval ( ntermh, ab2h, t )
    abram2 = fval / onerpi + x * ( x * x * log ( x ) * hval - gval )

  else

    v = three * ( ( x / two ) ** ( two / three ) )
    t = ( six / v - half ) - half
    asval = cheval ( nterma, ab2as, t )
    asln = log ( asval / rt3bpi ) + log ( v / three ) - v

    if ( asln < lnxmin ) then
      abram2 = zero
    else
      abram2 = exp ( asln )
    end if

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
!      ABRAM2(x) = Integral ( 0 <= t < infinity ) t^2 * exp ( -t^2 - x / t ) dt
!
!    The data was reported by McLeod.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function airy_ai_int ( xvalue )

!*****************************************************************************80
!
!! AIRY_AI_INT calculates the integral of the Airy function Ai.
!
!  Discussion:
!
!    The function is defined by:
!
!      AIRY_AI_INT(x) = Integral ( 0 <= t <= x ) Ai(t) dt
!
!    The program uses Chebyshev expansions, the coefficients of which
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) AIRY_AI_INT, the value of the function.
!
  implicit none

  real ( kind = 8 ) aaint1(0:25)
  real ( kind = 8 ) aaint2(0:21)
  real ( kind = 8 ) aaint3(0:40)
  real ( kind = 8 ) aaint4(0:17)
  real ( kind = 8 ) aaint5(0:17)
  real ( kind = 8 ) airy_ai_int
  real ( kind = 8 ) airzer
  real ( kind = 8 ) arg
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ) forty1
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ) fr996
  real ( kind = 8 ) gval
  real ( kind = 8 ) hval
  real ( kind = 8 ) nine
  real ( kind = 8 ) ninhun
  integer, parameter :: nterm1 = 22
  integer, parameter :: nterm2 = 17
  integer, parameter :: nterm3 = 37
  integer nterm4
  integer nterm5
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) piby4
  real ( kind = 8 ) pitim6
  real ( kind = 8 ) rt2b3p
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xhigh1
  real ( kind = 8 ) xlow1
  real ( kind = 8 ) xneg1
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00
  real ( kind = 8 ) z

  data aaint1(0)/  0.37713517694683695526d0/
  data aaint1(1)/ -0.13318868432407947431d0/
  data aaint1(2)/  0.3152497374782884809d-1/
  data aaint1(3)/ -0.318543076436574077d-2/
  data aaint1(4)/ -0.87398764698621915d-3/
  data aaint1(5)/  0.46699497655396971d-3/
  data aaint1(6)/ -0.9544936738983692d-4/
  data aaint1(7)/  0.542705687156716d-5/
  data aaint1(8)/  0.239496406252188d-5/
  data aaint1(9)/ -0.75690270205649d-6/
  data aaint1(10)/ 0.9050138584518d-7/
  data aaint1(11)/ 0.320529456043d-8/
  data aaint1(12)/-0.303825536444d-8/
  data aaint1(13)/ 0.48900118596d-9/
  data aaint1(14)/-0.1839820572d-10/
  data aaint1(15)/-0.711247519d-11/
  data aaint1(16)/ 0.151774419d-11/
  data aaint1(17)/-0.10801922d-12/
  data aaint1(18)/-0.963542d-14/
  data aaint1(19)/ 0.313425d-14/
  data aaint1(20)/-0.29446d-15/
  data aaint1(21)/-0.477d-17/
  data aaint1(22)/ 0.461d-17/
  data aaint1(23)/-0.53d-18/
  data aaint1(24)/ 0.1d-19/
  data aaint1(25)/ 0.1d-19/
  data aaint2(0)/  1.92002524081984009769d0/
  data aaint2(1)/ -0.4220049417256287021d-1/
  data aaint2(2)/ -0.239457722965939223d-2/
  data aaint2(3)/ -0.19564070483352971d-3/
  data aaint2(4)/ -0.1547252891056112d-4/
  data aaint2(5)/ -0.140490186137889d-5/
  data aaint2(6)/ -0.12128014271367d-6/
  data aaint2(7)/ -0.1179186050192d-7/
  data aaint2(8)/ -0.104315578788d-8/
  data aaint2(9)/ -0.10908209293d-9/
  data aaint2(10)/-0.929633045d-11/
  data aaint2(11)/-0.110946520d-11/
  data aaint2(12)/-0.7816483d-13/
  data aaint2(13)/-0.1319661d-13/
  data aaint2(14)/-0.36823d-15/
  data aaint2(15)/-0.21505d-15/
  data aaint2(16)/ 0.1238d-16/
  data aaint2(17)/-0.557d-17/
  data aaint2(18)/ 0.84d-18/
  data aaint2(19)/-0.21d-18/
  data aaint2(20)/ 0.4d-19/
  data aaint2(21)/-0.1d-19/
  data aaint3(0)/  0.47985893264791052053d0/
  data aaint3(1)/ -0.19272375126169608863d0/
  data aaint3(2)/  0.2051154129525428189d-1/
  data aaint3(3)/  0.6332000070732488786d-1/
  data aaint3(4)/ -0.5093322261845754082d-1/
  data aaint3(5)/  0.1284424078661663016d-1/
  data aaint3(6)/  0.2760137088989479413d-1/
  data aaint3(7)/ -0.1547066673866649507d-1/
  data aaint3(8)/ -0.1496864655389316026d-1/
  data aaint3(9)/  0.336617614173574541d-2/
  data aaint3(10)/ 0.530851163518892985d-2/
  data aaint3(11)/ 0.41371226458555081d-3/
  data aaint3(12)/-0.102490579926726266d-2/
  data aaint3(13)/-0.32508221672025853d-3/
  data aaint3(14)/ 0.8608660957169213d-4/
  data aaint3(15)/ 0.6671367298120775d-4/
  data aaint3(16)/ 0.449205999318095d-5/
  data aaint3(17)/-0.670427230958249d-5/
  data aaint3(18)/-0.196636570085009d-5/
  data aaint3(19)/ 0.22229677407226d-6/
  data aaint3(20)/ 0.22332222949137d-6/
  data aaint3(21)/ 0.2803313766457d-7/
  data aaint3(22)/-0.1155651663619d-7/
  data aaint3(23)/-0.433069821736d-8/
  data aaint3(24)/-0.6227777938d-10/
  data aaint3(25)/ 0.26432664903d-9/
  data aaint3(26)/ 0.5333881114d-10/
  data aaint3(27)/-0.522957269d-11/
  data aaint3(28)/-0.382229283d-11/
  data aaint3(29)/-0.40958233d-12/
  data aaint3(30)/ 0.11515622d-12/
  data aaint3(31)/ 0.3875766d-13/
  data aaint3(32)/ 0.140283d-14/
  data aaint3(33)/-0.141526d-14/
  data aaint3(34)/-0.28746d-15/
  data aaint3(35)/ 0.923d-17/
  data aaint3(36)/ 0.1224d-16/
  data aaint3(37)/ 0.157d-17/
  data aaint3(38)/-0.19d-18/
  data aaint3(39)/-0.8d-19/
  data aaint3(40)/-0.1d-19/
  data aaint4/1.99653305828522730048d0, &
             -0.187541177605417759d-2, &
             -0.15377536280305750d-3, &
             -0.1283112967682349d-4, &
             -0.108128481964162d-5, &
             -0.9182131174057d-7, &
             -0.784160590960d-8, &
             -0.67292453878d-9, &
             -0.5796325198d-10, &
             -0.501040991d-11, &
             -0.43420222d-12, &
             -0.3774305d-13, &
             -0.328473d-14, &
             -0.28700d-15, &
             -0.2502d-16, &
             -0.220d-17, &
             -0.19d-18, &
             -0.2d-19/
  data aaint5/1.13024602034465716133d0, &
             -0.464718064639872334d-2, &
             -0.35137413382693203d-3, &
             -0.2768117872545185d-4, &
             -0.222057452558107d-5, &
             -0.18089142365974d-6, &
             -0.1487613383373d-7, &
             -0.123515388168d-8, &
             -0.10310104257d-9, &
             -0.867493013d-11, &
             -0.73080054d-12, &
             -0.6223561d-13, &
             -0.525128d-14, &
             -0.45677d-15, &
             -0.3748d-16, &
             -0.356d-17, &
             -0.23d-18, &
             -0.4d-19/
  data nine,forty1/ 9.0d0, 41.0d0/
  data ninhun,fr996/ 900.0d0, 4996.0d0 /
  data piby4/0.78539816339744830962d0/
  data pitim6/18.84955592153875943078d0/
  data rt2b3p/0.46065886596178063902d0/
  data airzer/0.35502805388781723926d0/
!
!   Machine-dependant constants (suitable for IEEE machines)
!
  data nterm4,nterm5/15,15/
  data xlow1,xhigh1,xneg1/2.22045d-16,14.480884d0,-2.727134d10/

  x = xvalue

  if ( x < xneg1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AIRY_AI_INT - Fatal error!'
    write ( *, '(a)' ) '  X too negative for accurate computation.'
    airy_ai_int = -two / three
    return
  else if ( x < -eight ) then
    z = - ( x + x ) * sqrt ( -x ) / three
    arg = z + piby4
    temp = nine * z * z
    t = ( fr996 - temp ) / ( ninhun + temp )
    gval = cheval ( nterm4, aaint4, t )
    hval = cheval ( nterm5, aaint5, t )
    temp = gval * cos ( arg ) + hval * sin ( arg ) / z
    airy_ai_int = rt2b3p * temp / sqrt ( z ) - two / three
  else if ( x <= -xlow1 )then
    t = -x / four - one
    airy_ai_int = x * cheval ( nterm3, aaint3, t )
  else if ( x < xlow1 ) then
    airy_ai_int = airzer * x
  else if ( x <= four ) then
    t = x / two - one
    airy_ai_int = cheval ( nterm1, aaint1, t ) * x
  else if ( x <= xhigh1 ) then
    z = ( x + x ) * sqrt ( x ) / three
    temp = three * z
    t = ( forty1 - temp ) / ( nine + temp )
    temp = exp ( -z ) * cheval ( nterm2, aaint2, t ) / sqrt ( pitim6 * z )
    airy_ai_int = one / three - temp
  else
    airy_ai_int = one / three
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
!      AIRY_AI_INT(x) = Integral ( 0 <= t <= x ) Ai(t) dt
!
!    The data was reported by McLeod.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function airy_bi_int ( xvalue )

!*****************************************************************************80
!
!! AIRY_BI_INT calculates the integral of the Airy function Bi.
!
!  Discussion:
!
!    The function is defined by:
!
!      AIRY_BI_INT(x) = Integral ( 0 <= t <= x ) Bi(t) dt
!
!    The program uses Chebyshev expansions, the coefficients of which
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) AIRY_BI_INT, the value of the function.
!
  implicit none

  real ( kind = 8 ) abint1(0:36)
  real ( kind = 8 ) abint2(0:37)
  real ( kind = 8 ) abint3(0:37)
  real ( kind = 8 ) abint4(0:20)
  real ( kind = 8 ) abint5(0:20)
  real ( kind = 8 ) airy_bi_int
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  integer, parameter :: nterm1 = 33
  integer, parameter :: nterm2 = 30
  integer, parameter :: nterm3 = 34
  integer nterm4
  integer nterm5
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) arg,birzer,f1,f2,nine,ninhun, &
       onept5,piby4,rt2b3p,sixten,seven,t,temp, &
       thr644,xlow1,xhigh1,xmax,xneg1, &
       z
  data abint1(0)/  0.38683352445038543350d0/
  data abint1(1)/ -0.8823213550888908821d-1/
  data abint1(2)/  0.21463937440355429239d0/
  data abint1(3)/ -0.4205347375891315126d-1/
  data abint1(4)/  0.5932422547496086771d-1/
  data abint1(5)/ -0.840787081124270210d-2/
  data abint1(6)/  0.871824772778487955d-2/
  data abint1(7)/ -0.12191600199613455d-3/
  data abint1(8)/  0.44024821786023234d-3/
  data abint1(9)/  0.27894686666386678d-3/
  data abint1(10)/-0.7052804689785537d-4/
  data abint1(11)/ 0.5901080066770100d-4/
  data abint1(12)/-0.1370862587982142d-4/
  data abint1(13)/ 0.505962573749073d-5/
  data abint1(14)/-0.51598837766735d-6/
  data abint1(15)/ 0.397511312349d-8/
  data abint1(16)/ 0.9524985978055d-7/
  data abint1(17)/-0.3681435887321d-7/
  data abint1(18)/ 0.1248391688136d-7/
  data abint1(19)/-0.249097619137d-8/
  data abint1(20)/ 0.31775245551d-9/
  data abint1(21)/ 0.5434365270d-10/
  data abint1(22)/-0.4024566915d-10/
  data abint1(23)/ 0.1393855527d-10/
  data abint1(24)/-0.303817509d-11/
  data abint1(25)/ 0.40809511d-12/
  data abint1(26)/ 0.1634116d-13/
  data abint1(27)/-0.2683809d-13/
  data abint1(28)/ 0.896641d-14/
  data abint1(29)/-0.183089d-14/
  data abint1(30)/ 0.21333d-15/
  data abint1(31)/ 0.1108d-16/
  data abint1(32)/-0.1276d-16/
  data abint1(33)/ 0.363d-17/
  data abint1(34)/-0.62d-18/
  data abint1(35)/ 0.5d-19/
  data abint1(36)/ 0.1d-19/
  data abint2(0)/  2.04122078602516135181d0/
  data abint2(1)/  0.2124133918621221230d-1/
  data abint2(2)/  0.66617599766706276d-3/
  data abint2(3)/  0.3842047982808254d-4/
  data abint2(4)/  0.362310366020439d-5/
  data abint2(5)/  0.50351990115074d-6/
  data abint2(6)/  0.7961648702253d-7/
  data abint2(7)/  0.717808442336d-8/
  data abint2(8)/ -0.267770159104d-8/
  data abint2(9)/ -0.168489514699d-8/
  data abint2(10)/-0.36811757255d-9/
  data abint2(11)/ 0.4757128727d-10/
  data abint2(12)/ 0.5263621945d-10/
  data abint2(13)/ 0.778973500d-11/
  data abint2(14)/-0.460546143d-11/
  data abint2(15)/-0.183433736d-11/
  data abint2(16)/ 0.32191249d-12/
  data abint2(17)/ 0.29352060d-12/
  data abint2(18)/-0.1657935d-13/
  data abint2(19)/-0.4483808d-13/
  data abint2(20)/ 0.27907d-15/
  data abint2(21)/ 0.711921d-14/
  data abint2(22)/-0.1042d-16/
  data abint2(23)/-0.119591d-14/
  data abint2(24)/ 0.4606d-16/
  data abint2(25)/ 0.20884d-15/
  data abint2(26)/-0.2416d-16/
  data abint2(27)/-0.3638d-16/
  data abint2(28)/ 0.863d-17/
  data abint2(29)/ 0.591d-17/
  data abint2(30)/-0.256d-17/
  data abint2(31)/-0.77d-18/
  data abint2(32)/ 0.66d-18/
  data abint2(33)/ 0.3d-19/
  data abint2(34)/-0.15d-18/
  data abint2(35)/ 0.2d-19/
  data abint2(36)/ 0.3d-19/
  data abint2(37)/-0.1d-19/
  data abint3(0)/  0.31076961598640349251d0/
  data abint3(1)/ -0.27528845887452542718d0/
  data abint3(2)/  0.17355965706136543928d0/
  data abint3(3)/ -0.5544017909492843130d-1/
  data abint3(4)/ -0.2251265478295950941d-1/
  data abint3(5)/  0.4107347447812521894d-1/
  data abint3(6)/  0.984761275464262480d-2/
  data abint3(7)/ -0.1555618141666041932d-1/
  data abint3(8)/ -0.560871870730279234d-2/
  data abint3(9)/  0.246017783322230475d-2/
  data abint3(10)/ 0.165740392292336978d-2/
  data abint3(11)/-0.3277587501435402d-4/
  data abint3(12)/-0.24434680860514925d-3/
  data abint3(13)/-0.5035305196152321d-4/
  data abint3(14)/ 0.1630264722247854d-4/
  data abint3(15)/ 0.851914057780934d-5/
  data abint3(16)/ 0.29790363004664d-6/
  data abint3(17)/-0.64389707896401d-6/
  data abint3(18)/-0.15046988145803d-6/
  data abint3(19)/ 0.1587013535823d-7/
  data abint3(20)/ 0.1276766299622d-7/
  data abint3(21)/ 0.140578534199d-8/
  data abint3(22)/-0.46564739741d-9/
  data abint3(23)/-0.15682748791d-9/
  data abint3(24)/-0.403893560d-11/
  data abint3(25)/ 0.666708192d-11/
  data abint3(26)/ 0.128869380d-11/
  data abint3(27)/-0.6968663d-13/
  data abint3(28)/-0.6254319d-13/
  data abint3(29)/-0.718392d-14/
  data abint3(30)/ 0.115296d-14/
  data abint3(31)/ 0.42276d-15/
  data abint3(32)/ 0.2493d-16/
  data abint3(33)/-0.971d-17/
  data abint3(34)/-0.216d-17/
  data abint3(35)/-0.2d-19/
  data abint3(36)/ 0.6d-19/
  data abint3(37)/ 0.1d-19/
  data abint4(0)/  1.99507959313352047614d0/
  data abint4(1)/ -0.273736375970692738d-2/
  data abint4(2)/ -0.30897113081285850d-3/
  data abint4(3)/ -0.3550101982798577d-4/
  data abint4(4)/ -0.412179271520133d-5/
  data abint4(5)/ -0.48235892316833d-6/
  data abint4(6)/ -0.5678730727927d-7/
  data abint4(7)/ -0.671874810365d-8/
  data abint4(8)/ -0.79811649857d-9/
  data abint4(9)/ -0.9514271478d-10/
  data abint4(10)/-0.1137468966d-10/
  data abint4(11)/-0.136359969d-11/
  data abint4(12)/-0.16381418d-12/
  data abint4(13)/-0.1972575d-13/
  data abint4(14)/-0.237844d-14/
  data abint4(15)/-0.28752d-15/
  data abint4(16)/-0.3475d-16/
  data abint4(17)/-0.422d-17/
  data abint4(18)/-0.51d-18/
  data abint4(19)/-0.6d-19/
  data abint4(20)/-0.1d-19/
  data abint5(0)/  1.12672081961782566017d0/
  data abint5(1)/ -0.671405567525561198d-2/
  data abint5(2)/ -0.69812918017832969d-3/
  data abint5(3)/ -0.7561689886425276d-4/
  data abint5(4)/ -0.834985574510207d-5/
  data abint5(5)/ -0.93630298232480d-6/
  data abint5(6)/ -0.10608556296250d-6/
  data abint5(7)/ -0.1213128916741d-7/
  data abint5(8)/ -0.139631129765d-8/
  data abint5(9)/ -0.16178918054d-9/
  data abint5(10)/-0.1882307907d-10/
  data abint5(11)/-0.220272985d-11/
  data abint5(12)/-0.25816189d-12/
  data abint5(13)/-0.3047964d-13/
  data abint5(14)/-0.358370d-14/
  data abint5(15)/-0.42831d-15/
  data abint5(16)/-0.4993d-16/
  data abint5(17)/-0.617d-17/
  data abint5(18)/-0.68d-18/
  data abint5(19)/-0.10d-18/
  data abint5(20)/-0.1d-19/
  data onept5/ 1.5d0 /
  data seven/ 7.0d0 /
  data nine,sixten/ 9.0d0 , 16.0d0 /
  data ninhun,thr644/900.0d0 , 3644.0d0 /
  data piby4/0.78539816339744830962d0/
  data rt2b3p/0.46065886596178063902d0/
  data birzer/0.61492662744600073515d0/
!
!   Machine-dependent parameters (suitable for IEEE machines)
!
  data nterm4,nterm5/17,17/
  data xlow1,xhigh1/2.22044604d-16,104.587632d0/
  data xneg1,xmax/-2.727134d10,1.79d308/

  x = xvalue

  if ( x < xneg1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AIRY_BI_INT - Warning!'
    write ( *, '(a)' ) '  Argument is too negative for accurate computation.'
    airy_bi_int = zero
  else if ( x < -seven ) then
    z = - ( x + x ) * sqrt ( -x ) / three
    arg = z + piby4
    temp = nine * z * z
    t = ( thr644 - temp ) / ( ninhun + temp )
    f1 = cheval ( nterm4, abint4, t ) * sin ( arg )
    f2 = cheval ( nterm5, abint5, t ) * cos ( arg ) / z
    airy_bi_int = ( f2 - f1 ) * rt2b3p / sqrt ( z )
  else if ( x <= -xlow1 ) then
    t = - ( x + x ) / seven - one
    airy_bi_int = x * cheval ( nterm3, abint3, t )
  else if ( x < xlow1 ) then
    airy_bi_int = birzer * x
  else if ( x <= eight ) then
    t = x / four - one
    airy_bi_int = x * exp ( onept5 * x ) * cheval ( nterm1, abint1, t )
  else if ( x <= xhigh1 ) then
    t = sixten * sqrt ( eight / x ) / x - one
    z = ( x + x ) * sqrt ( x ) / three
    temp = rt2b3p * cheval ( nterm2, abint2, t ) / sqrt ( z )
    temp = z + log ( temp )
    airy_bi_int = exp ( temp )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AIRY_BI_INT - Warning!'
    write ( *, '(a)' ) '  Argument is too large for accurate computation.'
    airy_bi_int = xmax
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
!      AIRY_BI_INT(x) = Integral ( 0 <= t <= x ) Bi(t) dt
!
!    The data was reported by McLeod.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function airy_gi ( xvalue )

!*****************************************************************************80
!
!! AIRY_GI computes the modified Airy function Gi(x).
!
!  Discussion:
!
!    The function is defined by:
!
!      AIRY_GI(x) = Integral ( 0 <= t < infinity ) sin ( x*t+t^3/3) dt / pi
!
!    The approximation uses Chebyshev expansions with the coefficients
!    given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) AIRY_GI, the value of the function.
!
  implicit none

  real ( kind = 8 ) airy_gi
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: four = 4.0D+00
  integer, parameter :: nterm1 = 28
  integer, parameter :: nterm2 = 23
  integer, parameter :: nterm3 = 39
  integer nterm4
  integer nterm5
  integer nterm6
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) argip1(0:30),argip2(0:29),argin1(0:42), &
       arbin1(0:10),arbin2(0:11),arhin1(0:15), &
       bi,cheb1,cheb2,cosz,five,five14, &
       gizero,minate,nine,onebpi,one76,one024,piby4, &
       rtpiin,seven,seven2,sinz,t,temp,twelhu,twent8, &
       xcube,xhigh1,xhigh2,xhigh3,xlow1,xminus, &
       zeta
  data argip1(0)/  0.26585770795022745082d0/
  data argip1(1)/ -0.10500333097501922907d0/
  data argip1(2)/  0.841347475328454492d-2/
  data argip1(3)/  0.2021067387813439541d-1/
  data argip1(4)/ -0.1559576113863552234d-1/
  data argip1(5)/  0.564342939043256481d-2/
  data argip1(6)/ -0.59776844826655809d-3/
  data argip1(7)/ -0.42833850264867728d-3/
  data argip1(8)/  0.22605662380909027d-3/
  data argip1(9)/ -0.3608332945592260d-4/
  data argip1(10)/-0.785518988788901d-5/
  data argip1(11)/ 0.473252480746370d-5/
  data argip1(12)/-0.59743513977694d-6/
  data argip1(13)/-0.15917609165602d-6/
  data argip1(14)/ 0.6336129065570d-7/
  data argip1(15)/-0.276090232648d-8/
  data argip1(16)/-0.256064154085d-8/
  data argip1(17)/ 0.47798676856d-9/
  data argip1(18)/ 0.4488131863d-10/
  data argip1(19)/-0.2346508882d-10/
  data argip1(20)/ 0.76839085d-12/
  data argip1(21)/ 0.73227985d-12/
  data argip1(22)/-0.8513687d-13/
  data argip1(23)/-0.1630201d-13/
  data argip1(24)/ 0.356769d-14/
  data argip1(25)/ 0.25001d-15/
  data argip1(26)/-0.10859d-15/
  data argip1(27)/-0.158d-17/
  data argip1(28)/ 0.275d-17/
  data argip1(29)/-0.5d-19/
  data argip1(30)/-0.6d-19/
  data argip2(0)/  2.00473712275801486391d0/
  data argip2(1)/  0.294184139364406724d-2/
  data argip2(2)/  0.71369249006340167d-3/
  data argip2(3)/  0.17526563430502267d-3/
  data argip2(4)/  0.4359182094029882d-4/
  data argip2(5)/  0.1092626947604307d-4/
  data argip2(6)/  0.272382418399029d-5/
  data argip2(7)/  0.66230900947687d-6/
  data argip2(8)/  0.15425323370315d-6/
  data argip2(9)/  0.3418465242306d-7/
  data argip2(10)/ 0.728157724894d-8/
  data argip2(11)/ 0.151588525452d-8/
  data argip2(12)/ 0.30940048039d-9/
  data argip2(13)/ 0.6149672614d-10/
  data argip2(14)/ 0.1202877045d-10/
  data argip2(15)/ 0.233690586d-11/
  data argip2(16)/ 0.43778068d-12/
  data argip2(17)/ 0.7996447d-13/
  data argip2(18)/ 0.1494075d-13/
  data argip2(19)/ 0.246790d-14/
  data argip2(20)/ 0.37672d-15/
  data argip2(21)/ 0.7701d-16/
  data argip2(22)/ 0.354d-17/
  data argip2(23)/-0.49d-18/
  data argip2(24)/ 0.62d-18/
  data argip2(25)/-0.40d-18/
  data argip2(26)/-0.1d-19/
  data argip2(27)/ 0.2d-19/
  data argip2(28)/-0.3d-19/
  data argip2(29)/ 0.1d-19/
  data argin1(0)/ -0.20118965056732089130d0/
  data argin1(1)/ -0.7244175303324530499d-1/
  data argin1(2)/  0.4505018923894780120d-1/
  data argin1(3)/ -0.24221371122078791099d0/
  data argin1(4)/  0.2717884964361678294d-1/
  data argin1(5)/ -0.5729321004818179697d-1/
  data argin1(6)/ -0.18382107860337763587d0/
  data argin1(7)/  0.7751546082149475511d-1/
  data argin1(8)/  0.18386564733927560387d0/
  data argin1(9)/  0.2921504250185567173d-1/
  data argin1(10)/-0.6142294846788018811d-1/
  data argin1(11)/-0.2999312505794616238d-1/
  data argin1(12)/ 0.585937118327706636d-2/
  data argin1(13)/ 0.822221658497402529d-2/
  data argin1(14)/ 0.132579817166846893d-2/
  data argin1(15)/-0.96248310766565126d-3/
  data argin1(16)/-0.45065515998211807d-3/
  data argin1(17)/ 0.772423474325474d-5/
  data argin1(18)/ 0.5481874134758052d-4/
  data argin1(19)/ 0.1245898039742876d-4/
  data argin1(20)/-0.246196891092083d-5/
  data argin1(21)/-0.169154183545285d-5/
  data argin1(22)/-0.16769153169442d-6/
  data argin1(23)/ 0.9636509337672d-7/
  data argin1(24)/ 0.3253314928030d-7/
  data argin1(25)/ 0.5091804231d-10/
  data argin1(26)/-0.209180453553d-8/
  data argin1(27)/-0.41237387870d-9/
  data argin1(28)/ 0.4163338253d-10/
  data argin1(29)/ 0.3032532117d-10/
  data argin1(30)/ 0.340580529d-11/
  data argin1(31)/-0.88444592d-12/
  data argin1(32)/-0.31639612d-12/
  data argin1(33)/-0.1505076d-13/
  data argin1(34)/ 0.1104148d-13/
  data argin1(35)/ 0.246508d-14/
  data argin1(36)/-0.3107d-16/
  data argin1(37)/-0.9851d-16/
  data argin1(38)/-0.1453d-16/
  data argin1(39)/ 0.118d-17/
  data argin1(40)/ 0.67d-18/
  data argin1(41)/ 0.6d-19/
  data argin1(42)/-0.1d-19/
  data arbin1/1.99983763583586155980d0, &
             -0.8104660923669418d-4, &
              0.13475665984689d-6, &
             -0.70855847143d-9, &
              0.748184187d-11, &
             -0.12902774d-12, &
              0.322504d-14, &
             -0.10809d-15, &
              0.460d-17, &
             -0.24d-18, &
              0.1d-19/
  data arbin2/0.13872356453879120276d0, &
             -0.8239286225558228d-4, &
              0.26720919509866d-6, &
             -0.207423685368d-8, &
              0.2873392593d-10, &
             -0.60873521d-12, &
              0.1792489d-13, &
             -0.68760d-15, &
              0.3280d-16, &
             -0.188d-17, &
              0.13d-18, &
             -0.1d-19/
  data arhin1/1.99647720399779650525d0, &
             -0.187563779407173213d-2, &
             -0.12186470897787339d-3, &
             -0.814021609659287d-5, &
             -0.55050925953537d-6, &
             -0.3763008043303d-7, &
             -0.258858362365d-8, &
             -0.17931829265d-9, &
             -0.1245916873d-10, &
             -0.87171247d-12, &
             -0.6084943d-13, &
             -0.431178d-14, &
             -0.29787d-15, &
             -0.2210d-16, &
             -0.136d-17, &
             -0.14d-18/
  data five,seven,minate/ 5.0d0, 7.0d0 , -8.0d0 /
  data nine,twent8,seven2/ 9.0d0, 28.0d0 , 72.0d0 /
  data one76,five14/ 176.0d0 , 514.0d0 /
  data one024,twelhu/ 1024.0d0, 1200.0d0 /
  data gizero/0.20497554248200024505d0/
  data onebpi/0.31830988618379067154d0/
  data piby4/0.78539816339744830962d0/
  data rtpiin/0.56418958354775628695d0/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data nterm4,nterm5,nterm6/9,10,14/
  data xlow1,xhigh1/2.22045d-16,208063.8307d0/
  data xhigh2,xhigh3/0.14274d308,-2097152.0d0/

  x = xvalue

  if ( x < -xhigh1 * xhigh1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AIRY_GI - Fatal error!'
    write ( *, '(a)' ) '  Argument too negative for accurate computation.'
    airy_gi = zero
  else if ( x <= xhigh3 ) then
    xminus = -x
    t = xminus * sqrt ( xminus )
    zeta = ( t + t ) / three
    temp = rtpiin / sqrt ( sqrt ( xminus ) )
    cosz = cos ( zeta + piby4 )
    sinz = sin ( zeta + piby4 ) / zeta
    xcube = x * x * x
    bi = ( cosz + sinz * five / seven2 ) * temp
    t = ( xcube + twelhu ) / ( one76 - xcube )
    airy_gi = bi + cheval ( nterm6, arhin1, t ) * onebpi / x
  else if ( x < minate ) then
    xminus = -x
    t = xminus * sqrt ( xminus )
    zeta = ( t + t ) / three
    temp = rtpiin / sqrt ( sqrt ( xminus ) )
    cosz = cos ( zeta + piby4 )
    sinz = sin ( zeta + piby4 ) / zeta
    xcube = x * x * x
    t = - ( one024 / ( xcube ) + one )
    cheb1 = cheval ( nterm4, arbin1, t )
    cheb2 = cheval ( nterm5, arbin2, t )
    bi = ( cosz * cheb1 + sinz * cheb2 ) * temp
    t = ( xcube + twelhu ) / ( one76 - xcube )
    airy_gi = bi + cheval ( nterm6, arhin1, t ) * onebpi / x
  else if ( x <= -xlow1 ) then
    t = -( x + four ) / four
    airy_gi = cheval ( nterm3, argin1, t )
  else if ( x < xlow1 ) then
    airy_gi = gizero
  else if ( x <= seven ) then
    t = ( nine * x - twent8 ) / ( x + twent8 )
    airy_gi = cheval ( nterm1, argip1, t )
  else if ( x <= xhigh1 ) then
    xcube = x * x * x
    t = ( twelhu - xcube ) / ( five14 + xcube )
    airy_gi = onebpi * cheval ( nterm2, argip2, t ) / x
  else if ( x <= xhigh2 ) then
    airy_gi = onebpi / x
  else
    airy_gi = zero
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
!      AIRY_GI(x) = Integral ( 0 <= t < infinity ) sin ( x*t+t^3/3) dt / pi
!
!    The data was reported by McLeod.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function airy_hi ( xvalue )

!*****************************************************************************80
!
!! AIRY_HI computes the modified Airy function Hi(x).
!
!  Discussion:
!
!    The function is defined by:
!
!      AIRY_HI(x) = Integral ( 0 <= t < infinity ) exp(x*t-t^3/3) dt / pi
!
!    The approximation uses Chebyshev expansions with the coefficients
!    given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) AIRY_HI, the value of the function.
!
  implicit none

  real ( kind = 8 ) airy_hi
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: four = 4.0D+00
  integer, parameter :: nterm1 = 29
  integer, parameter :: nterm2 = 17
  integer, parameter :: nterm3 = 22
  integer nterm4
  integer nterm5
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) arhip(0:31),arbip(0:23),argip1(0:29), &
       arhin1(0:21),arhin2(0:15), &
       bi,five14,gi,hizero,lnrtpi, &
       minate,onebpi,one76,seven,t,temp, &
       thre43,twelhu,twelve,xcube, &
       xhigh1,xlow1,xmax,xneg1,xneg2, &
       zeta
  data arhip(0)/ 1.24013562561762831114d0/
  data arhip(1)/ 0.64856341973926535804d0/
  data arhip(2)/ 0.55236252592114903246d0/
  data arhip(3)/ 0.20975122073857566794d0/
  data arhip(4)/ 0.12025669118052373568d0/
  data arhip(5)/ 0.3768224931095393785d-1/
  data arhip(6)/ 0.1651088671548071651d-1/
  data arhip(7)/ 0.455922755211570993d-2/
  data arhip(8)/ 0.161828480477635013d-2/
  data arhip(9)/ 0.40841282508126663d-3/
  data arhip(10)/0.12196479721394051d-3/
  data arhip(11)/0.2865064098657610d-4/
  data arhip(12)/0.742221556424344d-5/
  data arhip(13)/0.163536231932831d-5/
  data arhip(14)/0.37713908188749d-6/
  data arhip(15)/0.7815800336008d-7/
  data arhip(16)/0.1638447121370d-7/
  data arhip(17)/0.319857665992d-8/
  data arhip(18)/0.61933905307d-9/
  data arhip(19)/0.11411161191d-9/
  data arhip(20)/0.2064923454d-10/
  data arhip(21)/0.360018664d-11/
  data arhip(22)/0.61401849d-12/
  data arhip(23)/0.10162125d-12/
  data arhip(24)/0.1643701d-13/
  data arhip(25)/0.259084d-14/
  data arhip(26)/0.39931d-15/
  data arhip(27)/0.6014d-16/
  data arhip(28)/0.886d-17/
  data arhip(29)/0.128d-17/
  data arhip(30)/0.18d-18/
  data arhip(31)/0.3d-19/
  data arbip(0)/  2.00582138209759064905d0/
  data arbip(1)/  0.294478449170441549d-2/
  data arbip(2)/  0.3489754514775355d-4/
  data arbip(3)/  0.83389733374343d-6/
  data arbip(4)/  0.3136215471813d-7/
  data arbip(5)/  0.167865306015d-8/
  data arbip(6)/  0.12217934059d-9/
  data arbip(7)/  0.1191584139d-10/
  data arbip(8)/  0.154142553d-11/
  data arbip(9)/  0.24844455d-12/
  data arbip(10)/ 0.4213012d-13/
  data arbip(11)/ 0.505293d-14/
  data arbip(12)/-0.60032d-15/
  data arbip(13)/-0.65474d-15/
  data arbip(14)/-0.22364d-15/
  data arbip(15)/-0.3015d-16/
  data arbip(16)/ 0.959d-17/
  data arbip(17)/ 0.616d-17/
  data arbip(18)/ 0.97d-18/
  data arbip(19)/-0.37d-18/
  data arbip(20)/-0.21d-18/
  data arbip(21)/-0.1d-19/
  data arbip(22)/ 0.2d-19/
  data arbip(23)/ 0.1d-19/
  data argip1(0)/  2.00473712275801486391d0/
  data argip1(1)/  0.294184139364406724d-2/
  data argip1(2)/  0.71369249006340167d-3/
  data argip1(3)/  0.17526563430502267d-3/
  data argip1(4)/  0.4359182094029882d-4/
  data argip1(5)/  0.1092626947604307d-4/
  data argip1(6)/  0.272382418399029d-5/
  data argip1(7)/  0.66230900947687d-6/
  data argip1(8)/  0.15425323370315d-6/
  data argip1(9)/  0.3418465242306d-7/
  data argip1(10)/ 0.728157724894d-8/
  data argip1(11)/ 0.151588525452d-8/
  data argip1(12)/ 0.30940048039d-9/
  data argip1(13)/ 0.6149672614d-10/
  data argip1(14)/ 0.1202877045d-10/
  data argip1(15)/ 0.233690586d-11/
  data argip1(16)/ 0.43778068d-12/
  data argip1(17)/ 0.7996447d-13/
  data argip1(18)/ 0.1494075d-13/
  data argip1(19)/ 0.246790d-14/
  data argip1(20)/ 0.37672d-15/
  data argip1(21)/ 0.7701d-16/
  data argip1(22)/ 0.354d-17/
  data argip1(23)/-0.49d-18/
  data argip1(24)/ 0.62d-18/
  data argip1(25)/-0.40d-18/
  data argip1(26)/-0.1d-19/
  data argip1(27)/ 0.2d-19/
  data argip1(28)/-0.3d-19/
  data argip1(29)/ 0.1d-19/
  data arhin1(0)/  0.31481017206423404116d0/
  data arhin1(1)/ -0.16414499216588964341d0/
  data arhin1(2)/  0.6176651597730913071d-1/
  data arhin1(3)/ -0.1971881185935933028d-1/
  data arhin1(4)/  0.536902830023331343d-2/
  data arhin1(5)/ -0.124977068439663038d-2/
  data arhin1(6)/  0.24835515596994933d-3/
  data arhin1(7)/ -0.4187024096746630d-4/
  data arhin1(8)/  0.590945437979124d-5/
  data arhin1(9)/ -0.68063541184345d-6/
  data arhin1(10)/ 0.6072897629164d-7/
  data arhin1(11)/-0.367130349242d-8/
  data arhin1(12)/ 0.7078017552d-10/
  data arhin1(13)/ 0.1187894334d-10/
  data arhin1(14)/-0.120898723d-11/
  data arhin1(15)/ 0.1189656d-13/
  data arhin1(16)/ 0.594128d-14/
  data arhin1(17)/-0.32257d-15/
  data arhin1(18)/-0.2290d-16/
  data arhin1(19)/ 0.253d-17/
  data arhin1(20)/ 0.9d-19/
  data arhin1(21)/-0.2d-19/
  data arhin2/1.99647720399779650525d0, &
             -0.187563779407173213d-2, &
             -0.12186470897787339d-3, &
             -0.814021609659287d-5, &
             -0.55050925953537d-6, &
             -0.3763008043303d-7, &
             -0.258858362365d-8, &
             -0.17931829265d-9, &
             -0.1245916873d-10, &
             -0.87171247d-12, &
             -0.6084943d-13, &
             -0.431178d-14, &
             -0.29787d-15, &
             -0.2210d-16, &
             -0.136d-17, &
             -0.14d-18/
  data seven/ 7.0d0 /
  data minate,twelve,one76/ -8.0d0 , 12.0d0 , 176.0d0 /
  data thre43,five14,twelhu/ 343.0d0, 514.0d0, 1200.0d0 /
  data hizero/0.40995108496400049010d0/
  data lnrtpi/0.57236494292470008707d0/
  data onebpi/0.31830988618379067154d0/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data nterm4,nterm5/19,14/
  data xlow1,xhigh1/2.220446d-16,104.4175d0/
  data xneg1,xneg2,xmax/-0.14274d308,-208063.831d0,1.79d308/

  x = xvalue

  if ( xhigh1 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AIRY_HI - Fatal error!'
    write ( *, '(a)' ) '  Argument too large.'
    airy_hi = xmax
    return
  end if
!
!  Code for x < 0.0
!
  if ( x < zero ) then

    if ( x < minate ) then

      if ( x < xneg1 ) then
        airy_hi = zero
      else
        if ( x < xneg2 ) then
          temp = one
          airy_hi = - temp * onebpi / x
        else
          xcube = x * x * x
          t = ( xcube + twelhu ) / ( one76 - xcube )
          temp = cheval ( nterm5, arhin2, t )
          airy_hi = - temp * onebpi / x
        end if

      end if
    else
      if ( -xlow1 < x ) then
        airy_hi = hizero
      else
        t = ( four * x + twelve ) / ( x - twelve )
        airy_hi = cheval ( nterm4, arhin1, t )
      end if

    end if
!
!   Code for x >= 0.0
!
  else

    if ( x <= seven ) then
      if ( x < xlow1 ) then
        airy_hi = hizero
      else
        t = ( x + x ) / seven - one
        temp = ( x + x + x ) / two
        airy_hi = exp ( temp ) * cheval ( nterm1, arhip, t )
      end if
    else
      xcube = x * x * x
      temp = sqrt ( xcube )
      zeta = ( temp + temp ) / three
      t = two * ( sqrt ( thre43 / xcube ) ) - one
      temp = cheval ( nterm2, arbip, t )
      temp = zeta + log ( temp ) - log ( x ) / four - lnrtpi
      bi = exp ( temp )
      t = ( twelhu - xcube ) / ( xcube + five14 )
      gi = cheval ( nterm3, argip1, t ) * onebpi / x
      airy_hi = bi - gi
    end if

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
!      AIRY_HI(x) =
!        Integral ( 0 <= t < infinity ) exp ( x * t - t^3 / 3 ) dt / pi
!
!    The data was reported by McLeod.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function arctan_int ( xvalue )

!*****************************************************************************80
!
!! ARCTAN_INT calculates the inverse tangent integral.
!
!  Discussion:
!
!    The function is defined by:
!
!      ARCTAN_INT(x) = Integral ( 0 <= t <= x ) arctan ( t ) / t dt
!
!    The approximation uses Chebyshev series with the coefficients
!    given to an accuracy of 20D.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    24 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) ARCTAN_INT, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 0:22 ) :: atnina = (/ &
     1.91040361296235937512d0, &
    -0.4176351437656746940d-1, &
     0.275392550786367434d-2, &
    -0.25051809526248881d-3, &
     0.2666981285121171d-4, &
    -0.311890514107001d-5, &
     0.38833853132249d-6, &
    -0.5057274584964d-7, &
     0.681225282949d-8, &
    -0.94212561654d-9, &
     0.13307878816d-9, &
    -0.1912678075d-10, &
     0.278912620d-11, &
    -0.41174820d-12, &
     0.6142987d-13, &
    -0.924929d-14, &
     0.140387d-14, &
    -0.21460d-15, &
     0.3301d-16, &
    -0.511d-17, &
     0.79d-18, &
    -0.12d-18, &
     0.2d-19 /)
  real ( kind = 8 ) arctan_int
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ind
  integer, parameter :: nterms = 19
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: twobpi = 0.63661977236758134308D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xlow = 7.4505806D-09
  real ( kind = 8 ), parameter :: xupper = 4.5036D+15
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  ind = 1
  x = xvalue

  if ( x < zero ) then
    x = -x
    ind = -1
  end if

  if ( x < xlow ) then
    arctan_int = x
  else if ( x <= one ) then
    t = x * x
    t =  ( t - half ) + ( t - half )
    arctan_int = x * cheval ( nterms, atnina, t )
  else if ( x <= xupper ) then
    t = one / ( x * x )
    t = ( t - half ) + ( t - half )
    arctan_int = log ( x ) / twobpi + cheval ( nterms, atnina, t ) / x
  else
    arctan_int = log ( x ) / twobpi
  end if

  if ( ind < 0 ) then
    arctan_int = -arctan_int
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
!      ARCTAN_INT(x) = Integral ( 0 <= t <= x ) arctan ( t ) / t dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function bessel_i0_int ( xvalue )

!*****************************************************************************80
!
!! BESSEL_I0_INT computes the integral of the modified Bessel function I0(X).
!
!  Discussion:
!
!    The function is defined by:
!
!      I0_INT(x) = Integral ( 0 <= t <= x ) I0(t) dt
!
!    The program uses Chebyshev expansions, the coefficients of
!    which are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) BESSEL_I0_INT, the value of the function.
!
  implicit none

  real ( kind = 8 ) ari01(0:28)
  real ( kind = 8 ) ari0a(0:33)
  real ( kind = 8 ), parameter :: ateen = 18.0D+00
  real ( kind = 8 ) bessel_i0_int
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ind
  real ( kind = 8 ), parameter :: lnr2pi = 0.91893853320467274178D+00
  integer, parameter :: nterm1 = 25
  integer, parameter :: nterm2 = 27
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: thirt6 = 36.0D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xhigh = 713.758339D+00
  real ( kind = 8 ), parameter :: xlow = 0.5161914D-07
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  data ari01(0)/  0.41227906926781516801d0/
  data ari01(1)/ -0.34336345150081519562d0/
  data ari01(2)/  0.22667588715751242585d0/
  data ari01(3)/ -0.12608164718742260032d0/
  data ari01(4)/  0.6012484628777990271d-1/
  data ari01(5)/ -0.2480120462913358248d-1/
  data ari01(6)/  0.892773389565563897d-2/
  data ari01(7)/ -0.283253729936696605d-2/
  data ari01(8)/  0.79891339041712994d-3/
  data ari01(9)/ -0.20053933660964890d-3/
  data ari01(10)/ 0.4416816783014313d-4/
  data ari01(11)/-0.822377042246068d-5/
  data ari01(12)/ 0.120059794219015d-5/
  data ari01(13)/-0.11350865004889d-6/
  data ari01(14)/ 0.69606014466d-9/
  data ari01(15)/ 0.180622772836d-8/
  data ari01(16)/-0.26039481370d-9/
  data ari01(17)/-0.166188103d-11/
  data ari01(18)/ 0.510500232d-11/
  data ari01(19)/-0.41515879d-12/
  data ari01(20)/-0.7368138d-13/
  data ari01(21)/ 0.1279323d-13/
  data ari01(22)/ 0.103247d-14/
  data ari01(23)/-0.30379d-15/
  data ari01(24)/-0.1789d-16/
  data ari01(25)/ 0.673d-17/
  data ari01(26)/ 0.44d-18/
  data ari01(27)/-0.14d-18/
  data ari01(28)/-0.1d-19/
  data ari0a(0)/  2.03739654571143287070d0/
  data ari0a(1)/  0.1917631647503310248d-1/
  data ari0a(2)/  0.49923334519288147d-3/
  data ari0a(3)/  0.2263187103659815d-4/
  data ari0a(4)/  0.158682108285561d-5/
  data ari0a(5)/  0.16507855636318d-6/
  data ari0a(6)/  0.2385058373640d-7/
  data ari0a(7)/  0.392985182304d-8/
  data ari0a(8)/  0.46042714199d-9/
  data ari0a(9)/ -0.7072558172d-10/
  data ari0a(10)/-0.6747183961d-10/
  data ari0a(11)/-0.2026962001d-10/
  data ari0a(12)/-0.87320338d-12/
  data ari0a(13)/ 0.175520014d-11/
  data ari0a(14)/ 0.60383944d-12/
  data ari0a(15)/-0.3977983d-13/
  data ari0a(16)/-0.8049048d-13/
  data ari0a(17)/-0.1158955d-13/
  data ari0a(18)/ 0.827318d-14/
  data ari0a(19)/ 0.282290d-14/
  data ari0a(20)/-0.77667d-15/
  data ari0a(21)/-0.48731d-15/
  data ari0a(22)/ 0.7279d-16/
  data ari0a(23)/ 0.7873d-16/
  data ari0a(24)/-0.785d-17/
  data ari0a(25)/-0.1281d-16/
  data ari0a(26)/ 0.121d-17/
  data ari0a(27)/ 0.214d-17/
  data ari0a(28)/-0.27d-18/
  data ari0a(29)/-0.36d-18/
  data ari0a(30)/ 0.7d-19/
  data ari0a(31)/ 0.6d-19/
  data ari0a(32)/-0.2d-19/
  data ari0a(33)/-0.1d-19/

  ind = 1
  x = xvalue

  if ( xvalue < zero ) then
    ind = -1
    x = -x
  end if

  if ( x < xlow ) then
    bessel_i0_int = x
  else if ( x <= ateen ) then
    t = ( three * x - ateen ) / ( x + ateen )
    bessel_i0_int = x * exp ( x ) * cheval ( nterm1, ari01, t )
  else if ( x <= xhigh ) then
    t = ( thirt6 / x - half ) - half
    temp = x - half * log ( x ) - lnr2pi + log ( cheval ( nterm2, ari0a, t ))
    bessel_i0_int = exp ( temp )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BESSEL_I0_INT - Fatal error!'
    write ( *, '(a)' ) '  Argument magnitude too large.'
    bessel_i0_int = exp ( xhigh - lnr2pi - half * log ( xhigh ) )
  end if

  if ( ind == -1 ) then
    bessel_i0_int = -bessel_i0_int
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
!      I0_INT(x) = Integral ( 0 <= t <= x ) I0(t) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function bessel_j0_int ( xvalue )

!*****************************************************************************80
!
!! BESSEL_J0_INT calculates the integral of the Bessel function J0.
!
!  Discussion:
!
!    The function is defined by:
!
!      J0_INT(x) = Integral ( 0 <= t <= x ) J0(t) dt
!
!    The code uses Chebyshev expansions whose coefficients are
!    given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) BESSEL_J0_INT, the value of the function.
!
  implicit none

  real ( kind = 8 ) bessel_j0_int
  real ( kind = 8 ) cheval
  integer ind
  integer, parameter :: nterm1 = 22
  integer, parameter :: nterm2 = 18
  integer, parameter :: nterm3 = 16
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) arj01(0:23),arj0a1(0:21),arj0a2(0:18), &
       five12,one28,pib41,pib411,pib412, &
       pib42,rt2bpi,sixten,t,temp,xhigh,xlow, &
       xmpi4
  data sixten/ 16.0d0 /
  data one28,five12/ 128.0d0 , 512d0 /
  data rt2bpi/0.79788456080286535588d0/
  data pib411,pib412/ 201.0d0 , 256.0d0/
  data pib42/0.24191339744830961566d-3/
  data arj01(0)/  0.38179279321690173518d0/
  data arj01(1)/ -0.21275636350505321870d0/
  data arj01(2)/  0.16754213407215794187d0/
  data arj01(3)/ -0.12853209772196398954d0/
  data arj01(4)/  0.10114405455778847013d0/
  data arj01(5)/ -0.9100795343201568859d-1/
  data arj01(6)/  0.6401345264656873103d-1/
  data arj01(7)/ -0.3066963029926754312d-1/
  data arj01(8)/  0.1030836525325064201d-1/
  data arj01(9)/ -0.255670650399956918d-2/
  data arj01(10)/ 0.48832755805798304d-3/
  data arj01(11)/-0.7424935126036077d-4/
  data arj01(12)/ 0.922260563730861d-5/
  data arj01(13)/-0.95522828307083d-6/
  data arj01(14)/ 0.8388355845986d-7/
  data arj01(15)/-0.633184488858d-8/
  data arj01(16)/ 0.41560504221d-9/
  data arj01(17)/-0.2395529307d-10/
  data arj01(18)/ 0.122286885d-11/
  data arj01(19)/-0.5569711d-13/
  data arj01(20)/ 0.227820d-14/
  data arj01(21)/-0.8417d-16/
  data arj01(22)/ 0.282d-17/
  data arj01(23)/-0.9d-19/
  data arj0a1(0)/  1.24030133037518970827d0/
  data arj0a1(1)/ -0.478125353632280693d-2/
  data arj0a1(2)/  0.6613148891706678d-4/
  data arj0a1(3)/ -0.186042740486349d-5/
  data arj0a1(4)/  0.8362735565080d-7/
  data arj0a1(5)/ -0.525857036731d-8/
  data arj0a1(6)/  0.42606363251d-9/
  data arj0a1(7)/ -0.4211761024d-10/
  data arj0a1(8)/  0.488946426d-11/
  data arj0a1(9)/ -0.64834929d-12/
  data arj0a1(10)/ 0.9617234d-13/
  data arj0a1(11)/-0.1570367d-13/
  data arj0a1(12)/ 0.278712d-14/
  data arj0a1(13)/-0.53222d-15/
  data arj0a1(14)/ 0.10844d-15/
  data arj0a1(15)/-0.2342d-16/
  data arj0a1(16)/ 0.533d-17/
  data arj0a1(17)/-0.127d-17/
  data arj0a1(18)/ 0.32d-18/
  data arj0a1(19)/-0.8d-19/
  data arj0a1(20)/ 0.2d-19/
  data arj0a1(21)/-0.1d-19/
  data arj0a2(0)/  1.99616096301341675339d0/
  data arj0a2(1)/ -0.190379819246668161d-2/
  data arj0a2(2)/  0.1539710927044226d-4/
  data arj0a2(3)/ -0.31145088328103d-6/
  data arj0a2(4)/  0.1110850971321d-7/
  data arj0a2(5)/ -0.58666787123d-9/
  data arj0a2(6)/  0.4139926949d-10/
  data arj0a2(7)/ -0.365398763d-11/
  data arj0a2(8)/  0.38557568d-12/
  data arj0a2(9)/ -0.4709800d-13/
  data arj0a2(10)/ 0.650220d-14/
  data arj0a2(11)/-0.99624d-15/
  data arj0a2(12)/ 0.16700d-15/
  data arj0a2(13)/-0.3028d-16/
  data arj0a2(14)/ 0.589d-17/
  data arj0a2(15)/-0.122d-17/
  data arj0a2(16)/ 0.27d-18/
  data arj0a2(17)/-0.6d-19/
  data arj0a2(18)/ 0.1d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow,xhigh/3.650024d-8,9.0072d15/

  x = xvalue
  ind = 1

  if ( x < zero ) then
    x = -x
    ind = -1
  end if

  if ( x < xlow ) then
    bessel_j0_int = x
  else if ( x <= sixten ) then
    t = x * x / one28 - one
    bessel_j0_int = x * cheval ( nterm1, arj01, t )
  else if ( x <= xhigh ) then
    t = five12 / ( x * x ) - one
    pib41 = pib411 / pib412
    xmpi4 = ( x - pib41 ) - pib42
    temp = cos ( xmpi4 ) * cheval ( nterm2, arj0a1, t ) / x
    temp = temp - sin ( xmpi4) * cheval ( nterm3, arj0a2, t )
    bessel_j0_int = one - rt2bpi * temp / sqrt ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BESSEL_J0_INT - Fatal error!'
    write ( *, '(a)' ) '  Argument magnitude too large.'
    bessel_j0_int = one
  end if

  if ( ind == -1 ) then
    bessel_j0_int = -bessel_j0_int
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
!      J0_INT(x) = Integral ( 0 <= t <= x ) J0(t) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function bessel_k0_int ( xvalue )

!*****************************************************************************80
!
!! BESSEL_K0_INT calculates the integral of the modified Bessel function K0(X).
!
!  Discussion:
!
!    The function is defined by:
!
!      K0_INT(x) = Integral ( 0 <= t <= x ) K0(t) dt
!
!    The code uses Chebyshev expansions, whose coefficients are
!    given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) BESSEL_K0_INT, the value of the function.
!
  implicit none

  real ( kind = 8 ) bessel_k0_int
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer, parameter :: nterm1 = 14
  integer, parameter :: nterm2 = 14
  integer, parameter :: nterm3 = 23
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) ak0in1(0:15),ak0in2(0:15),ak0ina(0:27), &
       const1,const2,eightn,fval, &
       piby2,rt2bpi,t,temp,twelve, &
       xhigh,xlow
  data twelve,eightn / 12.0d0 , 18.0d0 /
  data const1/1.11593151565841244881d0/
  data const2/-0.11593151565841244881d0/
  data piby2/1.57079632679489661923d0/
  data rt2bpi/0.79788456080286535588d0/
  data ak0in1/16.79702714464710959477d0, &
               9.79134687676889407070d0, &
               2.80501316044337939300d0, &
               0.45615620531888502068d0, &
               0.4716224457074760784d-1, &
               0.335265148269698289d-2, &
               0.17335181193874727d-3, &
               0.679951889364702d-5, &
               0.20900268359924d-6, &
               0.516603846976d-8, &
               0.10485708331d-9, &
               0.177829320d-11, &
               0.2556844d-13, &
               0.31557d-15, &
               0.338d-17, &
               0.3d-19/
  data ak0in2/10.76266558227809174077d0, &
               5.62333479849997511550d0, &
               1.43543664879290867158d0, &
               0.21250410143743896043d0, &
               0.2036537393100009554d-1, &
               0.136023584095623632d-2, &
               0.6675388699209093d-4, &
               0.250430035707337d-5, &
               0.7406423741728d-7, &
               0.176974704314d-8, &
               0.3485775254d-10, &
               0.57544785d-12, &
               0.807481d-14, &
               0.9747d-16, &
               0.102d-17, &
               0.1d-19/
  data ak0ina(0)/  1.91172065445060453895d0/
  data ak0ina(1)/ -0.4183064565769581085d-1/
  data ak0ina(2)/  0.213352508068147486d-2/
  data ak0ina(3)/ -0.15859497284504181d-3/
  data ak0ina(4)/  0.1497624699858351d-4/
  data ak0ina(5)/ -0.167955955322241d-5/
  data ak0ina(6)/  0.21495472478804d-6/
  data ak0ina(7)/ -0.3058356654790d-7/
  data ak0ina(8)/  0.474946413343d-8/
  data ak0ina(9)/ -0.79424660432d-9/
  data ak0ina(10)/ 0.14156555325d-9/
  data ak0ina(11)/-0.2667825359d-10/
  data ak0ina(12)/ 0.528149717d-11/
  data ak0ina(13)/-0.109263199d-11/
  data ak0ina(14)/ 0.23518838d-12/
  data ak0ina(15)/-0.5247991d-13/
  data ak0ina(16)/ 0.1210191d-13/
  data ak0ina(17)/-0.287632d-14/
  data ak0ina(18)/ 0.70297d-15/
  data ak0ina(19)/-0.17631d-15/
  data ak0ina(20)/ 0.4530d-16/
  data ak0ina(21)/-0.1190d-16/
  data ak0ina(22)/ 0.319d-17/
  data ak0ina(23)/-0.87d-18/
  data ak0ina(24)/ 0.24d-18/
  data ak0ina(25)/-0.7d-19/
  data ak0ina(26)/ 0.2d-19/
  data ak0ina(27)/-0.1d-19/
!
!   Machine-dependent values (suitable for IEEE machines)
!
  data xlow,xhigh/4.47034836d-8,36.0436534d0/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BESSEL_K0_INT - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    bessel_k0_int = zero
  else if ( x == zero ) then
    bessel_k0_int = zero
  else if ( x < xlow ) then
    bessel_k0_int = x * ( const1 - log ( x ) )
  else if ( x <= six ) then
    t = ( ( x * x ) / eightn - half ) - half
    fval = ( const2 + log ( x ) ) * cheval ( nterm2, ak0in2, t )
    bessel_k0_int = x * ( cheval ( nterm1, ak0in1, t ) - fval )
  else if ( x < xhigh ) then
    fval = piby2
    t = ( twelve / x - half ) - half
    temp = exp ( -x ) * cheval ( nterm3, ak0ina, t )
    fval = fval - temp / ( sqrt ( x ) * rt2bpi )
    bessel_k0_int = fval
  else
    bessel_k0_int = piby2
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
!      K0_INT(x) = Integral ( 0 <= t <= x ) K0(t) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function bessel_y0_int ( xvalue )

!*****************************************************************************80
!
!! BESSEL_Y0_INT calculates the integral of the Bessel function Y0.
!
!  Discussion:
!
!    The function is defined by:
!
!      Y0_INT(x) = Integral ( 0 <= t <= x ) Y0(t) dt
!
!    The code uses Chebyshev expansions whose coefficients are
!    given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    23 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) BESSEL_Y0_INT, the value of the function.
!
  implicit none

  real ( kind = 8 ) bessel_y0_int
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: nine = 9.0D+00
  integer, parameter :: nterm1 = 22
  integer, parameter :: nterm2 = 22
  integer, parameter :: nterm3 = 17
  integer, parameter :: nterm4 = 15
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: sixten = 16.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) arj01(0:23),ary01(0:24),ary0a1(0:21), &
       ary0a2(0:18),five12,gal2m1,gamln2, &
       one28,pib41,pib411,pib412, &
       pib42,rt2bpi,t,temp,twobpi,xhigh, &
       xlow,xmpi4
  data one28,five12/ 128.0d0 , 512.0d0 /
  data rt2bpi/0.79788456080286535588d0/
  data pib411,pib412/ 201.0d0, 256.0d0/
  data pib42/0.24191339744830961566d-3/
  data twobpi/0.63661977236758134308d0/
  data gal2m1/-1.11593151565841244881d0/
  data gamln2/-0.11593151565841244881d0/
  data arj01(0)/  0.38179279321690173518d0/
  data arj01(1)/ -0.21275636350505321870d0/
  data arj01(2)/  0.16754213407215794187d0/
  data arj01(3)/ -0.12853209772196398954d0/
  data arj01(4)/  0.10114405455778847013d0/
  data arj01(5)/ -0.9100795343201568859d-1/
  data arj01(6)/  0.6401345264656873103d-1/
  data arj01(7)/ -0.3066963029926754312d-1/
  data arj01(8)/  0.1030836525325064201d-1/
  data arj01(9)/ -0.255670650399956918d-2/
  data arj01(10)/ 0.48832755805798304d-3/
  data arj01(11)/-0.7424935126036077d-4/
  data arj01(12)/ 0.922260563730861d-5/
  data arj01(13)/-0.95522828307083d-6/
  data arj01(14)/ 0.8388355845986d-7/
  data arj01(15)/-0.633184488858d-8/
  data arj01(16)/ 0.41560504221d-9/
  data arj01(17)/-0.2395529307d-10/
  data arj01(18)/ 0.122286885d-11/
  data arj01(19)/-0.5569711d-13/
  data arj01(20)/ 0.227820d-14/
  data arj01(21)/-0.8417d-16/
  data arj01(22)/ 0.282d-17/
  data arj01(23)/-0.9d-19/
  data ary01(0)/  0.54492696302724365490d0/
  data ary01(1)/ -0.14957323588684782157d0/
  data ary01(2)/  0.11085634486254842337d0/
  data ary01(3)/ -0.9495330018683777109d-1/
  data ary01(4)/  0.6820817786991456963d-1/
  data ary01(5)/ -0.10324653383368200408d0/
  data ary01(6)/  0.10625703287534425491d0/
  data ary01(7)/ -0.6258367679961681990d-1/
  data ary01(8)/  0.2385645760338293285d-1/
  data ary01(9)/ -0.644864913015404481d-2/
  data ary01(10)/ 0.131287082891002331d-2/
  data ary01(11)/-0.20988088174989640d-3/
  data ary01(12)/ 0.2716042484138347d-4/
  data ary01(13)/-0.291199114014694d-5/
  data ary01(14)/ 0.26344333093795d-6/
  data ary01(15)/-0.2041172069780d-7/
  data ary01(16)/ 0.137124781317d-8/
  data ary01(17)/-0.8070680792d-10/
  data ary01(18)/ 0.419883057d-11/
  data ary01(19)/-0.19459104d-12/
  data ary01(20)/ 0.808782d-14/
  data ary01(21)/-0.30329d-15/
  data ary01(22)/ 0.1032d-16/
  data ary01(23)/-0.32d-18/
  data ary01(24)/ 0.1d-19/
  data ary0a1(0)/  1.24030133037518970827d0/
  data ary0a1(1)/ -0.478125353632280693d-2/
  data ary0a1(2)/  0.6613148891706678d-4/
  data ary0a1(3)/ -0.186042740486349d-5/
  data ary0a1(4)/  0.8362735565080d-7/
  data ary0a1(5)/ -0.525857036731d-8/
  data ary0a1(6)/  0.42606363251d-9/
  data ary0a1(7)/ -0.4211761024d-10/
  data ary0a1(8)/  0.488946426d-11/
  data ary0a1(9)/ -0.64834929d-12/
  data ary0a1(10)/ 0.9617234d-13/
  data ary0a1(11)/-0.1570367d-13/
  data ary0a1(12)/ 0.278712d-14/
  data ary0a1(13)/-0.53222d-15/
  data ary0a1(14)/ 0.10844d-15/
  data ary0a1(15)/-0.2342d-16/
  data ary0a1(16)/ 0.533d-17/
  data ary0a1(17)/-0.127d-17/
  data ary0a1(18)/ 0.32d-18/
  data ary0a1(19)/-0.8d-19/
  data ary0a1(20)/ 0.2d-19/
  data ary0a1(21)/-0.1d-19/
  data ary0a2(0)/  1.99616096301341675339d0/
  data ary0a2(1)/ -0.190379819246668161d-2/
  data ary0a2(2)/  0.1539710927044226d-4/
  data ary0a2(3)/ -0.31145088328103d-6/
  data ary0a2(4)/  0.1110850971321d-7/
  data ary0a2(5)/ -0.58666787123d-9/
  data ary0a2(6)/  0.4139926949d-10/
  data ary0a2(7)/ -0.365398763d-11/
  data ary0a2(8)/  0.38557568d-12/
  data ary0a2(9)/ -0.4709800d-13/
  data ary0a2(10)/ 0.650220d-14/
  data ary0a2(11)/-0.99624d-15/
  data ary0a2(12)/ 0.16700d-15/
  data ary0a2(13)/-0.3028d-16/
  data ary0a2(14)/ 0.589d-17/
  data ary0a2(15)/-0.122d-17/
  data ary0a2(16)/ 0.27d-18/
  data ary0a2(17)/-0.6d-19/
  data ary0a2(18)/ 0.1d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow,xhigh/3.16101364d-8,9.007199256d15/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BESSEL_Y0_INT - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    bessel_y0_int = zero
  else if ( x == zero ) then
    bessel_y0_int = zero
  else if ( x < xlow ) then
    bessel_y0_int = ( log ( x ) + gal2m1 ) * twobpi * x
  else if ( x <= sixten ) then
    t = x * x / one28 - one
    temp = ( log ( x ) + gamln2 ) * cheval ( nterm1, arj01, t )
    temp = temp - cheval ( nterm2, ary01, t )
    bessel_y0_int = twobpi * x * temp
  else if ( x <= xhigh ) then
    t = five12 / ( x * x ) - one
    pib41 = pib411 / pib412
    xmpi4 = ( x - pib41 ) - pib42
    temp = sin ( xmpi4 ) * cheval ( nterm3, ary0a1, t ) / x
    temp = temp + cos ( xmpi4 ) * cheval ( nterm4, ary0a2, t )
    bessel_y0_int = - rt2bpi * temp / sqrt ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BESSEL_Y0_INT - Fatal error!'
    write ( *, '(a)' ) '  Argument too large.'
    bessel_y0_int = zero
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
!      Y0_INT(x) = Integral ( 0 <= t <= x ) Y0(t) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
!
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function cheval ( n, a, t )

!*****************************************************************************80
!
!! CHEVAL evaluates a Chebyshev series.
!
!  Discussion:
!
!    This function evaluates a Chebyshev series, using the
!    Clenshaw method with Reinsch modification, as analysed
!    in the paper by Oliver.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    J Oliver,
!    An error analysis of the modified Clenshaw method for
!    evaluating Chebyshev and Fourier series,
!    Journal of the IMA,
!    Volume 20, 1977, pages 379-391.
!
!  Parameters:
!
!    Input, integer N, the number of terms in the sequence.
!
!    Input, real ( kind = 8 ) A(0:N), the coefficients of the Chebyshev series.
!
!    Input, real ( kind = 8 ) T, the value at which the series is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CHEVAL, the value of the Chebyshev series at T.
!
  implicit none

  integer n

  real ( kind = 8 ) :: a(0:n)
  real ( kind = 8 ) :: cheval
  real ( kind = 8 ) :: d1
  real ( kind = 8 ) :: d2
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer i
  real ( kind = 8 ) :: t
  real ( kind = 8 ), parameter :: test = 0.6D+00
  real ( kind = 8 ) :: tt
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) :: u0
  real ( kind = 8 ) :: u1
  real ( kind = 8 ) :: u2
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  u1 = zero
!
!  T <= -0.6, Reinsch modification.
!
  if ( t <= -test ) then

    d1 = zero
    tt = ( t + half ) + half
    tt = tt + tt

    do i = n, 0, -1
      d2 = d1
      u2 = u1
      d1 = tt * u2 + a(i) - d2
      u1 = d1 - u2
    end do

    cheval = ( d1 - d2 ) / two
!
!  -0.6 < T < 0.6, Standard Clenshaw method.
!
  else if ( t < test ) then

    u0 = zero
    tt = t + t

    do i = n, 0, -1
      u2 = u1
      u1 = u0
      u0 = tt * u1 + a(i) - u2
    end do

    cheval = ( u0 - u2 ) / two
!
!  0.6 <= T, Reinsch modification.
!
  else

    d1 = zero
    tt = ( t - half ) - half
    tt = tt + tt

    do i = n, 0, -1
      d2 = d1
      u2 = u1
      d1 = tt * u2 + a(i) + d2
      u1 = d1 + u2
    end do

    cheval = ( d1 + d2 ) / two

  end if

  return
end
function clausen ( xvalue )

!*****************************************************************************80
!
!! CLAUSEN calculates Clausen's integral.
!
!  Discussion:
!
!    The function is defined by:
!
!      CLAUSEN(x) = Integral ( 0 <= t <= x ) -ln ( 2 * sin ( t / 2 ) ) dt
!
!    The code uses Chebyshev expansions with the coefficients
!    given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) CLAUSEN, the value of the function.
!
  implicit none

  real ( kind = 8 ) aclaus(0:15)
  real ( kind = 8 ) cheval
  real ( kind = 8 ) clausen
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer indx
  integer, parameter :: nterms = 13
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00
  real ( kind = 8 ), parameter :: pisq = 9.8696044010893586188D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) twopi,twopia,twopib,xhigh,xsmall
  data twopi/6.2831853071795864769d0/
  data twopia,twopib/6.28125d0, 0.19353071795864769253d-2/
  data aclaus/2.14269436376668844709d0, &
              0.7233242812212579245d-1, &
              0.101642475021151164d-2, &
              0.3245250328531645d-4, &
              0.133315187571472d-5, &
              0.6213240591653d-7, &
              0.313004135337d-8, &
              0.16635723056d-9, &
              0.919659293d-11, &
              0.52400462d-12, &
              0.3058040d-13, &
              0.181969d-14, &
              0.11004d-15, &
              0.675d-17, &
              0.42d-18, &
              0.3d-19/
!
!  Set machine-dependent constants (suitable for IEEE machines)
!
  data xsmall,xhigh/2.3406689d-8,4.5036d15/

  x = xvalue

  if ( xhigh < abs ( x ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLAUSEN - Warning!'
    write ( *, '(a)' ) '  Argument magnitude too large for accurate computation.'
    clausen = zero
    return
  end if

  indx = 1
  if ( x < zero ) then
    x = -x
    indx = -1
  end if
!
!  Argument reduced using simulated extra precision
!
  if ( twopi < x ) then
    t = aint ( x / twopi )
    x =  ( x - t * twopia ) - t * twopib
  end if

  if ( pi < x ) then
    x = ( twopia - x ) + twopib
    indx = -indx
  end if

  if ( x == zero ) then
    clausen = zero
  else if ( x < xsmall ) then
    clausen = x * ( one - log ( x ) )
  else
    t = ( x * x ) / pisq - half
    t = t + t
    if ( one < t ) then
      t = one
    end if

    clausen = x * cheval ( nterms, aclaus, t ) - x * log ( x )

  end if

  if ( indx < 0 ) then
    clausen = -clausen
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
!      CLAUSEN(x) = Integral ( 0 <= t <= x ) -ln ( 2 * sin ( t / 2 ) ) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function debye1 ( xvalue )

!*****************************************************************************80
!
!! DEBYE1 calculates the Debye function of order 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE1(x) = 1 / x * Integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) DEBYE1, the value of the function.
!
  implicit none

  real ( kind = 8 ) adeb1(0:18)
  real ( kind = 8 ) cheval
  real ( kind = 8 ) debye1
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer i
  integer nexp
  integer, parameter :: nterms = 15
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: quart = 0.25D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) debinf,expmx, &
       nine,rk,sum1,t,thirt6,xk,xlim,xlow, &
       xupper
  data nine,thirt6 /9.0d0, 36.0d0 /
  data debinf/0.60792710185402662866d0/
  data adeb1/2.40065971903814101941d0, &
             0.19372130421893600885d0, &
            -0.623291245548957703d-2, &
             0.35111747702064800d-3, &
            -0.2282224667012310d-4, &
             0.158054678750300d-5, &
            -0.11353781970719d-6, &
             0.835833611875d-8, &
            -0.62644247872d-9, &
             0.4760334890d-10, &
            -0.365741540d-11, &
             0.28354310d-12, &
            -0.2214729d-13, &
             0.174092d-14, &
            -0.13759d-15, &
             0.1093d-16, &
            -0.87d-18, &
             0.7d-19, &
            -0.1d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow,xupper,xlim/0.298023d-7,35.35051d0,708.39642d0/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DEBYE1 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    debye1 = zero
  else if ( x < xlow ) then
    debye1 = ( ( x - nine ) * x + thirt6 ) / thirt6
  else if ( x <= four ) then
    t = ( ( x * x / eight ) - half ) - half
    debye1 = cheval ( nterms, adeb1, t ) - quart * x
  else

    debye1 = one / ( x * debinf )
    if ( x < xlim ) then
      expmx = exp ( -x )
      if ( xupper < x ) then
        debye1 = debye1 - expmx * ( one + one / x )
      else
        sum1 = zero
        rk = aint ( xlim / x )
        nexp = int ( rk )
        xk = rk * x
        do i = nexp, 1, -1
          t = ( one + one / xk ) / rk
          sum1 = sum1 * expmx + t
          rk = rk - one
          xk = xk - x
        end do
        debye1 = debye1 - sum1 * expmx
      end if
    end if
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
!      DEBYE1(x) = 1 / x * Integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
!
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function debye2 ( xvalue )

!*****************************************************************************80
!
!! DEBYE2 calculates the Debye function of order 2.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE2(x) = 2 / x^2 * Integral ( 0 <= t <= x ) t^2 / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    24 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) DEBYE2, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ) debye2
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer i
  integer nexp
  integer, parameter :: nterms = 17
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) adeb2(0:18),debinf,expmx, &
       rk,sum1,t,twent4,xk,xlim1, &
       xlim2,xlow,xupper
  data twent4/24.0d0/
  data debinf/4.80822761263837714160d0/
  data adeb2/2.59438102325707702826d0, &
             0.28633572045307198337d0, &
            -0.1020626561580467129d-1, &
             0.60491097753468435d-3, &
            -0.4052576589502104d-4, &
             0.286338263288107d-5, &
            -0.20863943030651d-6, &
             0.1552378758264d-7, &
            -0.117312800866d-8, &
             0.8973585888d-10, &
            -0.693176137d-11, &
             0.53980568d-12, &
            -0.4232405d-13, &
             0.333778d-14, &
            -0.26455d-15, &
             0.2106d-16, &
            -0.168d-17, &
             0.13d-18, &
            -0.1d-19/
!
!   Machine-dependent constants
!
  data xlow,xupper/0.298023d-7,35.35051d0/
  data xlim1,xlim2/708.39642d0,2.1572317d154/
!
  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DEBYE2 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    debye2 = zero
  else if ( x < xlow ) then
    debye2 = ( ( x - eight ) * x + twent4 ) / twent4
  else if ( x <= four ) then
    t = ( ( x * x / eight ) - half ) - half
    debye2 = cheval ( nterms, adeb2, t ) - x / three
  else if ( x <= xupper ) then

    expmx = exp ( -x )
    sum1 = zero
    rk = aint ( xlim1 / x )
    nexp = int ( rk )
    xk = rk * x
    do i = nexp, 1, -1
      t =  ( one + two / xk + two / ( xk * xk ) ) / rk
      sum1 = sum1 * expmx + t
      rk = rk - one
      xk = xk - x
    end do
    debye2 = debinf / ( x * x ) - two * sum1 * expmx

  else if ( x < xlim1 ) then

    expmx = exp ( -x )
    sum1 = ( ( x + two ) * x + two ) / ( x * x )
    debye2 = debinf / ( x * x ) - two * sum1 * expmx

  else if ( x <= xlim2 ) then
    debye2 = debinf / ( x * x )
  else
    debye2 = zero
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
!      DEBYE2(x) = 2 / x^2 * Integral ( 0 <= t <= x ) t^2 / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
!
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function debye3 ( xvalue )

!*****************************************************************************80
!
!! DEBYE3 calculates the Debye function of order 3.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE3(x) = 3 / x^3 * Integral ( 0 <= t <= x ) t^3 / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) DEBYE3, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ) debye3
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer i
  integer nexp
  integer, parameter :: nterms = 16
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) adeb3(0:18),debinf,expmx, &
       pt375,rk,sevp5,sum1,t,twenty, &
       xk,xki,xlim1,xlim2,xlow,xupper
  data pt375/0.375d0/
  data sevp5,twenty/7.5d0 , 20.0d0/
  data debinf/0.51329911273421675946d-1/
  data adeb3/2.70773706832744094526d0, &
             0.34006813521109175100d0, &
            -0.1294515018444086863d-1, &
             0.79637553801738164d-3, &
            -0.5463600095908238d-4, &
             0.392430195988049d-5, &
            -0.28940328235386d-6, &
             0.2173176139625d-7, &
            -0.165420999498d-8, &
             0.12727961892d-9, &
            -0.987963459d-11, &
             0.77250740d-12, &
            -0.6077972d-13, &
             0.480759d-14, &
            -0.38204d-15, &
             0.3048d-16, &
            -0.244d-17, &
             0.20d-18, &
            -0.2d-19/
!
!   Machine-dependent constants
!
  data xlow,xupper/0.298023d-7,35.35051d0/
  data xlim1,xlim2/708.39642d0,0.9487163d103/
!
  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DEBYE3 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    debye3 = zero
    return
  end if

  if ( x < xlow ) then
    debye3 = ( ( x - sevp5 ) * x + twenty ) / twenty
  else if ( x <= 4 ) then
    t = ( ( x * x / eight ) - half ) - half
    debye3 = cheval ( nterms, adeb3, t ) - pt375 * x
  else
!
!   Code for x > 4.0
!
     if ( xlim2 < x ) then
        debye3 = zero
     else
        debye3 = one / ( debinf * x * x * x )
        if ( x < xlim1 ) then
           expmx = exp ( -x )
           if ( xupper < x ) then
              sum1 = ((( x + three ) * x + six ) * x + six ) / ( x * x * x )
           else
              sum1 = zero
              rk = aint ( xlim1 / x )
              nexp = int ( rk )
              xk = rk * x
              do i = nexp, 1, -1
                 xki = one / xk
                 t =  ((( six * xki + six ) * xki + three ) * xki + one ) / rk
                 sum1 = sum1 * expmx + t
                 rk = rk - one
                 xk = xk - x
              end do
           end if
           debye3 = debye3 - three * sum1 * expmx
        end if
     end if
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
!      DEBYE3(x) = 3 / x^3 * Integral ( 0 <= t <= x ) t^3 / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
!
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function debye4 ( xvalue )

!*****************************************************************************80
!
!! DEBYE4 calculates the Debye function of order 4.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE4(x) = 4 / x^4 * Integral ( 0 <= t <= x ) t^4 / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) DEBYE4, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ) debye4
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four =  4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer i
  integer nexp
  integer, parameter :: nterms = 16
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) adeb4(0:18),debinf,eightn,expmx, &
       five,forty5,rk,sum1,t,twelve,twent4, &
       twopt5,xk,xki,xlim1,xlim2,xlow,xupper
  data twopt5,five/2.5d0, 5.0d0/
  data twelve,eightn/ 12.0d0 , 18.0d0/
  data twent4,forty5 /24.0d0 , 45.0d0 /
  data debinf/99.54506449376351292781d0/
  data adeb4/2.78186941502052346008d0, &
             0.37497678352689286364d0, &
            -0.1494090739903158326d-1, &
             0.94567981143704274d-3, &
            -0.6613291613893255d-4, &
             0.481563298214449d-5, &
            -0.35880839587593d-6, &
             0.2716011874160d-7, &
            -0.208070991223d-8, &
             0.16093838692d-9, &
            -0.1254709791d-10, &
             0.98472647d-12, &
            -0.7772369d-13, &
             0.616483d-14, &
            -0.49107d-15, &
             0.3927d-16, &
            -0.315d-17, &
             0.25d-18, &
            -0.2d-19/
!
!   Machine-dependent constants
!
  data xlow,xupper/0.298023d-7,35.35051d0/
  data xlim1,xlim2/708.39642d0,2.5826924d77/

  x = xvalue
!
!   Check XVALUE >= 0.0
!
  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DEBYE4 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    debye4 = zero
    return
  end if

  if ( x < xlow ) then
    debye4 = ( ( twopt5 * x - eightn ) * x + forty5 ) / forty5
  else if ( x <= four ) then
    t = ( ( x * x / eight ) - half ) - half
    debye4 = cheval ( nterms, adeb4, t ) - ( x + x ) / five
  else
!
!   Code for x > 4.0
!
     if ( xlim2 < x ) then
        debye4 = zero
     else
        t = x * x
        debye4 = ( debinf / t ) / t
        if ( x < xlim1 ) then
           expmx = exp ( -x )
           if ( xupper < x ) then
              sum1 = ( ( ( ( x + four ) * x + twelve ) * x + &
                    twent4 ) * x + twent4 ) / ( x * x * x * x )
           else
              sum1 = zero
              rk = aint ( xlim1 / x )
              nexp = int ( rk )
              xk = rk * x
              do i = nexp, 1, -1
                 xki = one / xk
                 t =  ( ( ( ( twent4 * xki + twent4 ) * xki + &
                      twelve ) * xki + four ) * xki + one ) / rk
                 sum1 = sum1 * expmx + t
                 rk = rk - one
                 xk = xk - x
              end do
           end if
           debye4 = debye4 - four * sum1 * expmx
        end if
     end if
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
!      DEBYE4(x) = 4 / x^4 * Integral ( 0 <= t <= x ) t^4 / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
!
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function exp3_int ( xvalue )

!*****************************************************************************80
!
!! EXP3_INT calculates the integral of exp(-t^3).
!
!  Discussion:
!
!    The function is defined by:
!
!      EXP3_INT(x) = Integral ( 0 <= t <= x ) exp ( -t^3 ) dt
!
!    The code uses Chebyshev expansions, whose coefficients are
!    given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) EXP3_INT, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ) exp3_int
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer, parameter :: nterm1 = 22
  integer, parameter :: nterm2 = 20
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) aexp3(0:24),aexp3a(0:24), &
         funinf,sixten,t, &
         xlow,xupper
  data sixten /16.0d0 /
  data funinf/0.89297951156924921122d0/
  data aexp3(0)/  1.26919841422112601434d0/
  data aexp3(1)/ -0.24884644638414098226d0/
  data aexp3(2)/  0.8052622071723104125d-1/
  data aexp3(3)/ -0.2577273325196832934d-1/
  data aexp3(4)/  0.759987887307377429d-2/
  data aexp3(5)/ -0.203069558194040510d-2/
  data aexp3(6)/  0.49083458669932917d-3/
  data aexp3(7)/ -0.10768223914202077d-3/
  data aexp3(8)/  0.2155172626428984d-4/
  data aexp3(9)/ -0.395670513738429d-5/
  data aexp3(10)/ 0.66992409338956d-6/
  data aexp3(11)/-0.10513218080703d-6/
  data aexp3(12)/ 0.1536258019825d-7/
  data aexp3(13)/-0.209909603636d-8/
  data aexp3(14)/ 0.26921095381d-9/
  data aexp3(15)/-0.3251952422d-10/
  data aexp3(16)/ 0.371148157d-11/
  data aexp3(17)/-0.40136518d-12/
  data aexp3(18)/ 0.4123346d-13/
  data aexp3(19)/-0.403375d-14/
  data aexp3(20)/ 0.37658d-15/
  data aexp3(21)/-0.3362d-16/
  data aexp3(22)/ 0.288d-17/
  data aexp3(23)/-0.24d-18/
  data aexp3(24)/ 0.2d-19/
  data aexp3a(0)/  1.92704649550682737293d0/
  data aexp3a(1)/ -0.3492935652048138054d-1/
  data aexp3a(2)/  0.145033837189830093d-2/
  data aexp3a(3)/ -0.8925336718327903d-4/
  data aexp3a(4)/  0.705423921911838d-5/
  data aexp3a(5)/ -0.66717274547611d-6/
  data aexp3a(6)/  0.7242675899824d-7/
  data aexp3a(7)/ -0.878258256056d-8/
  data aexp3a(8)/  0.116722344278d-8/
  data aexp3a(9)/ -0.16766312812d-9/
  data aexp3a(10)/ 0.2575501577d-10/
  data aexp3a(11)/-0.419578881d-11/
  data aexp3a(12)/ 0.72010412d-12/
  data aexp3a(13)/-0.12949055d-12/
  data aexp3a(14)/ 0.2428703d-13/
  data aexp3a(15)/-0.473311d-14/
  data aexp3a(16)/ 0.95531d-15/
  data aexp3a(17)/-0.19914d-15/
  data aexp3a(18)/ 0.4277d-16/
  data aexp3a(19)/-0.944d-17/
  data aexp3a(20)/ 0.214d-17/
  data aexp3a(21)/-0.50d-18/
  data aexp3a(22)/ 0.12d-18/
  data aexp3a(23)/-0.3d-19/
  data aexp3a(24)/ 0.1d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow,xupper/0.762939d-5,3.3243018d0/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXP3_INT - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    exp3_int = zero
  else if ( x < xlow ) then
    exp3_int = x
  else if ( x <= two ) then
    t = ( ( x * x * x / four ) - half ) - half
    exp3_int = x * cheval ( nterm1, aexp3, t )
  else if ( x <= xupper ) then
    t = ( ( sixten / ( x * x * x ) ) - half ) - half
    t = cheval ( nterm2, aexp3a, t )
    t = t * exp ( -x * x * x ) / ( three * x * x )
    exp3_int = funinf - t
  else
    exp3_int = funinf
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
!      EXP3_INT(x) = Integral ( 0 <= t <= x ) exp ( -t^3 ) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function goodwin ( xvalue )

!*****************************************************************************80
!
!! GOODWIN calculates the integral of exp(-t^2/(t+x)).
!
!  Discussion:
!
!    The function is defined by:
!
!      GOODWIN(x) = Integral ( 0 <= t < infinity ) exp ( -t^2 ) / ( t + x ) dt
!
!    The code uses Chebyshev expansions whose coefficients are
!    given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) GOODWIN, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ) goodwin
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer, parameter :: nterm1 = 26
  integer, parameter :: nterm2 = 20
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) agost(0:28),agosta(0:23), &
       fval,gamby2,rtpib2, &
       t,xhigh,xlow
  data gamby2/0.28860783245076643030d0/
  data rtpib2/0.88622692545275801365d0/
  data agost(0)/  0.63106560560398446247d0/
  data agost(1)/  0.25051737793216708827d0/
  data agost(2)/ -0.28466205979018940757d0/
  data agost(3)/  0.8761587523948623552d-1/
  data agost(4)/  0.682602267221252724d-2/
  data agost(5)/ -0.1081129544192254677d-1/
  data agost(6)/  0.169101244117152176d-2/
  data agost(7)/  0.50272984622615186d-3/
  data agost(8)/ -0.18576687204100084d-3/
  data agost(9)/ -0.428703674168474d-5/
  data agost(10)/ 0.1009598903202905d-4/
  data agost(11)/-0.86529913517382d-6/
  data agost(12)/-0.34983874320734d-6/
  data agost(13)/ 0.6483278683494d-7/
  data agost(14)/ 0.757592498583d-8/
  data agost(15)/-0.277935424362d-8/
  data agost(16)/-0.4830235135d-10/
  data agost(17)/ 0.8663221283d-10/
  data agost(18)/-0.394339687d-11/
  data agost(19)/-0.209529625d-11/
  data agost(20)/ 0.21501759d-12/
  data agost(21)/ 0.3959015d-13/
  data agost(22)/-0.692279d-14/
  data agost(23)/-0.54829d-15/
  data agost(24)/ 0.17108d-15/
  data agost(25)/ 0.376d-17/
  data agost(26)/-0.349d-17/
  data agost(27)/ 0.7d-19/
  data agost(28)/ 0.6d-19/
  data agosta(0)/  1.81775467984718758767d0/
  data agosta(1)/ -0.9921146570744097467d-1/
  data agosta(2)/ -0.894058645254819243d-2/
  data agosta(3)/ -0.94955331277726785d-3/
  data agosta(4)/ -0.10971379966759665d-3/
  data agosta(5)/ -0.1346694539578590d-4/
  data agosta(6)/ -0.172749274308265d-5/
  data agosta(7)/ -0.22931380199498d-6/
  data agosta(8)/ -0.3127844178918d-7/
  data agosta(9)/ -0.436197973671d-8/
  data agosta(10)/-0.61958464743d-9/
  data agosta(11)/-0.8937991276d-10/
  data agosta(12)/-0.1306511094d-10/
  data agosta(13)/-0.193166876d-11/
  data agosta(14)/-0.28844270d-12/
  data agosta(15)/-0.4344796d-13/
  data agosta(16)/-0.659518d-14/
  data agosta(17)/-0.100801d-14/
  data agosta(18)/-0.15502d-15/
  data agosta(19)/-0.2397d-16/
  data agosta(20)/-0.373d-17/
  data agosta(21)/-0.58d-18/
  data agosta(22)/-0.9d-19/
  data agosta(23)/-0.1d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow,xhigh/1.11022303d-16,1.80144d16/

  x = xvalue

  if ( x <= zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GOODWIN - Fatal error!'
    write ( *, '(a)' ) '  Argument X <= 0.'
    goodwin = zero
  else if ( x < xlow ) then
    goodwin = - gamby2 - log ( x )
  else if ( x <= two ) then
    t = ( x - half ) - half
    goodwin = cheval ( nterm1, agost, t ) - exp ( -x * x ) * log ( x )
  else if ( x <= xhigh ) then
    fval = rtpib2 / x
    t = ( six - x ) / ( two + x )
    goodwin = fval * cheval ( nterm2, agosta, t )
  else
    goodwin = rtpib2 / x
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
!      GOODWIN(x) = Integral ( 0 <= t < infinity ) exp ( -t^2 ) / ( t + x ) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function i0ml0 ( xvalue )

!*****************************************************************************80
!
!! I0ML0 calculates the difference between the Bessel I0 and Struve L0 functions.
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
!    The code uses Chebyshev expansions with the coefficients
!    given to an accuracy of 20D.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) I0ML0, the value of the function.
!
  implicit none

  real ( kind = 8 ) ai0l0(0:23)
  real ( kind = 8 ) ai0l0a(0:23)
  real ( kind = 8 ) cheval
  real ( kind = 8 ) i0ml0
  integer, parameter :: nterm1 = 21
  integer, parameter :: nterm2 = 21
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) atehun, &
       forty,sixten,t,twobpi,two88,xhigh, &
       xlow,xsq
  data sixten/ 16.0d0 /
  data forty / 40.0d0 /
  data two88,atehun/ 288.0d0, 800.0d0 /
  data twobpi/0.63661977236758134308d0/
  data ai0l0(0)/  0.52468736791485599138d0/
  data ai0l0(1)/ -0.35612460699650586196d0/
  data ai0l0(2)/  0.20487202864009927687d0/
  data ai0l0(3)/ -0.10418640520402693629d0/
  data ai0l0(4)/  0.4634211095548429228d-1/
  data ai0l0(5)/ -0.1790587192403498630d-1/
  data ai0l0(6)/  0.597968695481143177d-2/
  data ai0l0(7)/ -0.171777547693565429d-2/
  data ai0l0(8)/  0.42204654469171422d-3/
  data ai0l0(9)/ -0.8796178522094125d-4/
  data ai0l0(10)/ 0.1535434234869223d-4/
  data ai0l0(11)/-0.219780769584743d-5/
  data ai0l0(12)/ 0.24820683936666d-6/
  data ai0l0(13)/-0.2032706035607d-7/
  data ai0l0(14)/ 0.90984198421d-9/
  data ai0l0(15)/ 0.2561793929d-10/
  data ai0l0(16)/-0.710609790d-11/
  data ai0l0(17)/ 0.32716960d-12/
  data ai0l0(18)/ 0.2300215d-13/
  data ai0l0(19)/-0.292109d-14/
  data ai0l0(20)/-0.3566d-16/
  data ai0l0(21)/ 0.1832d-16/
  data ai0l0(22)/-0.10d-18/
  data ai0l0(23)/-0.11d-18/
  data ai0l0a(0)/ 2.00326510241160643125d0/
  data ai0l0a(1)/ 0.195206851576492081d-2/
  data ai0l0a(2)/ 0.38239523569908328d-3/
  data ai0l0a(3)/ 0.7534280817054436d-4/
  data ai0l0a(4)/ 0.1495957655897078d-4/
  data ai0l0a(5)/ 0.299940531210557d-5/
  data ai0l0a(6)/ 0.60769604822459d-6/
  data ai0l0a(7)/ 0.12399495544506d-6/
  data ai0l0a(8)/ 0.2523262552649d-7/
  data ai0l0a(9)/ 0.504634857332d-8/
  data ai0l0a(10)/0.97913236230d-9/
  data ai0l0a(11)/0.18389115241d-9/
  data ai0l0a(12)/0.3376309278d-10/
  data ai0l0a(13)/0.611179703d-11/
  data ai0l0a(14)/0.108472972d-11/
  data ai0l0a(15)/0.18861271d-12/
  data ai0l0a(16)/0.3280345d-13/
  data ai0l0a(17)/0.565647d-14/
  data ai0l0a(18)/0.93300d-15/
  data ai0l0a(19)/0.15881d-15/
  data ai0l0a(20)/0.2791d-16/
  data ai0l0a(21)/0.389d-17/
  data ai0l0a(22)/0.70d-18/
  data ai0l0a(23)/0.16d-18/
!
!   MACHINE-DEPENDENT CONSTANTS (suitable for IEEE-arithmetic machines)
!
  data xlow,xhigh/1.11022303d-16,1.8981253d9/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I0ML0 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    i0ml0 = zero
  else if ( x < xlow ) then
    i0ml0 = one
  else if ( x <= sixten ) then
    t = ( six * x - forty ) / ( x + forty )
    i0ml0 = cheval ( nterm1, ai0l0, t )
  else if ( x <= xhigh ) then
    xsq = x * x
    t = ( atehun - xsq ) / ( two88 + xsq )
    i0ml0 = cheval ( nterm2, ai0l0a, t ) * twobpi / x
  else
    i0ml0 = twobpi / x
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function i1ml1 ( xvalue )

!*****************************************************************************80
!
!! I1ML1 calculates the difference between the Bessel I1 and Struve L1 functions.
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
!    The code uses Chebyshev expansions with the coefficients
!    given to an accuracy of 20D.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!    0 <= XVALUE is required.
!
!    Output, real ( kind = 8 ) I1ML1, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ) i1ml1
  integer, parameter :: nterm1 = 20
  integer, parameter :: nterm2 = 22
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) ai1l1(0:23),ai1l1a(0:25),atehun, &
       forty,sixten,t,twobpi,two88, &
       xhigh,xlow,xsq
  data sixten,forty/ 16.0d0 , 40.0d0 /
  data two88,atehun/ 288.0d0 , 800.0d0 /
  data twobpi/0.63661977236758134308d0/
  data ai1l1(0)/  0.67536369062350576137d0/
  data ai1l1(1)/ -0.38134971097266559040d0/
  data ai1l1(2)/  0.17452170775133943559d0/
  data ai1l1(3)/ -0.7062105887235025061d-1/
  data ai1l1(4)/  0.2517341413558803702d-1/
  data ai1l1(5)/ -0.787098561606423321d-2/
  data ai1l1(6)/  0.214814368651922006d-2/
  data ai1l1(7)/ -0.50862199717906236d-3/
  data ai1l1(8)/  0.10362608280442330d-3/
  data ai1l1(9)/ -0.1795447212057247d-4/
  data ai1l1(10)/ 0.259788274515414d-5/
  data ai1l1(11)/-0.30442406324667d-6/
  data ai1l1(12)/ 0.2720239894766d-7/
  data ai1l1(13)/-0.158126144190d-8/
  data ai1l1(14)/ 0.1816209172d-10/
  data ai1l1(15)/ 0.647967659d-11/
  data ai1l1(16)/-0.54113290d-12/
  data ai1l1(17)/-0.308311d-14/
  data ai1l1(18)/ 0.305638d-14/
  data ai1l1(19)/-0.9717d-16/
  data ai1l1(20)/-0.1422d-16/
  data ai1l1(21)/ 0.84d-18/
  data ai1l1(22)/ 0.7d-19/
  data ai1l1(23)/-0.1d-19/
  data ai1l1a(0)/  1.99679361896789136501d0/
  data ai1l1a(1)/ -0.190663261409686132d-2/
  data ai1l1a(2)/ -0.36094622410174481d-3/
  data ai1l1a(3)/ -0.6841847304599820d-4/
  data ai1l1a(4)/ -0.1299008228509426d-4/
  data ai1l1a(5)/ -0.247152188705765d-5/
  data ai1l1a(6)/ -0.47147839691972d-6/
  data ai1l1a(7)/ -0.9020819982592d-7/
  data ai1l1a(8)/ -0.1730458637504d-7/
  data ai1l1a(9)/ -0.332323670159d-8/
  data ai1l1a(10)/-0.63736421735d-9/
  data ai1l1a(11)/-0.12180239756d-9/
  data ai1l1a(12)/-0.2317346832d-10/
  data ai1l1a(13)/-0.439068833d-11/
  data ai1l1a(14)/-0.82847110d-12/
  data ai1l1a(15)/-0.15562249d-12/
  data ai1l1a(16)/-0.2913112d-13/
  data ai1l1a(17)/-0.543965d-14/
  data ai1l1a(18)/-0.101177d-14/
  data ai1l1a(19)/-0.18767d-15/
  data ai1l1a(20)/-0.3484d-16/
  data ai1l1a(21)/-0.643d-17/
  data ai1l1a(22)/-0.118d-17/
  data ai1l1a(23)/-0.22d-18/
  data ai1l1a(24)/-0.4d-19/
  data ai1l1a(25)/-0.1d-19/
!
!   MACHINE-DEPENDENT CONSTANTS (suitable for IEEE machines)
!
  data xlow,xhigh/2.22044605d-16,1.8981253d9/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1ML1 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    i1ml1 = zero
  else if ( x < xlow ) then
    i1ml1 = x / two
  else if ( x <= sixten ) then
    t = ( six * x - forty ) / ( x + forty )
    i1ml1 = cheval ( nterm1, ai1l1, t ) * x / two
  else if ( x <= xhigh ) then
    xsq = x * x
    t = ( atehun - xsq ) / ( two88 + xsq )
    i1ml1 = cheval ( nterm2, ai1l1a, t ) * twobpi
  else
    i1ml1 = twobpi
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function lobachevsky ( xvalue )

!*****************************************************************************80
!
!! LOBACHEVSKY calculates the Lobachevsky function.
!
!  Discussion:
!
!    The function is defined by:
!
!      LOBACHEVSKY(x) = Integral ( 0 <= t <= x ) -ln ( abs ( cos ( t ) ) dt
!
!    The code uses Chebyshev expansions whose coefficients are given
!    to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) LOBACHEVSKY, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer indpi2
  integer indsgn
  integer npi
  integer, parameter :: nterm1 = 13
  integer, parameter :: nterm2 = 9
  real ( kind = 8 ) lobachevsky
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) arlob1(0:15),arlob2(0:10), &
       fval,fval1,lbpb21,lbpb22,lobpia,lobpib, &
       lobpi1,lobpi2,pi,piby2,piby21,piby22,piby4,pi1, &
       pi11,pi12,pi2,t,tcon,xcub,xhigh,xlow1, &
       xlow2,xlow3,xr
  data lobpia,lobpib/ 1115.0d0 , 512.0d0 /
  data lobpi2/-1.48284696397869499311d-4/
  data lbpb22/-7.41423481989347496556d-5/
  data pi11,pi12/ 201.0d0 , 64.0d0 /
  data pi2/9.67653589793238462643d-4/
  data piby22/4.83826794896619231322d-4/
  data tcon/3.24227787655480868620d0/
  data arlob1/0.34464884953481300507d0, &
              0.584198357190277669d-2, &
              0.19175029694600330d-3, &
              0.787251606456769d-5, &
              0.36507477415804d-6, &
              0.1830287272680d-7, &
              0.96890333005d-9, &
              0.5339055444d-10, &
              0.303408025d-11, &
              0.17667875d-12, &
              0.1049393d-13, &
              0.63359d-15, &
              0.3878d-16, &
              0.240d-17, &
              0.15d-18, &
              0.1d-19/
  data arlob2/2.03459418036132851087d0, &
              0.1735185882027407681d-1, &
              0.5516280426090521d-4, &
              0.39781646276598d-6, &
              0.369018028918d-8, &
              0.3880409214d-10, &
              0.44069698d-12, &
              0.527674d-14, &
              0.6568d-16, &
              0.84d-18, &
              0.1d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow1,xlow2/5.11091385d-103,4.71216091d-8/
  data xlow3,xhigh/6.32202727d-8,4.5035996d15/

  x = abs ( xvalue )
  indsgn = 1
  if ( xvalue < zero ) then
    indsgn = -1
  end if

  if ( xhigh < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOBACHEVSKY - Fatal error!'
    write ( *, '(a)' ) '  Argument magnitude too large.'
    lobachevsky = zero
    return
  end if
!
!  Reduce argument to [0,pi]
!
  pi1 = pi11 / pi12
  pi = pi1 + pi2
  piby2 = pi / two
  piby21 = pi1 / two
  piby4 = piby2 / two
  npi = int ( x / pi )
  xr = ( x - npi * pi1 ) - npi * pi2
!
!  Reduce argument to [0,pi/2]
!
  indpi2 = 0
  if ( piby2 < xr ) then
    indpi2 = 1
    xr = ( pi1 - xr ) + pi2
  end if
!
!  Code for argument in [0,pi/4]
!
  if ( xr <= piby4 ) then
     if ( xr < xlow1 ) then
        fval = zero
     else
        xcub = xr * xr * xr
        if ( xr < xlow2 ) then
          fval = xcub / six
        else
          t = ( tcon * xr * xr - half ) - half
          fval = xcub * cheval ( nterm1, arlob1, t )
        end if
     end if
  else
!
!  Code for argument in [pi/4,pi/2]
!
     xr = ( piby21 - xr ) + piby22
     if ( xr == zero ) then
        fval1 = zero
     else
        if ( xr < xlow3 ) then
          fval1 = xr * ( one - log ( xr ) )
        else
          t = ( tcon * xr * xr - half ) - half
          fval1 = xr * ( cheval ( nterm2, arlob2, t ) - log ( xr ) )
        end if
     end if
     lbpb21 = lobpia / ( lobpib + lobpib )
     fval = ( lbpb21 - fval1 ) + lbpb22
  end if

  lobpi1 = lobpia / lobpib
!
!  Compute value for argument in [pi/2,pi]
!
  if ( indpi2 == 1 ) then
     fval = ( lobpi1 - fval ) + lobpi2
  end if

  if ( npi <= 0 ) then
    lobachevsky = fval
  else
    lobachevsky = ( fval + npi * lobpi2 ) + npi * lobpi1
  end if

  if ( indsgn == -1 ) then
    lobachevsky = -lobachevsky
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
!      LOBACHEVSKY(x) = Integral ( 0 <= t <= x ) -ln ( abs ( cos ( t ) ) dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function stromgen ( xvalue )

!*****************************************************************************80
!
!! STROMGEN calculates Stromgen's integral.
!
!  Discussion:
!
!    The function is defined by:
!
!      STROMGEN(X) = Integral ( 0 <= t <= X ) t^7 * exp(2*t) / (exp(t)-1)^3 dt
!
!    The code uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) STROMGEN, the value of the function.
!
  implicit none

  real ( kind = 8 ) astrom(0:26)
  real ( kind = 8 ) cheval
  real ( kind = 8 ) epngln
  real ( kind = 8 ) epsln
  real ( kind = 8 ) f15bp4
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer k1
  integer k2
  integer, parameter :: nterms = 23
  integer numexp
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) one5ln
  real ( kind = 8 ) pi4b3
  real ( kind = 8 ) rk
  real ( kind = 8 ) stromgen
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) seven,sumexp,sum2,t,valinf,xhigh, &
       xk,xk1,xlow0,xlow1
  data seven/ 7.0d0 /
  data one5ln/ 0.4055d0 /
  data f15bp4/0.38497433455066256959d-1 /
  data pi4b3/1.29878788045336582982d2 /
  data valinf/196.51956920868988261257d0/
  data astrom(0)/  0.56556120872539155290d0/
  data astrom(1)/  0.4555731969101785525d-1/
  data astrom(2)/ -0.4039535875936869170d-1/
  data astrom(3)/ -0.133390572021486815d-2/
  data astrom(4)/  0.185862506250538030d-2/
  data astrom(5)/ -0.4685555868053659d-4/
  data astrom(6)/ -0.6343475643422949d-4/
  data astrom(7)/  0.572548708143200d-5/
  data astrom(8)/  0.159352812216822d-5/
  data astrom(9)/ -0.28884328431036d-6/
  data astrom(10)/-0.2446633604801d-7/
  data astrom(11)/ 0.1007250382374d-7/
  data astrom(12)/-0.12482986104d-9/
  data astrom(13)/-0.26300625283d-9/
  data astrom(14)/ 0.2490407578d-10/
  data astrom(15)/ 0.485454902d-11/
  data astrom(16)/-0.105378913d-11/
  data astrom(17)/-0.3604417d-13/
  data astrom(18)/ 0.2992078d-13/
  data astrom(19)/-0.163971d-14/
  data astrom(20)/-0.61061d-15/
  data astrom(21)/ 0.9335d-16/
  data astrom(22)/ 0.709d-17/
  data astrom(23)/-0.291d-17/
  data astrom(24)/ 0.8d-19/
  data astrom(25)/ 0.6d-19/
  data astrom(26)/-0.1d-19/
!
!  Machine-dependent constants
!
  data xlow0,xlow1/7.80293d-62,2.22045d-16/
  data epsln,epngln/-36.0436534d0,-36.7368006d0/
  data xhigh/3.1525197d16/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STROMGEN - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    stromgen = zero
    return
  end if

  if ( x < xlow0 ) then
    stromgen = zero
  else if ( x < xlow1 ) then
    stromgen = x**5 / pi4b3
  else if ( x <= four ) then
    t = ( ( x / two ) - half ) - half
    stromgen = x**5 * cheval ( nterms, astrom, t ) * f15bp4
  else
!
!  Code for x > 4.0
!
    if ( xhigh < x ) then
      sumexp = one
    else
      numexp = int ( epsln / ( one5ln - x ) ) + 1
      if ( 1 < numexp ) then
        t = exp ( -x )
      else
        t = one
      end if
      rk = zero
      do k1 = 1, numexp
        rk = rk + one
      end do
      sumexp = zero
      do k1 = 1, numexp
        sum2 = one
        xk = one / ( rk * x )
        xk1 = one
        do k2 = 1, 7
          sum2 = sum2 * xk1 * xk + one
          xk1 = xk1 + one
        end do
        sum2 = sum2 * ( rk + one ) / two
        sumexp = sumexp * t + sum2
        rk = rk - one
      end do

    end if

    t = seven * log ( x ) - x + log ( sumexp )

    if ( t < epngln ) then
      stromgen = valinf
    else
      stromgen = valinf - exp ( t ) * f15bp4
    end if

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
!      STROMGEN(X) = Integral ( 0 <= t <= X ) t^7 * exp(2*t) / (exp(t)-1)^3 dt
!
!    The data was reported by McLeod.
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
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function struve_h0 ( xvalue )

!*****************************************************************************80
!
!! STRUVE_H0 calculates the Struve function of order 0.
!
!  Discussion:
!
!    The function is defined by:
!
!      HO(x) = (2/pi) Integral ( 0 <= t <= pi/2 ) sin ( x * cos ( t ) ) dt
!
!    H0 also satisfies the second-order equation
!
!      x*D(Df) + Df + x * f = 2 * x / pi
!
!    The code uses Chebyshev expansions whose coefficients are
!    given to 20D.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) STRUVE_H0, the value of the function.
!
  implicit none

  real ( kind = 8 ) arrh0(0:19)
  real ( kind = 8 ) arrh0a(0:20)
  real ( kind = 8 ) ay0asp(0:12)
  real ( kind = 8 ) ay0asq(0:13)
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer indsgn
  integer, parameter :: nterm1 = 18
  integer, parameter :: nterm2 = 18
  integer, parameter :: nterm3 = 11
  integer, parameter :: nterm4 = 11
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) struve_h0
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) eleven,h0as, &
       piby4,rt2bpi,sixtp5,t,thr2p5,twenty, &
       twobpi,two62,xhigh,xlow,xmp4,xsq, &
       y0p,y0q,y0val
  data eleven/ 11.0d0/
  data twenty /20.0d0 /
  data sixtp5,two62,thr2p5/60.5d0, 262.0d0, 302.5d0/
  data piby4/0.78539816339744830962d0/
  data rt2bpi/0.79788456080286535588d0/
  data twobpi/0.63661977236758134308d0/
  data arrh0(0)/  0.28696487399013225740d0/
  data arrh0(1)/ -0.25405332681618352305d0/
  data arrh0(2)/  0.20774026739323894439d0/
  data arrh0(3)/ -0.20364029560386585140d0/
  data arrh0(4)/  0.12888469086866186016d0/
  data arrh0(5)/ -0.4825632815622261202d-1/
  data arrh0(6)/  0.1168629347569001242d-1/
  data arrh0(7)/ -0.198118135642418416d-2/
  data arrh0(8)/  0.24899138512421286d-3/
  data arrh0(9)/ -0.2418827913785950d-4/
  data arrh0(10)/ 0.187437547993431d-5/
  data arrh0(11)/-0.11873346074362d-6/
  data arrh0(12)/ 0.626984943346d-8/
  data arrh0(13)/-0.28045546793d-9/
  data arrh0(14)/ 0.1076941205d-10/
  data arrh0(15)/-0.35904793d-12/
  data arrh0(16)/ 0.1049447d-13/
  data arrh0(17)/-0.27119d-15/
  data arrh0(18)/ 0.624d-17/
  data arrh0(19)/-0.13d-18/
  data arrh0a(0)/  1.99291885751992305515d0/
  data arrh0a(1)/ -0.384232668701456887d-2/
  data arrh0a(2)/ -0.32871993712353050d-3/
  data arrh0a(3)/ -0.2941181203703409d-4/
  data arrh0a(4)/ -0.267315351987066d-5/
  data arrh0a(5)/ -0.24681031075013d-6/
  data arrh0a(6)/ -0.2295014861143d-7/
  data arrh0a(7)/ -0.215682231833d-8/
  data arrh0a(8)/ -0.20303506483d-9/
  data arrh0a(9)/ -0.1934575509d-10/
  data arrh0a(10)/-0.182773144d-11/
  data arrh0a(11)/-0.17768424d-12/
  data arrh0a(12)/-0.1643296d-13/
  data arrh0a(13)/-0.171569d-14/
  data arrh0a(14)/-0.13368d-15/
  data arrh0a(15)/-0.2077d-16/
  data arrh0a(16)/ 0.2d-19/
  data arrh0a(17)/-0.55d-18/
  data arrh0a(18)/ 0.10d-18/
  data arrh0a(19)/-0.4d-19/
  data arrh0a(20)/ 0.1d-19/
  data ay0asp/1.99944639402398271568d0, &
             -0.28650778647031958d-3, &
             -0.1005072797437620d-4, &
             -0.35835941002463d-6, &
             -0.1287965120531d-7, &
             -0.46609486636d-9, &
             -0.1693769454d-10, &
             -0.61852269d-12, &
             -0.2261841d-13, &
             -0.83268d-15, &
             -0.3042d-16, &
             -0.115d-17, &
             -0.4d-19/
  data ay0asq/1.99542681386828604092d0, &
             -0.236013192867514472d-2, &
             -0.7601538908502966d-4, &
             -0.256108871456343d-5, &
             -0.8750292185106d-7, &
             -0.304304212159d-8, &
             -0.10621428314d-9, &
             -0.377371479d-11, &
             -0.13213687d-12, &
             -0.488621d-14, &
             -0.15809d-15, &
             -0.762d-17, &
             -0.3d-19, &
             -0.3d-19/
!
!   MACHINE-DEPENDENT CONSTANTS (Suitable for IEEE-arithmetic machines)
!
  data xlow,xhigh/3.1610136d-8,4.50359963d15/

  x = xvalue
  indsgn = 1

  if ( x < zero ) then
    x = -x
    indsgn = -1
  end if

  if ( x < xlow ) then
    struve_h0 = twobpi * x
  else if ( x <= eleven ) then
    t = ( ( x * x ) / sixtp5 - half ) - half
    struve_h0 = twobpi * x * cheval ( nterm1, arrh0, t )
  else if ( x <= xhigh ) then
    xsq = x * x
    t = ( two62 - xsq ) / ( twenty + xsq )
    y0p = cheval ( nterm3, ay0asp, t )
    y0q = cheval ( nterm4, ay0asq, t ) / ( eight * x )
    xmp4 = x - piby4
    y0val = y0p * sin ( xmp4 ) - y0q * cos ( xmp4 )
    y0val = y0val * rt2bpi / sqrt ( x )
    t = ( thr2p5 - xsq ) / ( sixtp5 + xsq )
    h0as = twobpi * cheval ( nterm2, arrh0a, t ) / x
    struve_h0 = y0val + h0as
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRUVE_H0 - Fatal error!'
    write ( *, '(a)' ) '  Argument magnitude too large.'
    struve_h0 = zero
  end if

  if ( indsgn == -1 ) then
    struve_h0 = -struve_h0
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
!      HO(x) = (2/pi) * Integral ( 0 <= t <= pi/2 ) sin ( x * cos ( t ) ) dt
!
!    In Mathematica, the function can be evaluated by:
!
!      StruveH[0,x]
!
!    The data was reported by McLeod.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
!
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function struve_h1 ( xvalue )

!*****************************************************************************80
!
!! STRUVE_H1 calculates the Struve function of order 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      H1(x) = 2*x/pi * Integral ( 0 <= t <= pi/2 )
!        sin ( x * cos ( t ) )^2 * sin ( t ) dt
!
!    H1 also satisfies the second-order differential equation
!
!      x^2 * D^2 f  +  x * Df  +  (x^2 - 1)f  =  2x^2 / pi
!
!    The code uses Chebyshev expansions with the coefficients
!    given to 20D.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) STRUVE_H1, the value of the function.
!
  implicit none

  real ( kind = 8 ) arrh1(0:17)
  real ( kind = 8 ) arrh1a(0:21)
  real ( kind = 8 ) ay1asp(0:14)
  real ( kind = 8 ) ay1asq(0:15)
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer, parameter :: nterm1 = 15
  integer, parameter :: nterm2 = 17
  integer, parameter :: nterm3 = 12
  integer, parameter :: nterm4 = 13
  real ( kind = 8 ) struve_h1
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) fortp5, &
       h1as,nine,one82,rt2bpi,t,thpby4, &
       twenty,twobpi,tw02p5,xhigh,xlow1,xlow2, &
       xm3p4,xsq,y1p,y1q,y1val
  data nine / 9.0d0 /
  data twenty / 20.0d0 /
  data fortp5,one82,tw02p5/40.5d0, 182.0d0 , 202.5d0/
  data rt2bpi/0.79788456080286535588d0/
  data thpby4/2.35619449019234492885d0/
  data twobpi/0.63661977236758134308d0/
  data arrh1/0.17319061083675439319d0, &
            -0.12606917591352672005d0, &
             0.7908576160495357500d-1, &
            -0.3196493222321870820d-1, &
             0.808040581404918834d-2, &
            -0.136000820693074148d-2, &
             0.16227148619889471d-3, &
            -0.1442352451485929d-4, &
             0.99219525734072d-6, &
            -0.5441628049180d-7, &
             0.243631662563d-8, &
            -0.9077071338d-10, &
             0.285926585d-11, &
            -0.7716975d-13, &
             0.180489d-14, &
            -0.3694d-16, &
             0.67d-18, &
            -0.1d-19/
  data arrh1a(0)/  2.01083504951473379407d0/
  data arrh1a(1)/  0.592218610036099903d-2/
  data arrh1a(2)/  0.55274322698414130d-3/
  data arrh1a(3)/  0.5269873856311036d-4/
  data arrh1a(4)/  0.506374522140969d-5/
  data arrh1a(5)/  0.49028736420678d-6/
  data arrh1a(6)/  0.4763540023525d-7/
  data arrh1a(7)/  0.465258652283d-8/
  data arrh1a(8)/  0.45465166081d-9/
  data arrh1a(9)/  0.4472462193d-10/
  data arrh1a(10)/ 0.437308292d-11/
  data arrh1a(11)/ 0.43568368d-12/
  data arrh1a(12)/ 0.4182190d-13/
  data arrh1a(13)/ 0.441044d-14/
  data arrh1a(14)/ 0.36391d-15/
  data arrh1a(15)/ 0.5558d-16/
  data arrh1a(16)/-0.4d-19/
  data arrh1a(17)/ 0.163d-17/
  data arrh1a(18)/-0.34d-18/
  data arrh1a(19)/ 0.13d-18/
  data arrh1a(20)/-0.4d-19/
  data arrh1a(21)/ 0.1d-19/
  data ay1asp/2.00135240045889396402d0, &
              0.71104241596461938d-3, &
              0.3665977028232449d-4, &
              0.191301568657728d-5, &
              0.10046911389777d-6, &
              0.530401742538d-8, &
              0.28100886176d-9, &
              0.1493886051d-10, &
              0.79578420d-12, &
              0.4252363d-13, &
              0.227195d-14, &
              0.12216d-15, &
              0.650d-17, &
              0.36d-18, &
              0.2d-19/
  data ay1asq/5.99065109477888189116d0, &
             -0.489593262336579635d-2, &
             -0.23238321307070626d-3, &
             -0.1144734723857679d-4, &
             -0.57169926189106d-6, &
             -0.2895516716917d-7, &
             -0.147513345636d-8, &
             -0.7596537378d-10, &
             -0.390658184d-11, &
             -0.20464654d-12, &
             -0.1042636d-13, &
             -0.57702d-15, &
             -0.2550d-16, &
             -0.210d-17, &
              0.2d-19, &
             -0.2d-19/
!
!   MACHINE-DEPENDENT CONSTANTS (Suitable for IEEE-arithmetic machines)
!
  data xlow1,xlow2/2.23750222d-154,4.08085106d-8/
  data xhigh/4.50359963d15/
!
  x = abs ( xvalue )

  if ( x < xlow1 ) then
    struve_h1 = zero
  else if ( x < xlow2 ) then
    xsq = x * x
    struve_h1 = twobpi * xsq
  else if ( x <= nine ) then
    xsq = x * x
    t = ( xsq / fortp5 - half ) - half
    struve_h1 = twobpi * xsq * cheval ( nterm1, arrh1, t )
  else if ( x <= xhigh ) then
    xsq = x * x
    t = ( one82 - xsq ) / ( twenty + xsq )
    y1p = cheval ( nterm3, ay1asp, t )
    y1q = cheval ( nterm4, ay1asq, t ) / ( eight * x)
    xm3p4 = x - thpby4
    y1val = y1p * sin ( xm3p4 ) + y1q * cos ( xm3p4 )
    y1val = y1val * rt2bpi / sqrt ( x )
    t = ( tw02p5 - xsq ) / ( fortp5 + xsq )
    h1as = twobpi * cheval ( nterm2, arrh1a, t )
    struve_h1 = y1val + h1as
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRUVE_H1 - Fatal error!'
    write ( *, '(a)' ) '  Argument magnitude too large.'
    struve_h1 = zero
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
!      H1(x) = 2*x/pi * Integral ( 0 <= t <= pi/2 )
!        sin ( x * cos ( t ) )^2 * sin ( t ) dt
!
!    In Mathematica, the function can be evaluated by:
!
!      StruveH[1,x]
!
!    The data was reported by McLeod.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function struve_l0 ( xvalue )

!*****************************************************************************80
!
!! STRUVE_L0 calculates the modified Struve function of order 0.
!
!  Discussion:
!
!    This function calculates the modified Struve function of
!    order 0, denoted L0(x), defined as the solution of the
!    second-order equation
!
!      x*D(Df) + Df - x*f  =  2x/pi
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) STRUVE_L0, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: four = 4.0D+00
  integer indsgn
  integer, parameter :: nterm1 = 25
  integer, parameter :: nterm2 = 14
  integer, parameter :: nterm3 = 21
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) struve_l0
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) arl0(0:27),arl0as(0:15),ai0ml0(0:23), &
       atehun,ch1,ch2,lnr2pi, &
       sixten,t,test,twent4,twent8,twobpi,two88, &
       xhigh1,xhigh2,xlow,xmax,xsq
  data sixten/16.0d0/
  data twent4,twent8/24.0d0 , 28.0d0 /
  data two88,atehun/288.0d0 , 800.0d0/
  data lnr2pi/0.91893853320467274178d0/
  data twobpi/0.63661977236758134308d0/
  data arl0(0)/  0.42127458349979924863d0/
  data arl0(1)/ -0.33859536391220612188d0/
  data arl0(2)/  0.21898994812710716064d0/
  data arl0(3)/ -0.12349482820713185712d0/
  data arl0(4)/  0.6214209793866958440d-1/
  data arl0(5)/ -0.2817806028109547545d-1/
  data arl0(6)/  0.1157419676638091209d-1/
  data arl0(7)/ -0.431658574306921179d-2/
  data arl0(8)/  0.146142349907298329d-2/
  data arl0(9)/ -0.44794211805461478d-3/
  data arl0(10)/ 0.12364746105943761d-3/
  data arl0(11)/-0.3049028334797044d-4/
  data arl0(12)/ 0.663941401521146d-5/
  data arl0(13)/-0.125538357703889d-5/
  data arl0(14)/ 0.20073446451228d-6/
  data arl0(15)/-0.2588260170637d-7/
  data arl0(16)/ 0.241143742758d-8/
  data arl0(17)/-0.10159674352d-9/
  data arl0(18)/-0.1202430736d-10/
  data arl0(19)/ 0.262906137d-11/
  data arl0(20)/-0.15313190d-12/
  data arl0(21)/-0.1574760d-13/
  data arl0(22)/ 0.315635d-14/
  data arl0(23)/-0.4096d-16/
  data arl0(24)/-0.3620d-16/
  data arl0(25)/ 0.239d-17/
  data arl0(26)/ 0.36d-18/
  data arl0(27)/-0.4d-19/
  data arl0as(0)/  2.00861308235605888600d0/
  data arl0as(1)/  0.403737966500438470d-2/
  data arl0as(2)/ -0.25199480286580267d-3/
  data arl0as(3)/  0.1605736682811176d-4/
  data arl0as(4)/ -0.103692182473444d-5/
  data arl0as(5)/  0.6765578876305d-7/
  data arl0as(6)/ -0.444999906756d-8/
  data arl0as(7)/  0.29468889228d-9/
  data arl0as(8)/ -0.1962180522d-10/
  data arl0as(9)/  0.131330306d-11/
  data arl0as(10)/-0.8819190d-13/
  data arl0as(11)/ 0.595376d-14/
  data arl0as(12)/-0.40389d-15/
  data arl0as(13)/ 0.2651d-16/
  data arl0as(14)/-0.208d-17/
  data arl0as(15)/ 0.11d-18/
  data ai0ml0(0)/ 2.00326510241160643125d0/
  data ai0ml0(1)/ 0.195206851576492081d-2/
  data ai0ml0(2)/ 0.38239523569908328d-3/
  data ai0ml0(3)/ 0.7534280817054436d-4/
  data ai0ml0(4)/ 0.1495957655897078d-4/
  data ai0ml0(5)/ 0.299940531210557d-5/
  data ai0ml0(6)/ 0.60769604822459d-6/
  data ai0ml0(7)/ 0.12399495544506d-6/
  data ai0ml0(8)/ 0.2523262552649d-7/
  data ai0ml0(9)/ 0.504634857332d-8/
  data ai0ml0(10)/0.97913236230d-9/
  data ai0ml0(11)/0.18389115241d-9/
  data ai0ml0(12)/0.3376309278d-10/
  data ai0ml0(13)/0.611179703d-11/
  data ai0ml0(14)/0.108472972d-11/
  data ai0ml0(15)/0.18861271d-12/
  data ai0ml0(16)/0.3280345d-13/
  data ai0ml0(17)/0.565647d-14/
  data ai0ml0(18)/0.93300d-15/
  data ai0ml0(19)/0.15881d-15/
  data ai0ml0(20)/0.2791d-16/
  data ai0ml0(21)/0.389d-17/
  data ai0ml0(22)/0.70d-18/
  data ai0ml0(23)/0.16d-18/
!
!   MACHINE-DEPENDENT VALUES (Suitable for IEEE-arithmetic machines)
!
  data xlow,xmax/4.4703484d-8,1.797693d308/
  data xhigh1,xhigh2/5.1982303d8,2.5220158d17/

  x = xvalue
  indsgn = 1
  if ( x < zero ) then
    x = -x
    indsgn = -1
  end if

  if ( x < xlow ) then
    struve_l0 = twobpi * x
  else if ( x <= sixten ) then
    t = ( four * x - twent4 ) / ( x + twent4 )
    struve_l0 = twobpi * x * cheval ( nterm1, arl0, t ) * exp ( x )
  else
!
!   Code for |xvalue| > 16
!
     if ( xhigh2 < x ) then
        ch1 = one
     else
        t = ( x - twent8 ) / ( four - x )
        ch1 = cheval ( nterm2, arl0as, t )
     end if

     if ( xhigh1 < x ) then
        ch2 = one
     else
        xsq = x * x
        t = ( atehun - xsq ) / ( two88 + xsq )
        ch2 = cheval ( nterm3, ai0ml0, t )
     end if

     test = log ( ch1 ) - lnr2pi - log ( x ) / two + x

     if ( log ( xmax ) < test ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRUVE_L0 - Fatal error!'
        write ( *, '(a)' ) '  Argument would cause overflow.'
        struve_l0 = xmax
     else
        struve_l0 = exp ( test ) - twobpi * ch2 / x
     end if

  end if

  if ( indsgn == -1 ) then
    struve_l0 = -struve_l0
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
!  Modified:
!
!    03 September 2004
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
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function struve_l1 ( xvalue )

!*****************************************************************************80
!
!! STRUVE_L1 calculates the modified Struve function of order 1.
!
!  Discussion:
!
!    This function calculates the modified Struve function of
!    order 1, denoted L1(x), defined as the solution of
!
!      x*x*D(Df) + x*Df - (x*x+1)f = 2 * x * x / pi
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) STRUVE_L1, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: four = 4.0D+00
  integer, parameter :: nterm1 = 24
  integer, parameter :: nterm2 = 13
  integer, parameter :: nterm3 = 22
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: sixten = 16.0D+00
  real ( kind = 8 ) struve_l1
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) arl1(0:26),arl1as(0:16),ai1ml1(0:25), &
       atehun,ch1,ch2,lnr2pi, &
       pi3by2,t,test,thirty,twent4, &
       twobpi,two88,xhigh1,xhigh2,xlow1,xlow2, &
       xmax,xsq
  data twent4,thirty/24.0d0, 30.0d0/
  data two88,atehun/288.0d0, 800.0d0/
  data lnr2pi/0.91893853320467274178d0/
  data pi3by2/4.71238898038468985769d0/
  data twobpi/0.63661977236758134308d0/
  data arl1(0)/  0.38996027351229538208d0/
  data arl1(1)/ -0.33658096101975749366d0/
  data arl1(2)/  0.23012467912501645616d0/
  data arl1(3)/ -0.13121594007960832327d0/
  data arl1(4)/  0.6425922289912846518d-1/
  data arl1(5)/ -0.2750032950616635833d-1/
  data arl1(6)/  0.1040234148637208871d-1/
  data arl1(7)/ -0.350532294936388080d-2/
  data arl1(8)/  0.105748498421439717d-2/
  data arl1(9)/ -0.28609426403666558d-3/
  data arl1(10)/ 0.6925708785942208d-4/
  data arl1(11)/-0.1489693951122717d-4/
  data arl1(12)/ 0.281035582597128d-5/
  data arl1(13)/-0.45503879297776d-6/
  data arl1(14)/ 0.6090171561770d-7/
  data arl1(15)/-0.623543724808d-8/
  data arl1(16)/ 0.38430012067d-9/
  data arl1(17)/ 0.790543916d-11/
  data arl1(18)/-0.489824083d-11/
  data arl1(19)/ 0.46356884d-12/
  data arl1(20)/ 0.684205d-14/
  data arl1(21)/-0.569748d-14/
  data arl1(22)/ 0.35324d-15/
  data arl1(23)/ 0.4244d-16/
  data arl1(24)/-0.644d-17/
  data arl1(25)/-0.21d-18/
  data arl1(26)/ 0.9d-19/
  data arl1as(0)/  1.97540378441652356868d0/
  data arl1as(1)/ -0.1195130555088294181d-1/
  data arl1as(2)/  0.33639485269196046d-3/
  data arl1as(3)/ -0.1009115655481549d-4/
  data arl1as(4)/  0.30638951321998d-6/
  data arl1as(5)/ -0.953704370396d-8/
  data arl1as(6)/  0.29524735558d-9/
  data arl1as(7)/ -0.951078318d-11/
  data arl1as(8)/  0.28203667d-12/
  data arl1as(9)/ -0.1134175d-13/
  data arl1as(10)/ 0.147d-17/
  data arl1as(11)/-0.6232d-16/
  data arl1as(12)/-0.751d-17/
  data arl1as(13)/-0.17d-18/
  data arl1as(14)/ 0.51d-18/
  data arl1as(15)/ 0.23d-18/
  data arl1as(16)/ 0.5d-19/
  data ai1ml1(0)/  1.99679361896789136501d0/
  data ai1ml1(1)/ -0.190663261409686132d-2/
  data ai1ml1(2)/ -0.36094622410174481d-3/
  data ai1ml1(3)/ -0.6841847304599820d-4/
  data ai1ml1(4)/ -0.1299008228509426d-4/
  data ai1ml1(5)/ -0.247152188705765d-5/
  data ai1ml1(6)/ -0.47147839691972d-6/
  data ai1ml1(7)/ -0.9020819982592d-7/
  data ai1ml1(8)/ -0.1730458637504d-7/
  data ai1ml1(9)/ -0.332323670159d-8/
  data ai1ml1(10)/-0.63736421735d-9/
  data ai1ml1(11)/-0.12180239756d-9/
  data ai1ml1(12)/-0.2317346832d-10/
  data ai1ml1(13)/-0.439068833d-11/
  data ai1ml1(14)/-0.82847110d-12/
  data ai1ml1(15)/-0.15562249d-12/
  data ai1ml1(16)/-0.2913112d-13/
  data ai1ml1(17)/-0.543965d-14/
  data ai1ml1(18)/-0.101177d-14/
  data ai1ml1(19)/-0.18767d-15/
  data ai1ml1(20)/-0.3484d-16/
  data ai1ml1(21)/-0.643d-17/
  data ai1ml1(22)/-0.118d-17/
  data ai1ml1(23)/-0.22d-18/
  data ai1ml1(24)/-0.4d-19/
  data ai1ml1(25)/-0.1d-19/
!
!   MACHINE-DEPENDENT VALUES (Suitable for IEEE-arithmetic machines)
!
  data xlow1,xlow2,xmax/5.7711949d-8,3.3354714d-154,1.797693d308/
  data xhigh1,xhigh2/5.19823025d8,2.7021597d17/

  x = abs ( xvalue )

  if ( x <= xlow2 ) then
    struve_l1 = zero
  else if ( x < xlow1 ) then
    xsq = x * x
    struve_l1 = xsq / pi3by2
  else if ( x <= sixten ) then
    xsq = x * x
    t = ( four * x - twent4 ) / ( x + twent4 )
    struve_l1 = xsq * cheval ( nterm1, arl1, t ) * exp ( x ) / pi3by2
  else

    if ( xhigh2 < x ) then
      ch1 = one
    else
      t = ( x - thirty ) / ( two - x )
      ch1 = cheval ( nterm2, arl1as, t )
    end if

    if ( xhigh1 < x ) then
      ch2 = one
    else
      xsq = x * x
      t = ( atehun - xsq ) / ( two88 + xsq )
      ch2 = cheval ( nterm3, ai1ml1, t )
    end if

    test = log ( ch1 ) - lnr2pi - log ( x ) / two + x

    if ( log ( xmax ) < test ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STRUVE_L1 - Fatal error!'
      write ( *, '(a)' ) '  Argument would cause overflow.'
      struve_l1 = xmax
    else
      struve_l1 = exp ( test ) - twobpi * ch2
    end if

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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
           0.80950410749865126939d-06, &
           0.20724649092571514607d-03, &
           0.33191834066894516744d-02, &
           0.53942182623522663292d-01, &
           0.22676438105580863683d+00, &
           0.11027597873677158176d+01, &
           0.91692778117386847344d+01, &
           0.15541656652426660966d+03, &
           0.26703582852084829694d+04, &
           0.86505880175304633906d+06, &
           0.11026046613094942620d+07, &
           0.22846209494153934787d+07, &
           0.42454972750111979449d+08, &
           0.48869614587997695539d+09, &
           0.56578651292431051863d+10, &
           0.76853203893832108948d+12, &
           0.14707396163259352103d+17, &
           0.29030785901035567967d+21, &
           0.58447515883904682813d+25, &
           0.11929750788892311875d+30 /)
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function synch1 ( xvalue )

!*****************************************************************************80
!
!! SYNCH1 calculates the synchrotron radiation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      SYNCH1(x) = x * Integral ( x <= t < infinity ) K(5/3)(t) dt
!
!    where K(5/3) is a modified Bessel function of order 5/3.
!
!    The code uses Chebyshev expansions, the coefficients of which
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) SYNCH1, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer, parameter :: nterm1 = 12
  integer, parameter :: nterm2 = 10
  integer, parameter :: nterm3 = 21
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) synch1
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) async1(0:13),async2(0:11),asynca(0:24), &
       cheb1,cheb2,conlow, &
       lnrtp2,pibrt3,t,twelve,xhigh1, &
       xhigh2,xlow,xpowth
  data twelve/ 12.0d0 /
  data conlow/2.14952824153447863671d0/
  data pibrt3/1.81379936423421785059d0/
  data lnrtp2/0.22579135264472743236d0/
  data async1/30.36468298250107627340d0, &
              17.07939527740839457449d0, &
               4.56013213354507288887d0, &
               0.54928124673041997963d0, &
               0.3729760750693011724d-1, &
               0.161362430201041242d-2, &
               0.4819167721203707d-4, &
               0.105124252889384d-5, &
               0.1746385046697d-7, &
               0.22815486544d-9, &
               0.240443082d-11, &
               0.2086588d-13, &
               0.15167d-15, &
               0.94d-18/
  data async2/0.44907216235326608443d0, &
              0.8983536779941872179d-1, &
              0.810445737721512894d-2, &
              0.42617169910891619d-3, &
              0.1476096312707460d-4, &
              0.36286336153998d-6, &
              0.666348074984d-8, &
              0.9490771655d-10, &
              0.107912491d-11, &
              0.1002201d-13, &
              0.7745d-16, &
              0.51d-18/
  data asynca(0)/ 2.13293051613550009848d0/
  data asynca(1)/ 0.7413528649542002401d-1/
  data asynca(2)/ 0.869680999099641978d-2/
  data asynca(3)/ 0.117038262487756921d-2/
  data asynca(4)/ 0.16451057986191915d-3/
  data asynca(5)/ 0.2402010214206403d-4/
  data asynca(6)/ 0.358277563893885d-5/
  data asynca(7)/ 0.54477476269837d-6/
  data asynca(8)/ 0.8388028561957d-7/
  data asynca(9)/ 0.1306988268416d-7/
  data asynca(10)/0.205309907144d-8/
  data asynca(11)/0.32518753688d-9/
  data asynca(12)/0.5179140412d-10/
  data asynca(13)/0.830029881d-11/
  data asynca(14)/0.133527277d-11/
  data asynca(15)/0.21591498d-12/
  data asynca(16)/0.3499673d-13/
  data asynca(17)/0.569942d-14/
  data asynca(18)/0.92906d-15/
  data asynca(19)/0.15222d-15/
  data asynca(20)/0.2491d-16/
  data asynca(21)/0.411d-17/
  data asynca(22)/0.67d-18/
  data asynca(23)/0.11d-18/
  data asynca(24)/0.2d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow/2.98023224d-8/
  data xhigh1,xhigh2/809.595907d0,-708.396418d0/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SYNCH1 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    synch1 = zero
  else if ( x < xlow ) then
    xpowth = x ** ( one / three )
    synch1 = conlow * xpowth
  else if ( x <= four ) then
    xpowth = x ** ( one / three )
    t = ( x * x / eight - half ) - half
    cheb1 = cheval ( nterm1, async1, t )
    cheb2 = cheval ( nterm2, async2, t )
    t = xpowth * cheb1 - xpowth**11 * cheb2
    synch1 = t - pibrt3 * x
  else if ( x <= xhigh1 ) then
    t = ( twelve - x ) / ( x + four )
    cheb1 = cheval ( nterm3, asynca, t )
    t = lnrtp2 - x + log ( sqrt ( x ) * cheb1 )
    if ( t < xhigh2 ) then
      synch1 = zero
    else
      synch1 = exp ( t )
    end if
  else
    synch1 = zero
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
!      SYNCH1(x) = x * Integral ( x <= t < infinity ) K(5/3)(t) dt
!
!    where K(5/3) is a modified Bessel function of order 5/3.
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function synch2 ( xvalue )

!*****************************************************************************80
!
!! SYNCH2 calculates the synchrotron radiation function.
!
!  Discussion:
!
!    The function is defined by:
!
!      SYNCH2(x) = x * K(2/3)(x)
!
!    where K(2/3) is a modified Bessel function of order 2/3.
!
!    The code uses Chebyshev expansions, the coefficients of which
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) SYNCH2, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer, parameter :: nterm1 = 13
  integer, parameter :: nterm2 = 12
  integer, parameter :: nterm3 = 16
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) synch2
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) asyn21(0:14),asyn22(0:13),asyn2a(0:18), &
       cheb1,cheb2,conlow, &
       lnrtp2,t,ten,xhigh1, &
       xhigh2,xlow,xpowth
  data ten/ 10.0d0 /
  data conlow/1.07476412076723931836d0/
  data lnrtp2/0.22579135264472743236d0/
  data asyn21/38.61783992384308548014d0, &
              23.03771559496373459697d0, &
               5.38024998683357059676d0, &
               0.61567938069957107760d0, &
               0.4066880046688955843d-1, &
               0.172962745526484141d-2, &
               0.5106125883657699d-4, &
               0.110459595022012d-5, &
               0.1823553020649d-7, &
               0.23707698034d-9, &
               0.248872963d-11, &
               0.2152868d-13, &
               0.15607d-15, &
               0.96d-18, &
               0.1d-19/
  data asyn22/7.90631482706608042875d0, &
              3.13534636128534256841d0, &
              0.48548794774537145380d0, &
              0.3948166758272372337d-1, &
              0.196616223348088022d-2, &
              0.6590789322930420d-4, &
              0.158575613498559d-5, &
              0.2868653011233d-7, &
              0.40412023595d-9, &
              0.455684443d-11, &
              0.4204590d-13, &
              0.32326d-15, &
              0.210d-17, &
              0.1d-19/
  data asyn2a/2.02033709417071360032d0, &
              0.1095623712180740443d-1, &
              0.85423847301146755d-3, &
              0.7234302421328222d-4, &
              0.631244279626992d-5, &
              0.56481931411744d-6, &
              0.5128324801375d-7, &
              0.471965329145d-8, &
              0.43807442143d-9, &
              0.4102681493d-10, &
              0.386230721d-11, &
              0.36613228d-12, &
              0.3480232d-13, &
              0.333010d-14, &
              0.31856d-15, &
              0.3074d-16, &
              0.295d-17, &
              0.29d-18, &
              0.3d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow/2.98023224d-8/
  data xhigh1,xhigh2/809.595907d0,-708.396418d0/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SYNCH2 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    synch2 = zero
  else if ( x < xlow ) then
    xpowth = x ** ( one / three )
    synch2 = conlow * xpowth
  else if ( x <= four ) then
    xpowth = x ** ( one / three )
    t = ( x * x / eight - half ) - half
    cheb1 = cheval ( nterm1, asyn21, t )
    cheb2 = cheval ( nterm2, asyn22, t )
    synch2 = xpowth * cheb1 - xpowth**5 * cheb2
  else if ( x <= xhigh1 ) then
    t = ( ten - x ) / ( x + two )
    cheb1 = cheval ( nterm3, asyn2a, t )
    t = lnrtp2 - x + log ( sqrt ( x ) * cheb1 )
    if ( t < xhigh2 ) then
      synch2 = zero
    else
      synch2 = exp ( t )
    end if
  else
    synch2 = zero
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
function tran02 ( xvalue )

!*****************************************************************************80
!
!! TRAN02 calculates the transport integral of order 2.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN02(x) = Integral ( 0 <= t <= x ) t^2 exp(t) / ( exp(t) - 1 )^2 dt
!
!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) TRAN02, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer k1
  integer k2
  integer, parameter :: nterms = 17
  integer numexp
  integer, parameter :: numjn = 2
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tran02
  real ( kind = 8 ), parameter :: valinf = 0.32898681336964528729D+01
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) atran(0:19),rk, &
       rnumjn,sumexp,sum2,t,xhigh1,xhigh2, &
       xhigh3,xk,xk1,xlow1
  data rnumjn/ 2.0d0 /

  data atran/1.67176044643453850301d0, &
            -0.14773535994679448986d0, &
             0.1482138199469363384d-1, &
            -0.141953303263056126d-2, &
             0.13065413244157083d-3, &
            -0.1171557958675790d-4, &
             0.103334984457557d-5, &
            -0.9019113042227d-7, &
             0.781771698331d-8, &
            -0.67445656840d-9, &
             0.5799463945d-10, &
            -0.497476185d-11, &
             0.42596097d-12, &
            -0.3642189d-13, &
             0.311086d-14, &
            -0.26547d-15, &
             0.2264d-16, &
            -0.193d-17, &
             0.16d-18, &
            -0.1d-19/
!
!  Machine-dependent constants
!
  data xlow1/2.98023224d-8/
  data xhigh1,xhigh3/36.04365668d0,-36.73680056d0/
  data xhigh2/9.00719925d15/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAN02 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    tran02 = zero
  else if ( x < xlow1 ) then
    tran02 =  ( x ** ( numjn - 1 ) ) / ( rnumjn - one )
  else if ( x <= four ) then
    t = ( ( ( x * x ) / eight ) - half ) - half
    tran02 = ( x ** ( numjn - 1 ) ) * cheval ( nterms, atran, t )
  else

     if ( xhigh2 < x ) then
        sumexp = one
     else
        if ( x <= xhigh1 ) then
           numexp = int ( xhigh1 / x ) + 1
           t = exp ( -x )
        else
           numexp = 1
           t = one
        end if
        rk = zero
        do k1 = 1, numexp
           rk = rk + one
        end do
        sumexp = zero
        do k1 = 1, numexp
           sum2 = one
           xk = one / ( rk * x )
           xk1 = one
           do k2 = 1, numjn
              sum2 = sum2 * xk1 * xk + one
              xk1 = xk1 + one
           end do
           sumexp = sumexp * t + sum2
           rk = rk - one
        end do
     end if

     t = rnumjn * log ( x ) - x + log ( sumexp )
     if ( t < xhigh3 ) then
        tran02 = valinf
     else
        tran02 = valinf - exp ( t )
     end if

  end if

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
!      TRAN02(x) = Integral ( 0 <= t <= x ) t^2 exp(t) / ( exp(t) - 1 )^2 dt
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tran03 ( xvalue )

!*****************************************************************************80
!
!! TRAN03 calculates the transport integral of order 3.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN03(x) = Integral ( 0 <= t <= x ) t^3 * exp(t) / ( exp(t) - 1 )^2 dt
!
!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) TRAN03, the value of the function.
!
  implicit none

  real ( kind = 8 ) atran(0:19)
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer k1
  integer k2
  integer, parameter :: nterms = 17
  integer numexp
  integer numjn
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tran03
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) rk, &
       rnumjn,sumexp,sum2,t,valinf,xhigh1,xhigh2, &
       xhigh3,xk,xk1,xlow1,xlow2
  data numjn,rnumjn/ 3 , 3.0d0 /
  data valinf/0.72123414189575657124d1/
  data atran/0.76201254324387200657d0, &
            -0.10567438770505853250d0, &
             0.1197780848196578097d-1, &
            -0.121440152036983073d-2, &
             0.11550997693928547d-3, &
            -0.1058159921244229d-4, &
             0.94746633853018d-6, &
            -0.8362212128581d-7, &
             0.731090992775d-8, &
            -0.63505947788d-9, &
             0.5491182819d-10, &
            -0.473213954d-11, &
             0.40676948d-12, &
            -0.3489706d-13, &
             0.298923d-14, &
            -0.25574d-15, &
             0.2186d-16, &
            -0.187d-17, &
             0.16d-18, &
            -0.1d-19/
!
!  Machine-dependent constants
!
  data xlow1,xlow2/2.98023224d-8,2.10953733d-154/
  data xhigh1,xhigh3/36.04365668d0,-36.73680056d0/
  data xhigh2/1.35107988d16/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAN03 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    tran03 = zero
  else if ( x < xlow2 ) then
    tran03 = zero
  else if ( x < xlow1 ) then
    tran03 = ( x**( numjn - 1 ) ) / ( rnumjn - one )
  else if ( x <= four ) then
    t = ( ( ( x*x ) / eight ) - half ) - half
    tran03 = ( x**( numjn - 1 ) ) * cheval ( nterms, atran, t )
  else

     if ( xhigh2 < x ) then
        sumexp = one
     else
        if ( x <= xhigh1 ) then
           numexp = int ( xhigh1 / x ) + 1
           t = exp ( -x )
        else
           numexp = 1
           t = one
        end if
        rk = zero
        do k1 = 1, numexp
           rk = rk + one
        end do
        sumexp = zero
        do k1 = 1, numexp
           sum2 = one
           xk = one / ( rk * x )
           xk1 = one
           do k2 = 1, numjn
              sum2 = sum2 * xk1 * xk + one
              xk1 = xk1 + one
           end do
           sumexp = sumexp * t + sum2
           rk = rk - one
        end do
     end if

     t = rnumjn * log ( x ) - x + log ( sumexp )

     if ( t < xhigh3 ) then
        tran03 = valinf
     else
        tran03 = valinf - exp ( t )
     end if

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
!      TRAN03(x) = Integral ( 0 <= t <= x ) t^3 * exp(t) / ( exp(t) - 1 )^2 dt
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tran04 ( xvalue )

!*****************************************************************************80
!
!! TRAN04 calculates the transport integral of order 4.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN04(x) = Integral ( 0 <= t <= x ) t^4 * exp(t) / ( exp(t) - 1 )^2 dt
!
!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) TRAN04, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer k1
  integer k2
  integer, parameter :: nterms = 17
  integer numexp
  integer numjn
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tran04
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) atran(0:19),rk, &
       rnumjn,sumexp,sum2,t,valinf,xhigh1,xhigh2, &
       xhigh3,xk,xk1,xlow1,xlow2
  data numjn,rnumjn/ 4 , 4.0d0 /
  data valinf/0.25975757609067316596d2/
  data atran/0.48075709946151105786d0, &
            -0.8175378810321083956d-1, &
             0.1002700665975162973d-1, &
            -0.105993393598201507d-2, &
             0.10345062450304053d-3, &
            -0.964427054858991d-5, &
             0.87455444085147d-6, &
            -0.7793212079811d-7, &
             0.686498861410d-8, &
            -0.59995710764d-9, &
             0.5213662413d-10, &
            -0.451183819d-11, &
             0.38921592d-12, &
            -0.3349360d-13, &
             0.287667d-14, &
            -0.24668d-15, &
             0.2113d-16, &
            -0.181d-17, &
             0.15d-18, &
            -0.1d-19/
!
!  Machine-dependent constants
!
  data xlow1,xlow2/2.98023224d-8,4.05653502d-103/
  data xhigh1,xhigh3/36.04365668d0,-36.73680056d0/
  data xhigh2/1.80143985d16/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAN04 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    tran04 = zero
    return
  end if
!
!   Code for x < =  4.0
!
  if ( x <= four ) then
     if ( x < xlow2 ) then
        tran04 = zero
     else
        if ( x < xlow1 ) then
           tran04 =  ( x ** ( numjn-1 ) ) / ( rnumjn - one )
        else
           t = ( ( ( x * x ) / eight ) - half ) - half
           tran04 = ( x ** ( numjn-1 ) ) * cheval ( nterms, atran, t )
        end if
     end if
  else
!
!  Code for x > 4.0
!
     if ( xhigh2 < x ) then
        sumexp = one
     else
        if ( x <= xhigh1 ) then
           numexp = int ( xhigh1 / x ) + 1
           t = exp ( -x )
        else
           numexp = 1
           t = one
        end if
        rk = zero
        do k1 = 1, numexp
           rk = rk + one
        end do
        sumexp = zero
        do k1 = 1, numexp
           sum2 = one
           xk = one / ( rk * x )
           xk1 = one
           do k2 = 1, numjn
              sum2 = sum2 * xk1 * xk + one
              xk1 = xk1 + one
           end do
           sumexp = sumexp * t + sum2
           rk = rk - one
        end do
     end if

     t = rnumjn * log ( x ) - x + log ( sumexp )

     if ( t < xhigh3 ) then
        tran04 = valinf
     else
        tran04 = valinf - exp ( t )
     end if

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
!      TRAN04(x) = Integral ( 0 <= t <= x ) t^4 * exp(t) / ( exp(t) - 1 )^2 dt
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tran05 ( xvalue )

!*****************************************************************************80
!
!! TRAN05 calculates the transport integral of order 5.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN05(x) = Integral ( 0 <= t <= x ) t^5 * exp(t) / ( exp(t) - 1 )^2 dt
!
!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) TRAN05, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer k1
  integer k2
  integer, parameter :: nterms = 17
  integer numexp
  integer numjn
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tran05
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) atran(0:19),rk, &
       rnumjn,sumexp,sum2,t,valinf,xhigh1,xhigh2, &
       xhigh3,xk,xk1,xlow1,xlow2
  data numjn,rnumjn/ 5 , 5.0d0 /
  data valinf/0.12443133061720439116d3/
  data atran/0.34777777713391078928d0, &
            -0.6645698897605042801d-1, &
             0.861107265688330882d-2, &
            -0.93966822237555384d-3, &
             0.9363248060815134d-4, &
            -0.885713193408328d-5, &
             0.81191498914503d-6, &
            -0.7295765423277d-7, &
             0.646971455045d-8, &
            -0.56849028255d-9, &
             0.4962559787d-10, &
            -0.431093996d-11, &
             0.37310094d-12, &
            -0.3219769d-13, &
             0.277220d-14, &
            -0.23824d-15, &
             0.2044d-16, &
            -0.175d-17, &
             0.15d-18, &
            -0.1d-19/
!
!  Machine-dependent constants
!
  data xlow1,xlow2/2.98023224d-8,1.72723372d-77/
  data xhigh1,xhigh3/36.04365668d0,-36.73680056d0/
  data xhigh2/2.25179981d16/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAN05 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    tran05 = zero
    return
  end if
!
!   Code for x < =  4.0
!
  if ( x <= four ) then
     if ( x < xlow2 ) then
        tran05 = zero
     else
        if ( x < xlow1 ) then
           tran05 =  ( x ** ( numjn - 1 ) ) / ( rnumjn - one )
        else
           t = ( ( ( x * x ) / eight ) - half ) - half
           tran05 = ( x ** ( numjn-1 ) ) * cheval ( nterms, atran, t )
        end if
     end if
  else
!
!  Code for x > 4.0
!
     if ( xhigh2 < x ) then
        sumexp = one
     else
        if ( x <= xhigh1 ) then
           numexp = int ( xhigh1 / x )  + 1
           t = exp ( -x )
        else
           numexp = 1
           t = one
        end if
        rk = zero
        do k1 = 1, numexp
           rk = rk + one
        end do
        sumexp = zero
        do k1 = 1, numexp
           sum2 = one
           xk = one / ( rk * x )
           xk1 = one
           do k2 = 1, numjn
              sum2 = sum2 * xk1 * xk + one
              xk1 = xk1 + one
           end do
           sumexp = sumexp * t + sum2
           rk = rk - one
        end do
     end if
     t = rnumjn * log ( x ) - x + log ( sumexp )
     if ( t < xhigh3 ) then
        tran05 = valinf
     else
        tran05 = valinf - exp ( t )
     end if
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
!      TRAN05(x) = Integral ( 0 <= t <= x ) t^5 * exp(t) / ( exp(t) - 1 )^2 dt
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tran06 ( xvalue )

!*****************************************************************************80
!
!! TRAN06 calculates the transport integral of order 6.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN06(x) = Integral ( 0 <= t <= x ) t^6 * exp(t) / ( exp(t) - 1 )^2 dt
!
!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) TRAN06, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer k1
  integer k2
  integer, parameter :: nterms = 17
  integer numexp
  integer numjn
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tran06
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) atran(0:19),rk, &
       rnumjn,sumexp,sum2,t,valinf,xhigh1,xhigh2, &
       xhigh3,xk,xk1,xlow1,xlow2
  data numjn,rnumjn/ 6 , 6.0d0 /
  data valinf/0.73248700462880338059d3/
  data atran/0.27127335397840008227d0, &
            -0.5588610553191453393d-1, &
             0.753919513290083056d-2, &
            -0.84351138579211219d-3, &
             0.8549098079676702d-4, &
            -0.818715493293098d-5, &
             0.75754240427986d-6, &
            -0.6857306541831d-7, &
             0.611700376031d-8, &
            -0.54012707024d-9, &
             0.4734306435d-10, &
            -0.412701055d-11, &
             0.35825603d-12, &
            -0.3099752d-13, &
             0.267501d-14, &
            -0.23036d-15, &
             0.1980d-16, &
            -0.170d-17, &
             0.15d-18, &
            -0.1d-19/
!
!  Machine-dependent constants
!
  data xlow1,xlow2/2.98023224d-8,4.06689432d-62/
  data xhigh1,xhigh3/36.04365668d0,-36.73680056d0/
  data xhigh2/2.70215977d16/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAN06 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    tran06 = zero
    return
  end if
!
!   Code for x < =  4 .0
!
  if ( x <= four ) then
     if ( x < xlow2 ) then
        tran06 = zero
     else
        if ( x < xlow1 ) then
           tran06 =  ( x ** ( numjn-1 ) ) / ( rnumjn - one )
        else
           t =  ( ( ( x * x ) / eight ) - half ) - half
           tran06 = ( x ** ( numjn-1 )  ) * cheval ( nterms, atran, t )
        end if
     end if
  else
!
!  Code for x > 4 .0
!
     if ( xhigh2 < x ) then
        sumexp = one
     else
        if ( x <= xhigh1 ) then
           numexp = int ( xhigh1 / x ) + 1
           t = exp ( - x )
        else
           numexp = 1
           t = one
        end if
        rk = zero
        do k1 = 1, numexp
           rk = rk + one
        end do
        sumexp = zero
        do k1 = 1, numexp
           sum2 = one
           xk = one / ( rk * x )
           xk1 = one
           do k2 = 1, numjn
              sum2 = sum2 * xk1 * xk + one
              xk1 = xk1 + one
           end do
           sumexp = sumexp * t + sum2
           rk = rk - one
        end do
     end if
     t = rnumjn * log ( x ) - x + log ( sumexp )
     if ( t < xhigh3 ) then
        tran06 = valinf
     else
        tran06 = valinf - exp ( t )
     end if
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
!      TRAN06(x) = Integral ( 0 <= t <= x ) t^6 * exp(t) / ( exp(t) - 1 )^2 dt
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tran07 ( xvalue )

!*****************************************************************************80
!
!! TRAN07 calculates the transport integral of order 7.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN07(x) = Integral ( 0 <= t <= x ) t^7 * exp(t) / ( exp(t) - 1 )^2 dt
!
!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) TRAN07, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer k1
  integer k2
  integer, parameter :: nterms = 17
  integer numexp
  integer numjn
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tran07
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) atran(0:19),rk, &
       rnumjn,sumexp,sum2,t,valinf,xhigh1,xhigh2, &
       xhigh3,xk,xk1,xlow1,xlow2
  data numjn,rnumjn/ 7 , 7.0d0/
  data valinf/0.50820803580048910473d4/
  data atran/0.22189250734010404423d0, &
            -0.4816751061177993694d-1, &
             0.670092448103153629d-2, &
            -0.76495183443082557d-3, &
             0.7863485592348690d-4, &
            -0.761025180887504d-5, &
             0.70991696299917d-6, &
            -0.6468025624903d-7, &
             0.580039233960d-8, &
            -0.51443370149d-9, &
             0.4525944183d-10, &
            -0.395800363d-11, &
             0.34453785d-12, &
            -0.2988292d-13, &
             0.258434d-14, &
            -0.22297d-15, &
             0.1920d-16, &
            -0.165d-17, &
             0.14d-18, &
            -0.1d-19/
!
!  Machine-dependent constants
!
  data xlow1,xlow2/2.98023224d-8,7.14906557d-52/
  data xhigh1,xhigh3/36.04365668d0,-36.73680056d0/
  data xhigh2/3.15251973d16/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAN07 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    tran07 = zero
    return
  end if
!
!   Code for x <= 4.0
!
  if ( x <= four ) then
     if ( x < xlow2 ) then
        tran07 = zero
     else
        if ( x < xlow1 ) then
           tran07 = ( x**(numjn-1) ) / ( rnumjn - one )
        else
           t = ( ( ( x * x ) / eight ) - half ) - half
           tran07 = ( x**(numjn-1) ) * cheval ( nterms, atran, t )
        end if
     end if
  else
!
!  Code for x > 4.0
!
     if ( xhigh2 < x ) then
        sumexp = one
     else
        if ( x <= xhigh1 ) then
           numexp = int ( xhigh1 / x ) + 1
           t = exp ( -x )
        else
           numexp = 1
           t = one
        end if
        rk = zero
        do k1 = 1, numexp
           rk = rk + one
        end do
        sumexp = zero
        do k1 = 1, numexp
           sum2 = one
           xk = one / ( rk * x )
           xk1 = one
           do k2 = 1, numjn
              sum2 = sum2 * xk1 * xk + one
              xk1 = xk1 + one
           end do
           sumexp = sumexp * t + sum2
           rk = rk - one
        end do
     end if

     t = rnumjn * log ( x ) - x + log ( sumexp )

     if ( t < xhigh3 ) then
        tran07 = valinf
     else
        tran07 = valinf - exp ( t )
     end if

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
!      TRAN07(x) = Integral ( 0 <= t <= x ) t^7 * exp(t) / ( exp(t) - 1 )^2 dt
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tran08 ( xvalue )

!*****************************************************************************80
!
!! TRAN08 calculates the transport integral of order 8.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN08(x) = Integral ( 0 <= t <= x ) t^8 * exp(t) / ( exp(t) - 1 )^2 dt
!
!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) TRAN08, the value of the function.
!
  implicit none

  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer k1
  integer k2
  integer, parameter :: nterms = 17
  integer numexp
  integer numjn
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tran08
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) atran(0:19),rk, &
       rnumjn,sumexp,sum2,t,valinf,xhigh1,xhigh2, &
       xhigh3,xk,xk1,xlow1,xlow2
  data numjn,rnumjn/ 8, 8.0d0 /
  data valinf/0.40484399001901115764d5/
  data atran/0.18750695774043719233d0, &
            -0.4229527646093673337d-1, &
             0.602814856929065592d-2, &
            -0.69961054811814776d-3, &
             0.7278482421298789d-4, &
            -0.710846250050067d-5, &
             0.66786706890115d-6, &
            -0.6120157501844d-7, &
             0.551465264474d-8, &
            -0.49105307052d-9, &
             0.4335000869d-10, &
            -0.380218700d-11, &
             0.33182369d-12, &
            -0.2884512d-13, &
             0.249958d-14, &
            -0.21605d-15, &
             0.1863d-16, &
            -0.160d-17, &
             0.14d-18, &
            -0.1d-19/
!
!  Machine-dependent constants
!
  data xlow1,xlow2/2.98023224d-8,1.48029723d-44/
  data xhigh1,xhigh3/36.04365668d0,-36.73680056d0/
  data xhigh2/3.6028797d16/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAN08 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    tran08 = zero
    return
  end if
!
!   Code for x < =  4.0
!
  if ( x <= four ) then
     if ( x < xlow2 ) then
        tran08 = zero
     else
        if ( x < xlow1 ) then
           tran08 = ( x ** ( numjn - 1 ) ) / ( rnumjn - one )
        else
           t = ( ( ( x * x ) / eight ) - half ) - half
           tran08 = ( x ** ( numjn - 1 ) ) * cheval ( nterms, atran, t )
        end if
     end if
  else
!
!  Code for x > 4.0
!
     if ( xhigh2 < x ) then
        sumexp = one
     else
        if ( x <= xhigh1 ) then
           numexp = int ( xhigh1 / x ) + 1
           t = exp ( - x )
        else
           numexp = 1
           t = one
        end if
        rk = zero
        do k1 = 1, numexp
           rk = rk + one
        end do
        sumexp = zero
        do k1 = 1, numexp
           sum2 = one
           xk = one / ( rk * x )
           xk1 = one
           do k2 = 1, numjn
              sum2 = sum2 * xk1 * xk + one
              xk1 = xk1 + one
           end do
           sumexp = sumexp * t + sum2
           rk = rk - one
        end do
     end if
     t = rnumjn * log ( x ) - x + log ( sumexp )
     if ( t < xhigh3 ) then
        tran08 = valinf
     else
        tran08 = valinf - exp ( t )
     end if
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
!      TRAN08(x) = Integral ( 0 <= t <= x ) t^8 * exp(t) / ( exp(t) - 1 )^2 dt
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tran09 ( xvalue )

!*****************************************************************************80
!
!! TRAN09 calculates the transport integral of order 9.
!
!  Discussion:
!
!    The function is defined by:
!
!      TRAN09(x) = Integral ( 0 <= t <= x ) t^9 * exp(t) / ( exp(t) - 1 )^2 dt
!
!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) TRAN09, the value of the function.
!
  implicit none

  real ( kind = 8 ) atran(0:19)
  real ( kind = 8 ) cheval
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer k1
  integer k2
  integer, parameter :: nterms = 17
  integer numexp
  integer, parameter :: numjn = 9
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) rk
  real ( kind = 8 ), parameter :: rnumjn = 9.0D+00
  real ( kind = 8 ) sumexp
  real ( kind = 8 ) sum2
  real ( kind = 8 ) t
  real ( kind = 8 ) tran09
  real ( kind = 8 ), parameter :: valinf = 0.36360880558872871397d6
  real ( kind = 8 ) x
  real ( kind = 8 ) xhigh1
  real ( kind = 8 ) xhigh2
  real ( kind = 8 ) xhigh3
  real ( kind = 8 ) xk
  real ( kind = 8 ) xk1
  real ( kind = 8 ) xlow1
  real ( kind = 8 ) xlow2
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  data atran/0.16224049991949846835d0, &
            -0.3768351452195937773d-1, &
             0.547669715917719770d-2, &
            -0.64443945009449521d-3, &
             0.6773645285280983d-4, &
            -0.666813497582042d-5, &
             0.63047560019047d-6, &
            -0.5807478663611d-7, &
             0.525551305123d-8, &
            -0.46968861761d-9, &
             0.4159395065d-10, &
            -0.365808491d-11, &
             0.32000794d-12, &
            -0.2787651d-13, &
             0.242017d-14, &
            -0.20953d-15, &
             0.1810d-16, &
            -0.156d-17, &
             0.13d-18, &
            -0.1d-19/
!
!  Machine-dependent constants (for IEEE machines)
!
  data xlow1,xlow2/2.98023224d-8,4.5321503d-39/
  data xhigh1,xhigh3/36.04365668d0,-36.73680056d0/
  data xhigh2/4.05323966d16/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAN09 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    tran09 = zero
    return
  end if
!
!   Code for x < =  4.0
!
  if ( x <= four ) then
     if ( x < xlow2 ) then
        tran09 = zero
     else
        if ( x < xlow1 ) then
           tran09 = ( x ** ( numjn - 1 ) ) / ( rnumjn - one )
        else
           t = ( ( ( x * x ) / eight ) - half ) - half
           tran09 = ( x ** ( numjn - 1 ) ) * cheval ( nterms, atran, t )
        end if
     end if
  else
!
!  Code for x > 4.0
!
     if ( xhigh2 < x ) then
        sumexp = one
     else
        if ( x <= xhigh1 ) then
           numexp = int ( xhigh1 / x ) + 1
           t = exp ( -x )
        else
           numexp = 1
           t = one
        end if
        rk = zero
        do k1 = 1, numexp
           rk = rk + one
        end do
        sumexp = zero
        do k1 = 1, numexp
           sum2 = one
           xk = one / ( rk * x )
           xk1 = one
           do k2 = 1, numjn
              sum2 = sum2 * xk1 * xk + one
              xk1 = xk1 + one
           end do
           sumexp = sumexp * t + sum2
           rk = rk - one
        end do
     end if

     t = rnumjn * log ( x ) - x + log ( sumexp )

     if ( t < xhigh3 ) then
        tran09 = valinf
     else
        tran09 = valinf - exp ( t )
     end if

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
!      TRAN09(x) = Integral ( 0 <= t <= x ) t^9 * exp(t) / ( exp(t) - 1 )^2 dt
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
!    US Department of Commerce, 1964.
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
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
  integer n_data
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
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
