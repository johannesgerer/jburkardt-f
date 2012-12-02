subroutine testdt ( k, n, x, y )

!*****************************************************************************80
!
!! TESTDT returns one of five sets of test nodes.
!
!  Discussion:
!
!    This routine returns one of five sets of nodes used
!    for testing scattered data fitting methods.  All five sets
!    approximately cover the unit square [0,1] X [0,1]:  the
!    convex hulls of sets 1 and 3 extend slightly outside the
!    square but do not completely cover it, those of sets 2 and
!    5 coincide with the unit square, and the convex hull of
!    set 4 is a large subset of the unit square.
!
!  Modified:
!
!    26 January 2012
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Richard Franke,
!    A Critical Comparison of Some Methods for Interpolation of Scattered Data,
!    Naval Postgraduate School Technical Report,
!    NPS-53-79-003, 1979.
!
!    Robert Renka, Ron Brown,
!    Algorithm 792: Accuracy tests for ACM algorithms for interpolation of
!    scattered data in the plane,
!    ACM Transactions on Mathematical Software,
!    Volume 25, Number 1, March 1999, pages 78-94.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, chooses the data set:
!    1, Franke's 100 node set;
!    2, Franke's 33 node set;
!    3, Lawson's 25 node set;
!    4, Random 100 node set;
!    5, Gridded 81 node set.
!
!    Output, integer ( kind = 4 ) N, the size of the data vectors.
!
!    Output, real ( kind = 8 ), X(N), Y(N), the coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) x(*)
  real ( kind = 8 ), dimension ( 100 ) :: x1 = (/ &
    0.0227035,  0.0539888,  0.0217008,  0.0175129,  0.0019029, &
   -0.0509685,  0.0395408, -0.0487061,  0.0315828, -0.0418785, &
    0.1324189,  0.1090271,  0.1254439,  0.0934540,  0.0767578, &
    0.1451874,  0.0626494,  0.1452734,  0.0958668,  0.0695559, &
    0.2645602,  0.2391645,  0.2088990,  0.2767329,  0.1714726, &
    0.2266781,  0.1909212,  0.1867647,  0.2304634,  0.2426219, &
    0.3663168,  0.3857662,  0.3832392,  0.3179087,  0.3466321, &
    0.3776591,  0.3873159,  0.3812917,  0.3795364,  0.2803515, &
    0.4149771,  0.4277679,  0.4200010,  0.4663631,  0.4855658, &
    0.4092026,  0.4792578,  0.4812279,  0.3977761,  0.4027321, &
    0.5848691,  0.5730076,  0.6063893,  0.5013894,  0.5741311, &
    0.6106955,  0.5990105,  0.5380621,  0.6096967,  0.5026188, &
    0.6616928,  0.6427836,  0.6396475,  0.6703963,  0.7001181, &
    0.6333590,  0.6908947,  0.6895638,  0.6718889,  0.6837675, &
    0.7736939,  0.7635332,  0.7410424,  0.8258981,  0.7306034, &
    0.8086609,  0.8214531,  0.7290640,  0.8076643,  0.8170951, &
    0.8424572,  0.8684053,  0.8366923,  0.9418461,  0.8478122, &
    0.8599583,  0.9175700,  0.8596328,  0.9279871,  0.8512805, &
    1.0449820,  0.9670631,  0.9857884,  0.9676313,  1.0129299, &
    0.9657040,  1.0019855,  1.0359297,  1.0414677,  0.9471506 /)
  real ( kind = 8 ), dimension ( 33 ) :: x2 = (/ &
    0.05,  0.00,  0.00,  0.00,  0.10, &
    0.10,  0.15,  0.20,  0.25,  0.30, &
    0.35,  0.50,  0.50,  0.55,  0.60, &
    0.60,  0.60,  0.65,  0.70,  0.70, &
    0.70,  0.75,  0.75,  0.75,  0.80, &
    0.80,  0.85,  0.90,  0.90,  0.95, &
    1.00,  1.00,  1.00 /)
  real ( kind = 8 ), dimension ( 25 ) :: x3 = (/ &
    0.13750,   0.91250,   0.71250,   0.22500,  -0.05000, &
    0.47500,   0.05000,   0.45000,   1.08750,   0.53750, &
   -0.03750,   0.18750,   0.71250,   0.85000,   0.70000, &
    0.27500,   0.45000,   0.81250,   0.45000,   1.00000, &
    0.50000,   0.18750,   0.58750,   1.05000,   0.10000 /)
  real ( kind = 8 ), dimension ( 100 ) :: x4 = (/ &
    0.0096326,  0.0216348,  0.0298360,  0.0417447,  0.0470462, &
    0.0562965,  0.0646857,  0.0740377,  0.0873907,  0.0934832, &
    0.1032216,  0.1110176,  0.1181193,  0.1251704,  0.1327330, &
    0.1439536,  0.1564861,  0.1651043,  0.1786039,  0.1886405, &
    0.2016706,  0.2099886,  0.2147003,  0.2204141,  0.2343715, &
    0.2409660,  0.2527740,  0.2570839,  0.2733365,  0.2853833, &
    0.2901755,  0.2964854,  0.3019725,  0.3125695,  0.3307163, &
    0.3378504,  0.3439061,  0.3529922,  0.3635507,  0.3766172, &
    0.3822429,  0.3869838,  0.3973137,  0.4170708,  0.4255588, &
    0.4299218,  0.4372839,  0.4705033,  0.4736655,  0.4879299, &
    0.4940260,  0.5055324,  0.5162593,  0.5219219,  0.5348529, &
    0.5483213,  0.5569571,  0.5638611,  0.5784908,  0.5863950, &
    0.5929148,  0.5987839,  0.6117561,  0.6252296,  0.6331381, &
    0.6399048,  0.6488972,  0.6558537,  0.6677405,  0.6814074, &
    0.6887812,  0.6940896,  0.7061687,  0.7160957,  0.7317445, &
    0.7370798,  0.7462030,  0.7566957,  0.7699998,  0.7879347, &
    0.7944014,  0.8164468,  0.8192794,  0.8368405,  0.8500993, &
    0.8588255,  0.8646496,  0.8792329,  0.8837536,  0.8900077, &
    0.8969894,  0.9044917,  0.9083947,  0.9203972,  0.9347906, &
    0.9434519,  0.9490328,  0.9569571,  0.9772067,  0.9983493 /)
  real ( kind = 8 ), dimension ( 81 ) :: x5 = (/ &
    0.125,  0.000,  0.000,  0.000,  0.000, &
    0.000,  0.000,  0.000,  0.000,  0.000, &
    0.125,  0.125,  0.125,  0.125,  0.125, &
    0.125,  0.125,  0.125,  0.250,  0.250, &
    0.250,  0.250,  0.250,  0.250,  0.250, &
    0.250,  0.250,  0.375,  0.375,  0.375, &
    0.375,  0.375,  0.375,  0.375,  0.375, &
    0.375,  0.500,  0.500,  0.500,  0.500, &
    0.500,  0.500,  0.500,  0.500,  0.500, &
    0.625,  0.625,  0.625,  0.625,  0.625, &
    0.625,  0.625,  0.625,  0.625,  0.750, &
    0.750,  0.750,  0.750,  0.750,  0.750, &
    0.750,  0.750,  0.750,  0.875,  0.875, &
    0.875,  0.875,  0.875,  0.875,  0.875, &
    0.875,  0.875,  1.000,  1.000,  1.000, &
    1.000,  1.000,  1.000,  1.000,  1.000, &
    1.000 /)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ), dimension ( 100 ) :: y1 = (/ &
  -0.0310206,   0.1586742,   0.2576924,   0.3414014,   0.4943596, &
   0.5782854,   0.6993418,   0.7470194,   0.9107649,   0.9962890, &
   0.0501330,   0.0918555,   0.2592973,   0.3381592,   0.4171125, &
   0.5615563,   0.6552235,   0.7524066,   0.9146523,   0.9632421, &
   0.0292939,   0.0602303,   0.2668783,   0.3696044,   0.4801738, &
   0.5940595,   0.6878797,   0.8185576,   0.9046507,   0.9805412, &
   0.0396955,   0.0684484,   0.2389548,   0.3124129,   0.4902989, &
   0.5199303,   0.6445227,   0.8203789,   0.8938079,   0.9711719, &
  -0.0284618,   0.1560965,   0.2262471,   0.3175094,   0.3891417, &
   0.5084949,   0.6324247,   0.7511007,   0.8489712,   0.9978728, &
  -0.0271948,   0.1272430,   0.2709269,   0.3477728,   0.4259422, &
   0.6084711,   0.6733781,   0.7235242,   0.9242411,   1.0308762, &
   0.0255959,   0.0707835,   0.2008336,   0.3259843,   0.4890704, &
   0.5096324,   0.6697880,   0.7759569,   0.9366096,   1.0064516, &
   0.0285374,   0.1021403,   0.1936581,   0.3235775,   0.4714228, &
   0.6091595,   0.6685053,   0.8022808,   0.8476790,   1.0512371, &
   0.0380499,   0.0902048,   0.2083092,   0.3318491,   0.4335632, &
   0.5910139,   0.6307383,   0.8144841,   0.9042310,   0.9696030, &
  -0.0120900,   0.1334114,   0.2695844,   0.3795281,   0.4396054, &
   0.5044425,   0.6941519,   0.7459923,   0.8682081,   0.9801409 /)
  real ( kind = 8 ), dimension ( 33 ) :: y2 = (/ &
    0.45,  0.50,  1.00,  0.00,  0.15, &
    0.75,  0.30,  0.10,  0.20,  0.35, &
    0.85,  0.00,  1.00,  0.95,  0.25, &
    0.65,  0.85,  0.70,  0.20,  0.65, &
    0.90,  0.10,  0.35,  0.85,  0.40, &
    0.65,  0.25,  0.35,  0.80,  0.90, &
    0.00,  0.50,  1.00 /)
  real ( kind = 8 ), dimension ( 25 ) :: y3 = (/ &
    0.97500,   0.98750,   0.76250,   0.83750,   0.41250, &
    0.63750,  -0.05000,   1.03750,   0.55000,   0.80000, &
    0.75000,   0.57500,   0.55000,   0.43750,   0.31250, &
    0.42500,   0.28750,   0.18750,  -0.03750,   0.26250, &
    0.46250,   0.26250,   0.12500,  -0.06125,   0.11250 /)
  real ( kind = 8 ), dimension ( 100 ) :: y4 = (/ &
    0.3083158,  0.2450434,  0.8613847,  0.0977864,  0.3648355, &
    0.7156339,  0.5311312,  0.9755672,  0.1781117,  0.5452797, &
    0.1603881,  0.7837139,  0.9982015,  0.6910589,  0.1049580, &
    0.8184662,  0.7086405,  0.4456593,  0.1178342,  0.3189021, &
    0.9668446,  0.7571834,  0.2016598,  0.3232444,  0.4368583, &
    0.8907869,  0.0647260,  0.5692618,  0.2947027,  0.4332426, &
    0.3347464,  0.7436284,  0.1066265,  0.8845357,  0.5158730, &
    0.9425637,  0.4799701,  0.1783069,  0.1146760,  0.8225797, &
    0.2270688,  0.4073598,  0.8875080,  0.7631616,  0.9972804, &
    0.4959884,  0.3410421,  0.2498120,  0.6409007,  0.1058690, &
    0.5411969,  0.0089792,  0.8784268,  0.5515874,  0.4038952, &
    0.1654023,  0.2965158,  0.3660356,  0.0366554,  0.9502420, &
    0.2638101,  0.9277386,  0.5377694,  0.7374676,  0.4674627, &
    0.9186109,  0.0416884,  0.1291029,  0.6763676,  0.8444238, &
    0.3273328,  0.1893879,  0.0645923,  0.0180147,  0.8904992, &
    0.4160648,  0.4688995,  0.2174508,  0.5734231,  0.8853319, &
    0.8018436,  0.6388941,  0.8931002,  0.1000558,  0.2789506, &
    0.9082948,  0.3259159,  0.8318747,  0.0508513,  0.9708450, &
    0.5120548,  0.2859716,  0.9581641,  0.6183429,  0.3779934, &
    0.4010423,  0.9478657,  0.7425486,  0.8883287,  0.5496750 /)
  real ( kind = 8 ), dimension ( 81 ) :: y5 = (/ &
    0.000,  0.125,  0.250,  0.375,  0.500, &
    0.625,  0.750,  0.875,  1.000,  0.000, &
    0.125,  0.250,  0.375,  0.500,  0.625, &
    0.750,  0.875,  1.000,  0.000,  0.125, &
    0.250,  0.375,  0.500,  0.625,  0.750, &
    0.875,  1.000,  0.000,  0.125,  0.250, &
    0.375,  0.500,  0.625,  0.750,  0.875, &
    1.000,  0.000,  0.125,  0.250,  0.375, &
    0.500,  0.625,  0.750,  0.875,  1.000, &
    0.000,  0.125,  0.250,  0.375,  0.500, &
    0.625,  0.750,  0.875,  1.000,  0.000, &
    0.125,  0.250,  0.375,  0.500,  0.625, &
    0.750,  0.875,  1.000,  0.000,  0.125, &
    0.250,  0.375,  0.500,  0.625,  0.750, &
    0.875,  1.000,  0.000,  0.125,  0.250, &
    0.375,  0.500,  0.625,  0.750,  0.875, &
    1.000 /)

  if ( k == 1 ) then
    n = 100
    x(1:n) = x1(1:n)
    y(1:n) = y1(1:n)
  else if ( k == 2 ) then
    n = 33
    x(1:n) = x2(1:n)
    y(1:n) = y2(1:n)
  else if ( k == 3 ) then
    n = 25
    x(1:n) = x3(1:n)
    y(1:n) = y3(1:n)
  else if ( k == 4 ) then
    n = 100
    x(1:n) = x4(1:n)
    y(1:n) = y4(1:n)
  else if ( k == 5 ) then
    n = 81 
    x(1:n) = x5(1:n)
    y(1:n) = y5(1:n)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTDT - Fatal error!'
    write ( *, '(a)' ) '  Input K does not satisfy 1 <= K <= 5.'
    stop
  end if

  return
end
subroutine tstfn1 ( k, x, y, flag, f, fx, fy )

!*********************************************************************72
!
!! TSTFN1 computes one of ten bivariate test functions.
!
!  Discussion:
!
!    This subroutine computes the value and, optionally, the
!    first partial derivatives of one of ten bivariate test
!    functions.  The first six functions were chosen by Richard
!    Franke to test interpolation software (See the reference
!    below).  The last four functions represent more challenging 
!    surface fitting problems.
!
!    Note that function 6 is only defined inside a circle of
!    radius 8/9 centered at (.5,.5).  Thus, if (X-.5)**2 +
!    (Y-.5)**2 .GE. 64/81, the value (and partials if FLAG=1)
!    are set to 0 for this function.  Also, the first partial
!    derivatives of function 10 are not defined at (.5,.5) --
!    again, zeros are returned.
!
!  Modified:
!
!    26 January 2012
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Richard Franke,
!    A Critical Comparison of Some Methods for Interpolation of Scattered Data,
!    Naval Postgraduate School Technical Report,
!    NPS-53-79-003, 1979.
!
!    Robert Renka, Ron Brown,
!    Algorithm 792: Accuracy tests for ACM algorithms for interpolation of
!    scattered data in the plane,
!    ACM Transactions on Mathematical Software,
!    Volume 25, Number 1, March 1999, pages 78-94.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, chooses the function:
!     1, Exponential;
!     2, Cliff;
!     3, Saddle;
!     4, Gentle;
!     5, Steep;
!     6, Sphere;
!     7, Trig;
!     8, Synergistic Gaussian;
!     9, Cloverleaf Asymmetric Peak/Valley;
!    10, Cosine Peak.
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point where the
!    function is to be evaluated.
!
!    Input, integer ( kind = 4 ) FLAG, indicates whether derivatives are desired.
!    0, only a function value is required;
!    1, both the function and its first partial derivatives are required.
!
!    Output, real ( kind = 8 ) F, the value of the function at (X,Y).
!
!    Output, real ( kind = 8 ) FX, FY, the first partial derivatives of the
!    function at (X,Y) if 1 <= FLAG.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) flag
  real ( kind = 8 ) fx
  real ( kind = 8 ) fy
  integer ( kind = 4 ) k
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Exponential:
!
  if ( k == 1 ) then

    f = 0.75 * exp ( - ( ( 9.0 * x - 2.0 )**2 + ( 9.0 * y - 2.0 )**2 ) / 4.0 ) + &
        0.75 * exp ( - ( ( 9.0 * x + 1.0 )**2 ) / 49.0 - ( 9.0 * y + 1.0 ) / 10.0 ) + &
         0.5 * exp ( - ( ( 9.0 * x - 7.0 )**2 + ( 9.0 * y - 3.0 )**2 ) / 4.0 ) - &
         0.2 * exp ( - ( 9.0 * x - 4.0 )**2 - ( 9.0 * y - 7.0 )**2 )

    if ( 1 <= flag ) then

      t1 = exp ( - ( ( 9.0 * x - 2.0 )**2 + ( 9.0 * y - 2.0 )**2 ) / 4.0 )
      t2 = exp ( - ( ( 9.0 * x + 1.0 )**2 ) / 49.0 - ( 9.0 * y + 1.0 ) / 10.0 )
      t3 = exp ( - ( ( 9.0 * x - 7.0 )**2 + ( 9.0 * y - 3.0 )**2 ) / 4.0 )
      t4 = exp ( - ( 9.0 * x - 4.0 )**2 - ( 9.0 * y - 7.0 )**2 )

      fx = - 3.375 * ( 9.0 * x - 2.0 ) * t1 - ( 27.0 / 98.0 ) * ( 9.0 * x + 1.0 ) * t2 &
        - 2.25 * ( 9.0 * x - 7.0 ) * t3 + 3.6 * ( 9.0 * x - 4.0 ) * t4

      fy = - 3.375 * ( 9.0 * y - 2.0 ) * t1 - 0.675 * t2 &
        - 2.25 * ( 9.0 * y - 3.0 ) * t3 + 3.6 * ( 9.0 * y - 7.0 ) * t4

    end if
!
!  Cliff:
!
  else if ( k == 2 ) then

    f = ( tanh ( 9.0 * ( y - x ) ) + 1.0 ) / 9.0

    if ( 1 <= flag ) then
      t1 = 18.0 * ( y - x )
      fx = - 4.0 / ( exp ( t1 ) + 2.0 + exp ( - t1 ) )
      fy = - fx
    end if
!
!  Saddle:
!
  else if ( k == 3 ) then

    f = ( 1.25 + cos ( 5.4 * y ) ) / ( 6.0 + 6.0 * ( 3.0 * x - 1.0 )**2 )

    if ( 1 <= flag ) then
      t1 = 5.4 * y
      t2 = 1.0 + ( 3.0 * x - 1.0 )**2
      fx = - ( 3.0 * x - 1.0 ) * ( 1.25 + cos ( t1 ) ) / ( t2 * t2 )
      fy = - 0.9 * sin ( t1 ) / t2
    end if
!
!  Gentle:
!
  else if ( k == 4 ) then

    f = exp ( - 5.0625 * ( ( x - 0.5 )**2 + ( y - 0.5 )**2 ) ) / 3.0

    if ( 1 <= flag ) then
      t1 = x - 0.5
      t2 = y - 0.5
      t3 = - 3.375 * exp ( - 5.0625 * ( t1 * t1 + t2 * t2 ) )
      fx = t1 * t3
      fy = t2 * t3
    end if
!
!  Steep:
!
  else if ( k == 5 ) then

    f = exp ( - 20.25 * ( ( x - 0.5 )**2 + ( y - 0.5 )**2 ) ) / 3.0

    if ( 1 <= flag ) then
      t1 = x - 0.5
      t2 = y - 0.5
      t3 = - 13.5 * exp ( - 20.25 * ( t1 * t1 + t2 * t2 ) )
      fx = t1 * t3
      fy = t2 * t3
    end if
!
!  Sphere:
!
  else if ( k == 6 ) then

    t4 = 64.0 - 81.0 * ( ( x - 0.5 )**2 + ( y - 0.5 )**2 )

    if ( 0.0 <= t4 ) then
      f = sqrt ( t4 ) / 9.0 - 0.5
    else
      f = 0.0
    end if

    if ( 1 <= flag ) then
      t1 = x - 0.5
      t2 = y - 0.5
      if ( 0.0 < t4 ) then
        t3 = - 9.0 / sqrt ( t4 )
        fx = t1 * t3
        fy = t2 * t3
      else
        fx = 0.0
        fy = 0.0
      end if

    end if
!
!  Trig:
!
  else if ( k == 7 ) then

    f = 2.0 * cos ( 10.0 * x ) * sin ( 10.0 * y ) + sin ( 10.0 * x * y )

    if ( 1 <= flag ) then
      t1 = 10.0 * x
      t2 = 10.0 * y
      t3 = 10.0 * cos ( 10.0 * x * y )
      fx = - 20.0 * sin ( t1 ) * sin ( t2 ) + t3 * y
      fy = 20.0 * cos ( t1 ) * cos ( t2 ) + t3 * x
    end if
!
!  gaussx(1,.5,.1) + gaussy(.75,.5,.1) + gaussx(1,.5,.1)*
!  gaussy(.75,.5,.1), where gaussx(a,b,c) is the gaussian
!  function of x with amplitude a, center (mean) b, and
!  width (standard deviation) c.
!
  else if ( k == 8 ) then

    t1 = 5.0 - 10.0 * x
    t2 = 5.0 - 10.0 * y
    t3 = exp ( - 0.5 * t1 * t1 )
    t4 = exp ( - 0.5 * t2 * t2 )
    f = t3 + 0.75 * t4 * ( 1.0 + t3 )

    if ( 1 <= flag ) then
      fx = t1 * t3 * ( 10.0 + 7.5 * t4 )
      fy = t2 * t4 * ( 7.5 + 7.5 * t3 )
    end if
!
!  Cloverleaf asymmetric hill/valley:
!
  else if ( k == 9 ) then

    t1 = exp ( ( 10.0 - 20.0 * x ) / 3.0 )
    t2 = exp ( ( 10.0 - 20.0 * y ) / 3.0 )
    t3 = 1.0 / ( 1.0 + t1 )
    t4 = 1.0 / ( 1.0 + t2 )
    f = ( ( 20.0 / 3.0 )**3 * t1 * t2 )**2 * ( t3 * t4 )**5 * &
      ( t1 - 2.0 * t3 ) * ( t2 - 2.0 * t4 )

    if ( 1 <= flag ) then

      fx = ( ( 20.0 / 3.0 ) * t1 )**2 * ( ( 20.0 / 3.0 ) * t3 )**5 * &
        ( 2.0 * t1 - 3.0 * t3 - 5.0 + 12.0 * t3 * t3 ) * t2 * t2 * t4**5 * &
        ( t2 - 2.0 * t4 )

      fy = ( ( 20.0 / 3.0 ) * t1 )**2 * ( ( 20.0 / 3.0 ) * t3 )**5 * &
        ( 2.0 * t2 - 3.0 * t4 - 5.0 + 12.0 * t4 * t4 ) * t2 * t2 * t4**5 * &
        ( t1 - 2.0 * t3 )

    end if
!
!  Cosine peak:
!
  else if ( k == 10 ) then

    t1 = sqrt ( ( 80.0 * x - 40.0 )**2 + ( 90.0 * y - 45.0 )**2 )
    t2 = exp ( - 0.04 * t1 )
    t3 = cos ( 0.15 * t1 )
    f = t2 * t3

    if ( 1 <= flag ) then

      if ( t1 == 0.0 ) then
        fx = 0.0
        fy = 0.0
      else
        t4 = sin ( 0.15 * t1 )
        fx = - t2 * ( 12.0 * t4 + 3.2 * t3 ) * ( 80.0 * x - 40.0 ) / t1
        fy = - t2 * ( 13.5 * t4 + 3.6 * t3 ) * ( 90.0 * y - 45.0 ) / t1
      end if

    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TSTFN1 - Fatal error!'
    write ( *, '(a)' ) '  Input K does not satisfy 1 <= K <= 10.'
    stop

  end if

  return
  end
  subroutine tstfn2 ( k, x, y, flag, f, fx, fy, fxx, fxy, fyy )

!*********************************************************************72
!
!! TSTFN2 evaluates one of ten bivariate test functions and derivatives.
!
!  Discussion:
!
!    This subroutine computes the value and, optionally, the
!    first and/or second partial derivatives of one of ten
!    bivariate test functions.  The first six functions were
!    chosen by Richard Franke to test interpolation software.
!    The last four functions represent more challenging surface fitting problems.
!
!    Note that function 6 is only defined inside a circle of
!    radius 8/9 centered at (.5,.5).  Thus, if (X-.5)**2 +
!    (Y-.5)**2 .GE. 64/81, the value (and partials if IFLAG=1)
!    are set to 0 for this function.  Also, the first partial
!    derivatives of function 10 are not defined at (.5,.5) --
!    again, zeros are returned.
!
!  Modified:
!
!    26 January 2012
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Richard Franke,
!    A Critical Comparison of Some Methods for Interpolation of Scattered Data,
!    Naval Postgraduate School Technical Report,
!    NPS-53-79-003, 1979.
!
!    Robert Renka, Ron Brown,
!    Algorithm 792: Accuracy tests for ACM algorithms for interpolation of
!    scattered data in the plane,
!    ACM Transactions on Mathematical Software,
!    Volume 25, Number 1, March 1999, pages 78-94.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, chooses the function:
!     1, Exponential;
!     2, Cliff;
!     3, Saddle;
!     4, Gentle;
!     5, Steep;
!     6, Sphere;
!     7, Trig;
!     8, Synergistic Gaussian;
!     9, Cloverleaf Asymmetric Peak/Valley;
!    10, Cosine Peak.
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point where the
!    function is to be evaluated.
!
!    Input, integer ( kind = 4 ) FLAG, indicates whether derivatives are desired.
!    0, only a function value is required;
!    1, the function and its first partial derivatives are required.
!    2, the function, first and second partial derivatives are required.
!
!    Output, real ( kind = 8 ) F, the value of the function at (X,Y).
!
!    Output, real ( kind = 8 ) FX, FY, the first partial derivatives of the
!    function at (X,Y) if 1 <= FLAG.
!
!    Output, real ( kind = 8 ) FXX, FXY, FYY, the second partial derivatives
!    of the function at (X,Y) if 2 <= FLAG.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) flag
  real ( kind = 8 ) fx
  real ( kind = 8 ) fxx
  real ( kind = 8 ) fxy
  real ( kind = 8 ) fy
  real ( kind = 8 ) fyy
  integer ( kind = 4 ) k
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) t5
  real ( kind = 8 ) t6
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Exponential:
!
  if ( k == 1 ) then

    f = 0.75 * exp ( - ( ( 9.0 * x - 2.0 )**2 + ( 9.0 * y - 2.0 )**2 ) / 4.0 ) + &
        0.75 * exp ( - ( ( 9.0 * x + 1.0 )**2 ) / 49.0 - ( 9.0 * y + 1.0 ) / 10.0 ) + &
         0.5 * exp ( - ( ( 9.0 * x - 7.0 )**2 + ( 9.0 * y - 3.0 )**2 ) / 4.0 ) - &
         0.2 * exp ( - ( 9.0 * x - 4.0 )**2 - ( 9.0 * y - 7.0 )**2 )

    if ( 1 <= flag ) then

      t1 = exp ( - ( ( 9.0 * x - 2.0 )**2 + ( 9.0 * y - 2.0 )**2 ) / 4.0 )
      t2 = exp ( - ( ( 9.0 * x + 1.0 )**2 ) / 49.0 - ( 9.0 * y + 1.0 ) / 10.0 )
      t3 = exp ( - ( ( 9.0 * x - 7.0 )**2 + ( 9.0 * y - 3.0 )**2 ) / 4.0 )
      t4 = exp ( - ( 9.0 * x - 4.0 )**2 - ( 9.0 * y - 7.0 )**2 )

      fx = - 3.375 * ( 9.0 * x - 2.0 ) * t1 - ( 27.0 / 98.0 ) * ( 9.0 * x + 1.0 ) * t2 &
        - 2.25 * ( 9.0 * x - 7.0 ) * t3 + 3.6 * ( 9.0 * x - 4.0 ) * t4

      fy = - 3.375 * ( 9.0 * y - 2.0 ) * t1 - 0.675 * t2 &
        - 2.25 * ( 9.0 * y - 3.0 ) * t3 + 3.6 * ( 9.0 * y - 7.0 ) * t4

    end if

    if ( 2 <= flag ) then

      fxx = 15.1875 * ( ( 9.0 * x - 2.0 )**2 - 2.0 ) * t1 + &
        60.75 * ( ( 9.0 * x + 1.0 )**2 - 24.5 ) * t2 + &
        10.125 * ( ( 9.0 * x - 7.0 )**2 - 2.0 ) * t3 - &
        64.8 * ( ( 9.0 * x - 4.0 )**2 - 0.5 ) * t4

      fxy = 15.1875 * ( 9.0 * x - 2.0 ) * ( 9.0 * y - 2.0 ) * t1 + &
        ( 243.0 / 980.0 ) * ( 9.0 * x + 1.0 ) * t2 + &
        10.125 * ( 9.0 * x - 7.0 ) * ( 9.0 * y - 3.0 ) * t3 - &
        64.8 * ( 9.0 * x - 4.0 ) * ( 9.0 * y - 7.0 ) * t4

      fyy = 15.1875 * ( ( 9.0 * y - 2.0 )**2 - 2.0 ) * t1 + &
        0.6075 * t2 + &
        10.125 * ( ( 9.0 * y - 3.0 )**2 - 2.0 ) * t3 - &
        64.8 * ( ( 9.0 * y - 7.0 )**2 - 0.5 ) * t4

    end if
!
!  Cliff:
!
  else if ( k == 2 ) then

    f = ( tanh ( 9.0 * ( y - x ) ) + 1.0 ) / 9.0

    if ( 1 <= flag ) then
      t1 = 18.0 * ( y - x )
      fx = - 4.0 / ( exp ( t1 ) + 2.0 + exp ( - t1 ) )
      fy = - fx
    end if

    if ( 2 <= flag ) then
      fxx = 18.0 * tanh ( 0.5 * t1 ) * fx
      fxy = - fxx
      fyy = fxx
    end if
!
!  Saddle:
!
  else if ( k == 3 ) then

    f = ( 1.25 + cos ( 5.4 * y ) ) / ( 6.0 + 6.0 * ( 3.0 * x - 1.0 )**2 )

    if ( 1 <= flag ) then
      t1 = 5.4 * y
      t2 = 1.0 + ( 3.0 * x - 1.0 )**2
      fx = - ( 3.0 * x - 1.0 ) * ( 1.25 + cos ( t1 ) ) / ( t2 * t2 )
      fy = - 0.9 * sin ( t1 ) / t2
    end if

    if ( 2 <= flag ) then
      fxx = 3.0 * ( 1.25 + cos ( t1 ) ) * ( 3.0 * t2 - 4.0 ) / ( t2**3 )
      fxy = 5.4 * ( 3.0 * x - 1.0 ) * sin ( t1 ) / ( t2 * t2 )
      fyy = - 4.86 * cos ( t1 ) / t2
    end if
!
!  Gentle:
!
  else if ( k == 4 ) then

    f = exp ( - 5.0625 * ( ( x - 0.5 )**2 + ( y - 0.5 )**2 ) ) / 3.0

    if ( 1 <= flag ) then
      t1 = x - 0.5
      t2 = y - 0.5
      t3 = - 3.375 * exp ( - 5.0625 * ( t1 * t1 + t2 * t2 ) )
      fx = t1 * t3
      fy = t2 * t3
    end if

    if ( 2 <= flag ) then
      fxx = ( 1.0 - 10.125 * t1 * t1 ) * t3
      fxy = - 10.125 * t1 * t2 * t3
      fyy = ( 1.0 - 10.125 * t2 * t2 ) * t3
    end if
!
!  Steep:
!
  else if ( k == 5 ) then

    f = exp ( - 20.25 * ( ( x - 0.5 )**2 + ( y - 0.5 )**2 ) ) / 3.0

    if ( 1 <= flag ) then
      t1 = x - 0.5
      t2 = y - 0.5
      t3 = - 13.5 * exp ( - 20.25 * ( t1 * t1 + t2 * t2 ) )
      fx = t1 * t3
      fy = t2 * t3
    end if

    if ( 2 <= flag ) then
      fxx = ( 1.0 - 40.5 * t1 * t1 ) * t3
      fxy = - 40.5 * t1 * t2 * t3
      fyy = ( 1.0 - 40.5 * t2 * t2 ) * t3
    end if
!
!  Sphere:
!
  else if ( k == 6 ) then

    t4 = 64.0 - 81.0 * ( ( x - 0.5 )**2 + ( y - 0.5 )**2 )

    if ( 0.0 <= t4 ) then
      f = sqrt ( t4 ) / 9.0 - 0.5
    else
      f = 0.0
    end if

    if ( 1 <= flag ) then
      t1 = x - 0.5
      t2 = y - 0.5
      if ( 0.0 < t4 ) then
        t3 = - 9.0 / sqrt ( t4 )
        fx = t1 * t3
        fy = t2 * t3
      else
        fx = 0.0
        fy = 0.0
      end if

    end if

    if ( 2 <= flag ) then
      fxx = ( 1.0 + fx * fx ) * t3
      fxy = fx * fy
      fyy = ( 1.0 + fy * fy ) * t3
    end if
!
!  Trig:
!
  else if ( k == 7 ) then

    f = 2.0 * cos ( 10.0 * x ) * sin ( 10.0 * y ) + sin ( 10.0 * x * y )

    if ( 1 <= flag ) then
      t1 = 10.0 * x
      t2 = 10.0 * y
      t3 = 10.0 * cos ( 10.0 * x * y )
      fx = - 20.0 * sin ( t1 ) * sin ( t2 ) + t3 * y
      fy = 20.0 * cos ( t1 ) * cos ( t2 ) + t3 * x
    end if

    if ( 2 <= flag ) then
      t4 = 100.0 * sin ( 10.0 * x * y )
      fxx = - 200.0 * cos ( t1 ) * sin ( t2 ) - t4 * y * y
      fxy = - 200.0 * sin ( t1 ) * cos ( t2 ) + t3 - t4 * x * y
      fyy = - 200.0 * cos ( t1 ) * sin ( t2 ) - t4 * x * x
    end if
!
!  gaussx(1,.5,.1) + gaussy(.75,.5,.1) + gaussx(1,.5,.1)*
!  gaussy(.75,.5,.1), where gaussx(a,b,c) is the gaussian
!  function of x with amplitude a, center (mean) b, and
!  width (standard deviation) c.
!
  else if ( k == 8 ) then

    t1 = 5.0 - 10.0 * x
    t2 = 5.0 - 10.0 * y
    t3 = exp ( - 0.5 * t1 * t1 )
    t4 = exp ( - 0.5 * t2 * t2 )
    f = t3 + 0.75 * t4 * ( 1.0 + t3 )

    if ( 1 <= flag ) then
      fx = t1 * t3 * ( 10.0 + 7.5 * t4 )
      fy = t2 * t4 * ( 7.5 + 7.5 * t3 )
    end if

    if ( 2 <= flag ) then
      fxx = t3 * ( t1 * t1 - 1.0 ) * ( 100.0 + 75.0 * t4 )
      fxy = 75.0 * t1 * t2 * t3 * t4
      fyy = t4 * ( t2 * t2 - 1.0 ) * ( 75.0 + 75.0 * t3 )
    end if
!
!  Cloverleaf asymmetric hill/valley:
!
  else if ( k == 9 ) then

    t1 = exp ( ( 10.0 - 20.0 * x ) / 3.0 )
    t2 = exp ( ( 10.0 - 20.0 * y ) / 3.0 )
    t3 = 1.0 / ( 1.0 + t1 )
    t4 = 1.0 / ( 1.0 + t2 )
    f = ( ( 20.0 / 3.0 )**3 * t1 * t2 )**2 * ( t3 * t4 )**5 * &
      ( t1 - 2.0 * t3 ) * ( t2 - 2.0 * t4 )

    if ( 1 <= flag ) then

      fx = ( ( 20.0 / 3.0 ) * t1 )**2 * ( ( 20.0 / 3.0 ) * t3 )**5 * &
        ( 2.0 * t1 - 3.0 * t3 - 5.0 + 12.0 * t3 * t3 ) * t2 * t2 * t4**5 * &
        ( t2 - 2.0 * t4 )

      fy = ( ( 20.0 / 3.0 ) * t1 )**2 * ( ( 20.0 / 3.0 ) * t3 )**5 * &
        ( 2.0 * t2 - 3.0 * t4 - 5.0 + 12.0 * t4 * t4 ) * t2 * t2 * t4**5 * &
        ( t1 - 2.0 * t3 )

    end if

    if ( 2 <= flag ) then

      t5 = 20.0 / 3.0
      t6 = ( t5 * t1 * t2 )**2 * ( t5 * t3 * t4 )**5

      fxx = t5 * t6 * ( t2 - 2.0 * t4 ) * &
        ( ( ( - 84.0 * t3 + 78.0 ) * t3 + 23.0 ) * t3 + 4.0 * t1 - 25.0 )

      fxy = t5 * t6 * &
        ( ( 12.0 * t4 - 3.0 ) * t4 + 2.0 * t2 - 5.0 ) * &
        ( ( 12.0 * t3 - 3.0 ) * t3 + 2.0 * t1 - 5.0 )

      fyy = t5 * t6 * ( t1 - 2.0 * t3 ) * &
        ( ( ( - 84.0 * t4 + 78.0 ) * t4 + 23.0 ) * t4 + 4.0 * t2 - 25.0 )

    end if
!
!  Cosine peak:
!
  else if ( k == 10 ) then

    t1 = sqrt ( ( 80.0 * x - 40.0 )**2 + ( 90.0 * y - 45.0 )**2 )
    t2 = exp ( - 0.04 * t1 )
    t3 = cos ( 0.15 * t1 )
    f = t2 * t3

    if ( 1 <= flag ) then

      if ( t1 == 0.0 ) then
        fx = 0.0
        fy = 0.0
      else
        t4 = sin ( 0.15 * t1 )
        fx = - t2 * ( 12.0 * t4 + 3.2 * t3 ) * ( 80.0 * x - 40.0 ) / t1
        fy = - t2 * ( 13.5 * t4 + 3.6 * t3 ) * ( 90.0 * y - 45.0 ) / t1
      end if

    end if

    if ( 2 <= flag ) then
      if ( t1 == 0.0 ) then
        fxx = 0.0
        fxy = 0.0
        fyy = 0.0
      else
        t5 = t2 / ( t1**3 )
        fxx = t5 * ( t1 * ( 76.8 * t4 - 133.76 * t3 ) * ( 80.0 * x - 40.0 )**2 - &
          ( 960.0 * t4 + 256.0 * t3 ) * ( 90.0 * y - 45.0 )**2 )

        fxy = t5 * ( t1 * ( 86.4 * t4 - 150.48 * t3 ) + 1080.0 * t4 + &
          288.0 * t3 ) * ( 80.0 * x - 40.0 ) * ( 90.0 * y - 45.0 )

        fyy = t5 * ( t1 * ( 97.2 * t4 - 169.29 * t3 ) * ( 90.0 * y - 45.0 )**2 - &
          ( 1215.0 * t4 + 324.0 * t3 ) * ( 80.0 * x - 40.0 )**2 )
      end if
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TSTFN2 - Fatal error!'
    write ( *, '(a)' ) '  Input K does not satisfy 1 <= K <= 10.'
    stop

  end if

  return
end
