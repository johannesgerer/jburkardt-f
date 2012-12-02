subroutine angle_shift ( alpha, beta, gamma )

!*****************************************************************************80
!
!! ANGLE_SHIFT shifts angle ALPHA to lie between BETA and BETA+2PI.
!
!  Discussion:
!
!    The input angle ALPHA is shifted by multiples of 2 * PI to lie
!    between BETA and BETA+2*PI.
!
!    The resulting angle GAMMA has all the same trigonometric function
!    values as ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the angle to be shifted.
!
!    Input, real ( kind = 8 ) BETA, defines the lower endpoint of
!    the angle range.
!
!    Output, real ( kind = 8 ), GAMMA, the shifted angle.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  if ( alpha < beta ) then
    gamma = beta - mod ( beta - alpha, 2.0D+00 * pi ) + 2.0D+00 * pi
  else
    gamma = beta + mod ( alpha - beta, 2.0D+00 * pi )
  end if

  return
end
subroutine arccos_cordic ( t, n, theta )

!*****************************************************************************80
!
!! ARCCOS_CORDIC returns the arccosine of an angle using the CORDIC method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Michel Muller,
!    Elementary Functions: Algorithms and Implementation,
!    Second Edition,
!    Birkhaeuser, 2006,
!    ISBN13: 978-0-8176-4372-0,
!    LC: QA331.M866.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the cosine of an angle.  -1 <= T <= 1.
!
!    Input, integer ( kind = 4 ) N, the number of iterations to take.
!    A value of 10 is low.  Good accuracy is achieved with 20 or more
!    iterations.
!
!    Output, real ( kind = 8 ) THETA, an angle whose cosine is T.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ANGLES(60) = arctan ( (1/2)^(0:59) );
!
  implicit none

  integer ( kind = 4 ), parameter :: angles_length = 60

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( angles_length ) :: angles = (/ &
    7.8539816339744830962D-01, &
    4.6364760900080611621D-01, &
    2.4497866312686415417D-01, &
    1.2435499454676143503D-01, &
    6.2418809995957348474D-02, &
    3.1239833430268276254D-02, &
    1.5623728620476830803D-02, &
    7.8123410601011112965D-03, &
    3.9062301319669718276D-03, &
    1.9531225164788186851D-03, &
    9.7656218955931943040D-04, &
    4.8828121119489827547D-04, &
    2.4414062014936176402D-04, &
    1.2207031189367020424D-04, &
    6.1035156174208775022D-05, &
    3.0517578115526096862D-05, &
    1.5258789061315762107D-05, &
    7.6293945311019702634D-06, &
    3.8146972656064962829D-06, &
    1.9073486328101870354D-06, &
    9.5367431640596087942D-07, &
    4.7683715820308885993D-07, &
    2.3841857910155798249D-07, &
    1.1920928955078068531D-07, &
    5.9604644775390554414D-08, &
    2.9802322387695303677D-08, &
    1.4901161193847655147D-08, &
    7.4505805969238279871D-09, &
    3.7252902984619140453D-09, &
    1.8626451492309570291D-09, &
    9.3132257461547851536D-10, &
    4.6566128730773925778D-10, &
    2.3283064365386962890D-10, &
    1.1641532182693481445D-10, &
    5.8207660913467407226D-11, &
    2.9103830456733703613D-11, &
    1.4551915228366851807D-11, &
    7.2759576141834259033D-12, &
    3.6379788070917129517D-12, &
    1.8189894035458564758D-12, &
    9.0949470177292823792D-13, &
    4.5474735088646411896D-13, &
    2.2737367544323205948D-13, &
    1.1368683772161602974D-13, &
    5.6843418860808014870D-14, &
    2.8421709430404007435D-14, &
    1.4210854715202003717D-14, &
    7.1054273576010018587D-15, &
    3.5527136788005009294D-15, &
    1.7763568394002504647D-15, &
    8.8817841970012523234D-16, &
    4.4408920985006261617D-16, &
    2.2204460492503130808D-16, &
    1.1102230246251565404D-16, &
    5.5511151231257827021D-17, &
    2.7755575615628913511D-17, &
    1.3877787807814456755D-17, &
    6.9388939039072283776D-18, &
    3.4694469519536141888D-18, &
    1.7347234759768070944D-18 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) poweroftwo
  real ( kind = 8 ) r(2,2)
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sign_z2
  real ( kind = 8 ) t
  real ( kind = 8 ) t_copy
  real ( kind = 8 ) theta
  real ( kind = 8 ) z(2)

  if ( t < -1.0D+00 .or. 1.0D+00 < t ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ARCCOS_CORDIC - Fatal error!'
    write ( *, '(a)' ) '  1.0 < |T|.'
    stop
  end if

  theta = 0.0D+00
  z = (/ 1.0D+00, 0.0D+00 /)
  poweroftwo = 1.0D+00
  r(1,1) = 1.0D+00
  r(2,2) = 1.0D+00
  t_copy = t

  do j = 1, n

    if ( z(2) < 0.0D+00 ) then
      sign_z2 = -1.0D+00
    else
      sign_z2 = 1.0D+00
    end if

    if ( t_copy <= z(1) ) then
      sigma = + sign_z2
    else
      sigma = - sign_z2
    end if

    if ( j <= angles_length ) then
      angle = angles(j)
    else
      angle = angle / 2.0D+00
    end if

    r(1,2) = - sigma * poweroftwo
    r(2,1) = + sigma * poweroftwo

    z = matmul ( r, matmul ( r, z ) )

    theta = theta + 2.0D+00 * sigma * angle

    t_copy = t_copy + t_copy * poweroftwo * poweroftwo

    poweroftwo = poweroftwo / 2.0D+00

  end do

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
subroutine arcsin_cordic ( t, n, theta )

!*****************************************************************************80
!
!! ARCSIN_CORDIC returns the arcsine of an angle using the CORDIC method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Michel Muller,
!    Elementary Functions: Algorithms and Implementation,
!    Second Edition,
!    Birkhaeuser, 2006,
!    ISBN13: 978-0-8176-4372-0,
!    LC: QA331.M866.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the sine of an angle.  -1 <= T <= 1.
!
!    Input, integer ( kind = 4 ) N, the number of iterations to take.
!    A value of 10 is low.  Good accuracy is achieved with 20 or more
!    iterations.
!
!    Output, real ( kind = 8 ) THETA, an angle whose sine is T.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ANGLES(60) = arctan ( (1/2)^(0:59) );
!
  implicit none

  integer ( kind = 4 ), parameter :: angles_length = 60

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( angles_length ) :: angles = (/ &
    7.8539816339744830962D-01, &
    4.6364760900080611621D-01, &
    2.4497866312686415417D-01, &
    1.2435499454676143503D-01, &
    6.2418809995957348474D-02, &
    3.1239833430268276254D-02, &
    1.5623728620476830803D-02, &
    7.8123410601011112965D-03, &
    3.9062301319669718276D-03, &
    1.9531225164788186851D-03, &
    9.7656218955931943040D-04, &
    4.8828121119489827547D-04, &
    2.4414062014936176402D-04, &
    1.2207031189367020424D-04, &
    6.1035156174208775022D-05, &
    3.0517578115526096862D-05, &
    1.5258789061315762107D-05, &
    7.6293945311019702634D-06, &
    3.8146972656064962829D-06, &
    1.9073486328101870354D-06, &
    9.5367431640596087942D-07, &
    4.7683715820308885993D-07, &
    2.3841857910155798249D-07, &
    1.1920928955078068531D-07, &
    5.9604644775390554414D-08, &
    2.9802322387695303677D-08, &
    1.4901161193847655147D-08, &
    7.4505805969238279871D-09, &
    3.7252902984619140453D-09, &
    1.8626451492309570291D-09, &
    9.3132257461547851536D-10, &
    4.6566128730773925778D-10, &
    2.3283064365386962890D-10, &
    1.1641532182693481445D-10, &
    5.8207660913467407226D-11, &
    2.9103830456733703613D-11, &
    1.4551915228366851807D-11, &
    7.2759576141834259033D-12, &
    3.6379788070917129517D-12, &
    1.8189894035458564758D-12, &
    9.0949470177292823792D-13, &
    4.5474735088646411896D-13, &
    2.2737367544323205948D-13, &
    1.1368683772161602974D-13, &
    5.6843418860808014870D-14, &
    2.8421709430404007435D-14, &
    1.4210854715202003717D-14, &
    7.1054273576010018587D-15, &
    3.5527136788005009294D-15, &
    1.7763568394002504647D-15, &
    8.8817841970012523234D-16, &
    4.4408920985006261617D-16, &
    2.2204460492503130808D-16, &
    1.1102230246251565404D-16, &
    5.5511151231257827021D-17, &
    2.7755575615628913511D-17, &
    1.3877787807814456755D-17, &
    6.9388939039072283776D-18, &
    3.4694469519536141888D-18, &
    1.7347234759768070944D-18 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) poweroftwo
  real ( kind = 8 ) r(2,2)
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sign_z1
  real ( kind = 8 ) t
  real ( kind = 8 ) t_copy
  real ( kind = 8 ) theta
  real ( kind = 8 ) z(2)

  if ( t < -1.0D+00 .or. 1.0D+00 < t ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ARCSIN_CORDIC - Fatal error!'
    write ( *, '(a)' ) '  1.0 < |T|.'
    stop
  end if

  theta = 0.0D+00
  z = (/ 1.0D+00, 0.0D+00 /)
  poweroftwo = 1.0D+00
  r(1,1) = 1.0D+00
  r(2,2) = 1.0D+00
  t_copy = t

  do j = 1, n

    if ( z(1) < 0.0D+00 ) then
      sign_z1 = -1.0D+00
    else
      sign_z1 = 1.0D+00
    end if

    if ( z(2) <= t_copy ) then
      sigma = + sign_z1
    else
      sigma = - sign_z1
    end if

    if ( j <= angles_length ) then
      angle = angles(j)
    else
      angle = angle / 2.0D+00
    end if

    r(1,2) = - sigma * poweroftwo
    r(2,1) = + sigma * poweroftwo

    z = matmul ( r, matmul ( r, z ) )

    theta = theta + 2.0D+00 * sigma * angle

    t_copy = t_copy + t_copy * poweroftwo * poweroftwo

    poweroftwo = poweroftwo / 2.0D+00

  end do

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
subroutine arctan_cordic ( x, y, n, theta )

!*****************************************************************************80
!
!! ARCTAN_CORDIC returns the arctangent of an angle using the CORDIC method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Michel Muller,
!    Elementary Functions: Algorithms and Implementation,
!    Second Edition,
!    Birkhaeuser, 2006,
!    ISBN13: 978-0-8176-4372-0,
!    LC: QA331.M866.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, define the tangent of an angle as Y/X.
!
!    Input, integer ( kind = 4 )  N, the number of iterations to take.
!    A value of 10 is low.  Good accuracy is achieved with 20 or more
!    iterations.
!
!    Output, real ( kind = 8 ) THETA, the angle whose tangent is Y/X.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ANGLES(60) = arctan ( (1/2)^(0:59) );
!
  implicit none

  integer ( kind = 4 ), parameter :: angles_length = 60

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( angles_length ) :: angles = (/ &
    7.8539816339744830962D-01, &
    4.6364760900080611621D-01, &
    2.4497866312686415417D-01, &
    1.2435499454676143503D-01, &
    6.2418809995957348474D-02, &
    3.1239833430268276254D-02, &
    1.5623728620476830803D-02, &
    7.8123410601011112965D-03, &
    3.9062301319669718276D-03, &
    1.9531225164788186851D-03, &
    9.7656218955931943040D-04, &
    4.8828121119489827547D-04, &
    2.4414062014936176402D-04, &
    1.2207031189367020424D-04, &
    6.1035156174208775022D-05, &
    3.0517578115526096862D-05, &
    1.5258789061315762107D-05, &
    7.6293945311019702634D-06, &
    3.8146972656064962829D-06, &
    1.9073486328101870354D-06, &
    9.5367431640596087942D-07, &
    4.7683715820308885993D-07, &
    2.3841857910155798249D-07, &
    1.1920928955078068531D-07, &
    5.9604644775390554414D-08, &
    2.9802322387695303677D-08, &
    1.4901161193847655147D-08, &
    7.4505805969238279871D-09, &
    3.7252902984619140453D-09, &
    1.8626451492309570291D-09, &
    9.3132257461547851536D-10, &
    4.6566128730773925778D-10, &
    2.3283064365386962890D-10, &
    1.1641532182693481445D-10, &
    5.8207660913467407226D-11, &
    2.9103830456733703613D-11, &
    1.4551915228366851807D-11, &
    7.2759576141834259033D-12, &
    3.6379788070917129517D-12, &
    1.8189894035458564758D-12, &
    9.0949470177292823792D-13, &
    4.5474735088646411896D-13, &
    2.2737367544323205948D-13, &
    1.1368683772161602974D-13, &
    5.6843418860808014870D-14, &
    2.8421709430404007435D-14, &
    1.4210854715202003717D-14, &
    7.1054273576010018587D-15, &
    3.5527136788005009294D-15, &
    1.7763568394002504647D-15, &
    8.8817841970012523234D-16, &
    4.4408920985006261617D-16, &
    2.2204460492503130808D-16, &
    1.1102230246251565404D-16, &
    5.5511151231257827021D-17, &
    2.7755575615628913511D-17, &
    1.3877787807814456755D-17, &
    6.9388939039072283776D-18, &
    3.4694469519536141888D-18, &
    1.7347234759768070944D-18 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) poweroftwo
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sign_factor
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  x1 = x
  y1 = y
!
!  Account for signs.
!
  if ( x1 < 0.0D+00 .and. y1 < 0.0D+00 ) then
    x1 = -x1
    y1 = -y1
  end if

  if ( x1 < 0.0D+00 ) then
    x1 = -x1
    sign_factor = -1.0D+00
  else if ( y1 < 0.0D+00 ) then
    y1 = -y1
    sign_factor = -1.0D+00
  else
    sign_factor = +1.0D+00
  end if

  theta = 0.0D+00
  poweroftwo = 1.0D+00

  do j = 1, n

    if ( y1 <= 0.0D+00 ) then
      sigma = +1.0D+00
    else
      sigma = -1.0D+00
    end if

    if ( j <= angles_length ) then
      angle = angles(j)
    else
      angle = angle / 2.0D+00
    end if

    x2 =                      x1 - sigma * poweroftwo * y1
    y2 = sigma * poweroftwo * x1 +                      y1
    theta  = theta - sigma * angle

    x1 = x2
    y1 = y2

    poweroftwo = poweroftwo / 2.0D+00

  end do

  theta = sign_factor * theta

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
subroutine cbrt_cordic ( x, n, y )

!*****************************************************************************80
!
!! CBRT_CORDIC returns the cube root of a value using the CORDIC method.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose square root is desired.
!
!    Input, integer ( kind = 4 ) N, the number of iterations to take.
!    This is essentially the number of binary digits of accuracy, and
!    might go as high as 53.
!
!    Output, real ( kind = 8 ) Y, the approximate cube root of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 4 ) poweroftwo
  real ( kind = 8 ) x
  real ( kind = 8 ) x_mag
  real ( kind = 8 ) y

  x_mag = abs ( x )

  if ( x == 0.0D+00 ) then
    y = 0.0D+00
    return
  end if

  if ( x_mag == 1.0D+00 ) then
    y = x
    return
  end if

  poweroftwo = 1.0D+00

  if ( x_mag < 1.0D+00 ) then

    do while ( x_mag <= poweroftwo * poweroftwo * poweroftwo )
      poweroftwo = poweroftwo / 2.0D+00
    end do

    y = poweroftwo

  else if ( 1.0D+00 < x_mag ) then

    do while ( poweroftwo * poweroftwo * poweroftwo <= x_mag )
      poweroftwo = 2.0D+00 * poweroftwo
    end do

    y = poweroftwo / 2.0D+00

  end if

  do i = 1, n
    poweroftwo = poweroftwo / 2.0D+00
    if ( ( y + poweroftwo ) * ( y + poweroftwo ) * ( y + poweroftwo ) &
      <= x_mag ) then
      y = y + poweroftwo
    end if
  end do

  if ( x < 0.0D+00 ) then
    y = -y
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
!    Output real ( kind = 8 ) FX, the value of the function.
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
subroutine cossin_cordic ( beta, n, c, s )

!*****************************************************************************80
!
!! COSSIN_CORDIC returns the sine and cosine using the CORDIC method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2007
!
!  Author:
!
!    Based on MATLAB code in a Wikipedia article.
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BETA, the angle (in radians).
!
!    Input, integer ( kind = 4 ) N, the number of iterations to take.
!    A value of 10 is low.  Good accuracy is achieved with 20 or more
!    iterations.
!
!    Output, real ( kind = 8 ) C, S, the cosine and sine of the angle.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ANGLES(60) = arctan ( (1/2)^(0:59) );
!
!    Local, real ( kind = 8 ) KPROD(33).
!    KPROD(j) = product ( 0 <= i <= j ) K(i)
!    where K(i) = 1 / sqrt ( 1 + (1/2)^(2i) ).
!
  implicit none

  integer ( kind = 4 ), parameter :: angles_length = 60
  integer ( kind = 4 ), parameter :: kprod_length = 33

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( angles_length ) :: angles = (/ &
    7.8539816339744830962D-01, &
    4.6364760900080611621D-01, &
    2.4497866312686415417D-01, &
    1.2435499454676143503D-01, &
    6.2418809995957348474D-02, &
    3.1239833430268276254D-02, &
    1.5623728620476830803D-02, &
    7.8123410601011112965D-03, &
    3.9062301319669718276D-03, &
    1.9531225164788186851D-03, &
    9.7656218955931943040D-04, &
    4.8828121119489827547D-04, &
    2.4414062014936176402D-04, &
    1.2207031189367020424D-04, &
    6.1035156174208775022D-05, &
    3.0517578115526096862D-05, &
    1.5258789061315762107D-05, &
    7.6293945311019702634D-06, &
    3.8146972656064962829D-06, &
    1.9073486328101870354D-06, &
    9.5367431640596087942D-07, &
    4.7683715820308885993D-07, &
    2.3841857910155798249D-07, &
    1.1920928955078068531D-07, &
    5.9604644775390554414D-08, &
    2.9802322387695303677D-08, &
    1.4901161193847655147D-08, &
    7.4505805969238279871D-09, &
    3.7252902984619140453D-09, &
    1.8626451492309570291D-09, &
    9.3132257461547851536D-10, &
    4.6566128730773925778D-10, &
    2.3283064365386962890D-10, &
    1.1641532182693481445D-10, &
    5.8207660913467407226D-11, &
    2.9103830456733703613D-11, &
    1.4551915228366851807D-11, &
    7.2759576141834259033D-12, &
    3.6379788070917129517D-12, &
    1.8189894035458564758D-12, &
    9.0949470177292823792D-13, &
    4.5474735088646411896D-13, &
    2.2737367544323205948D-13, &
    1.1368683772161602974D-13, &
    5.6843418860808014870D-14, &
    2.8421709430404007435D-14, &
    1.4210854715202003717D-14, &
    7.1054273576010018587D-15, &
    3.5527136788005009294D-15, &
    1.7763568394002504647D-15, &
    8.8817841970012523234D-16, &
    4.4408920985006261617D-16, &
    2.2204460492503130808D-16, &
    1.1102230246251565404D-16, &
    5.5511151231257827021D-17, &
    2.7755575615628913511D-17, &
    1.3877787807814456755D-17, &
    6.9388939039072283776D-18, &
    3.4694469519536141888D-18, &
    1.7347234759768070944D-18 /)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) factor
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( kprod_length ) :: kprod = (/ &
    0.70710678118654752440D+00, &
    0.63245553203367586640D+00, &
    0.61357199107789634961D+00, &
    0.60883391251775242102D+00, &
    0.60764825625616820093D+00, &
    0.60735177014129595905D+00, &
    0.60727764409352599905D+00, &
    0.60725911229889273006D+00, &
    0.60725447933256232972D+00, &
    0.60725332108987516334D+00, &
    0.60725303152913433540D+00, &
    0.60725295913894481363D+00, &
    0.60725294104139716351D+00, &
    0.60725293651701023413D+00, &
    0.60725293538591350073D+00, &
    0.60725293510313931731D+00, &
    0.60725293503244577146D+00, &
    0.60725293501477238499D+00, &
    0.60725293501035403837D+00, &
    0.60725293500924945172D+00, &
    0.60725293500897330506D+00, &
    0.60725293500890426839D+00, &
    0.60725293500888700922D+00, &
    0.60725293500888269443D+00, &
    0.60725293500888161574D+00, &
    0.60725293500888134606D+00, &
    0.60725293500888127864D+00, &
    0.60725293500888126179D+00, &
    0.60725293500888125757D+00, &
    0.60725293500888125652D+00, &
    0.60725293500888125626D+00, &
    0.60725293500888125619D+00, &
    0.60725293500888125617D+00 /)
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) poweroftwo
  real ( kind = 8 ) s
  real ( kind = 8 ) s2
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sign_factor
  real ( kind = 8 ) theta
!
!  Shift angle to interval [-pi,pi].
!
  call angle_shift ( beta, -pi, theta )
!
!  Shift angle to interval [-pi/2,pi/2] and account for signs.
!
  if ( theta < - 0.5D+00 * pi ) then
    theta = theta + pi
    sign_factor = -1.0D+00
  else if ( 0.5D+00 * pi < theta ) then
    theta = theta - pi
    sign_factor = -1.0D+00
  else
    sign_factor = +1.0D+00
  end if

  c = 1.0D+00
  s = 0.0D+00

  poweroftwo = 1.0D+00
  angle = angles(1)

  do j = 1, n

    if ( theta < 0.0D+00 ) then
      sigma = -1.0D+00
    else
      sigma = 1.0D+00
    end if

    factor = sigma * poweroftwo

    c2 =          c - factor * s
    s2 = factor * c +          s

    c = c2
    s = s2
!
!  Update the remaining angle.
!
    theta = theta - sigma * angle

    poweroftwo = poweroftwo / 2.0D+00
!
!  Update the angle from table, or eventually by just dividing by two.
!
    if ( angles_length < j + 1 ) then
      angle = angle / 2.0D+00
    else
      angle = angles(j+1)
    end if

  end do
!
!  Adjust length of output vector to be [cos(beta), sin(beta)]
!
!  KPROD is essentially constant after a certain point, so if N is
!  large, just take the last available value.
!
  if ( 0 < n ) then
    c = c * kprod ( min ( n, kprod_length ) )
    s = s * kprod ( min ( n, kprod_length ) )
  end if
!
!  Adjust for possible sign change because angle was originally
!  not in quadrant 1 or 4.
!
  c = sign_factor * c
  s = sign_factor * s

  return
end
subroutine exp_cordic ( x, n, fx )

!*****************************************************************************80
!
!! EXP_CORDIC evaluates the exponential function using the CORDIC method.
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
!    Frederick Ruckdeschel,
!    BASIC Scientific Subroutines,
!    Volume II,
!    McGraw-Hill, 1980,
!    ISBN: 0-07-054202-3,
!    LC: QA76.95.R82.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) N, the number of steps to take.
!
!    Output, real ( kind = 8 ) FX, the exponential of X.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(1:25) = exp ( (1/2)^(1:25) );
!
  implicit none

  integer ( kind = 4 ), parameter :: a_length = 25
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( a_length ) :: a = (/ &
    1.648721270700128D+00, &
    1.284025416687742D+00, &
    1.133148453066826D+00, &
    1.064494458917859D+00, &
    1.031743407499103D+00, &
    1.015747708586686D+00, &
    1.007843097206488D+00, &
    1.003913889338348D+00, &
    1.001955033591003D+00, &
    1.000977039492417D+00, &
    1.000488400478694D+00, &
    1.000244170429748D+00, &
    1.000122077763384D+00, &
    1.000061037018933D+00, &
    1.000030518043791D+00, &
    1.0000152589054785D+00, &
    1.0000076294236351D+00, &
    1.0000038147045416D+00, &
    1.0000019073504518D+00, &
    1.0000009536747712D+00, &
    1.0000004768372719D+00, &
    1.0000002384186075D+00, &
    1.0000001192092967D+00, &
    1.0000000596046466D+00, &
    1.0000000298023228D+00 /)
  real ( kind = 8 ) ai
  real ( kind = 8 ), parameter :: e = 2.718281828459045D+00
  real ( kind = 8 ) fx
  integer ( kind = 4 ) i
  real ( kind = 8 ) poweroftwo
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x
  integer ( kind = 4 ) x_int
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  x_int = floor ( x )
!
!  Determine the weights.
!
  poweroftwo = 0.5D+00
  z = x - real ( x_int, kind = 8 )

  do i = 1, n
    w(i) = 0.0D+00
    if ( poweroftwo < z ) then
      w(i) = 1.0D+00
      z = z - poweroftwo
    end if
    poweroftwo = poweroftwo / 2.0D+00
  end do
!
!  Calculate products.
!
  fx = 1.0D+00

  do i = 1, n

    if ( i <= a_length ) then
      ai = a(i)
    else
      ai = 1.0D+00 + ( ai - 1.0D+00 ) / 2.0D+00
    end if

    if ( 0.0D+00 < w(i) ) then
      fx = fx * ai
    end if

  end do
!
!  Perform residual multiplication.
!
  fx = fx                     &
    * ( 1.0D+00 + z           &
    * ( 1.0D+00 + z / 2.0D+00 &
    * ( 1.0D+00 + z / 3.0D+00 &
    * ( 1.0D+00 + z / 4.0D+00 ))))
!
!  Account for factor EXP(X_INT).
!
  if ( x_int < 0 ) then

    do i = 1, -x_int
      fx = fx / e
    end do

  else

    do i = 1, x_int
      fx = fx * e
    end do

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

  integer ( kind = 4 ), parameter :: n_max = 19

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.000045399929762484851536D+00, &
    0.0067379469990854670966D+00, &
    0.36787944117144232160D+00, &
    1.0000000000000000000D+00, &
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
    4.8516519540979027797D+08 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     -10.0D+00, &
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
       3.1415926535897932385D+00, &
       5.0D+00, &
      10.0D+00, &
      20.0D+00 /)

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
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
subroutine ln_cordic ( x, n, fx )

!*****************************************************************************80
!
!! LN_CORDIC evaluates the natural logarithm using the CORDIC method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Frederick Ruckdeschel,
!    BASIC Scientific Subroutines,
!    Volume II,
!    McGraw-Hill, 1980,
!    ISBN: 0-07-054202-3,
!    LC: QA76.95.R82.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) N, the number of steps to take.
!
!    Output, real ( kind = 8 ) FX, the natural logarithm of X.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(1:25) = exp ( (1/2)^(1:25) );
!
  implicit none

  integer ( kind = 4 ), parameter :: a_length = 25
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( a_length ) :: a = (/ &
    1.648721270700128D+00, &
    1.284025416687742D+00, &
    1.133148453066826D+00, &
    1.064494458917859D+00, &
    1.031743407499103D+00, &
    1.015747708586686D+00, &
    1.007843097206488D+00, &
    1.003913889338348D+00, &
    1.001955033591003D+00, &
    1.000977039492417D+00, &
    1.000488400478694D+00, &
    1.000244170429748D+00, &
    1.000122077763384D+00, &
    1.000061037018933D+00, &
    1.000030518043791D+00, &
    1.0000152589054785D+00, &
    1.0000076294236351D+00, &
    1.0000038147045416D+00, &
    1.0000019073504518D+00, &
    1.0000009536747712D+00, &
    1.0000004768372719D+00, &
    1.0000002384186075D+00, &
    1.0000001192092967D+00, &
    1.0000000596046466D+00, &
    1.0000000298023228D+00 /)
  real ( kind = 8 ) ai
  real ( kind = 8 ), parameter :: e = 2.718281828459045D+00
  real ( kind = 8 ) fx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) poweroftwo
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x
  real ( kind = 8 ) x_copy

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LN_CORDIC - Fatal error!'
    write ( *, '(a)' ) '  Input argument X <= 0.0'
    stop
  end if

  x_copy = x

  k = 0

  do while ( e <= x_copy )
    k = k + 1
    x_copy = x_copy / e
  end do

  do while ( x_copy < 1.0D+00 )
    k = k - 1
    x_copy = x_copy * e
  end do
!
!  Determine the weights.
!
  do i = 1, n

    w(i) = 0.0D+00

    if ( i <= a_length ) then
      ai = a(i)
    else
      ai = 1.0D+00 + ( ai - 1.0D+00 ) / 2.0D+00
    end if

    if ( ai < x_copy ) then
      w(i) = 1.0D+00
      x_copy = x_copy / ai
    end if

  end do

  x_copy = x_copy - 1.0D+00

  x_copy = x_copy &
    * ( 1.0D+00 - ( x_copy / 2.0D+00 ) &
    * ( 1.0D+00 + ( x_copy / 3.0D+00 ) &
    * ( 1.0D+00 -   x_copy / 4.0D+00 )))
!
!  Assemble
!
  poweroftwo = 0.5D+00

  do i = 1, n
    x_copy = x_copy + w(i) * poweroftwo
    poweroftwo = poweroftwo / 2.0D+00
  end do

  fx = real ( k, kind = 8 ) + x_copy

  return
end
subroutine ln_values ( n_data, x, fx )

!*****************************************************************************80
!
!! LN_VALUES returns some values of the natural logarithm function.
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
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
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
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) i4_huge
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
    seed = seed + i4_huge ( )
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

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
subroutine sqrt_cordic ( x, n, y )

!*****************************************************************************80
!
!! SQRT_CORDIC returns the square root of a value using the CORDIC method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose square root is desired.
!
!    Input, integer ( kind = 4 ) N, the number of iterations to take.
!    This is essentially the number of binary digits of accuracy, and
!    might go as high as 53.
!
!    Output, real ( kind = 8 ) Y, the approximate square root of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 4 ) poweroftwo
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SQRT_CORDIC - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop
  end if

  if ( x == 0.0D+00 ) then
    y = 0.0D+00
    return
  end if

  if ( x == 1.0D+00 ) then
    y = 1.0D+00
    return
  end if

  poweroftwo = 1.0D+00

  if ( x < 1.0D+00 ) then

    do while ( x <= poweroftwo * poweroftwo )
      poweroftwo = poweroftwo / 2.0D+00
    end do

    y = poweroftwo

  else if ( 1.0D+00 < x ) then

    do while ( poweroftwo * poweroftwo <= x )
      poweroftwo = 2.0D+00 * poweroftwo
    end do

    y = poweroftwo / 2.0D+00

  end if

  do i = 1, n
    poweroftwo = poweroftwo / 2.0D+00
    if ( ( y + poweroftwo ) * ( y + poweroftwo ) <= x ) then
      y = y + poweroftwo
    end if
  end do

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
subroutine tan_cordic ( beta, n, t )

!*****************************************************************************80
!
!! TAN_CORDIC returns the tangent of an angle using the CORDIC method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BETA, the angle (in radians).
!
!    Input, integer ( kind = 4 ) N, the number of iterations to take.
!    A value of 10 is low.  Good accuracy is achieved with 20 or more
!    iterations.
!
!    Output, real ( kind = 8 ) T, the tangent of the angle.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ANGLES(60) = arctan ( (1/2)^(0:59) );
!
  implicit none

  integer ( kind = 4 ), parameter :: angles_length = 60

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( angles_length ) :: angles = (/ &
    7.8539816339744830962D-01, &
    4.6364760900080611621D-01, &
    2.4497866312686415417D-01, &
    1.2435499454676143503D-01, &
    6.2418809995957348474D-02, &
    3.1239833430268276254D-02, &
    1.5623728620476830803D-02, &
    7.8123410601011112965D-03, &
    3.9062301319669718276D-03, &
    1.9531225164788186851D-03, &
    9.7656218955931943040D-04, &
    4.8828121119489827547D-04, &
    2.4414062014936176402D-04, &
    1.2207031189367020424D-04, &
    6.1035156174208775022D-05, &
    3.0517578115526096862D-05, &
    1.5258789061315762107D-05, &
    7.6293945311019702634D-06, &
    3.8146972656064962829D-06, &
    1.9073486328101870354D-06, &
    9.5367431640596087942D-07, &
    4.7683715820308885993D-07, &
    2.3841857910155798249D-07, &
    1.1920928955078068531D-07, &
    5.9604644775390554414D-08, &
    2.9802322387695303677D-08, &
    1.4901161193847655147D-08, &
    7.4505805969238279871D-09, &
    3.7252902984619140453D-09, &
    1.8626451492309570291D-09, &
    9.3132257461547851536D-10, &
    4.6566128730773925778D-10, &
    2.3283064365386962890D-10, &
    1.1641532182693481445D-10, &
    5.8207660913467407226D-11, &
    2.9103830456733703613D-11, &
    1.4551915228366851807D-11, &
    7.2759576141834259033D-12, &
    3.6379788070917129517D-12, &
    1.8189894035458564758D-12, &
    9.0949470177292823792D-13, &
    4.5474735088646411896D-13, &
    2.2737367544323205948D-13, &
    1.1368683772161602974D-13, &
    5.6843418860808014870D-14, &
    2.8421709430404007435D-14, &
    1.4210854715202003717D-14, &
    7.1054273576010018587D-15, &
    3.5527136788005009294D-15, &
    1.7763568394002504647D-15, &
    8.8817841970012523234D-16, &
    4.4408920985006261617D-16, &
    2.2204460492503130808D-16, &
    1.1102230246251565404D-16, &
    5.5511151231257827021D-17, &
    2.7755575615628913511D-17, &
    1.3877787807814456755D-17, &
    6.9388939039072283776D-18, &
    3.4694469519536141888D-18, &
    1.7347234759768070944D-18 /)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) factor
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) poweroftwo
  real ( kind = 8 ) s
  real ( kind = 8 ) s2
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t
  real ( kind = 8 ) theta
!
!  Shift angle to interval [-pi,pi].
!
  call angle_shift ( beta, -pi, theta )
!
!  Shift angle to interval [-pi/2,pi/2].
!
  if ( theta < - 0.5D+00 * pi ) then
    theta = theta + pi
  else if ( 0.5D+00 * pi < theta ) then
    theta = theta - pi
  end if

  c = 1.0D+00
  s = 0.0D+00

  poweroftwo = 1.0D+00
  angle = angles(1)

  do j = 1, n

    if ( theta < 0.0D+00 ) then
      sigma = -1.0D+00
    else
      sigma = +1.0D+00
    end if

    factor = sigma * poweroftwo

    c2 =          c - factor * s
    s2 = factor * c +          s

    c = c2
    s = s2
!
!  Update the remaining angle.
!
    theta = theta - sigma * angle

    poweroftwo = poweroftwo / 2.0D+00
!
!  Update the angle from table, or eventually by just dividing by two.
!
    if ( angles_length < j + 1 ) then
      angle = angle / 2.0D+00
    else
      angle = angles(j+1)
    end if

  end do

  t = s / c

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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
