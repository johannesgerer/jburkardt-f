subroutine p00_data ( prob, dim_num, data_num, p_data )

!*****************************************************************************80
!
!! P00_DATA returns the data for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the dependent
!    variables.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Output, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) p_data(dim_num,data_num)
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_data ( dim_num, data_num, p_data )
  else if ( prob == 2 ) then
    call p02_data ( dim_num, data_num, p_data )
  else if ( prob == 3 ) then
    call p03_data ( dim_num, data_num, p_data )
  else if ( prob == 4 ) then
    call p04_data ( dim_num, data_num, p_data )
  else if ( prob == 5 ) then
    call p05_data ( dim_num, data_num, p_data )
  else if ( prob == 6 ) then
    call p06_data ( dim_num, data_num, p_data )
  else if ( prob == 7 ) then
    call p07_data ( dim_num, data_num, p_data )
  else if ( prob == 8 ) then
    call p08_data ( dim_num, data_num, p_data )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_DATA - Fatal error!'
    write ( *, '(a)' ) '  Unexpected input value of PROB.'
    stop
  end if

  return
end
subroutine p00_data_num ( prob, data_num )

!*****************************************************************************80
!
!! P00_DATA_NUM returns the number of data points for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Output, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_data_num ( data_num )
  else if ( prob == 2 ) then
    call p02_data_num ( data_num )
  else if ( prob == 3 ) then
    call p03_data_num ( data_num )
  else if ( prob == 4 ) then
    call p04_data_num ( data_num )
  else if ( prob == 5 ) then
    call p05_data_num ( data_num )
  else if ( prob == 6 ) then
    call p06_data_num ( data_num )
  else if ( prob == 7 ) then
    call p07_data_num ( data_num )
  else if ( prob == 8 ) then
    call p08_data_num ( data_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_DATA_NUM - Fatal error!'
    write ( *, '(a)' ) '  Unexpected input value of PROB.'
    stop
  end if

  return
end
subroutine p00_dim_num ( prob, dim_num )

!*****************************************************************************80
!
!! P00_DIM_NUM returns the spatial dimension for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the
!    dependent variables.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_dim_num ( dim_num )
  else if ( prob == 2 ) then
    call p02_dim_num ( dim_num )
  else if ( prob == 3 ) then
    call p03_dim_num ( dim_num )
  else if ( prob == 4 ) then
    call p04_dim_num ( dim_num )
  else if ( prob == 5 ) then
    call p05_dim_num ( dim_num )
  else if ( prob == 6 ) then
    call p06_dim_num ( dim_num )
  else if ( prob == 7 ) then
    call p07_dim_num ( dim_num )
  else if ( prob == 8 ) then
    call p08_dim_num ( dim_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_DIM_NUM - Fatal error!'
    write ( *, '(a)' ) '  Unexpected input value of PROB.'
    stop
  end if

  return
end
subroutine p00_prob_num ( prob_num )

!*****************************************************************************80
!
!! P00_PROB_NUM returns the number of test problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROB_NUM, the number of test problems.
!
  implicit none

  integer ( kind = 4 ) prob_num

  prob_num = 8

  return
end
subroutine p00_story ( prob )

!*****************************************************************************80
!
!! P00_STORY prints the "story" for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2012
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

  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_story ( )
  else if ( prob == 2 ) then
    call p02_story ( )
  else if ( prob == 3 ) then
    call p03_story ( )
  else if ( prob == 4 ) then
    call p04_story ( )
  else if ( prob == 5 ) then
    call p05_story ( )
  else if ( prob == 6 ) then
    call p06_story ( )
  else if ( prob == 7 ) then
    call p07_story ( )
  else if ( prob == 8 ) then
    call p08_story ( )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_STORY - Fatal error!'
    write ( *, '(a)' ) '  Unexpected input value of PROB.'
    stop
  end if

  return
end
subroutine p01_data ( dim_num, data_num, p_data )

!*****************************************************************************80
!
!! P01_DATA returns the data for problem p01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the dependent
!    variables.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Output, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) p_data(dim_num,data_num)

  p_data = reshape ( (/ &
    0.0D+00, 4.0D+00, &
    1.0D+00, 5.0D+00, &
    2.0D+00, 6.0D+00, &
    4.0D+00, 6.0D+00, &
    5.0D+00, 5.0D+00, &
    6.0D+00, 3.0D+00, &
    7.0D+00, 1.0D+00, &
    8.0D+00, 1.0D+00, &
    9.0D+00, 1.0D+00, &
   10.0D+00, 3.0D+00, &
   11.0D+00, 4.0D+00, &
   12.0D+00, 4.0D+00, &
   13.0D+00, 3.0D+00, &
   14.0D+00, 3.0D+00, &
   15.0D+00, 4.0D+00, &
   16.0D+00, 4.0D+00, &
   17.0D+00, 3.0D+00, &
   18.0D+00, 0.0D+00  &
  /), (/ 2, 18 /) )

 return
end
subroutine p01_data_num ( data_num )

!*****************************************************************************80
!
!! P01_DATA_NUM returns the number of data points for problem p01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
  implicit none

  integer ( kind = 4 ) data_num

  data_num = 18

  return
end
subroutine p01_dim_num ( dim_num )

!*****************************************************************************80
!
!! P01_DIM_NUM returns the spatial dimension for problem p01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the 
!    dependent variables.
!
  implicit none

  integer ( kind = 4 ) dim_num

  dim_num = 2

  return
end
subroutine p01_story ( )

!*****************************************************************************80
!
!! P01_STORY prints the "story" for problem p01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hans-Joerg Wenz,
!    Interpolation of Curve Data by Blended Generalized Circles,
!    Computer Aided Geometric Design,
!    Volume 13, Number 8, November 1996, pages 673-680.
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  This example is due to Hans-Joerg Wenz.'
  write ( *, '(a)' ) &
    '  It is an example of good data, which is dense enough in areas'
  write ( *, '(a)' ) &
    '  where the expected curvature of the interpolant is large.'
  write ( *, '(a)' ) &
    '  Good results can be expected with almost any reasonable'
  write ( *, '(a)' ) &
    '  interpolation method.'

  return
end
subroutine p02_data ( dim_num, data_num, p_data )

!*****************************************************************************80
!
!! P02_DATA returns the data for problem p02.
!
!  Discussion:
!
!    Two pairs of identical X values have now been slightly separated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the dependent
!    variables.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Output, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) p_data(dim_num,data_num)

  p_data = reshape ( (/ &
     0.00D+00,   0.00D+00, &
     1.34D+00,   5.00D+00, &
     5.00D+00,   8.66D+00, &
    10.00D+00,  10.00D+00, &
    10.60D+00,  10.40D+00, &
    10.70D+00,  12.00D+00, &
    10.705D+00, 28.60D+00, &
    10.80D+00,  30.20D+00, &
    11.40D+00,  30.60D+00, &
    19.60D+00,  30.60D+00, &
    20.20D+00,  30.20D+00, &
    20.295D+00, 28.60D+00, &
    20.30D+00,  12.00D+00, &
    20.40D+00,  10.40D+00, &
    21.00D+00,  10.00D+00, &
    26.00D+00,   8.66D+00, &
    29.66D+00,   5.00D+00, &
    31.00D+00,   0.00D+00  &
  /), (/ dim_num, data_num /) )

  return
end
subroutine p02_data_num ( data_num )

!*****************************************************************************80
!
!! P02_DATA_NUM returns the number of data points for problem p02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
  implicit none

  integer ( kind = 4 ) data_num

  data_num = 18

  return
end
subroutine p02_dim_num ( dim_num )

!*****************************************************************************80
!
!! P02_DIM_NUM returns the spatial dimension for problem p02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the 
!    dependent variables.
!
  implicit none

  integer ( kind = 4 ) dim_num

  dim_num = 2

  return
end
subroutine p02_story ( )

!*****************************************************************************80
!
!! P02_STORY prints the "story" for problem p02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ETY Lee,
!    Choosing Nodes in Parametric Curve Interpolation,
!    Computer-Aided Design,
!    Volume 21, Number 6, July/August 1989, pages 363-370.
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  This example is due to ETY Lee of Boeing.'
  write ( *, '(a)' ) &
    '  Data near the corners is more dense than in regions of small curvature.'
  write ( *, '(a)' ) &
    '  A local interpolation method will produce a more plausible'
  write ( *, '(a)' ) &
    '  interpolant than a nonlocal interpolation method, such as'
  write ( *, '(a)' ) &
    '  cubic splines.'

  return
end
subroutine p03_data ( dim_num, data_num, p_data )

!*****************************************************************************80
!
!! P03_DATA returns the data for problem p03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the dependent
!    variables.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Output, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) p_data(dim_num,data_num)

  p_data = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     2.0D+00, 10.0D+00, &
     3.0D+00, 10.0D+00, &
     5.0D+00, 10.0D+00, &
     6.0D+00, 10.0D+00, &
     8.0D+00, 10.0D+00, &
     9.0D+00, 10.5D+00, &
    11.0D+00, 15.0D+00, &
    12.0D+00, 50.0D+00, &
    14.0D+00, 60.0D+00, &
    15.0D+00, 85.0D+00  &
  /), (/ dim_num, data_num /) )

  return
end
subroutine p03_data_num ( data_num )

!*****************************************************************************80
!
!! P03_DATA_NUM returns the number of data points for problem p03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
  implicit none

  integer ( kind = 4 ) data_num

  data_num = 11

  return
end
subroutine p03_dim_num ( dim_num )

!*****************************************************************************80
!
!! P03_DIM_NUM returns the spatial dimension for problem p03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the 
!    dependent variables.
!
  implicit none

  integer ( kind = 4 ) dim_num

  dim_num = 2

  return
end
subroutine p03_story ( )

!*****************************************************************************80
!
!! P03_STORY prints the "story" for problem p03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
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
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  This example is due to Fred Fritsch and Ralph Carlson.'
  write ( *, '(a)' ) &
    '  This data can cause problems for interpolation methods.'
  write ( *, '(a)' ) &
    '  There are sudden changes in direction, and at the same time,'
  write ( *, '(a)' ) &
    '  sparsely-placed data.  This can cause an interpolant to overshoot'
  write ( *, '(a)' ) &
    '  the data in a way that seems implausible.'

  return
end
subroutine p04_data ( dim_num, data_num, p_data )

!*****************************************************************************80
!
!! P04_DATA returns the data for problem p04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the dependent
!    variables.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Output, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) p_data(dim_num,data_num)

  p_data = reshape ( (/ &
    0.00D+00,  0.00D+00, &
    0.05D+00,  0.70D+00, &
    0.10D+00,  1.00D+00, &
    0.20D+00,  1.00D+00, &
    0.80D+00,  0.30D+00, &
    0.85D+00,  0.05D+00, &
    0.90D+00,  0.10D+00, &
    1.00D+00,  1.00D+00  &
  /), (/ dim_num, data_num /) )

  return
end
subroutine p04_data_num ( data_num )

!*****************************************************************************80
!
!! P04_DATA_NUM returns the number of data points for problem p04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
  implicit none

  integer ( kind = 4 ) data_num

  data_num = 8

  return
end
subroutine p04_dim_num ( dim_num )

!*****************************************************************************80
!
!! P04_DIM_NUM returns the spatial dimension for problem p04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the 
!    dependent variables.
!
  implicit none

  integer ( kind = 4 ) dim_num

  dim_num = 2

  return
end
subroutine p04_story ( )

!*****************************************************************************80
!
!! P04_STORY prints the "story" for problem p04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Larry Irvine, Samuel Marin, Philip Smith,
!    Constrained Interpolation and Smoothing,
!    Constructive Approximation,
!    Volume 2, Number 1, December 1986, pages 129-151.
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  This example is due to Larry Irvine, Samuel Marin and Philip Smith.'
  write ( *, '(a)' ) &
    '  This data can cause problems for interpolation methods.'
  write ( *, '(a)' ) &
    '  There are sudden changes in direction, and at the same time,'
  write ( *, '(a)' ) &
    '  sparsely-placed data.  This can cause an interpolant to overshoot'
  write ( *, '(a)' ) &
    '  the data in a way that seems implausible.'

  return
end
subroutine p05_data ( dim_num, data_num, p_data )

!*****************************************************************************80
!
!! P05_DATA returns the data for problem p05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the dependent
!    variables.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Output, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) p_data(dim_num,data_num)

  p_data = reshape ( (/ &
    0.00D+00, 0.00D+00, &
    0.10D+00, 0.90D+00, &
    0.20D+00, 0.95D+00, &
    0.30D+00, 0.90D+00, &
    0.40D+00, 0.10D+00, &
    0.50D+00, 0.05D+00, &
    0.60D+00, 0.05D+00, &
    0.80D+00, 0.20D+00, &
    1.00D+00, 1.00D+00  &
  /), (/ dim_num, data_num /) )

  return
end
subroutine p05_data_num ( data_num )

!*****************************************************************************80
!
!! P05_DATA_NUM returns the number of data points for problem p05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
  implicit none

  integer ( kind = 4 ) data_num

  data_num = 9

  return
end
subroutine p05_dim_num ( dim_num )

!*****************************************************************************80
!
!! P05_DIM_NUM returns the spatial dimension for problem p05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the 
!    dependent variables.
!
  implicit none

  integer ( kind = 4 ) dim_num

  dim_num = 2

  return
end
subroutine p05_story ( )

!*****************************************************************************80
!
!! P05_STORY prints the "story" for problem p05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Larry Irvine, Samuel Marin, Philip Smith,
!    Constrained Interpolation and Smoothing,
!    Constructive Approximation,
!    Volume 2, Number 1, December 1986, pages 129-151.
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  This example is due to Larry Irvine, Samuel Marin and Philip Smith.'
  write ( *, '(a)' ) &
    '  This data can cause problems for interpolation methods.'
  write ( *, '(a)' ) &
    '  There are sudden changes in direction, and at the same time,'
  write ( *, '(a)' ) &
    '  sparsely-placed data.  This can cause an interpolant to overshoot'
  write ( *, '(a)' ) &
    '  the data in a way that seems implausible.'

  return
end
subroutine p06_data ( dim_num, data_num, p_data )

!*****************************************************************************80
!
!! P06_DATA returns the data for problem p06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the dependent
!    variables.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Output, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) p_data(dim_num,data_num)

  p_data = reshape ( (/ &
   595.0D+00, 0.644D+00, &
   605.0D+00, 0.622D+00, &
   615.0D+00, 0.638D+00, &
   625.0D+00, 0.649D+00, &
   635.0D+00, 0.652D+00, &
   645.0D+00, 0.639D+00, &
   655.0D+00, 0.646D+00, &
   665.0D+00, 0.657D+00, &
   675.0D+00, 0.652D+00, &
   685.0D+00, 0.655D+00, &
   695.0D+00, 0.644D+00, &
   705.0D+00, 0.663D+00, &
   715.0D+00, 0.663D+00, &
   725.0D+00, 0.668D+00, &
   735.0D+00, 0.676D+00, &
   745.0D+00, 0.676D+00, &
   755.0D+00, 0.686D+00, &
   765.0D+00, 0.679D+00, &
   775.0D+00, 0.678D+00, &
   785.0D+00, 0.683D+00, &
   795.0D+00, 0.694D+00, &
   805.0D+00, 0.699D+00, &
   815.0D+00, 0.710D+00, &
   825.0D+00, 0.730D+00, &
   835.0D+00, 0.763D+00, &
   845.0D+00, 0.812D+00, &
   855.0D+00, 0.907D+00, &
   865.0D+00, 1.044D+00, &
   875.0D+00, 1.336D+00, &
   885.0D+00, 1.881D+00, &
   895.0D+00, 2.169D+00, &
   905.0D+00, 2.075D+00, &
   915.0D+00, 1.598D+00, &
   925.0D+00, 1.211D+00, &
   935.0D+00, 0.916D+00, &
   945.0D+00, 0.746D+00, &
   955.0D+00, 0.672D+00, &
   965.0D+00, 0.627D+00, &
   975.0D+00, 0.615D+00, &
   985.0D+00, 0.607D+00, &
   995.0D+00, 0.606D+00, &
  1005.0D+00, 0.609D+00, &
  1015.0D+00, 0.603D+00, &
  1025.0D+00, 0.601D+00, &
  1035.0D+00, 0.603D+00, &
  1045.0D+00, 0.601D+00, &
  1055.0D+00, 0.611D+00, &
  1065.0D+00, 0.601D+00, &
  1075.0D+00, 0.608D+00 /), (/ 2, 49 /) )

  return
end
subroutine p06_data_num ( data_num )

!*****************************************************************************80
!
!! P06_DATA_NUM returns the number of data points for problem p06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
  implicit none

  integer ( kind = 4 ) data_num

  data_num = 49

  return
end
subroutine p06_dim_num ( dim_num )

!*****************************************************************************80
!
!! P06_DIM_NUM returns the spatial dimension for problem p06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the 
!    dependent variables.
!
  implicit none

  integer ( kind = 4 ) dim_num

  dim_num = 2

  return
end
subroutine p06_story ( )

!*****************************************************************************80
!
!! P06_STORY prints the "story" for problem p06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl DeBoor, John Rice,
!    Least-squares cubic spline approximation II - variable knots. 
!    Technical Report CSD TR 21, 
!    Purdue University, Lafayette, Indiana, 1968.
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data is due to deBoor and Rice.'
  write ( *, '(a)' ) &
    '  The data represents a temperature dependent property of titanium.'
  write ( *, '(a)' ) &
    '  The data has been used extensively as an example in spline'
  write ( *, '(a)' ) &
    '  approximation with variably-spaced knots.'
  write ( *, '(a)' ) '  DeBoor considers two sets of knots:'
  write ( *, '(a)' ) '  (595,675,755,835,915,995,1075)'
  write ( *, '(a)' ) '  and'
  write ( *, '(a)' ) '  (595,725,850,910,975,1040,1075).'

  return
end
subroutine p07_data ( dim_num, data_num, p_data )

!*****************************************************************************80
!
!! P07_DATA returns the data for problem p07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the dependent
!    variables.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Output, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) p_data(dim_num,data_num)

  p_data = reshape ( (/ &
   0.0D+00, 1.0D+00, &
   1.0D+00, 2.0D+00, &
   4.0D+00, 2.0D+00, &
   5.0D+00, 1.0D+00 /), (/ 2, 4 /) )

  return
end
subroutine p07_data_num ( data_num )

!*****************************************************************************80
!
!! P07_DATA_NUM returns the number of data points for problem p07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
  implicit none

  integer ( kind = 4 ) data_num

  data_num = 4

  return
end
subroutine p07_dim_num ( dim_num )

!*****************************************************************************80
!
!! P07_DIM_NUM returns the spatial dimension for problem p07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the 
!    dependent variables.
!
  implicit none

  integer ( kind = 4 ) dim_num

  dim_num = 2

  return
end
subroutine p07_story ( )

!*****************************************************************************80
!
!! P07_STORY prints the "story" for problem p07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2012
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This data is a simple symmetric set of 4 points,'
  write ( *, '(a)' ) '  for which it is interesting to develop the Shepard'
  write ( *, '(a)' ) '  interpolants for varying values of the exponent p.'

  return
end
subroutine p08_data ( dim_num, data_num, p_data )

!*****************************************************************************80
!
!! P08_DATA returns the data for problem p07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the dependent
!    variables.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Output, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) p_data(dim_num,data_num)

  p_data = reshape ( (/ &
   -1.0D+00, 1.00D+00, &
   -0.8D+00, 0.64D+00, &
   -0.6D+00, 0.36D+00, &
   -0.4D+00, 0.16D+00, &
   -0.2D+00, 0.04D+00, &
    0.0D+00, 0.00D+00, &
    0.2D+00, 0.04D+00, &
    0.20001D+00, 0.05D+00, &
    0.4D+00, 0.16D+00, &
    0.6D+00, 0.36D+00, &
    0.8D+00, 0.64D+00, &
    1.0D+00, 1.00D+00 /), (/ 2, 12 /) )

  return
end
subroutine p08_data_num ( data_num )

!*****************************************************************************80
!
!! P08_DATA_NUM returns the number of data points for problem p07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
  implicit none

  integer ( kind = 4 ) data_num

  data_num = 12

  return
end
subroutine p08_dim_num ( dim_num )

!*****************************************************************************80
!
!! P08_DIM_NUM returns the spatial dimension for problem p07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the 
!    dependent variables.
!
  implicit none

  integer ( kind = 4 ) dim_num

  dim_num = 2

  return
end
subroutine p08_story ( )

!*****************************************************************************80
!
!! P08_STORY prints the "story" for problem p07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2012
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This is equally spaced data for y = x^2,'
  write ( *, '(a)' ) '  except for one extra point whose x value is'
  write ( *, '(a)' ) '  close to another, but whose y value is not so close.'
  write ( *, '(a)' ) '  A small disagreement in nearby data can be disaster.'

  return
end
