subroutine f00_f0 ( fi, n, x, y, f )

!*****************************************************************************80
!
!! F00_F0 returns the value of any function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FI, the index of the function.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) fi
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  if ( fi == 1 ) then
    call f01_f0 ( n, x, y, f )
  else if ( fi == 2 ) then
    call f02_f0 ( n, x, y, f )
  else if ( fi == 3 ) then
    call f03_f0 ( n, x, y, f )
  else if ( fi == 4 ) then
    call f04_f0 ( n, x, y, f )
  else if ( fi == 5 ) then
    call f05_f0 ( n, x, y, f )
  else if ( fi == 6 ) then
    call f06_f0 ( n, x, y, f )
  else if ( fi == 7 ) then
    call f07_f0 ( n, x, y, f )
  else if ( fi == 8 ) then
    call f08_f0 ( n, x, y, f )
  else if ( fi == 9 ) then
    call f09_f0 ( n, x, y, f )
  else if ( fi == 10 ) then
    call f10_f0 ( n, x, y, f )
  else if ( fi == 11 ) then
    call f11_f0 ( n, x, y, f )
  else if ( fi == 12 ) then
    call f12_f0 ( n, x, y, f )
  else if ( fi == 13 ) then
    call f13_f0 ( n, x, y, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F00_F0 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal function index FI = ', fi
    stop
  end if

  return
end
subroutine f00_f1 ( fi, n, x, y, fx, fy )

!*****************************************************************************80
!
!! F00_F1 returns first derivatives of any function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FI, the index of the function.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the first derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) fi
  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  if ( fi == 1 ) then
    call f01_f1 ( n, x, y, fx, fy )
  else if ( fi == 2 ) then
    call f02_f1 ( n, x, y, fx, fy )
  else if ( fi == 3 ) then
    call f03_f1 ( n, x, y, fx, fy )
  else if ( fi == 4 ) then
    call f04_f1 ( n, x, y, fx, fy )
  else if ( fi == 5 ) then
    call f05_f1 ( n, x, y, fx, fy )
  else if ( fi == 6 ) then
    call f06_f1 ( n, x, y, fx, fy )
  else if ( fi == 7 ) then
    call f07_f1 ( n, x, y, fx, fy )
  else if ( fi == 8 ) then
    call f08_f1 ( n, x, y, fx, fy )
  else if ( fi == 9 ) then
    call f09_f1 ( n, x, y, fx, fy )
  else if ( fi == 10 ) then
    call f10_f1 ( n, x, y, fx, fy )
  else if ( fi == 11 ) then
    call f11_f1 ( n, x, y, fx, fy )
  else if ( fi == 12 ) then
    call f12_f1 ( n, x, y, fx, fy )
  else if ( fi == 13 ) then
    call f13_f1 ( n, x, y, fx, fy )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F00_F1 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal function index FI = ', fi
    stop
  end if

  return
end
subroutine f00_f2 ( fi, n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F00_F2 returns second derivatives of any function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FI, the index of the function.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), the second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) fi
  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  if ( fi == 1 ) then
    call f01_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 2 ) then
    call f02_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 3 ) then
    call f03_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 4 ) then
    call f04_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 5 ) then
    call f05_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 6 ) then
    call f06_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 7 ) then
    call f07_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 8 ) then
    call f08_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 9 ) then
    call f09_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 10 ) then
    call f10_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 11 ) then
    call f11_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 12 ) then
    call f12_f2 ( n, x, y, fxx, fxy, fyy )
  else if ( fi == 13 ) then
    call f13_f2 ( n, x, y, fxx, fxy, fyy )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F00_F2 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal function index FI = ', fi
    stop
  end if

  return
end
subroutine f00_num ( fun_num )

!*****************************************************************************80
!
!! F00_NUM returns the number of test functions available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!   Output, integer ( kind = 4 ) FUN_NUM, the number of test functions.
!
  implicit none

  integer ( kind = 4 ) fun_num

  fun_num = 13

  return
end
subroutine f00_title ( fi, ft )

!*****************************************************************************80
!
!! F00_TITLE returns the title for any function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FI, the index of the function.
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  integer ( kind = 4 ) fi
  character ( len = * ) ft

  if ( fi == 1 ) then
    call f01_title ( ft )
  else if ( fi == 2 ) then
    call f02_title ( ft )
  else if ( fi == 3 ) then
    call f03_title ( ft )
  else if ( fi == 4 ) then
    call f04_title ( ft )
  else if ( fi == 5 ) then
    call f05_title ( ft )
  else if ( fi == 6 ) then
    call f06_title ( ft )
  else if ( fi == 7 ) then
    call f07_title ( ft )
  else if ( fi == 8 ) then
    call f08_title ( ft )
  else if ( fi == 9 ) then
    call f09_title ( ft )
  else if ( fi == 10 ) then
    call f10_title ( ft )
  else if ( fi == 11 ) then
    call f11_title ( ft )
  else if ( fi == 12 ) then
    call f12_title ( ft )
  else if ( fi == 13 ) then
    call f13_title ( ft )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal function index FI = ', fi
    stop
  end if

  return
end
subroutine f01_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F01_F0 returns the value of function 1.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  f(1:n) = &
      0.75 * exp ( - ( ( 9.0 * x(1:n) - 2.0 )**2            &
                     + ( 9.0 * y(1:n) - 2.0 )**2 ) / 4.0 )  &
    + 0.75 * exp ( - ( ( 9.0 * x(1:n) + 1.0 )**2 ) / 49.0   &
                     - ( 9.0 * y(1:n) + 1.0 ) / 10.0 )      &
    + 0.5  * exp ( - ( ( 9.0 * x(1:n) - 7.0 )**2            &
                     + ( 9.0 * y(1:n) - 3.0 )**2 ) / 4.0 )  &
    - 0.2  * exp ( - (   9.0 * x(1:n) - 4.0 )**2            &
                     - ( 9.0 * y(1:n) - 7.0 )**2 )

  return
end
subroutine f01_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F01_F1 returns first derivatives of function 1.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = exp ( - ( ( 9.0 * x(1:n) - 2.0 )**2 &
                    + ( 9.0 * y(1:n) - 2.0 )**2 ) / 4.0 )
  t2(1:n) = exp ( - ( ( 9.0 * x(1:n) + 1.0 )**2 ) / 49.0 &
                    - ( 9.0 * y(1:n) + 1.0 ) / 10.0 )
  t3(1:n) = exp ( - ( ( 9.0 * x(1:n) - 7.0 )**2 &
                    + ( 9.0 * y(1:n) - 3.0 )**2 ) / 4.0 )
  t4(1:n) = exp ( -   ( 9.0 * x(1:n) - 4.0 )**2 &
                    - ( 9.0 * y(1:n) - 7.0 )**2 )

  fx(1:n) = &
    - 3.375           * ( 9.0 * x(1:n) - 2.0 ) * t1 &
    - ( 27.0 / 98.0 ) * ( 9.0 * x(1:n) + 1.0 ) * t2 &
    - 2.25            * ( 9.0 * x(1:n) - 7.0 ) * t3 &
    + 3.6             * ( 9.0 * x(1:n) - 4.0 ) * t4

  fy(1:n) = &
    - 3.375 * ( 9.0 * y(1:n) - 2.0 ) * t1 &
    - 0.675                          * t2 &
    - 2.25  * ( 9.0 * y(1:n) - 3.0 ) * t3 &
    + 3.6   * ( 9.0 * y(1:n) - 7.0 ) * t4

  return
end
subroutine f01_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F01_F2 returns second derivatives of function 1.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = exp ( - ( ( 9.0 * x(1:n) - 2.0 )**2 &
                    + ( 9.0 * y(1:n) - 2.0 )**2 ) / 4.0 )
  t2(1:n) = exp ( - ( ( 9.0 * x(1:n) + 1.0 )**2 ) / 49.0 &
                    - ( 9.0 * y(1:n) + 1.0 ) / 10.0 )
  t3(1:n) = exp ( - ( ( 9.0 * x(1:n) - 7.0 )**2 &
                    + ( 9.0 * y(1:n) - 3.0 )**2 ) / 4.0 )
  t4(1:n) = exp ( -   ( 9.0 * x(1:n) - 4.0 )**2 &
                    - ( 9.0 * y(1:n) - 7.0 )**2 )

  fxx(1:n) = &
      15.1875 * ( ( 9.0 * x(1:n) - 2.0 )**2 - 2.0 )  * t1(1:n) &
    + 60.75   * ( ( 9.0 * x(1:n) + 1.0 )**2 - 24.5 ) * t2(1:n) &
    + 10.125  * ( ( 9.0 * x(1:n) - 7.0 )**2 - 2.0 )  * t3(1:n) &
    - 64.8    * ( ( 9.0 * x(1:n) - 4.0 )**2 - 0.5 )  * t4(1:n)

  fxy(1:n) = &
      15.1875 * ( 9.0 * x(1:n) - 2.0 ) * ( 9.0 * y(1:n) - 2.0 ) * t1(1:n) &
    + ( 243.0 / 980.0 ) * ( 9.0 * x(1:n) + 1.0 ) * t2(1:n) &
    + 10.125 * ( 9.0 * x(1:n) - 7.0 ) * ( 9.0 * y(1:n) - 3.0 ) * t3(1:n) &
    - 64.8 * ( 9.0 * x(1:n) - 4.0 ) * ( 9.0 * y(1:n) - 7.0 ) * t4(1:n)

  fyy(1:n) = &
      15.1875 * ( ( 9.0 * y(1:n) - 2.0 )**2 - 2.0 ) * t1(1:n) &
    + 0.6075  *                                       t2(1:n) &
    + 10.125  * ( ( 9.0 * y(1:n) - 3.0 )**2 - 2.0 ) * t3(1:n) &
    - 64.8    * ( ( 9.0 * y(1:n) - 7.0 )**2 - 0.5 ) * t4(1:n)

  return
end
subroutine f01_title ( ft )

!*****************************************************************************80
!
!! F01_TITLE returns the title for function 1.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Exponential'

  return
end
subroutine f02_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F02_F0 returns the value of function 2.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  f(1:n) = ( tanh ( 9.0 * ( y(1:n) - x(1:n) ) ) + 1.0 ) / 9.0

  return
end
subroutine f02_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F02_F1 returns first derivatives of function 2.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = 18.0 * ( y(1:n) - x(1:n) )
  fx(1:n) = - 4.0 / ( exp ( t1(1:n) ) + 2.0 + exp ( - t1(1:n) ) )
  fy(1:n) = - fx(1:n)

  return
end
subroutine f02_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F02_F2 returns second derivatives of function 2.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = 18.0 * ( y(1:n) - x(1:n) )

  fxx(1:n) = 18.0 * tanh ( 0.5 * t1(1:n) ) &
    * ( tanh ( 9.0 * ( y(1:n) - x(1:n) ) ) + 1.0 ) / 9.0
  fxy(1:n) = - fxx(1:n)
  fyy(1:n) = fxx(1:n)

  return
end
subroutine f02_title ( ft )

!*****************************************************************************80
!
!! F02_TITLE returns the title for function 2.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Cliff'

  return
end
subroutine f03_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F03_F0 returns the value of function 3.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  f(1:n) = ( 1.25 + cos ( 5.4 * y(1:n) ) ) &
    / ( 6.0 + 6.0 * ( 3.0 * x(1:n) - 1.0 )**2 )

  return
end
subroutine f03_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F03_F1 returns first derivatives of function 3.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = 5.4 * y(1:n)
  t2(1:n) = 1.0 + ( 3.0 * x(1:n) - 1.0 )**2
  fx(1:n) = - ( 3.0 * x(1:n) - 1.0 ) &
    * ( 1.25 + cos ( t1(1:n) ) ) / ( t2(1:n) * t2(1:n) )
  fy(1:n) = - 0.9 * sin ( t1(1:n) ) / t2(1:n)

  return
end
subroutine f03_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F03_F2 returns second derivatives of function 3.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = 5.4 * y(1:n)
  t2(1:n) = 1.0 + ( 3.0 * x(1:n) - 1.0 )**2

  fxx(1:n) = 3.0 * ( 1.25 + cos ( t1(1:n) ) ) * ( 3.0 * t2(1:n) - 4.0 ) &
    / ( t2(1:n)**3 )
  fxy(1:n) = 5.4 * ( 3.0 * x(1:n) - 1.0 ) * sin ( t1(1:n) ) &
    / ( t2(1:n) * t2(1:n) )
  fyy(1:n) = - 4.86 * cos ( t1(1:n) ) / t2(1:n)

  return
end
subroutine f03_title ( ft )

!*****************************************************************************80
!
!! F03_TITLE returns the title for function 3.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Saddle'

  return
end
subroutine f04_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F04_F0 returns the value of function 4.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  f(1:n) = exp ( - 5.0625 * ( ( x(1:n) - 0.5 )**2 &
                            + ( y(1:n) - 0.5 )**2 ) ) / 3.0

  return
end
subroutine f04_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F04_F1 returns first derivatives of function 4.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = x(1:n) - 0.5
  t2(1:n) = y(1:n) - 0.5
  t3(1:n) = - 3.375 &
    * exp ( - 5.0625 * ( t1(1:n) * t1(1:n) + t2(1:n) * t2(1:n) ) )
  fx(1:n) = t1(1:n) * t3(1:n)
  fy(1:n) = t2(1:n) * t3(1:n)

  return
end
subroutine f04_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F04_F2 returns second derivatives of function 4.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = x(1:n) - 0.5
  t2(1:n) = y(1:n) - 0.5
  t3(1:n) = - 3.375 &
    * exp ( - 5.0625 * ( t1(1:n) * t1(1:n) + t2(1:n) * t2(1:n) ) )

  fxx(1:n) = ( 1.0 - 10.125 * t1(1:n) * t1(1:n) ) * t3(1:n)
  fxy(1:n) = - 10.125 * t1(1:n) * t2(1:n) * t3(1:n)
  fyy(1:n) = ( 1.0 - 10.125 * t2(1:n) * t2(1:n) ) * t3(1:n)

  return
end
subroutine f04_title ( ft )

!*****************************************************************************80
!
!! F04_TITLE returns the title for function 4.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Gentle'

  return
end
subroutine f05_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F05_F0 returns the value of function 5.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  f(1:n) = exp ( - 20.25 * ( ( x(1:n) - 0.5 )**2 + ( y(1:n) - 0.5 )**2 ) ) / 3.0

  return
end
subroutine f05_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F05_F1 returns first derivatives of function 5.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = x(1:n) - 0.5
  t2(1:n) = y(1:n) - 0.5
  t3(1:n) = - 13.5 * exp ( - 20.25 * ( t1(1:n) * t1(1:n) + t2(1:n) * t2(1:n) ) )
  fx(1:n) = t1(1:n) * t3(1:n)
  fy(1:n) = t2(1:n) * t3(1:n)

  return
end
subroutine f05_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F05_F2 returns second derivatives of function 5.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = x(1:n) - 0.5
  t2(1:n) = y(1:n) - 0.5
  t3(1:n) = - 13.5 * exp ( - 20.25 * ( t1(1:n) * t1(1:n) + t2(1:n) * t2(1:n) ) )

  fxx(1:n) = ( 1.0 - 40.5 * t1(1:n) * t1(1:n) ) * t3(1:n)
  fxy(1:n) = - 40.5 * t1(1:n) * t2(1:n) * t3(1:n)
  fyy(1:n) = ( 1.0 - 40.5 * t2(1:n) * t2(1:n) ) * t3(1:n)

  return
end
subroutine f05_title ( ft )

!*****************************************************************************80
!
!! F05_TITLE returns the title for function 5.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Steep'

  return
end
subroutine f06_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F06_F0 returns the value of function 6.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t4(1:n) = 64.0 - 81.0 * ( ( x(1:n) - 0.5 )**2 + ( y(1:n) - 0.5 )**2 )

  do i = 1, n

    if ( 0.0 <= t4(i) ) then
      f(i) = sqrt ( t4(i) ) / 9.0 - 0.5
    else
      f(i) = 0.0
    end if

  end do

  return
end
subroutine f06_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F06_F1 returns first derivatives of function 6.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t4(1:n) = 64.0 - 81.0 * ( ( x(1:n) - 0.5 )**2 + ( y(1:n) - 0.5 )**2 )

  do i = 1, n

    if ( 0.0 < t4(i) ) then
      t1 = x(i) - 0.5
      t2 = y(i) - 0.5
      t3 = - 9.0 / sqrt ( t4(i) )
      fx(i) = t1 * t3
      fy(i) = t2 * t3
    else
      fx(i) = 0.0
      fy(i) = 0.0
    end if

  end do

  return
end
subroutine f06_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F06_F2 returns second derivatives of function 6.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t4(1:n) = 64.0 - 81.0 * ( ( x(1:n) - 0.5 )**2 + ( y(1:n) - 0.5 )**2 )

  do i = 1, n

    if ( 0.0 < t4(i) ) then
      t1 = x(i) - 0.5
      t2 = y(i) - 0.5
      t3 = - 9.0 / sqrt ( t4(i) )
      fxx(i) = ( 1.0 + t1 * t3 * t1 * t3 ) * t3
      fxy(i) = t1 * t3 * t2 * t3
      fyy(i) = ( 1.0 + t2 * t3 * t2 * t3 ) * t3
    else
      fxx(i) = 0.0
      fxy(i) = 0.0
      fyy(i) = 0.0
    end if

  end do

  return
end
subroutine f06_title ( ft )

!*****************************************************************************80
!
!! F06_TITLE returns the title for function 6.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Sphere'

  return
end
subroutine f07_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F07_F0 returns the value of function 7.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  f(1:n) = 2.0 * cos ( 10.0 * x(1:n) ) * sin ( 10.0 * y(1:n) ) &
    + sin ( 10.0 * x(1:n) * y(1:n) )

  return
end
subroutine f07_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F07_F1 returns first derivatives of function 7.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = 10.0 * x(1:n)
  t2(1:n) = 10.0 * y(1:n)
  t3(1:n) = 10.0 * cos ( 10.0 * x(1:n) * y(1:n) )
  fx(1:n) = - 20.0 * sin ( t1(1:n) ) * sin ( t2(1:n) ) + t3(1:n) * y(1:n)
  fy(1:n) =   20.0 * cos ( t1(1:n) ) * cos ( t2(1:n) ) + t3(1:n) * x(1:n)

  return
end
subroutine f07_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F07_F2 returns second derivatives of function 7.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = 10.0 * x(1:n)
  t2(1:n) = 10.0 * y(1:n)
  t3(1:n) = 10.0 * cos ( 10.0 * x(1:n) * y(1:n) )
  t4(1:n) = 100.0 * sin ( 10.0 * x(1:n) * y(1:n) )

  fxx(1:n) = - 200.0 * cos ( t1(1:n) ) * sin ( t2(1:n) ) &
    - t4(1:n) * y(1:n) * y(1:n)
  fxy(1:n) = - 200.0 * sin ( t1(1:n) ) * cos ( t2(1:n) ) &
    + t3(1:n) - t4(1:n) * x(1:n) * y(1:n)
  fyy(1:n) = - 200.0 * cos ( t1(1:n) ) * sin ( t2(1:n) ) &
    - t4(1:n) * x(1:n) * x(1:n)

  return
end
subroutine f07_title ( ft )

!*****************************************************************************80
!
!! F07_TITLE returns the title for function 7.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Trig'

  return
end
subroutine f08_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F08_F0 returns the value of function 8.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = 5.0 - 10.0 * x(1:n)
  t2(1:n) = 5.0 - 10.0 * y(1:n)
  t3(1:n) = exp ( - 0.5 * t1(1:n) * t1(1:n) )
  t4(1:n) = exp ( - 0.5 * t2(1:n) * t2(1:n) )
  f(1:n) = t3(1:n) + 0.75 * t4(1:n) * ( 1.0 + t3(1:n) )

  return
end
subroutine f08_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F08_F1 returns first derivatives of function 8.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = 5.0 - 10.0 * x(1:n)
  t2(1:n) = 5.0 - 10.0 * y(1:n)
  t3(1:n) = exp ( - 0.5 * t1(1:n) * t1(1:n) )
  t4(1:n) = exp ( - 0.5 * t2(1:n) * t2(1:n) )

  fx(1:n) = t1(1:n) * t3(1:n) * ( 10.0 + 7.5 * t4(1:n) )
  fy(1:n) = t2(1:n) * t4(1:n) * ( 7.5 + 7.5 * t3(1:n) )

  return
end
subroutine f08_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F08_F2 returns second derivatives of function 8.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = 5.0 - 10.0 * x(1:n)
  t2(1:n) = 5.0 - 10.0 * y(1:n)
  t3(1:n) = exp ( - 0.5 * t1(1:n) * t1(1:n) )
  t4(1:n) = exp ( - 0.5 * t2(1:n) * t2(1:n) )

  fxx(1:n) = t3(1:n) * ( t1(1:n) * t1(1:n) - 1.0 ) * ( 100.0 + 75.0 * t4(1:n) )
  fxy(1:n) = 75.0 * t1(1:n) * t2(1:n) * t3(1:n) * t4(1:n)
  fyy(1:n) = t4(1:n) * ( t2(1:n) * t2(1:n) - 1.0 ) * ( 75.0 + 75.0 * t3(1:n) )

  return
end
subroutine f08_title ( ft )

!*****************************************************************************80
!
!! F08_TITLE returns the title for function 8.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Synergistic Gaussian'

  return
end
subroutine f09_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F09_F0 returns the value of function 9.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = exp ( ( 10.0 - 20.0 * x(1:n) ) / 3.0 )
  t2(1:n) = exp ( ( 10.0 - 20.0 * y(1:n) ) / 3.0 )
  t3(1:n) = 1.0 / ( 1.0 + t1(1:n) )
  t4(1:n) = 1.0 / ( 1.0 + t2(1:n) )
  f(1:n) = ( ( 20.0 / 3.0 )**3 * t1(1:n) * t2(1:n) )**2 &
    * ( t3(1:n) * t4(1:n) )**5 &
    * ( t1(1:n) - 2.0 * t3(1:n) ) * ( t2(1:n) - 2.0 * t4(1:n) )

  return
end
subroutine f09_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F09_F1 returns first derivatives of function 9.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = exp ( ( 10.0 - 20.0 * x(1:n) ) / 3.0 )
  t2(1:n) = exp ( ( 10.0 - 20.0 * y(1:n) ) / 3.0 )
  t3(1:n) = 1.0 / ( 1.0 + t1(1:n) )
  t4(1:n) = 1.0 / ( 1.0 + t2(1:n) )

  fx(1:n) = ( ( 20.0 / 3.0 ) * t1(1:n) )**2 * ( ( 20.0 / 3.0 ) * t3(1:n) )**5 &
    * ( 2.0 * t1(1:n) - 3.0 * t3(1:n) - 5.0 + 12.0 * t3(1:n) * t3(1:n) ) &
    * t2(1:n) * t2(1:n) * t4(1:n)**5 &
    * ( t2(1:n) - 2.0 * t4(1:n) )

  fy(1:n) = ( ( 20.0 / 3.0 ) * t1(1:n) )**2 * ( ( 20.0 / 3.0 ) * t3(1:n) )**5 &
    * ( 2.0 * t2(1:n) - 3.0 * t4(1:n) - 5.0 + 12.0 * t4(1:n) * t4(1:n) ) &
    * t2(1:n) * t2(1:n) * t4(1:n)**5 &
    * ( t1(1:n) - 2.0 * t3(1:n) )

  return
end
subroutine f09_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F09_F2 returns second derivatives of function 9.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) t4(n)
  real ( kind = 8 ) t5
  real ( kind = 8 ) t6(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = exp ( ( 10.0 - 20.0 * x(1:n) ) / 3.0 )
  t2(1:n) = exp ( ( 10.0 - 20.0 * y(1:n) ) / 3.0 )
  t3(1:n) = 1.0 / ( 1.0 + t1(1:n) )
  t4(1:n) = 1.0 / ( 1.0 + t2(1:n) )
  t5 = 20.0 / 3.0
  t6(1:n) = ( t5 * t1(1:n) * t2(1:n) )**2 * ( t5 * t3(1:n) * t4(1:n) )**5

  fxx(1:n) = t5 * t6(1:n) * ( t2(1:n) - 2.0 * t4(1:n) ) &
    * ( ( ( - 84.0 * t3(1:n) + 78.0 ) * t3(1:n) + 23.0 ) * t3(1:n) &
    + 4.0 * t1(1:n) - 25.0 )

  fxy(1:n) = t5 * t6(1:n) &
    * ( ( 12.0 * t4(1:n) - 3.0 ) * t4(1:n) + 2.0 * t2(1:n) - 5.0 ) &
    * ( ( 12.0 * t3(1:n) - 3.0 ) * t3(1:n) + 2.0 * t1(1:n) - 5.0 )

  fyy(1:n) = t5 * t6(1:n) * ( t1(1:n) - 2.0 * t3(1:n) ) &
    * ( ( ( - 84.0 * t4(1:n) + 78.0 ) * t4(1:n) + 23.0 ) * t4(1:n) &
    + 4.0 * t2(1:n) - 25.0 )

  return
end
subroutine f09_title ( ft )

!*****************************************************************************80
!
!! F09_TITLE returns the title for function 9.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Cloverleaf Asymmetric Peak/Valley'

  return
end
subroutine f10_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F10_F0 returns the value of function f10.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2(n)
  real ( kind = 8 ) t3(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = sqrt ( ( 80.0 * x(1:n) - 40.0 )**2 + ( 90.0 * y(1:n) - 45.0 )**2 )
  t2(1:n) = exp ( - 0.04 * t1(1:n) )
  t3(1:n) = cos ( 0.15 * t1(1:n) )

  f(1:n) = t2(1:n) * t3(1:n)

  return
end
subroutine f10_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F10_F1 returns first derivatives of function f10.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = sqrt ( ( 80.0 * x(1:n) - 40.0 )**2 + ( 90.0 * y(1:n) - 45.0 )**2 )

  do i = 1, n

    if ( t1(i) == 0.0 ) then
      fx(i) = 0.0
      fy(i) = 0.0
    else
      t2 = exp ( - 0.04 * t1(i) )
      t3 = cos ( 0.15 * t1(i) )
      t4 = sin ( 0.15 * t1(i) )
      fx(i) = - t2 * ( 12.0 * t4 + 3.2 * t3 ) * ( 80.0 * x(i) - 40.0 ) / t1(i)
      fy(i) = - t2 * ( 13.5 * t4 + 3.6 * t3 ) * ( 90.0 * y(i) - 45.0 ) / t1(i)
    end if

  end do

  return
end
subroutine f10_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F10_F2 returns second derivatives of function f10.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t1(n)
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) t5
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  t1(1:n) = sqrt ( ( 80.0 * x(1:n) - 40.0 )**2 + ( 90.0 * y(1:n) - 45.0 )**2 )

  do i = 1, n

    if ( t1(i) == 0.0 ) then
      fxx(i) = 0.0
      fxy(i) = 0.0
      fyy(i) = 0.0
    else
      t2 = exp ( - 0.04 * t1(i) )
      t3 = cos ( 0.15 * t1(i) )
      t4 = sin ( 0.15 * t1(i) )
      t5 = t2 / t1(i)**3

      fxx(i) = t5 * ( t1(i) * ( 76.8 * t4 - 133.76 * t3 ) &
        * ( 80.0 * x(i) - 40.0 )**2 &
        - ( 960.0 * t4 + 256.0 * t3 ) * ( 90.0 * y(i) - 45.0 )**2 )

      fxy(i) = t5 * ( t1(i) * ( 86.4 * t4 - 150.48 * t3 ) + 1080.0 * t4 + &
        288.0 * t3 ) * ( 80.0 * x(i) - 40.0 ) * ( 90.0 * y(i) - 45.0 )

      fyy(i) = t5 * ( t1(i) * ( 97.2 * t4 - 169.29 * t3 ) &
        * ( 90.0 * y(i) - 45.0 )**2 &
        - ( 1215.0 * t4 + 324.0 * t3 ) * ( 80.0 * x(i) - 40.0 )**2 )

    end if

  end do

  return
end
subroutine f10_title ( ft )

!*****************************************************************************80
!
!! F10_TITLE returns the title for function f10.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Cosine Peak'

  return
end
subroutine f11_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F11_F0 returns the value of function f11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  f(1:n) = x(1:n) * ( y(1:n) + 1.0D+00 )

  return
end
subroutine f11_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F11_F1 returns first derivatives of function f11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  fx(1:n) = y(1:n) + 1.0D+00
  fy(1:n) = x(1:n)

  return
end
subroutine f11_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F11_F2 returns second derivatives of function f11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  fxx(1:n) = 0.0D+00
  fxy(1:n) = 1.0D+00
  fyy(1:n) = 0.0D+00

  return
end
subroutine f11_title ( ft )

!*****************************************************************************80
!
!! F11_TITLE returns the title for function f11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Bilinear function'

  return
end
subroutine f12_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F12_F0 returns the value of function f12.
!
!  Discussion:
!
!    This is an example from Vicente Romero.
!
!      R = sqrt ( X^2 + Y^2 )
!      T = atan ( Y / X )
!      F(X,Y) = ( 0.8 * R + 0.35 * sin ( 2.4 * pi * R / sqrt ( 2 )  ) )
!               * 1.5 * sin ( 1.3 * T )
!
!    The mean and standard deviation of the function over the interval
!    are approximately:
!
!      mu    = 0.581608
!      sigma = 0.343208
!
!    Since the interpolation interval is the unit square, this means the
!    integral of the function over the interval can also be estimated as
!
!      I = 0.581608
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Vicente Romero, John Burkardt, Max Gunzburger, Janet Peterson,
!    Initial Evaluation of Centroidal Voronoi Tessellation for
!    Statistical Sampling and Function Integration,
!    Fourth International Symposium on Uncertainty Modeling and Analysis,
!    (ISUMA) 2003.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  r(1:n) = sqrt ( x(1:n)**2 + y(1:n)**2 )
  t(1:n) = atan2 ( y(1:n), x(1:n) )

  f(1:n) = 1.5D+00 * ( 0.8D+00 * r(1:n) &
    + 0.35D+00 * sin ( 2.4D+00 * pi * r(1:n) / sqrt ( 2.0D+00 ) ) ) &
    * sin ( 1.3D+00 * t(1:n) )

  return
end
subroutine f12_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F12_F1 returns first derivatives of function f12.
!
!  Discussion:
!
!    Currently, the derivative information is of no interest to me.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  fx(1:n) = 0.0D+00
  fy(1:n) = 0.0D+00

  return
end
subroutine f12_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F12_F2 returns second derivatives of function f12.
!
!  Discussion:
!
!    Currently, the derivative information is of no interest to me.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  fxx(1:n) = 0.0D+00
  fxy(1:n) = 0.0D+00
  fyy(1:n) = 0.0D+00

  return
end
subroutine f12_title ( ft )

!*****************************************************************************80
!
!! F12_TITLE returns the title for function f12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Vicente Romero function'

  return
end
subroutine f13_f0 ( n, x, y, f )

!*****************************************************************************80
!
!! F13_F0 returns the value of function f13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  f(1:n) = 1.0D+00 / ( ( 10.0D+00 * x(1:n) - 5.0D+00 )**2 &
                     + ( 10.0D+00 * y(1:n) - 5.0D+00 )**2 + 1.0D+00 )

  return
end
subroutine f13_f1 ( n, x, y, fx, fy )

!*****************************************************************************80
!
!! F13_F1 returns first derivatives of function f13.
!
!  Discussion:
!
!    Currently, the derivative information is of no interest to me.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FX(N), FY(N), the derivative values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) fy(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)  

  fx(1:n) = 0.0D+00
  fy(1:n) = 0.0D+00

  return
end
subroutine f13_f2 ( n, x, y, fxx, fxy, fyy )

!*****************************************************************************80
!
!! F13_F2 returns second derivatives of function f13.
!
!  Discussion:
!
!    Currently, the derivative information is of no interest to me.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evalution points.
!
!    Output, real ( kind = 8 ) FXX(N), FXY(N), FYY(N), second derivatives.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fxx(n)
  real ( kind = 8 ) fxy(n)
  real ( kind = 8 ) fyy(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  fxx(1:n) = 0.0D+00
  fxy(1:n) = 0.0D+00
  fyy(1:n) = 0.0D+00

  return
end
subroutine f13_title ( ft )

!*****************************************************************************80
!
!! F13_TITLE returns the title for function f13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FT, the function title.
!
  implicit none

  character ( len = * ) ft

  ft = 'Rescaled Runge function'

  return
end
subroutine g00_num ( grid_num )

!*****************************************************************************80
!
!! G00_NUM returns the number of grids available.
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
!    John Burkardt
!
!  Parameters:
!
!   Output, integer ( kind = 4 ) GRID_NUM, the number of grids.
!
  implicit none

  integer ( kind = 4 ) grid_num

  grid_num = 5

  return
end
subroutine g00_size ( gi, gn )

!*****************************************************************************80
!
!! G00_SIZE returns the size for any grid.
!
!  Discussion:
!
!    The "grid size" is simply the number of points in the grid.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GI, the index of the grid.
!
!    Output, integer ( kind = 4 ) GN, the grid size.
!
  implicit none

  integer ( kind = 4 ) gi
  integer ( kind = 4 ) gn

  if ( gi == 1 ) then
    call g01_size ( gn )
  else if ( gi == 2 ) then
    call g02_size ( gn )
  else if ( gi == 3 ) then
    call g03_size ( gn )
  else if ( gi == 4 ) then
    call g04_size ( gn )
  else if ( gi == 5 ) then
    call g05_size ( gn )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G00_SIZE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal grid index GI = ', gi
    stop
  end if

  return
end
subroutine g00_title ( gi, gt )

!*****************************************************************************80
!
!! G00_TITLE returns the title for any grid.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GI, the index of the grid.
!
!    Output, character ( len = * ) GT, the grid title.
!
  implicit none

  integer ( kind = 4 ) gi
  character ( len = * ) gt

  if ( gi == 1 ) then
    call g01_title ( gt )
  else if ( gi == 2 ) then
    call g02_title ( gt )
  else if ( gi == 3 ) then
    call g03_title ( gt )
  else if ( gi == 4 ) then
    call g04_title ( gt )
  else if ( gi == 5 ) then
    call g05_title ( gt )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal grid index GI = ', gi
    stop
  end if

  return
end
subroutine g00_xy ( gi, gn, gx, gy )

!*****************************************************************************80
!
!! G00_XY returns the grid points for any grid.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GI, the index of the grid.
!
!    Input, integer ( kind = 4 ) GN, the grid size.
!
!    Output, real ( kind = 8 ) GX(GN), GY(GN), the grid coordinates.
!
  implicit none

  integer ( kind = 4 ) gn

  integer ( kind = 4 ) gi
  real ( kind = 8 ) gx(gn)
  real ( kind = 8 ) gy(gn)

  if ( gi == 1 ) then
    call g01_xy ( gn, gx, gy )
  else if ( gi == 2 ) then
    call g02_xy ( gn, gx, gy )
  else if ( gi == 3 ) then
    call g03_xy ( gn, gx, gy )
  else if ( gi == 4 ) then
    call g04_xy ( gn, gx, gy )
  else if ( gi == 5 ) then
    call g05_xy ( gn, gx, gy )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G00_XY - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal grid index GI = ', gi
    stop
  end if

  return
end
subroutine g01_size ( gn )

!*****************************************************************************80
!
!! G01_SIZE returns the size for grid 1.
!
!  Discussion:
!
!    The "grid size" is simply the number of points in the grid.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) GN, the grid size.
!
  implicit none

  integer ( kind = 4 ) gn

  gn = 100

  return
end
subroutine g01_title ( gt )

!*****************************************************************************80
!
!! G01_TITLE returns the title for grid 1.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) GT, the grid title.
!
  implicit none

  character ( len = * ) gt

  gt = 'Franke''s 100 node set'

  return
end
subroutine g01_xy ( gn, gx, gy )

!*****************************************************************************80
!
!! G01_XY returns the grid points for grid 1.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GN, the grid size.
!
!    Output, real ( kind = 8 ) GX(GN), GY(GN), the grid coordinates.
!
  implicit none

  integer ( kind = 4 ) gn

  real ( kind = 8 ) gx(gn)
  real ( kind = 8 ) gy(gn)
  real ( kind = 8 ), dimension ( 100 ) :: x = (/ &
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
  real ( kind = 8 ), dimension ( 100 ) :: y = (/ &
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

  gx(1:gn) = x(1:gn)
  gy(1:gn) = y(1:gn)

  return
end
subroutine g02_size ( gn )

!*****************************************************************************80
!
!! G02_SIZE returns the size for grid 2.
!
!  Discussion:
!
!    The "grid size" is simply the number of points in the grid.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) GN, the grid size.
!
  implicit none

  integer ( kind = 4 ) gn

  gn = 33

  return
end
subroutine g02_title ( gt )

!*****************************************************************************80
!
!! G02_TITLE returns the title for grid 2.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) GT, the grid title.
!
  implicit none

  character ( len = * ) gt

  gt = 'Franke''s 33 node set'

  return
end
subroutine g02_xy ( gn, gx, gy )

!*****************************************************************************80
!
!! G02_XY returns the grid points for grid 2.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GN, the grid size.
!
!    Output, real ( kind = 8 ) GX(GN), GY(GN), the grid coordinates.
!
  implicit none

  integer ( kind = 4 ) gn

  real ( kind = 8 ) gx(gn)
  real ( kind = 8 ) gy(gn)
  real ( kind = 8 ), dimension ( 33 ) :: x = (/ &
    0.05,  0.00,  0.00,  0.00,  0.10, &
    0.10,  0.15,  0.20,  0.25,  0.30, &
    0.35,  0.50,  0.50,  0.55,  0.60, &
    0.60,  0.60,  0.65,  0.70,  0.70, &
    0.70,  0.75,  0.75,  0.75,  0.80, &
    0.80,  0.85,  0.90,  0.90,  0.95, &
    1.00,  1.00,  1.00 /)
  real ( kind = 8 ), dimension ( 33 ) :: y = (/ &
    0.45,  0.50,  1.00,  0.00,  0.15, &
    0.75,  0.30,  0.10,  0.20,  0.35, &
    0.85,  0.00,  1.00,  0.95,  0.25, &
    0.65,  0.85,  0.70,  0.20,  0.65, &
    0.90,  0.10,  0.35,  0.85,  0.40, &
    0.65,  0.25,  0.35,  0.80,  0.90, &
    0.00,  0.50,  1.00 /)

  gx(1:gn) = x(1:gn)
  gy(1:gn) = y(1:gn)

  return
end
subroutine g03_size ( gn )

!*****************************************************************************80
!
!! G03_SIZE returns the size for grid 3.
!
!  Discussion:
!
!    The "grid size" is simply the number of points in the grid.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) GN, the grid size.
!
  implicit none

  integer ( kind = 4 ) gn

  gn = 25

  return
end
subroutine g03_title ( gt )

!*****************************************************************************80
!
!! G03_TITLE returns the title for grid 3.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) GT, the grid title.
!
  implicit none

  character ( len = * ) gt

  gt = 'Lawson''s 25 node set'

  return
end
subroutine g03_xy ( gn, gx, gy )

!*****************************************************************************80
!
!! G03_XY returns the grid points for grid 3.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GN, the grid size.
!
!    Output, real ( kind = 8 ) GX(GN), GY(GN), the grid coordinates.
!
  implicit none

  integer ( kind = 4 ) gn

  real ( kind = 8 ) gx(gn)
  real ( kind = 8 ) gy(gn)
  real ( kind = 8 ), dimension ( 25 ) :: x = (/ &
    0.13750,   0.91250,   0.71250,   0.22500,  -0.05000, &
    0.47500,   0.05000,   0.45000,   1.08750,   0.53750, &
   -0.03750,   0.18750,   0.71250,   0.85000,   0.70000, &
    0.27500,   0.45000,   0.81250,   0.45000,   1.00000, &
    0.50000,   0.18750,   0.58750,   1.05000,   0.10000 /)
  real ( kind = 8 ), dimension ( 25 ) :: y = (/ &
    0.97500,   0.98750,   0.76250,   0.83750,   0.41250, &
    0.63750,  -0.05000,   1.03750,   0.55000,   0.80000, &
    0.75000,   0.57500,   0.55000,   0.43750,   0.31250, &
    0.42500,   0.28750,   0.18750,  -0.03750,   0.26250, &
    0.46250,   0.26250,   0.12500,  -0.06125,   0.11250 /)

  gx(1:gn) = x(1:gn)
  gy(1:gn) = y(1:gn)

  return
end
subroutine g04_size ( gn )

!*****************************************************************************80
!
!! G04_SIZE returns the size for grid 4.
!
!  Discussion:
!
!    The "grid size" is simply the number of points in the grid.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) GN, the grid size.
!
  implicit none

  integer ( kind = 4 ) gn

  gn = 100

  return
end
subroutine g04_title ( gt )

!*****************************************************************************80
!
!! G04_TITLE returns the title for grid 4.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) GT, the grid title.
!
  implicit none

  character ( len = * ) gt

  gt = 'Random 100 node set'

  return
end
subroutine g04_xy ( gn, gx, gy )

!*****************************************************************************80
!
!! G04_XY returns the grid points for grid 4.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GN, the grid size.
!
!    Output, real ( kind = 8 ) GX(GN), GY(GN), the grid coordinates.
!
  implicit none

  integer ( kind = 4 ) gn

  real ( kind = 8 ) gx(gn)
  real ( kind = 8 ) gy(gn)
  real ( kind = 8 ), dimension ( 100 ) :: x = (/ &
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
  real ( kind = 8 ), dimension ( 100 ) :: y = (/ &
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

  gx(1:gn) = x(1:gn)
  gy(1:gn) = y(1:gn)

  return
end
subroutine g05_size ( gn )

!*****************************************************************************80
!
!! G05_SIZE returns the size for grid 5.
!
!  Discussion:
!
!    The "grid size" is simply the number of points in the grid.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) GN, the grid size.
!
  implicit none

  integer ( kind = 4 ) gn

  gn = 81

  return
end
subroutine g05_title ( gt )

!*****************************************************************************80
!
!! G05_TITLE returns the title for grid 5.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) GT, the grid title.
!
  implicit none

  character ( len = * ) gt

  gt = 'Gridded 81 node set'

  return
end
subroutine g05_xy ( gn, gx, gy )

!*****************************************************************************80
!
!! G05_XY returns the grid points for grid 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GN, the grid size.
!
!    Output, real ( kind = 8 ) GX(GN), GY(GN), the grid coordinates.
!
  implicit none

  integer ( kind = 4 ) gn

  real ( kind = 8 ) gx(gn)
  real ( kind = 8 ) gy(gn)
  real ( kind = 8 ), dimension ( 81 ) :: x = (/ &
    0.000,  0.000,  0.000,  0.000,  0.000, &
    0.000,  0.000,  0.000,  0.000,  0.125, &
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
  real ( kind = 8 ), dimension ( 81 ) :: y = (/ &
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

  gx(1:gn) = x(1:gn)
  gy(1:gn) = y(1:gn)

  return
end
