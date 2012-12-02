subroutine p00_f ( problem, n, x, f )

!*****************************************************************************80
!
!! P00_F evaluates the objective function for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Evelyn Beale,
!    On an Iterative Method for Finding a Local Minimum of a Function
!    of More than One Variable,
!    Technical Report 25, 
!    Statistical Techniques Research Group,
!    Princeton University, 1958.
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    John Dennis, David Gay, Phuong Vu,
!    A new nonlinear equations test problem,
!    Technical Report 83-16,
!    Mathematical Sciences Department,
!    Rice University (1983 - revised 1985).
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization 
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    Noel deVilliers, David Glasser,
!    A continuation method for nonlinear regression,
!    SIAM Journal on Numerical Analysis,
!    Volume 18, 1981, pages 1139-1154.
!
!    Chris Fraley,
!    Solution of nonlinear least-squares problems,
!    Technical Report STAN-CS-1165,
!    Computer Science Department,
!    Stanford University, 1987.
!
!    Chris Fraley,
!    Software performance on nonlinear least-squares problems,
!    Technical Report SOL 88-17,
!    Systems Optimization Laboratory,
!    Department of Operations Research,
!    Stanford University, 1988.
!
!    A Leon,
!    A Comparison of Eight Known Optimizing Procedures,
!    in Recent Advances in Optimization Techniques,
!    edited by Abraham Lavi, Thomas Vogl,
!    Wiley, 1966.
!
!    JJ McKeown,
!    Specialized versus general-purpose algorithms for functions that are sums
!    of squared terms,
!    Mathematical Programming,
!    Volume 9, 1975a, pages 57-68.
!
!    JJ McKeown,
!    On algorithms for sums of squares problems,
!    in Towards Global Optimization,
!    L. C. W. Dixon and G. Szego (eds.),
!    North-Holland, 1975, pages 229-257.
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    Algorithm 566 - Testing unconstrained optimization software,
!    ACM Transactions on Mathematical Software,
!    Volume 7, 1981, pages 17-41.
!
!    Michael Powell,
!    An Efficient Method for Finding the Minimum of a Function of
!    Several Variables Without Calculating Derivatives,
!    Computer Journal, 
!    Volume 7, Number 2, pages 155-162, 1964.
!  
!    DE Salane,
!    A continuation approach for solving large residual nonlinear least squares
!    problems,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 8, 1987, pages 655-671.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_f ( n, x, f )
  else if ( problem == 2 ) then
    call p02_f ( n, x, f )
  else if ( problem == 3 ) then
    call p03_f ( n, x, f )
  else if ( problem == 4 ) then
    call p04_f ( n, x, f )
  else if ( problem == 5 ) then
    call p05_f ( n, x, f )
  else if ( problem == 6 ) then
    call p06_f ( n, x, f )
  else if ( problem == 7 ) then
    call p07_f ( n, x, f )
  else if ( problem == 8 ) then
    call p08_f ( n, x, f )
  else if ( problem == 9 ) then
    call p09_f ( n, x, f )
  else if ( problem == 10 ) then
    call p10_f ( n, x, f )
  else if ( problem == 11 ) then
    call p11_f ( n, x, f )
  else if ( problem == 12 ) then
    call p12_f ( n, x, f )
  else if ( problem == 13 ) then
    call p13_f ( n, x, f )
  else if ( problem == 14 ) then
    call p14_f ( n, x, f )
  else if ( problem == 15 ) then
    call p15_f ( n, x, f )
  else if ( problem == 16 ) then
    call p16_f ( n, x, f )
  else if ( problem == 17 ) then
    call p17_f ( n, x, f )
  else if ( problem == 18 ) then
    call p18_f ( n, x, f )
  else if ( problem == 19 ) then
    call p19_f ( n, x, f )
  else if ( problem == 20 ) then
    call p20_f ( n, x, f )
  else if ( problem == 21 ) then
    call p21_f ( n, x, f )
  else if ( problem == 22 ) then
    call p22_f ( n, x, f )
  else if ( problem == 23 ) then
    call p23_f ( n, x, f )
  else if ( problem == 24 ) then
    call p24_f ( n, x, f )
  else if ( problem == 25 ) then
    call p25_f ( n, x, f )
  else if ( problem == 26 ) then
    call p26_f ( n, x, f )
  else if ( problem == 27 ) then
    call p27_f ( n, x, f )
  else if ( problem == 28 ) then
    call p28_f ( n, x, f )
  else if ( problem == 29 ) then
    call p29_f ( n, x, f )
  else if ( problem == 30 ) then
    call p30_f ( n, x, f )
  else if ( problem == 31 ) then
    call p31_f ( n, x, f )
  else if ( problem == 32 ) then
    call p32_f ( n, x, f )
  else if ( problem == 33 ) then
    call p33_f ( n, x, f )
  else if ( problem == 34 ) then
    call p34_f ( n, x, f )
  else if ( problem == 35 ) then
    call p35_f ( n, x, f )
  else if ( problem == 36 ) then
    call p36_f ( n, x, f )
  else if ( problem == 37 ) then
    call p37_f ( n, x, f )
  else if ( problem == 38 ) then
    call p38_f ( n, x, f )
  else if ( problem == 39 ) then
    call p39_f ( n, x, f )
  else if ( problem == 40 ) then
    call p40_f ( n, x, f )
  else if ( problem == 41 ) then
    call p41_f ( n, x, f )
  else if ( problem == 42 ) then
    call p42_f ( n, x, f )
  else if ( problem == 43 ) then
    call p43_f ( n, x, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_g ( problem, n, x, g )

!*****************************************************************************80
!
!! P00_G evaluates the gradient for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_g ( n, x, g )
  else if ( problem == 2 ) then
    call p02_g ( n, x, g )
  else if ( problem == 3 ) then
    call p03_g ( n, x, g )
  else if ( problem == 4 ) then
    call p04_g ( n, x, g )
  else if ( problem == 5 ) then
    call p05_g ( n, x, g )
  else if ( problem == 6 ) then
    call p06_g ( n, x, g )
  else if ( problem == 7 ) then
    call p07_g ( n, x, g )
  else if ( problem == 8 ) then
    call p08_g ( n, x, g )
  else if ( problem == 9 ) then
    call p09_g ( n, x, g )
  else if ( problem == 10 ) then
    call p10_g ( n, x, g )
  else if ( problem == 11 ) then
    call p11_g ( n, x, g )
  else if ( problem == 12 ) then
    call p12_g ( n, x, g )
  else if ( problem == 13 ) then
    call p13_g ( n, x, g )
  else if ( problem == 14 ) then
    call p14_g ( n, x, g )
  else if ( problem == 15 ) then
    call p15_g ( n, x, g )
  else if ( problem == 16 ) then
    call p16_g ( n, x, g )
  else if ( problem == 17 ) then
    call p17_g ( n, x, g )
  else if ( problem == 18 ) then
    call p18_g ( n, x, g )
  else if ( problem == 19 ) then
    call p19_g ( n, x, g )
  else if ( problem == 20 ) then
    call p20_g ( n, x, g )
  else if ( problem == 21 ) then
    call p21_g ( n, x, g )
  else if ( problem == 22 ) then
    call p22_g ( n, x, g )
  else if ( problem == 23 ) then
    call p23_g ( n, x, g )
  else if ( problem == 24 ) then
    call p24_g ( n, x, g )
  else if ( problem == 25 ) then
    call p25_g ( n, x, g )
  else if ( problem == 26 ) then
    call p26_g ( n, x, g )
  else if ( problem == 27 ) then
    call p27_g ( n, x, g )
  else if ( problem == 28 ) then
    call p28_g ( n, x, g )
  else if ( problem == 29 ) then
    call p29_g ( n, x, g )
  else if ( problem == 30 ) then
    call p30_g ( n, x, g )
  else if ( problem == 31 ) then
    call p31_g ( n, x, g )
  else if ( problem == 32 ) then
    call p32_g ( n, x, g )
  else if ( problem == 33 ) then
    call p33_g ( n, x, g )
  else if ( problem == 34 ) then
    call p34_g ( n, x, g )
  else if ( problem == 35 ) then
    call p35_g ( n, x, g )
  else if ( problem == 36 ) then
    call p36_g ( n, x, g )
  else if ( problem == 37 ) then
    call p37_g ( n, x, g )
  else if ( problem == 38 ) then
    call p38_g ( n, x, g )
  else if ( problem == 39 ) then
    call p39_g ( n, x, g )
  else if ( problem == 40 ) then
    call p40_g ( n, x, g )
  else if ( problem == 41 ) then
    call p41_g ( n, x, g )
  else if ( problem == 42 ) then
    call p42_g ( n, x, g )
  else if ( problem == 43 ) then
    call p43_g ( n, x, g )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_G - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_gdif ( problem, n, x, gdif )

!*****************************************************************************80
!
!! P00_GDIF approximates the gradient via finite differences.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the gradient 
!    is to be approximated.
!
!    Output, real ( kind = 8 ) GDIF(N), the approximated gradient vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dx
  real ( kind = 8 ) eps
  real ( kind = 8 ) fminus
  real ( kind = 8 ) fplus
  real ( kind = 8 ) gdif(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xi

  eps = ( epsilon ( eps ) )**0.33D+00

  do i = 1, n

    if ( 0.0D+00 <= x(i) ) then
      dx = eps * ( x(i) + 1.0D+00 )
    else
      dx = eps * ( x(i) - 1.0D+00 )
    end if

    xi = x(i)
    x(i) = xi + dx
    call p00_f ( problem, n, x, fplus )

    x(i) = xi - dx
    call p00_f ( problem, n, x, fminus )

    gdif(i) = ( fplus - fminus ) / ( 2.0D+00 * dx )

    x(i) = xi

  end do

  return
end
subroutine p00_h ( problem, n, x, h )

!*****************************************************************************80
!
!! P00_H evaluates the Hessian for any problem.
!
!  Discussion:
!
!    H(I,J) = d2 F(X) / dX(I)dX(J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_h ( n, x, h )
  else if ( problem == 2 ) then
    call p02_h ( n, x, h )
  else if ( problem == 3 ) then
    call p03_h ( n, x, h )
  else if ( problem == 4 ) then
    call p04_h ( n, x, h )
  else if ( problem == 5 ) then
    call p05_h ( n, x, h )
  else if ( problem == 6 ) then
    call p06_h ( n, x, h )
  else if ( problem == 7 ) then
    call p07_h ( n, x, h )
  else if ( problem == 8 ) then
    call p08_h ( n, x, h )
  else if ( problem == 9 ) then
    call p09_h ( n, x, h )
  else if ( problem == 10 ) then
    call p10_h ( n, x, h )
  else if ( problem == 11 ) then
    call p11_h ( n, x, h )
  else if ( problem == 12 ) then
    call p12_h ( n, x, h )
  else if ( problem == 13 ) then
    call p13_h ( n, x, h )
  else if ( problem == 14 ) then
    call p14_h ( n, x, h )
  else if ( problem == 15 ) then
    call p15_h ( n, x, h )
  else if ( problem == 16 ) then
    call p16_h ( n, x, h )
  else if ( problem == 17 ) then
    call p17_h ( n, x, h )
  else if ( problem == 18 ) then
    call p18_h ( n, x, h )
  else if ( problem == 19 ) then
    call p19_h ( n, x, h )
  else if ( problem == 20 ) then
    call p20_h ( n, x, h )
  else if ( problem == 21 ) then
    call p21_h ( n, x, h )
  else if ( problem == 22 ) then
    call p22_h ( n, x, h )
  else if ( problem == 23 ) then
    call p23_h ( n, x, h )
  else if ( problem == 24 ) then
    call p24_h ( n, x, h )
  else if ( problem == 25 ) then
    call p25_h ( n, x, h )
  else if ( problem == 26 ) then
    call p26_h ( n, x, h )
  else if ( problem == 27 ) then
    call p27_h ( n, x, h )
  else if ( problem == 28 ) then
    call p28_h ( n, x, h )
  else if ( problem == 29 ) then
    call p29_h ( n, x, h )
  else if ( problem == 30 ) then
    call p30_h ( n, x, h )
  else if ( problem == 31 ) then
    call p31_h ( n, x, h )
  else if ( problem == 32 ) then
    call p32_h ( n, x, h )
  else if ( problem == 33 ) then
    call p33_h ( n, x, h )
  else if ( problem == 34 ) then
    call p34_h ( n, x, h )
  else if ( problem == 35 ) then
    call p35_h ( n, x, h )
  else if ( problem == 36 ) then
    call p36_h ( n, x, h )
  else if ( problem == 37 ) then
    call p37_h ( n, x, h )
  else if ( problem == 38 ) then
    call p38_h ( n, x, h )
  else if ( problem == 39 ) then
    call p39_h ( n, x, h )
  else if ( problem == 40 ) then
    call p40_h ( n, x, h )
  else if ( problem == 41 ) then
    call p41_h ( n, x, h )
  else if ( problem == 42 ) then
    call p42_h ( n, x, h )
  else if ( problem == 43 ) then
    call p43_h ( n, x, h )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_H - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_hdif ( problem, n, x, hdif )

!*****************************************************************************80
!
!! P00_HDIF approximates the Hessian via finite differences.
!
!  Discussion:
!
!    The values returned by this routine will be only approximate.
!    In some cases, they will be so poor that they are useless.
!    This is particularly true for expressions in which a term like
!    X**6 is dominant.  For such terms, a small deviation in the argument
!    may hardly show up.  Using a LARGER value of EPS may sometimes help
!    in these cases.
!
!    However, one of the best applications of this routine is for
!    checking your own Hessian calculations, since as Heraclitus
!    said, you'll never get the same result twice when you differentiate
!    a complicated expression by hand.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) HDIF(N,N), the approximated N by N 
!    Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) eps
  real ( kind = 8 ) f00
  real ( kind = 8 ) fmm
  real ( kind = 8 ) fmp
  real ( kind = 8 ) fpm
  real ( kind = 8 ) fpp
  real ( kind = 8 ) hdif(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) problem
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xi
  real ( kind = 8 ) xj
!
!  Choose the stepsizes.
!
  eps = ( epsilon ( eps ) )**0.33D+00

  do i = 1, n
    s(i) = eps * max ( abs ( x(i) ), 1.0D+00 )
  end do
!
!  Calculate the diagonal elements.
!
  do i = 1, n

    xi = x(i)

    call p00_f ( problem, n, x, f00 )

    x(i) = xi + s(i)
    call p00_f ( problem, n, x, fpp )

    x(i) = xi - s(i)
    call p00_f ( problem, n, x, fmm )

    hdif(i,i) = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s(i) / s(i)

    x(i) = xi

  end do
!
!  Calculate the off-diagonal elements.
!
  do i = 1, n

    xi = x(i)

    do j = i+1, n

      xj = x(j)

      x(i) = xi + s(i)
      x(j) = xj + s(j)
      call p00_f ( problem, n, x, fpp )

      x(i) = xi + s(i)
      x(j) = xj - s(j)
      call p00_f ( problem, n, x, fpm )

      x(i) = xi - s(i)
      x(j) = xj + s(j)
      call p00_f ( problem, n, x, fmp )

      x(i) = xi - s(i)
      x(j) = xj - s(j)
      call p00_f ( problem, n, x, fmm )

      hdif(j,i) = ( ( fpp - fpm ) + ( fmm - fmp ) ) / ( 4.0D+00 * s(i) * s(j) )

      hdif(i,j) = hdif(j,i)

      x(j) = xj

    end do

    x(i) = xi

  end do

  return
end
subroutine p00_problem_num ( problem_num )

!*****************************************************************************80
!
!! P00_PROBLEM_NUM returns the number of problems available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!   Output, integer ( kind = 4 ) PROBLEM_NUM, the number of problems available.
!
  implicit none

  integer ( kind = 4 ) problem_num

  problem_num = 43

  return
end
subroutine p00_n ( problem, n )

!*****************************************************************************80
!
!! P00_N returns the number of variables for any problem.
!
!  Discussion:
!
!    Some of the problems in this set have only a single appropriate
!    size.  Others can take any value for N.  Others, alas, can take
!    SOME values of N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the number of the problem.
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_n ( n )
  else if ( problem == 2 ) then
    call p02_n ( n )
  else if ( problem == 3 ) then
    call p03_n ( n )
  else if ( problem == 4 ) then
    call p04_n ( n )
  else if ( problem == 5 ) then
    call p05_n ( n )
  else if ( problem == 6 ) then
    call p06_n ( n )
  else if ( problem == 7 ) then
    call p07_n ( n )
  else if ( problem == 8 ) then
    call p08_n ( n )
  else if ( problem == 9 ) then
    call p09_n ( n )
  else if ( problem == 10 ) then
    call p10_n ( n )
  else if ( problem == 11 ) then
    call p11_n ( n )
  else if ( problem == 12 ) then
    call p12_n ( n )
  else if ( problem == 13 ) then
    call p13_n ( n )
  else if ( problem == 14 ) then
    call p14_n ( n )
  else if ( problem == 15 ) then
    call p15_n ( n )
  else if ( problem == 16 ) then
    call p16_n ( n )
  else if ( problem == 17 ) then
    call p17_n ( n )
  else if ( problem == 18 ) then
    call p18_n ( n )
  else if ( problem == 19 ) then
    call p19_n ( n )
  else if ( problem == 20 ) then
    call p20_n ( n )
  else if ( problem == 21 ) then
    call p21_n ( n )
  else if ( problem == 22 ) then
    call p22_n ( n )
  else if ( problem == 23 ) then
    call p23_n ( n )
  else if ( problem == 24 ) then
    call p24_n ( n )
  else if ( problem == 25 ) then
    call p25_n ( n )
  else if ( problem == 26 ) then
    call p26_n ( n )
  else if ( problem == 27 ) then
    call p27_n ( n )
  else if ( problem == 28 ) then
    call p28_n ( n )
  else if ( problem == 29 ) then
    call p29_n ( n )
  else if ( problem == 30 ) then
    call p30_n ( n )
  else if ( problem == 31 ) then
    call p31_n ( n )
  else if ( problem == 32 ) then
    call p32_n ( n )
  else if ( problem == 33 ) then
    call p33_n ( n )
  else if ( problem == 34 ) then
    call p34_n ( n )
  else if ( problem == 35 ) then
    call p35_n ( n )
  else if ( problem == 36 ) then
    call p36_n ( n )
  else if ( problem == 37 ) then
    call p37_n ( n )
  else if ( problem == 38 ) then
    call p38_n ( n )
  else if ( problem == 39 ) then
    call p39_n ( n )
  else if ( problem == 40 ) then
    call p40_n ( n )
  else if ( problem == 41 ) then
    call p41_n ( n )
  else if ( problem == 42 ) then
    call p42_n ( n )
  else if ( problem == 43 ) then
    call p43_n ( n )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_N - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_sol ( problem, n, know, x )

!*****************************************************************************80
!
!! P00_SOL returns the solution for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_sol ( n, know, x )
  else if ( problem == 2 ) then
    call p02_sol ( n, know, x )
  else if ( problem == 3 ) then
    call p03_sol ( n, know, x )
  else if ( problem == 4 ) then
    call p04_sol ( n, know, x )
  else if ( problem == 5 ) then
    call p05_sol ( n, know, x )
  else if ( problem == 6 ) then
    call p06_sol ( n, know, x )
  else if ( problem == 7 ) then
    call p07_sol ( n, know, x )
  else if ( problem == 8 ) then
    call p08_sol ( n, know, x )
  else if ( problem == 9 ) then
    call p09_sol ( n, know, x )
  else if ( problem == 10 ) then
    call p10_sol ( n, know, x )
  else if ( problem == 11 ) then
    call p11_sol ( n, know, x )
  else if ( problem == 12 ) then
    call p12_sol ( n, know, x )
  else if ( problem == 13 ) then
    call p13_sol ( n, know, x )
  else if ( problem == 14 ) then
    call p14_sol ( n, know, x )
  else if ( problem == 15 ) then
    call p15_sol ( n, know, x )
  else if ( problem == 16 ) then
    call p16_sol ( n, know, x )
  else if ( problem == 17 ) then
    call p17_sol ( n, know, x )
  else if ( problem == 18 ) then
    call p18_sol ( n, know, x )
  else if ( problem == 19 ) then
    call p19_sol ( n, know, x )
  else if ( problem == 20 ) then
    call p20_sol ( n, know, x )
  else if ( problem == 21 ) then
    call p21_sol ( n, know, x )
  else if ( problem == 22 ) then
    call p22_sol ( n, know, x )
  else if ( problem == 23 ) then
    call p23_sol ( n, know, x )
  else if ( problem == 24 ) then
    call p24_sol ( n, know, x )
  else if ( problem == 25 ) then
    call p25_sol ( n, know, x )
  else if ( problem == 26 ) then
    call p26_sol ( n, know, x )
  else if ( problem == 27 ) then
    call p27_sol ( n, know, x )
  else if ( problem == 28 ) then
    call p28_sol ( n, know, x )
  else if ( problem == 29 ) then
    call p29_sol ( n, know, x )
  else if ( problem == 30 ) then
    call p30_sol ( n, know, x )
  else if ( problem == 31 ) then
    call p31_sol ( n, know, x )
  else if ( problem == 32 ) then
    call p32_sol ( n, know, x )
  else if ( problem == 33 ) then
    call p33_sol ( n, know, x )
  else if ( problem == 34 ) then
    call p34_sol ( n, know, x )
  else if ( problem == 35 ) then
    call p35_sol ( n, know, x )
  else if ( problem == 36 ) then
    call p36_sol ( n, know, x )
  else if ( problem == 37 ) then
    call p37_sol ( n, know, x )
  else if ( problem == 38 ) then
    call p38_sol ( n, know, x )
  else if ( problem == 39 ) then
    call p39_sol ( n, know, x )
  else if ( problem == 40 ) then
    call p40_sol ( n, know, x )
  else if ( problem == 41 ) then
    call p41_sol ( n, know, x )
  else if ( problem == 42 ) then
    call p42_sol ( n, know, x )
  else if ( problem == 43 ) then
    call p43_sol ( n, know, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SOL - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_start ( problem, n, x )

!*****************************************************************************80
!
!! P00_START returns a starting point for optimization for any problem.
!
!  Discussion:
!
!    The point returned by this routine does not produce an optimal
!    value of the objective function.  Instead, it is "reasonably far"
!    from an optimizing point, so that an optimization program has
!    a proper workout.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the number of the problem.
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_start ( n, x )
  else if ( problem == 2 ) then
    call p02_start ( n, x )
  else if ( problem == 3 ) then
    call p03_start ( n, x )
  else if ( problem == 4 ) then
    call p04_start ( n, x )
  else if ( problem == 5 ) then
    call p05_start ( n, x )
  else if ( problem == 6 ) then
    call p06_start ( n, x )
  else if ( problem == 7 ) then
    call p07_start ( n, x )
  else if ( problem == 8 ) then
    call p08_start ( n, x )
  else if ( problem == 9 ) then
    call p09_start ( n, x )
  else if ( problem == 10 ) then
    call p10_start ( n, x )
  else if ( problem == 11 ) then
    call p11_start ( n, x )
  else if ( problem == 12 ) then
    call p12_start ( n, x )
  else if ( problem == 13 ) then
    call p13_start ( n, x )
  else if ( problem == 14 ) then
    call p14_start ( n, x )
  else if ( problem == 15 ) then
    call p15_start ( n, x )
  else if ( problem == 16 ) then
    call p16_start ( n, x )
  else if ( problem == 17 ) then
    call p17_start ( n, x )
  else if ( problem == 18 ) then
    call p18_start ( n, x )
  else if ( problem == 19 ) then
    call p19_start ( n, x )
  else if ( problem == 20 ) then
    call p20_start ( n, x )
  else if ( problem == 21 ) then
    call p21_start ( n, x )
  else if ( problem == 22 ) then
    call p22_start ( n, x )
  else if ( problem == 23 ) then
    call p23_start ( n, x )
  else if ( problem == 24 ) then
    call p24_start ( n, x )
  else if ( problem == 25 ) then
    call p25_start ( n, x )
  else if ( problem == 26 ) then
    call p26_start ( n, x )
  else if ( problem == 27 ) then
    call p27_start ( n, x )
  else if ( problem == 28 ) then
    call p28_start ( n, x )
  else if ( problem == 29 ) then
    call p29_start ( n, x )
  else if ( problem == 30 ) then
    call p30_start ( n, x )
  else if ( problem == 31 ) then
    call p31_start ( n, x )
  else if ( problem == 32 ) then
    call p32_start ( n, x )
  else if ( problem == 33 ) then
    call p33_start ( n, x )
  else if ( problem == 34 ) then
    call p34_start ( n, x )
  else if ( problem == 35 ) then
    call p35_start ( n, x )
  else if ( problem == 36 ) then
    call p36_start ( n, x )
  else if ( problem == 37 ) then
    call p37_start ( n, x )
  else if ( problem == 38 ) then
    call p38_start ( n, x )
  else if ( problem == 39 ) then
    call p39_start ( n, x )
  else if ( problem == 40 ) then
    call p40_start ( n, x )
  else if ( problem == 41 ) then
    call p41_start ( n, x )
  else if ( problem == 42 ) then
    call p42_start ( n, x )
  else if ( problem == 43 ) then
    call p43_start ( n, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_START - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of nPROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_title ( problem, title )

!*****************************************************************************80
!
!! P00_TITLE returns a title for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the number of the problem.
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer ( kind = 4 ) problem
  character ( len = * ) title

  if ( problem == 1 ) then
    call p01_title ( title )
  else if ( problem == 2 ) then
    call p02_title ( title )
  else if ( problem == 3 ) then
    call p03_title ( title )
  else if ( problem == 4 ) then
    call p04_title ( title )
  else if ( problem == 5 ) then
    call p05_title ( title )
  else if ( problem == 6 ) then
    call p06_title ( title )
  else if ( problem == 7 ) then
    call p07_title ( title )
  else if ( problem == 8 ) then
    call p08_title ( title )
  else if ( problem == 9 ) then
    call p09_title ( title )
  else if ( problem == 10 ) then
    call p10_title ( title )
  else if ( problem == 11 ) then
    call p11_title ( title )
  else if ( problem == 12 ) then
    call p12_title ( title )
  else if ( problem == 13 ) then
    call p13_title ( title )
  else if ( problem == 14 ) then
    call p14_title ( title )
  else if ( problem == 15 ) then
    call p15_title ( title )
  else if ( problem == 16 ) then
    call p16_title ( title )
  else if ( problem == 17 ) then
    call p17_title ( title )
  else if ( problem == 18 ) then
    call p18_title ( title )
  else if ( problem == 19 ) then
    call p19_title ( title )
  else if ( problem == 20 ) then
    call p20_title ( title )
  else if ( problem == 21 ) then
    call p21_title ( title )
  else if ( problem == 22 ) then
    call p22_title ( title )
  else if ( problem == 23 ) then
    call p23_title ( title )
  else if ( problem == 24 ) then
    call p24_title ( title )
  else if ( problem == 25 ) then
    call p25_title ( title )
  else if ( problem == 26 ) then
    call p26_title ( title )
  else if ( problem == 27 ) then
    call p27_title ( title )
  else if ( problem == 28 ) then
    call p28_title ( title )
  else if ( problem == 29 ) then
    call p29_title ( title )
  else if ( problem == 30 ) then
    call p30_title ( title )
  else if ( problem == 31 ) then
    call p31_title ( title )
  else if ( problem == 32 ) then
    call p32_title ( title )
  else if ( problem == 33 ) then
    call p33_title ( title )
  else if ( problem == 34 ) then
    call p34_title ( title )
  else if ( problem == 35 ) then
    call p35_title ( title )
  else if ( problem == 36 ) then
    call p36_title ( title )
  else if ( problem == 37 ) then
    call p37_title ( title )
  else if ( problem == 38 ) then
    call p38_title ( title )
  else if ( problem == 39 ) then
    call p39_title ( title )
  else if ( problem == 40 ) then
    call p40_title ( title )
  else if ( problem == 41 ) then
    call p41_title ( title )
  else if ( problem == 42 ) then
    call p42_title ( title )
  else if ( problem == 43 ) then
    call p43_title ( title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p01_f ( n, x, f )

!*****************************************************************************80
!
!! P01_F evaluates the objective function for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)

  call p01_th ( x, th )

  f = 100.0D+00 * ( x(3) - 10.0D+00 * th )**2 &
    + 100.0D+00 * ( sqrt ( x(1) * x(1) + x(2) * x(2) ) - 1.0D+00 )**2 &
    + x(3) * x(3)

  return
end
subroutine p01_g ( n, x, g )

!*****************************************************************************80
!
!! P01_G evaluates the gradient for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) s1
  real ( kind = 8 ) t
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)

  call p01_th ( x, th )

  r = sqrt ( x(1) * x(1) + x(2) * x(2) )
  t = x(3) - 10.0D+00 * th
  s1 = 5.0D+00 * t / ( pi * r * r )

  g(1) = 200.0D+00 * ( x(1) - x(1) / r + x(2) * s1 )
  g(2) = 200.0D+00 * ( x(2) - x(2) / r - x(1) * s1 )
  g(3) = 2.0D+00 * ( 100.0D+00 * t + x(3) )

  return
end
subroutine p01_h ( n, x, h )

!*****************************************************************************80
!
!! P01_H evaluates the Hessian for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)

  call p01_th ( x, th )

  h(1,1) = 200.0D+00 - 200.0D+00 * x(2)**2 &
         * ( 1.0D+00 / sqrt ( x(1)**2 + x(2)**2 )**3 &
         - 25.0D+00 / ( pi * ( x(1)**2 + x(2)**2 ) )**2 ) &
         - 2000.0D+00 * x(1) * x(2) * ( x(3) - 10.0D+00 * th ) &
         / ( pi * ( x(1)**2 + x(2)**2 )**2 )

  h(1,2) = 200.0D+00 * x(1) * x(2) / sqrt ( x(1)**2 + x(2)**2 )**3 &
         + 1000.0D+00 / ( pi * ( x(1)**2 + x(2)**2 )**2 ) &
         * ( ( x(3) - 10.0D+00 * th ) * ( x(1)**2 - x(2)**2 ) &
         - 5.0D+00 * x(1) * x(2) / pi )

  h(1,3) = 1000.0D+00 * x(2) / ( pi * ( x(1)**2 + x(2)**2 ) )

  h(2,1) = 200.0D+00 * x(1) * x(2) / sqrt ( x(1)**2 + x(2)**2 )**3 &
         + 1000.0D+00 / ( pi * ( x(1)**2 + x(2)**2 )**2 ) &
         * ( ( x(3) - 10.0D+00 * th ) * ( x(1)**2 - x(2)**2 ) &
         - 5.0D+00 * x(1) * x(2) / pi )

  h(2,2) = 200.0D+00 - 200.0D+00 * x(1)**2 &
         * ( 1.0D+00 / sqrt ( x(1)**2 + x(2)**2 )**3 &
         - 25.0D+00 / ( pi * ( x(1)**2 + x(2)**2 ) )**2 ) &
         + 2000.0D+00 * x(1) * x(2) * ( x(3) - 10.0D+00 * th ) &
         / ( pi * ( x(1)**2 + x(2)**2 )**2 )

  h(2,3) = - 1000.0D+00 * x(1) / ( pi * ( x(1)**2 + x(2)**2 ) )

  h(3,1) = 1000.0D+00 * x(2) / ( pi * ( x(1)**2 + x(2)**2 ) )
  h(3,2) = - 1000.0D+00 * x(1) / ( pi * ( x(1)**2 + x(2)**2 ) )
  h(3,3) = 202.0D+00

  return
end
subroutine p01_n ( n )

!*****************************************************************************80
!
!! P01_N returns the number of variables for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p01_sol ( n, know, x )

!*****************************************************************************80
!
!! P01_SOL returns the solution for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    x(1:3) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
    know = know + 1
  else
    know = 0
  end if

  return
end
subroutine p01_start ( n, x )

!*****************************************************************************80
!
!! P01_START returns a starting point for optimization for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:3) = (/ -1.0D+00, 0.0D+00, 0.0D+00 /)

  return
end
subroutine p01_th ( x, th )

!*****************************************************************************80
!
!! P01_TH evaluates a term used by problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(2), the values of the variables.
!
!    Output, real ( kind = 8 ) TH, the value of the term.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) th
  real ( kind = 8 ) x(2)

  if ( 0.0D+00 < x(1) ) then
    th = 0.5D+00 * atan ( x(2) / x(1) ) / pi
  else if ( x(1) < 0.0D+00 ) then
    th = 0.5D+00 * atan ( x(2) / x(1) ) / pi + 0.5D+00
  else if ( 0.0D+00 < x(2) ) then
    th = 0.25D+00
  else if ( x(2) < 0.0D+00 ) then
    th = - 0.25D+00
  else
    th = 0.0D+00
  end if

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns a title for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Fletcher-Powell helical valley function.'

  return
end
subroutine p02_f ( n, x, f )

!*****************************************************************************80
!
!! P02_F evaluates the objective function for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) f
  real ( kind = 8 ) fi
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  f = 0.0D+00

  do i = 1, 13

    c = - real ( i, kind = 8 ) / 10.0D+00

    fi = x(3)     * exp ( c * x(1) )         - x(4) * exp ( c * x(2) ) &
       + x(6)     * exp ( c * x(5) )         -        exp ( c ) &
       + 5.0D+00  * exp ( 10.0D+00 * c ) - 3.0D+00  * exp ( 4.0D+00 * c )

    f = f + fi * fi

  end do

  return
end
subroutine p02_g ( n, x, g )

!*****************************************************************************80
!
!! P02_G evaluates the gradient for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) dfdx1
  real ( kind = 8 ) dfdx2
  real ( kind = 8 ) dfdx3
  real ( kind = 8 ) dfdx4
  real ( kind = 8 ) dfdx5
  real ( kind = 8 ) dfdx6
  real ( kind = 8 ) fi
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  do i = 1, 13

    c = - real ( i, kind = 8 ) / 10.0D+00

    fi =        x(3) * exp ( c * x(1) )         - x(4) * exp ( c * x(2) ) &
          +     x(6) * exp ( c * x(5) )         -        exp ( c ) &
          +  5.0D+00 * exp ( 10.0D+00 * c ) -  3.0D+00 * exp ( 4.0D+00 * c )

    dfdx1 =     c    * x(3) * exp ( c * x(1) ) 
    dfdx2 =   - c    * x(4) * exp ( c * x(2) )
    dfdx3 =                   exp ( c * x(1) )
    dfdx4 =   -               exp ( c * x(2) )
    dfdx5 =     c    * x(6) * exp ( c * x(5) )
    dfdx6 =                   exp ( c * x(5) )

    g(1) = g(1) + 2.0D+00 * fi * dfdx1
    g(2) = g(2) + 2.0D+00 * fi * dfdx2
    g(3) = g(3) + 2.0D+00 * fi * dfdx3
    g(4) = g(4) + 2.0D+00 * fi * dfdx4
    g(5) = g(5) + 2.0D+00 * fi * dfdx5
    g(6) = g(6) + 2.0D+00 * fi * dfdx6

  end do

  return
end
subroutine p02_h ( n, x, h )

!*****************************************************************************80
!
!! P02_H evaluates the Hessian for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) dfdx1
  real ( kind = 8 ) dfdx2
  real ( kind = 8 ) dfdx3
  real ( kind = 8 ) dfdx4
  real ( kind = 8 ) dfdx5
  real ( kind = 8 ) dfdx6
  real ( kind = 8 ) d2fdx11
  real ( kind = 8 ) d2fdx13
  real ( kind = 8 ) d2fdx22
  real ( kind = 8 ) d2fdx24
  real ( kind = 8 ) d2fdx31
  real ( kind = 8 ) d2fdx42
  real ( kind = 8 ) d2fdx55
  real ( kind = 8 ) d2fdx56
  real ( kind = 8 ) d2fdx65
  real ( kind = 8 ) fi
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do i = 1, 13

    c = - real ( i, kind = 8 ) / 10.0D+00

    fi =   x(3) * exp ( c * x(1) )  - x(4) * exp ( c * x(2) ) &
          + x(6) * exp ( c * x(5) ) -        exp ( c ) &
          + 5.0D+00 * exp ( 10.0D+00 * c ) -  3.0D+00 * exp ( 4.0D+00 * c )

    dfdx1 =     c     * x(3) * exp ( c * x(1) ) 
    d2fdx11 =   c * c * x(3) * exp ( c * x(1) ) 
    d2fdx13 =   c     *        exp ( c * x(1) ) 
    dfdx2 =   - c     * x(4) * exp ( c * x(2) )
    d2fdx22 = - c * c * x(4) * exp ( c * x(2) )
    d2fdx24 = - c     *        exp ( c * x(2) )
    dfdx3 =                    exp ( c * x(1) )
    d2fdx31 =   c     *        exp ( c * x(1) ) 
    dfdx4 =   -                exp ( c * x(2) )
    d2fdx42 = - c     *        exp ( c * x(2) )
    dfdx5 =     c     * x(6) * exp ( c * x(5) )
    d2fdx55 =   c * c * x(6) * exp ( c * x(5) )
    d2fdx56 =   c     *        exp ( c * x(5) )
    dfdx6 =                    exp ( c * x(5) )
    d2fdx65 =   c     *        exp ( c * x(5) )

    h(1,1) = h(1,1) + 2.0D+00 * dfdx1 * dfdx1 + 2.0D+00 * fi * d2fdx11
    h(1,2) = h(1,2) + 2.0D+00 * dfdx2 * dfdx1 
    h(1,3) = h(1,3) + 2.0D+00 * dfdx3 * dfdx1 + 2.0D+00 * fi * d2fdx13
    h(1,4) = h(1,4) + 2.0D+00 * dfdx4 * dfdx1 
    h(1,5) = h(1,5) + 2.0D+00 * dfdx5 * dfdx1 
    h(1,6) = h(1,6) + 2.0D+00 * dfdx6 * dfdx1 

    h(2,1) = h(2,1) + 2.0D+00 * dfdx1 * dfdx2 
    h(2,2) = h(2,2) + 2.0D+00 * dfdx2 * dfdx2 + 2.0D+00 * fi * d2fdx22
    h(2,3) = h(2,3) + 2.0D+00 * dfdx3 * dfdx2 
    h(2,4) = h(2,4) + 2.0D+00 * dfdx4 * dfdx2 + 2.0D+00 * fi * d2fdx24
    h(2,5) = h(2,5) + 2.0D+00 * dfdx5 * dfdx2 
    h(2,6) = h(2,6) + 2.0D+00 * dfdx6 * dfdx2 

    h(3,1) = h(3,1) + 2.0D+00 * dfdx1 * dfdx3 + 2.0D+00 * fi * d2fdx31
    h(3,2) = h(3,2) + 2.0D+00 * dfdx2 * dfdx3 
    h(3,3) = h(3,3) + 2.0D+00 * dfdx3 * dfdx3 
    h(3,4) = h(3,4) + 2.0D+00 * dfdx4 * dfdx3
    h(3,5) = h(3,5) + 2.0D+00 * dfdx5 * dfdx3 
    h(3,6) = h(3,6) + 2.0D+00 * dfdx6 * dfdx3 

    h(4,1) = h(4,1) + 2.0D+00 * dfdx1 * dfdx4 
    h(4,2) = h(4,2) + 2.0D+00 * dfdx2 * dfdx4 + 2.0D+00 * fi * d2fdx42
    h(4,3) = h(4,3) + 2.0D+00 * dfdx3 * dfdx4 
    h(4,4) = h(4,4) + 2.0D+00 * dfdx4 * dfdx4 
    h(4,5) = h(4,5) + 2.0D+00 * dfdx5 * dfdx4 
    h(4,6) = h(4,6) + 2.0D+00 * dfdx6 * dfdx4 

    h(5,1) = h(5,1) + 2.0D+00 * dfdx1 * dfdx5 
    h(5,2) = h(5,2) + 2.0D+00 * dfdx2 * dfdx5 
    h(5,3) = h(5,3) + 2.0D+00 * dfdx3 * dfdx5 
    h(5,4) = h(5,4) + 2.0D+00 * dfdx4 * dfdx5 
    h(5,5) = h(5,5) + 2.0D+00 * dfdx5 * dfdx5 + 2.0D+00 * fi * d2fdx55
    h(5,6) = h(5,6) + 2.0D+00 * dfdx6 * dfdx5 + 2.0D+00 * fi * d2fdx56

    h(6,1) = h(6,1) + 2.0D+00 * dfdx1 * dfdx6 
    h(6,2) = h(6,2) + 2.0D+00 * dfdx2 * dfdx6 
    h(6,3) = h(6,3) + 2.0D+00 * dfdx3 * dfdx6 
    h(6,4) = h(6,4) + 2.0D+00 * dfdx4 * dfdx6
    h(6,5) = h(6,5) + 2.0D+00 * dfdx5 * dfdx6 + 2.0D+00 * fi * d2fdx65
    h(6,6) = h(6,6) + 2.0D+00 * dfdx6 * dfdx6

  end do

  return
end
subroutine p02_n ( n )

!*****************************************************************************80
!
!! P02_N returns the number of variables for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 6

  return
end
subroutine p02_sol ( n, know, x )

!*****************************************************************************80
!
!! P02_SOL returns the solution for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:6) = (/ 1.0D+00, 10.0D+00, 1.0D+00, 5.0D+00, 4.0D+00, 3.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p02_start ( n, x )

!*****************************************************************************80
!
!! P02_START returns a starting point for optimization for problem 2.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:6) = (/ 1.0D+00, 2.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns a title for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Biggs EXP6 function.'

  return
end
subroutine p03_f ( n, x, f )

!*****************************************************************************80
!
!! P03_F evaluates the objective function for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(15)

  call p03_yvec ( y ) 

  f = 0.0D+00

  do i = 1, 15

    t = x(1) * exp ( - 0.5D+00 * x(2) * &
      ( 3.5D+00 - 0.5D+00 * real ( i - 1, kind = 8 ) - x(3) )**2 ) - y(i)

    f = f + t * t

  end do

  return
end
subroutine p03_g ( n, x, g )

!*****************************************************************************80
!
!! P03_G evaluates the gradient for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(15)

  call p03_yvec ( y )

  g(1:n) = 0.0D+00

  do i = 1, 15

    d1 = 0.5D+00 * real ( i - 1, kind = 8 )
    d2 = 3.5D+00 - d1 - x(3)
    arg = - 0.5D+00 * x(2) * d2 * d2
    t = x(1) * exp ( arg ) - y(i)

    g(1) = g(1) + 2.0D+00 * exp ( arg ) * t
    g(2) = g(2) - x(1) * exp ( arg ) * t * d2 * d2
    g(3) = g(3) + 2.0D+00 * x(1) * x(2) * exp ( arg ) * t * d2

  end do

  return
end
subroutine p03_yvec ( y )

!*****************************************************************************80
!
!! P03_YVEC is an auxilliary routine for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) Y(15), data values needed for the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) y(15)

  y(1:15) = (/ 0.0009D+00, 0.0044D+00, 0.0175D+00, 0.0540D+00, 0.1295D+00, &
               0.2420D+00, 0.3521D+00, 0.3989D+00, 0.3521D+00, 0.2420D+00, &
               0.1295D+00, 0.0540D+00, 0.0175D+00, 0.0044D+00, 0.0009D+00 /)

  return
end
subroutine p03_h ( n, x, h )

!*****************************************************************************80
!
!! P03_H evaluates the Hessian for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(15)

  call p03_yvec ( y )

  h(1:n,1:n) = 0.0D+00

  do i = 1, 15

    d1 = 0.5D+00 * real ( i - 1, kind = 8 )
    d2 = 3.5D+00 - d1 - x(3)
    arg = 0.5D+00 * x(2) * d2 * d2
    r = exp ( - arg )
    t = x(1) * r - y(i)
    t1 = 2.0D+00 * x(1) * r - y(i)

    h(1,1) = h(1,1) + r * r
    h(2,2) = h(2,2) + r * t1 * d2**4
    h(3,3) = h(3,3) + r * ( x(2) * t1 * d2 * d2 - t )
    h(2,1) = h(2,1) - r * t1 * d2 * d2
    h(3,1) = h(3,1) + d2 * r * t1
    h(3,2) = h(3,2) + d2 * r * ( t - arg * t1 )

  end do

  h(1,1) = 2.0D+00 * h(1,1)
  h(2,2) = 0.5D+00 * x(1) * h(2,2)
  h(3,3) = 2.0D+00 * x(1) * x(2) * h(3,3)
  h(3,1) = 2.0D+00 * x(2) * h(3,1)
  h(3,2) = 2.0D+00 * x(1) * h(3,2)

  h(1,2) = h(2,1)
  h(1,3) = h(3,1)
  h(2,3) = h(3,2)

  return
end
subroutine p03_n ( n )

!*****************************************************************************80
!
!! P03_N returns the number of variables for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p03_sol ( n, know, x )

!*****************************************************************************80
!
!! P03_SOL returns the solution for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  know = 0
  x(1:n) = 0.0D+00

  return
end
subroutine p03_start ( n, x )

!*****************************************************************************80
!
!! P03_START returns a starting point for optimization for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:3) = (/ 0.4D+00, 1.0D+00, 0.0D+00 /)

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns a title for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Gaussian function.'

  return
end
subroutine p04_f ( n, x, f )

!*****************************************************************************80
!
!! P04_F evaluates the objective function for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) x(n)

  f1 = 10000.0D+00 * x(1) * x(2) - 1.0D+00
  f2 = exp ( - x(1) ) + exp ( - x(2) ) - 1.0001D+00

  f = f1 * f1 + f2 * f2

  return
end
subroutine p04_g ( n, x, g )

!*****************************************************************************80
!
!! P04_G evaluates the gradient for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df1dx2
  real ( kind = 8 ) df2dx1
  real ( kind = 8 ) df2dx2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  f1 = 10000.0D+00 * x(1) * x(2) - 1.0D+00
  df1dx1 = 10000.0D+00 * x(2)
  df1dx2 = 10000.0D+00 * x(1)

  f2 = exp ( - x(1) ) + exp ( - x(2) ) - 1.0001D+00
  df2dx1 = - exp ( - x(1) )
  df2dx2 = - exp ( - x(2) )

  g(1) = 2.0D+00 * f1 * df1dx1 + 2.0D+00 * f2 * df2dx1
  g(2) = 2.0D+00 * f1 * df1dx2 + 2.0D+00 * f2 * df2dx2

  return
end
subroutine p04_h ( n, x, h )

!*****************************************************************************80
!
!! P04_H evaluates the Hessian for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d2f1dx12
  real ( kind = 8 ) d2f1dx21
  real ( kind = 8 ) d2f2dx11
  real ( kind = 8 ) d2f2dx22
  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df1dx2
  real ( kind = 8 ) df2dx1
  real ( kind = 8 ) df2dx2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  f1 = 10000.0D+00 * x(1) * x(2) - 1.0D+00
  df1dx1 = 10000.0D+00 * x(2)
  df1dx2 = 10000.0D+00 * x(1)
  d2f1dx12 = 10000.0D+00
  d2f1dx21 = 10000.0D+00

  f2 = exp ( - x(1) ) + exp ( - x(2) ) - 1.0001D+00
  df2dx1 = - exp ( - x(1) )
  df2dx2 = - exp ( - x(2) )
  d2f2dx11 = exp ( - x(1) )
  d2f2dx22 = exp ( - x(2) )

  h(1,1) = 2.0D+00 * df1dx1 * df1dx1 &
    + 2.0D+00 * f2 * d2f2dx11 + 2.0D+00 * df2dx1 * df2dx1

  h(1,2) = 2.0D+00 * f1 * d2f1dx12 + 2.0D+00 * df1dx2 * df1dx1 &
    + 2.0D+00 * df2dx2 * df2dx1

  h(2,1) = 2.0D+00 * f1 * d2f1dx21 + 2.0D+00 * df1dx2 * df1dx1 &
    + 2.0D+00 * df2dx2 * df2dx1

  h(2,2) = 2.0D+00 * df1dx2 * df1dx2 &
    + 2.0D+00 * f2 * d2f2dx22 + 2.0D+00 * df2dx2 * df2dx2

  return
end
subroutine p04_n ( n )

!*****************************************************************************80
!
!! P04_N returns the number of variables for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p04_sol ( n, know, x )

!*****************************************************************************80
!
!! P04_SOL returns the solution for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:2) = (/ 1.098159D-05, 9.106146D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p04_start ( n, x )

!*****************************************************************************80
!
!! P04_START returns a starting point for optimization for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ 0.0D+00, 1.0D+00 /)

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns a title for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Powell badly scaled function.'

  return
end
subroutine p05_f ( n, x, f )

!*****************************************************************************80
!
!! P05_F evaluates the objective function for problem 5.
!
!  Discussion:
!
!    The function is formed by the sum of squares of 10 separate terms.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) f
  real ( kind = 8 ) fi
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  f = 0.0D+00

  do i = 1, 10

    c = - real ( i, kind = 8 ) / 10.0D+00

    fi = exp ( c * x(1) ) - exp ( c * x(2) ) - x(3) * &
      ( exp ( c ) - exp ( 10.0D+00 * c ) )
   
    f = f + fi * fi

  end do

  return
end
subroutine p05_g ( n, x, g )

!*****************************************************************************80
!
!! P05_G evaluates the gradient for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) dfidx1
  real ( kind = 8 ) dfidx2
  real ( kind = 8 ) dfidx3
  real ( kind = 8 ) fi
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  do i = 1, 10

    c = - real ( i, kind = 8 ) / 10.0D+00

    fi = exp ( c * x(1) ) - exp ( c * x(2) ) &
      - x(3) * ( exp ( c ) - exp ( 10.0D+00 * c ) )

    dfidx1 =   c * exp ( c * x(1) )
    dfidx2 = - c * exp ( c * x(2) )
    dfidx3 = - ( exp ( c ) - exp ( 10.0D+00 * c ) )

    g(1) = g(1) + 2.0D+00 * fi * dfidx1
    g(2) = g(2) + 2.0D+00 * fi * dfidx2
    g(3) = g(3) + 2.0D+00 * fi * dfidx3

  end do

  return
end
subroutine p05_h ( n, x, h )

!*****************************************************************************80
!
!! P05_H evaluates the Hessian for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) d2fidx11
  real ( kind = 8 ) d2fidx22
  real ( kind = 8 ) dfidx1
  real ( kind = 8 ) dfidx2
  real ( kind = 8 ) dfidx3
  real ( kind = 8 ) fi
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do i = 1, 10

    c = - real ( i, kind = 8 ) / 10.0D+00

    fi = exp ( c * x(1) ) - exp ( c * x(2) ) &
      - x(3) * ( exp ( c ) - exp ( 10.0D+00 * c ) )

    dfidx1 =     c     * exp ( c * x(1) )
    d2fidx11 =   c * c * exp ( c * x(1) )
    dfidx2 =   - c     * exp ( c * x(2) )
    d2fidx22 = - c * c * exp ( c * x(2) )
    dfidx3 =  - ( exp ( c ) - exp ( 10.0D+00 * c ) )

    h(1,1) = h(1,1) + 2.0D+00 * dfidx1 * dfidx1 + 2.0D+00 * fi * d2fidx11 
    h(1,2) = h(1,2) + 2.0D+00 * dfidx1 * dfidx2
    h(1,3) = h(1,3) + 2.0D+00 * dfidx1 * dfidx3

    h(2,1) = h(2,1) + 2.0D+00 * dfidx2 * dfidx1
    h(2,2) = h(2,2) + 2.0D+00 * dfidx2 * dfidx2 + 2.0D+00 * fi * d2fidx22 
    h(2,3) = h(2,3) + 2.0D+00 * dfidx2 * dfidx3

    h(3,1) = h(3,1) + 2.0D+00 * dfidx3 * dfidx1
    h(3,2) = h(3,2) + 2.0D+00 * dfidx3 * dfidx2
    h(3,3) = h(3,3) + 2.0D+00 * dfidx3 * dfidx3

  end do

  return
end
subroutine p05_n ( n )

!*****************************************************************************80
!
!! P05_N returns the number of variables for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p05_sol ( n, know, x )

!*****************************************************************************80
!
!! P05_SOL returns the solution for problem 5.
!
!  Discussion:
!
!    The function has a minimum of 0 at (1,10,1) and also for
!    any point of the form (x,x,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:3) = (/ 1.0D+00, 10.0D+00, 1.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p05_start ( n, x )

!*****************************************************************************80
!
!! P05_START returns a starting point for optimization for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:3) = (/ 0.0D+00, 10.0D+00, 5.0D+00 /)

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns a title for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Box 3-dimensional function.'

  return
end
subroutine p06_f ( n, x, f )

!*****************************************************************************80
!
!! P06_F evaluates the objective function for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  f1 = 0.0D+00
  do i = 1, n
    f1 = f1 + real ( i, kind = 8 ) * ( x(i) - 1.0D+00 )
  end do

  f2 = 0.0D+00
  do i = 1, n
    f2 = f2 + ( x(i) - 1.0D+00 )**2
  end do

  f = f1 * f1 * ( 1.0D+00 + f1 * f1 ) + f2

  return
end
subroutine p06_g ( n, x, g )

!*****************************************************************************80
!
!! P06_G evaluates the gradient for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) df1dxi
  real ( kind = 8 ) df2dxi
  real ( kind = 8 ) f1
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  f1 = 0.0D+00
  do i = 1, n
    f1 = f1 + real ( i, kind = 8 ) * ( x(i) - 1.0D+00 )
  end do

  do i = 1, n
    df1dxi = real ( i, kind = 8 )
    df2dxi = 2.0D+00 * ( x(i) - 1.0D+00 )
    g(i) = ( 2.0D+00 * f1 + 4.0D+00 * f1**3 ) * df1dxi + df2dxi
  end do

  return
end
subroutine p06_h ( n, x, h )

!*****************************************************************************80
!
!! P06_H evaluates the Hessian for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d2f2dxii
  real ( kind = 8 ) df1dxi
  real ( kind = 8 ) df1dxj
  real ( kind = 8 ) f1
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  f1 = 0.0D+00
  do i = 1, n
    f1 = f1 + real ( i, kind = 8 ) * ( x(i) - 1.0D+00 )
  end do

  do i = 1, n
    df1dxi = real ( i, kind = 8 )
    d2f2dxii = 2.0D+00
    do j = 1, n
      df1dxj = real ( j, kind = 8 )
      h(i,j) = ( 2.0D+00 + 12.0D+00 * f1 * f1 ) * df1dxi * df1dxj
    end do
    h(i,i) = h(i,i) + d2f2dxii
  end do

  return
end
subroutine p06_n ( n )

!*****************************************************************************80
!
!! P06_N returns the number of variables for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p06_sol ( n, know, x )

!*****************************************************************************80
!
!! P06_SOL returns the solution for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = 1.0D+00
  else
    know = 0
  end if

  return
end
subroutine p06_start ( n, x )

!*****************************************************************************80
!
!! P06_START returns a starting point for optimization for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = real ( n - i, kind = 8 ) / real ( n, kind = 8 )
  end do

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns a title for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The variably dimensioned function.'

  return
end
subroutine p07_f ( n, x, f )

!*****************************************************************************80
!
!! P07_F evaluates the objective function for problem 7.
!
!  Discussion:
!
!    For N = 9, the problem of minimizing the Watson function is
!    very ill conditioned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) x(n)

  f = 0.0D+00
  do i = 1, 29

    s1 = 0.0D+00
    d = 1.0D+00
    do j = 2, n
      s1 = s1 + real ( j - 1, kind = 8 ) * d * x(j)
      d = d * real ( i, kind = 8 ) / 29.0D+00
    end do

    s2 = 0.0D+00
    d = 1.0D+00
    do j = 1, n
      s2 = s2 + d * x(j)
      d = d * real ( i, kind = 8 ) / 29.0D+00
    end do

    f = f + ( s1 - s2 * s2 - 1.0D+00 )**2

  end do

  f = f + x(1) * x(1) + ( x(2) - x(1) * x(1) - 1.0D+00 )**2

  return
end
subroutine p07_g ( n, x, g )

!*****************************************************************************80
!
!! P07_G evaluates the gradient for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) s3
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  do i = 1, 29

    d1 = real ( i, kind = 8 ) / 29.0D+00
    s1 = 0.0D+00
    d2 = 1.0D+00
    do j = 2, n
      s1 = s1 + real ( j - 1, kind = 8 ) * d2 * x(j)
      d2 = d2 * real ( i, kind = 8 ) / 29.0D+00
    end do

    s2 = 0.0D+00
    d2 = 1.0D+00
    do j = 1, n
      s2 = s2 + d2 * x(j)
      d2 = d2 * real ( i, kind = 8 ) / 29.0D+00
    end do

    t = s1 - s2 * s2 - 1.0D+00
    s3 = 2.0D+00 * s2 * real ( i, kind = 8 ) / 29.0D+00
    d2 = 2.0D+00 / d1

    do j = 1, n
      g(j) = g(j) + d2 * ( real ( j - 1, kind = 8 ) - s3 ) * t
      d2 = d2 * real ( i, kind = 8 ) / 29.0D+00
    end do

  end do

  t1 = x(2) - x(1) * x(1) - 1.0D+00

  g(1) = g(1) + 2.0D+00 * x(1) - 4.0D+00 * x(1) * t1
  g(2) = g(2) + 2.0D+00 * t1

  return
end
subroutine p07_h ( n, x, h )

!*****************************************************************************80
!
!! P07_H evaluates the Hessian for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) s3
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t3
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do i = 1, 29

    d1 = real ( i, kind = 8 ) / 29.0D+00
    d2 = 1.0D+00
    s1 = 0.0D+00
    s2 = x(1)
    do j = 2, n
      s1 = s1 + real ( j - 1, kind = 8 ) * d2 * x(j)
      d2 = d1 * d2
      s2 = s2 + d2 * x(j)
    end do
    t = 2.0D+00 * ( s1 - s2 * s2 - 1.0D+00 ) * d1 * d1
    s3 = 2.0D+00 * d1 * s2
    d2 = 1.0D+00 / d1
    do j = 1, n
      t1 = real ( j - 1, kind = 8 ) - s3
      h(j,j) = h(j,j) + 2.0D+00 * ( t1 * t1 - t ) * d2 * d2
      d3 = 1.0D+00 / d1
      do k = 1, j-1
        h(j,k) = h(j,k) &
          + 2.0D+00 * ( t1 * ( real ( k - 1, kind = 8 ) - s3 ) - t ) * d2 * d3
        d3 = d1 * d3
      end do
      d2 = d1 * d2
    end do

  end do

  t3 = x(2) - x(1) * x(1) - 1.0D+00
  h(1,1) = h(1,1) + 2.0D+00 - 4.0D+00 * ( t3 - 2.0D+00 * x(1) * x(1) )
  h(2,2) = h(2,2) + 2.0D+00
  h(2,1) = h(2,1) - 4.0D+00 * x(1)

  do i = 1, n
    h(i,i+1:n) = h(i+1:n,i)
  end do

  return
end
subroutine p07_n ( n )

!*****************************************************************************80
!
!! P07_N returns the number of variables for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 2

  return
end
subroutine p07_sol ( n, know, x )

!*****************************************************************************80
!
!! P07_SOL returns the solution for problem 7.
!
!  Discussion:
!
!    The values of the approximate solutions are taken from Brent.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 .and. n == 6 ) then
    know = 1
    x(1:n) = (/ -0.015725D+00, 1.012435D+00, -0.232992D+00, 1.260430D+00, &
      -1.513729D+00, 0.992996D+00 /)
  else if ( know == 0 .and. n == 9 ) then
    know = 1
    x(1:n) = (/ -0.000015D+00, 0.999790D+00, 0.014764D+00, 0.146342D+00, &
      1.000821D+00, -2.617731D+00, 4.104403D+00, -3.143612D+00, 1.052627D+00 /)
  else
    know = 0
    x(1:n) = 0.0D+00
  end if

  return
end
subroutine p07_start ( n, x )

!*****************************************************************************80
!
!! P07_START returns a starting point for optimization for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = 0.0D+00

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns a title for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Watson function.'

  return
end
subroutine p08_f ( n, x, f )

!*****************************************************************************80
!
!! P08_F evaluates the objective function for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: ap = 0.00001D+00
  real ( kind = 8 ) f
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) x(n)

  t1 = - 0.25D+00 + sum ( x(1:n)**2 )

  t2 = sum ( ( x(1:n) - 1.0D+00 )**2 )

  f = ap * t2 + t1 * t1

  return
end
subroutine p08_g ( n, x, g )

!*****************************************************************************80
!
!! P08_G evaluates the gradient for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: ap = 0.00001D+00
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) t1
  real ( kind = 8 ) x(n)

  t1 = - 0.25D+00 + sum ( x(1:n)**2 )

  g(1:n) = 2.0D+00 * ap * ( x(1:n) - 1.0D+00 ) + 4.0D+00 * x(1:n) * t1

  return
end
subroutine p08_h ( n, x, h )

!*****************************************************************************80
!
!! P08_H evaluates the Hessian for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: ap = 0.00001D+00
  real ( kind = 8 ) d1
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t1
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)

  t1 = - 0.25D+00 + sum ( x(1:n)**2 )

  d1 = 2.0D+00 * ap
  th = 4.0D+00 * t1

  do i = 1, n
    h(i,i) = d1 + th + 8.0D+00 * x(i)**2
    do j = 1, i-1
      h(i,j) = 8.0D+00 * x(i) * x(j)
    end do
  end do

  do i = 1, n
    do j = i+1, n
      h(i,j) = h(j,i)
    end do
  end do

  return
end
subroutine p08_n ( n )

!*****************************************************************************80
!
!! P08_N returns the number of variables for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p08_sol ( n, know, x )

!*****************************************************************************80
!
!! P08_SOL returns the solution for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  know = 0
  x(1:n) = 0.0D+00

  return
end
subroutine p08_start ( n, x )

!*****************************************************************************80
!
!! P08_START returns a starting point for optimization for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns a title for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Penalty Function #1.'

  return
end
subroutine p09_f ( n, x, f )

!*****************************************************************************80
!
!! P09_F evaluates the objective function for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: ap = 0.00001D+00
  real ( kind = 8 ) d2
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) s3
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) x(n)

  t1 = -1.0D+00
  t2 = 0.0D+00
  t3 = 0.0D+00
  d2 = 1.0D+00
  s2 = 0.0D+00
  do j = 1, n
    t1 = t1 + real ( n - j + 1, kind = 8 ) * x(j)**2
    s1 = exp ( x(j) / 10.0D+00 )
    if ( 1 < j ) then
      s3 = s1 + s2 - d2 * ( exp ( 0.1D+00 ) + 1.0D+00 )
      t2 = t2 + s3 * s3
      t3 = t3 + ( s1 - 1.0D+00 / exp ( 0.1D+00 ) )**2
    end if
    s2 = s1
    d2 = d2 * exp ( 0.1D+00 )
  end do

  f = ap * ( t2 + t3 ) + t1 * t1 + ( x(1) - 0.2D+00 )**2

  return
end
subroutine p09_g ( n, x, g )

!*****************************************************************************80
!
!! P09_G evaluates the gradient for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: ap = 0.00001D+00
  real ( kind = 8 ) d2
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) s3
  real ( kind = 8 ) t1
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)

  t1 = -1.0D+00
  do j = 1, n
    t1 = t1 + real ( n - j + 1, kind = 8 ) * x(j)**2
  end do

  d2 = 1.0D+00
  th = 4.0D+00 * t1
  s2 = 0.0D+00
  do j = 1, n
    g(j) = real ( n - j + 1, kind = 8 ) * x(j) * th
    s1 = exp ( x(j) / 10.0D+00 )
    if ( 1 < j ) then
      s3 = s1 + s2 - d2 * ( exp ( 0.1D+00 ) + 1.0D+00 )
      g(j) = g(j) + ap * s1 * ( s3 + s1 - 1.0D+00 / exp ( 0.1D+00 ) ) / 5.0D+00
      g(j-1) = g(j-1) + ap * s2 * s3 / 5.0D+00
    end if
    s2 = s1
    d2 = d2 * exp ( 0.1D+00 )
  end do

  g(1) = g(1) + 2.0D+00 * ( x(1) - 0.2D+00 )

  return
end
subroutine p09_h ( n, x, h )

!*****************************************************************************80
!
!! P09_H evaluates the Hessian for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: ap = 0.00001D+00
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) s3
  real ( kind = 8 ) t1
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)

  t1 = - 1.0D+00
  do j = 1, n
    t1 = t1 + real ( n - j + 1, kind = 8 ) * x(j) * x(j)
  end do

  d1 = exp ( 0.1D+00 )
  d2 = 1.0D+00
  s2 = 0.0D+00
  th = 4.0D+00 * t1

  do j = 1, n

    h(j,j) = 8.0D+00 * ( real ( n - j + 1, kind = 8 ) * x(j) )**2 &
      + real ( n - j + 1, kind = 8 ) * th

    s1 = exp ( x(j) / 10.0D+00 )

    if ( 1 < j ) then

      s3 = s1 + s2 - d2 * ( d1 + 1.0D+00 )
      h(j,j) = h(j,j) + ap * s1 * ( s3 + s1 - 1.0D+00 &
        / d1 + 2.0D+00 * s1 ) / 50.0D+00
      h(j-1,j-1) = h(j-1,j-1) + ap * s2 * ( s2 + s3 ) / 50.0D+00
      do k = 1, j-1
        t1 = exp ( real ( k, kind = 8 ) / 10.0D+00 )
        h(j,k) = 8.0D+00 * real ( n - j + 1, kind = 8 ) &
          * real ( n - k + 1, kind = 8 ) * x(j) * x(k)
      end do

      h(j,j-1) = h(j,j-1) + ap * s1 * s2 / 50.0D+00

    end if

    s2 = s1
    d2 = d1 * d2

  end do

  h(1,1) = h(1,1) + 2.0D+00

  do i = 1, n
    h(i,i+1:n) = h(i+1:n,i)
  end do

  return
end
subroutine p09_n ( n )

!*****************************************************************************80
!
!! P09_N returns the number of variables for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p09_sol ( n, know, x )

!*****************************************************************************80
!
!! P09_SOL returns the solution for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  know = 0
  x(1:n) = 0.0D+00

  return
end
subroutine p09_start ( n, x )

!*****************************************************************************80
!
!! P09_START returns a starting point for optimization for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = 0.5D+00

  return
end
subroutine p09_title ( title )

!*****************************************************************************80
!
!! P09_TITLE returns a title for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Penalty Function #2.'

  return
end
subroutine p10_f ( n, x, f )

!*****************************************************************************80
!
!! P10_F evaluates the objective function for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) x(n)

  f = ( x(1) - 1000000.0D+00 )**2 &
    + ( x(2) - 0.000002D+00 )**2 &
    + ( x(1) * x(2) - 2.0D+00 )**2

  return
end
subroutine p10_g ( n, x, g )

!*****************************************************************************80
!
!! P10_G evaluates the gradient for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  g(1) = 2.0D+00 * x(1) - 2000000.0D+00 + 2.0D+00 * x(1) * x(2) * x(2) &
    - 4.0D+00 * x(2)

  g(2) = 2.0D+00 * x(2) - 0.000004 + 2.0D+00 * x(1)**2 * x(2) - 4.0D+00 * x(1)

  return
end
subroutine p10_h ( n, x, h )

!*****************************************************************************80
!
!! P10_H evaluates the Hessian for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1,1) = 2.0D+00 + 2.0D+00 * x(2) * x(2)
  h(1,2) = 4.0D+00 * x(1) * x(2) - 4.0D+00

  h(2,1) = 4.0D+00 * x(1) * x(2) - 4.0D+00
  h(2,2) = 2.0D+00 + 2.0D+00 * x(1) * x(1)

  return
end
subroutine p10_n ( n )

!*****************************************************************************80
!
!! P10_N returns the number of variables for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p10_sol ( n, know, x )

!*****************************************************************************80
!
!! P10_SOL returns the solution for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:2) = (/ 1.0D+06, 2.0D-06 /)
  else
    know = 0
  end if

  return
end
subroutine p10_start ( n, x )

!*****************************************************************************80
!
!! P10_START returns a starting point for optimization for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ 1.0D+00, 1.0D+00 /)

  return
end
subroutine p10_title ( title )

!*****************************************************************************80
!
!! P10_TITLE returns a title for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Brown Badly Scaled Function.'

  return
end
subroutine p11_f ( n, x, f )

!*****************************************************************************80
!
!! P11_F evaluates the objective function for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  f = 0.0D+00

  do i = 1, 20

    c = real ( i, kind = 8 ) / 5.0D+00
    f1 = x(1) + c * x(2) - exp ( c )
    f2 = x(3) + sin ( c ) * x(4) - cos ( c )

    f = f + f1**4 + 2.0D+00 * f1 * f1 * f2 * f2 + f2**4

  end do

  return
end
subroutine p11_g ( n, x, g )

!*****************************************************************************80
!
!! P11_G evaluates the gradient for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df1dx2
  real ( kind = 8 ) df2dx3
  real ( kind = 8 ) df2dx4
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  do i = 1, 20

    c = real ( i, kind = 8 ) / 5.0D+00

    f1 = x(1) + c * x(2) - exp ( c )
    f2 = x(3) + sin ( c ) * x(4) - cos ( c )

    df1dx1 = 1.0D+00
    df1dx2 = c
    df2dx3 = 1.0D+00
    df2dx4 = sin ( c )

    g(1) = g(1) + 4.0D+00 * ( f1**3 * df1dx1 + f1 * f2 * f2 * df1dx1 )
    g(2) = g(2) + 4.0D+00 * ( f1**3 * df1dx2 + f1 * f2 * f2 * df1dx2 )
    g(3) = g(3) + 4.0D+00 * ( f1 * f1 * f2 * df2dx3 + f2**3 * df2dx3 )
    g(4) = g(4) + 4.0D+00 * ( f1 * f1 * f2 * df2dx4 + f2**3 * df2dx4 )

  end do

  return
end
subroutine p11_h ( n, x, h )

!*****************************************************************************80
!
!! P11_H evaluates the Hessian for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df1dx2
  real ( kind = 8 ) df2dx3
  real ( kind = 8 ) df2dx4
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do i = 1, 20

    c = real ( i, kind = 8 ) / 5.0D+00

    f1 = x(1) + c * x(2) - exp ( c )
    f2 = x(3) + sin ( c ) * x(4) - cos ( c )

    df1dx1 = 1.0D+00
    df1dx2 = c
    df2dx3 = 1.0D+00
    df2dx4 = sin ( c )

    h(1,1) = h(1,1) + 12.0D+00 * f1 * f1 * df1dx1 * df1dx1 &
                    +  4.0D+00 * f2 * f2 * df1dx1 * df1dx1
    h(1,2) = h(1,2) + 12.0D+00 * f1 * f1 * df1dx1 * df1dx2 &
                    +  4.0D+00 * f2 * f2 * df1dx1 * df1dx2
    h(1,3) = h(1,3) +  8.0D+00 * f1 * f2 * df1dx1 * df2dx3
    h(1,4) = h(1,4) +  8.0D+00 * f1 * f2 * df1dx1 * df2dx4

    h(2,1) = h(2,1) + 12.0D+00 * f1 * f1 * df1dx2 * df1dx1 &
                    +  4.0D+00 * f2 * f2 * df1dx2 * df1dx1
    h(2,2) = h(2,2) + 12.0D+00 * f1 * f1 * df1dx2 * df1dx2 &
                    +  4.0D+00 * f2 * f2 * df1dx2 * df1dx1
    h(2,3) = h(2,3) +  8.0D+00 * f1 * f2 * df1dx2 * df2dx3
    h(2,4) = h(2,4) +  8.0D+00 * f1 * f2 * df1dx2 * df2dx4

    h(3,1) = h(3,1) +  8.0D+00 * f1 * f2 * df2dx3 * df1dx1
    h(3,2) = h(3,2) +  8.0D+00 * f1 * f2 * df2dx3 * df1dx2
    h(3,3) = h(3,3) +  4.0D+00 * f1 * f1 * df2dx3 * df2dx3 &
                    + 12.0D+00 * f2 * f2 * df2dx3 * df2dx3
    h(3,4) = h(3,4) +  4.0D+00 * f1 * f1 * df2dx4 * df2dx3 &
                    + 12.0D+00 * f2 * f2 * df2dx3 * df2dx4

    h(4,1) = h(4,1) +  8.0D+00 * f1 * f2 * df2dx4 * df1dx1
    h(4,2) = h(4,2) +  8.0D+00 * f1 * f2 * df2dx4 * df1dx2
    h(4,3) = h(4,3) +  4.0D+00 * f1 * f1 * df2dx3 * df2dx4 &
                    + 12.0D+00 * f2 * f2 * df2dx4 * df2dx3
    h(4,4) = h(4,4) +  4.0D+00 * f1 * f1 * df2dx4 * df2dx4 &
                    + 12.0D+00 * f2 * f2 * df2dx4 * df2dx4

  end do

  return
end
subroutine p11_n ( n )

!*****************************************************************************80
!
!! P11_N returns the number of variables for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 4

  return
end
subroutine p11_sol ( n, know, x )

!*****************************************************************************80
!
!! P11_SOL returns the solution for problem 11.
!
!  Discussion:
!
!    A local minimizer is approximately known.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ -11.5844D+00, 13.1999D+00, -0.406200D+00, 0.240998D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p11_start ( n, x )

!*****************************************************************************80
!
!! P11_START returns a starting point for optimization for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:4) =  (/ 25.0D+00, 5.0D+00, -5.0D+00, -1.0D+00 /)

  return
end
subroutine p11_title ( title )

!*****************************************************************************80
!
!! P11_TITLE returns a title for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Brown and Dennis Function.'

  return
end
subroutine p12_f ( n, x, f )

!*****************************************************************************80
!
!! P12_F evaluates the objective function for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  f = 0.0D+00
  do i = 1, 99
    arg = real ( i, kind = 8 ) / 100.0D+00
    r = abs ( ( - 50.0D+00 * log ( arg ) )**( 2.0D+00 / 3.0D+00 ) &
      + 25.0D+00 - x(2) )

    t = exp ( - r**x(3) / x(1) ) - arg

    f = f + t * t

  end do

  return
end
subroutine p12_g ( n, x, g )

!*****************************************************************************80
!
!! P12_G evaluates the gradient for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  do i = 1, 99

    arg = real ( i, kind = 8 ) / 100.0D+00
    r = abs ( ( - 50.0D+00 * log ( arg ) )**( 2.0D+00 / 3.0D+00 ) &
      + 25.0D+00 - x(2) )
    t1 = r**x(3) / x(1)
    t2 = exp ( - t1 )

    t = exp ( - r**x(3) / x(1) ) - arg

    g(1) = g(1) + 2.0D+00 * t * t1 * t2 / x(1)
    g(2) = g(2) + 2.0D+00 * t * t1 * t2 * x(3) / r
    g(3) = g(3) - 2.0D+00 * t * t1 * t2 * log ( r )

  end do

  return
end
subroutine p12_h ( n, x, h )

!*****************************************************************************80
!
!! P12_H evaluates the Hessian for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) d1
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) logr
  real ( kind = 8 ) r
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  d1 = 2.0D+00 / 3.0D+00

  do i = 1, 99

    arg = real ( i, kind = 8 ) / 100.0D+00
    r = ( - 50.0D+00 * log ( arg ) )**d1 + 25.0D+00 - x(2)
    t1 = abs ( r )**x(3) / x(1)
    t2 = exp ( - t1 )
    t3 = t1 * t2 * ( t1 * t2 + ( t1 - 1.0D+00 ) * ( t2 - arg ) )
    t = t1 * t2 * ( t2 - arg )
    logr = log ( abs ( r ) )

    h(1,1) = h(1,1) + 2.0D+00 * t3 - 2.0D+00 * t
    h(1,2) = h(1,2) + 2.0D+00 * t3 / r
    h(1,3) = h(1,3) - 2.0D+00 * t3 * logr

    h(2,1) = h(2,1) + 2.0D+00 * t3 / r
    h(2,2) = h(2,2) + 2.0D+00 * ( t + x(3) * t3 ) / r / r
    h(2,3) = h(2,3) + 2.0D+00 * ( t - x(3) * t3 * logr ) / r

    h(3,1) = h(3,1) - 2.0D+00 * t3 * logr
    h(3,2) = h(3,2) + 2.0D+00 * ( t - x(3) * t3 * logr ) / r
    h(3,3) = h(3,3) + 2.0D+00 * t3 * logr * logr

  end do

  h(1,1) = ( h(1,1) / x(1) ) / x(1)
  h(1,2) = h(1,2) * x(3) / x(1)
  h(1,3) = h(1,3) / x(1)

  h(2,1) = h(2,1) * x(3) / x(1)
  h(2,2) = h(2,2) * x(3)

  h(3,1) = h(3,1) / x(1)

  return
end
subroutine p12_n ( n )

!*****************************************************************************80
!
!! P12_N returns the number of variables for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p12_sol ( n, know, x )

!*****************************************************************************80
!
!! P12_SOL returns the solution for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:3) = (/ 50.0D+00, 25.0D+00, 1.5D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p12_start ( n, x )

!*****************************************************************************80
!
!! P12_START returns a starting point for optimization for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:3) = (/ 40.0D+00, 20.0D+00, 1.20D+00 /)

  return
end
subroutine p12_title ( title )

!*****************************************************************************80
!
!! P12_TITLE returns a title for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Gulf R&D Function.'

  return
end
subroutine p13_f ( n, x, f )

!*****************************************************************************80
!
!! P13_F evaluates the objective function for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) s1
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  s1 = sum ( cos ( x(1:n) ) )

  f = 0.0D+00
  do j = 1, n
    t = real ( n + j, kind = 8 ) - sin ( x(j) ) &
      - s1 - real ( j, kind = 8 ) * cos ( x(j) )
    f = f + t * t
  end do

  return
end
subroutine p13_g ( n, x, g )

!*****************************************************************************80
!
!! P13_G evaluates the gradient for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  s1 = sum ( cos ( x(1:n) ) )

  s2 = 0.0D+00
  do j = 1, n
    t = real ( n + j, kind = 8 ) - sin ( x(j) ) &
      - s1 - real ( j, kind = 8 ) * cos ( x(j) )
    s2 = s2 + t
    g(j) = ( real ( j, kind = 8 ) * sin ( x(j) ) - cos ( x(j) ) ) * t
  end do

  do j = 1, n
    g(j) = 2.0D+00 * ( g(j) + sin ( x(j) ) * s2 )
  end do

  return
end
subroutine p13_h ( n, x, h )

!*****************************************************************************80
!
!! P13_H evaluates the Hessian for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) t
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  s1 = sum ( cos ( x(1:n) ) )

  do j = 1, n
    h(j,j) = sin ( x(j) )
  end do

  s2 = 0.0D+00

  do j = 1, n

    th = cos ( x(j) )
    t = real ( n + j, kind = 8 ) - h(j,j) - s1 - real ( j, kind = 8 ) * th
    s2 = s2 + t
    do k = 1, j-1
      h(j,k) = 2.0D+00 * ( sin ( x(k) ) &
        * ( real ( n + j + k, kind = 8 ) * h(j,j) - th ) - &
        h(j,j) * cos ( x(k) ) )
    end do

    h(j,j) = real ( j * ( j + 2 ) + n, kind = 8 ) * h(j,j) * h(j,j) + th * &
      ( th - real ( 2 * j + 2, kind = 8 ) * h(j,j) ) &
      + t * ( real ( j, kind = 8 ) * th + h(j,j) )

  end do

  do j = 1, n
    h(j,j) = 2.0D+00 * ( h(j,j) + cos ( x(j) ) * s2 )
  end do

  do i = 1, n
    h(i,i+1:n) = h(i+1:n,i)
  end do

  return
end
subroutine p13_n ( n )

!*****************************************************************************80
!
!! P13_N returns the number of variables for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p13_sol ( n, know, x )

!*****************************************************************************80
!
!! P13_SOL returns the solution for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  know = 0
  x(1:n) = 0.0D+00

  return
end
subroutine p13_start ( n, x )

!*****************************************************************************80
!
!! P13_START returns a starting point for optimization for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = 1.0D+00 / real ( n, kind = 8 )

  return
end
subroutine p13_title ( title )

!*****************************************************************************80
!
!! P13_TITLE returns a title for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Trigonometric Function.'

  return
end
subroutine p14_f ( n, x, f )

!*****************************************************************************80
!
!! P14_F evaluates the objective function for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  f = 0.0D+00
  do j = 1, n
    if ( mod ( j, 2 ) == 1 ) then
      f = f + ( 1.0D+00 - x(j) )**2
    else
      f = f + 100.0D+00 * ( x(j) - x(j-1) * x(j-1) )**2
    end if
  end do

  return
end
subroutine p14_g ( n, x, g )

!*****************************************************************************80
!
!! P14_G evaluates the gradient for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  do j = 1, n

    if ( mod ( j, 2 ) == 1 ) then
      g(j) = - 2.0D+00 * ( 1.0D+00 - x(j) )
    else
      g(j) = 200.0D+00 * ( x(j) - x(j-1) * x(j-1) )
      g(j-1) = g(j-1) - 400.0D+00 * x(j-1) * ( x(j) - x(j-1) * x(j-1) )
    end if

  end do

  return
end
subroutine p14_h ( n, x, h )

!*****************************************************************************80
!
!! P14_H evaluates the Hessian for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do j = 1, n

    if ( mod ( j, 2 ) == 1 ) then
      h(j,j) = 2.0D+00
    else
      h(j,j) = 200.0D+00
      h(j,j-1) = - 400.0D+00 * x(j-1)
      h(j-1,j) = h(j-1,j) - 400.0D+00 * x(j-1) 
      h(j-1,j-1) = h(j-1,j-1) - 400.0D+00 * ( x(j) - 3.0D+00 * x(j-1) * x(j-1) )
    end if

  end do

  return
end
subroutine p14_n ( n )

!*****************************************************************************80
!
!! P14_N returns the number of variables for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p14_sol ( n, know, x )

!*****************************************************************************80
!
!! P14_SOL returns the solution for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = 1.0D+00
  else
    know = 0
  end if

  return
end
subroutine p14_start ( n, x )

!*****************************************************************************80
!
!! P14_START returns a starting point for optimization for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( mod ( i, 2 ) == 1 ) then
      x(i) = - 1.2D+00
    else
      x(i) = 1.0D+00
    end if
  end do

  return
end
subroutine p14_title ( title )

!*****************************************************************************80
!
!! P14_TITLE returns a title for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Extended Rosenbrock parabolic valley Function.'

  return
end
subroutine p15_f ( n, x, f )

!*****************************************************************************80
!
!! P15_F evaluates the objective function for problem 15.
!
!  Discussion:
!
!    The Hessian matrix is doubly singular at the minimizer,
!    suggesting that most optimization routines will experience
!    a severe slowdown in convergence.
!
!    The problem is usually only defined for N being a multiple of 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xjp1
  real ( kind = 8 ) xjp2
  real ( kind = 8 ) xjp3

  f = 0.0D+00

  do j = 1, n, 4

    if ( j + 1 <= n ) then
      xjp1 = x(j+1)
    else
      xjp1 = 0.0D+00
    end if

    if ( j + 2 <= n ) then
      xjp2 = x(j+2)
    else
      xjp2 = 0.0D+00
    end if

    if ( j + 3 <= n ) then
      xjp3 = x(j+3)
    else
      xjp3 = 0.0D+00
    end if
 
    f1 = x(j) + 10.0D+00 * xjp1

    if ( j + 1 <= n ) then
      f2 = xjp2 - xjp3
    else
      f2 = 0.0D+00
    end if

    if ( j + 2 <= n ) then
      f3 = xjp1 - 2.0D+00 * xjp2
    else
      f3 = 0.0D+00
    end if

    if ( j + 3 <= n ) then
      f4 = x(j) - xjp3
    else
      f4 = 0.0D+00
    end if

    f = f +            f1 * f1 &
          +  5.0D+00 * f2 * f2 &
          +            f3 * f3 * f3 * f3 &
          + 10.0D+00 * f4 * f4 * f4 * f4

  end do

  return
end
subroutine p15_g ( n, x, g )

!*****************************************************************************80
!
!! P15_G evaluates the gradient for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) df1dxj
  real ( kind = 8 ) df1dxjp1
  real ( kind = 8 ) df2dxjp2
  real ( kind = 8 ) df2dxjp3
  real ( kind = 8 ) df3dxjp1
  real ( kind = 8 ) df3dxjp2
  real ( kind = 8 ) df4dxj
  real ( kind = 8 ) df4dxjp3
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xjp1
  real ( kind = 8 ) xjp2
  real ( kind = 8 ) xjp3

  do j = 1, n, 4

    if ( j + 1 <= n ) then
      xjp1 = x(j+1)
    else
      xjp1 = 0.0D+00
    end if

    if ( j + 2 <= n ) then
      xjp2 = x(j+2)
    else
      xjp2 = 0.0D+00
    end if

    if ( j + 3 <= n ) then
      xjp3 = x(j+3)
    else
      xjp3 = 0.0D+00
    end if

    f1 = x(j) + 10.0D+00 * xjp1
    df1dxj = 1.0D+00
    df1dxjp1 = 10.0D+00

    if ( j + 1 <= n ) then
      f2 = xjp2 - xjp3
      df2dxjp2 = 1.0D+00
      df2dxjp3 = -1.0D+00
    else
      f2 = 0.0D+00
      df2dxjp2 = 0.0D+00
      df2dxjp3 = 0.0D+00
    end if

    if ( j + 2 <= n ) then
      f3 = xjp1 - 2.0D+00 * xjp2
      df3dxjp1 = 1.0D+00
      df3dxjp2 = -2.0D+00
    else
      f3 = 0.0D+00
      df3dxjp1 = 0.0D+00
      df3dxjp2 = 0.0D+00
    end if

    if ( j + 3 <= n ) then
      f4 = x(j) - xjp3
      df4dxj = 1.0D+00
      df4dxjp3 = -1.0D+00
    else
      f4 = 0.0D+00
      df4dxj = 0.0D+00
      df4dxjp3 = 0.0D+00
    end if

    g(j) = 2.0D+00 * f1 * df1dxj + 10.0D+00 * 4.0D+00 * f4**3 * df4dxj

    if ( j+1 <= n ) then
      g(j+1) = 2.0D+00 * f1 * df1dxjp1 + 4.0D+00 * f3**3 * df3dxjp1
    end if

    if ( j+2 <= n ) then
      g(j+2) = 2.0D+00 * 5.0D+00 * f2 * df2dxjp2 + 4.0D+00 * f3**3 * df3dxjp2
    end if

    if ( j+3 <= n ) then
      g(j+3) = 2.0D+00 * 5.0D+00 * f2 * df2dxjp3 &
        + 10.0D+00 * 4.0D+00 * f4**3 * df4dxjp3
    end if

  end do

  return
end
subroutine p15_h ( n, x, h )

!*****************************************************************************80
!
!! P15_H evaluates the Hessian for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) df1dxj
  real ( kind = 8 ) df1dxjp1
  real ( kind = 8 ) df2dxjp2
  real ( kind = 8 ) df2dxjp3
  real ( kind = 8 ) df3dxjp1
  real ( kind = 8 ) df3dxjp2
  real ( kind = 8 ) df4dxj
  real ( kind = 8 ) df4dxjp3
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xjp1
  real ( kind = 8 ) xjp2
  real ( kind = 8 ) xjp3

  h(1:n,1:n) = 0.0D+00

  do j = 1, n, 4

    if ( j + 1 <= n ) then
      xjp1 = x(j+1)
    else
      xjp1 = 0.0D+00
    end if

    if ( j + 2 <= n ) then
      xjp2 = x(j+2)
    else
      xjp2 = 0.0D+00
    end if

    if ( j + 3 <= n ) then
      xjp3 = x(j+3)
    else
      xjp3 = 0.0D+00
    end if

    f1 = x(j) + 10.0D+00 * xjp1
    df1dxj = 1.0D+00
    df1dxjp1 = 10.0D+00

    if ( j + 1 <= n ) then
      f2 = xjp2 - xjp3
      df2dxjp2 = 1.0D+00
      df2dxjp3 = -1.0D+00
    else
      f2 = 0.0D+00
      df2dxjp2 = 0.0D+00
      df2dxjp3 = 0.0D+00
    end if

    if ( j + 2 <= n ) then
      f3 = xjp1 - 2.0D+00 * xjp2
      df3dxjp1 = 1.0D+00
      df3dxjp2 = -2.0D+00
    else
      f3 = 0.0D+00
      df3dxjp1 = 0.0D+00
      df3dxjp2 = 0.0D+00
    end if

    if ( j + 3 <= n ) then
      f4 = x(j) - xjp3
      df4dxj = 1.0D+00
      df4dxjp3 = -1.0D+00
    else
      f4 = 0.0D+00
      df4dxj = 0.0D+00
      df4dxjp3 = 0.0D+00
    end if

    h(j  ,j  ) = 2.0D+00 * df1dxj * df1dxj &
               + 120.0D+00 * f4 * f4 * df4dxj * df4dxj

    if ( j + 1 <= n ) then
      h(j  ,j+1) = 2.0D+00 * df1dxj * df1dxjp1
      h(j+1,j  ) = 2.0D+00 * df1dxj * df1dxjp1
      h(j+1,j+1) = 2.0D+00 * df1dxjp1 * df1dxjp1 &
                 + 12.0D+00 * f3 * f3 * df3dxjp1 * df3dxjp1
    end if

    if ( j + 2 <= n ) then
      h(j  ,j+2) = 0.0D+00
      h(j+2,j  ) = 0.0D+00
      h(j+1,j+2) = 12.0D+00 * f3 * f3 * df3dxjp2 * df3dxjp1
      h(j+2,j+1) = 12.0D+00 * f3 * f3 * df3dxjp1 * df3dxjp2
      h(j+2,j+2) = 10.0D+00 * df2dxjp2 * df2dxjp2 &
                 + 12.0D+00 * f3 * f3 * df3dxjp2 * df3dxjp2
    end if

    if ( j + 3 <= n ) then
      h(j  ,j+3) = 120.0D+00 * f4 * f4 * df4dxj * df4dxjp3
      h(j+3,j  ) = 120.0D+00 * f4 * f4 * df4dxj * df4dxjp3
      h(j+1,j+3) = 0.0D+00
      h(j+3,j+1) = 0.0D+00
      h(j+2,j+3) = 10.0D+00 * df2dxjp3 * df2dxjp2
      h(j+3,j+2) = 10.0D+00 * df2dxjp2 * df2dxjp3
      h(j+3,j+3) = 10.0D+00 * df2dxjp3 * df2dxjp3 &
                 + 120.0D+00 * f4 * f4 * df4dxjp3 * df4dxjp3
    end if

  end do

  return
end
subroutine p15_n ( n )

!*****************************************************************************80
!
!! P15_N returns the number of variables for problem 15.
!
!  Discussion:
!
!    The number of variables may be any multiple of 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 4

  return
end
subroutine p15_sol ( n, know, x )

!*****************************************************************************80
!
!! P15_SOL returns the solution for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p15_start ( n, x )

!*****************************************************************************80
!
!! P15_START returns a starting point for optimization for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( mod ( i, 4 ) == 1 ) then
      x(i) = 3.0D+00
    else if ( mod ( i, 4 ) == 2 ) then
      x(i) = - 1.0D+00
    else if ( mod ( i, 4 ) == 3 ) then
      x(i) = 0.0D+00
    else
      x(i) = 1.0D+00
    end if

  end do

  return
end
subroutine p15_title ( title )

!*****************************************************************************80
!
!! P15_TITLE returns a title for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Extended Powell Singular Quartic Function.'

  return
end
subroutine p16_f ( n, x, f )

!*****************************************************************************80
!
!! P16_F evaluates the objective function for problem 16.
!
!  Discussion:
!
!    This function has a valley approaching the line X(2) = 1.
!
!    The function has a global minimizer:
!
!      X(*) = ( 3.0, 0.5 ), F(X*) = 0.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Evelyn Beale,
!    On an Iterative Method for Finding a Local Minimum of a Function
!    of More than One Variable,
!    Technical Report 25, 
!    Statistical Techniques Research Group,
!    Princeton University, 1958.
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) x(n)

  f1 = 1.5D+00   - x(1) * ( 1.0D+00 - x(2)    )
  f2 = 2.25D+00  - x(1) * ( 1.0D+00 - x(2) * x(2) )
  f3 = 2.625D+00 - x(1) * ( 1.0D+00 - x(2) * x(2) * x(2) )

  f = f1 * f1 + f2 * f2 + f3 * f3

  return
end
subroutine p16_g ( n, x, g )

!*****************************************************************************80
!
!! P16_G evaluates the gradient for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df1dx2
  real ( kind = 8 ) df2dx1
  real ( kind = 8 ) df2dx2
  real ( kind = 8 ) df3dx1
  real ( kind = 8 ) df3dx2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  f1 = 1.5D+00   - x(1) * ( 1.0D+00 - x(2)    )
  f2 = 2.25D+00  - x(1) * ( 1.0D+00 - x(2) * x(2) )
  f3 = 2.625D+00 - x(1) * ( 1.0D+00 - x(2) * x(2) * x(2) )

  df1dx1 = - ( 1.0D+00 - x(2) )
  df1dx2 = x(1)
  df2dx1 = - ( 1.0D+00 - x(2) * x(2) )
  df2dx2 = 2.0D+00 * x(1) * x(2)
  df3dx1 = - ( 1.0D+00 - x(2) * x(2) * x(2) )
  df3dx2 = 3.0D+00 * x(1) * x(2) * x(2)

  g(1) = 2.0D+00 * ( f1 * df1dx1 + f2 * df2dx1 + f3 * df3dx1 )
  g(2) = 2.0D+00 * ( f1 * df1dx2 + f2 * df2dx2 + f3 * df3dx2 )

  return
end
subroutine p16_h ( n, x, h )

!*****************************************************************************80
!
!! P16_H evaluates the Hessian for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d2f1dx12
  real ( kind = 8 ) d2f1dx21
  real ( kind = 8 ) d2f2dx12
  real ( kind = 8 ) d2f2dx21
  real ( kind = 8 ) d2f2dx22
  real ( kind = 8 ) d2f3dx12
  real ( kind = 8 ) d2f3dx21
  real ( kind = 8 ) d2f3dx22
  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df1dx2
  real ( kind = 8 ) df2dx1
  real ( kind = 8 ) df2dx2
  real ( kind = 8 ) df3dx1
  real ( kind = 8 ) df3dx2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  f1 = 1.5D+00   - x(1) * ( 1.0D+00 - x(2)    )
  f2 = 2.25D+00  - x(1) * ( 1.0D+00 - x(2) * x(2) )
  f3 = 2.625D+00 - x(1) * ( 1.0D+00 - x(2) * x(2) * x(2) )

  df1dx1 = - ( 1.0D+00 - x(2) )
  df1dx2 = x(1)
  df2dx1 = - ( 1.0D+00 - x(2) * x(2) )
  df2dx2 = 2.0D+00 * x(1) * x(2)
  df3dx1 = - ( 1.0D+00 - x(2) * x(2) * x(2) )
  df3dx2 = 3.0D+00 * x(1) * x(2) * x(2)

  d2f1dx12 = 1.0D+00
  d2f1dx21 = 1.0D+00
 
  d2f2dx12 = 2.0D+00 * x(2)
  d2f2dx21 = 2.0D+00 * x(2)
  d2f2dx22 = 2.0D+00 * x(1)

  d2f3dx12 = 3.0D+00 * x(2) * x(2)
  d2f3dx21 = 3.0D+00 * x(2) * x(2)
  d2f3dx22 = 6.0D+00 * x(1) * x(2)

  h(1,1) = 2.0D+00 * ( df1dx1 * df1dx1 &
                     + df2dx1 * df2dx1 &
                     + df3dx1 * df3dx1 )

  h(1,2) = 2.0D+00 * ( df1dx2 * df1dx1 + f1 * d2f1dx12 &
                     + df2dx2 * df2dx1 + f2 * d2f2dx12 &
                     + df3dx2 * df3dx1 + f3 * d2f3dx12 )

  h(2,1) = 2.0D+00 * ( df1dx1 * df1dx2 + f1 * d2f1dx21 &
                     + df2dx1 * df2dx2 + f2 * d2f2dx21 &
                     + df3dx1 * df3dx2 + f3 * d2f3dx21 )

  h(2,2) = 2.0D+00 * ( df1dx2 * df1dx2 &
                     + df2dx2 * df2dx2 + f2 * d2f2dx22 &
                     + df3dx2 * df3dx2 + f3 * d2f3dx22 )

  return
end
subroutine p16_n ( n )

!*****************************************************************************80
!
!! P16_N returns the number of variables for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p16_sol ( n, know, x )

!*****************************************************************************80
!
!! P16_SOL returns the solution for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:2) = (/ 3.0D+00, 0.5D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p16_start ( n, x )

!*****************************************************************************80
!
!! P16_START returns a starting point for optimization for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ 1.0D+00, 1.0D+00 /)

  return
end
subroutine p16_title ( title )

!*****************************************************************************80
!
!! P16_TITLE returns a title for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Beale Function.'

  return
end
subroutine p17_f ( n, x, f )

!*****************************************************************************80
!
!! P17_F evaluates the objective function for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  real ( kind = 8 ) f5
  real ( kind = 8 ) f6
  real ( kind = 8 ) x(n)

  f1 = x(2) - x(1) * x(1)
  f2 = 1.0D+00 - x(1)
  f3 = x(4) - x(3) * x(3)
  f4 = 1.0D+00 - x(3)
  f5 = x(2) + x(4) - 2.0D+00
  f6 = x(2) - x(4)

  f = 100.0D+00 * f1 * f1 &
    +             f2 * f2 &
    +  90.0D+00 * f3 * f3 &
    +             f4 * f4 &
    +  10.0D+00 * f5 * f5 &
    +   0.1D+00 * f6 * f6

  return
end
subroutine p17_g ( n, x, g )

!*****************************************************************************80
!
!! P17_G evaluates the gradient for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df1dx2
  real ( kind = 8 ) df2dx1
  real ( kind = 8 ) df3dx3
  real ( kind = 8 ) df3dx4
  real ( kind = 8 ) df4dx3
  real ( kind = 8 ) df5dx2
  real ( kind = 8 ) df5dx4
  real ( kind = 8 ) df6dx2
  real ( kind = 8 ) df6dx4
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  real ( kind = 8 ) f5
  real ( kind = 8 ) f6
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  f1 = x(2) - x(1) * x(1)
  f2 = 1.0D+00 - x(1)
  f3 = x(4) - x(3)**2
  f4 = 1.0D+00 - x(3)
  f5 = x(2) + x(4) - 2.0D+00
  f6 = x(2) - x(4)

  df1dx1 = - 2.0D+00 * x(1)
  df1dx2 =   1.0D+00
  df2dx1 = - 1.0D+00
  df3dx3 = - 2.0D+00 * x(3)
  df3dx4 =   1.0D+00
  df4dx3 = - 1.0D+00
  df5dx2 =   1.0D+00
  df5dx4 =   1.0D+00
  df6dx2 =   1.0D+00
  df6dx4 = - 1.0D+00

  g(1) = 2.0D+00 * ( 100.0D+00 * f1 * df1dx1 + f2 * df2dx1 )

  g(2) = 2.0D+00 * ( 100.0D+00 * f1 * df1dx2 + 10.0D+00 * f5 * df5dx2 &
    +  0.1D+00 * f6 * df6dx2 )

  g(3) = 2.0D+00 * ( 90.0D+00 * f3 * df3dx3 + f4 * df4dx3 )

  g(4) = 2.0D+00 * ( 90.0D+00 * f3 * df3dx4 + 10.0D+00 * f5 * df5dx4 &
    + 0.1D+00 * f6 * df6dx4 )

  return
end
subroutine p17_h ( n, x, h )

!*****************************************************************************80
!
!! P17_H evaluates the Hessian for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d2f1dx11
  real ( kind = 8 ) d2f3dx33
  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df1dx2
  real ( kind = 8 ) df2dx1
  real ( kind = 8 ) df3dx3
  real ( kind = 8 ) df3dx4
  real ( kind = 8 ) df4dx3
  real ( kind = 8 ) df5dx2
  real ( kind = 8 ) df5dx4
  real ( kind = 8 ) df6dx2
  real ( kind = 8 ) df6dx4
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  real ( kind = 8 ) f5
  real ( kind = 8 ) f6
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  f1 = x(2) - x(1)**2
  f2 = 1.0D+00 - x(1)
  f3 = x(4) - x(3)**2
  f4 = 1.0D+00 - x(3)
  f5 = x(2) + x(4) - 2.0D+00
  f6 = x(2) - x(4)

  df1dx1 = - 2.0D+00 * x(1)
  df1dx2 =   1.0D+00
  df2dx1 = - 1.0D+00
  df3dx3 = - 2.0D+00 * x(3)
  df3dx4 =   1.0D+00
  df4dx3 = - 1.0D+00
  df5dx2 =   1.0D+00
  df5dx4 =   1.0D+00
  df6dx2 =   1.0D+00
  df6dx4 = - 1.0D+00

  d2f1dx11 = - 2.0D+00
  d2f3dx33 = - 2.0D+00

  h(1,1) = 200.0D+00 * df1dx1 * df1dx1 + 200.0D+00 * f1 * d2f1dx11 &
         + 2.0D+00 * df2dx1 * df2dx1
  h(1,2) = 200.0D+00 * df1dx2 * df1dx1  
  h(1,3) = 0.0D+00
  h(1,4) = 0.0D+00

  h(2,1) = 200.0D+00 * df1dx1 * df1dx2
  h(2,2) = 200.0D+00 * df1dx2 * df1dx2 + 20.0D+00 * df5dx2 * df5dx2 &
         + 0.2D+00 * df6dx2 * df6dx2 
  h(2,3) = 0.0D+00
  h(2,4) = 20.0D+00 * df5dx4 * df5dx2 + 0.2D+00 * df6dx4 * df6dx2

  h(3,1) = 0.0D+00
  h(3,2) = 0.0D+00
  h(3,3) = 180.0D+00 * df3dx3 * df3dx3 + 180.0D+00 * f3 * d2f3dx33 &
         + 2.0D+00 * df4dx3 * df4dx3
  h(3,4) = 180.0D+00 * df3dx4 * df3dx3

  h(4,1) = 0.0D+00
  h(4,2) = 20.0D+00 * df5dx2 * df5dx4 + 0.2D+00 * df6dx2 * df6dx4
  h(4,3) = 180.0D+00 * df3dx3 * df3dx4
  h(4,4) = 180.0D+00 * df3dx4 * df3dx4 + 20.0D+00 * df5dx4 * df5dx4 &
         + 0.2D+00 * df6dx4 * df6dx4

  return
end
subroutine p17_n ( n )

!*****************************************************************************80
!
!! P17_N returns the number of variables for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 4

  return
end
subroutine p17_sol ( n, know, x )

!*****************************************************************************80
!
!! P17_SOL returns the solution for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:4) = (/ 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p17_start ( n, x )

!*****************************************************************************80
!
!! P17_START returns a starting point for optimization for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:4) = (/ -3.0D+00, -1.0D+00, -3.0D+00, -1.0D+00 /)

  return
end
subroutine p17_title ( title )

!*****************************************************************************80
!
!! P17_TITLE returns a title for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Wood Function.'

  return
end
subroutine p18_f ( n, x, f )

!*****************************************************************************80
!
!! P18_F evaluates the objective function for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) fvec(n)
  real ( kind = 8 ) x(n)
!
!  Compute FVEC.
!
  call p18_fvec ( n, x, fvec )
!
!  Compute F.
!
  f = sum ( fvec(1:n)**2 )

  return
end
subroutine p18_fvec ( n, x, fvec )

!*****************************************************************************80
!
!! P18_FVEC is an auxilliary routine for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) FVEC(N), an auxilliary vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)

  fvec(1:n) = 0.0D+00

  do j = 1, n
    t1 = 1.0D+00
    t2 = 2.0D+00 * x(j) - 1.0D+00
    t = 2.0D+00 * t2
    do i = 1, n
      fvec(i) = fvec(i) + t2
      th = t * t2 - t1
      t1 = t2
      t2 = th
    end do
  end do

  do i = 1, n
    fvec(i) = fvec(i) / real ( n, kind = 8 )
    if ( mod ( i, 2 ) == 0 ) then
      fvec(i) = fvec(i) + 1.0D+00 / ( real ( i, kind = 8 )**2 - 1.0D+00 )
    end if
  end do

  return
end
subroutine p18_g ( n, x, g )

!*****************************************************************************80
!
!! P18_G evaluates the gradient for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fvec(n)
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)
!
!  Compute FVEC.
!
  call p18_fvec ( n, x, fvec )
!
!  Compute G.
!
  do j = 1, n
    g(j) = 0.0D+00
    t1 = 1.0D+00
    t2 = 2.0D+00 * x(j) - 1.0D+00
    t = 2.0D+00 * t2
    s1 = 0.0D+00
    s2 = 2.0D+00
    do i = 1, n
      g(j) = g(j) + fvec(i) * s2
      th = 4.0D+00 * t2 + t * s2 - s1
      s1 = s2
      s2 = th
      th = t * t2 - t1
      t1 = t2
      t2 = th
    end do
  end do

  g(1:n) = 2.0D+00 * g(1:n) / real ( n, kind = 8 )

  return
end
subroutine p18_h ( n, x, h )

!*****************************************************************************80
!
!! P18_H evaluates the Hessian for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) fvec(n)
  real ( kind = 8 ) gvec(n)
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) ss1
  real ( kind = 8 ) ss2
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) th
  real ( kind = 8 ) tt
  real ( kind = 8 ) tth
  real ( kind = 8 ) tt1
  real ( kind = 8 ) tt2
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  call p18_fvec ( n, x, fvec )

  d1 = 1.0D+00 / real ( n, kind = 8 )
  d2 = 2.0D+00 * d1

  do j = 1, n

    h(j,j) = 4.0D+00 * d1
    t1 = 1.0D+00
    t2 = 2.0D+00 * x(j) - 1.0D+00
    t = 2.0D+00 * t2
    s1 = 0.0D+00
    s2 = 2.0D+00
    p1 = 0.0D+00
    p2 = 0.0D+00
    gvec(1) = s2

    do i = 2, n
      th = 4.0D+00 * t2 + t * s2 - s1
      s1 = s2
      s2 = th
      th = t * t2 - t1
      t1 = t2
      t2 = th
      th = 8.0D+00 * s1 + t * p2 - p1
      p1 = p2
      p2 = th
      gvec(i) = s2
      h(j,j) = h(j,j) + fvec(i) * th + d1 * s2 * s2
    end do

    h(j,j) = d2 * h(j,j)

    do k = 1, j-1

      h(j,k) = 0.0
      tt1 = 1.0D+00
      tt2 = 2.0D+00 * x(k) - 1.0D+00
      tt = 2.0D+00 * tt2
      ss1 = 0.0D+00
      ss2 = 2.0D+00

      do i = 1, n
        h(j,k) = h(j,k) + ss2 * gvec(i)
        tth = 4.0D+00 * tt2 + tt * ss2 - ss1
        ss1 = ss2
        ss2 = tth
        tth = tt * tt2 - tt1
        tt1 = tt2
        tt2 = tth
      end do

      h(j,k) = d2 * d1 * h(j,k)

    end do

  end do

  do i = 1, n
    h(i,i+1:n) = h(i+1:n,i)
  end do

  return
end
subroutine p18_n ( n )

!*****************************************************************************80
!
!! P18_N returns the number of variables for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p18_sol ( n, know, x )

!*****************************************************************************80
!
!! P18_SOL returns the solution for problem 18.
!
!  Discussion:
!
!    The solution values are taken from Brent.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 .and. n == 2 ) then
    know = 1
    x(1:2) = (/ 0.2113249D+00, 0.7886751D+00 /)
  else if ( know == 0 .and. n == 4 ) then
    know = 1
    x(1:4) = (/ 0.1026728D+00, 0.4062037D+00, 0.5937963D+00, 0.8973272D+00 /)
  else if ( know == 0 .and. n == 6 ) then
    know = 1
    x(1:6) = (/ 0.066877D+00, 0.288741D+00, 0.366682D+00, 0.633318D+00, &
      0.711259D+00, 0.933123D+00 /)
  else if ( know == 0 .and. n == 8 ) then
    know = 1
    x(1:8) = (/ 0.043153D+00, 0.193091D+00, 0.266329D+00, 0.500000D+00, &
      0.500000D+00, 0.733671D+00, 0.806910D+00, 0.956847D+00 /)
  else
    know = 0
    x(1:n) = 0.0D+00
  end if

  return
end
subroutine p18_start ( n, x )

!*****************************************************************************80
!
!! P18_START returns a starting point for optimization for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = real ( i, kind = 8 ) / real ( n + 1, kind = 8 )
  end do

  return
end
subroutine p18_title ( title )

!*****************************************************************************80
!
!! P18_TITLE returns a title for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Chebyquad Function'

  return
end
subroutine p19_f ( n, x, f )

!*****************************************************************************80
!
!! P19_F evaluates the objective function for problem 19.
!
!  Discussion:
!
!    The function is similar to Rosenbrock's.  The "valley" follows
!    the curve X(2) = X(1)**3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    A Leon,
!    A Comparison of Eight Known Optimizing Procedures,
!    in Recent Advances in Optimization Techniques,
!    edited by Abraham Lavi, Thomas Vogl,
!    Wiley, 1966.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) x(n)

  f1 = x(2) - x(1) * x(1) * x(1)
  f2 = 1.0D+00 - x(1)

  f = 100.0D+00 * f1 * f1 &
    +             f2 * f2

  return
end
subroutine p19_g ( n, x, g )

!*****************************************************************************80
!
!! P19_G evaluates the gradient for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  g(1) = - 600.0D+00 * ( x(2) - x(1)**3 ) * x(1) * x(1) &
    - 2.0D+00 * ( 1.0D+00 - x(1) )

  g(2) = 200.0D+00 * ( x(2) - x(1)**3 )

  return
end
subroutine p19_h ( n, x, h )

!*****************************************************************************80
!
!! P19_H evaluates the Hessian for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1,1) = - 1200.0D+00 * x(1) * x(2) + 3000.0D+00 * x(1)**4 + 2.0D+00
  h(1,2) = - 600.0D+00 * x(1) * x(1)

  h(2,1) = - 600.0D+00 * x(1) * x(1)
  h(2,2) = 200.0D+00

  return
end
subroutine p19_n ( n )

!*****************************************************************************80
!
!! P19_N returns the number of variables for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p19_sol ( n, know, x )

!*****************************************************************************80
!
!! P19_SOL returns the solution for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:2) = (/ 1.0D+00, 1.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p19_start ( n, x )

!*****************************************************************************80
!
!! P19_START returns a starting point for optimization for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ -1.2D+00, -1.0D+00 /)

  return
end
subroutine p19_title ( title )

!*****************************************************************************80
!
!! P19_TITLE returns a title for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Leon cubic valley function'

  return
end
subroutine p20_f ( n, x, f )

!*****************************************************************************80
!
!! P20_F evaluates the objective function for problem 20.
!
!  Discussion:
!
!    The function has the form 
!      f = x'*A*x - 2*x(1)
!    where A is the (-1,2,-1) tridiagonal matrix, except that A(1,1) is 1.  
!    The minimum value of F(X) is -N, which occurs for
!      x = ( n, n-1, ..., 2, 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  f = x(1) * x(1) + 2.0D+00 * sum ( x(2:n)**2 )

  do i = 1, n-1
    f = f - 2.0D+00 * x(i) * x(i+1)
  end do

  f = f - 2.0D+00 * x(1)

  return
end
subroutine p20_g ( n, x, g )

!*****************************************************************************80
!
!! P20_G evaluates the gradient for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( i == 1 ) then
      g(i) = x(i)
    else
      g(i) = 2.0D+00 * x(i)
    end if

    if ( 1 < i ) then
      g(i) = g(i) - x(i-1)
    end if

    if ( i < n ) then
      g(i) = g(i) - x(i+1)
    end if

  end do

  g(1) = g(1) - 2.0D+00

  return
end
subroutine p20_h ( n, x, h )

!*****************************************************************************80
!
!! P20_H evaluates the Hessian for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  do i = 1, n
    do j = 1, n

      if ( i == j ) then
        if ( i == 1 ) then
          h(i,j) = 2.0D+00
        else
          h(i,j) = 4.0D+00
        end if
      else if ( i == j + 1 .or. i == j - 1 ) then
        h(i,j) = - 2.0D+00
      else
        h(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine p20_n ( n )

!*****************************************************************************80
!
!! P20_N returns the number of variables for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p20_sol ( n, know, x )

!*****************************************************************************80
!
!! P20_SOL returns the solution for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    do i = 1, n
      x(i) = real ( n + 1 - i, kind = 8 )
    end do
  else
    know = 0
  end if

  return
end
subroutine p20_start ( n, x )

!*****************************************************************************80
!
!! P20_START returns a starting point for optimization for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = 0.0D+00

  return
end
subroutine p20_title ( title )

!*****************************************************************************80
!
!! P20_TITLE returns a title for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Gregory and Karney Tridiagonal Matrix Function'

  return
end
subroutine p21_f ( n, x, f )

!*****************************************************************************80
!
!! P21_F evaluates the objective function for problem 21.
!
!  Discussion:
!
!    The function has the form 
!      f = x'*A*x
!    where A is the Hilbert matrix.  The minimum value
!    of F(X) is 0, which occurs for
!      x = ( 0, 0, ..., 0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  f = 0.0D+00

  do i = 1, n
    do j = 1, n
      f = f + x(i) * x(j) / real ( i + j - 1, kind = 8 )
    end do
  end do

  return
end
subroutine p21_g ( n, x, g )

!*****************************************************************************80
!
!! P21_G evaluates the gradient for problem 21.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  do i = 1, n

    g(i) = 0.0D+00
    do j = 1, n
      g(i) = g(i) + 2.0D+00 * x(j) / real ( i + j - 1, kind = 8 )
    end do

  end do

  return
end
subroutine p21_h ( n, x, h )

!*****************************************************************************80
!
!! P21_H evaluates the Hessian for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  do i = 1, n
    do j = 1, n
      h(i,j) = 2.0D+00 / real ( i + j - 1, kind = 8 )
    end do
  end do

  return
end
subroutine p21_n ( n )

!*****************************************************************************80
!
!! P21_N returns the number of variables for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p21_sol ( n, know, x )

!*****************************************************************************80
!
!! P21_SOL returns the solution for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p21_start ( n, x )

!*****************************************************************************80
!
!! P21_START returns a starting point for optimization for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = 1.0D+00

  return
end
subroutine p21_title ( title )

!*****************************************************************************80
!
!! P21_TITLE returns a title for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Hilbert Matrix Function F = x''Ax'

  return
end
subroutine p22_f ( n, x, f )

!*****************************************************************************80
!
!! P22_F evaluates the objective function for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) x(n)

  f = sum ( x(1:n)**2 )

  return
end
subroutine p22_g ( n, x, g )

!*****************************************************************************80
!
!! P22_G evaluates the gradient for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  g(1:n) = 2.0D+00 * x(1:n)

  return
end
subroutine p22_h ( n, x, h )

!*****************************************************************************80
!
!! P22_H evaluates the Hessian for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do i = 1, n
    h(i,i) = 2.0D+00
  end do

  return
end
subroutine p22_n ( n )

!*****************************************************************************80
!
!! P22_N returns the number of variables for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p22_sol ( n, know, x )

!*****************************************************************************80
!
!! P22_SOL returns the solution for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p22_start ( n, x )

!*****************************************************************************80
!
!! P22_START returns a starting point for optimization for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: xhi = 5.12D+00
  real ( kind = 8 ), parameter :: xlo = -5.12D+00

  call r8vec_even ( xlo, xhi, n, x )

  return
end
subroutine p22_title ( title )

!*****************************************************************************80
!
!! P22_TITLE returns a title for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The De Jong Function F1'

  return
end
subroutine p23_f ( n, x, f )

!*****************************************************************************80
!
!! P23_F evaluates the objective function for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) x(n)

  f = 100.0D+00 * ( x(1) * x(1) - x(2) )**2 + ( 1.0D+00 - x(1) )**2

  return
end
subroutine p23_g ( n, x, g )

!*****************************************************************************80
!
!! P23_G evaluates the gradient for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  g(1) = 400.0D+00 * ( x(1) * x(1) - x(2) ) * x(1) - 2.0D+00 + 2.0D+00 * x(1)
  g(2) = - 200.0D+00 * ( x(1) * x(1) - x(2) )

  return
end
subroutine p23_h ( n, x, h )

!*****************************************************************************80
!
!! P23_H evaluates the Hessian for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1,1) = 1200.0D+00 * x(1) * x(1) - 400.0D+00 * x(2) + 2.0D+00
  h(1,2) = - 400.0D+00 * x(1)
  h(2,1) = -400.0D+00 * x(1)
  h(2,2) = 200.0D+00

  return
end
subroutine p23_n ( n )

!*****************************************************************************80
!
!! P23_N returns the number of variables for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p23_sol ( n, know, x )

!*****************************************************************************80
!
!! P23_SOL returns the solution for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = 1.0D+00
  else
    know = 0
  end if

  return
end
subroutine p23_start ( n, x )

!*****************************************************************************80
!
!! P23_START returns a starting point for optimization for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: xhi =  2.048D+00
  real ( kind = 8 ), parameter :: xlo = -2.048D+00

  call r8vec_even ( xlo, xhi, n, x )

  return
end
subroutine p23_title ( title )

!*****************************************************************************80
!
!! P23_TITLE returns a title for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The De Jong Function F2'

  return
end
subroutine p24_f ( n, x, f )

!*****************************************************************************80
!
!! P24_F evaluates the objective function for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) x(n)

  f = real ( sum ( int ( x(1:n) ) ), kind = 8 )

  return
end
subroutine p24_g ( n, x, g )

!*****************************************************************************80
!
!! P24_G evaluates the gradient for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  return
end
subroutine p24_h ( n, x, h )

!*****************************************************************************80
!
!! P24_H evaluates the Hessian for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  return
end
subroutine p24_n ( n )

!*****************************************************************************80
!
!! P24_N returns the number of variables for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 5

  return
end
subroutine p24_sol ( n, know, x )

!*****************************************************************************80
!
!! P24_SOL returns the solution for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = -5.0D+00
  else
    know = 0
  end if

  return
end
subroutine p24_start ( n, x )

!*****************************************************************************80
!
!! P24_START returns a starting point for optimization for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: xhi = 5.12D+00
  real ( kind = 8 ), parameter :: xlo = -5.12D+00

  call r8vec_even ( xlo, xhi, n, x )

  return
end
subroutine p24_title ( title )

!*****************************************************************************80
!
!! P24_TITLE returns a title for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The De Jong Function F3, (discontinuous)'

  return
end
subroutine p25_f ( n, x, f )

!*****************************************************************************80
!
!! P25_F evaluates the objective function for problem 25.
!
!  Discussion:
!
!    The function includes Gaussian noise, multiplied by a parameter P.
!
!    If P is zero, then the function is a proper function, and evaluating
!    it twice with the same argument will yield the same results.
!    Moreover, P25_G and P25_H are the correct gradient and hessian routines.
!
!    If P is nonzero, then evaluating the function at the same argument
!    twice will generally yield two distinct values; this means the function
!    is not even a well defined mathematical function, let alone continuous;
!    the gradient and hessian are not correct.  And yet, at least for small
!    values of P, it may be possible to approximate the minimizer of the 
!    (underlying well-defined ) function.
!
!    The value of the parameter P is by default 1.  The user can manipulate
!    this value by calling P25_P_GET or P25_P_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) gauss
  integer ( kind = 4 ) i
  real ( kind = 8 ) p
  real ( kind = 8 ) x(n)

  call p25_p_get ( p )

  call normal_01_sample ( gauss )

  f = p * gauss
  do i = 1, n
    f = f + real ( i, kind = 8 ) * x(i)**4
  end do

  return
end
subroutine p25_g ( n, x, g )

!*****************************************************************************80
!
!! P25_G evaluates the gradient for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    g(i) = real ( i, kind = 8 ) * 4.0D+00 * x(i)**3
  end do

  return
end
subroutine p25_h ( n, x, h )

!*****************************************************************************80
!
!! P25_H evaluates the Hessian for problem 25.
!
!  Discussion:
!
!    Note that, for P = 0, the Hessian matrix should be diagonal.
!    However, if it is estimated using finite differences, off diagonal
!    terms may appear.  This occurs when the argument increment dX is
!    too small to be significant in terms such as X^6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do i = 1, n
    h(i,i) = real ( i, kind = 8 ) * 12.0D+00 * x(i)**2
  end do

  return
end
subroutine p25_n ( n )

!*****************************************************************************80
!
!! P25_N returns the number of variables for problem 25.
!
!  Discussion:
!
!    The function is actually well defined for any positive value of N.
!    The value given here is that specified in the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 30

  return
end
subroutine p25_p_get ( p )

!*****************************************************************************80
!
!! P25_P_GET gets the value of a parameter for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P, the value of the parameter that multiplies the
!    Gaussian noise.
!
  implicit none

  real ( kind = 8 ) p

  call p25_p_val ( 'GET', p )

  return
end
subroutine p25_p_set ( p )

!*****************************************************************************80
!
!! P25_P_SET sets the value of a parameter for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the parameter that multiplies the 
!    Gaussian noise.
!
  implicit none

  real ( kind = 8 ) p

  call p25_p_val ( 'SET', p )

  return
end
subroutine p25_p_val ( action, p )

!*****************************************************************************80
!
!! P25_P_VAL sets or gets the value of a parameter for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.  
!    'S', set the interval value to the input value of P.
!    'G', set the output value of P to the internal value.
!
!    Input/output, real ( kind = 8 ) P, the value of the parameter that
!    multiplies the Gaussian noise.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ) p
  real ( kind = 8 ), save :: p_save = 1.0D+00

  if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then
    p_save = p
  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then
    p = p_save
  end if

  return
end
subroutine p25_sol ( n, know, x )

!*****************************************************************************80
!
!! P25_SOL returns the solution for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p25_start ( n, x )

!*****************************************************************************80
!
!! P25_START returns a starting point for optimization for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: xhi = 1.28D+00
  real ( kind = 8 ), parameter :: xlo = -1.28D+00

  call r8vec_even ( xlo, xhi, n, x )

  return
end
subroutine p25_title ( title )

!*****************************************************************************80
!
!! P25_TITLE returns a title for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The De Jong Function F4 (with Gaussian noise)'

  return
end
subroutine p26_f ( n, x, f )

!*****************************************************************************80
!
!! P26_F evaluates the objective function for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1
  integer ( kind = 4 ) a2
  real ( kind = 8 ) f
  real ( kind = 8 ) fi
  real ( kind = 8 ) fj
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ), parameter :: jroot = 5
  integer ( kind = 4 ), parameter :: k = 500
  real ( kind = 8 ) x(n)

  fi = real ( k, kind = 8 )

  do j = 1, jroot * jroot

    j1 = mod ( j - 1, jroot ) + 1
    a1 = - 32 + j1 * 16

    j2 = ( j - 1 ) / jroot
    a2 = - 32 + j2 * 16

    fj = real ( j, kind = 8 ) + ( x(1) - real ( a1, kind = 8 ) )**6 &
      + ( x(2) - real ( a2, kind = 8 ) )**6

    fi = fi + 1.0D+00 / fj

  end do

  f = 1.0D+00 / fi

  return
end
subroutine p26_g ( n, x, g )

!*****************************************************************************80
!
!! P26_G evaluates the gradient for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ), parameter :: jroot = 5
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1
  integer ( kind = 4 ) a2
  real ( kind = 8 ) dfidx1
  real ( kind = 8 ) dfidx2
  real ( kind = 8 ) fi
  real ( kind = 8 ) fj
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ), parameter :: k = 500
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  fi = real ( k, kind = 8 )
  dfidx1 = 0.0D+00
  dfidx2 = 0.0D+00

  do j = 1, jroot * jroot

    j1 = mod ( j-1, jroot ) + 1
    a1 = -32 + j1 * 16

    j2 = ( j - 1 ) / jroot
    a2 = -32 + j2 * 16

    fj = ( real ( j, kind = 8 ) + ( x(1) - real ( a1, kind = 8 ) )**6 &
      + ( x(2) - real ( a2, kind = 8 ) )**6 )

    fi = fi + 1.0D+00 / fj

    dfidx1 = dfidx1 - ( 1.0D+00 / fj**2 ) * 6.0D+00 &
      * ( x(1) - real ( a1, kind = 8 ) )**5
    dfidx2 = dfidx2 - ( 1.0D+00 / fj**2 ) * 6.0D+00 &
      * ( x(2) - real ( a2, kind = 8 ) )**5

  end do

  g(1) = - ( 1.0D+00 / fi**2 ) * dfidx1
  g(2) = - ( 1.0D+00 / fi**2 ) * dfidx2

  return
end
subroutine p26_h ( n, x, h )

!*****************************************************************************80
!
!! P26_H evaluates the Hessian for problem 26.
!
!  Discussion:
!
!    I haven't properly written this routine yet.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  return
end
subroutine p26_n ( n )

!*****************************************************************************80
!
!! P26_N returns the number of variables for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p26_sol ( n, know, x )

!*****************************************************************************80
!
!! P26_SOL returns the solution for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none
  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ -32.0D+00, -32.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p26_start ( n, x )

!*****************************************************************************80
!
!! P26_START returns a starting point for optimization for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ -32.01D+00, -32.02D+00 /)

  return
end
subroutine p26_title ( title )

!*****************************************************************************80
!
!! P26_TITLE returns a title for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The De Jong Function F5'

  return
end
subroutine p27_f ( n, x, f )

!*****************************************************************************80
!
!! P27_F evaluates the objective function for problem 27.
!
!  Discussion:
!
!    F can be regarded as a function of R = SQRT ( X(1)^2 + X(2)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) f
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)

  r = sqrt ( x(1)**2 + x(2)**2 )

  a = ( 1.0D+00 + 0.001D+00 * r**2 )**( -2 )

  b = ( sin ( r ) )**2 - 0.5D+00

  f = 0.5D+00 + a * b

  return
end
subroutine p27_g ( n, x, g )

!*****************************************************************************80
!
!! P27_G evaluates the gradient for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) ar
  real ( kind = 8 ) b
  real ( kind = 8 ) br
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) r
  real ( kind = 8 ) rx1
  real ( kind = 8 ) rx2
  real ( kind = 8 ) x(n)

  r = sqrt ( x(1)**2 + x(2)**2 )

  if ( r == 0.0D+00 ) then
    g(1) = 0.0D+00
    g(2) = 0.0D+00
    return
  end if

  rx1 = x(1) / r
  rx2 = x(2) / r

  a = ( 1.0D+00 + 0.001D+00 * r**2 )**( -2 )
  ar = - 0.004D+00 * r * ( 1.0D+00 + 0.001D+00 * r**2 )**( -3 )

  b = ( sin ( r ) )**2 - 0.5D+00
  br = sin ( 2.0D+00 * r )

  g(1) = ( ar * b + a * br ) * rx1
  g(2) = ( ar * b + a * br ) * rx2

  return
end
subroutine p27_h ( n, x, h )

!*****************************************************************************80
!
!! P27_H evaluates the Hessian for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) ar
  real ( kind = 8 ) arr
  real ( kind = 8 ) b
  real ( kind = 8 ) br
  real ( kind = 8 ) brr
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) r
  real ( kind = 8 ) rx1
  real ( kind = 8 ) rx1x1
  real ( kind = 8 ) rx1x2
  real ( kind = 8 ) rx2
  real ( kind = 8 ) rx2x1
  real ( kind = 8 ) rx2x2
  real ( kind = 8 ) x(n)

  r = sqrt ( x(1)**2 + x(2)**2 )

  if ( r == 0.0D+00 ) then
    h(1,1) = 1.0D+00
    h(1,2) = 0.0D+00
    h(2,1) = 0.0D+00
    h(2,2) = 1.0D+00
    return
  end if

  rx1 = x(1) / r
  rx2 = x(2) / r

  rx1x1 = x(2)**2 / r**3
  rx1x2 = - x(1) * x(2) / r**3
  rx2x1 = - x(1) * x(2) / r**3
  rx2x2 = x(1)**2 / r**3

  a = ( 1.0D+00 + 0.001D+00 * r**2 )**( -2 )
  ar = - 0.004D+00 * r * ( 1.0D+00 + 0.001D+00 * r**2 )**( -3 )
  arr = - 0.004D+00 * ( 1.0D+00 + 0.001D+00 * r**2 )**( -3 ) &
    + 0.000024D+00 * r * ( 1.0D+00 + 0.001D+00 * r**2 )**( -4 )

  b = ( sin ( r ) )**2 - 0.5D+00
  br = sin ( 2.0D+00 * r )
  brr = 2.0D+00 * cos ( 2.0D+00 * r )

  h(1,1) = ( arr * b + 2.0D+00 * ar * br + a * brr ) * rx1 * rx1 &
    + ( ar * b + a * br ) * rx1x1

  h(1,2) = ( arr * b + 2.0D+00 * ar * br + a * brr ) * rx1 * rx2 &
    + ( ar * b + a * br ) * rx1x2

  h(2,1) = ( arr * b + 2.0D+00 * ar * br + a * brr ) * rx2 * rx1 &
    + ( ar * b + a * br ) * rx2x1

  h(2,2) = ( arr * b + 2.0D+00 * ar * br + a * brr ) * rx2 * rx2 &
    + ( ar * b + a * br ) * rx2x2

  return
end
subroutine p27_n ( n )

!*****************************************************************************80
!
!! P27_N returns the number of variables for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p27_sol ( n, know, x )

!*****************************************************************************80
!
!! P27_SOL returns the solution for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 0.0D+00, 0.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p27_start ( n, x )

!*****************************************************************************80
!
!! P27_START returns a starting point for optimization for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ -5.0D+00, +10.0D+00 /)

  return
end
subroutine p27_title ( title )

!*****************************************************************************80
!
!! P27_TITLE returns a title for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Schaffer Function F6'

  return
end
subroutine p28_f ( n, x, f )

!*****************************************************************************80
!
!! P28_F evaluates the objective function for problem 28.
!
!  Discussion:
!
!    Note that F is a function of R^2 = X(1)^2 + X(2)^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)

  r = sqrt ( x(1)**2 + x(2)**2 )

  f = sqrt ( r ) * ( 1.0D+00 + ( sin ( 50.0D+00 * r**0.2D+00 ) )**2 )

  return
end
subroutine p28_g ( n, x, g )

!*****************************************************************************80
!
!! P28_G evaluates the gradient for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) ar
  real ( kind = 8 ) b
  real ( kind = 8 ) br
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) r
  real ( kind = 8 ) rx1
  real ( kind = 8 ) rx2
  real ( kind = 8 ) x(n)

  r = sqrt ( x(1)**2 + x(2)**2 )

  if ( r == 0.0D+00 ) then
    g(1) = 0.0D+00
    g(2) = 0.0D+00
    return
  end if

  a = sqrt ( r )
  ar = 0.5D+00 / sqrt ( r )

  b = 1.0D+00 + ( sin ( 50.0D+00 * r**0.2D+00 ) )**2
  br = 10.0D+00 * sin ( 100.0D+00 * r**0.2D+00 ) * r**(-0.8D+00)

  rx1 = x(1) / r
  rx2 = x(2) / r

  g(1) = ( ar * b + a * br ) * rx1
  g(2) = ( ar * b + a * br ) * rx2

  return
end
subroutine p28_h ( n, x, h )

!*****************************************************************************80
!
!! P28_H evaluates the Hessian for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) ar
  real ( kind = 8 ) arr
  real ( kind = 8 ) b
  real ( kind = 8 ) br
  real ( kind = 8 ) brr
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) r
  real ( kind = 8 ) rx1
  real ( kind = 8 ) rx1x1
  real ( kind = 8 ) rx1x2
  real ( kind = 8 ) rx2
  real ( kind = 8 ) rx2x1
  real ( kind = 8 ) rx2x2
  real ( kind = 8 ) x(n)

  r = sqrt ( x(1)**2 + x(2)**2 )

  rx1 = x(1) / r
  rx2 = x(2) / r

  rx1x1 = x(2)**2 / r**3
  rx1x2 = - x(1) * x(2) / r**3
  rx2x1 = - x(1) * x(2) / r**3
  rx2x2 = x(1)**2 / r**3
!
!  F = A * B
!  dFdX1 = ( Ar * B + A * Br ) * Rx1
!  d2FdX1dX1 = ( Arr * B + Ar * Br ) * Rx1**2 + ( Ar * B + A * Br ) * Rx1x1
!  etc
!
  a = sqrt ( r )
  ar = 0.5D+00 / sqrt ( r )
  arr = - 0.25D+00 / sqrt ( r**3 )

  b = 1.0D+00 + ( sin ( 50.0D+00 * r**0.2D+00 ) )**2
  br = 10.0D+00 * sin ( 100.0D+00 * r**0.2D+00 ) * r**(-0.8D+00)
  brr = 200.0D+00 * cos ( 100.0D+00 * r**0.2D+00 ) * r**(-1.6D+00) &
    - 10.0D+00 * sin ( 100.0D+00 * r**0.2D+00 ) * 0.8D+00 * r**(-1.8D+00)

  h(1,1) = ( arr * b + 2.0D+00 * ar * br + a * brr ) * rx1 * rx1 &
    + ( ar * b + a * br ) * rx1x1

  h(1,2) = ( arr * b + 2.0D+00 * ar * br + a * brr ) * rx1 * rx2 &
    + ( ar * b + a * br ) * rx1x2

  h(2,1) = ( arr * b + 2.0D+00 * ar * br + a * brr ) * rx2 * rx1 &
    + ( ar * b + a * br ) * rx2x1

  h(2,2) = ( arr * b + 2.0D+00 * ar * br + a * brr ) * rx2 * rx2 &
    + ( ar * b + a * br ) * rx2x2

  return
end
subroutine p28_n ( n )

!*****************************************************************************80
!
!! P28_N returns the number of variables for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p28_sol ( n, know, x )

!*****************************************************************************80
!
!! P28_SOL returns the solution for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 0.0D+00, 0.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p28_start ( n, x )

!*****************************************************************************80
!
!! P28_START returns a starting point for optimization for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ -5.0D+00, +10.0D+00 /)

  return
end
subroutine p28_title ( title )

!*****************************************************************************80
!
!! P28_TITLE returns a title for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Schaffer Function F7'

  return
end
subroutine p29_f ( n, x, f )

!*****************************************************************************80
!
!! P29_F evaluates the objective function for problem 29.
!
!  Discussion:
!
!    Note that F is a polynomial in X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) x(n)

  a = x(1) + x(2) + 1.0D+00

  b = 19.0D+00 - 14.0D+00 * x(1) + 3.0D+00 * x(1) * x(1) - 14.0D+00 * x(2) &
    + 6.0D+00 * x(1) * x(2) + 3.0D+00 * x(2) * x(2)

  c = 2.0D+00 * x(1) - 3.0D+00 * x(2)

  d = 18.0D+00 - 32.0D+00 * x(1) + 12.0D+00 * x(1) * x(1) + 48.0D+00 * x(2) &
    - 36.0D+00 * x(1) * x(2) + 27.0D+00 * x(2) * x(2)

  f = ( 1.0D+00 + a * a * b ) * ( 30.0D+00 + c * c * d )

  return
end
subroutine p29_g ( n, x, g )

!*****************************************************************************80
!
!! P29_G evaluates the gradient for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dbdx1
  real ( kind = 8 ) dbdx2
  real ( kind = 8 ) dddx1
  real ( kind = 8 ) dddx2
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  a = x(1) + x(2) + 1.0D+00

  b = 19.0D+00 - 14.0D+00 * x(1) + 3.0D+00 * x(1)**2 - 14.0D+00 * x(2) &
    + 6.0D+00 * x(1) * x(2) + 3.0D+00 * x(2)**2

  dbdx1 = - 14.0D+00 + 6.0D+00 * x(1) + 6.0D+00 * x(2)
  dbdx2 = - 14.0D+00 + 6.0D+00 * x(1) + 6.0D+00 * x(2)

  c = 2.0D+00 * x(1) - 3.0D+00 * x(2)

  d = 18.0D+00 - 32.0D+00 * x(1) + 12.0D+00 * x(1)**2 + 48.0D+00 * x(2) &
    - 36.0D+00 * x(1) * x(2) + 27.0D+00 * x(2)**2
  dddx1 = - 32.0D+00 + 24.0D+00 * x(1) - 36.0D+00 * x(2)
  dddx2 = 48.0D+00 - 36.0D+00 * x(1) + 54.0D+00 * x(2)

  g(1) = ( 1.0D+00 + a**2 * b ) * ( 4.0D+00 * c * d + c**2 * dddx1 ) &
        + ( 2.0D+00 * a * b + a**2 * dbdx1 ) * ( 30.0D+00 + c**2 * d )

  g(2) = ( 1.0D+00 + a**2 * b ) * ( -6.0D+00 * c * d + c**2 * dddx2 ) &
        + ( 2.0D+00 * a * b + a**2 * dbdx2 ) * ( 30.0D+00 + c**2 * d )

  return
end
subroutine p29_h ( n, x, h )

!*****************************************************************************80
!
!! P29_H evaluates the Hessian for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x(n)

  a = x(1) + x(2) + 1.0D+00

  b = 19.0D+00 - 14.0D+00 * x(1) + 3.0D+00 * x(1)**2 - 14.0D+00 * x(2) &
    + 6.0D+00 * x(1) * x(2) + 3.0D+00 * x(2)**2

  e = - 14.0D+00 + 6.0D+00 * x(1) + 6.0D+00 * x(2)

  c = 2.0D+00 * x(1) - 3.0D+00 * x(2)

  d = 18.0D+00 - 32.0D+00 * x(1) + 12.0D+00 * x(1)**2 + 48.0D+00 * x(2) &
    - 36.0D+00 * x(1) * x(2) + 27.0D+00 * x(2)**2

  r = - 32.0D+00 + 24.0D+00 * x(1) - 36.0D+00 * x(2)
  s = 48.0D+00 - 36.0D+00 * x(1) + 54.0D+00 * x(2)

  h(1,1) = 2.0D+00 * ( 2.0D+00 * a * b + a**2 * e ) &
    * ( 4.0D+00 * c * d + c**2 * r ) + ( 1.0D+00 + a**2 * b ) &
    * ( 8.0D+00 * d + 4.0D+00 * c * r + 4.0D+00 * c * r + 24.0D+00 * c**2 ) &
    + ( 2.0D+00 * b + 2.0D+00 * a * e + 2.0D+00 * a * e + 6.0D+00 * a**2 ) &
    * ( 30.0D+00 + c**2 * d )

  h(1,2) = ( 2.0D+00 * a * b + a**2 * e ) &
    * ( -2.0D+00 * c * d + c**2 * ( r + s ) ) &
    + ( 1.0D+00 + a**2 * b ) &
    * ( -12.0D+00 * d + 4.0D+00 * c * s -6.0D+00 * c * r - 36.0D+00 * c**2 ) &
    + ( 2.0D+00 * b + 4.0D+00 * a * e + 6.0D+00 * a**2 ) &
    * ( 30.0D+00 + c**2 * d )

  h(2,1) = ( 2.0D+00 * a * b + a**2 * e ) &
    * ( -2.0D+00 * c * d + c**2 * ( s + r ) ) &
    + ( 1.0D+00 + a**2 * b ) &
    * ( -12.0D+00 * d - 6.0D+00 * c * r + 4.0D+00 * c * s - 36.0D+00 * c**2 ) &
    + ( 2.0D+00 * b + 4.0D+00 * a * e + 6.0D+00 * a**2 ) &
    * ( 30.0D+00 + c**2 * d )

  h(2,2) = 2.0D+00 * ( 2.0D+00 * a * b + a**2 * e ) &
    * ( -6.0D+00 * c * d + c**2 * s ) &
    + ( 1.0D+00 + a**2 * b ) &
    * ( 18.0D+00 * d - 6.0D+00 * c * s - 6.0D+00 * c * s + 54.0D+00 * c**2 ) &
    + ( 2.0D+00 * b + 2.0D+00 * a * e + 2.0D+00 * a * e + 6.0D+00 * a**2 ) &
    * ( 30.0D+00 + c**2 * d )

  return
end
subroutine p29_n ( n )

!*****************************************************************************80
!
!! P29_N returns the number of variables for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p29_sol ( n, know, x )

!*****************************************************************************80
!
!! P29_SOL returns the solution for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 0.0D+00, -1.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p29_start ( n, x )

!*****************************************************************************80
!
!! P29_START returns a starting point for optimization for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ -0.5D+00, +0.25D+00 /)

  return
end
subroutine p29_title ( title )

!*****************************************************************************80
!
!! P29_TITLE returns a title for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Goldstein Price Polynomial'

  return
end
subroutine p30_f ( n, x, f )

!*****************************************************************************80
!
!! P30_F evaluates the objective function for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: a = 1.0D+00
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: d = 6.0D+00
  real ( kind = 8 ), parameter :: e = 10.0D+00
  real ( kind = 8 ) f
  real ( kind = 8 ) ff
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  b = 5.1D+00 / ( 4.0D+00 * pi**2 )
  c = 5.0D+00 / pi
  ff = 1.0D+00 / ( 8.0D+00 * pi )

  f = a * ( x(2) - b * x(1)**2 + c * x(1) - d )**2 &
    + e * ( 1.0D+00 - ff ) * cos ( x(1) ) + e

  return
end
subroutine p30_g ( n, x, g )

!*****************************************************************************80
!
!! P30_G evaluates the gradient for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: a = 1.0D+00
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: d = 6.0D+00
  real ( kind = 8 ), parameter :: e = 10.0D+00
  real ( kind = 8 ) ff
  real ( kind = 8 ) g(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  b = 5.1D+00 / ( 4.0D+00 * pi**2 )
  c = 5.0D+00 / pi
  ff = 1.0D+00 / ( 8.0D+00 * pi )

  g(1) = 2.0D+00 * a * ( x(2) - b * x(1)**2 + c * x(1) - d ) &
    * ( - 2.0D+00 * b * x(1) + c ) - e * ( 1.0D+00 - ff ) * sin ( x(1) )

  g(2) = 2.0D+00 * a * ( x(2) - b * x(1)**2 + c * x(1) - d )

  return
end
subroutine p30_h ( n, x, h )

!*****************************************************************************80
!
!! P30_H evaluates the Hessian for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: a = 1.0D+00
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: d = 6.0D+00
  real ( kind = 8 ), parameter :: e = 10.0D+00
  real ( kind = 8 ) ff
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  b = 5.1D+00 / ( 4.0D+00 * pi**2 )
  c = 5.0D+00 / pi
  ff = 1.0D+00 / ( 8.0D+00 * pi )

  h(1,1) = 2.0D+00 * a * ( - 2.0D+00 * b * x(1) + c ) &
    * ( - 2.0D+00 * b * x(1) + c ) &
    - 4.0D+00 * a * b * ( x(2) - b * x(1)**2 + c * x(1) - d ) &
    - e * ( 1.0D+00 - ff ) * cos ( x(1) )

  h(1,2) = 2.0D+00 * a * ( - 2.0D+00 * b * x(1) + c )

  h(2,1) = 2.0D+00 * a * ( - 2.0D+00 * b * x(1) + c )
  h(2,2) = 2.0D+00 * a

  return
end
subroutine p30_n ( n )

!*****************************************************************************80
!
!! P30_N returns the number of variables for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p30_sol ( n, know, x )

!*****************************************************************************80
!
!! P30_SOL returns the solution for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    x(1:n) = (/ -pi, 12.275D+00 /)
    know = know + 1
  else if ( know == 1 ) then
    x(1:n) = (/  pi,  2.275D+00 /)
    know = know + 1
  else if ( know == 2 ) then
    x(1:n) = (/ 9.42478D+00, 2.475D+00 /)
    know = know + 1
  else
    know = 0
  end if

  return
end
subroutine p30_start ( n, x )

!*****************************************************************************80
!
!! P30_START returns a starting point for optimization for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ -1.0D+00, 1.0D+00 /)

  return
end
subroutine p30_title ( title )

!*****************************************************************************80
!
!! P30_TITLE returns a title for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Branin RCOS Function'

  return
end
subroutine p31_f ( n, x, f )

!*****************************************************************************80
!
!! P31_F evaluates the objective function for problem 31.
!
!  Discussion:
!
!    The minimal function value is -10.15320.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( 4, m ) :: a = reshape ( &
    (/ 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
       1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
       8.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
       6.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
       3.0D+00, 7.0D+00, 3.0D+00, 7.0D+00 /), (/ 4, m /) )
  real ( kind = 8 ), save, dimension ( m ) :: c = &
    (/ 0.1D+00, 0.2D+00, 0.2D+00, 0.4D+00, 0.6D+00 /)
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  f = 0.0D+00
  do j = 1, m
    f = f - 1.0D+00 / ( c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 ) )
  end do

  return
end
subroutine p31_g ( n, x, g )

!*****************************************************************************80
!
!! P31_G evaluates the gradient for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( 4, m ) :: a = reshape ( &
    (/ 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
       1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
       8.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
       6.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
       3.0D+00, 7.0D+00, 3.0D+00, 7.0D+00 /), (/ 4, m /) )
  real ( kind = 8 ), save, dimension ( m ) :: c = &
    (/ 0.1D+00, 0.2D+00, 0.2D+00, 0.4D+00, 0.6D+00 /)
  real ( kind = 8 ) d
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  do k = 1, n
    do j = 1, m
      d = c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 )
      g(k) = g(k) + ( 2.0D+00 * ( x(k) - a(k,j) ) ) / d**2
    end do
  end do

  return
end
subroutine p31_h ( n, x, h )

!*****************************************************************************80
!
!! P31_H evaluates the Hessian for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( 4, m ) :: a = reshape ( &
    (/ 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
       1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
       8.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
       6.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
       3.0D+00, 7.0D+00, 3.0D+00, 7.0D+00 /), (/ 4, m /) )
  real ( kind = 8 ), save, dimension ( m ) :: c = &
    (/ 0.1D+00, 0.2D+00, 0.2D+00, 0.4D+00, 0.6D+00 /)
  real ( kind = 8 ) d
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do ii = 1, n
    do jj = 1, n
      do j = 1, m
        d = c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 )
        h(ii,jj) = h(ii,jj) &
          - 8.0D+00 * ( x(ii) - a(ii,j) ) * ( x(jj) - a(jj,j) ) / d**3
        if ( ii == jj ) then
          h(ii,jj) = h(ii,jj) + 2.0D+00 / d**2
        end if
      end do
    end do
  end do

  return
end
subroutine p31_n ( n )

!*****************************************************************************80
!
!! P31_N returns the number of variables for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 4

  return
end
subroutine p31_sol ( n, know, x )

!*****************************************************************************80
!
!! P31_SOL returns the solution for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p31_start ( n, x )

!*****************************************************************************80
!
!! P31_START returns a starting point for optimization for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:4) = (/ 1.0D+00, 3.0D+00, 5.0D+00, 6.0D+00 /)

  return
end
subroutine p31_title ( title )

!*****************************************************************************80
!
!! P31_TITLE returns a title for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Shekel SQRN5 Function'

  return
end
subroutine p32_f ( n, x, f )

!*****************************************************************************80
!
!! P32_F evaluates the objective function for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( 4, m ) :: a = reshape ( &
    (/ 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
       1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
       8.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
       6.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
       3.0D+00, 7.0D+00, 3.0D+00, 7.0D+00, &
       2.0D+00, 9.0D+00, 2.0D+00, 9.0D+00, &
       5.0D+00, 5.0D+00, 3.0D+00, 3.0D+00 /), (/ 4, m /) )
  real ( kind = 8 ), save, dimension ( m ) :: c = &
    (/ 0.1D+00, 0.2D+00, 0.2D+00, 0.4D+00, 0.6D+00, 0.6D+00, 0.3D+00 /)
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  f = 0.0D+00
  do j = 1, m
    f = f - 1.0D+00 / ( c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 ) )
  end do

  return
end
subroutine p32_g ( n, x, g )

!*****************************************************************************80
!
!! P32_G evaluates the gradient for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( 4, m ) :: a = reshape ( &
    (/ 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
       1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
       8.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
       6.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
       3.0D+00, 7.0D+00, 3.0D+00, 7.0D+00, &
       2.0D+00, 9.0D+00, 2.0D+00, 9.0D+00, &
       5.0D+00, 5.0D+00, 3.0D+00, 3.0D+00 /), (/ 4, m /) )

  real ( kind = 8 ), save, dimension ( m ) :: c = &
    (/ 0.1D+00, 0.2D+00, 0.2D+00, 0.4D+00, 0.6D+00, 0.6D+00, 0.3D+00 /)
  real ( kind = 8 ) d
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  do k = 1, n
    do j = 1, m
      d = c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 )
      g(k) = g(k) + ( 2.0D+00 * ( x(k) - a(k,j) ) ) / d**2
    end do
  end do

  return
end
subroutine p32_h ( n, x, h )

!*****************************************************************************80
!
!! P32_H evaluates the Hessian for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( 4, m ) :: a = reshape ( &
    (/ 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
       1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
       8.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
       6.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
       3.0D+00, 7.0D+00, 3.0D+00, 7.0D+00, &
       2.0D+00, 9.0D+00, 2.0D+00, 9.0D+00, &
       5.0D+00, 5.0D+00, 3.0D+00, 3.0D+00 /), (/ 4, m /) )

  real ( kind = 8 ), save, dimension ( m ) :: c = &
    (/ 0.1D+00, 0.2D+00, 0.2D+00, 0.4D+00, 0.6D+00, 0.6D+00, 0.3D+00 /)
  real ( kind = 8 ) d
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do ii = 1, n
    do jj = 1, n
      do j = 1, m
        d = c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 )
        h(ii,jj) = h(ii,jj) &
          - 8.0D+00 * ( x(ii) - a(ii,j) ) * ( x(jj) - a(jj,j) ) / d**3
        if ( ii == jj ) then
          h(ii,jj) = h(ii,jj) + 2.0D+00 / d**2
        end if
      end do
    end do
  end do

  return
end
subroutine p32_n ( n )

!*****************************************************************************80
!
!! P32_N returns the number of variables for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 4

  return
end
subroutine p32_sol ( n, know, x )

!*****************************************************************************80
!
!! P32_SOL returns the solution for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p32_start ( n, x )

!*****************************************************************************80
!
!! P32_START returns a starting point for optimization for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:4) = (/ 1.0D+00, 3.0D+00, 5.0D+00, 6.0D+00 /)

  return
end
subroutine p32_title ( title )

!*****************************************************************************80
!
!! P32_TITLE returns a title for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Shekel SQRN7 Function'

  return
end
subroutine p33_f ( n, x, f )

!*****************************************************************************80
!
!! P33_F evaluates the objective function for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( 4, m ) :: a = reshape ( &
    (/ 4.0, 4.0, 4.0, 4.0, &
       1.0, 1.0, 1.0, 1.0, &
       8.0, 8.0, 8.0, 8.0, &
       6.0, 6.0, 6.0, 6.0, &
       3.0, 7.0, 3.0, 7.0, &
       2.0, 9.0, 2.0, 9.0, &
       5.0, 5.0, 3.0, 3.0, &
       8.0, 1.0, 8.0, 1.0, &
       6.0, 2.0, 6.0, 2.0, &
       7.0, 3.6, 7.0, 3.6 /), (/ 4, m /) )

  real ( kind = 8 ), save, dimension ( m ) :: c = &
    (/ 0.1D+00, 0.2D+00, 0.2D+00, 0.4D+00, 0.6D+00, &
       0.6D+00, 0.3D+00, 0.7D+00, 0.5D+00, 0.5D+00 /)
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  f = 0.0D+00
  do j = 1, m
    f = f - 1.0D+00 / ( c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 ) )
  end do

  return
end
subroutine p33_g ( n, x, g )

!*****************************************************************************80
!
!! P33_G evaluates the gradient for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( 4, m ) :: a = reshape ( &
    (/ 4.0, 4.0, 4.0, 4.0, &
       1.0, 1.0, 1.0, 1.0, &
       8.0, 8.0, 8.0, 8.0, &
       6.0, 6.0, 6.0, 6.0, &
       3.0, 7.0, 3.0, 7.0, &
       2.0, 9.0, 2.0, 9.0, &
       5.0, 5.0, 3.0, 3.0, &
       8.0, 1.0, 8.0, 1.0, &
       6.0, 2.0, 6.0, 2.0, &
       7.0, 3.6, 7.0, 3.6 /), (/ 4, m /) )

  real ( kind = 8 ), save, dimension ( m ) :: c = &
    (/ 0.1D+00, 0.2D+00, 0.2D+00, 0.4D+00, 0.6D+00, &
       0.6D+00, 0.3D+00, 0.7D+00, 0.5D+00, 0.5D+00 /)
  real ( kind = 8 ) d
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  g(1:n) = 0.0D+00

  do k = 1, n
    do j = 1, m
      d = c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 )
      g(k) = g(k) + ( 2.0D+00 * ( x(k) - a(k,j) ) ) / d**2
    end do
  end do

  return
end
subroutine p33_h ( n, x, h )

!*****************************************************************************80
!
!! P33_H evaluates the Hessian for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( 4, m ) :: a = reshape ( &
    (/ 4.0, 4.0, 4.0, 4.0, &
       1.0, 1.0, 1.0, 1.0, &
       8.0, 8.0, 8.0, 8.0, &
       6.0, 6.0, 6.0, 6.0, &
       3.0, 7.0, 3.0, 7.0, &
       2.0, 9.0, 2.0, 9.0, &
       5.0, 5.0, 3.0, 3.0, &
       8.0, 1.0, 8.0, 1.0, &
       6.0, 2.0, 6.0, 2.0, &
       7.0, 3.6, 7.0, 3.6 /), (/ 4, m /) )

  real ( kind = 8 ), save, dimension ( m ) :: c = &
    (/ 0.1D+00, 0.2D+00, 0.2D+00, 0.4D+00, 0.6D+00, &
       0.6D+00, 0.3D+00, 0.7D+00, 0.5D+00, 0.5D+00 /)
  real ( kind = 8 ) d
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  do ii = 1, n
    do jj = 1, n
      do j = 1, m
        d = c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 )
        h(ii,jj) = h(ii,jj) &
          - 8.0D+00 * ( x(ii) - a(ii,j) ) * ( x(jj) - a(jj,j) ) / d**3
        if ( ii == jj ) then
          h(ii,jj) = h(ii,jj) + 2.0D+00 / d**2
        end if
      end do
    end do
  end do

  return
end
subroutine p33_n ( n )

!*****************************************************************************80
!
!! P33_N returns the number of variables for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 4

  return
end
subroutine p33_sol ( n, know, x )

!*****************************************************************************80
!
!! P33_SOL returns the solution for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p33_start ( n, x )

!*****************************************************************************80
!
!! P33_START returns a starting point for optimization for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:4) = (/ 1.0D+00, 3.0D+00, 5.0D+00, 6.0D+00 /)

  return
end
subroutine p33_title ( title )

!*****************************************************************************80
!
!! P33_TITLE returns a title for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Shekel SQRN10 Function'

  return
end
subroutine p34_f ( n, x, f )

!*****************************************************************************80
!
!! P34_F evaluates the objective function for problem 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) x(n)

  f = ( 4.0D+00 - 2.1D+00 * x(1)**2 + x(1)**4 / 3.0D+00 ) * x(1)**2 &
    + x(1) * x(2) + 4.0D+00 * ( x(2)**2 - 1.0D+00 ) * x(2)**2

  return
end
subroutine p34_g ( n, x, g )

!*****************************************************************************80
!
!! P34_G evaluates the gradient for problem 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  g(1) = 8.0D+00 * x(1) - 8.4D+00 * x(1)**3 + 2.0D+00 * x(1)**5 + x(2)
  g(2) = x(1) - 8.0D+00 * x(2) + 16.0D+00 * x(2)**3

  return
end
subroutine p34_h ( n, x, h )

!*****************************************************************************80
!
!! P34_H evaluates the Hessian for problem 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1,1) = 8.0D+00 - 25.2D+00 * x(1)**2 + 10.0D+00 * x(1)**4
  h(1,2) = 1.0D+00
  h(2,1) = 1.0D+00
  h(2,2) = -8.0D+00 + 48.0D+00 * x(2)**2

  return
end
subroutine p34_n ( n )

!*****************************************************************************80
!
!! P34_N returns the number of variables for problem 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p34_sol ( n, know, x )

!*****************************************************************************80
!
!! P34_SOL returns the solution for problem 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    x(1:n) = (/ -0.0898D+00,  0.7126D+00 /)
    know = know + 1
  else if ( know == 1 ) then
    x(1:n) = (/  0.0898D+00, -0.7126D+00 /)
    know = know + 1
  else
    know = 0
  end if

  return
end
subroutine p34_start ( n, x )

!*****************************************************************************80
!
!! P34_START returns a starting point for optimization for problem 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ -1.5D+00, 0.5D+00 /)

  return
end
subroutine p34_title ( title )

!*****************************************************************************80
!
!! P34_TITLE returns a title for problem 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Six-Hump Camel-Back Polynomial'

  return
end
subroutine p35_f ( n, x, f )

!*****************************************************************************80
!
!! P35_F evaluates the objective function for problem 35.
!
!  Discussion:
!
!    For -10 <= X(I) <= 10, the function has 760 local minima, 18 of which
!    are global minima, with minimum value -186.73.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!    Bruno Shubert,
!    A sequential method seeking the global maximum of a function,
!    SIAM Journal on Numerical Analysis,
!    Volume 9, pages 379-388, 1972.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) k_r8
  real ( kind = 8 ) x(n)

  f = 1.0D+00

  do i = 1, n
    factor = 0.0D+00
    do k = 1, 5
      k_r8 = real ( k, kind = 8 )
      factor = factor + k_r8 * cos ( ( k_r8 + 1.0D+00 ) * x(1) + k_r8 )
    end do
    f = f * factor
  end do

  return
end
subroutine p35_g ( n, x, g )

!*****************************************************************************80
!
!! P35_G evaluates the gradient for problem 35.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df2dx2
  real ( kind = 8 ) factor1
  real ( kind = 8 ) factor2
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y

  factor1 = 0.0D+00
  df1dx1 = 0.0D+00
  do i = 1, 5
    y = real ( i, kind = 8 )
    factor1 = factor1 + y * cos ( ( y + 1.0D+00 ) * x(1) + y )
    df1dx1 = df1dx1 - y * ( y + 1.0D+00 ) * sin ( ( y + 1.0D+00 ) * x(1) + y )
  end do

  factor2 = 0.0D+00
  df2dx2 = 0.0D+00
  do i = 1, 5
    y = real ( i, kind = 8 )
    factor2 = factor2 + y * cos ( ( y + 1.0D+00 ) * x(2) + y )
    df2dx2 = df2dx2 - y * ( y + 1.0D+00 ) * sin ( ( y + 1.0D+00 ) * x(2) + y )
  end do

  g(1) = df1dx1 * factor2
  g(2) = factor1 * df2dx2

  return
end
subroutine p35_h ( n, x, h )

!*****************************************************************************80
!
!! P35_H evaluates the Hessian for problem 35.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) df1dx1
  real ( kind = 8 ) df1dx11
  real ( kind = 8 ) df2dx2
  real ( kind = 8 ) df2dx22
  real ( kind = 8 ) factor1
  real ( kind = 8 ) factor2
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y

  factor1 = 0.0D+00
  df1dx1 = 0.0D+00
  df1dx11 = 0.0D+00
  do i = 1, 5
    y = real ( i, kind = 8 )
    factor1 = factor1 + y * cos ( ( y + 1.0D+00 ) * x(1) + y )
    df1dx1 = df1dx1 - y * ( y + 1.0D+00 ) * sin ( ( y + 1.0D+00 ) * x(1) + y )
    df1dx11 = df1dx11 &
      - y * ( y + 1.0D+00 )**2 * cos ( ( y + 1.0D+00 ) * x(1) + y )
  end do

  factor2 = 0.0D+00
  df2dx2 = 0.0D+00
  df2dx22 = 0.0D+00
  do i = 1, 5
    y = real ( i, kind = 8 )
    factor2 = factor2 + y * cos ( ( y + 1.0D+00 ) * x(2) + y )
    df2dx2 = df2dx2 - y * ( y + 1.0D+00 ) * sin ( ( y + 1.0D+00 ) * x(2) + y )
    df2dx22 = df2dx22 &
      - y * ( y + 1.0D+00 )**2 * cos ( ( y + 1.0D+00 ) * x(2) + y )
  end do

  h(1,1) = df1dx11 * factor2
  h(1,2) = df1dx1 * df2dx2
  h(2,1) = df1dx1 * df2dx2
  h(2,2) = factor1 * df2dx22

  return
end
subroutine p35_n ( n )

!*****************************************************************************80
!
!! P35_N returns the number of variables for problem 35.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p35_sol ( n, know, x )

!*****************************************************************************80
!
!! P35_SOL returns the solution for problem 35.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 0.0D+00, 0.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p35_start ( n, x )

!*****************************************************************************80
!
!! P35_START returns a starting point for optimization for problem 35.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 0.5D+00, 1.0D+00 /)

  return
end
subroutine p35_title ( title )

!*****************************************************************************80
!
!! P35_TITLE returns a title for problem 35.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Shubert Function'

  return
end
subroutine p36_f ( n, x, f )

!*****************************************************************************80
!
!! P36_F evaluates the objective function for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) f
  real ( kind = 8 ) m1
  real ( kind = 8 ) m2
  real ( kind = 8 ) r11
  real ( kind = 8 ) r12
  real ( kind = 8 ) r21
  real ( kind = 8 ) r22
  real ( kind = 8 ) r8_aint
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  call p36_p_get ( b, m1, m2, r11, r12, r21, r22, seed )

  a1 = r8_aint ( abs ( x(1) - r11 ) ) + r8_aint ( abs ( x(2) - r21 ) )
  a2 = r8_aint ( abs ( x(1) - r12 ) ) + r8_aint ( abs ( x(2) - r22 ) )

  if ( x(1) <= b ) then
    if ( a1 == 0.0D+00 ) then
      f = r8_aint ( m1 )
    else
      f = r8_aint ( m1 * sin ( a1 ) / a1 )
    end if
  else
    if ( a2 == 0.0D+00 ) then
      f = r8_aint ( m2 )
    else
      f = r8_aint ( m2 * sin ( a2 ) / a2 )
    end if
  end if

  return
end
subroutine p36_g ( n, x, g )

!*****************************************************************************80
!
!! P36_G evaluates the gradient for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  g(1:n) = (/ 0.0D+00, 0.0D+00 /)

  return
end
subroutine p36_h ( n, x, h )

!*****************************************************************************80
!
!! P36_H evaluates the Hessian for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  return
end
subroutine p36_n ( n )

!*****************************************************************************80
!
!! P36_N returns the number of variables for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p36_p_get ( b, m1, m2, r11, r12, r21, r22, seed )

!*****************************************************************************80
!
!! P36_P_GET gets the values of the parameters for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) B, M1, M2, R11, R12, R21, R22, 
!    problem parameters.
!
!    Output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) m1
  real ( kind = 8 ) m2
  real ( kind = 8 ) r11
  real ( kind = 8 ) r12
  real ( kind = 8 ) r21
  real ( kind = 8 ) r22
  integer ( kind = 4 ) seed

  call p36_p_val ( 'GET', b, m1, m2, r11, r12, r21, r22, seed )

  return
end
subroutine p36_p_init ( seed )

!*****************************************************************************80
!
!! P36_P_INIT initializes parameters for problem 36.
!
!  Discussion:
!
!    This routine can be called to choose values for the parameters at random.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2007
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
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) m1
  real ( kind = 8 ) m2
  real ( kind = 8 ) r11
  real ( kind = 8 ) r12
  real ( kind = 8 ) r21
  real ( kind = 8 ) r22
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  m1  = r8_uniform ( 0.0D+00, 100.0D+00, seed, m1 )
  m2  = r8_uniform ( 0.0D+00, 100.0D+00, seed, m2 )
  b   = r8_uniform ( 0.0D+00, 10.0D+00,  seed, b )
  r11 = r8_uniform ( 0.0D+00, b,         seed, r11 )
  r12 = r8_uniform ( b,       10.0D+00,  seed, r12 )
  r21 = r8_uniform ( 0.0D+00, 10.0D+00,  seed, r21 )
  r22 = r8_uniform ( 0.0D+00, 10.0D+00,  seed, r22 )

  call p36_p_set ( b, m1, m2, r11, r12, r21, r22, seed )

  return
end
subroutine p36_p_set ( b, m1, m2, r11, r12, r21, r22, seed )

!*****************************************************************************80
!
!! P36_P_SET sets parameters for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) B, M1, M2, R11, R12, R21, R22, problem parameters.
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) m1
  real ( kind = 8 ) m2
  real ( kind = 8 ) r11
  real ( kind = 8 ) r12
  real ( kind = 8 ) r21
  real ( kind = 8 ) r22
  integer ( kind = 4 ) seed

  call p36_p_val ( 'SET', b, m1, m2, r11, r12, r21, r22, seed )

  return
end
subroutine p36_p_val ( action, b, m1, m2, r11, r12, r21, r22, seed )

!*****************************************************************************80
!
!! P36_P_VAL sets or gets parameters for problem 36.
!
!  Discussion:
!
!    If ACTION is 'SET', the parameters are input values, and set 
!    the internal values.
!
!    If ACTION is 'GET', the parametrs are output values, copied from
!    the internal values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, specifies the action desired.
!    'SET' sets the parameters;
!    'GET' gets the current values.
!
!    Input/output, real ( kind = 8 ) B, M1, M2, R11, R12, R21, R22, 
!    problem parameters.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ) b
  real ( kind = 8 ), save :: b_save = 0.0D+00
  real ( kind = 8 ) m1
  real ( kind = 8 ), save :: m1_save = 0.0D+00
  real ( kind = 8 ) m2
  real ( kind = 8 ), save :: m2_save = 0.0D+00
  real ( kind = 8 ) r11
  real ( kind = 8 ), save :: r11_save = 0.0D+00
  real ( kind = 8 ) r12
  real ( kind = 8 ), save :: r12_save = 0.0D+00
  real ( kind = 8 ) r21
  real ( kind = 8 ), save :: r21_save = 0.0D+00
  real ( kind = 8 ) r22
  real ( kind = 8 ), save :: r22_save = 0.0D+00
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_save

  if ( action == 'GET' .or. action == 'get' ) then

    b = b_save
    m1 = m1_save
    m2 = m2_save
    r11 = r11_save
    r12 = r12_save
    r21 = r21_save
    r22 = r22_save
    seed = seed_save

  else if ( action == 'SET' .or. action == 'set' ) then

    b_save = b
    m1_save = m1
    m2_save = m2
    r11_save = r11
    r12_save = r12
    r21_save = r21
    r22_save = r22
    seed_save = seed

  end if

  return
end
subroutine p36_sol ( n, know, x )

!*****************************************************************************80
!
!! P36_SOL returns the solution for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  integer ( kind = 4 ) know
  real ( kind = 8 ) m1
  real ( kind = 8 ) m2
  real ( kind = 8 ) r11
  real ( kind = 8 ) r12
  real ( kind = 8 ) r21
  real ( kind = 8 ) r22
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
  
    know = 1

    call p36_p_get ( b, m1, m2, r11, r12, r21, r22, seed )

    if ( m2 < m1 ) then
      x(1:2) = (/ r11, r21 /)
    else
      x(1:2) = (/ r12, r22 /)
    end if

  else

    know = 0

  end if

  return
end
subroutine p36_start ( n, x )

!*****************************************************************************80
!
!! P36_START returns a starting point for optimization for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 0.5D+00, 1.0D+00 /)

  return
end
subroutine p36_title ( title )

!*****************************************************************************80
!
!! P36_TITLE returns a title for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Stuckman Function'

  return
end
subroutine p37_f ( n, x, f )

!*****************************************************************************80
!
!! P37_F evaluates the objective function for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  arg = - ( x(1) - pi )**2 - ( x(2) - pi )**2
  f = - cos ( x(1) ) * cos ( x(2) ) * exp ( arg )

  return
end
subroutine p37_g ( n, x, g )

!*****************************************************************************80
!
!! P37_G evaluates the gradient for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) g(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  arg = - ( x(1) - pi )**2 - ( x(2) - pi )**2
 
  g(1) = ( sin ( x(1) ) * cos ( x(2) ) &
       + 2.0D+00 * cos ( x(1) ) * cos ( x(2) ) * ( x(1) - pi ) ) &
       * exp ( arg )

  g(2) = ( cos ( x(1) ) * sin ( x(2) ) &
       + 2.0D+00 * cos ( x(1) ) * cos ( x(2) ) * ( x(2) - pi ) ) &
       * exp ( arg )

  return
end
subroutine p37_h ( n, x, h )

!*****************************************************************************80
!
!! P37_H evaluates the Hessian for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) dargdx1
  real ( kind = 8 ) dargdx2
  real ( kind = 8 ) dfdx1
  real ( kind = 8 ) dfdx2
  real ( kind = 8 ) factor
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  arg = - ( x(1) - pi )**2 - ( x(2) - pi )**2
  dargdx1 = - 2.0D+00 * ( x(1) - pi )
  dargdx2 = - 2.0D+00 * ( x(2) - pi )

  factor = cos ( x(2) ) * ( sin ( x(1) ) - cos ( x(1) ) * dargdx1 )
  dfdx1 = cos ( x(2) ) * ( cos ( x(1) ) + sin ( x(1) ) * dargdx1 &
    + 2.0D+00 * cos ( x(1) ) ) 
  dfdx2 = - sin ( x(2) ) * ( sin ( x(1) ) - cos ( x(1) ) * dargdx1 )

  h(1,1) = ( dfdx1 + factor * dargdx1 ) * exp ( arg )
  h(1,2) = ( dfdx2 + factor * dargdx2 ) * exp ( arg )

  factor = cos ( x(1) ) * ( sin ( x(2) ) - cos ( x(2) ) * dargdx2 )
  dfdx1 = - sin ( x(1) ) * ( sin ( x(2) ) - cos ( x(2) ) * dargdx2 )
  dfdx2 = cos ( x(1) ) * ( cos ( x(2) ) + sin ( x(2) ) * dargdx2 &
    + 2.0D+00 * cos ( x(2) ) )

  h(2,1) = ( dfdx1 + factor * dargdx1 ) * exp ( arg )
  h(2,2) = ( dfdx2 + factor * dargdx2 ) * exp ( arg )

  return
end
subroutine p37_n ( n )

!*****************************************************************************80
!
!! P37_N returns the number of variables for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p37_sol ( n, know, x )

!*****************************************************************************80
!
!! P37_SOL returns the solution for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ pi, pi /)
  else
    know = 0
  end if

  return
end
subroutine p37_start ( n, x )

!*****************************************************************************80
!
!! P37_START returns a starting point for optimization for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 0.5D+00, 1.0D+00 /)

  return
end
subroutine p37_title ( title )

!*****************************************************************************80
!
!! P37_TITLE returns a title for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Easom Function'

  return
end
subroutine p38_f ( n, x, f )

!*****************************************************************************80
!
!! P38_F evaluates the objective function for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  f =           x(1) * x(1) - 0.3D+00 * cos ( 3.0D+00 * pi * x(1) ) &
    + 2.0D+00 * x(2) * x(2) - 0.4D+00 * cos ( 4.0D+00 * pi * x(2) ) &
    + 0.7D+00

  return
end
subroutine p38_g ( n, x, g )

!*****************************************************************************80
!
!! P38_G evaluates the gradient for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  g(1) = 2.0D+00 * x(1) + 0.9D+00 * pi * sin ( 3.0D+00 * pi * x(1) )
  g(2) = 4.0D+00 * x(2) + 1.6D+00 * pi * sin ( 4.0D+00 * pi * x(2) )

  return
end
subroutine p38_h ( n, x, h )

!*****************************************************************************80
!
!! P38_H evaluates the Hessian for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  h(1,1) = 2.0D+00 + 2.7D+00 * pi**2 * cos ( 3.0D+00 * pi * x(1) )
  h(1,2) = 0.0D+00

  h(2,1) = 0.0D+00
  h(2,2) = 4.0D+00 + 6.4D+00 * pi**2 * cos ( 4.0D+00 * pi * x(2) )

  return
end
subroutine p38_n ( n )

!*****************************************************************************80
!
!! P38_N returns the number of variables for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p38_sol ( n, know, x )

!*****************************************************************************80
!
!! P38_SOL returns the solution for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 0.0D+00, 0.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p38_start ( n, x )

!*****************************************************************************80
!
!! P38_START returns a starting point for optimization for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 0.5D+00, 1.0D+00 /)

  return
end
subroutine p38_title ( title )

!*****************************************************************************80
!
!! P38_TITLE returns a title for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Bohachevsky Function #1'

  return
end
subroutine p39_f ( n, x, f )

!*****************************************************************************80
!
!! P39_F evaluates the objective function for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  f = x(1) * x(1) + 2.0D+00 * x(2) * x(2) &
    - 0.3D+00 * cos ( 3.0D+00 * pi * x(1) ) &
    * cos ( 4.0D+00 * pi * x(2) ) + 0.3D+00

  return
end
subroutine p39_g ( n, x, g )

!*****************************************************************************80
!
!! P39_G evaluates the gradient for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  g(1) = 2.0D+00 * x(1) &
    + 0.9D+00 * pi * sin ( 3.0D+00 * pi * x(1) ) &
    * cos ( 4.0D+00 * pi * x(2) )

  g(2) = 4.0D+00 * x(2) &
    + 1.2D+00 * pi * cos ( 3.0D+00 * pi * x(1) ) &
    * sin ( 4.0D+00 * pi * x(2) )

  return
end
subroutine p39_h ( n, x, h )

!*****************************************************************************80
!
!! P39_H evaluates the Hessian for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  h(1,1) = 2.0D+00 + 2.7D+00 * pi**2 * cos ( 3.0D+00 * pi * x(1) ) &
    * cos ( 4.0D+00 * pi * x(2) )

  h(1,2) = - 3.6D+00 * pi**2 * sin ( 3.0D+00 * pi * x(1) ) &
    * sin ( 4.0D+00 * pi * x(2) )

  h(2,1) = - 3.6D+00 * pi**2 * sin ( 3.0D+00 * pi * x(1) ) &
    * sin ( 4.0D+00 * pi * x(2) )

  h(2,2) = 4.0D+00 + 4.8D+00 * pi**2 * cos ( 3.0D+00 * pi * x(1) ) & 
    * cos ( 4.0D+00 * pi * x(2) )

  return
end
subroutine p39_n ( n )

!*****************************************************************************80
!
!! P39_N returns the number of variables for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p39_sol ( n, know, x )

!*****************************************************************************80
!
!! P39_SOL returns the solution for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 0.0D+00, 0.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p39_start ( n, x )

!*****************************************************************************80
!
!! P39_START returns a starting point for optimization for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 0.6D+00, 1.3D+00 /)

  return
end
subroutine p39_title ( title )

!*****************************************************************************80
!
!! P39_TITLE returns a title for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Bohachevsky Function #2'

  return
end
subroutine p40_f ( n, x, f )

!*****************************************************************************80
!
!! P40_F evaluates the objective function for problem 40.
!
!  Discussion:
!
!    There is a typo in the reference.  I'm just guessing at the correction.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  f = x(1)**2 + 2.0D+00 * x(2)**2 &
    - 0.3D+00 * cos ( 3.0D+00 * pi * x(1) ) &
    + cos ( 4.0D+00 * pi * x(2) ) + 0.3D+00

  return
end
subroutine p40_g ( n, x, g )

!*****************************************************************************80
!
!! P40_G evaluates the gradient for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  g(1) = 2.0D+00 * x(1) + 0.9D+00 * pi * sin ( 3.0D+00 * pi * x(1) )
  g(2) = 4.0D+00 * x(2) - 4.0D+00 * pi * sin ( 4.0D+00 * pi * x(2) )

  return
end
subroutine p40_h ( n, x, h )

!*****************************************************************************80
!
!! P40_H evaluates the Hessian for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  h(1,1) = 2.0D+00 + 2.7D+00 * pi**2 * cos ( 3.0D+00 * pi * x(1) )
  h(1,2) = 0.0D+00
  h(2,1) = 0.0D+00
  h(2,2) = 4.0D+00 - 16.0D+00 * pi**2 * cos ( 4.0D+00 * pi * x(2) )

  return
end
subroutine p40_n ( n )

!*****************************************************************************80
!
!! P40_N returns the number of variables for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p40_sol ( n, know, x )

!*****************************************************************************80
!
!! P40_SOL returns the solution for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 0.0D+00, 0.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p40_start ( n, x )

!*****************************************************************************80
!
!! P40_START returns a starting point for optimization for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 0.5D+00, 1.0D+00 /)

  return
end
subroutine p40_title ( title )

!*****************************************************************************80
!
!! P40_TITLE returns a title for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Bohachevsky Function #3'

  return
end
subroutine p41_f ( n, x, f )

!*****************************************************************************80
!
!! P41_F evaluates the objective function for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) x(n)

  f = 100.0D+00 * ( x(2) - x(1)**2 )**2 &
    + ( 1.0D+00 - x(1) )**2 &
    + 90.0D+00 * ( x(4) - x(3)**2 )**2 &
    + ( 1.0D+00 - x(3) )**2 &
    + 10.1D+00 * ( ( x(2) - 1.0D+00 )**2 + ( x(4) - 1.0D+00 )**2 ) &
    + 19.8D+00 * ( x(2) - 1.0D+00 ) * ( x(4) - 1.0D+00 )

  return
end
subroutine p41_g ( n, x, g )

!*****************************************************************************80
!
!! P41_G evaluates the gradient for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  g(1) = 400.0D+00 * x(1)**3 - 400.0D+00 * x(2) * x(1) &
    + 2.0D+00 * x(1) - 2.0D+00

  g(2) = -200.0D+00 * x(1)**2 + 220.2D+00 * x(2) + 19.8D+00 * x(4) - 40.0D+00

  g(3) = -360.0D+00 * x(3) * x(4) + 360.0D+00 * x(3)**3 &
    + 2.0D+00 * x(3) - 2.0D+00

  g(4) = + 180.0D+00 * x(4) - 180.0D+00 * x(3)**2 + 20.2D+00 * x(4) &
    + 19.8D+00 * x(2) - 40.0D+00

  return
end
subroutine p41_h ( n, x, h )

!*****************************************************************************80
!
!! P41_H evaluates the Hessian for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  h(1,1) = 1200.0D+00 * x(1)**2 - 400.0D+00 * x(2) + 2.0D+00
  h(1,2) = - 400.0D+00 * x(1)

  h(2,1) = -400.0D+00 * x(1)
  h(2,2) = 220.2D+00
  h(2,4) = 19.8D+00

  h(3,3) = -360.0D+00 * x(4) + 1080.0D+00 * x(3)**2 + 2.0D+00
  h(3,4) = - 360.0D+00 * x(3)

  h(4,2) = 19.8D+00
  h(4,3) = - 360.0D+00 * x(3)
  h(4,4) = 200.2D+00

  return
end
subroutine p41_n ( n )

!*****************************************************************************80
!
!! P41_N returns the number of variables for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 4

  return
end
subroutine p41_sol ( n, know, x )

!*****************************************************************************80
!
!! P41_SOL returns the solution for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p41_start ( n, x )

!*****************************************************************************80
!
!! P41_START returns a starting point for optimization for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 0.5D+00, 1.0D+00,-0.5D+00,-1.0D+00 /)

  return
end
subroutine p41_title ( title )

!*****************************************************************************80
!
!! P41_TITLE returns a title for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Colville Polynomial'

  return
end
subroutine p42_f ( n, x, f )

!*****************************************************************************80
!
!! P42_F evaluates the objective function for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    MJD Powell,
!    An Efficient Method for Finding the Minimum of a Function of
!    Several Variables Without Calculating Derivatives,
!    Computer Journal, 
!    Volume 7, Number 2, pages 155-162, 1964.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) term
  real ( kind = 8 ) x(n)

  if ( x(2) == 0.0D+00 ) then
    term = 0.0D+00
  else
    arg = ( x(1) + 2.0D+00 * x(2) + x(3) ) / x(2)
    term = exp ( - arg**2 )
  end if

  f = 3.0D+00 &
    - 1.0D+00 / ( 1.0D+00 + ( x(1) - x(2) )**2 ) &
    - sin ( 0.5D+00 * pi * x(2) * x(3) ) &
    - term

  return
end
subroutine p42_g ( n, x, g )

!*****************************************************************************80
!
!! P42_G evaluates the gradient for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) g(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) term
  real ( kind = 8 ) x(n)

  g(1) = 2.0D+00 * ( x(1) - x(2) ) / ( 1.0D+00 + ( x(1) - x(2) )**2 )**2

  g(2) = 2.0D+00 * ( x(2) - x(1) ) / ( 1.0D+00 + ( x(1) - x(2) )**2 )**2 &
    - 0.5D+00 * pi * x(3) * cos ( 0.5D+00 * pi * x(2) * x(3) )

  g(3) = - 0.5D+00 * pi * x(2) * cos ( 0.5D+00 * pi * x(2) * x(3) )
 
  if ( x(2) /= 0.0D+00 ) then

    arg = ( x(1) + 2.0D+00 * x(2) + x(3) ) / x(2)
    term = exp ( - arg**2 )

    g(1) = g(1) + 2.0D+00 * term * ( x(1) + 2.0D+00 * x(2) + x(3) ) / x(2)**2
    g(2) = g(2) - 2.0D+00 * term * ( x(1) + 2.0D+00 * x(2) + x(3) ) &
      * ( x(1) + x(3) ) / x(2)**3
    g(3) = g(3) + 2.0D+00 * term * ( x(1) + 2.0D+00 * x(2) + x(3) ) / x(2)**2

  end if

  return
end
subroutine p42_h ( n, x, h )

!*****************************************************************************80
!
!! P42_H evaluates the Hessian for problem 42.
!
!  Discussion:
!
!    I haven't written this yet.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1:n,1:n) = 0.0D+00

  return
end
subroutine p42_n ( n )

!*****************************************************************************80
!
!! P42_N returns the number of variables for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p42_sol ( n, know, x )

!*****************************************************************************80
!
!! P42_SOL returns the solution for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:n) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p42_start ( n, x )

!*****************************************************************************80
!
!! P42_START returns a starting point for optimization for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 0.0D+00, 1.0D+00, 2.0D+00 /)

  return
end
subroutine p42_title ( title )

!*****************************************************************************80
!
!! P42_TITLE returns a title for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Powell 3D Function'

  return
end
subroutine p43_f ( n, x, f )

!*****************************************************************************80
!
!! P43_F evaluates the objective function for problem 43.
!
!  Discussion:
!
!    This function has 4 global minima:
!
!      X* = (  3,        2       ), F(X*) = 0.
!      X* = (  3.58439, -1.84813 ), F(X*) = 0.
!      X* = ( -3.77934, -3.28317 ), F(X*) = 0.
!      X* = ( -2.80512,  3.13134 ), F(X*) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Himmelblau,
!    Applied Nonlinear Programming,
!    McGraw Hill, 1972,
!    ISBN13: 978-0070289215,
!   LC: T57.8.H55.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) x(n)

 f = ( x(1)**2 + x(2) - 11.0D+00 )**2 &
    + ( x(1) + x(2)**2 - 7.0D+00 )**2

  return
end
subroutine p43_g ( n, x, g )

!*****************************************************************************80
!
!! P43_G evaluates the gradient for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) G(N), the gradient of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(n)
  real ( kind = 8 ) x(n)

  g(1) = 4.0D+00 * ( x(1)**2 + x(2) - 11.0D+00 ) * x(1) &
       + 2.0D+00 * ( x(1) + x(2)**2 - 7.0D+00 )

  g(2) = 2.0D+00 * ( x(1)**2 + x(2) - 11.0D+00 ) &
       + 4.0D+00 * ( x(1) + x(2)**2 - 7.0D+00 ) * x(2)

  return
end
subroutine p43_h ( n, x, h )

!*****************************************************************************80
!
!! P43_H evaluates the Hessian for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1,1) = 8.0D+00 * x(1) * x(1) &
         + 4.0D+00 * ( x(1) * x(1) + x(2) - 11.0D+00 ) &
         + 2.0D+00

  h(1,2) = 4.0D+00 * x(1) + 4.0D+00 * x(2)

  h(2,1) = 4.0D+00 * x(1) + 4.0D+00 * x(2)

  h(2,2) = 2.0D+00 &
         + 8.0D+00 * x(2) * x(2) &
         + 4.0D+00 * ( x(1) + x(2) * x(2) - 7.0D+00 )

  return
end
subroutine p43_n ( n )

!*****************************************************************************80
!
!! P43_N returns the number of variables for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.  If N is positive, N represents the
!    only legal value for N for this problem.  If N is
!    negative, the absolute value of N is the least legal
!    value of N, but other values are allowable.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p43_sol ( n, know, x )

!*****************************************************************************80
!
!! P43_SOL returns the solution for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the problem.  This value
!    is only needed for those problems with variable N.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(N), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(n)

  if ( know == 0 ) then
    know = 1
    x(1:2) = (/ 3.0D+00, 2.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p43_start ( n, x )

!*****************************************************************************80
!
!! P43_START returns a starting point for optimization for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables X.
!
!    Output, real ( kind = 8 ) X(N), a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ -1.3D+00, 2.7D+00 /)

  return
end
subroutine p43_title ( title )

!*****************************************************************************80
!
!! P43_TITLE returns a title for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The Himmelblau function.'

  return
end
subroutine normal_01_sample ( x )

!*****************************************************************************80
!
!! NORMAL_01_SAMPLE samples the standard Normal PDF.
!
!  Discussion:
!
!    The standard normal distribution has mean 0 and standard
!    deviation 1.
!
!  Method:
!
!    The Box-Muller method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ), save :: iset = -1
  real ( kind = 8 ), parameter :: PI = &
    3.14159265358979323846264338327950288419716939937510D+00
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xsave = 0.0D+00

  if ( iset == -1 ) then
    call random_seed ( )
    iset = 0
  end if

  if ( iset == 0 ) then

    call random_number ( harvest = v1 )

    if ( v1 <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  V1 <= 0.'
      write ( *, '(a,g14.6)' ) '  V1 = ', v1
      stop
    end if

    call random_number ( harvest = v2 )

    if ( v2 <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  V2 <= 0.'
      write ( *, '(a,g14.6)' ) '  V2 = ', v2
      stop
    end if

    x = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * PI * v2 )

    xsave = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * PI * v2 )

    iset = 1

  else

    x = xsave
    iset = 0

  end if

  return
end
function r8_aint ( x )

!****************************************************************************80
!
!! R8_AINT truncates an R8 argument to an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2011
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) VALUE, the truncated version of X.
!
  implicit none

  real ( kind = 8 ) r8_aint
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    value = - int ( abs ( x ) )
  else
    value =   int ( abs ( x ) )
  end if

  r8_aint = value

  return
end
function r8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

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
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer ( kind = 4 ) arithmetic never requires more than 32 bits,
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
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8vec_even ( alo, ahi, n, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns N evenly spaced double precision values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) ahi
  real ( kind = 8 ) alo
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
