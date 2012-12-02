subroutine box_nd ( func, dim_num, order, xtab, weight, result, eval_num )

!*****************************************************************************80
!
!! BOX_ND estimates a multidimensional integral using a product rule.
!
!  Discussion:
!
!    The routine creates a DIM_NUM-dimensional product rule from a 1D rule
!    supplied by the user.  The routine is fairly inflexible.  If
!    you supply a rule for integration from -1 to 1, then your product
!    box must be a product of DIM_NUM copies of the interval [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, a routine which evaluates
!    the function to be integrated, of the form:
!      function func ( dim_num, x )
!      integer ( kind = 4 ) dim_num
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(dim_num)
!      func = ...
!      return
!      end
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ORDER, the number of points used in the 1D rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), the abscissas of the 1D rule.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the 1D rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
!    Output, integer ( kind = 4 ) EVAL_NUM, the number of function evaluations.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order

  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) indx(dim_num)
  integer ( kind = 4 ) k
  real ( kind = 8 ) result
  real ( kind = 8 ) w
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) xtab(order)

  eval_num = 0

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BOX_ND - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BOX_ND - Fatal error!'
    write ( *, '(a)' ) '  ORDER < 1'
    write ( *, '(a,i8)' ) '  ORDER = ', order
    stop
  end if

  k = 0
  result = 0.0D+00

  do

    call tuple_next ( 1, order, dim_num, k, indx )

    if ( k == 0  ) then
      exit
    end if

    w = product ( weight(indx(1:dim_num)) )

    x(1:dim_num) = xtab(indx(1:dim_num))

    result = result + w * func ( dim_num, x )
    eval_num = eval_num + 1

  end do

  return
end
subroutine monte_carlo_nd ( func, dim_num, a, b, eval_num, seed, result )

!*****************************************************************************80
!
!! MONTE_CARLO_ND estimates a multidimensional integral using Monte Carlo.
!
!  Discussion:
!
!    Unlike the other routines, this routine requires the user to specify
!    the number of function evaluations as an INPUT quantity.
!
!    No attempt at error estimation is made.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, a routine which evaluates
!    the function to be integrated, of the form:
!      function func ( dim_num, x )
!      integer ( kind = 4 ) dim_num
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(dim_num)
!      func = ...
!      return
!      end
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the integration limits.
!
!    Input, integer ( kind = 4 ) EVAL_NUM, the number of function evaluations.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  real ( kind = 8 ) volume
  real ( kind = 8 ) x(dim_num)

  result = 0.0D+00

  do i = 1, eval_num

    call r8vec_uniform_01 ( dim_num, seed, x )

    result = result + func ( dim_num, x )

  end do

  volume = product ( b(1:dim_num) - a(1:dim_num) )

  result = result * volume / real ( eval_num, kind = 8 )

  return
end
subroutine p5_nd ( func, dim_num, a, b, result, eval_num )

!*****************************************************************************80
!
!! P5_ND estimates a multidimensional integral with a formula of exactness 5.
!
!  Discussion:
!
!    The routine uses a method which is exact for polynomials of total
!    degree 5 or less.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, a routine which evaluates
!    the function to be integrated, of the form:
!      function func ( dim_num, x )
!      integer ( kind = 4 ) dim_num
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(dim_num)
!      func = ...
!      return
!      end
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the integration limits.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
!    Output, integer ( kind = 4 ) EVAL_NUM, the number of function evaluations.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) a0
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) a3
  real ( kind = 8 ) a4
  real ( kind = 8 ) a5
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ) en
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) result
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) sum3
  real ( kind = 8 ) volume
  real ( kind = 8 ) work(dim_num)

  eval_num = 0

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P5_ND - Fatal error!'
    write ( *, '(a,i8)' ) '  DIM_NUM < 1, DIM_NUM = ', dim_num
    stop
  end if

  a2 = 25.0D+00 / 324.0D+00
  a3 = sqrt ( 0.6D+00 )
  en = real ( dim_num, kind = 8 )
  a0 = ( 25.0D+00 * en * en - 115.0D+00 * en + 162.0D+00 ) / 162.0D+00
  a1 = ( 70.0D+00 - 25.0D+00 * en ) / 162.0D+00

  volume = product ( b(1:dim_num) - a(1:dim_num) )
  work(1:dim_num) = 0.5D+00 * ( a(1:dim_num) + b(1:dim_num) )

  result = 0.0D+00
  if ( volume == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P5_ND - Warning!'
    write ( *, '(a)' ) '  Volume = 0, integral = 0.'
    return
  end if

  sum1 = a0 * func ( dim_num, work )
  eval_num = eval_num + 1

  sum2 = 0.0D+00
  sum3 = 0.0D+00

  do i = 1, dim_num

    work(i) = 0.5D+00 * ( ( a(i) + b(i) ) + a3 * ( b(i) - a(i) ) )
    sum2 = sum2 + func ( dim_num, work )
    eval_num = eval_num + 1

    work(i) = 0.5D+00 * ( ( a(i) + b(i) ) - a3 * ( b(i) - a(i) ) )
    sum2 = sum2 + func ( dim_num, work )
    eval_num = eval_num + 1

    work(i) = 0.5D+00 * ( a(i) + b(i) )

  end do

  if ( 1 < dim_num ) then

    a4 = a3

    do

      do i = 1, dim_num - 1

        work(i) = 0.5D+00 * ( ( a(i) + b(i) ) + a4 * ( b(i) - a(i) ) )
        a5 = a3

        do

          do j = i + 1, dim_num
            work(j) = 0.5D+00 * ( ( a(j) + b(j) ) + a5 * ( b(j) - a(j) ) )
            sum3 = sum3 + func ( dim_num, work )
            eval_num = eval_num + 1
            work(j) = 0.5D+00 * ( a(j) + b(j) )
          end do

          a5 = -a5

          if ( 0.0D+00 <= a5 ) then
            exit
          end if

        end do

        work(i) = 0.5D+00 * ( a(i) + b(i) )

      end do

      a4 = -a4

      if ( 0.0D+00 <= a4 ) then
        exit
      end if

    end do

  end if

  result = volume * ( sum1 + a1 * sum2 + a2 * sum3 )

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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
      seed = seed + huge ( seed )
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine romberg_nd ( func, a, b, dim_num, sub_num, it_max, tol, result, &
  ind, eval_num )

!*****************************************************************************80
!
!! ROMBERG_ND estimates a multidimensional integral using Romberg integration.
!
!  Discussion:
!
!    The routine uses a Romberg method based on the midpoint rule.
!
!    In the reference, this routine is called "NDIMRI".
!
!    Thanks to Barak Bringoltz for pointing out problems in a previous
!    FORTRAN90 implementation of this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2006
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, a routine which evaluates
!    the function to be integrated, of the form:
!      function func ( dim_num, x )
!      integer ( kind = 4 ) dim_num
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(dim_num)
!      func = ...
!      return
!      end
!
!    Input, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the integration limits.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SUB_NUM(DIM_NUM), the number of subintervals into
!    which the I-th integration interval (A(I), B(I)) is
!    initially subdivided.  SUB_NUM(I) must be greater than 0.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations to
!    be performed.  The number of function evaluations on
!    iteration J is at least J**DIM_NUM, which grows very rapidly.
!    IT_MAX should be small!
!
!    Input, real ( kind = 8 ) TOL, an error tolerance for the approximation
!    of the integral.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
!    Output, integer ( kind = 4 ) IND, error return flag.
!    IND = -1 if the error tolerance could not be achieved.
!    IND = 1 if the error tolerance was achieved.
!
!    Output, integer ( kind = 4 ) EVAL_NUM, the number of function evaluations.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) IWORK(DIM_NUM), a pointer used to generate all the
!    points X in the product region.
!
!    Local, integer ( kind = 4 ) IWORK2(IT_MAX), a counter of the number of points
!    used at each step of the Romberg iteration.
!
!    Local, integer ( kind = 4 ) SUB_NUM2(DIM_NUM), the number of subintervals used
!    in each direction, a refinement of the user's input SUB_NUM.
!
!    Local, real ( kind = 8 ) TABLE(IT_MAX), the difference table.
!
!    Local, real ( kind = 8 ) X(DIM_NUM), an evaluation point.
!
  implicit none

  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ) en
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ) factor
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iwork(dim_num)
  integer ( kind = 4 ) iwork2(it_max)
  integer ( kind = 4 ) kdim
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) sub_num(dim_num)
  integer ( kind = 4 ) sub_num2(dim_num)
  real ( kind = 8 ) result
  real ( kind = 8 ) result_old
  real ( kind = 8 ) rnderr
  real ( kind = 8 ) submid
  real ( kind = 8 ) sum1
  real ( kind = 8 ) weight
  real ( kind = 8 ) table(it_max)
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(dim_num)

  eval_num = 0

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROMBERG_ND - Fatal error!'
    write ( *, '(a,i8)' ) '  DIM_NUM is less than 1.  DIM_NUM = ', dim_num
    stop
  end if

  if ( it_max < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROMBERG_ND - Fatal error!'
    write ( *, '(a,i8)' ) '  IT_MAX is less than 1.  IT_MAX = ', it_max
    stop
  end if

  do i = 1, dim_num
    if ( sub_num(i) <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROMBERG_ND - Fatal error!'
      write ( *, '(a)' ) '  SUB_NUM(I) is less than 1.'
      write ( *, '(a,i8)' ) '  for I = ', i
      write ( *, '(a,i8)' ) '  SUB_NUM(I) = ', sub_num(i)
      stop
    end if
  end do

  ind = 0
  rnderr = epsilon ( 1.0D+00 )
  iwork2(1) = 1
  sub_num2(1:dim_num) = sub_num(1:dim_num)

  if ( 1 < it_max ) then
    iwork2(2) = 2
  end if

  it = 1

  do

    sum1 = 0.0D+00

    weight = product ( ( b(1:dim_num) - a(1:dim_num) ) &
      / real ( sub_num2(1:dim_num), kind = 8 ) )
!
!  Generate every point X in the product region, and evaluate F(X).
!
    iwork(1:dim_num) = 1

    do

      x(1:dim_num) = &
       ( real ( 2 * sub_num2(1:dim_num) - 2 * iwork(1:dim_num) + 1, kind = 8 ) &
       * a(1:dim_num)   &
       + real (                         + 2 * iwork(1:dim_num) - 1, kind = 8 ) &
       * b(1:dim_num) ) &
       / real ( 2 * sub_num2(1:dim_num),                            kind = 8 )

      sum1 = sum1 + func ( dim_num, x )
      eval_num = eval_num + 1

      kdim = dim_num

      do while ( 0 < kdim )

        if ( iwork(kdim) < sub_num2(kdim) ) then
          iwork(kdim) = iwork(kdim) + 1
          exit
        end if

        iwork(kdim) = 1

        kdim = kdim - 1

      end do

      if ( kdim == 0 ) then
        exit
      end if

    end do
!
!  Done with summing.
!
    table(it) = weight * sum1

    if ( it <= 1 ) then

      result = table(1)
      result_old = result

      if ( it_max <= it ) then
        ind = 1
        exit
      end if

      it = it + 1

      sub_num2(1:dim_num) = iwork2(it) * sub_num2(1:dim_num)

      cycle

    end if
!
!  Compute the difference table for Richardson extrapolation.
!
    do ll = 2, it
      i = it + 1 - ll
      factor = real ( iwork2(i)**2, kind = 8 ) &
             / real ( iwork2(it)**2 - iwork2(i)**2, kind = 8 )
      table(i) = table(i+1) + ( table(i+1) - table(i) ) * factor
    end do

    result = table(1)
!
!  Terminate successfully if the estimated error is acceptable.
!
    if ( abs ( result - result_old ) <= abs ( result * ( tol + rnderr ) ) ) then
      ind = 1
      exit
    end if
!
!  Terminate unsuccessfully if the iteration limit has been reached.
!
    if ( it_max <= it ) then
      ind = -1
      exit
    end if
!
!  Prepare for another step.
!
    result_old = result

    it = it + 1

    iwork2(it) = int ( 1.5D+00 * real ( iwork2(it-1), kind = 8 ) )

    sub_num2(1:dim_num) = int ( 1.5D+00 * real ( sub_num2(1:dim_num), kind = 8 ) )

  end do

  return
end
subroutine sample_nd ( func, k1, k2, dim_num, est1, err1, dev1, est2, &
  err2, dev2, eval_num )

!*****************************************************************************80
!
!! SAMPLE_ND estimates a multidimensional integral using sampling.
!
!  Discussion:
!
!    This routine computes two sequences of integral estimates, EST1
!    and EST2, for indices K going from K1 to K2.  These estimates are
!    produced by the generation of 'random' abscissas in the region.
!    The process can become very expensive if high accuracy is needed.
!
!    The total number of function evaluations is
!    4*(K1**DIM_NUM+(K1+1)**DIM_NUM+...+(K2-1)**DIM_NUM+K2**DIM_NUM), and K2
!    should be chosen so as to make this quantity reasonable.
!    In most situations, EST2(K) are much better estimates than EST1(K).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2007
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, a routine which evaluates
!    the function to be integrated, of the form:
!      function func ( dim_num, x )
!      integer ( kind = 4 ) dim_num
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(dim_num)
!      func = ...
!      return
!      end
!
!    Input, integer ( kind = 4 ) K1, the beginning index for the iteration.
!    1 <= K1 <= K2.
!
!    Input, integer ( kind = 4 ) K2, the final index for the iteration.  K1 <= K2.
!    Increasing K2 increases the accuracy of the calculation,
!    but vastly increases the work and running time of the code.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM <= 10.
!
!    Output, real ( kind = 8 ) EST1(K2).  Entries K1 through K2 contain
!    successively better estimates of the integral.
!
!    Output, real ( kind = 8 ) ERR1(K2).  Entries K1 through K2 contain
!    the corresponding estimates of the integration errors.
!
!    Output, real ( kind = 8 ) DEV1(K2).  Entries K1 through K2 contain
!    estimates of the reliability of the the integration.
!    If consecutive values DEV1(K) and DEV1(K+1) do not differ
!    by more than 10 percent, then ERR1(K) can be taken as
!    a reliable upper bound on the difference between EST1(K)
!    and the true value of the integral.
!
!    Output, real ( kind = 8 ) EST2(K2).  Entries K2 through K2 contain
!    successively better estimates of the integral.
!
!    Output, real ( kind = 8 ) ERR2(K2).  Entries K2 through K2 contain
!    the corresponding estimates of the integration errors.
!
!    Output, real ( kind = 8 ) DEV2(K2).  Entries K2 through K2 contain
!    estimates of the reliability of the the integration.
!    If consecutive values DEV2(K) and DEV2(K+2) do not differ
!    by more than 10 percent, then ERR2(K) can be taken as
!    a reliable upper bound on the difference between EST2(K)
!    and the true value of the integral.
!
!    Output, integer ( kind = 4 ) EVAL_NUM, the number of function evaluations.
!
  implicit none

  integer ( kind = 4 ) k2
  integer ( kind = 4 ), parameter :: dim_max = 10
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) ak
  real ( kind = 8 ) ak1
  real ( kind = 8 ) akn
  real ( kind = 8 ), save, dimension ( dim_max ) :: al =  (/ &
       0.4142135623730950D+00, &
       0.7320508075688773D+00, &
       0.2360679774997897D+00, &
       0.6457513110645906D+00, &
       0.3166247903553998D+00, &
       0.6055512754639893D+00, &
       0.1231056256176605D+00, &
       0.3589989435406736D+00, &
       0.7958315233127195D+00, &
       0.3851648071345040D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ) be(dim_num)
  real ( kind = 8 ) bk
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dev1(k2)
  real ( kind = 8 ) dev2(k2)
  real ( kind = 8 ) dex(dim_num)
  real ( kind = 8 ) err1(k2)
  real ( kind = 8 ) err2(k2)
  real ( kind = 8 ) est1(k2)
  real ( kind = 8 ) est2(k2)
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) g
  real ( kind = 8 ) ga(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) key
  logical more
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ) p4(dim_num)
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) t
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4

  eval_num = 0
!
!  Check input
!
  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SAMPLE_ND - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM must be at least 1,'
    write ( *, '(a,i8)' ) '  but DIM_NUM = ', dim_num
    stop
  end if

  if ( dim_max < dim_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SAMPLE_ND - Fatal error!'
    write ( *, '(a,i8)' ) '  DIM_NUM must be no more than DIM_MAX = ', dim_max
    write ( *, '(a,i8)' ) '  but DIM_NUM = ', dim_num
    stop
  end if

  if ( k1 < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SAMPLE_ND - Fatal error!'
    write ( *, '(a,i8)' ) '  K1 must be at least 1, but K1 = ', k1
    stop
  end if

  if ( k2 < k1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SAMPLE_ND - Fatal error!'
    write ( *, '(a)' ) '  K1 may not be greater than K2, but '
    write ( *, '(a,i8)' ) '  K1 = ', k1
    write ( *, '(a,i8)' ) '  K2 = ', k2
    stop
  end if

  be(1:dim_num) = al(1:dim_num)
  ga(1:dim_num) = al(1:dim_num)
  dex(1:dim_num) = 0.0D+00

  do k = k1, k2

    ak = real ( k, kind = 8 )
    key = 0
    ak1 = ak - 1.1D+00
    s1 = 0.0D+00
    d1 = 0.0D+00
    s2 = 0.0D+00
    d2 = 0.0D+00
    akn = ak**dim_num
    t = sqrt ( ak**dim_num ) * ak
    bk = 1.0D+00 / ak

    do

      key = key + 1

      if ( key /= 1 ) then

        key = key - 1
        more = .false.

        do j = 1, dim_num

          if ( dex(j) <= ak1 ) then
            dex(j) = dex(j) + 1.0D+00
            more = .true.
            exit
          end if

          dex(j) = 0.0D+00

        end do

        if ( .not. more ) then
          exit
        end if

      end if

      do i = 1, dim_num

        b = be(i) + al(i)
        if ( 1.0D+00 < b ) then
          b = b - 1.0D+00
        end if

        g = ga(i) + b
        if ( 1.0D+00 < g ) then
          g = g - 1.0D+00
        end if

        be(i) = b + al(i)
        if ( 1.0D+00 < be(i) ) then
          be(i) = be(i) - 1.0D+00
        end if

        ga(i) = be(i) + g
        if ( 1.0D+00 < ga(i) ) then
          ga(i) = ga(i) - 1.0D+00
        end if

        p1(i) = ( dex(i) + g ) * bk
        p2(i) = ( dex(i) + 1.0D+00 - g ) * bk
        p3(i) = ( dex(i) + ga(i) ) * bk
        p4(i) = ( dex(i) + 1.0D+00 - ga(i) ) * bk

      end do

      y1 = func ( dim_num, p1 )
      eval_num = eval_num + 1
!
!  There may be an error in the next two lines,
!  but oddly enough, that is how the original reads
!
      y3 = func ( dim_num, p2 )
      eval_num = eval_num + 1
      y2 = func ( dim_num, p3 )
      eval_num = eval_num + 1
      y4 = func ( dim_num, p4 )
      eval_num = eval_num + 1

      s1 = s1 + y1 + y2
      d1 = d1 + ( y1 - y2 )**2
      s2 = s2 + y3 + y4
      d2 = d2 + ( y1 + y3 - y2 - y4 )**2

    end do

    est1(k) = 0.5D+00 * s1 / akn
    err1(k) = 1.5D+00 * sqrt ( d1 ) / akn
    dev1(k) = err1(k) * t
    est2(k) = 0.25D+00 * ( s1 + s2 ) / akn
    err2(k) = 0.75D+00 * sqrt ( d2 ) / akn
    dev2(k) = err2(k) * t * ak

  end do

  return
end
subroutine sum2_nd ( func, xtab, weight, order, dim_num, result, eval_num )

!*****************************************************************************80
!
!! SUM2_ND estimates a multidimensional integral using a product rule.
!
!  Discussion:
!
!    The routine uses a product rule supplied by the user.
!
!    The region may be a product of any combination of finite,
!    semi-infinite, or infinite intervals.
!
!    For each factor in the region, it is assumed that an integration
!    rule is given, and hence, the region is defined implicitly by
!    the integration rule chosen.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, a routine which evaluates
!    the function to be integrated, of the form:
!      function func ( dim_num, x )
!      integer ( kind = 4 ) dim_num
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(dim_num)
!      func = ...
!      return
!      end
!
!    Input, real ( kind = 8 ) XTAB(DIM_NUM,ORDER_MAX).  XTAB(I,J) is the
!    I-th abscissa of the J-th rule.
!
!    Input, real ( kind = 8 ) WEIGHT(DIM_NUM,ORDER_MAX).  WEIGHT(I,J) is the
!    I-th weight for the J-th rule.
!
!    Input, integer ( kind = 4 ) ORDER(DIM_NUM).  ORDER(I) is the number of
!    abscissas to be used in the J-th rule.  ORDER(I) must be
!    greater than 0 and less than or equal to ORDER_MAX.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
!    Output, integer ( kind = 4 ) EVAL_NUM, the number of function evaluations.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwork(dim_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) order(dim_num)
  real ( kind = 8 ) result
  real ( kind = 8 ) w1
  real ( kind = 8 ) weight(dim_num,*)
  real ( kind = 8 ) work(dim_num)
  real ( kind = 8 ) xtab(dim_num,*)
!
!  Default values.
!
  result = 0.0D+00
  eval_num = 0

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUM2_ND - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  do i = 1, dim_num

    if ( order(i) < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUM2_ND - Fatal error!'
      write ( *, '(a)' ) '  ORDER(I) < 1.'
      write ( *, '(a,i8)' ) '  For I = ', i
      write ( *, '(a,i8)' ) '  ORDER(I) = ', order(i)
      stop
    end if

  end do

  iwork(1:dim_num) = 1

  do

    k = 1

    w1 = 1.0D+00
    do i = 1, dim_num
      m1 = iwork(i)
      work(i) = xtab(i,m1)
      w1 = w1 * weight(i,m1)
    end do

    result = result + w1 * func ( dim_num, work )
    eval_num = eval_num + 1

    do while ( iwork(k) == order(k) )

      iwork(k) = 1
      k = k + 1

      if ( dim_num < k ) then
        return
      end if

    end do

    iwork(k) = iwork(k) + 1

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
subroutine tuple_next ( m1, m2, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT computes the next element of a tuple space.
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between M1 and M2.  The elements are produced one at a time.
!    The first element is
!      (M1,M1,...,M1),
!    the second element is
!      (M1,M1,...,M1+1),
!    and the last element is
!      (M2,M2,...,M2)
!    Intermediate elements are produced in lexicographic order.
!
!  Example:
!
!    N = 2, M1 = 1, M2 = 3
!
!    INPUT        OUTPUT
!    -------      -------
!    Rank  X      Rank   X
!    ----  ---    -----  ---
!    0     * *    1      1 1
!    1     1 1    2      1 2
!    2     1 2    3      1 3
!    3     1 3    4      2 1
!    4     2 1    5      2 2
!    5     2 2    6      2 3
!    6     2 3    7      3 1
!    7     3 1    8      3 2
!    8     3 2    9      3 3
!    9     3 3    0      0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M1, M2, the minimum and maximum entries.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input/output, integer ( kind = 4 ) RANK, counts the elements.
!    On first call, set RANK to 0.  Thereafter, the output value of RANK
!    will indicate the order of the element returned.  When there are no
!    more elements, RANK will be returned as 0.
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( m2 < m1 ) then
    rank = 0
    return
  end if

  if ( rank <= 0 ) then

    x(1:n) = m1
    rank = 1

  else

    rank = rank + 1
    i = n

    do

      if ( x(i) < m2 ) then
        x(i) = x(i) + 1
        exit
      end if

      x(i) = m1

      if ( i == 1 ) then
        rank = 0
        x(1:n) = m1
        exit
      end if

      i = i - 1

    end do

  end if

  return
end
