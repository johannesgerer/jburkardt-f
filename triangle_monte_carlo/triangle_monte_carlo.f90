subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

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
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
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
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine reference_to_physical_t3 ( t, n, ref, phy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 3 physical triangle and a point 
!    (XSI,ETA) in the reference triangle, the routine computes the value 
!    of the corresponding image point (X,Y) in physical space.
!
!    This routine is also appropriate for an order 4 triangle,
!    as long as the fourth node is the centroid of the triangle.
!
!    This routine may also be appropriate for an order 6
!    triangle, if the mapping between reference and physical space
!    is linear.  This implies, in particular, that the sides of the
!    image triangle are straight and that the "midside" nodes in the
!    physical triangle are halfway along the sides of
!    the physical triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the coordinates of the vertices.  
!    The vertices are assumed to be the images of (0,0), (1,0) and 
!    (0,1) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) REF(2,N), points in the reference element.
!
!    Output, real ( kind = 8 ) PHY(2,N), corresponding points in the
!    physical element.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) phy(2,n)
  real ( kind = 8 ) ref(2,n)
  real ( kind = 8 ) t(2,3)

  do i = 1, 2
    phy(i,1:n) = t(i,1) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) ) &
               + t(i,2) *             ref(1,1:n)                &
               + t(i,3) *                          ref(2,1:n)
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
subroutine triangle_area ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA computes the area of a triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) AREA, the absolute area of the triangle.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) t(2,3)

  area = 0.5D+00 * abs ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end
subroutine triangle_integrand_01 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TRIANGLE_INTEGRAND_01 evaluates 1 integrand function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input, real ( kind = 8 ) P(2,P_NUM), the evaluation points.
!
!    Input, integer ( kind = 4 ) F_NUM, the number of integrands.
!
!    Output, real ( kind = 8 ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer ( kind = 4 ) f_num
  integer ( kind = 4 ) p_num

  real ( kind = 8 ) fp(f_num,p_num)
  real ( kind = 8 ) p(2,p_num)

  fp(1,1:p_num) = 1.0D+00

  return
end
subroutine triangle_integrand_02 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TRIANGLE_INTEGRAND_02 evaluates 2 integrand functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input, real ( kind = 8 ) P(2,P_NUM), the evaluation points.
!
!    Input, integer ( kind = 4 ) F_NUM, the number of integrands.
!
!    Output, real ( kind = 8 ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer ( kind = 4 ) f_num
  integer ( kind = 4 ) p_num

  real ( kind = 8 ) fp(f_num,p_num)
  real ( kind = 8 ) p(2,p_num)

  fp(1,1:p_num) = p(1,1:p_num)
  fp(2,1:p_num) = p(2,1:p_num)

  return
end
subroutine triangle_integrand_03 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TRIANGLE_INTEGRAND_03 evaluates 3 integrand functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input, real ( kind = 8 ) P(2,P_NUM), the evaluation points.
!
!    Input, integer ( kind = 4 ) F_NUM, the number of integrands.
!
!    Output, real ( kind = 8 ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer ( kind = 4 ) f_num
  integer ( kind = 4 ) p_num

  real ( kind = 8 ) fp(f_num,p_num)
  real ( kind = 8 ) p(2,p_num)

  fp(1,1:p_num) = p(1,1:p_num) * p(1,1:p_num)
  fp(2,1:p_num) = p(1,1:p_num) * p(2,1:p_num)
  fp(3,1:p_num) = p(2,1:p_num) * p(2,1:p_num)

  return
end
subroutine triangle_integrand_04 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TRIANGLE_INTEGRAND_04 evaluates 4 integrand functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input, real ( kind = 8 ) P(2,P_NUM), the evaluation points.
!
!    Input, integer ( kind = 4 ) F_NUM, the number of integrands.
!
!    Output, real ( kind = 8 ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer ( kind = 4 ) f_num
  integer ( kind = 4 ) p_num

  real ( kind = 8 ) fp(f_num,p_num)
  real ( kind = 8 ) p(2,p_num)

  fp(1,1:p_num) = p(1,1:p_num) * p(1,1:p_num) * p(1,1:p_num)
  fp(2,1:p_num) = p(1,1:p_num) * p(1,1:p_num) * p(2,1:p_num)
  fp(3,1:p_num) = p(1,1:p_num) * p(2,1:p_num) * p(2,1:p_num)
  fp(4,1:p_num) = p(2,1:p_num) * p(2,1:p_num) * p(2,1:p_num)

  return
end
subroutine triangle_integrand_05 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TRIANGLE_INTEGRAND_05 evaluates 5 integrand functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input, real ( kind = 8 ) P(2,P_NUM), the evaluation points.
!
!    Input, integer ( kind = 4 ) F_NUM, the number of integrands.
!
!    Output, real ( kind = 8 ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer ( kind = 4 ) f_num
  integer ( kind = 4 ) p_num

  real ( kind = 8 ) fp(f_num,p_num)
  real ( kind = 8 ) p(2,p_num)

  fp(1,1:p_num) = p(1,1:p_num)**4
  fp(2,1:p_num) = p(1,1:p_num)**3 * p(2,1:p_num)
  fp(3,1:p_num) = p(1,1:p_num)**2 * p(2,1:p_num)**2
  fp(4,1:p_num) = p(1,1:p_num)    * p(2,1:p_num)**3
  fp(5,1:p_num) =                   p(2,1:p_num)**4
  
  return
end
subroutine triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample, &
  triangle_integrand, seed, result )

!*****************************************************************************80
!
!! TRIANGLE_MONTE_CARLO applies the Monte Carlo rule to integrate a function.
!
!  Discussion:
!
!    The function f(x,y) is to be integrated over a triangle T.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Input, integer ( kind = 4 ) P_NUM, the number of sample points.
!
!    Input, integer ( kind = 4 ) F_NUM, the number of functions to integrate.
!
!    Input, external TRIANGLE_UNIT_SAMPLE, the sampling routine.
!
!    Input, external TRIANGLE_INTEGRAND, the integrand routine.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) RESULT(F_NUM), the approximate integrals.
!
  implicit none

  integer ( kind = 4 ) f_num
  integer ( kind = 4 ) p_num

  real ( kind = 8 ) area
  real ( kind = 8 ) fp(f_num,p_num)
  integer ( kind = 4 ) i
  real ( kind = 8 ) p(2,p_num)
  real ( kind = 8 ) p2(2,p_num)
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(2,3)
  external triangle_sample
  external triangle_integrand

  call triangle_area ( t, area )

  call triangle_unit_sample ( p_num, seed, p )

  call reference_to_physical_t3 ( t, p_num, p, p2 )

  call triangle_integrand ( p_num, p2, f_num, fp )

  do i = 1, f_num
    result(i) = area * sum ( fp(i,1:p_num) ) / real ( p_num, kind = 8 )
  end do

  return
end
subroutine triangle_unit_sample_01 ( p_num, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SAMPLE_01 selects points from the unit triangle.
!
!  Discussion:
!
!    The unit triangle has vertices (1,0), (0,1), (0,0).
!
!    Any point in the unit triangle CAN be chosen by this algorithm.
!
!    However, the points that are chosen tend to be clustered near
!    the centroid.
!
!    This routine is supplied as an example of "bad" sampling.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) P(2,P_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) p_num

  real ( kind = 8 ) e(3)
  real ( kind = 8 ) e_sum
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(2,p_num)
  integer ( kind = 4 ) seed

  do j = 1, p_num

    call r8vec_uniform_01 ( 3, seed, e )

    e_sum = sum ( e(1:3) )

    e(1:3) = e(1:3) / e_sum
!
!  We may take the values E(1:3) as being the barycentric
!  coordinates of the point.
!
    p(1:2,j) = e(1:2)

  end do

  return
end
subroutine triangle_unit_sample_02 ( p_num, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SAMPLE_02 selects points from the unit triangle.
!
!  Discussion:
!
!    The unit triangle has vertices (1,0), (0,1), (0,0).
!
!    The sampling is uniform.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) P(2,P_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) p_num

  integer ( kind = 4 ) j
  real ( kind = 8 ) r(2)
  real ( kind = 8 ) p(2,p_num)
  integer ( kind = 4 ) seed
!
!  Generate the points using barycentric coordinates.
!
  do j = 1, p_num

    call r8vec_uniform_01 ( 2, seed, r )

    if ( 1.0D+00 < sum ( r(1:2) ) ) then
      r(1:2) = 1.0D+00 - r(1:2)
    end if

    p(1:2,j) = r(1:2)

  end do

  return
end
subroutine triangle_unit_sample_03 ( p_num, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SAMPLE_03 selects points from the unit triangle.
!
!  Discussion:
!
!    The unit triangle has vertices (1,0), (0,1), (0,0).
!
!    This routine uses Turk's rule #1.
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
!    Greg Turk,
!    Generating Random Points in a Triangle,
!    in Graphics Gems,
!    edited by Andrew Glassner,
!    AP Professional, 1990, pages 24-28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) P(2,P_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) p_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(2,p_num)
  real ( kind = 8 ) r(2)
  integer ( kind = 4 ) seed
!
!  Generate the points using Turk's rule 1.
!
  do j = 1, p_num

    call r8vec_uniform_01 ( 2, seed, r )

    a = 1.0D+00            - sqrt ( r(2) )
    b = ( 1.0D+00 - r(1) ) * sqrt ( r(2) )
    c =             r(1)   * sqrt ( r(2) )

    p(1,j) = a
    p(2,j) = b

  end do

  return
end
subroutine triangle_unit_sample_04 ( p_num, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SAMPLE_04 selects points from the unit triangle.
!
!  Discussion:
!
!    The unit triangle has vertices (1,0), (0,1), (0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity 
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) P(2,P_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) p_num

  real ( kind = 8 ) e(3)
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(2,p_num)
  integer ( kind = 4 ) seed
!
!  The construction begins by sampling DIM_NUM+1 points from the
!  exponential distribution with parameter 1.
!
  do j = 1, p_num

    call r8vec_uniform_01 ( 3, seed, e )

    e(1:3) = - log ( e(1:3) )

    p(1:2,j) = e(1:2) / sum ( e(1:3) )

  end do

  return
end
