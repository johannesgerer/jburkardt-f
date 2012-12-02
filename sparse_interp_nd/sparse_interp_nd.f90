subroutine cc_compute_points ( n, points )

!*****************************************************************************80
!
!! CC_COMPUTE_POINTS: abscissas of a Clenshaw Curtis rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) POINTS(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) points(n)

  if ( n < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CC_COMPUTE_POINTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  else if ( n == 1 ) then

    points(1) = 0.0D+00

  else

    do i = 1, n
      points(i) = cos ( real ( n - i, kind = 8 ) * pi &
                      / real ( n - 1, kind = 8 ) )
    end do

    points(1) = -1.0D+00
    if ( mod ( n, 2 ) == 1 ) then
      points((n+1)/2) = 0.0D+00
    end if
    points(n) = +1.0D+00

  end if

  return
end
subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose compositions are desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the composition.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer ( kind = 4 ) H, T, two internal parameters needed
!    for the computation.  The user should allocate space for these in the
!    calling program, include them in the calling sequence, but never alter
!    them!
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) h
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else
!
!  If the first entry A(1) is positive, then set H to zero,
!  so that when we increment H, it points to A(1); we will decrement A(1) by 1
!  and increment A(2).
!
    if ( 1 < t ) then
      h = 0
    end if
!
!  Otherwise, A(1) is 0.  Then by H + 1 is the entry we incremented last time.
!  Set H = H + 1, zero A(H), adding all but one of its value to A(1),
!  and incrementing A(H+1) by 1.
!
    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )

  return
end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K) as an I4.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end
function i4_mop ( i )

!*****************************************************************************80
!
!! I4_MOP returns the I-th power of -1 as an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, integer ( kind = 4 ) I4_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_mop

  if ( mod ( i, 2 ) == 0 ) then
    i4_mop = 1
  else
    i4_mop = -1
  end if

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb ) 

!*****************************************************************************80
!
!! LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the interpolation nodes.
!
!    Input, integer ( kind = 4 ) NI, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XI(NI), the evaluation points.
!
!    Output, real ( kind = 8 ) LB(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lb(ni,nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  
  do i = 1, ni
    do j = 1, nd
      lb(i,j) = product ( ( xi(i) - xd(1:j-1)  ) / ( xd(j) - xd(1:j-1)  ) ) &
              * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
    end do
  end do

  return
end
subroutine lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Output, real ( kind = 8 ) XD(M,ND), the points at which data was sampled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xd(m,nd)
!
!  Compute the data points.
!
  xd(1:m,1:nd) = 0.0D+00
  do i = 1, m
    call order_from_level_135 ( ind(i), n )
    allocate ( x_1d(1:n) )
    call cc_compute_points ( n, x_1d )
    x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                          + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
    call r8vec_direct_product ( i, n, x_1d, m, nd, xd )
    deallocate ( x_1d )
  end do

  return
end
subroutine lagrange_interp_nd_size2 ( m, ind, nd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Output, integer ( kind = 4 ) ND, the number of points in the product grid.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
!
!  Determine the number of data points.
!
  nd = 1
  do i = 1, m
    call order_from_level_135 ( ind(i), n )
    nd = nd * n
  end do

  return
end
subroutine lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zi )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!
!    Input, integer ( kind = 4 ) NI, the number of points at which the 
!    interpolant is to be evaluated.
!
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xi(m,ni)
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)

  do j = 1, ni

    w(1:nd) = 1.0D+00

    do i = 1, m
      call order_from_level_135 ( ind(i), n )
      allocate ( x_1d(1:n) )
      allocate ( value(1:n) )
      call cc_compute_points ( n, x_1d )
      x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                            + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
      call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
      call r8vec_direct_product2 ( i, n, value, m, nd, w )
      deallocate ( value )
      deallocate ( x_1d )
    end do

    zi(j) = dot_product ( w, zd )

  end do

  return
end
subroutine order_from_level_135 ( l, n )

!*****************************************************************************80
!
!! ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.
!
!  Discussion:
!
!    Clenshaw Curtis rules, and some others, often use the following
!    scheme:
!
!    L: 0  1  2  3   4   5
!    N: 1  3  5  9  17  33 ... 2^L+1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the level, which should be 0 or greater.
!
!    Output, integer ( kind = 4 ) N, the order.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  if ( l < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ORDER_FROM_LEVEL_135 - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of L!'
    stop
  else if ( l == 0 ) then
    n = 1
  else
    n = ( 2 ** l ) + 1
  end if

  return
end
subroutine smolyak_coefficients ( l_max, m, c, w )

!*****************************************************************************80
!
!! SMOLYAK_COEFFICIENTS returns the Smolyak coefficients and counts.
!
!  Discussion:
!
!    The Smolyak sparse interpolant can be written as:
!
!      A(L,M)(X) = sum ( L-M+1 <= |L| <= L_max ) 
!        C(|L|) * g(l1)(x1) * g(l2)(x2) * ... * g(lm)(xm).
!
!    where:
!
!    * L=(l1,l2,...,lm) is a vector of M nonnegative integers;
!    * |L| is the sum of the entries of L;
!    * X=(x1,x2,...,xm) is an M-dimensional point in a product space;
!    * g(i)(xj) is the i-th 1-d interpolation function in dimension j;
!
!    Note that:
!
!    * W(|L|) will represent the number of distinct interpolants for which
!      the sublevel, or sum of the L vector entries, is |L|;
!
!    * the coefficients C and counts W will be zero for sublevels 
!      0 through L_MAX - M (and MATLAB indices 1 through L_MAX-M+1).
!
!    * it will be the case that W' * C = 1, essentially because the interpolant
!      to the identity function must be the identity function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L_MAX, the (maximum) level.
!    0 <= L_MAX.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!    1 <= M.
!
!    Output, integer ( kind = 4 ) C(0:L_MAX), the coefficients for objects 
!    at sublevels 0 through L_MAX.
!
!    Output, integer ( kind = 4 ) W(0:L_MAX), the number of objects at 
!    sublevels 0 through L_MAX.
!
  implicit none

  integer ( kind = 4 ) l_max

  integer ( kind = 4 ) c(0:l_max)
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) i4_mop
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_min
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w(0:l_max)

  l_min = max ( l_max - m + 1, 0 )

  c(0:l_min-1) = 0
  do l = l_min, l_max
    c(l) = i4_mop ( l_max - l ) * i4_choose ( m - 1, l_max - l )
  end do

  w(0:l_min-1) = 0
  do l = l_min, l_max
    w(l) = i4_choose ( l + m - 1, m - 1 )
  end do

  return
end
