subroutine bandwidth_mesh ( element_order, element_num, element_node, &
  ml, mu, m )

!*****************************************************************************80
!
!! BANDWIDTH_MESH: bandwidth of finite element mesh.
!
!  Discussion:
!
!    The quantity computed here is the "geometric" bandwidth determined
!    by the finite element mesh alone.
!
!    If a single finite element variable is associated with each node
!    of the mesh, and if the nodes and variables are numbered in the
!    same way, then the geometric bandwidth is the same as the bandwidth
!    of a typical finite element matrix.
!
!    The bandwidth M is defined in terms of the lower and upper bandwidths:
!
!      M = ML + 1 + MU
!
!    where 
!
!      ML = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but earlier column,
!
!      MU = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but later column.
!
!    Because the finite element node adjacency relationship is symmetric,
!    we are guaranteed that ML = MU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM); 
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Output, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths of 
!    the matrix.
!
!    Output, integer ( kind = 4 ) M, the bandwidth of the matrix.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) global_i
  integer ( kind = 4 ) global_j
  integer ( kind = 4 ) local_i
  integer ( kind = 4 ) local_j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu

  ml = 0
  mu = 0
  
  do element = 1, element_num

    do local_i = 1, element_order
      global_i = element_node(local_i,element)

      do local_j = 1, element_order
        global_j = element_node(local_j,element)

        mu = max ( mu, global_j - global_i )
        ml = max ( ml, global_i - global_j )

      end do
    end do
  end do

  m = ml + 1 + mu

  return
end
subroutine legendre_com ( norder, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_COM computes abscissas and weights for Gauss-Legendre quadrature.
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 1998
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NORDER, the order of the rule.
!    NORDER must be greater than 0.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) norder

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2pn
  real ( kind = 8 ) d3pn
  real ( kind = 8 ) d4pn
  real ( kind = 8 ) dp
  real ( kind = 8 ) dpn
  real ( kind = 8 ) e1
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iback
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp1mi
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) nmove
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00 
  real ( kind = 8 ) pk
  real ( kind = 8 ) pkm1
  real ( kind = 8 ) pkp1
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x0
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) xtemp
  real ( kind = 8 ) weight(norder)

  if ( norder < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of NORDER = ', norder
    stop
  end if
 
  e1 = real ( norder * ( norder + 1 ), kind = 8 )
 
  m = ( norder + 1 ) / 2
 
  do i = 1, ( norder + 1 ) / 2
 
    mp1mi = m + 1 - i
    t = pi * real ( 4 * i - 1,      kind = 8 ) &
           / real ( 4 * norder + 2, kind = 8 )
    x0 = cos(t) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 / &
      real ( norder, kind = 8 ) ) / real ( 8 * norder * norder, kind = 8 ) )
 
    pkm1 = 1.0D+00
    pk = x0

    do k = 2, norder
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do
 
    d1 = real ( norder, kind = 8 ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0 * x0 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 * x0 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) / &
      ( 1.0D+00 - x0 * x0 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) &
      / ( 1.0D+00 - x0 * x0 )

    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn &
      / ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )

    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )

    h = h - p / dp
 
    xtemp = x0 + h

    xtab(mp1mi) = xtemp
 
    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

    weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp * xtemp ) / ( fx * fx )
 
  end do
 
  if ( mod ( norder, 2 ) == 1 ) then
    xtab(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = int ( ( norder + 1 ) / 2 )
  ncopy = norder - nmove

  do i = 1, nmove
    iback = norder + 1 - i
    xtab(iback) = xtab(iback-ncopy)
    weight(iback) = weight(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, norder - nmove
    xtab(i) = - xtab(norder+1-i)
    weight(i) = weight(norder+1-i)
  end do
 
  return
end
subroutine local_basis_1d ( order, node_x, x, phi )

!*****************************************************************************80
!
!! LOCAL_BASIS_1D evaluates the basis functions in an element.
!
!  Discussion:
!
!    PHI(I)(X) = product ( J ~= I ) ( X         - NODE_X(I) ) 
!                                 / ( NODE_X(J) - NODE_X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the element.
!    0 <= ORDER.  ORDER = 1 means piecewise linear.
!
!    Input, real ( kind = 8 ) NODE_X(ORDER), the element nodes.  
!    These must be distinct.  Basis function I is 1 when X = NODE_X(I) 
!    and 0 when X is equal to any other node.
!
!    Input, real ( kind = 8 ) X, the point at which the basis functions are to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) PHI(ORDER), the basis functions.
!
  implicit none

  integer ( kind = 4 ) order

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) node_x(order)
  real ( kind = 8 ) phi(order)
  real ( kind = 8 ) x

  phi(1:order) = 1.0D+00

  do i = 1, order
    do j = 1, order
      if ( j /= i ) then
        phi(j) = ( phi(j) * ( x - node_x(i) ) ) / ( node_x(j) - node_x(i) )
      end if
    end do
  end do

  return
end
subroutine local_basis_prime_1d ( order, node_x, x, dphidx )

!*****************************************************************************80
!
!! LOCAL_BASIS_PRIME_1D evaluates the basis function derivatives in an element.
!
!  Discussion:
!
!    PHI(I)(X) = product ( J ~= I ) ( X - NODE_X(I) ) 
!                                 / ( NODE_X(J) - NODE_X(I) )
!
!    dPHIdx(I)(X) = sum ( J ~= I ) ( 1 / ( NODE_X(J) - NODE_X(I) ) *
!      product ( K ~= ( J, I ) ) ( X - NODE_X(I) ) / ( NODE_X(J) - NODE_X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the element.
!    0 <= ORDER.  ORDER = 1 means piecewise linear.
!
!    Input, real ( kind = 8 ) NODE_X(ORDER), the element nodes.  
!    These must be distinct.  Basis function I is 1 when X = NODE_X(I) 
!    and 0 when X is equal to any other node.
!
!    Input, real ( kind = 8 ) X, the point at which the basis functions are to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) PHI(ORDER), the basis functions.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) dphidx(order)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) node_x(order)
  real ( kind = 8 ) term
  real ( kind = 8 ) x

  dphidx(1:order) = 0.0D+00

  do i = 1, order
    do j = 1, order
      if ( j /= i ) then
        term = 1.0D+00 / ( node_x(j) - node_x(i) )
        do k = 1, order
          if ( k /= i .and. k /= j ) then
            term = term * ( x - node_x(i) ) / ( node_x(k) - node_x(i) )
          end if
        end do
        dphidx(i) = dphidx(i) + term
      end if
    end do
  end do

  return
end
subroutine local_fem_1d ( order, node_x, node_v, sample_num, sample_x, &
  sample_v )

!*****************************************************************************80
!
!! LOCAL_FEM_1D evaluates a local finite element function.
!
!  Discussion:
!
!    A local finite element function is a finite element function
!    defined over a single element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the element.
!    0 <= ORDER.  ORDER = 1 means piecewise linear.
!
!    Input, real ( kind = 8 ) NODE_X(ORDER), the element nodes.  
!    These must be distinct.  Basis function I is 1 when X = NODE_X(I) and 0 
!    when X is equal to any other node.
!
!    Input, real ( kind = 8 ) NODE_V(ORDER), the value of the finite element 
!    function at each node.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sample points.
!
!    Input, real ( kind = 8 ) SAMPLE_X(SAMPLE_NUM), the sample points at which 
!    the local finite element function is to be evaluated.
!
!    Output, real ( kind = 8 ) SAMPLE_V(SAMPLE_NUM), the values of the local 
!    finite element basis functions.
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) sample_num

  real ( kind = 8 ) node_v(order)
  real ( kind = 8 ) node_x(order)
  real ( kind = 8 ) phi(order)
  integer ( kind = 4 ) sample
  real ( kind = 8 ) sample_v(sample_num)
  real ( kind = 8 ) sample_x(sample_num)
  real ( kind = 8 ) x

  sample_v(1:sample_num) = 0.0D+00

  do sample = 1, sample_num

    x = sample_x(sample)
    call local_basis_1d ( order, node_x, x, phi )
    sample_v(sample) = dot_product ( node_v(1:order), phi(1:order) )

  end do

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
!    For now, the input quantity SEED is an integer variable.
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
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
    seed = seed + i4_huge
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
