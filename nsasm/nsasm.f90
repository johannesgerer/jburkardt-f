subroutine assemble ( p, t, u, nt, np, Pr, Ir, Jc, L, mu )

!*****************************************************************************80
!
!! ASSEMBLE assembles the local stiffness and residual into global arrays.
!
!  Modified:
!
!    23 January 2011
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P[2][NP], the coordinates of the mesh nodes.
!
!    Input, integer T[6][NT], the indices of the nodes in each element.
!
!    Input, real ( kind = 8 ) U[2*NP+NP0+NE], the current solution vector.
!
!    Input, integer ( kind = 4 ) NT, the number of elements.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input/output, real ( kind = 8 ) Pr[NNZ], the values of the nonzero entries
!    of the sparse matrix.
!
!    Input, integer ( kind = 4 ) Ir[NNZ], the row indices of the nonzero entries
!    of the sparse matrix.
!
!    Input, integer ( kind = 4 ) Jc[N+1], points to the first element of each
!    column.
!
!    Input/output, real ( kind = 8 ) L[2*NP+NP0+NE], the residual vector.
!
!    Input, real ( kind = 8 ) MU, the kinematic viscosity.
!
  implicit none

  integer ( kind = 4 ), parameter :: ngp = 7
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nt

  real    ( kind = 8 ) gp(3,ngp)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ir(*)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) Jc(*)
  integer ( kind = 4 ) j2
  real    ( kind = 8 ) L(*)
  real    ( kind = 8 ) lK(15,15)
  real    ( kind = 8 ) lL(15)
  real    ( kind = 8 ) mu
  real    ( kind = 8 ) p(2,np)
  real    ( kind = 8 ) pr(*)
  integer ( kind = 4 ) ri
  integer ( kind = 4 ) rj
  real    ( kind = 8 ) sh1(3,ngp)
  real    ( kind = 8 ) sh1r(3,ngp)
  real    ( kind = 8 ) sh1s(3,ngp)
  real    ( kind = 8 ) sh2(6,ngp)
  real    ( kind = 8 ) sh2r(6,ngp)
  real    ( kind = 8 ) sh2s(6,ngp)
  integer ( kind = 4 ) t(6,nt)
  real    ( kind = 8 ) u(*)
  real    ( kind = 8 ) w(ngp)
!
!  Get the quadrature rule.
!
  call quad_rule ( ngp, w, gp )
!
!  Initialize the shape functions.
!
  call init_shape ( ngp, gp, sh1, sh1r, sh1s, sh2, sh2r, sh2s )
!
!  For each element, determine the local matrix and right hand side.
!
  do it = 1, nt

    call localKL ( ngp, w, sh1, sh1r, sh1s, sh2, sh2r, sh2s, p, t(1,it), &
      u, np, mu, lK, lL )
!
!  Add the local right hand side and matrix to the global data.
!
    do i = 1, 15

      i2 = mod ( i - 1, 6 ) + 1
      ri = t(i2,it) + np * ( i / 6 )
      L(ri) = L(ri) + lL(i)

      do j = 1, 15
        j2 = mod ( j - 1, 6 ) + 1
        rj = t(j2,it) + np * ( j / 6 ) - 1
        call sparse_set ( Pr, Ir, Jc, ri, rj, lK(i,j) )
      end do
    end do

  end do

  return
end
subroutine assemble_constr ( e, u, np, ne, n0, Pr, Ir, Jc, L )

!*****************************************************************************80
!
!! ASSEMBLE_CONSTR assembles the constraints.
!
!  Modified:
!
!    23 January 2011
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) E(3,NE), the node, variable type (0,1,2) and
!    value of each constraint.
!
!    Input, real ( kind = 8 ) U(2*NP+NP0+NE), the current solution vector.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, integer ( kind = 4 ) NE, the number of constraints.
!
!    Input, integer ( kind = 4 ) N0, the number of variables.
!
!    Input/output, real ( kind = 8 ) Pr(NNZ), the values of the nonzero entries
!    of the sparse matrix.
!
!    Input, integer ( kind = 4 ) Ir(NNZ), the row indices of the nonzero entries
!    of the sparse matrix.
!
!    Input, integer ( kind = 4 ) Jc(N+1), points to the first element of each
!    column.
!
!    Input/output, real ( kind = 8 ) L(2*NP+NP0+NE), the residual vector.
!    On output, the residual vector has been updated to include the constraints.
!
  implicit none

  integer ( kind = 4 ) ne

  real    ( kind = 8 ) e(3,ne)
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) Ir(*)
  integer ( kind = 4 ) Jc(*)
  real    ( kind = 8 ) l(*)
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) np
  real    ( kind = 8 ) Pr(*)
  integer ( kind = 4 ) ri
  integer ( kind = 4 ) rj
  real    ( kind = 8 ) u(*)
  real    ( kind = 8 ) value

  value = 1.0D+00

  do ie = 1, ne

    ri = n0 + ie
    rj = int ( e(2,ie) ) * np + int ( e(1,ie) )

    call sparse_set ( Pr, Ir, Jc, ri, rj, value )
    call sparse_set ( Pr, Ir, Jc, rj, ri, value )

    L(rj) = L(rj) + u(ri)
    L(ri) = u(rj) - e(3,ie)

  end do

  return
end
subroutine init_shape ( ngp, gp, sh1, sh1r, sh1s, sh2, sh2r, sh2s )

!*****************************************************************************80
!
!! INIT_SHAPE evaluates the shape functions at the quadrature points.
!
!  Modified:
!
!    21 January 2011
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NGP, the number of Gauss quadrature points.
!
!    Input, real ( kind = 8 ) GP(3,NGP), the Gauss quadrature points.
!
!    Output, real ( kind = 8 ) SH1(3,NGP), SH1R(3,NGP), SH1S(3,NGP), the
!    P1 shape function, and R and S derivatives, evaluated at the
!    Gauss points.
!
!    Output, real ( kind = 8 ) SH2(6,NGP), SH2R(6,NGP), SH2S(6,NGP), the
!    P2 shape function, and R and S derivatives, evaluated at the
!    Gauss points.
!
  implicit none

  integer ( kind = 4 ) ngp

  real    ( kind = 8 ) gp(3,ngp)
  integer ( kind = 4 ) j
  real    ( kind = 8 ) sh1(3,ngp)
  real    ( kind = 8 ) sh1r(3,ngp)
  real    ( kind = 8 ) sh1s(3,ngp)
  real    ( kind = 8 ) sh2(6,ngp)
  real    ( kind = 8 ) sh2r(6,ngp)
  real    ( kind = 8 ) sh2s(6,ngp)

  do j = 1, ngp

    sh1(1,j) = gp(3,j)
    sh1(2,j) = gp(1,j)
    sh1(3,j) = gp(2,j)

    sh1r(1,j) = -1.0D+00
    sh1r(2,j) =  1.0D+00
    sh1r(3,j) =  0.0D+00

    sh1s(1,j) = -1.0D+00
    sh1s(2,j) =  0.0D+00
    sh1s(3,j) =  1.0D+00

    sh2(1,j) = 1.0D+00 - 3.0D+00 * gp(1,j) - 3.0D+00 * gp(2,j) &
      + 2.0D+00 * gp(1,j) * gp(1,j) + 4.0D+00 * gp(1,j) * gp(2,j) &
      + 2.0D+00 * gp(2,j) * gp(2,j)
    sh2(2,j) = - gp(1,j) + 2.0D+00 * gp(1,j) * gp(1,j)
    sh2(3,j) = - gp(2,j) + 2.0D+00 * gp(2,j) * gp(2,j)
    sh2(4,j) = 4.0D+00 * gp(1,j) * gp(2,j)
    sh2(5,j) = 4.0D+00 * gp(2,j) - 4.0D+00 * gp(1,j) * gp(2,j) &
      - 4.0D+00 * gp(2,j) * gp(2,j)
    sh2(6,j) = 4.0D+00 * gp(1,j) - 4.0D+00 * gp(1,j) * gp(2,j) &
      - 4.0D+00 * gp(1,j) * gp(1,j)

    sh2r(1,j) = - 3.0D+00 + 4.0D+00 * gp(1,j) + 4.0D+00 * gp(2,j)
    sh2r(2,j) = - 1.0D+00 + 4.0D+00 * gp(1,j)
    sh2r(3,j) = 0.0D+00
    sh2r(4,j) = 4.0D+00 * gp(2,j)
    sh2r(5,j) = - 4.0D+00 * gp(2,j)
    sh2r(6,j) = 4.0D+00 - 8.0D+00 * gp(1,j) - 4.0D+00 * gp(2,j)

    sh2s(1,j) = - 3.0D+00 + 4.0D+00 * gp(1,j) + 4.0D+00 * gp(2,j)
    sh2s(2,j) =   0.0D+00
    sh2s(3,j) = - 1.0D+00 + 4.0D+00 * gp(2,j)
    sh2s(4,j) =   4.0D+00 * gp(1,j)
    sh2s(5,j) =   4.0D+00 - 8.0D+00 * gp(2,j) - 4.0D+00 * gp(1,j)
    sh2s(6,j) = - 4.0D+00 * gp(1,j)

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an I4VEC into an descending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
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
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) t

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  t    = a(1)
  a(1) = a(n)
  a(n) = t
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    t    = a(1)
    a(1) = a(n1)
    a(n1) = t

  end do

  return
end
subroutine localKL ( ngp, w, sh1, sh1r, sh1s, sh2, sh2r, sh2s, p, tt, &
  u0, np, mu, lK, lL )

!*****************************************************************************80
!
!! LOCALKL assembles the local stiffness matrix and residual.
!
!  Modified:
!
!    22 January 2011
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NGP, the number of quadrature points.
!
!    Input, real ( kind = 8 ) W(NGP), the quadrature weights.
!
!    Input, real ( kind = 8 ) SH1(3,NGP), SH1R(3,NGP), SH1S(3,NGP), the
!    P1 shape function, and R and S derivatives, evaluated at the
!    Gauss points.
!
!    Input, real ( kind = 8 ) SH2(6,NGP), SH2R(6,NGP), SH2S(6,NGP), the
!    P2 shape function, and R and S derivatives, evaluated at the
!    Gauss points.
!
!    Input, real ( kind = 8 ) P(2,NP), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) TT(6), the nodes for the current element.
!
!    Input, real ( kind = 8 ) U0(2*NP+NP0+NE), called "U" elsewhere in the
!    program, but renamed "U0" here because we want to use "U" for the
!    horizontal velocity.  U0 contains the current finite element solution.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, real ( kind = 8 ) MU, the kinematic viscosity.
!
!    Output, real ( kind = 8 ) lK(15,15), the local stiffness matrix.
!
!    Output, real ( kind = 8 ) lL(15), the local residual vector.
!
  implicit none

  integer ( kind = 4 ) ngp
  integer ( kind = 4 ) np

  integer ( kind = 4 ) i
  integer ( kind = 4 ) igp
  integer ( kind = 4 ) j
  real    ( kind = 8 ) Jdet
  real    ( kind = 8 ) Jinv(2,2)
  real    ( kind = 8 ) lK(15,15)
  real    ( kind = 8 ) lL(15)
  real    ( kind = 8 ) mu
  real    ( kind = 8 ) mul
  real    ( kind = 8 ) p(3,np)
  real    ( kind = 8 ) px
  real    ( kind = 8 ) py
  real    ( kind = 8 ) sh1(3,ngp)
  real    ( kind = 8 ) sh1r(3,ngp)
  real    ( kind = 8 ) sh1s(3,ngp)
  real    ( kind = 8 ) sh2(6,ngp)
  real    ( kind = 8 ) sh2r(6,ngp)
  real    ( kind = 8 ) sh2s(6,ngp)
  real    ( kind = 8 ) sh1x(3)
  real    ( kind = 8 ) sh1y(3)
  real    ( kind = 8 ) sh2x(6)
  real    ( kind = 8 ) sh2y(6)
  integer ( kind = 4 ) tt(6)
  real    ( kind = 8 ) u
  real    ( kind = 8 ) u0(*)
  real    ( kind = 8 ) ux
  real    ( kind = 8 ) uy
  real    ( kind = 8 ) v
  real    ( kind = 8 ) vx
  real    ( kind = 8 ) vy
  real    ( kind = 8 ) w(ngp)
  real    ( kind = 8 ) xr
  real    ( kind = 8 ) xs
  real    ( kind = 8 ) yr
  real    ( kind = 8 ) ys
!
!  Zero out lK and lL.
!
  lK(1:15,1:15) = 0.0D+00
  lL(1:15) = 0.0D+00

  do igp = 1, ngp
!
!  Compute the Jacobian.
!
    xr = 0.0D+00
    xs = 0.0D+00
    yr = 0.0D+00
    ys = 0.0D+00

    do i = 1, 6
      xr = xr + sh2r(i,igp) * p(1,tt(i))
      xs = xs + sh2s(i,igp) * p(1,tt(i))
      yr = yr + sh2r(i,igp) * p(2,tt(i))
      ys = ys + sh2s(i,igp) * p(2,tt(i))
    end do

    Jdet = xr * ys - xs * yr

    Jinv(1,1) =  ys / Jdet
    Jinv(1,2) = -xs / Jdet
    Jinv(2,1) = -yr / Jdet
    Jinv(2,2) =  xr / Jdet
!
!  Set the X and Y derivatives of the shape functions.
!
    do i = 1, 3
      sh1x(i) = sh1r(i,igp) * Jinv(1,1) + sh1s(i,igp) * Jinv(2,1)
      sh1y(i) = sh1r(i,igp) * Jinv(1,2) + sh1s(i,igp) * Jinv(2,2)
    end do

    do i = 1, 6
      sh2x(i) = sh2r(i,igp) * Jinv(1,1) + sh2s(i,igp) * Jinv(2,1)
      sh2y(i) = sh2r(i,igp) * Jinv(1,2) + sh2s(i,igp) * Jinv(2,2)
    end do
!
!  Solution and derivatives.
!
    u  = 0.0D+00
    ux = 0.0D+00
    uy = 0.0D+00

    v  = 0.0D+00
    vx = 0.0D+00
    vy = 0.0D+00

    do i = 1, 6
      u  = u  + sh2(i,igp) * u0(tt(i))
      ux = ux + sh2x(i)    * u0(tt(i))
      uy = uy + sh2y(i)    * u0(tt(i))

      v  = v  + sh2(i,igp) * u0(np+tt(i))
      vx = vx + sh2x(i)    * u0(np+tt(i))
      vy = vy + sh2y(i)    * u0(np+tt(i))
    end do

    px = 0.0D+00
    py = 0.0D+00
    do i = 1, 3
      px = px + sh1x(i) * u0(2*np+tt(i))
      py = py + sh1y(i) * u0(2*np+tt(i))
    end do
!
!  Local K.
!
    mul = w(igp) * Jdet / 2.0D+00

    do i = 1, 6

      lL(i)   = lL(i)   + mul * ( &
        ( u * ux + v * uy + px ) * sh2(i,igp) &
        + mu * ( ux * sh2x(i) + uy * sh2y(i) ) )

      lL(6+i) = lL(6+i) + mul * ( &
        ( u * vx + v * vy + py ) * sh2(i,igp) &
        + mu * ( vx * sh2x(i) + vy * sh2y(i) ) )

      do j = 1, 6

        lK(i,j) = lK(i,j) &
          + mu * ( sh2x(i) * sh2x(j) + sh2y(i) * sh2y(j) ) * mul

        lK(6+i,6+j) = lK(6+i,6+j) &
          + mu * ( sh2x(i) * sh2x(j) + sh2y(i) * sh2y(j) ) * mul

        lK(i,j) = lK(i,j) &
          + ( u  * sh2(i,igp) * sh2x(j) &
          +   v  * sh2(i,igp) * sh2y(j) ) * mul

        lK(6+i,6+j) = lK(6+i,6+j) &
          + ( u  * sh2(i,igp) * sh2x(j) &
          +   v  * sh2(i,igp) * sh2y(j) ) * mul

        lK(i,j) = lK(i,j) &
          +  ( ux * sh2(i,igp) * sh2(j,igp) ) * mul

        lK(i,6+j) = lK(i,6+j) &
          + ( uy * sh2(i,igp) * sh2(j,igp) ) * mul

        lK(6+i,j) = lK(6+i,j) &
          + ( vx * sh2(i,igp) * sh2(j,igp) ) * mul

        lK(6+i,6+j) = lK(6+i,6+j) &
          + ( vy * sh2(i,igp) * sh2(j,igp) ) * mul

      end do

      do j = 1, 3
        lK(i,12+j)   = lK(i,12+j)   + ( sh2(i,igp) * sh1x(j) ) * mul
        lK(6+i,12+j) = lK(6+i,12+j) + ( sh2(i,igp) * sh1y(j) ) * mul
      end do

    end do
!
!  Local L.
!
    do i = 1, 3
      lL(12+i) = lL(12+i) + ( ux + vy ) * sh1(i,igp) * mul
      do j = 1, 6
        lK(12+i,j) =   lK(12+i,j)   + ( sh1(i,igp) * sh2x(j) ) * mul
        lK(12+i,6+j) = lK(12+i,6+j) + ( sh1(i,igp) * sh2y(j) ) * mul
      end do
    end do

  end do

  return
end
subroutine quad_rule ( ngp, w, gp )

!*****************************************************************************80
!
!! QUAD_RULE returns the points and weights of a quadrature rule.
!
!  Discussion:
!
!    At the moment, only a 7-point rule is available.
!
!  Modified:
!
!    22 January 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NGP, the number of quadrature points.
!
!    Output, real ( kind = 8 ) W(NGP), the quadrature weights.
!
!    Output, real ( kind = 8 ) GP(3,NGP), the quadrature points.
!
  implicit none

  integer ( kind = 4 ) ngp

  real ( kind = 8 ) gp(3,ngp)
  real ( kind = 8 ), save :: gp_7(3,7) = reshape ( (/ &
    0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333334D+00, &
    0.059715871789770D+00, 0.470142064105115D+00, 0.470142064105115D+00, &
    0.470142064105115D+00, 0.059715871789770D+00, 0.470142064105115D+00, &
    0.470142064105115D+00, 0.470142064105115D+00, 0.059715871789770D+00, &
    0.797426985353087D+00, 0.101286507323457D+00, 0.101286507323457D+00, &
    0.101286507323457D+00, 0.797426985353087D+00, 0.101286507323457D+00, &
    0.101286507323457D+00, 0.101286507323457D+00, 0.797426985353087D+00 /), &
  (/ 3, 7 /) )
  real ( kind = 8 ) w(ngp)
  real ( kind = 8 ), save :: w_7(7) = (/ &
    0.225000000000000D+00, &
    0.132394152788506D+00, &
    0.132394152788506D+00, &
    0.132394152788506D+00, &
    0.125939180544827D+00, &
    0.125939180544827D+00, &
    0.125939180544827D+00 /)

  if ( ngp == 7 ) then
    w(1:ngp)      = w_7(1:ngp)
    gp(1:3,1:ngp) = gp_7(1:3,1:ngp)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUAD_RULE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected input value of NGP.'
    stop
  end if

  return
end
subroutine r8sp_print_some ( m, n, nz_num, row, col, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8SP_PRINT_SOME prints some of an R8SP matrix.
!
!  Discussion:
!
!    This version of R8SP_PRINT_SOME has been specifically modified to allow,
!    and correctly handle, the case in which a single matrix location
!    A(I,J) is referenced more than once by the sparse matrix structure.
!    In such cases, the routine prints out the sum of all the values.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) nz_num

  real    ( kind = 8 ) a(nz_num)
  real    ( kind = 8 ) aij(incx)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '
    write ( *, '(''  Col:  '',5(i7,7x))' ) ( j, j = j2lo, j2hi )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      aij(1:inc) = 0.0D+00
!
!  Is matrix entry K actually the value of A(I,J), with J2LO <= J <= J2HI?
!  Because MATLAB seems to allow for multiple (I,J,A) entries, we have
!  to sum up what we find.
!
      do k = 1, nz_num

        if ( i == row(k) .and. &
             j2lo <= col(k) .and. &
             col(k) <= j2hi ) then

          j2 = col(k) - j2lo + 1
          aij(j2) = aij(j2) + a(k)

        end if

      end do

      if ( any ( aij(1:inc) /= 0.0D+00 ) ) then
        write ( *, '(i5,1x,5g14.6)' ) i, aij(1:inc)
      end if

    end do

  end do

  return
end
subroutine sparse_set ( pr, ir, jc, ri, rj, val )

!*****************************************************************************80
!
!! SPARSE_SET increments an entry of the sparse matrix.
!
!  Discussion:
!
!    We know RJ, the column in which the entry occurs.  We know that entries
!    in column RJ occur in positions K1 = Jc(RJ) through K2 = JC(RJ+1)-1.
!    We may assume that the corresponding entries Ir(K1) through Ir(K2)
!    are ascending sorted and that one of these entries is equal to RI.
!
!    We now simply use binary search on Ir to locate the index K for which
!    Ir(K) = RI, and then increment Pr(K).
!
!  Modified:
!
!    21 January 2011
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) PR(NNZ), the values of the nonzero
!    entries of the sparse matrix.
!
!    Input, integer ( kind = 4 ) IR(NNZ), the row indices of the nonzero
!    entries of the sparse matrix.
!
!    Input, integer ( kind = 4 ) JC(N+1), points to the location in IR and
!    PR of the first element of each column.
!
!    Input, integer ( kind = 4 ) RI, RJ, the row and column index of the
!    matrix entry that is to be incremented.
!
!    Input, real ( kind = 8 ) VAL, the amount to be added to the matrix entry.
!
  implicit none

  integer ( kind = 4 ) cr
  integer ( kind = 4 ) ir(*)
  integer ( kind = 4 ) jc(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real    ( kind = 8 ) pr(*)
  integer ( kind = 4 ) ri
  integer ( kind = 4 ) rj
  real    ( kind = 8 ) val

  cr = ri
!
!  Set the range K1 <= K2 of entries in Ir to be searched.
!
  k1 = jc(rj)
  k2 = jc(rj+1) - 1
!
!  We seek the index K so that IR(K) = RI.
!
  do

    k = ( k1 + k2 ) / 2

    if ( cr < ir(k) ) then
      k2 = k - 1
    else
      k1 = k + 1
    end if

    if ( cr == ir(k) ) then
      exit
    end if

    if ( k2 < k1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSE_SET - Fatal error!'
      write ( *, '(a)' ) '  Could not locate the sparse matrix storage'
      write ( *, '(a)' ) '  index K for logical matrix entry (RI,RJ).'
      stop
    end if

  end do
!
!  Increment the matrix entry.
!
  pr(k) = pr(k) + val

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
