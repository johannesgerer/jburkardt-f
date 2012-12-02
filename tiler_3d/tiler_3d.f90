program main

!*****************************************************************************80
!
!! MAIN is the main program for TILER_3D.
!
!  Discussion:
!
!    TILER_3D illustrates the use of 3D blending.
!
!    This main program works in a 3D rectangular space we can think
!    of as "UVW" space.  We are interested in the space of data values
!    bounded by (umin,vmin,wmin) and (umax,vmax,wmax), and we plan to
!    divide this up, using indices (I,J,K), into ni by nj by nk sub-boxes.
!
!    The code below considers each sub-box indexed by (I,J,K) and determines
!    the values (u0,v0,w0) and (u1,v1,w1) that characterize its corners.
!    The coordinates of this sub-box and the coordinates of the big box
!    are then passed to SUB_BOX_TILER_3D.
!
!    The picture would be a LOT more interesting if the boundary were
!    a bit more wiggly, there were more sub-boxes, and the object in
!    each sub-box had more parts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon, Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, pages 461-477, 1973.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: ni = 3
  integer ( kind = 4 ), parameter :: nj = 3
  integer ( kind = 4 ), parameter :: nk = 3
  real ( kind = 8 ) u0
  real ( kind = 8 ) u1
  real ( kind = 8 ), parameter :: umax =   30.0D+00
  real ( kind = 8 ), parameter :: umin =  150.0D+00
  real ( kind = 8 ) v0
  real ( kind = 8 ) v1
  real ( kind = 8 ), parameter :: vmax =    5.0D+00
  real ( kind = 8 ), parameter :: vmin =    1.0D+00
  real ( kind = 8 ) w0
  real ( kind = 8 ) w1
  real ( kind = 8 ), parameter :: wmax = 30.0D+00
  real ( kind = 8 ), parameter :: wmin = -30.0D+00

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TILER_3D'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A simple example of transfinite'
  write ( *, '(a)' ) '  interpolation in 3D.'
!
!  Write the first line of the output file, which is the number
!  of triangles in the output 3D shape.  (In our simple example,
!  each sub-box will contain a tetrahedron.)
!
  open ( unit = 1, file = 'tiler_3d.tri', status = 'replace' )

  write ( 1, '(i12)' ) 4 * ni * nj * nk
!
!  Consider items with index (I,*,*):
!
  do i = 1, ni

    u0 = ( real ( ni - i + 1, kind = 8 ) * umin    &
         + real (      i - 1, kind = 8 ) * umax )  &
         / real ( ni,         kind = 8 )

    u1 = ( real ( ni - i, kind = 8 ) * umin   &
         + real (      i, kind = 8 ) * umax ) &
         / real ( ni,     kind = 8 )
!
!  Consider items with index (I,J,*):
!
    do j = 1, nj

      v0 = ( real ( nj - j + 1, kind = 8 ) * vmin   &
           + real (      j - 1, kind = 8 ) * vmax ) &
           / real ( nj,         kind = 8 )

      v1 = ( real ( nj - j, kind = 8 ) * vmin   &
           + real (      j, kind = 8 ) * vmax ) &
           / real ( nj,     kind = 8 )
!
!  Consider items with index (I,J,K):
!
      do k = 1, nk

        w0 = ( real ( nk - k + 1, kind = 8 ) * wmin   &
             + real (      k - 1, kind = 8 ) * wmax ) &
             / real ( nk,         kind = 8 )

        w1 = ( real ( nk - k, kind = 8 ) * wmin   &
             + real (      k, kind = 8 ) * wmax ) &
             / real ( nk,     kind = 8 )
!
!  Fill sub-box (I,J,K) with the 3-D "tile".
!
        call sub_box_tiler_3d ( umin, vmin, wmin, umax, vmax, wmax, &
          u0, v0, w0, u1, v1, w1 )

      end do

    end do

  end do

  close ( unit = 1 )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TILER_3D:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, u, v, w, x, y, z )

!*****************************************************************************80
!
!! BOUNDARY_3D returns the (X,Y,Z) coordinates of a point (U,V,W).
!
!  Discussion:
!
!    In this example, a single formula describes the mapping
!
!      (U,V,W) => (X,Y,Z).
!
!    It is more common (and more interesting) for the formula
!    to depend on which face of the boundary is being considered.
!
!    This routine is only called for points (U,V,W) where one of the
!    values U, V and W is "extreme", so there are generally six cases
!    to consider, for general boundaries.  The coding has been set
!    up with this in mind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, pages 461-477, 1973.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) UMIN, VMIN, WMIN, UMAX, VMAX, WMAX, the (U,V,W) coordinates
!    of the two opposite corners of the big box.
!
!    Input, real ( kind = 8 ) U, V, W, the (U,V,W) coordinates of a point in the big box.
!
!    Output, real ( kind = 8 ) X, Y, Z, the (X,Y,Z) coordinates of the point.
!
  implicit none

  real ( kind = 8 ), parameter :: deg2rad = 3.14159265D+00 / 180.0D+00
  real ( kind = 8 ) psi
  real ( kind = 8 ) theta
  real ( kind = 8 ) u
  real ( kind = 8 ) umax
  real ( kind = 8 ) umin
  real ( kind = 8 ) v
  real ( kind = 8 ) vmax
  real ( kind = 8 ) vmin
  real ( kind = 8 ) w
  real ( kind = 8 ) wmax
  real ( kind = 8 ) wmin
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  theta = u * deg2rad
  psi = w * deg2rad

       if ( u == umin ) then
    x = v * cos ( theta ) * cos ( psi )
    y = v * sin ( theta ) * cos ( psi )
    z = v                 * sin ( psi )
  else if ( u == umax ) then
    x = v * cos ( theta ) * cos ( psi )
    y = v * sin ( theta ) * cos ( psi )
    z = v                 * sin ( psi )
  else if ( v == vmin ) then
    x = v * cos ( theta ) * cos ( psi )
    y = v * sin ( theta ) * cos ( psi )
    z = v                 * sin ( psi )
  else if ( v == vmax ) then
    x = v * cos ( theta ) * cos ( psi )
    y = v * sin ( theta ) * cos ( psi )
    z = v                 * sin ( psi )
  else if ( w == wmin ) then
    x = v * cos ( theta ) * cos ( psi )
    y = v * sin ( theta ) * cos ( psi )
    z = v                 * sin ( psi )
  else if ( w == wmax ) then
    x = v * cos ( theta ) * cos ( psi )
    y = v * sin ( theta ) * cos ( psi )
    z = v                 * sin ( psi )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BOUNDARY_3D - Fatal error!'
    write ( *, '(a)' ) '  The input point is not on the boundary.'
    stop
  end if

  return
end
subroutine sub_box_tiler_3d ( umin, vmin, wmin, umax, vmax, wmax, &
  u0, v0, w0, u1, v1, w1 )

!*****************************************************************************80
!
!! SUB_BOX_TILER_3D "tiles" a 3D sub-box with a given pattern.
!
!  Discussion:
!
!    This routine knows the (U,V,W) coordinates of the big box and the
!    sub box, and knows the shape of the object to be place in the sub-box.
!    It uses transfinite interpolation to put the shape in the box.
!    This requires that, for each point of the shape to be mapped, the
!    (X,Y,Z) coordinates be evaluated for points on the surface of
!    the big box, namely, at 8 corners, 12 edges, and 6 faces.  These
!    values are then blended to give a sensible (X,Y,Z) coordinate for
!    the point belonging to the shape.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, pages 461-477, 1973.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) UMIN, VMIN, WMIN, UMAX, VMAX, WMAX, the (U,V,W) coordinates
!    of the two opposite corners of the big box.
!
!    Input, real ( kind = 8 ) U0, V0, W0, U1, V1, W1, the (U,V,W) coordinates of the
!    two oppositie corners of the sub-box.
!
  implicit none

  integer ( kind = 4 ), parameter :: npoint = 4

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) r_tab(npoint)
  real ( kind = 8 ), parameter :: r0 = 0.0D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ) s
  real ( kind = 8 ) s_tab(npoint)
  real ( kind = 8 ), parameter :: s0 = 0.0D+00
  real ( kind = 8 ), parameter :: s1 = 1.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) t_tab(npoint)
  real ( kind = 8 ), parameter :: t0 = 0.0D+00
  real ( kind = 8 ), parameter :: t1 = 1.0D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) u0
  real ( kind = 8 ) u1
  real ( kind = 8 ) umax
  real ( kind = 8 ) umin
  real ( kind = 8 ) v
  real ( kind = 8 ) v0
  real ( kind = 8 ) v1
  real ( kind = 8 ) vmax
  real ( kind = 8 ) vmin
  real ( kind = 8 ) w
  real ( kind = 8 ) w0
  real ( kind = 8 ) w1
  real ( kind = 8 ) wmax
  real ( kind = 8 ) wmin
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) x000
  real ( kind = 8 ) x00t
  real ( kind = 8 ) x001
  real ( kind = 8 ) x010
  real ( kind = 8 ) x011
  real ( kind = 8 ) x01t
  real ( kind = 8 ) x0s0
  real ( kind = 8 ) x0s1
  real ( kind = 8 ) x0st
  real ( kind = 8 ) x100
  real ( kind = 8 ) x101
  real ( kind = 8 ) x10t
  real ( kind = 8 ) x110
  real ( kind = 8 ) x111
  real ( kind = 8 ) x11t
  real ( kind = 8 ) x1s0
  real ( kind = 8 ) x1s1
  real ( kind = 8 ) x1st
  real ( kind = 8 ) xr00
  real ( kind = 8 ) xr01
  real ( kind = 8 ) xr0t
  real ( kind = 8 ) xr10
  real ( kind = 8 ) xr11
  real ( kind = 8 ) xr1t
  real ( kind = 8 ) xrs0
  real ( kind = 8 ) xrs1
  real ( kind = 8 ) y(npoint)
  real ( kind = 8 ) y000
  real ( kind = 8 ) y00t
  real ( kind = 8 ) y001
  real ( kind = 8 ) y010
  real ( kind = 8 ) y011
  real ( kind = 8 ) y01t
  real ( kind = 8 ) y0s0
  real ( kind = 8 ) y0s1
  real ( kind = 8 ) y0st
  real ( kind = 8 ) y100
  real ( kind = 8 ) y101
  real ( kind = 8 ) y10t
  real ( kind = 8 ) y110
  real ( kind = 8 ) y111
  real ( kind = 8 ) y11t
  real ( kind = 8 ) y1s0
  real ( kind = 8 ) y1s1
  real ( kind = 8 ) y1st
  real ( kind = 8 ) yr00
  real ( kind = 8 ) yr01
  real ( kind = 8 ) yr0t
  real ( kind = 8 ) yr10
  real ( kind = 8 ) yr11
  real ( kind = 8 ) yr1t
  real ( kind = 8 ) yrs0
  real ( kind = 8 ) yrs1
  real ( kind = 8 ) z(npoint)
  real ( kind = 8 ) z000
  real ( kind = 8 ) z00t
  real ( kind = 8 ) z001
  real ( kind = 8 ) z010
  real ( kind = 8 ) z011
  real ( kind = 8 ) z01t
  real ( kind = 8 ) z0s0
  real ( kind = 8 ) z0s1
  real ( kind = 8 ) z0st
  real ( kind = 8 ) z100
  real ( kind = 8 ) z101
  real ( kind = 8 ) z10t
  real ( kind = 8 ) z110
  real ( kind = 8 ) z111
  real ( kind = 8 ) z11t
  real ( kind = 8 ) z1s0
  real ( kind = 8 ) z1s1
  real ( kind = 8 ) z1st
  real ( kind = 8 ) zr00
  real ( kind = 8 ) zr01
  real ( kind = 8 ) zr0t
  real ( kind = 8 ) zr10
  real ( kind = 8 ) zr11
  real ( kind = 8 ) zr1t
  real ( kind = 8 ) zrs0
  real ( kind = 8 ) zrs1
!
!  Here are the (R,S,T) coordinates of the tetrahedron that we can think
!  of as the "tile" we need to place in each sub-box.
!
!  The (R,S,T) coordinate system is assumed to range from 0 to 1.
!
  r_tab(1) =  0.2D+00
  s_tab(1) =  0.2D+00
  t_tab(1) =  0.2D+00

  r_tab(2) =  0.8D+00
  s_tab(2) =  0.2D+00
  t_tab(2) =  0.2D+00

  r_tab(3) =  0.5D+00
  s_tab(3) =  0.8D+00
  t_tab(3) =  0.2D+00

  r_tab(4) =  0.5D+00
  s_tab(4) =  0.5D+00
  t_tab(4) =  0.8D+00
!
!  Compute the (X,Y,Z) coordinates of the corners of the (U,V,W) box.
!  These really only need to be computed once ever.
!
  call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umin, vmin, wmin, &
    x000, y000, z000 )

  call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umax, vmin, wmin, &
    x100, y100, z100 )

  call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umin, vmax, wmin, &
    x010, y010, z010 )

  call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umax, vmax, wmin, &
    x110, y110, z110 )

  call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umin, vmin, wmax, &
    x001, y001, z001 )

  call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umax, vmin, wmax, &
    x101, y101, z101 )

  call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umin, vmax, wmax, &
    x011, y011, z011 )

  call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umax, vmax, wmax, &
    x111, y111, z111 )
!
!  Now figure out the (X,Y,Z) coordinates of the tile point with
!  given (R,S,T) coordinates.  This depends on the positions of
!  all sorts of points on the corners, edges, and faces of the big box.
!
  do i = 1, npoint
!
!  Get the (R,S,T) coordinates of point I.
!
    r = r_tab(i)
    s = s_tab(i)
    t = t_tab(i)
!
!  Get the corresponding point (U,V,W) in the rectangular space.
!
    u = ( ( r1 - r      ) * u0   &
        + (      r - r0 ) * u1 ) &
        / ( r1     - r0 )

    v = ( ( s1 - s      ) * v0   &
        + (      s - s0 ) * v1 ) &
        / ( s1     - s0 )

    w = ( ( t1 - t      ) * w0   &
        + (      t - t0 ) * w1 ) &
        / ( t1     - t0 )
!
!  Evaluate (X,Y,Z) on the 12 edges "near" (U,V,W).
!
    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umin, vmin, w, &
      x00t, y00t, z00t )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umin, vmax, w, &
      x01t, y01t, z01t )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umax, vmin, w, &
      x10t, y10t, z10t )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umax, vmax, w, &
      x11t, y11t, z11t )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umin, v, wmin, &
      x0s0, y0s0, z0s0 )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umin, v, wmax, &
      x0s1, y0s1, z0s1 )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umax, v, wmin, &
      x1s0, y1s0, z1s0 )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umax, v, wmax, &
      x1s1, y1s1, z1s1 )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, u, vmin, wmin, &
      xr00, yr00, zr00 )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, u, vmin, wmax, &
      xr01, yr01, zr01 )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, u, vmax, wmin, &
      xr10, yr10, zr10 )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, u, vmax, wmax, &
      xr11, yr11, zr11 )
!
!  Evaluate (X,Y,Z) on the six faces near (U,V,W).
!
    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umin, v, w, &
      x0st, y0st, z0st )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, umax, v, w, &
      x1st, y1st, z1st )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, u, vmin, w, &
      xr0t, yr0t, zr0t )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, u, vmax, w, &
      xr1t, yr1t, zr1t )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, u, v, wmin, &
      xrs0, yrs0, zrs0 )

    call boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, u, v, wmax, &
      xrs1, yrs1, zrs1 )
!
!  Now figure out the location of point
!    I => (R,S,T) => (U,V,W) => (X(I),Y(I),Z(I)).
!
    x(i) = ( ( umax - u    ) * ( vmax - v    ) * ( wmax - w    ) * x000 &
           - ( umax - u    ) * ( vmax - v    ) * ( wmax - wmin ) * x00t &
           + ( umax - u    ) * ( vmax - v    ) * ( w    - wmin ) * x001 &
           - ( umax - u    ) * ( vmax - vmin ) * ( wmax - w    ) * x0s0 &
           + ( umax - u    ) * ( vmax - vmin ) * ( wmax - wmin ) * x0st &
           - ( umax - u    ) * ( vmax - vmin ) * ( w    - wmin ) * x0s1 &
           + ( umax - u    ) * ( v    - vmin ) * ( wmax - w    ) * x010 &
           - ( umax - u    ) * ( v    - vmin ) * ( wmax - wmin ) * x01t &
           + ( umax - u    ) * ( v    - vmin ) * ( w    - wmin ) * x011 &
           - ( umax - umin ) * ( vmax - v    ) * ( wmax - w    ) * xr00 &
           + ( umax - umin ) * ( vmax - v    ) * ( wmax - wmin ) * xr0t &
           - ( umax - umin ) * ( vmax - v    ) * ( w    - wmin ) * xr01 &
           + ( umax - umin ) * ( vmax - vmin ) * ( wmax - w    ) * xrs0 &
           + ( umax - umin ) * ( vmax - vmin ) * ( w    - wmin ) * xrs1 &
           - ( umax - umin ) * ( v    - vmin ) * ( wmax - w    ) * xr10 &
           + ( umax - umin ) * ( v    - vmin ) * ( wmax - wmin ) * xr1t &
           - ( umax - umin ) * ( v    - vmin ) * ( w    - wmin ) * xr11 &
           + ( u    - umin ) * ( vmax - v    ) * ( wmax - w    ) * x100 &
           - ( u    - umin ) * ( vmax - v    ) * ( wmax - wmin ) * x10t &
           + ( u    - umin ) * ( vmax - v    ) * ( w    - wmin ) * x101 &
           - ( u    - umin ) * ( vmax - vmin ) * ( wmax - w    ) * x1s0 &
           + ( u    - umin ) * ( vmax - vmin ) * ( wmax - wmin ) * x1st &
           - ( u    - umin ) * ( vmax - vmin ) * ( w    - wmin ) * x1s1 &
           + ( u    - umin ) * ( v    - vmin ) * ( wmax - w    ) * x110 &
           - ( u    - umin ) * ( v    - vmin ) * ( wmax - wmin ) * x11t &
           + ( u    - umin ) * ( v    - vmin ) * ( w    - wmin ) * x111 ) &
         / ( ( umax - umin ) * ( vmax - vmin ) * ( wmax - wmin ) )

    y(i) = ( ( umax - u    ) * ( vmax - v    ) * ( wmax - w    ) * y000 &
           - ( umax - u    ) * ( vmax - v    ) * ( wmax - wmin ) * y00t &
           + ( umax - u    ) * ( vmax - v    ) * ( w    - wmin ) * y001 &
           - ( umax - u    ) * ( vmax - vmin ) * ( wmax - w    ) * y0s0 &
           + ( umax - u    ) * ( vmax - vmin ) * ( wmax - wmin ) * y0st &
           - ( umax - u    ) * ( vmax - vmin ) * ( w    - wmin ) * y0s1 &
           + ( umax - u    ) * ( v    - vmin ) * ( wmax - w    ) * y010 &
           - ( umax - u    ) * ( v    - vmin ) * ( wmax - wmin ) * y01t &
           + ( umax - u    ) * ( v    - vmin ) * ( w    - wmin ) * y011 &
           - ( umax - umin ) * ( vmax - v    ) * ( wmax - w    ) * yr00 &
           + ( umax - umin ) * ( vmax - v    ) * ( wmax - wmin ) * yr0t &
           - ( umax - umin ) * ( vmax - v    ) * ( w    - wmin ) * yr01 &
           + ( umax - umin ) * ( vmax - vmin ) * ( wmax - w    ) * yrs0 &
           + ( umax - umin ) * ( vmax - vmin ) * ( w    - wmin ) * yrs1 &
           - ( umax - umin ) * ( v    - vmin ) * ( wmax - w    ) * yr10 &
           + ( umax - umin ) * ( v    - vmin ) * ( wmax - wmin ) * yr1t &
           - ( umax - umin ) * ( v    - vmin ) * ( w    - wmin ) * yr11 &
           + ( u    - umin ) * ( vmax - v    ) * ( wmax - w    ) * y100 &
           - ( u    - umin ) * ( vmax - v    ) * ( wmax - wmin ) * y10t &
           + ( u    - umin ) * ( vmax - v    ) * ( w    - wmin ) * y101 &
           - ( u    - umin ) * ( vmax - vmin ) * ( wmax - w    ) * y1s0 &
           + ( u    - umin ) * ( vmax - vmin ) * ( wmax - wmin ) * y1st &
           - ( u    - umin ) * ( vmax - vmin ) * ( w    - wmin ) * y1s1 &
           + ( u    - umin ) * ( v    - vmin ) * ( wmax - w    ) * y110 &
           - ( u    - umin ) * ( v    - vmin ) * ( wmax - wmin ) * y11t &
           + ( u    - umin ) * ( v    - vmin ) * ( w    - wmin ) * y111 ) &
         / ( ( umax - umin ) * ( vmax - vmin ) * ( wmax - wmin ) )

    z(i) = ( ( umax - u    ) * ( vmax - v    ) * ( wmax - w    ) * z000 &
           - ( umax - u    ) * ( vmax - v    ) * ( wmax - wmin ) * z00t &
           + ( umax - u    ) * ( vmax - v    ) * ( w    - wmin ) * z001 &
           - ( umax - u    ) * ( vmax - vmin ) * ( wmax - w    ) * z0s0 &
           + ( umax - u    ) * ( vmax - vmin ) * ( wmax - wmin ) * z0st &
           - ( umax - u    ) * ( vmax - vmin ) * ( w    - wmin ) * z0s1 &
           + ( umax - u    ) * ( v    - vmin ) * ( wmax - w    ) * z010 &
           - ( umax - u    ) * ( v    - vmin ) * ( wmax - wmin ) * z01t &
           + ( umax - u    ) * ( v    - vmin ) * ( w    - wmin ) * z011 &
           - ( umax - umin ) * ( vmax - v    ) * ( wmax - w    ) * zr00 &
           + ( umax - umin ) * ( vmax - v    ) * ( wmax - wmin ) * zr0t &
           - ( umax - umin ) * ( vmax - v    ) * ( w    - wmin ) * zr01 &
           + ( umax - umin ) * ( vmax - vmin ) * ( wmax - w    ) * zrs0 &
           + ( umax - umin ) * ( vmax - vmin ) * ( w    - wmin ) * zrs1 &
           - ( umax - umin ) * ( v    - vmin ) * ( wmax - w    ) * zr10 &
           + ( umax - umin ) * ( v    - vmin ) * ( wmax - wmin ) * zr1t &
           - ( umax - umin ) * ( v    - vmin ) * ( w    - wmin ) * zr11 &
           + ( u    - umin ) * ( vmax - v    ) * ( wmax - w    ) * z100 &
           - ( u    - umin ) * ( vmax - v    ) * ( wmax - wmin ) * z10t &
           + ( u    - umin ) * ( vmax - v    ) * ( w    - wmin ) * z101 &
           - ( u    - umin ) * ( vmax - vmin ) * ( wmax - w    ) * z1s0 &
           + ( u    - umin ) * ( vmax - vmin ) * ( wmax - wmin ) * z1st &
           - ( u    - umin ) * ( vmax - vmin ) * ( w    - wmin ) * z1s1 &
           + ( u    - umin ) * ( v    - vmin ) * ( wmax - w    ) * z110 &
           - ( u    - umin ) * ( v    - vmin ) * ( wmax - wmin ) * z11t &
           + ( u    - umin ) * ( v    - vmin ) * ( w    - wmin ) * z111 ) &
         / ( ( umax - umin ) * ( vmax - vmin ) * ( wmax - wmin ) )

  end do

  write ( 1, '(6f12.6)' ) x(1), y(1), z(1), 0.0D+00, 0.0D+00, 0.0D+00
  write ( 1, '(6f12.6)' ) x(2), y(2), z(2), 0.0D+00, 0.0D+00, 0.0D+00
  write ( 1, '(6f12.6)' ) x(3), y(3), z(3), 0.0D+00, 0.0D+00, 0.0D+00

  write ( 1, '(6f12.6)' ) x(2), y(2), z(2), 0.0D+00, 0.0D+00, 0.0D+00
  write ( 1, '(6f12.6)' ) x(1), y(1), z(1), 0.0D+00, 0.0D+00, 0.0D+00
  write ( 1, '(6f12.6)' ) x(4), y(4), z(4), 0.0D+00, 0.0D+00, 0.0D+00

  write ( 1, '(6f12.6)' ) x(3), y(3), z(3), 0.0D+00, 0.0D+00, 0.0D+00
  write ( 1, '(6f12.6)' ) x(4), y(4), z(4), 0.0D+00, 0.0D+00, 0.0D+00
  write ( 1, '(6f12.6)' ) x(1), y(1), z(1), 0.0D+00, 0.0D+00, 0.0D+00

  write ( 1, '(6f12.6)' ) x(4), y(4), z(4), 0.0D+00, 0.0D+00, 0.0D+00
  write ( 1, '(6f12.6)' ) x(3), y(3), z(3), 0.0D+00, 0.0D+00, 0.0D+00
  write ( 1, '(6f12.6)' ) x(2), y(2), z(2), 0.0D+00, 0.0D+00, 0.0D+00

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
