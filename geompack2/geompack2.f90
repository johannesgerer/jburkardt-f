function angle ( xa, ya, xb, yb, xc, yc )

!*****************************************************************************80
!
!! ANGLE computes the interior angle at a vertex defined by 3 points.
!
!  Discussion:
!
!    ANGLE computes the interior angle, in radians, at vertex
!    (XB,YB) of the chain formed by the directed edges from
!    (XA,YA) to (XB,YB) to (XC,YC).  The interior is to the
!    left of the two directed edges.
!
!  Modified:
!
!    17 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XA, YA, XB, YB, XC, YC, the coordinates of the 
!    vertices.
!
!    Output, real ( kind = 8 ) ANGLE, the interior angle formed by
!    the vertex, in radians, between 0 and 2*PI.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) ya
  real ( kind = 8 ) yb
  real ( kind = 8 ) yc

  x1 = xa - xb
  y1 = ya - yb
  x2 = xc - xb
  y2 = yc - yb

  t = sqrt ( ( x1 * x1 + y1 * y1 ) * ( x2 * x2 + y2 * y2 ) )

  if ( t == 0.0D+00 ) then
    angle = pi
    return
  end if

  t = ( x1 * x2 + y1 * y2 ) / t

  if ( t < -1.0D+00 ) then
    t = -1.0D+00
  else if ( 1.0D+00 < t ) then
    t = 1.0D+00
  end if

  angle = acos ( t )

  if ( x2 * y1 - y2 * x1 < 0.0D+00 ) then
    angle = 2.0D+00 * pi - angle
  end if

  return
end
function areapg ( nvrt, xc, yc )

!*****************************************************************************80
!
!! AREAPG computes twice the signed area of a simple polygon.
!
!  Modified:
!
!    13 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of 
!    the polygon.  N must be at least 3.
!
!    Input, real ( kind = 8 ) XC(NVRT), YC(NVRT), the X and Y coordinates 
!    of the vertices.
!
!    Output, real ( kind = 8 ) AREAPG, twice the signed area of the polygon,
!    which will be positive if the vertices were listed in counter clockwise
!    order, and negative otherwise.
!
  implicit none

  integer ( kind = 4 ) nvrt

  real ( kind = 8 ) areapg
  integer ( kind = 4 ) i
  real ( kind = 8 ) sum2
  real ( kind = 8 ) xc(nvrt)
  real ( kind = 8 ) yc(nvrt)

  sum2 = xc(1) * ( yc(2) - yc(nvrt) )

  do i = 2, nvrt-1
    sum2 = sum2 + xc(i) * ( yc(i+1) - yc(i-1) )
  end do

  sum2 = sum2 + xc(nvrt) * ( yc(1) - yc(nvrt-1) )

  areapg = sum2

  return
end
function areatr ( xa, ya, xb, yb, xc, yc )

!*****************************************************************************80
!
!! AREATR computes twice the signed area of a triangle.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XA, YA, XB, YB, XC, YC, the coordinates of the 
!    vertices.
!
!    Output, real ( kind = 8 ) AREATR, twice the signed area of the triangle.
!    This will be positive if the vertices are listed in counter clockwise 
!    order.
!
  implicit none

  real ( kind = 8 ) areatr 
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) ya
  real ( kind = 8 ) yb
  real ( kind = 8 ) yc

  areatr = ( xb - xa ) * ( yc - ya ) - ( xc - xa ) * ( yb - ya )

  return
end
subroutine bedgmv ( nvc, npolg, nvert, maxvc, h, vcl, hvl, pvl, vstart, vnum, &
  ierror )

!*****************************************************************************80
!
!! BEDGMV generates boundary edge mesh vertices.
!
!  Purpose: 
!
!    Generate mesh vertices on boundary of convex polygons
!    of decomposition with spacing determined by H array.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) NVC, the number of coordinates or positions used 
!    in VCL array.
!
!    Input, integer ( kind = 4 ) NPOLG, the number of polygons or positions used 
!    in HVL array.
!
!    Input, integer ( kind = 4 ) NVERT, the number of vertices or positions used in 
!    PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, real ( kind = 8 ) H(1:NPOLG), the spacing of mesh vertices for 
!    convex polygons.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:NPOLG, the head vertex list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:NVERT), the polygon vertex list.
!
!    Output, integer ( kind = 4 ) VSTART(1:NVERT), the start location in VCL for mesh 
!    vertices on each edge in PVL if there are any, else 0.
!
!    Output, integer ( kind = 4 ) VNUM(1:NVERT), the number of mesh vertices on interior 
!    of each edge in PVL; entry is negated if mesh vertices are listed in 
!    backward order in VCL.
!
!    Output, integer ( kind = 4 ) IERROR, is set to 3 on error.
!
  implicit none

  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvert

  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ), parameter :: edgv = 4
  real ( kind = 8 ) h(npolg)
  real ( kind = 8 ) hh
  integer ( kind = 4 ) hvl(npolg)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) leng
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,nvert)
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v
  integer ( kind = 4 ) vstart(nvert)
  integer ( kind = 4 ) vnum(nvert)
  real ( kind = 8 ) vcl(2,maxvc)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  ierror = 0
  vstart(1:nvert) = -1

  do k = 1, npolg

    i = hvl(k)

    do

      j = pvl(succ,i)

      if ( vstart(i) == -1 ) then

        u = pvl(loc,i)
        v = pvl(loc,j)
        x = vcl(1,u)
        y = vcl(2,u)
        leng = sqrt ( ( vcl(1,v) - x )**2 + ( vcl(2,v) - y )**2 )
        ia = pvl(edgv,i)

        if ( ia <= 0 ) then
          hh = h(k)
        else
          hh = sqrt ( h(k) * h(pvl(polg,ia)) )
        end if

        if ( hh == 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'BEDGMV - Fatal error!'
          write ( *, '(a)' ) '  HH = 0.'
          stop
        end if

        l = int ( leng / hh )

        if ( real ( l, kind = 8 ) / real ( 2 * l + 1, kind = 8 ) &
          < ( leng / hh ) - real ( l, kind = 8 ) ) then
          l = l + 1
        end if

        if ( l <= 1 ) then

          vstart(i) = 0
          vnum(i) = 0

        else

          dx = ( vcl(1,v) - x ) / real ( l, kind = 8 )
          dy = ( vcl(2,v) - y ) / real ( l, kind = 8 )
          l = l - 1

          if ( maxvc < nvc + l ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'BEDGMV - Fatal error!'
            write ( *, '(a)' ) '  MAXVC < NVC + L.'
            ierror = 3
            return
          end if

          vstart(i) = nvc + 1
          vnum(i) = l

          do m = 1, l
            x = x + dx
            y = y + dy
            nvc = nvc + 1
            vcl(1,nvc) = x
            vcl(2,nvc) = y
          end do

        end if

        if ( 0 < ia ) then
          vstart(ia) = vstart(i)
          vnum(ia) = -vnum(i)
        end if

      end if

      i = j

      if ( i == hvl(k) ) then
        exit
      end if

    end do

  end do

  return
end
subroutine bnsrt2 ( binexp, n, a, map, bin, iwk )

!*****************************************************************************80
!
!! BNSRT2 bin sorts N points in 2D into increasing bin order.
!
!  Purpose: 
!
!    Use a bin sort to obtain the permutation of N 2D
!    double precision points so that points are in increasing bin
!    order, where the N points are assigned to about N**BINEXP bins.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BINEXP, the exponent for the number of bins.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) A(2,*), the points to be binned.
!
!    Input/output, integer ( kind = 4 ) MAP(N); on input, the points of A with indices 
!    MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output, MAP has
!    been permuted so bin of MAP(1) <= bin of MAP(2) <= ... <= bin of MAP(N).
!
!    Workspace, integer BIN(N), used for bin numbers and permutation of 1 to N.
!
!    Workspace, integer IWK(N), used for copy of MAP array.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,*)
  integer ( kind = 4 ) bin(n)
  real ( kind = 8 ) binexp
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwk(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) nside
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  nside = int ( real ( n, kind = 8 )**( binexp / 2.0D+00 ) + 0.5D+00 )

  if ( nside <= 1 ) then
    return
  end if

  xmin = a(1,map(1))
  ymin = a(2,map(1))
  xmax = xmin
  ymax = ymin
  do i = 2, n
    j = map(i)
    xmin = min ( xmin, a(1,j) )
    xmax = max ( xmax, a(1,j) )
    ymin = min ( ymin, a(2,j) )
    ymax = max ( ymax, a(2,j) )
  end do

  dx = 1.0001D+00 * ( xmax - xmin ) / real ( nside, kind = 8 )
  dy = 1.0001D+00 * ( ymax - ymin ) / real ( nside, kind = 8 )

  if ( dx == 0.0D+00 ) then
    dx = 1.0D+00
  end if

  if ( dy == 0.0D+00 ) then
    dy = 1.0D+00
  end if

  do i = 1, n
    j = map(i)
    iwk(i) = j
    map(i) = i
    k = int ( ( a(1,j) - xmin ) / dx )
    l = int ( ( a(2,j) - ymin ) / dy )
    if ( mod ( k, 2 ) == 0 ) then
      bin(i) = k * nside + l
    else
      bin(i) = ( k + 1 ) * nside - l - 1
    end if
  end do

  call ihpsrt ( 1, n, 1, bin, map )

  bin(1:n) = map(1:n)

  do i = 1, n
    map(i) = iwk(bin(i))
  end do

  return
end
function cmcirc ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! CMCIRC determines whether a point lies within a circle through 3 points.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the point to
!    be tested.
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates of
!    three points that define a circle.
!
!    Output, integer ( kind = 4 ) CMCIRC, reports the test results:
!     2, if the three vertices are collinear,
!     1, if (X0,Y0) is inside the circle,
!     0, if (X0,Y0) is on the circle,
!    -1, if (X0,Y0) is outside the circle.
!
  real ( kind = 8 ) a11
  real ( kind = 8 ) a12
  real ( kind = 8 ) a21
  real ( kind = 8 ) a22
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  integer ( kind = 4 ) cmcirc
  real ( kind = 8 ) det
  real ( kind = 8 ) diff
  real ( kind = 8 ) rsq
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs
  real ( kind = 8 ) xc
  real ( kind = 8 ) yc
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  tol = 100.0D+00 * epsilon ( tol )
  cmcirc = 2
  a11 = x2 - x1
  a12 = y2 - y1
  a21 = x3 - x1
  a22 = y3 - y1
  tolabs = tol * max ( abs ( a11), abs ( a12), abs ( a21), abs ( a22) )
  det = a11 * a22 - a21 * a12

  if ( abs ( det ) <= tolabs ) then
    return
  end if

  b1 = a11 * a11 + a12 * a12
  b2 = a21 * a21 + a22 * a22
  det = 2.0D+00 * det
  xc = ( b1 * a22 - b2 * a12 ) / det
  yc = ( b2 * a11 - b1 * a21 ) / det
  rsq = xc * xc + yc * yc
  diff = ( ( x0 - x1 - xc)**2 + ( y0 - y1 - yc )**2 ) - rsq
  tolabs = tol * rsq

  if ( diff < - tolabs ) then
    cmcirc = 1
  else if ( tolabs < diff ) then
    cmcirc = -1
  else
    cmcirc = 0
  end if

  return
end
subroutine cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, &
  maxpv, maxiw, maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

!*****************************************************************************80
!
!! CVDEC2 decomposes a polygonal region into convex polygons.
!
!  Purpose: 
!
!    Decompose general polygonal region (which is decomposed
!    into simple polygons on input) into convex polygons using
!    vertex coordinate list, head vertex list, and polygon vertex
!    list data structures.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in radians 
!    used in controlling vertices to be considered as an endpoint of a
!    separator.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter in radians 
!    used in accepting separator(s).
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or positions 
!    used in VCL.
!
!    Input/output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions or 
!    positions used in HVL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of polygon vertices or positions 
!    used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array, 
!    should be greater than or equal to the number of vertex coordinates
!    required for decomposition.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM arrays, 
!    should be greater than or equal to the number of polygons required for
!    decomposition.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays; 
!    should be greater than or equal to the number of polygon vertices 
!    required for decomposition.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should be 
!    about 3 times maximum number of vertices in any polygon.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be 
!    about 5 times maximum number of vertices in any polygon.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), region numbers.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT),
!    the polygon vertex list and interior angles; see routine DSPGDC for 
!    more details.  Note that the data structures should be as output from 
!    routine SPDEC2.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  For abnormal return,
!    IERROR is set to 3, 4, 5, 6, 7, 206, 207, 208, 209, 210, or 212.
!
  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  integer ( kind = 4 ) hvl(maxhv)
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) piptol
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  real ( kind = 8 ) tol
  integer ( kind = 4 ) v
  real ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) w1
  integer ( kind = 4 ) w2
  real ( kind = 8 ) wk(maxwk)

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  For each reflex vertex, resolve it with one or two separators
!  and update VCL, HVL, PVL, IANG.
!
  piptol = pi + tol
  v = 1

  do

    if ( nvert < v ) then
      exit
    end if

    if ( piptol < iang(v) ) then

      call resvrt ( v, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
        maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      call insed2 ( v ,w1, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
        pvl, iang, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      if ( 0 < w2 ) then
        call insed2 ( v, w2, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
          pvl, iang, ierror )
      end if

      if ( ierror /= 0 ) then
        return
      end if

    end if

    v = v + 1

  end do

  return
end
subroutine cvdtri ( inter, ldv, nt, vcl, til, tedg, sptr, ierror )

!*****************************************************************************80
!
!! CVDTRI converts boundary triangles to Delaunay triangles.
!
!  Purpose: 
!
!    Convert triangles in strip near boundary of polygon
!    or inside polygon to Delaunay triangles.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical INTER, is .TRUE. if and only if there is at least 
!    one interior mesh vertex.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input, integer ( kind = 4 ) NT, the number of triangles in strip or polygon.
!
!    Input, VCL(1:2,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:NT), the triangle incidence list.
!
!    Input/output, integer ( kind = 4 ) TEDG(1:3,1:NT) - TEDG(J,I) refers to edge with 
!    vertices TIL(J:J+1,I) and contains index of merge edge or
!    a value greater than NT for edge of chains.
!
!    Workspace, SPTR(1:NT) - SPTR(I) = -1 if merge edge I is not in LOP stack,
!    else greater than or equal to 0 and pointer (index of SPTR) to next 
!    edge in stack (0 indicates bottom of stack).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return:
!    IERROR is set to 231.
!
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) nt

  integer ( kind = 4 ) e
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind(2)
  logical inter
  integer ( kind = 4 ) itr(2)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mxtr
  logical sflag
  integer ( kind = 4 ) sptr(nt)
  integer ( kind = 4 ) tedg(3,nt)
  integer ( kind = 4 ) til(3,nt)
  integer ( kind = 4 ) top
  real ( kind = 8 ) vcl(ldv,*)

  ierror = 0
  sflag = .true.
  sptr(1:nt) = -1

  do k = 1, nt

    mxtr = k + 1

    if ( k == nt ) then
      if ( .not. inter ) then
        return
      end if
      mxtr = nt
      sflag = .false.
    end if

    top = k
    sptr(k) = 0

    do

      e = top
      top = sptr(e)

      call fndtri ( e, mxtr, sflag, tedg, itr, ind, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      call lop ( itr, ind, k, top, ldv, vcl, til, tedg, sptr )

      if ( top <= 0 ) then
        exit
      end if

    end do

  end do

  return
end
function degrees_to_radians ( angle )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, an angle in degrees.
!
!    Output, real ( kind = 8 ) DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  real ( kind = 8 ) angle
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  degrees_to_radians = ( angle / 180.0D+00 ) * pi

  return
end
subroutine delaunay_print ( num_pts, xc, num_tri, nodtri, tnbr )

!*****************************************************************************80
!
!! DELAUNAY_PRINT prints out information defining a Delaunay triangulation.
!
!  Modified:
!
!    08 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUM_PTS, the number of points.
!
!    Input, real ( kind = 8 ) XC(2,NUM_PTS), the point coordinates.
!
!    Input, integer ( kind = 4 ) NUM_TRI, the number of triangles.
!
!    Input, integer ( kind = 4 ) NODTRI(3,NUM_TRI), the nodes that make up the triangles.
!
!    Input, integer ( kind = 4 ) TNBR(3,NUM_TRI), the triangle neighbors on each side.
!
  integer ( kind = 4 ) num_pts
  integer ( kind = 4 ) num_tri

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nodtri(3,num_tri)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t
  integer ( kind = 4 ) tnbr(3,num_tri)
  real ( kind = 8 ) xc(2,num_pts)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DELAUNAY_PRINT'
  write ( *, '(a)' ) '  Information defining a Delaunay triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points is ', num_pts

  call r8mat_print ( num_pts, num_pts, 2, transpose ( xc ), &
    '  Point coordinates (transpose of internal array)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of triangles is ', num_tri
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sets of three points are used as vertices of'
  write ( *, '(a)' ) '  the triangles.  For each triangle, the points'
  write ( *, '(a)' ) '  are listed in counterclockwise order.'

  call i4mat_print ( num_tri, num_tri, 3, transpose ( nodtri ), &
    '  Nodes that make up triangles (transpose of internal array)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  On each side of a given triangle, there is either'
  write ( *, '(a)' ) '  another triangle, or a piece of the convex hull.'
  write ( *, '(a)' ) '  For each triangle, we list the indices of the three'
  write ( *, '(a)' ) '  neighbors, or (if negative) the codes of the'
  write ( *, '(a)' ) '  segments of the convex hull.'

  call i4mat_print ( num_tri, num_tri, 3, transpose ( tnbr ), &
    '  Indices of neighboring triangles (transpose of internal array)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of boundary points (and segments) is ', &
    2 * num_pts - num_tri - 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The segments that make up the convex hull can be'
  write ( *, '(a)' ) '  determined from the negative entries of the triangle'
  write ( *, '(a)' ) '  neighbor list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  # Tri Side  N1  N2'
  write ( *, '(a)' ) ' '
  k = 0

  do i = 1, num_tri
    do j = 1, 3
      if ( tnbr(j,i) < 0 ) then
        s = - tnbr(j,i)
        t = s / 3
        s = mod ( s, 3 ) + 1
        k = k + 1
        n1 = nodtri(s,t)
        n2 = nodtri(i4_wrap(s+1,1,3),t)
        write ( *, '(5i4)' ) k, t, s, n1, n2
      end if
    end do
  end do

  return
end
subroutine dhpsrt ( k, n, lda, a, map )

!*****************************************************************************80
!
!! DHPSRT sorts points into lexicographic order using heap sort
!
!  Discussion:
!
!    This routine uses heapsort to obtain the permutation of N K-dimensional
!    points so that the points are in lexicographic increasing order.
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Modified:
!
!    19 February 2001
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the dimension of the points (for instance, 2
!    for points in the plane).
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in the calling
!    routine; LDA should be at least K.
!
!    Input, real ( kind = 8 ) A(LDA,N); A(I,J) contains the I-th coordinate
!    of point J.
!
!    Input/output, integer ( kind = 4 ) MAP(N).
!    On input, the points of A with indices MAP(1), MAP(2), ...,
!    MAP(N) are to be sorted.
!
!    On output, MAP contains a permutation of its input values, with the
!    property that, lexicographically,
!      A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(n)

  do i = n / 2, 1, -1
    call dsftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    call i4_swap ( map(1), map(i) )
    call dsftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! DIAEDG triangulates 4 points using the circumcircle criterion.
!
!  Diagram:
!
!    3---2      3---2
!    |\  |      |  /|
!    | \ |  or  | / |
!    |  \|      |/  |
!    0---1      0---1
!
!  Discussion:
!
!    Given four points, to be organized into two triangles, the routine
!    chooses 0--2 or 1--3 as the diagonal edge, based on the circumcircle
!    criterion.  
!
!    The points may be regarded as the vertices of a simple quadrilateral,
!    and should be listed in counterclockwise order.
!
!  Modified:
!
!    13 April 2004
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the points.
!
!    Output, integer ( kind = 4 ) DIAEDG, chooses a diagonal:
!    +1, choose edge 0--2;
!    -1, choose edge 1--3;
!     0, the four vertices are essentially cocircular.
!
  implicit none

  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  integer ( kind = 4 ) diaedg
  real ( kind = 8 ) dx10
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx30
  real ( kind = 8 ) dx32
  real ( kind = 8 ) dy10
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy30
  real ( kind = 8 ) dy32
  real ( kind = 8 ) s
  real ( kind = 8 ) tol
  real ( kind = 8 ) tola
  real ( kind = 8 ) tolb
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  tol = 100.0D+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32
!
!  If both angles 301 and 123 are acute, choose 1-3 as the edge.
!
  if ( tola < ca .and. tolb < cb ) then

    diaedg = -1
!
!  If both angles 301 and 123 are obtuse, choose 0-2 as the edge.
!
  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1
 
  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( tola < s ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end
subroutine diam2 ( nvrt, xc, yc, i1, i2, diamsq, ierror )

!*****************************************************************************80
!
!! DIAM2 finds the diameter of a convex polygon.
!
!  Purpose: 
!
!    Find the diameter of a convex polygon with vertices
!    given in counter clockwise order and with all interior angles < PI.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices.
!
!    Input, real ( kind = 8 ) XC(NVRT),YC(NVRT), the vertex coordinates in 
!    counter-clockwise order.
!
!    Output, integer ( kind = 4 ) I1, I2 , indices in XC, YC of the diameter edge; the
!    diameter is from (XC(I1),YC(I1)) to (XC(I2),YC(I2)).
!
!    Output, real ( kind = 8 ) DIAMSQ, the square of the diameter.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error was detected.
!    200, an error was detected.
!
  integer ( kind = 4 ) nvrt

  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) areatr
  real ( kind = 8 ) c1mtol
  real ( kind = 8 ) c1ptol
  real ( kind = 8 ) diamsq
  real ( kind = 8 ) dist
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) m
  real ( kind = 8 ) tol
  real ( kind = 8 ) xc(nvrt)
  real ( kind = 8 ) yc(nvrt)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Find the first vertex which is farthest from the edge connecting
!  vertices with indices NVRT, 1.
!
  c1mtol = 1.0D+00 - tol
  c1ptol = 1.0D+00 + tol
  j = nvrt
  jp1 = 1
  k = 2
  area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

  do

    area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k+1), yc(k+1) )

    if ( area2 <= area1 * c1ptol ) then
      exit
    end if

    area1 = area2
    k = k + 1

  end do

  m = k
  diamsq = 0.0D+00
!
!  Find diameter = maximum distance of antipodal pairs.
!
  area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

  do

    kp1 = k + 1
    if ( nvrt < kp1 ) then
      kp1 = 1
    end if

    area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(kp1), yc(kp1) )

    if ( area1 * c1ptol < area2 ) then
      k = k + 1
      area1 = area2
    else if ( area2 < area1 * c1mtol ) then
      j = jp1
      jp1 = j + 1
      area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )
    else
      k = k + 1
      j = jp1
      jp1 = j + 1
      area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )
    end if

    if ( m < j .or. nvrt < k ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIAM2 - Fatal error!'
      write ( *, '(a)' ) '  M < J or NVRT < K.'
      ierror = 200
      return
    end if

    dist = ( xc(j) - xc(k) )**2 + ( yc(j) - yc(k) )**2

    if ( diamsq < dist ) then
      diamsq = dist
      i1 = j
      i2 = k
    end if

    if ( j == m .and. k == nvrt ) then
      exit
    end if

  end do

  return
end
function dless ( k, p, q )

!*****************************************************************************80
!
!! DLESS determine whether P is lexicographically less than Q.
!
!  Discussion:
!
!    P and Q are K-dimensional points.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the spatial dimension of the points.
!
!    Input, real ( kind = 8 ) P(K), Q(K), the points to be compared.
!
!    Output, logical RLESS, is TRUE if P < Q, FALSE otherwise.
!
  integer ( kind = 4 ) k

  real ( kind = 8 ) cmax
  logical dless
  integer ( kind = 4 ) i
  real ( kind = 8 ) p(k)
  real ( kind = 8 ) q(k)
  real ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )

  do i = 1, k

    cmax = max ( abs ( p(i) ), abs ( q(i) ) )

    if ( tol * cmax < abs ( p(i) - q(i) ) .and. tol < cmax ) then

      if ( p(i) < q(i) ) then
        dless = .true.
      else
        dless = .false.
      end if

      return
    end if

  end do

  dless = .false.

  return
end
subroutine dsftdw ( l, u, k, lda, a, map )

!*****************************************************************************80
!
!! DSFTDW sifts A(*,MAP(L)) down a heap of size U.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, U, the lower and upper indices of part of the heap.
!
!    Input, integer ( kind = 4 ) K, the spatial dimension of the points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A in the calling routine.
!
!    Input, real ( kind = 8 ) A(LDA,N); A(I,J) contains the I-th coordinate
!    of point J.
!
!    Input/output, integer ( kind = 4 ) MAP(N).
!    On input, the points of A with indices MAP(1), MAP(2), ...,
!    MAP(N) are to be sorted.
!
!    On output, MAP contains a permutation of its input values, with the
!    property that, lexicographically,
!      A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
  integer ( kind = 4 ) lda

  real ( kind = 8 ) a(lda,*)
  logical dless
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) map(*)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) u

  i = l
  j = 2 * i
  t = map(i)

  do while ( j <= u )

    if ( j < u ) then
      if ( dless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( dless ( k, a(1,map(j)), a(1,t) ) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end
subroutine dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxho, nvc, npolg, &
  nvert, nhola, regnum, hvl, pvl, iang, holv, ierror )

!*****************************************************************************80
!
!! DSMCPR initializes the polygonal decomposition data structure.
!
!  Purpose: 
!
!    Initialize the polygonal decomposition data structure
!    given a multiply-connected polygonal region with 1 outer
!    boundary curve and 0 or more inner boundary curves of holes.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NHOLE, the number of holes in the region.
!
!    Input, integer ( kind = 4 ) NVBC(1:NHOLE+1), the number of vertices per boundary curve; 
!    first boundary curve is the outer boundary of the region.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinates of 
!    boundary curves in counter clockwise order; 
!    NVC = NVBC(1) + ... + NVBC(NHOLE+1); 
!    positions 1 to NVBC(1) of VCL contain the vertex coordinates of the
!    outer boundary in counter clockwise order; positions NVBC(1)+1 to
!    NVBC(1)+NVBC(2) contain the vertex coordinates of the
!    first hole boundary in counter clockwise order, etc.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM arrays, 
!    should be greater than or equal to NHOLE + 1.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays;
!    should be greater than or equal to NVC.
!
!    Input, integer ( kind = 4 ) MAXHO, the maximum size available for HOLV array; should be
!    greater than or equal to NHOLE*2.
!
!    Output, integer ( kind = 4 ) NVC, the number of vertex coordinates, set to sum of 
!    NVBC(I).
!
!    Output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions, set to 1.
!    For consistency with DSPGDC.
!
!    Output, integer ( kind = 4 ) NVERT, the number of vertices in PVL, set to NVC.
!    For consistency with DSPGDC.
!
!    Output, integer ( kind = 4 ) NHOLA, number of attached holes, set to 0.
!    For consistency with DSPGDC.
!
!    Output, integer ( kind = 4 ) REGNUM(1:1), region number of only subregion, set to 1
!    For consistency with DSPGDC.
!
!    Output, integer ( kind = 4 ) HVL(1:NHOLE+1), the head vertex list; first entry is the 
!    head vertex (index in PVL) of outer boundary curve; next
!    NHOLE entries contain the head vertex of a hole.
!
!    Output, integer ( kind = 4 ) PVL(1:4,1:NVC), IANG(1:NVC), the polygon vertex list and 
!    interior angles; vertices of outer boundary curve are in counter clockwise 
!    order followed by vertices of each hole in CW hole; vertices
!    of each polygon are in a circular linked list; see
!    routine DSPGDC for more details of this data structure.
!
!    Output, integer ( kind = 4 ) HOLV(1:NHOLE*2), the indices in PVL of top and bottom 
!    vertices of holes; first (last) NHOLE entries are for top (bottom)
!    vertices; top (bottom) vertices are sorted in decreasing
!    (increasing) lexicographic (y,x) order of coordinates.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    For abnormal return, IERROR is set to 2, 4, or 5.
!
  integer ( kind = 4 ) maxho
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) nhole

  real ( kind = 8 ) angle
  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) ivs
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lvp
  integer ( kind = 4 ) lvs
  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) nvs
  integer ( kind = 4 ) hvl(nhole+1)
  integer ( kind = 4 ) holv(maxho)
  integer ( kind = 4 ) nvbc(nhole+1)
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(1)
  integer ( kind = 4 ), parameter :: succ = 3
  real ( kind = 8 ) vcl(2,*)

  ierror = 0
  nvc = sum ( nvbc(1:nhole+1) )
  npolg = 1
  nvert = nvc
  nhola = 0
  regnum(1) = 1

  if ( maxhv < nhole + 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DSMCPR - Fatal error!'
    write ( *, '(a)' ) '  MAXHV < NHOLE + 1.'
    ierror = 4
    return
  end if

  if ( maxpv < nvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DSMCPR - Fatal error!'
    write ( *, '(a)' ) '  MAXPV < NVC.'
    ierror = 5
    return
  end if

  if ( maxho < 2 * nhole ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DSMCPR - Fatal error!'
    write ( *, '(a)' ) '  MAXHO < 2 * NHOLE.'
    ierror = 2
    return
  end if
!
!  Initialize the HVL and PVL arrays.
!
  hvl(1) = 1
  nv = nvbc(1)

  do i = 1, nv
    pvl(loc,i) = i
    pvl(polg,i) = 1
    pvl(succ,i) = i + 1
    pvl(edgv,i) = 0
  end do

  pvl(succ,nv) = 1

  do j = 1, nhole
    hvl(j+1) = nv + 1
    nvs = nv + nvbc(j+1)
    do i = nv+1, nvs
      pvl(loc,i) = i
      pvl(polg,i) = 1
      pvl(succ,i) = i - 1
      pvl(edgv,i) = 0
    end do
    pvl(succ,nv+1) = nvs
    nv = nvs
  end do
!
!  Initialize the IANG array.
!
  do i = 1, nhole+1

    j = hvl(i)
    lvp = pvl(loc,j)
    iv = pvl(succ,j)
    lv = pvl(loc,iv)

    do

      ivs = pvl(succ,iv)
      lvs = pvl(loc,ivs)
      iang(iv) = angle ( vcl(1,lvp), vcl(2,lvp), vcl(1,lv), vcl(2,lv), &
        vcl(1,lvs), vcl(2,lvs) )

      if ( iv == j ) then
        exit
      end if

      lvp = lv
      iv = ivs
      lv = lvs

    end do

  end do
!
!  Initialize the HOLV array.
!
  if ( 0 < nhole ) then
    call holvrt ( nhole, vcl, hvl(2), pvl, holv )
  end if

  return
end
subroutine dsmdf2 ( hflag, nvc, npolg, maxwk, vcl, hvl, pvl, iang, ivrt, &
  xivrt, widsq, edgval, vrtval, area, wk, ierror )

!*****************************************************************************80
!
!! DSMDF2 sets up a data structure for a heuristic mesh distribution.
!
!  Purpose: 
!
!    Set up the data structure for heuristic mesh distribution
!    function from data structure for convex polygon decomposition
!    if HFLAG is .TRUE., else set up only IVRT and XIVRT.
!
!    Also compute areas of convex polygons.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical HFLAG, set to .TRUE. if data structure is to be constructed,
!    .FALSE. if only IVRT, XIVRT, AREA are to be computed.
!
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates in VCL array.
!
!    Input, integer ( kind = 4 ) NPOLG, the number of polygonal subregions in HVL array.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    2 times maximum number of vertices in any polygon.
!
!    Input, VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), real ( kind = 8 ) IANG(1:*), the polygon
!    vertex list, interior angles.
!
!    Output, integer ( kind = 4 ) IVRT(1:*), the indices of polygon vertices in VCL, ordered 
!    by polygon; same size as PVL.  For heuristic MDF data structure.
!
!    Output, XIVRT(1:NPOLG+1), the pointer to first vertex of each polygon
!    in IVRT; vertices of polygon K are IVRT(I) for I from
!    XIVRT(K) to XIVRT(K+1)-1.  For heuristic MDF data structure.
!
!    Output, real ( kind = 8 ) WIDSQ(1:NPOLG), the square of width of convex 
!    polygons.  For heuristic MDF data structure.
!
!    Output, real ( kind = 8 ) EDGVAL(1:*), the value associated with each 
!    edge of decomposition; same size as PVL.  For heuristic MDF data structure.
!
!    Output, real ( kind = 8 ) VRTVAL(1:NVC), the value associated with each 
!    vertex of decomposition.  For heuristic MDF data structure.
!
!    Output, real ( kind = 8 ) AREA(1:NPOLG), the area of convex polygons.
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 7 or 201.
!
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvc

  real ( kind = 8 ) area(npolg)
  real ( kind = 8 ) areapg
  integer ( kind = 4 ), parameter :: edgv = 4
  real ( kind = 8 ) edgval(*)
  logical hflag
  integer ( kind = 4 ) hvl(npolg)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(*)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jl
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nvrt
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pimtol
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,*)
  real ( kind = 8 ) s
  integer ( kind = 4 ), parameter :: succ = 3
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(2,nvc)
  real ( kind = 8 ) vrtval(nvc)
  real ( kind = 8 ) widsq(npolg)
  real ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xc
  integer ( kind = 4 ) xivrt(npolg+1)
  integer ( kind = 4 ) yc

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Compute area and square of width of polygons.
!
  pimtol = pi - tol

  do k = 1, npolg

    nvrt = 0
    i = hvl(k)

    do

      if ( iang(i) < pimtol ) then
        nvrt = nvrt + 1
      end if

      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    if ( maxwk < 2 * nvrt ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DSMDF2 - Fatal error!'
      write ( *, '(a)' ) '  MAXWK < 2 * NVRT.'
      write ( *, '(a,i6)' ) '  NVRT = ', nvrt
      write ( *, '(a,i6)' ) '  MAXWK = ', maxwk
      ierror = 7
      return
    end if

    xc = 0

    do

      if ( iang(i) < pimtol ) then
        j = pvl(loc,i)
        xc = xc + 1
        wk(xc) = vcl(1,j)
        wk(xc+nvrt) = vcl(2,j)
      end if

      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    xc = 1
    yc = xc + nvrt
    area(k) = areapg ( nvrt, wk(xc), wk(yc) ) * 0.5D+00

    if ( hflag ) then

      call width2 ( nvrt, wk(xc), wk(yc), i, j, widsq(k), ierror )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DSMDF2 - Fatal error!'
        write ( *, '(a)' ) '  WIDTH2 returns an error condition.'
        return
      end if

    end if

  end do
!
!  Set up IVRT, XIVRT, EDGVAL, VRTVAL arrays.
!
  l = 1

  do k = 1, npolg

    xivrt(k) = l
    i = hvl(k)
    il = pvl(loc,i)

    do

      ivrt(l) = il
      j = pvl(succ,i)
      jl = pvl(loc,j)

      if ( hflag ) then

        s = min ( (vcl(1,jl) - vcl(1,il) )**2 + ( vcl(2,jl) - vcl(2,il) )**2, &
          widsq(k) )

        m = pvl(edgv,i)
        if ( 0 < m ) then
          s = min ( s, widsq(pvl(polg,m) ) )
        end if

        edgval(l) = s

      end if

      l = l + 1
      i = j
      il = jl

      if ( i == hvl(k) ) then
        exit
      end if

    end do

  end do

  xivrt(npolg+1) = l

  if ( .not. hflag ) then
    return
  end if

  vrtval(1:nvc) = 0.0D+00

  do k = 1, npolg

    j = xivrt(k+1) - 1
    l = j

    do i = xivrt(k),l

      il = ivrt(i)

      if ( vrtval(il) == 0.0D+00 ) then
        vrtval(il) = min ( edgval(i), edgval(j) )
      else
        vrtval(il) = min ( vrtval(il), edgval(i), edgval(j) )
      end if

      j = i

    end do

  end do

  return
end
subroutine dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv,  &
  maxho, npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, holv, htsiz,  &
  maxedg, ht, edge, map, ierror )

!*****************************************************************************80
!
!! DSPGDC initializes the polygonal decomposition data structure.
!
!  Purpose: 
!
!    Initialize the polygonal decomposition data structure
!    given an initial decomposition of a polygonal region which
!    may have holes and/or cut, separator, and hole interfaces.
!    Holes and hole interfaces must be simple polygons.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVC, the number of distinct vertex coordinates in region.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinates of 
!    boundary curves in arbitrary order.
!
!    Input, integer ( kind = 4 ) INCR, a positive integer greater than or equal to NVC, 
!    e.g. 10000, added to some elements of IVRT array.
!
!    Input, integer ( kind = 4 ) NCUR, the number of boundary curves (includes outer boundary
!    curves of subregions and boundary curves of holes
!    and hole interfaces).
!
!    Input, integer ( kind = 4 ) NVBC(1:NCUR), the number of vertices per boundary curve.
!
!    Input, integer ( kind = 4 ) ICUR(1:NCUR), indicates type and location of the curves:
!    ICUR(I) = 0 if Ith curve is outer boundary curve,
!    ICUR(I) = K if Ith curve is a hole and is inside
!      the subregion to the left of Kth curve,
!    ICUR(I) = -K if Ith curve is a hole interface and is
!      inside the subregion to the left of Kth curve.
!    K must be the index of an outer or hole interface
!    boundary curve (hole interfaces may be nested).
!    If the Ith curve is inside more than one subregion
!    due to nesting of hole interfaces, then the subregion
!    to the left of Kth curve must be the smallest
!    subregion containing the Ith curve.
!
!    Input, integer ( kind = 4 ) IVRT(1:NV), indices in VCL of vertices of boundary curves;
!    NV = NVBC(1) + ... + NVBC(NCUR); the vertices of each
!    boundary curve must be in counter clockwise order; the first NVBC(1)
!    positions of IVRT are used for the first curve; the
!    next NVBC(2) positions are used for second curve, etc.
!    If the Ith curve is the outer boundary of a subregion
!    determined from cut and separator interfaces, then the
!    elements of IVRT which correspond to this curve are used
!    both for an index in VCL and indicating the type of the
!    edge joining a vertex and its successor as follows.
!    Let J be in range of positions used for the Ith curve
!    and K be the index in VCL of the coordinates of a vertex
!    of the Ith curve. Consider the edge originating from this
!    vertex. IVRT(J) = -K if the edge is part of a cut or
!    separator interface (i.e. there is a subregion to right
!    of edge). IVRT(J) = K if the edge is part of the outer
!    boundary of the region (i.e. the unbounded exterior of
!    the region is to the right of edge). IVRT(J) = K + INCR
!    if the edge is part of the boundary of a hole (i.e.
!    there is a bounded area to the right of edge which is
!    not in the region. If the Ith curve is the boundary of
!    a hole or hole interface, then only IVRT(J) = K is used.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM arrays, 
!    should be greater than or equal to NCUR + (number of hole interfaces).
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays; 
!    should be greater than or equal to NVERT (see below).
!
!    Input, integer ( kind = 4 ) MAXHO, the maximum size available for HOLV array; should be
!    greater than or equal to NHOLE*2 + NHOLA (see below).
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT; should be a prime number 
!    which is about NSC/2 where NSC is number of separator and cut
!    interface edges.
!
!    Input, integer ( kind = 4 ) MAXEDG, the maximum size available for EDGE array; should 
!    be at least NSC.
!
!    Output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions, set to number 
!    of outer subregion boundaries plus number of hole interfaces.
!
!    Output, integer ( kind = 4 ) NVERT, the number of vertices in PVL, set to NV plus number
!    of vertices in holes and hole interfaces (< 2*NV).
!
!    Output, integer ( kind = 4 ) NHOLE, the number of holes and hole interfaces.
!
!    Output, integer ( kind = 4 ) NHOLA, the number of 'attached' holes; these holes are 
!    attached to the outer boundary of a subregion through vertices
!    or cut interfaces and have their edges in consecutive
!    order on the boundary (<= NV/4).
!
!    Output, integer ( kind = 4 ) REGNUM(1:NPOLG). region numbers to left of outer and hole
!    interface boundary curves, which are set to the indices
!    of ICUR or NVBC; this array may be useful in some
!    applications for identifying which original region a
!    subpolygon belongs to.
!
!    Output, HVL(1:NPOLG+NHOLE), the head vertex list; the first NPOLG
!    positions contain the head vertex (index in PVL) of an
!    outer or hole interface boundary curve in which the
!    vertices of the curve are in counter clockwise order in PVL; next
!    NHOLE positions contain the head vertex of a hole or
!    hole interface in which vertices are in CW order in PVL.
!
!    Output, PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT), the polygon 
!    vertex list and interior angles; contains the 5 'arrays' LOC, POLG, SUCC
!    EDGV, IANG (the first 4 are integer arrays, the last
!    is a double precision array); the vertices of each
!    polygon (except for holes) are stored in counter clockwise order in a
!    circular linked list. PVL(LOC,V) is the location in VCL
!    of the coordinates of 'vertex' (index) V. IANG(V) is
!    the interior angle at vertex V. PVL(POLG,V) is polygon
!    number (index of HVL) of subregion containing vertex V
!    (this entry is different from the polygon index only
!    for holes). PVL(SUCC,V) is index in PVL of successor
!    vertex of vertex V. PVL(EDGV,V) gives information about
!    the edge joining vertices V and its successor - if the
!    edge is part of 1 polygon then PVL(EDGV,V) = 0; if the
!    edge is common to 2 polygons then 0 < PVL(EDGV,V) and
!    is equal to the index in PVL of the successor vertex
!    as represented in the other polygon; i.e. in latter
!    case, PVL(LOC,PVL(EDGV,V)) = PVL(LOC,PVL(SUCC,V)) and
!    PVL(EDGV,PVL(EDGV,V)) = V.
!
!    Output, integer ( kind = 4 ) HOLV(1:NHOLE*2+NHOLA), indices in PVL of top or bottom 
!    vertex of holes; first (next) NHOLE entries are for top (bottom)
!    vertices of holes and hole interfaces, with top (bottom)
!    vertices sorted in decreasing (increasing) lexicographic
!    (y,x) order of coord; last NHOLA entries are for attached
!    holes; if bottom vertex of attached hole is a simple
!    vertex of boundary curve containing the hole then entry
!    contains index of bottom vertex otherwise entry contains
!    index of top vertex (which is simple).
!
!    Workspace, integer MAP(1:NCUR), used for mapping input boundary curve 
!    numbers to polygon numbers.
!
!    Workspace, HT(0:HTSIZ-1), EDGE(1:4,1:MAXEDG) - hash table and edge records
!    used to determine matching occurrences of separator or
!    cut interface edges by calling routine EDGHT.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 1, 2, 4, 5, 215, or 216.
!
  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) maxedg
  integer ( kind = 4 ) maxho
  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) ncur

  real ( kind = 8 ) angle
  integer ( kind = 4 ) edge(4,maxedg)
  integer ( kind = 4 ), parameter :: edgv = 4
  logical first
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) holv(maxho)
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) icur(ncur)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) ivs
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jend
  integer ( kind = 4 ) jstr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) kmin
  integer ( kind = 4 ) kpoly
  integer ( kind = 4 ) l
  integer ( kind = 4 ) last
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lvp
  integer ( kind = 4 ) lvs
  integer ( kind = 4 ) map(ncur)
  integer ( kind = 4 ) mpoly
  integer ( kind = 4 ) nh2
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) nholi
  integer ( kind = 4 ) nht
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvbc(ncur)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  integer ( kind = 4 ), parameter :: succ = 3
  real ( kind = 8 ) vcl(2,nvc)
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  ierror = 0
  nhola = 0
  nhole = 0
  nholi = 0
  nvert = 0

  do i = 1, ncur

    nvert = nvert + nvbc(i)

    if ( 0 < icur(i) ) then
      nhole = nhole + 1
    else if ( icur(i) < 0 ) then
      nholi = nholi + 1
      nvert = nvert + nvbc(i)
    end if

  end do

  npolg = ncur - nhole
  ipoly = 0
  iv = 0
  nv = 0
  hdfree = 0
  last = 0
  nht = 0

  ht(0:htsiz-1) = 0

  if ( maxhv < ncur + nholi ) then
    ierror = 4
    return
  else if ( maxpv < nvert ) then
    ierror = 5
    return
  else if ( maxho < ( nhole + nholi ) * 2 ) then
    ierror = 2
    return
  end if
!
!  Initialize REGNUM, HVL, PVL arrays for outer boundary curves.
!
  do i = 1, ncur

    if ( icur(i) /= 0 ) then
      map(i) = 0
      nv = nv + nvbc(i)
      cycle
    end if

    ipoly = ipoly + 1
    regnum(ipoly) = i
    hvl(ipoly) = iv + 1
    map(i) = ipoly
    jstr = nv + 1
    jend = nv + nvbc(i)

    do j = jstr, jend

      iv = iv + 1
      pvl(loc,iv) = abs ( ivrt(j) )
      pvl(polg,iv) = ipoly
      pvl(succ,iv) = iv + 1

      if ( 0 < ivrt(j) ) then
        pvl(edgv,iv) = 0
      else
!
!  The edge originating from current vertex is on a cut or
!  separator interface. Search in hash table for edge, and
!  insert or delete edge.  Set EDGV value if possible.
!
         lv = abs ( ivrt(j) )
         if ( incr < lv ) then
           lv = lv - incr
         end if

         if ( j < jend ) then
           lvs = abs ( ivrt(j+1) )
         else
           lvs = abs ( ivrt(jstr) )
         end if

         if ( incr < lvs ) then
           lvs = lvs - incr
         end if

         call edght ( lv, lvs, iv, nvc, htsiz, maxedg, hdfree, last, ht, &
           edge, ivs, ierror )

         if ( ierror /= 0 ) then
           return
         end if

         if ( 0 < ivs ) then
           pvl(edgv,iv) = ivs
           pvl(edgv,ivs) = iv
           nht = nht - 1
         else
           nht = nht + 1
         end if

       end if

     end do

     pvl(succ,iv) = hvl(ipoly)

     nv = nv + nvbc(i)

  end do

  if ( nht /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DSPGDC - Fatal error!'
    write ( *, '(a)' ) '  NHT /= 0.'
    ierror = 215
    return
  end if
!
!  Initialize REGNUM, HVL, PVL arrays for the hole interfaces.
!
  if ( nholi == 0 ) then
    go to 100
  end if

  do i = 1, ncur

    if ( icur(i) < 0 ) then
      ipoly = ipoly + 1
      map(i) = ipoly
    end if

  end do

  nv = 0

  do i = 1, ncur

    if ( icur(i) < 0 ) then

      ipoly = ipoly + 1
      kpoly = ipoly - nholi
      mpoly = map(-icur(i))
      regnum(kpoly) = i
      hvl(kpoly) = iv + 1
      hvl(ipoly) = iv + 2
      jstr = nv + 1
      jend = nv + nvbc(i)

      do j = jstr, jend
        iv = iv + 2
        pvl(loc,iv-1) = ivrt(j)
        pvl(polg,iv-1) = kpoly
        pvl(succ,iv-1) = iv + 1
        pvl(edgv,iv-1) = iv + 2
        pvl(loc,iv) = ivrt(j)
        pvl(polg,iv) = mpoly
        pvl(succ,iv) = iv - 2
        pvl(edgv,iv) = iv - 3
      end do

      pvl(succ,iv-1) = hvl(kpoly)
      pvl(edgv,iv-1) = hvl(ipoly)
      pvl(succ,hvl(ipoly)) = iv
      pvl(edgv,hvl(ipoly)) = iv - 1

    end if

    nv = nv + nvbc(i)
   
  end do
!
!  Initialize HVL, PVL arrays for the ordinary holes.
!
100 continue

  if ( nhole /= 0 ) then

    nv = 0

    do i = 1, ncur

      if ( 0 < icur(i) ) then

        ipoly = ipoly + 1
        mpoly = map(icur(i))
        hvl(ipoly) = iv + 1
        jstr = nv + 1
        jend = nv + nvbc(i)

        do j = jstr, jend
          iv = iv + 1
          pvl(loc,iv) = ivrt(j)
          pvl(polg,iv) = mpoly
          pvl(succ,iv) = iv - 1
          pvl(edgv,iv) = 0
        end do

        pvl(succ,hvl(ipoly)) = iv

      end if

      nv = nv + nvbc(i)

    end do

  end if
!
!  Determine bottom or top simple vertex of attached holes.
!
  nhole = nhole + nholi
  nh2 = nhole + nhole
  j1 = 0
  j2 = 0

  do i = 1, npolg-nholi

   j = hvl(i)

  150    continue

      if ( incr < pvl(loc,j) ) then
         j = pvl(succ,j)
         if ( j /= hvl(i) ) then
          go to 150
         else
          ierror = 216
          return
         end if
      end if

   first = .true.

  160    continue

      lv = pvl(loc,j)

      if ( 0 < j1 ) then
        if ( lv <= incr ) then
          j2 = j
        else if ( lv - incr == lvs ) then
          j2 = j
        else
          pvl(loc,j) = lv - incr
        end if
      else if ( incr < lv ) then
        j1 = j
        lvs = lv - incr
        pvl(loc,j) = lvs
      end if

      if ( 0 < j2 ) then
!
!  (Part of) hole starts at vertex J1 and ends at J2.
!
         if ( lv <= incr .and. lv /= lvs ) then
           go to 180
         end if

         k = j1

  170          continue

          if ( k == j1 ) then
            kmin = k
            kmax = k
            xmin = vcl(1,lvs)
            ymin = vcl(2,lvs)
            xmax = xmin
            ymax = ymin

          else

            l = pvl(loc,k)
            x = vcl(1,l)
            y = vcl(2,l)

            if ( y < ymin .or. y == ymin .and. x < xmin ) then
              kmin = k
              xmin = x
              ymin = y
            else if ( ymax < y .or. y == ymax .and. xmax < x ) then
              kmax = k
              xmax = x
              ymax = y
            end if

          end if

          k = pvl(succ,k)
         if ( k /= j2 ) then
           go to 170
         end if

         if ( kmin == j1 ) then
           kmin = kmax
         end if

         nhola = nhola + 1

         if ( maxho < nh2 + nhola ) then
          ierror = 2
          return
         end if

         holv(nh2+nhola) = kmin

180      continue

         j1 = 0
         j2 = 0

         if ( incr < lv ) then
           j1 = j
           pvl(loc,j) = lvs
         end if

      end if
      j = pvl(succ,j)

   if ( first ) then
     first = .false.
     jend = j
     go to 160
   else if ( j /= jend ) then
     go to 160
   end if

  end do
!
!  Initialize the IANG array.
!
  do i = 1, npolg+nhole

    j = hvl(i)
    lvp = pvl(loc,j)
    iv = pvl(succ,j)
    lv = pvl(loc,iv)

    do

      ivs = pvl(succ,iv)
      lvs = pvl(loc,ivs)
      iang(iv) = angle ( vcl(1,lvp), vcl(2,lvp), vcl(1,lv), vcl(2,lv), &
        vcl(1,lvs), vcl(2,lvs) )

      if ( iv == j ) then
        exit
      end if

      lvp = lv
      iv = ivs
      lv = lvs

    end do

  end do
!
!  Initialize HOLV array.
!
  if ( 0 < nhole ) then
    call holvrt ( nhole, vcl, hvl(npolg+1), pvl, holv )
  end if

  return
end
subroutine dtris2 ( npt, vcl, ind, ntri, til, tnbr, ierror )

!*****************************************************************************80
!
!! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!    Note that DTRIS2 or RTRIS2, the fundamental routine for constructing
!    the Delaunay triangulation, alters the input coordinate data by
!    sorting it.  This has caused me so many problems that I finally 
!    wrote a modified version of DTRIS2/RTRIS2 that undoes the sorting
!    before return.  In all other programs that use DTRIS2/RTRIS2, I
!    use the modified version, but I have left the original here in this
!    package.
!
!    On abnormal return, IERROR is set to 8, 224, or 225.
!
!  Modified:
!
!    22 July 2003
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPT, the number of vertices.
!
!    Input, real ( kind = 8 ) VCL(2,NPT), the coordinates of the vertices.
!
!    Input/output, integer ( kind = 4 ) IND(NPT), the indices in VCL of the vertices 
!    to be triangulated.  On output, IND has been permuted by the sort.
!
!    Output, integer ( kind = 4 ) NTRI, the number of triangles in the triangulation; 
!    NTRI is equal to 2*NPT - NB - 2, where NB is the number of boundary 
!    vertices.
!
!    Output, integer ( kind = 4 ) TIL(3,NTRI), the nodes that make up each triangle.
!    The elements are indices of VCL.  The vertices of the triangles are 
!    in counter clockwise order.
!
!    Output, integer ( kind = 4 ) TNBR(3,NTRI), the triangle neighbor list.
!    Positive elements are indices of TIL; negative elements are used for links
!    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
!    where I, J = triangle, edge index; TNBR(J,I) refers to
!    the neighbor along edge from vertex J to J+1 (mod 3).
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, nonzero if an error occurred.
!
  integer ( kind = 4 ) npt

  real ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind(npt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(npt)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,npt*2)
  integer ( kind = 4 ) tnbr(3,npt*2)
  real ( kind = 8 ) tol
  integer ( kind = 4 ) top
  real ( kind = 8 ) vcl(2,npt)

  maxst = npt

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  ierror = 0
!
!  Sort the vertices.
!
  call dhpsrt ( 2, npt, 2, vcl, ind )
!
!  Ensure that no two consecutive points are too close.
!
  m1 = ind(1)

  do i = 2, npt

    m = m1
    m1 = ind(i)

    k = 0
    do j = 1, 2

      cmax = max ( abs ( vcl(j,m) ), abs ( vcl(j,m1) ) )

      if ( tol * cmax < abs ( vcl(j,m) - vcl(j,m1) ) .and. tol < cmax ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      ierror = 224
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Consecutive points are too close.'
      return
    end if

  end do
!
!  Take the first two points, M1 and M2, and find a suitable non-collinear
!  third, M.  All points between M2 and M are very close to collinear
!  with M1 and M2.
!
  m1 = ind(1)
  m2 = ind(2)
  j = 3

  do

    if ( npt < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Could not find a noncollinear third point.'
      ierror = 225
      return
    end if

    m = ind(j)
    lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
      vcl(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1
  
  end do

  ntri = j - 2
!
!  Depending on the orientation of M1, M2, and M, set up the initial
!  triangle data.
!
  if ( lr == -1 ) then

    til(1,1) = m1
    til(2,1) = m2
    til(3,1) = m
    tnbr(3,1) = -3

    do i = 2, ntri

      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m1
      til(2,i) = m2
      til(3,i) = m
      tnbr(1,i-1) = -3 * i
      tnbr(2,i-1) = i
      tnbr(3,i) = i - 1

    end do

    tnbr(1,ntri) = -3 * ntri - 1
    tnbr(2,ntri) = -5
    ledg = 2
    ltri = ntri

  else

    til(1,1) = m2
    til(2,1) = m1
    til(3,1) = m
    tnbr(1,1) = -4

    do i = 2, ntri
      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m2
      til(2,i) = m1
      til(3,i) = m
      tnbr(3,i-1) = i
      tnbr(1,i) = -3 * i - 3
      tnbr(2,i) = i - 1
    end do

    tnbr(3,ntri) = -3 * ntri
    tnbr(2,1) = -3 * ntri - 2
    ledg = 2
    ltri = 1

  end if

  if ( msglvl == 4 ) then

    m2 = ind(1)
    write ( *, '(i7,4f15.7)' ) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)

    do i = 2, j-1
      m1 = m2
      m2 = ind(i)
      write ( *, '(i7,4f15.7)' ) 1,vcl(1,m1),vcl(2,m1),vcl(1,m2),vcl(2,m2)
      write ( *, '(i7,4f15.7)' ) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    end do

  end if
!
!  Insert vertices one at a time from outside the convex hull, determine
!  the visible boundary edges, and apply diagonal edge swaps until
!  the Delaunay triangulation of the vertices (so far) is obtained.
!
  top = 0

  do i = j+1, npt

    if ( msglvl == 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(i7)' ) i
    end if

    m = ind(i)
    m1 = til(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = til(ledg+1,ltri)
    else
      m2 = til(1,ltri)
    end if

    lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
      vcl(2,m2), 0.0D+00 )

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tnbr(ledg,ltri)
      rtri = l / 3
      redg = mod ( l, 3 ) + 1
    end if

    call vbedg ( vcl(1,m), vcl(2,m), vcl, til, tnbr, ltri, ledg, rtri, redg )

    n = ntri + 1
    l = -tnbr(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -tnbr(e,t)
      m2 = til(e,t)

      if ( e <= 2 ) then
        m1 = til(e+1,t)
      else
        m1 = til(1,t)
      end if

      ntri = ntri + 1
      tnbr(e,t) = ntri
      til(1,ntri) = m1
      til(2,ntri) = m2
      til(3,ntri) = m
      tnbr(1,ntri) = t
      tnbr(2,ntri) = ntri - 1
      tnbr(3,ntri) = ntri + 1
      top = top + 1

      if ( maxst < top ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Ran out of stack space.'
        write ( *, '(a,i6)' ) '  MAXST = ', maxst
        write ( *, '(a,i6)' ) '  TOP =   ', top
        ierror = 8
        return
      end if

      stack(top) = ntri

      if ( msglvl == 4 ) then
        write ( *, '(i7,4f15.7)' ) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
      end if

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    if ( msglvl == 4 ) then
      write ( *, '(i7,4f15.7)' ) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
    end if

    tnbr(ledg,ltri) = -3 * n - 1
    tnbr(2,n) = -3 * ntri - 2
    tnbr(3,ntri) = -l
    ltri = n
    ledg = 2

    call swapec ( m, top, maxst, ltri, ledg, vcl, til, tnbr, stack, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a,i6)' ) '  SWAPEC returned IERROR = ', ierror
      return
    end if

  end do

  if ( msglvl == 4 ) then
    write ( *, '(i7)' ) npt + 1
  end if

  return
end
subroutine dtriw2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierror )

!*****************************************************************************80
!
!! DTRIW2 constructs an incremental Delaunay triangulation in 2D.
!
!  Purpose: 
!
!    Construct the Delaunay triangulation of 2D vertices using
!    incremental approach and diagonal edge swaps. Vertices are
!    inserted one at a time in order given by IND array. The initial
!    triangles created due to a new vertex are obtained by a walk
!    through the triangulation until location of vertex is known.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPT, the number of 2D points (vertices).
!
!    Input, integer ( kind = 4 ) MAXST, the maximum size available for STACK array; should 
!    be about NPT to be safe, but MAX(10,2*LOG2(NPT)) usually enough.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) IND(1:NPT), indices in VCL of vertices to be triangulated;
!    vertices are inserted in order given by this array.
!
!    Output, integer ( kind = 4 ) NTRI, the number of triangles in triangulation; equal to
!    2*NPT - NB - 2 where NB = number of boundary vertices.
!
!    Output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list; elements 
!    are indices of VCL; vertices of triangles are in counter clockwise order.
!
!    Output, integer ( kind = 4 ) TNBR(1:3,1:NTRI), the triangle neighbor list; positive 
!    elements are indices of TIL; negative elements are used for links
!    of counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
!    where I, J = triangle, edge index; TNBR(J,I) refers to
!    the neighbor along edge from vertex J to J+1 (mod 3).
!
!    Workspace, integer STACK(1:MAXST), used for stack of triangles for which
!    circumcircle test must be made.
!
!    Output, integer ( kind = 4 ) IERROR, error flag. For abnormal return,
!    IERROR is set to 8, 224, 225, or 226.
!
  integer ( kind = 4 ) maxst
  integer ( kind = 4 ) npt

  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  real ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind(npt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(maxst)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,npt*2)
  integer ( kind = 4 ) tnbr(3,npt*2)
  integer ( kind = 4 ) top
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(2,*)

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Determine the initial triangle.
!
  m1 = ind(1)
  m2 = ind(2)

  do j = 1, 2
    cmax = max ( abs ( vcl(j,m1) ), abs ( vcl(j,m2) ) )
    if ( tol * cmax < abs ( vcl(j,m1) - vcl(j,m2) ) .and. tol < cmax ) then
      go to 20
    end if
  end do

  ierror = 224
  return

20 continue

  i3 = 3

30 continue

  if ( npt < i3 ) then
    ierror = 225
    return
  end if

  m = ind(i3)
  lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
    vcl(2,m2), 0.0D+00 )

  if ( lr == 0 ) then
    i3 = i3 + 1
    go to 30
  end if

  if ( i3 /= 3 ) then
    ind(i3) = ind(3)
    ind(3) = m
  end if

  ntri = 1

  if ( lr == -1 ) then
    til(1,1) = m1
    til(2,1) = m2
  else
    til(1,1) = m2
    til(2,1) = m1
  end if

  til(3,1) = m
  tnbr(1,1) = -4
  tnbr(2,1) = -5
  tnbr(3,1) = -3

  if ( msglvl == 4 ) then
    write ( *,600) 1,vcl(1,m1),vcl(2,m1),vcl(1,m2),vcl(2,m2)
    write ( *,600) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
  end if
!
!  Insert vertices one at a time from anywhere.
!  Walk through the triangulation to determine the location of the new vertex.
!  Apply diagonal edge swaps until Delaunay triangulation of vertices
!  (so far) is obtained.
!
  top = 0

  do i = 4, npt

    if ( msglvl == 4 ) then
      write ( *,600) i
    end if

    m = ind(i)
    rtri = ntri

    call walkt2 ( vcl(1,m), vcl(2,m), ntri, vcl, til, tnbr, rtri, redg, ierror )

    if ( redg == 0 ) then

      m1 = til(1,rtri)
      m2 = til(2,rtri)
      m3 = til(3,rtri)
      til(3,rtri) = m

      if ( 0 < tnbr(1,rtri) ) then
        top = 1
        stack(top) = rtri
      end if

      ntri = ntri + 1
      til(1,ntri) = m2
      til(2,ntri) = m3
      til(3,ntri) = m
      n = tnbr(2,rtri)
      tnbr(1,ntri) = n

      if ( 0 < n ) then
        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if
        top = top + 1
        stack(top) = ntri
      end if

      ntri = ntri + 1
      til(1,ntri) = m3
      til(2,ntri) = m1
      til(3,ntri) = m
      n = tnbr(3,rtri)
      tnbr(1,ntri) = n

      if ( 0 < n ) then
        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if
        top = top + 1
        stack(top) = ntri
      end if

      tnbr(2,rtri) = ntri - 1
      tnbr(3,rtri) = ntri
      tnbr(2,ntri-1) = ntri
      tnbr(3,ntri-1) = rtri
      tnbr(2,ntri) = rtri
      tnbr(3,ntri) = ntri - 1

      if ( tnbr(1,ntri-1) <= 0 ) then

        t = rtri
        e = 1

        do

          if ( tnbr(e,t) <= 0 ) then
            exit
          end if

          t = tnbr(e,t)

          if ( til(1,t) == m2 ) then
            e = 3
          else if ( til(2,t) == m2 ) then
            e = 1
          else
            e = 2
          end if

        end do

        tnbr(e,t) = -3 * ntri + 3

      end if

      if ( tnbr(1,ntri) <= 0 ) then

        t = ntri - 1
        e = 1

        do

          if ( tnbr(e,t) <= 0 ) then
            exit
          end if

          t = tnbr(e,t)
          if ( til(1,t) == m3 ) then
            e = 3
          else if ( til(2,t) == m3 ) then
            e = 1
          else
            e = 2
          end if

        end do

        tnbr(e,t) = -3 * ntri

      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
      end if

   else if ( redg < 0 ) then

      redg = -redg
      ltri = 0
      call vbedg ( vcl(1,m), vcl(2,m), vcl, til, tnbr, ltri, ledg, rtri, redg )
      n = ntri + 1
      l = -tnbr(ledg,ltri)

60    continue

         t = l / 3
         e = mod ( l, 3 ) + 1
         l = -tnbr(e,t)
         m2 = til(e,t)

         if ( e <= 2 ) then
           m1 = til(e+1,t)
         else
           m1 = til(1,t)
         end if

         ntri = ntri + 1
         tnbr(e,t) = ntri
         til(1,ntri) = m1
         til(2,ntri) = m2
         til(3,ntri) = m
         tnbr(1,ntri) = t
         tnbr(2,ntri) = ntri - 1
         tnbr(3,ntri) = ntri + 1
         top = top + 1

         if ( maxst < top ) then
           ierror = 8
           go to 100
         end if
 
        stack(top) = ntri

         if ( msglvl == 4 ) then
           write (*,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
         end if

      if ( t /= rtri .or. e /= redg ) then
        go to 60
      end if

      if ( msglvl == 4 ) then
        write (*,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
      end if

      tnbr(ledg,ltri) = -3*n - 1
      tnbr(2,n) = -3*ntri - 2
      tnbr(3,ntri) = -l

   else if ( redg <= 3 ) then

      m1 = til(redg,rtri)

      if ( redg == 1 ) then
        e = 2
        ep1 = 3
      else if ( redg == 2 ) then
        e = 3
        ep1 = 1
      else
        e = 1
        ep1 = 2
      end if

      m2 = til(e,rtri)
      til(e,rtri) = m
      m3 = til(ep1,rtri)

      if ( 0 < tnbr(ep1,rtri) ) then
        top = 1
        stack(top) = rtri
      end if

      ntri = ntri + 1
      til(1,ntri) = m
      til(2,ntri) = m2
      til(3,ntri) = m3
      n = tnbr(e,rtri)
      tnbr(2,ntri) = n
      tnbr(3,ntri) = rtri
      tnbr(e,rtri) = ntri

      if ( 0 < n ) then
        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if
        top = top + 1
        stack(top) = ntri
      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
      end if

      ltri = tnbr(redg,rtri)

      if ( ltri <= 0 ) then
        tnbr(1,ntri) = ltri
        tnbr(redg,rtri) = -3*ntri
        if ( tnbr(2,ntri) <= 0 ) then
          tnbr(1,ntri) = -3*ntri - 1
        end if

      else

        tnbr(1,ntri) = ntri + 1
        tnbr(redg,rtri) = ltri

        if ( til(1,ltri) == m2 ) then
          ledg = 1
          em1 = 2
          e = 3
        else if ( til(2,ltri) == m2 ) then
          ledg = 2
          em1 = 3
          e = 1
        else
          ledg = 3
          em1 = 1
          e = 2
        end if

        til(ledg,ltri) = m
        m3 = til(e,ltri)
        if ( 0 < tnbr(em1,ltri) ) then
          top = top + 1
          stack(top) = ltri
        end if
        ntri = ntri + 1
        til(1,ntri) = m2
        til(2,ntri) = m
        til(3,ntri) = m3
        tnbr(1,ntri) = ntri - 1
        tnbr(2,ntri) = ltri
        n = tnbr(e,ltri)
        tnbr(3,ntri) = n
        tnbr(e,ltri) = ntri

        if ( 0 < n ) then
          if ( tnbr(1,n) == ltri ) then
            tnbr(1,n) = ntri
          else if ( tnbr(2,n) == ltri ) then
            tnbr(2,n) = ntri
          else
            tnbr(3,n) = ntri
          end if
          top = top + 1
          stack(top) = ntri
        end if

        if ( msglvl == 4 ) then
          write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
        end if

        if ( tnbr(2,ntri-1) <= 0 ) then

          t = ntri
          e = 3

          do

            if ( tnbr(e,t) <= 0 ) then
              exit
            end if

            t = tnbr(e,t)
            if ( til(1,t) == m2 ) then
              e = 3
            else if ( til(2,t) == m2 ) then
              e = 1
            else
              e = 2
            end if

          end do

          tnbr(e,t) = -3 * ntri + 2

        end if

         if ( tnbr(3,ntri) <= 0 ) then

          t = ltri

          if ( ledg <= 2 ) then
             e = ledg + 1
          else
             e = 1
          end if

          do

            if ( tnbr(e,t) <= 0 ) then
              exit
            end if

            t = tnbr(e,t)
            if ( til(1,t) == m3 ) then
              e = 3
            else if ( til(2,t) == m3 ) then
              e = 1
            else
              e = 2
            end if

          end do

          tnbr(e,t) = -3 * ntri - 2

         end if

      end if

    else
      ierror = 224
      go to 100
    end if

    btri = 0
    bedg = 0

    call swapec ( m, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIW2 - Fatal error!'
      write ( *, '(a,i6)' ) '  SWAPEC returned error flag IERROR = ', ierror
      exit
    end if

  end do

100 continue

  if ( i3 /= 3 ) then
    t = ind(i3)
    ind(i3) = ind(3)
    ind(3) = t
  end if

  if ( msglvl == 4 ) then
    write ( *,600) npt+1
  end if

600 format (1x,i7,4f15.7)

  return
end
subroutine edght ( a, b, v, n, htsiz, maxedg, hdfree, last, ht, edge, w, &
  ierror )

!*****************************************************************************80
!
!! EDGHT searches a hash table for a record in EDGE containing key (A,B).
!
!  Purpose: 
!
!    Search in hash table HT for record in EDGE containing
!    key (A,B).
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, vertex indices, greater than 0, of edge (also 
!    key of hash table).
!
!    Input, integer ( kind = 4 ) V, value associated with edge.
!
!    Input, integer ( kind = 4 ) N, upper bound on A, B.
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT.
!
!    Input, integer ( kind = 4 ) MAXEDG, the maximum size available for EDGE array.
!
!    Input/output, integer ( kind = 4 ) HDFREE, head pointer to linked list of free entries 
!    of EDGE array due to deletions.  Before first call to this routine, HDFREE
!    should be set to 0.
!
!    Input/output, integer ( kind = 4 ) LAST, index of last entry used in EDGE array.  
!    Before first call to this routine, LAST should be set to 0.
!
!    Input/output, integer ( kind = 4 ) HT(0:HTSIZ-1), hash table of head pointers (direct 
!    chaining with ordered lists is used).  Before first call to this routine, 
!    entries of HT should be set to 0.  If key with A,B is found then this 
!    record is deleted from hash table, else record is inserted in hash table.
!
!    Input/output, integer ( kind = 4 ) EDGE(1:4,1:MAXEDG), entries of hash table records;
!      EDGE(1,I) = MIN(A,B); EDGE(2,I) = MAX(A,B);
!      EDGE(3,I) = V; EDGE(4,I) = link
!    If key with A,B is found then this record is deleted
!    from hash table, else record is inserted in hash table.
!
!    Output, integer ( kind = 4 ) W, EDGE(3,INDEX), where INDEX is index of record, if found;
!    else 0.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  For abnormal return,
!    IERROR is set to 1
!
  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) maxedg

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) bptr
  integer ( kind = 4 ) edge(4,maxedg)
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) last
  integer ( kind = 4 ) n
  integer ( kind = 4 ) newp
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w

  ierror = 0

  if ( a < b ) then
    aa = a
    bb = b
  else
    aa = b
    bb = a
  end if

  k = mod ( aa*n + bb, htsiz )
  bptr = -1
  ptr = ht(k)

10 continue

  if ( ptr /= 0 ) then

    if ( aa < edge(1,ptr) ) then

      go to 20

    else if ( edge(1,ptr) == aa ) then

      if ( bb < edge(2,ptr) ) then

        go to 20

      else if ( edge(2,ptr) == bb ) then

        if ( bptr == -1 ) then
          ht(k) = edge(4,ptr)
        else
          edge(4,bptr) = edge(4,ptr)
        end if

        edge(4,ptr) = hdfree
        hdfree = ptr
        w = edge(3,ptr)
        return

      end if

    end if

    bptr = ptr
    ptr = edge(4,ptr)
    go to 10

  end if

20 continue

  if ( 0 < hdfree ) then

    newp = hdfree
    hdfree = edge(4,hdfree)

  else

    last = last + 1
    newp = last

    if ( maxedg < last ) then
      ierror = 1
      return
    end if

  end if

  if ( bptr == -1 ) then
    ht(k) = newp
  else
    edge(4,bptr) = newp
  end if

  edge(1,newp) = aa
  edge(2,newp) = bb
  edge(3,newp) = v
  edge(4,newp) = ptr
  w = 0

  return
end
subroutine eqdis2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid,  &
  nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl,  &
  pvl, iang, area, psi, h, iwk, wk, ierror )

!*****************************************************************************80
!
!! EQDIS2 further subdivides convex polygons for mesh equidistribution.
!
!  Purpose: 
!
!    Further subdivide convex polygons so that an approximately
!    equidistributing triangular mesh can be constructed with
!    respect to a heuristic or a user-supplied mesh distribution
!    function (MDF), and determine triangle size for each polygon of
!    decomposition.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical HFLAG, is .TRUE. if heuristic mdf, .FALSE. if
!    user-supplied mdf.
!
!    Input, external UMDF, a user-supplied mesh distribution function 
!    of the form:
!
!      function umdf ( x, y )
!      double precision umdf
!      double precision x
!      double precision y
!
!    Input, real ( kind = 8 ) KAPPA, the mesh smoothness parameter in 
!    the interval [0.0,1.0].
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in radians
!    used to determine extra points as possible endpoints of separators.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter in radians
!    used in accepting separators.
!
!    Input, real ( kind = 8 ) DMIN, a parameter used to determine if variation 
!    of mdf in polygon is 'sufficiently high'.
!
!    Input, integer ( kind = 4 ) NMIN, a parameter used to determine if 'sufficiently large'
!    number of triangles in polygon.
!
!    Input, integer ( kind = 4 ) NTRID, the desired number of triangles in mesh.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or 
!    positions used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions or 
!    positions used in HVL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of polygon vertices or positions 
!    used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array, should 
!    be greater than or equal to the number of vertex coordinates required
!    for decomposition (approximately NVC + 2*NS where NS is expected number 
!    of new separators).
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM, AREA, 
!    PSI, H arrays; should be greater than or equal to the number of polygons
!    required for decomposition (approximately NPOLG + NS).
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays; 
!    should be greater than or equal to the number of polygon vertices 
!    required for decomposition (approximately NVERT + 5*NS).
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should 
!    be greater than or equal to 
!      MAX(2*NP, NVERT + NPOLG + 3*NVRT + INT(2*PI/ANGSPC))
!    where NVRT is maximum number of vertices in a convex
!    polygon of the (input) decomposition, NP is expected
!    value of NPOLG on output.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be greater than or equal to 
!      NVC + NVERT + 2*NPOLG + 3*(NVRT + INT(2*PI/ANGSPC)).
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, real ( kind = 8 ) HVL(1:NPOLG), head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT), 
!    the polygon vertex list and interior angles; see routine DSPGDC for more 
!    details.  Note that the data structures should be as output from routine
!    CVDEC2.
!
!    Output, real ( kind = 8 ) AREA(1:NPOLG), the area of convex polygons 
!    in the decomposition.
!
!    Output, real ( kind = 8 ) PSI(1:NPOLG), the smoothed mean mdf 
!    values in the convex polygons.
!
!    Output, real ( kind = 8 ) H(1:NPOLG), the triangle size for convex
!    polygons.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  For abnormal return,
!    IERROR is set to 3, 4, 5, 6, 7, 200, 201, or 222.
!
  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  real ( kind = 8 ) area(maxhv)
  real ( kind = 8 ) dmin
  integer ( kind = 4 ) edgval
  real ( kind = 8 ) h(maxhv)
  logical hflag
  integer ( kind = 4 ) hvl(maxhv)
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) iwk(maxiw)
  real ( kind = 8 ) kappa
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) ntrid
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  real ( kind = 8 ) psi(maxhv)
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  real ( kind = 8 ) umdf
  external umdf
  real ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vrtval
  integer ( kind = 4 ) widsq
  real ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xivrt

  ierror = 0
  ivrt = 1
  xivrt = ivrt + nvert
  m = xivrt + npolg

  if ( maxiw < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
    write ( *, '(a)' ) '  MAXIW < M'
    ierror = 6
    return
  end if

  widsq = 1

  if ( hflag ) then

    edgval = widsq + npolg
    vrtval = edgval + nvert
    n = npolg + nvert + nvc

    if ( maxwk < n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
      write ( *, '(a)' ) '  MAXWK < N'
      ierror = 7
      return
    end if

  else

    edgval = 1
    vrtval = 1
    n = 0

  end if

  call dsmdf2 ( hflag, nvc, npolg, maxwk-n, vcl, hvl, pvl, iang, iwk(ivrt), &
    iwk(xivrt), wk(widsq), wk(edgval), wk(vrtval), area, wk(n+1), ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
    write ( *, '(a)' ) '  DSMDF2 returns error condition.'
    return
  end if

  call mfdec2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid, nvc, &
    npolg, nvert, maxvc, maxhv, maxpv, maxiw-m, maxwk-n, vcl, regnum, hvl, &
    pvl, iang, iwk(ivrt), iwk(xivrt), wk(widsq), wk(edgval), wk(vrtval), &
    area, psi, iwk(m+1), wk(n+1), ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
    write ( *, '(a)' ) '  MFDEC2 returns error condition.'
    return
  end if

  if ( maxiw < 2 * npolg ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
    write ( *, '(a)' ) '  MAXIW < 2 * NPOLG.'
    ierror = 6
    return
  end if

  call trisiz ( ntrid, npolg, hvl, pvl, area, psi, h, iwk, iwk(npolg+1) )

  return
end
subroutine fndsep ( angac1, xr, yr, nvrt, xc, yc, ivis, theta, nv, iv,  &
  vcl, pvl, iang, angsep, i1, i2, wkang )

!*****************************************************************************80
!
!! FNDSEP finds separators to resolve a reflex vertex.
!
!  Purpose: 
!
!    Find 1 or 2 separators which can resolve a reflex vertex
!    (XR,YR) using a max-min angle criterion from list of vertices
!    in increasing polar angle with respect to the reflex vertex.
!
!    Preference is given to 1 separator.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGAC1, the angle tolerance parameter used 
!    for preference in accepting one separator.
!
!    Input, real ( kind = 8 ) XR, YR, the coordinates of reflex vertex.
!
!    Input, integer ( kind = 4 ) NVRT, (number of vertices) - 1.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex coordinates 
!    of possible endpoints of a separator.
!
!    Input, integer ( kind = 4 ) IVIS(0:NVRT), contains information about the vertices of
!    XC, YC arrays with respect to the polygon vertex list; if
!    0 < IVIS(I) then vertex (XC(I),YC(I)) has index IVIS(I)
!    in PVL; if IVIS(I) < 0 then vertex (XC(I),YC(I)) is on
!    the edge joining vertices with indices -IVIS(I) and
!    SUCC(-IVIS(I)) in PVL.
!
!    Input, real ( kind = 8 ) THETA(0:NVRT), the polar angles of vertices 
!    in increasing order; THETA(NVRT) is the interior angle of reflex vertex;
!    THETA(I), 0 <= I, is the polar angle of (XC(I),YC(I))
!    with respect to reflex vertex.
!
!    Input, integer ( kind = 4 ) NV, (number of vertices to be considered as endpoint of a
!    separator) - 1.
!
!    Input, integer ( kind = 4 ) IV(0:NV), the indices of vertices in XC, YC arrays to be
!    considered as endpoint of a separator; angle between
!    consecutive vertices is assumed to be < 180 degrees.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), real ( kind = 8 ) IANG(1:*), the polygon 
!    vertex list, interior angles.
!
!    Output, real ( kind = 8 ) ANGSEP, the minimum of the 4 or 7 angles at the
!    boundary resulting from 1 or 2 separators, respectively.
!
!    Output, integer ( kind = 4 ) I1, I2, the indices of endpoints of separators in XC, 
!    YC arrays; I2 = -1 if there is only one separator, else I1 < I2.
!
!    Workspace, real ( kind = 8 ) WKANG(0:NV).
!
  real ( kind = 8 ) ang
  real ( kind = 8 ) angac1
  real ( kind = 8 ) angsep
  real ( kind = 8 ) angsp2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real ( kind = 8 ) iang(*)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) iv(0:nv)
  integer ( kind = 4 ) ivis(0:nvrt)
  real ( kind = 8 ) minang
  integer ( kind = 4 ) p
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) pvl(4,*)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  real ( kind = 8 ) theta(0:nvrt)
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(2,*)
  real ( kind = 8 ) wkang(0:nv)
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) xr
  real ( kind = 8 ) yc(0:nvrt)
  real ( kind = 8 ) yr

  tol = 100.0D+00 * epsilon ( tol )
!
!  Determine the vertices in the inner cone - indices P to Q.
!
  i = 0
  p = -1
  phi = theta(nvrt) - pi + tol

  do while ( p < 0 )

    if ( phi <= theta(iv(i)) ) then
      p = i
    else
      i = i + 1
    end if

  end do

  i = nv
  q = -1
  phi = pi - tol

  do while ( q < 0 )

    if ( theta(iv(i)) <= phi ) then
      q = i
    else
      i = i - 1
    end if

  end do
!
!  Use the max-min angle criterion to find the best separator
!  in inner cone.
!
  angsep = 0.0

  do i = p, q

    k = iv(i)
    ang = minang ( xr, yr, xc(k), yc(k), ivis(k), theta(k), theta(nvrt), &
      vcl, pvl, iang )

    if ( angsep < ang ) then
      angsep = ang
      ii = iv(i)
    end if

  end do

  angsp2 = angsep
  if ( angac1 <= angsep ) then
    go to 110
  end if
!
!  If the best separator in inner cone is not 'good' enough,
!  use max-min angle criterion to try to find a better pair
!  of separators from the right and left cones.
!
  nr = 0
  nl = 0

  do r = 0, p - 1

    wkang(r) = 0.0D+00

    if ( angsep < theta(iv(r)) ) then

      k = iv(r)

      ang = minang ( xr, yr, xc(k), yc(k), ivis(k), theta(k), theta(nvrt), &
        vcl, pvl, iang )

      if ( angsep < ang ) then
        nr = nr + 1
        wkang(r) = ang
      end if

    end if

  end do

  if ( nr == 0 ) then
    go to 110
  end if

  phi = theta(nvrt) - angsep

  do l = q+1, nv

    wkang(l) = 0.0D+00

    if ( theta(iv(l)) < phi ) then

      k = iv(l)
      ang = minang ( xr, yr, xc(k), yc(k), ivis(k), theta(k), theta(nvrt), &
        vcl, pvl, iang )

      if ( angsep < ang ) then
        nl = nl + 1
        wkang(l) = ang
      end if

    end if

  end do

  if ( nl == 0 ) then
    go to 110
  end if
!
!  Check all possible pairs for the best pair of separators
!  in the right and left cones.
!
  m = nv

  do r = p-1, 0, -1

     if ( q < m .and. angsp2 < wkang(r) ) then

      phi = theta(iv(r))

80    continue

      if ( q < m .and. &
           ( wkang(m) <= angsp2 .or. &
             pi - tol < theta(iv(m)) - phi ) ) then
        m = m - 1
        go to 80
      end if

      do l = q+1, m

         if ( angsp2 < wkang(l) ) then

            ang = min ( theta(iv(l)) - phi, wkang(r), wkang(l) )

            if ( angsp2 < ang ) then
             angsp2 = ang
             i1 = iv(r)
             i2 = iv(l)
            end if

         end if

       end do

     end if

  end do
!
!  Choose 1 or 2 separators based on max-min angle criterion or
!  ANGAC1 parameter.
!
  110 continue

  if ( angsp2 <= angsep ) then
    i1 = ii
    i2 = -1
  else
    angsep = angsp2
  end if

  return
end
subroutine fndtri ( iedg, mxtr, sflag, tedg, itr, ind, ierror )

!*****************************************************************************80
!
!! FNDTRI finds two triangles containing a given edge.
!
!  Purpose: 
!
!    Find two triangles containing edge with index IEDG in array TEDG.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IEDG, the index of edge to be searched in TEDG.
!
!    Input, integer ( kind = 4 ) MXTR, the maximum index of triangle to be searched in TEDG.
!
!    Input, logical SFLAG, is .TRUE. if and only if the second triangle is to be 
!    searched from end of array.
!
!    Input, integer ( kind = 4 ) TEDG(1:3,1:MXTR), triangle edge indices; see routine CVDTRI.
!
!    Output, integer ( kind = 4 ) ITR(1:2), IND(1:2), indices such that IEDG =
!    TEDG(IND(1),ITR(1)) = TEDG(IND(2),ITR(2)).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 231.
!
  integer ( kind = 4 ) mxtr

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedg
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind(2)
  integer ( kind = 4 ) itr(2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical sflag
  integer ( kind = 4 ) tedg(3,mxtr)
!
!  Search from end of array TEDG.
!
  ierror = 0
  k = 1
  j = 1
  i = mxtr

10 continue

  do

    if ( tedg(j,i) == iedg ) then
      exit
    end if

    j = j + 1

    if ( 3 < j ) then
      j = 1
      i = i - 1
      if ( i <= 0 ) then
        ierror = 231
        return
      end if
    end if

  end do

  itr(k) = i
  ind(k) = j

  if ( k == 2 ) then
    return
  end if

  k = 2

  if ( sflag ) then

    j = 1
    i = i - 1

    if ( i <= 0 ) then
      ierror = 231
      return
    end if

    go to 10

  end if
!
!  Search from beginning of array TEDG for second triangle.
!
  j = 1
  i = 1

20 continue

  if ( itr(1) <= i ) then
    ierror = 231
    return
  end if

30 continue

  if ( tedg(j,i) /= iedg ) then
    j = j + 1
    if ( 3 < j ) then
      j = 1
      i = i + 1
      go to 20
    else
      go to 30
    end if
  end if

  itr(2) = i
  ind(2) = j

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine gtime ( time )

!*****************************************************************************80
!
!! GTIME gets the current CPU time in seconds.
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) TIME, the current reading of the CPU clock
!    in seconds.
!
  integer ( kind = 4 ) clock_count
  integer ( kind = 4 ) clock_max
  integer ( kind = 4 ) clock_rate
  real ( kind = 8 ) time

  call system_clock ( clock_count, clock_rate, clock_max )

  time = real ( clock_count, kind = 8 ) / real ( clock_rate, kind = 8 )

  return
end
subroutine holvrt ( nhole, vcl, hvl, pvl, holv )

!*****************************************************************************80
!
!! HOLVRT determines top and bottom vertices of holes in polygonal regions.
!
!  Purpose: 
!
!    Determine top and bottom vertices of holes in polygonal
!    regions, and sort top vertices in decreasing (y,x) order
!    and bottom vertices in increasing (y,x) order.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NHOLE, the number of holes in region(s).
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:NHOLE), the head vertex list; HVL(I) is index in 
!    PVL of head vertex of Ith hole.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), the polygon vertex list; see routine DSPGDC.
!
!    Output, integer ( kind = 4 ) HOLV(1:NHOLE*2), the indices in PVL of top and bottom 
!    vertices of holes; first (last) NHOLE entries are for top (bottom)
!    vertices; top (bottom) vertices are sorted in decreasing
!    (increasing) lexicographic (y,x) order of coordinates.
!
  integer ( kind = 4 ) nhole

  integer ( kind = 4 ) holv(nhole*2)
  integer ( kind = 4 ) hv
  integer ( kind = 4 ) hvl(nhole)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) nhp1
  integer ( kind = 4 ) pvl(4,*)
  integer ( kind = 4 ), parameter :: succ = 3
  real ( kind = 8 ) vcl(2,*)
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine top and bottom vertices of holes.
!
  do i = 1, nhole

    hv = hvl(i)
    iv = hv

    do

      lv = pvl(loc,iv)

      if ( iv == hv ) then

        imin = iv
        imax = iv
        xmin = vcl(1,lv)
        ymin = vcl(2,lv)
        xmax = xmin
        ymax = ymin

      else

        x = vcl(1,lv)
        y = vcl(2,lv)

        if ( y < ymin .or. y == ymin .and. x < xmin ) then
          imin = iv
          xmin = x
          ymin = y
        else if ( ymax < y .or. y == ymax .and. xmax < x ) then
          imax = iv
          xmax = x
          ymax = y
        end if

      end if

      iv = pvl(succ,iv)

      if ( iv == hv) then
        exit
      end if

    end do

    holv(i) = imax
    holv(i+nhole) = imin

  end do
!
!  Use linear insertion sort to sort the top vertices of the holes
!  in decreasing (y,x) order, then bottom vertices in increasing
!  (y,x) order.  It is assumed that NHOLE is small.
!
  do i = 2, nhole

    hv = holv(i)
    lv = pvl(loc,hv)
    x = vcl(1,lv)
    y = vcl(2,lv)
    j = i

30  continue

    iv = holv(j-1)
    lv = pvl(loc,iv)

    if ( vcl(2,lv) < y .or. ( y == vcl(2,lv) .and. vcl(1,lv) < x ) ) then
      holv(j) = iv
      j = j - 1
      if ( 1 < j ) then
        go to 30
      end if
    end if

    holv(j) = hv

  end do

  nhp1 = nhole + 1

  do i = nhp1+1, nhole+nhole

    hv = holv(i)
    lv = pvl(loc,hv)
    x = vcl(1,lv)
    y = vcl(2,lv)
    j = i

50  continue

    iv = holv(j-1)
    lv = pvl(loc,iv)

    if ( y < vcl(2,lv) .or. y == vcl(2,lv) .and. x < vcl(1,lv) ) then

      holv(j) = iv
      j = j - 1

      if ( nhp1 < j ) then
        go to 50
      end if

    end if

    holv(j) = hv

  end do

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
!
!  Formula:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!  Examples:
!
!        I     J     MOD  I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integer values.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an integer to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) wide

  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i4_wrap = ilo
  else
    i4_wrap = ilo + i4_modp ( ival-ilo, wide )
  end if

  return
end
subroutine ihpsrt ( k, n, lda, a, map )

!*****************************************************************************80
!
!! IHPSRT uses heapsort on integer points in K-dimension.
!
!  Purpose:
!
!    Use heapsort to obtain the permutation of N K-dimensional
!    integer points so that the points are in lexicographic
!    increasing order.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling 
!    routine; K <= LDA.
!
!    Input, integer ( kind = 4 ) A(1:K,1:*), the array of at least N K-dimensional 
!    integer points.
!
!    Input/output, integer ( kind = 4 ) MAP(1:N), the points of A with indices MAP(1), 
!    MAP(2), ..., MAP(N) are to be sorted.
!    On output, elements are permuted so that A(*,MAP(1)) <=
!    A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) t

  do i = n/2, 1, -1
    call isftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    t = map(1)
    map(1) = map(i)
    map(i) = t
    call isftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end
function iless ( k, p, q )

!*****************************************************************************80
!
!! ILESS determines whether a K-dimensional point P is lexically less than Q.
!
!  Purpose: 
!
!    Determine whether P is lexicographically less than Q.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) P(1:K), Q(1:K), two points to compare.
!
!    Output, logical ILESS, is .TRUE. if P < Q, .FALSE. otherwise.
!
  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  logical iless
  integer ( kind = 4 ) p(k)
  integer ( kind = 4 ) q(k)

  do i = 1, k

    if ( p(i) /= q(i) ) then

      if ( p(i) < q(i) ) then
        iless = .true.
      else
        iless = .false.
      end if

      return

    end if

  end do

  iless = .false.

  return
end
subroutine i4mat_print ( lda, m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Modified:
!
!    08 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) m
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 10
    jhi = min ( jlo + 9, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,10i7)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine insed2 ( v, w, npolg, nvert, maxhv, maxpv, vcl, regnum, &
  hvl, pvl, iang, ierror )

!*****************************************************************************80
!
!! INSED2 inserts an edge into the head and polygon vertex lists.
!
!  Purpose: 
!
!    Insert the edge joining vertices V, W into the head vertex
!    list and polygon vertex list data structures.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) V, W, indices in PVL of vertices which are the endpoints
!    of an edge to be added to decomposition.
!
!    Input, integer ( kind = 4 ) NPOLG, the number of positions used in HVL array.
!
!    Input, integer ( kind = 4 ) NVERT, the number of positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL array.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL array.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT), 
!    the polygon vertex list and interior angles.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 4 or 5.
!
  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxpv

  real ( kind = 8 ) angle
  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lw
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) v
  real ( kind = 8 ) vcl(2,*)
  integer ( kind = 4 ) vv
  integer ( kind = 4 ) w
  integer ( kind = 4 ) ww

  ierror = 0

  if ( maxhv <= npolg ) then
    ierror = 4
    return
  else if ( maxpv < nvert + 2 ) then
    ierror = 5
    return
  end if
!
!  Split linked list of vertices of the polygon containing vertices
!  V and W into two linked list of vertices of polygons with common
!  edge joining V and W.
!
  nvert = nvert + 2
  vv = nvert - 1
  ww = nvert
  lv = pvl(loc,v)
  lw = pvl(loc,w)
  pvl(loc,vv) = lv
  pvl(loc,ww) = lw
  pvl(polg,ww) = pvl(polg,v)
  pvl(succ,vv) = pvl(succ,v)
  pvl(succ,ww) = pvl(succ,w)
  pvl(succ,v) = ww
  pvl(succ,w) = vv
  pvl(edgv,vv) = pvl(edgv,v)
  pvl(edgv,ww) = pvl(edgv,w)
  pvl(edgv,v) = w
  pvl(edgv,w) = v

  if ( 0 < pvl(edgv,vv) ) then
    pvl(edgv,pvl(edgv,vv)) = vv
  end if

  if ( 0 < pvl(edgv,ww) ) then
    pvl(edgv,pvl(edgv,ww)) = ww
  end if

  l = pvl(loc,pvl(succ,vv))

  iang(vv) = angle ( vcl(1,lw), vcl(2,lw), vcl(1,lv), vcl(2,lv), vcl(1,l), &
    vcl(2,l) )

  iang(v) = iang(v) - iang(vv)
  l = pvl(loc,pvl(succ,ww))
  iang(ww) = angle ( vcl(1,lv), vcl(2,lv), vcl(1,lw), vcl(2,lw), vcl(1,l), &
    vcl(2,l) )

  iang(w) = iang(w) - iang(ww)
  npolg = npolg + 1
  i = vv

  do

    pvl(polg,i) = npolg
    i = pvl(succ,i)

    if ( i == vv ) then
      exit
    end if

  end do

  hvl(pvl(polg,v)) = v
  hvl(npolg) = vv
  regnum(npolg) = regnum(pvl(polg,v))

  if ( msglvl == 2 ) then
    write ( *, '(2i6,4f15.7)' ) v, w, vcl(1:2,lv), vcl(1:2,lw)
  end if

  return
end
subroutine insvr2 ( xi, yi, wp, nvc, nvert, maxvc, maxpv, vcl, pvl, &
  iang, w, ierror )

!*****************************************************************************80
!
!! INSVR2 inserts a point into the vertex coordinate and polygon vertex lists.
!
!  Purpose: 
!
!    Insert point (XI,YI) into vertex coordinate list and
!    polygon vertex list data structures.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XI, YI, the coordinates of point to be inserted.
!
!    Input, integer ( kind = 4 ) WP, the index of vertex in PVL which is to be the
!    predecessor vertex of the inserted vertex.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of positions used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL array.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT),
!    the polygon vertex list and interior angles.
!
!    Output, integer ( kind = 4 ) W, the index of inserted vertex in PVL.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 3 or 5.
!
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc

  integer ( kind = 4 ), parameter :: edgv = 4
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ), parameter :: succ = 3
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) w
  integer ( kind = 4 ) wp
  integer ( kind = 4 ) ws
  integer ( kind = 4 ) ww
  integer ( kind = 4 ) wwp
  integer ( kind = 4 ) wws
  real ( kind = 8 ) xi
  real ( kind = 8 ) yi

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( maxvc <= nvc ) then
    ierror = 3
    return
  end if

  if ( maxpv < nvert + 2 ) then
    ierror = 5
    return
  end if
!
!  Update linked list of vertices of polygon containing vertex WP.
!
  nvc = nvc + 1
  vcl(1,nvc) = xi
  vcl(2,nvc) = yi
  ws = pvl(succ,wp)
  nvert = nvert + 1
  w = nvert
  pvl(loc,w) = nvc
  pvl(polg,w) = pvl(polg,wp)
  pvl(succ,wp) = w
  pvl(succ,w) = ws
  iang(w) = pi
  pvl(edgv,w) = pvl(edgv,wp)
!
!  If edge containing (XI,YI) is shared by another polygon,
!  then also update linked list of vertices of that polygon.
!
  if ( 0 < pvl(edgv,wp) ) then
    wws = pvl(edgv,wp)
    wwp = pvl(succ,wws)
    nvert = nvert + 1
    ww = nvert
    pvl(loc,ww) = nvc
    pvl(polg,ww) = pvl(polg,wws)
    pvl(succ,wws) = ww
    pvl(succ,ww) = wwp
    iang(ww) = pi
    pvl(edgv,wp) = ww
    pvl(edgv,ww) = wp
    pvl(edgv,wws) = w
  end if

  return
end
subroutine intpg ( nvrt, xc, yc, ctrx, ctry, arpoly, hflag, umdf, wsq, nev, &
  ifv, listev, ivrt, edgval, vrtval, vcl, mdfint, mean, stdv, mdftr )

!*****************************************************************************80
!
!! INTPG integrates the mesh distribution function in a convex polygon.
!
!  Purpose: 
!
!    Compute integral of MDF2(X,Y) [heuristic mdf] or
!    UMDF(X,Y) [user-supplied mdf] in convex polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices in polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the coordinates of 
!    polygon vertices in counter clockwise order, translated so that 
!    centroid is at origin; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, real ( kind = 8 ) CTRX, CTRY, the coordinates of centroid 
!    before translation.
!
!    Input, real ( kind = 8 ) ARPOLY, the area of polygon.
!
!    Input, logical HFLAG, is .TRUE. if heuristic mdf, .FALSE. if 
!    user-supplied mdf.
!
!    Input, external UMDF, the name of the user supplied MDF routine, of 
!    the form:
!
!      function umdf ( x, y )
!      double precision umdf
!      double precision x
!      double precision y
!
!    Input, real ( kind = 8 ) WSQ, the square of width of original polygon 
!    of decomposition.
!
!    Input, integer ( kind = 4 ) NEV, integer IFV, integer LISTEV(1:NEV), output from 
!    routine PRMDF2.
!
!    Input, IVRT(1:*), real ( kind = 8 ) EDGVAL(1:*), 
!    double precision VRTVAL(1:*), arrays output from DSMDF2;
!    if .NOT. HFLAG then only first array exists.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Output, real ( kind = 8 ) MDFINT, the integral of MDF in polygon.
!
!    Output, real ( kind = 8 ) MEAN, the mean MDF value in polygon.
!
!    Output, real ( kind = 8 ) STDV, the standard deviation of MDF in polygon.
!
!    Output, real ( kind = 8 ) MDFTR(0:NVRT-1), the mean MDF value in each 
!    triangle of polygon;  triangles are determined by polygon vertices 
!    and centroid.
!
  integer ( kind = 4 ) nev
  integer ( kind = 4 ), parameter :: nqpt = 3
  integer ( kind = 4 ) nvrt

  real ( kind = 8 ) areatr
  real ( kind = 8 ) arpoly
  real ( kind = 8 ) ctrx
  real ( kind = 8 ) ctry
  real ( kind = 8 ) d
  real ( kind = 8 ) edgval(*)
  logical hflag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifv
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) listev(nev)
  integer ( kind = 4 ) m
  real ( kind = 8 ) mdfint
  real ( kind = 8 ) mdfsqi
  real ( kind = 8 ) mdftr(0:nvrt-1)
  real ( kind = 8 ) mean
  real ( kind = 8 ), save, dimension ( 3, nqpt ) :: qc = reshape ( (/ &
    0.6666666666666666D+00, 0.1666666666666667D+00, 0.1666666666666667D+00, &
    0.1666666666666667D+00, 0.6666666666666666D+00, 0.1666666666666667D+00, &
    0.1666666666666667D+00, 0.1666666666666667D+00, 0.6666666666666666D+00/), &
    (/ 3, nqpt /) )
  real ( kind = 8 ) s
  real ( kind = 8 ) stdv
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) umdf
  real ( kind = 8 ) val
  real ( kind = 8 ) vcl(2,*)
  real ( kind = 8 ) vrtval(*)
  real ( kind = 8 ) wsq
  real ( kind = 8 ), save, dimension ( nqpt ) :: wt = &
    (/ 0.3333333333333333D+00, 0.3333333333333333D+00, 0.3333333333333333D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) xx
  real ( kind = 8 ) y
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) yc(0:nvrt)
  real ( kind = 8 ) yy

  external umdf
!
!  NQPT is number of quadrature points for numerical integration in triangle.
!  WT(I) is weight of Ith quadrature point.
!  QC(1:3,I) are barycentric coordinates of Ith quadrature point.
!
  mdfint = 0.0D+00
  mdfsqi = 0.0D+00

  do l = 0, nvrt-1

    areatr = 0.5D+00 * ( xc(l) * yc(l+1) - xc(l+1) * yc(l) )
    sum1 = 0.0D+00
    sum2 = 0.0D+00

    do m = 1, nqpt

      xx = qc(1,m) * xc(l) + qc(2,m) * xc(l+1)
      yy = qc(1,m) * yc(l) + qc(2,m) * yc(l+1)
      if ( hflag ) then
!
!             VAL = MDF2(XX+CTRX,YY+CTRY,WSQ,NEV,IFV,LISTEV,IVRT, &
!               EDGVAL,VRTVAL,VCL)
!
!  Insert code for function MDF2 to reduce number of calls.
!
         x = xx + ctrx
         y = yy + ctry
           s = wsq
           do i = 1, nev
            k = listev(i)
            if ( k < 0 ) then
               k = -k
               d = ( vcl(1,k) - x )**2 + ( vcl(2,k) - y )**2
               d = max ( 0.25D+00 * d, vrtval(k) )
               s = min ( s, d )
            else
               kp1 = k + 1
               if ( i == nev .and. 0 < ifv ) then
                 kp1 = ifv
               end if
               j = ivrt(kp1)
               x0 = x - vcl(1,j)
               y0 = y - vcl(2,j)
               x1 = vcl(1,ivrt(k)) - vcl(1,j)
               y1 = vcl(2,ivrt(k)) - vcl(2,j)

               if ( x0 * x1 + y0 * y1 <= 0.0D+00 ) then
                 d = x0**2 + y0**2
               else
                 x0 = x0 - x1
                 y0 = y0 - y1
                 if ( 0.0D+00 <= x0 * x1 + y0 * y1 ) then
                   d = x0**2 + y0**2
                 else
                   d = ( x1 * y0 - y1 * x0 )**2 / ( x1**2 + y1**2 )
                 end if
               end if

               d = max ( 0.25D+00 * d, edgval(k) )
               s = min ( s, d )
            end if
           end do
           val = 1.0D+00 / s
      else
        val = umdf ( xx+ctrx, yy+ctry )
      end if

      temp = wt(m) * val
      sum1 = sum1 + temp
      sum2 = sum2 + temp * val

    end do

    mdftr(l) = sum1
    mdfint = mdfint + sum1 * areatr
    mdfsqi = mdfsqi + sum2 * areatr

  end do

  mean = mdfint / arpoly
  stdv = mdfsqi / arpoly - mean**2
  stdv = sqrt ( max ( stdv, 0.0D+00 ) )

  return
end
subroutine inttri ( nvrt, xc, yc, h, ibot, costh, sinth, ldv, nvc, ntri,  &
  maxvc, maxti, maxcw, vcl, til, ncw, cwalk, ierror )

!*****************************************************************************80
!
!! INTTRI generates triangles inside a convex polygon.
!
!  Purpose: 
!
!    Generate triangles inside convex polygon using quasi-uniform grid of
!    spacing H.  It is assumed that the diameter of the polygon is parallel 
!    to the Y axis.
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of
!    convex polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex coordinates 
!    in counter clockwise order; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, real ( kind = 8 ) H, the spacing of mesh vertices in polygon.
!
!    Input, integer ( kind = 4 ) IBOT, the index of bottom vertex; diameter contains vertices
!    (XC(0),YC(0)) and (XC(IBOT),YC(IBOT)).
!
!    Input, real ( kind = 8 ) COSTH, SINTH; COS(THETA), SIN(THETA) where 
!    THETA in [-PI,PI] is rotation angle to get diameter parallel to y-axis.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of coordinates or positions 
!    used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NTRI, the number of triangles or positions 
!    used in TIL.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXTI, the maximum size available for TIL array.
!
!    Input, integer ( kind = 4 ) MAXCW, the maximum size available for CWALK array; 
!    assumed to be at least 6*(1 + INT((YC(0) - YC(IBOT))/H)).
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list.
!
!    Output, integer ( kind = 4 ) NCW, the number of mesh vertices in closed walk, 
!    except NCW = 0 for 1 vertex.
!
!    Output, integer ( kind = 4 ) CWALK(0:NCW), indices in VCL of mesh vertices of closed
!    walk; CWALK(0) = CWALK(NCW)
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 9, or 10
!
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) maxcw
  integer ( kind = 4 ) maxti
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) nvrt

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) costh
  integer ( kind = 4 ) cwalk(0:maxcw)
  real ( kind = 8 ) cy
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) il
  integer ( kind = 4 ) im1l
  integer ( kind = 4 ) im1r
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lw
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncw
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) p
  integer ( kind = 4 ) r
  integer ( kind = 4 ) r0
  integer ( kind = 4 ) r1
  integer ( kind = 4 ) rw
  real ( kind = 8 ) sinth
  real ( kind = 8 ) sy
  integer ( kind = 4 ) til(3,maxti)
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(ldv,maxvc)
  real ( kind = 8 ) x
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) xj
  real ( kind = 8 ) xk
  real ( kind = 8 ) xl
  real ( kind = 8 ) xm1l
  real ( kind = 8 ) xm1r
  real ( kind = 8 ) xr
  real ( kind = 8 ) y
  real ( kind = 8 ) yc(0:nvrt)

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  n = int ( ( yc(0) - yc(ibot) ) / h )
  y = yc(0) - 0.5D+00 * ( yc(0) - yc(ibot ) - real ( n, kind = 8 ) * h )
  l = 0
  r = nvrt

  do i = 0, n
!
!  Determine left and right x-coordinates of polygon for
!  scan line with y-coordinate Y, and generate mesh vertices.
!
    do while ( y < yc(l+1) )
      l = l + 1
    end do

    do while ( y < yc(r-1) )
      r = r - 1
    end do

    xl = xc(l) + ( xc(l+1) - xc(l) ) * ( y - yc(l) ) / ( yc(l+1) - yc(l) )
    xr = xc(r) + ( xc(r-1) - xc(r) ) * ( y - yc(r) ) / ( yc(r-1) - yc(r) )
    m = int ( ( xr - xl ) / h )
    x = xl + 0.5D+00 * ( xr - xl - real ( m, kind = 8 ) * h )

    if ( maxvc < nvc + m + 1 ) then
      ierror = 3
      return
    end if

    cy = costh * y
    sy = sinth * y
    il = nvc + 1
    xl = x

    do j = 0, m
      nvc = nvc + 1
      vcl(1,nvc) = costh * x + sy
      vcl(2,nvc) = cy - sinth * x
      x = x + h
    end do

    ir = nvc
    xr = x - h

    if ( n == 0 ) then

      ncw = 0
      cwalk(0) = nvc
      return

    else if ( i == 0 ) then

      lw = 0
      cwalk(lw) = il
      rw = maxcw + 1

      do j = il, ir
        rw = rw - 1
        cwalk(rw) = j
      end do

      go to 100

    end if
!
!  Generate triangles between scan lines Y+H and Y.
!
    a = max ( xl, xm1l )
    b = min ( xr, xm1r )

    if ( xm1l == a ) then
      l0 = im1l
      x = ( xm1l - xl ) / h
      j = int(x + tol)
      if ( abs ( x - real ( j, kind = 8 ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      l1 = il + j
    else
      l1 = il
      x = ( xl - xm1l ) / h
      j = int ( x + tol )
      if ( abs ( x - real ( j, kind = 8 ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      l0 = im1l + j
    end if

    if ( xm1r == b ) then
      r0 = im1r
      x = ( xr - xm1r ) / h
      j = int ( x + tol )
      if ( abs ( x - real ( j, kind = 8 ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      r1 = ir - j
    else
      r1 = ir
      x = ( xm1r - xr ) / h
      j = int ( x + tol )
      if ( abs ( x - real ( j, kind = 8 ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      r0 = im1r - j
    end if

    if ( l0 < r0 .or. l1 < r1 ) then

      j = l0
      k = l1
      xj = xm1l + real ( j-im1l, kind = 8 ) * h
      xk = xl + real ( k - il, kind = 8 ) * h

      do

        if ( k < r1 .and. ( xk <= xj .or. j == r0 ) ) then
          p = k
          k = k + 1
          xk = xk + h
        else
          p = j
          j = j + 1
          xj = xj + h
        end if

        ntri = ntri + 1

        if ( maxti < ntri ) then
          ierror = 9
          return
        end if

        til(1,ntri) = j
        til(2,ntri) = p
        til(3,ntri) = k

        if ( r0 <= j .and. r1 <= k ) then
          exit
        end if

      end do

    end if
!
!  Generate paths of closed walk between scan lines Y+H and Y.
!
    if ( xm1l < xl ) then
      do j = im1l+1, l0
        lw = lw + 1
        cwalk(lw) = j
      end do
      lw = lw + 1
      cwalk(lw) = il
    else
      do j = l1, il, -1
        lw = lw + 1
        cwalk(lw) = j
      end do
    end if

    if ( xr < xm1r ) then
      do j = im1r-1, r0, -1
        rw = rw - 1
        cwalk(rw) = j
      end do
      rw = rw - 1
      cwalk(rw) = ir
    else
      do j = r1, ir
        rw = rw - 1
        cwalk(rw) = j
      end do
    end if

100 continue

    y = y - h
    im1l = il
    im1r = ir
    xm1l = xl
    xm1r = xr

  end do
!
!  Add last path of left walk and shift indices of right walk.
!
  if ( m == 0 ) then
    rw = rw + 1
  else
    do j = il+1, ir-1
      lw = lw + 1
      cwalk(lw) = j
    end do
  end if

  if ( rw <= lw ) then
    ierror = 10
    return
  end if

  do j = rw, maxcw
    lw = lw + 1
    cwalk(lw) = cwalk(j)
  end do

  ncw = lw

  return
end
subroutine isftdw ( l, u, k, lda, a, map )

!*****************************************************************************80
!
!! ISFTDW sifts A(*,MAP(L)) down a heap of size U.
!
!
!  Purpose: 
!
!    Sift A(*,MAP(L)) down a heap of size U.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, U, the lower and upper index of part of heap.
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, integer ( kind = 4 ) A(1:K,1:*); see routine IHPSRT.
!
!    Input/output, integer ( kind = 4 ) MAP(1:*); see routine IHPSRT.
!
  integer ( kind = 4 ) lda

  integer ( kind = 4 ) a(lda,*)
  integer ( kind = 4 ) i
  logical iless
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) map(*)
  integer ( kind = 4 ) u
  integer ( kind = 4 ) t

  i = l
  j = 2 * i
  t = map(i)

  do while ( j <= u )

    if ( j < u ) then
      if ( iless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( iless ( k, a(1,map(j)), a(1,t) ) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Modified:
!
!    09 November 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine jnhole ( itophv, angspc, angtol, nvc, nvert, maxvc, maxpv,  &
  maxiw, maxwk, vcl, hvl, pvl, iang, iwk, wk, ierror )

!*****************************************************************************80
!
!! JNHOLE joins a hole boundary to the boundary of a polygon.
!
!  Purpose: 
!
!    Join hole boundary to boundary of polygon containing hole
!    by finding a cut edge originating from the top vertex of hole
!    to a point on outer polygon boundary above it.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITOPHV, the index in PVL of top vertex of hole.
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter used in 
!    controlling the vertices to be considered as an endpoint of a separator.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter used 
!    in accepting separator(s).
!
!    Input/output, integer ( kind = 4 ) NVC, the number of positions used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should 
!    be about 3 times number of vertices in outer polygon.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should 
!    be about 5 times number of vertices in outer polygon.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC),the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:*), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT), 
!    the polygon vertex list and interior angles.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 5, 6, 7, 206 to 210, 212, 218, or 219.
!
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real ( kind = 8 ) angle
  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  real ( kind = 8 ) dy
  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ) hv
  integer ( kind = 4 ) hvl(*)
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilft
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) irgt
  integer ( kind = 4 ) itophv
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) ivs
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lw
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  real ( kind = 8 ) s
  real ( kind = 8 ) slft
  real ( kind = 8 ) srgt
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) succil
  integer ( kind = 4 ) succir
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vp
  integer ( kind = 4 ) vr
  integer ( kind = 4 ) vs
  integer ( kind = 4 ) vv
  integer ( kind = 4 ) w
  real ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) ww
  real ( kind = 8 ) xint
  real ( kind = 8 ) xlft
  real ( kind = 8 ) xrgt
  real ( kind = 8 ) xt
  real ( kind = 8 ) xv
  real ( kind = 8 ) xvs
  real ( kind = 8 ) ylft
  real ( kind = 8 ) yrgt
  real ( kind = 8 ) yt
  real ( kind = 8 ) ytmtol
  real ( kind = 8 ) ytptol
  real ( kind = 8 ) yv
  real ( kind = 8 ) yvs
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( maxvc < nvc + 3 ) then
    ierror = 3
    return
  end if

  if ( maxpv < nvert+5 ) then
    ierror = 5
    return
  end if
!
!  Determine 'closest' vertices on outer boundary which are to the
!  left and right of top vertex of hole and on the horizontal line
!  through top vertex. The two closest vertices must be on edges
!  which intersect the horizontal line and are partially above the
!  line. Ties are broken (in the case of a vertex on a cut edge)
!  by choosing the vertex on the edge of maximum or minimum dx/dy
!  slope depending on whether the vertex is to the left or right
!  of top vertex, respectively.
!
  ipoly = pvl(polg,itophv)
  lv = pvl(loc,itophv)
  xt = vcl(1,lv)
  yt = vcl(2,lv)
  dy = 0.0D+00
  hv = hvl(ipoly)
  iv = hv
  yv = vcl(2,pvl(loc,iv))

  do

    iv = pvl(succ,iv)
    yvs = vcl(2,pvl(loc,iv))
    dy = max ( dy, abs ( yvs - yv ) )
    yv = yvs

    if ( iv == hv ) then
      exit
    end if

  end do

  ytmtol = yt - tol * dy
  ytptol = yt + tol * dy
  ilft = 0
  irgt = 0
  xlft = 0.0D+00
  xrgt = 0.0D+00
  hv = hvl(ipoly)
  iv = hv
  lv = pvl(loc,iv)
  xv = vcl(1,lv)
  yv = vcl(2,lv)

20 continue

  ivs = pvl(succ,iv)
  lv = pvl(loc,ivs)
  xvs = vcl(1,lv)
  yvs = vcl(2,lv)

  if ( yv <= ytptol .and. ytptol < yvs ) then

    if ( ytmtol <= yv ) then
      if ( xt < xv ) then
        if ( xv < xrgt .or. irgt == 0 ) then
          irgt = iv
          xrgt = xv
          yrgt = yv
          srgt = (xvs - xv ) / ( yvs - yv )
        else if ( xv == xrgt ) then
          s = ( xvs - xv ) / ( yvs - yv )
          if ( s < srgt ) then
            irgt = iv
            yrgt = yv
            srgt = s
          end if
        end if
      end if
    else
      xint = ( yt - yv ) * ( xvs - xv ) / ( yvs - yv ) + xv
      if ( xt < xint ) then
        if ( xint < xrgt .or. irgt == 0 ) then
          irgt = iv
          xrgt = xint
          yrgt = yt
        end if
      end if
    end if

  else if ( ytptol < yv .and. yvs <= ytptol ) then

    if ( ytmtol <= yvs ) then
      if ( xvs < xt ) then
        if ( xlft < xvs .or. ilft == 0 ) then
          ilft = iv
          xlft = xvs
          ylft = yvs
          slft = ( xvs - xv ) / ( yvs - yv )
        else if ( xvs == xlft ) then
          s = ( xvs - xv ) / ( yvs - yv )
          if ( slft < s ) then
            ilft = iv
            ylft = yvs
            slft = s
          end if
        end if
      end if
    else
      xint = ( yt - yv ) * ( xvs - xv ) / ( yvs - yv ) + xv
      if ( xint < xt ) then
        if ( xlft < xint .or. ilft == 0 ) then
          ilft = iv
          xlft = xint
          ylft = yt
        end if
      end if
    end if
  end if

  iv = ivs
  xv = xvs
  yv = yvs

  if ( iv /= hv ) then
    go to 20
  end if

  if ( ilft == 0 .or. irgt  ==  0 ) then
   ierror = 218
   return
  end if
!
!  Temporarily modify PVL to pass the subregion 'above' top vertex
!  of hole to routine RESVRT. The top vertex is the reflex vertex
!  passed to RESVRT (in the temporary subregion, it has interior
!  angle PI). This causes one separator to be chosen by RESVRT
!  and its other endpoint is above the top vertex.
!
  succil = pvl(succ,ilft)
  succir = pvl(succ,irgt)
  vcl(1,nvc+2) = xlft
  vcl(2,nvc+2) = ylft
  vcl(1,nvc+3) = xrgt
  vcl(2,nvc+3) = yrgt
  vp = nvert + 3
  vr = nvert + 4
  vs = nvert + 5

  iang(vr) = angle ( xlft, ylft, xt, yt, xrgt, yrgt )

  if ( iang(vr) < pi - tol .or. pi + tol < iang(vr) ) then
    ierror = 219
    return
  end if

  pvl(loc,vp) = nvc + 2
  pvl(polg,vp) = ipoly
  pvl(succ,vp) = vr
  pvl(edgv,vp) = 0
  pvl(loc,vr) = pvl(loc,itophv)
  pvl(polg,vr) = ipoly
  pvl(succ,vr) = vs
  pvl(edgv,vr) = 0
  pvl(loc,vs) = nvc + 3
  pvl(polg,vs) = ipoly
  pvl(succ,vs) = succir
  pvl(edgv,vs) = pvl(edgv,irgt)
  pvl(succ,ilft) = vp
  lv = pvl(loc,ilft)
  iang(vp) = angle ( vcl(1,lv), vcl(2,lv), xlft, ylft, xt, yt )
  lv = pvl(loc,succir)
  iang(vs) = angle ( xt, yt, xrgt, yrgt, vcl(1,lv), vcl(2,lv) )
  w = 0

  call resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, maxwk, &
    vcl, pvl, iang, w, ww, iwk, wk, ierror )
!
!  Remove temporary modification to PVL. There are three cases
!  depending on where the endpoint of separator is located:
!  successor of closest vertex to the right of top vertex,
!  predecessor of closest vertex to the left of top vertex,
!  or neither of these.
!
  if ( pvl(succ,vs) == w ) then
    pvl(succ,ilft) = succil
    pvl(succ,irgt) = w
    pvl(edgv,irgt) = pvl(edgv,vs)
    if ( 0 < pvl(edgv,irgt) ) then
      pvl(edgv,pvl(edgv,irgt)) = irgt
    end if
  else if ( pvl(succ,ilft) == w ) then
    pvl(succ,w) = succil
  else
    pvl(succ,ilft) = succil
  end if

  if ( ierror /= 0 ) then
    return
  end if
!
!  Update PVL with cut edge, i.e. join linked lists of vertices
!  of the hole polygon and the outer boundary polygon into one
!  linked list of vertices by adding the cut edge from the top
!  vertex of hole to the vertex on the outer boundary.
!
  nvert = nvert + 2
  vv = nvert - 1
  ww = nvert
  lv = pvl(loc,itophv)
  lw = pvl(loc,w)
  pvl(loc,vv) = lv
  pvl(loc,ww) = lw
  pvl(polg,vv) = ipoly
  pvl(polg,ww) = ipoly
  pvl(succ,vv) = pvl(succ,itophv)
  pvl(succ,ww) = pvl(succ,w)
  pvl(succ,itophv) = ww
  pvl(succ,w) = vv
  pvl(edgv,vv) = pvl(edgv,itophv)
  pvl(edgv,ww) = pvl(edgv,w)
  pvl(edgv,itophv) = w
  pvl(edgv,w) = itophv

  if ( 0 < pvl(edgv,vv) ) then
    pvl(edgv,pvl(edgv,vv)) = vv
  end if

  if ( 0 < pvl(edgv,ww) ) then
    pvl(edgv,pvl(edgv,ww)) = ww
  end if

  l = pvl(loc,pvl(succ,vv))
  iang(vv) = angle ( vcl(1,lw), vcl(2,lw), vcl(1,lv), vcl(2,lv), vcl(1,l), &
    vcl(2,l) )
  iang(itophv) = iang(itophv) - iang(vv)
  l = pvl(loc,pvl(succ,ww))
  iang(ww) = angle ( vcl(1,lv), vcl(2,lv), vcl(1,lw), vcl(2,lw), vcl(1,l), &
    vcl(2,l) )
  iang(w) = iang(w) - iang(ww)

  if ( msglvl == 2 ) then
    write ( *,600) itophv,w,vcl(1,lv),vcl(2,lv),vcl(1,lw),vcl(2,lw)
  end if

600 format (1x,2i7,4f15.7)

  return
end
subroutine lop ( itr, ind, mxedg, top, ldv, vcl, til, tedg, sptr )

!*****************************************************************************80
!
!! LOP applies the local optimization procedure to two triangles.
!
!  Purpose: 
!
!    Apply local optimization procedure to two triangles
!    indicated by ITR(1) and ITR(2). This may result in swapping
!    diagonal edge of quadrilateral.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITR(1:2), the indices of triangles for LOP.
!
!    Input, integer ( kind = 4 ) IND(1:2), indices indicating common edge of triangles.
!
!    Input, integer ( kind = 4 ) MXEDG, the maximum index of edge to be considered for LOP.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of SPTR indicating top of stack.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input/output, integer ( kind = 4 ) TEDG(1:3,1:*), the triangle edge indices; 
!    see routine CVDTRI.
!
!    Input/output, integer ( kind = 4 ) SPTR(1:*), stack pointers; see routine CVDTRI.
!
  integer ( kind = 4 ) ldv

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) iedg
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ind(2)
  integer ( kind = 4 ) ind1m1
  integer ( kind = 4 ) ind1p1
  integer ( kind = 4 ) ind2m1
  integer ( kind = 4 ) ind2p1
  integer ( kind = 4 ) itr(2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mxedg
  integer ( kind = 4 ) top
  integer ( kind = 4 ) sptr(*)
  integer ( kind = 4 ) tedg(3,*)
  integer ( kind = 4 ) til(3,*)
  real ( kind = 8 ) vcl(ldv,*)
!
!  Common edge is BC, other two vertices are A and D.
!
  iedg = tedg(ind(1),itr(1))
  sptr(iedg) = -1

  ind1m1 = i4_wrap ( ind(1) - 1, 1, 3 )
  ind1p1 = i4_wrap ( ind(1) + 1, 1, 3 )
  ind2m1 = i4_wrap ( ind(2) - 1, 1, 3 )
  ind2p1 = i4_wrap ( ind(2) + 1, 1, 3 )

  b = til(ind(1),itr(1))
  c = til(ind1p1,itr(1))
  a = til(ind1m1,itr(1))
  d = til(ind2m1,itr(2))

  in = diaedg ( vcl(1,d), vcl(2,d), vcl(1,c), vcl(2,c), vcl(1,a), vcl(2,a), &
    vcl(1,b), vcl(2,b) )

  if ( in == 1 ) then
!
!  Check if four edges of quadrilateral should be put on LOP
!  stack, and swap edge BC for AD.
!
   i = tedg(ind1m1,itr(1))

   do j = 1, 4

      if ( j == 2 ) then
        i = tedg(ind1p1,itr(1))
      else if ( j == 3 ) then
        i = tedg(ind2m1,itr(2))
      else if ( j == 4 ) then
        i = tedg(ind2p1,itr(2))
      end if

      if ( i <= mxedg ) then
        if ( sptr(i) == -1 ) then
          sptr(i) = top
          top = i
        end if
      end if

    end do

    til(ind1p1,itr(1)) = d
    til(ind2p1,itr(2)) = a
    tedg(ind(1),itr(1)) = tedg(ind2p1,itr(2))
    tedg(ind(2),itr(2)) = tedg(ind1p1,itr(1))
    tedg(ind1p1,itr(1)) = iedg
    tedg(ind2p1,itr(2)) = iedg

  end if

  return
end
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines if a point is left of, right or, or on a directed line.
!
!  Discussion:
!
!    The directed line is parallel to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, real ( kind = 8 ) DV, the signed distance of the directed line 
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).  
!    DV is positive for a line to the left of the base line.
! 
!    Output, integer ( kind = 4 ) LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
  implicit none

  real ( kind = 8 ) dv
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxu
  real ( kind = 8 ) dy
  real ( kind = 8 ) dyu
  integer ( kind = 4 ) lrline
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2

  tol = 100.0D+00 * epsilon ( tol )

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
    abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else
    lrline = -1
  end if

  return
end
subroutine lufac ( a, lda, n, tol, ipvt, singlr )

!*****************************************************************************80
!
!! LUFAC computes the LU factorization of a matrix.
!
!  Purpose: 
!
!    Obtain LU factorization of matrix A, i.e. apply Gaussian
!    elimination with partial pivoting to A.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(1:N,1:N), the matrix.  On input, 
!    the N by N matrix to be factored.  On output, the upper triangular 
!    matrix U and multipliers of unit lower triangular matrix L (if matrix 
!    A is nonsingular).
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, integer ( kind = 4 ) N, the order of matrix A.
!
!    Input, double preicison TOL, the relative tolerance for detecting 
!    singularity of A.
!
!    Output, integer ( kind = 4 ) IPVT(1:N-1), the pivot indices.
!
!    Output, logical SINGLR, is .TRUE. if matrix is singular; this occurs 
!    when the magnitude of a pivot element is <= TOL*MAX(|A(I,J)|).
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) m
  logical singlr
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs

  if ( n < 1 ) then
    return
  end if

  singlr = .true.

  t = 0.0D+00
  do j = 1, n
    do i = 1, n
      t = max ( t, abs ( a(i,j) ) )
    end do
  end do

  tolabs = tol * t

  do k = 1, n-1

    kp1 = k + 1
    m = k

    do i = k+1, n
      if ( abs ( a(m,k) ) < abs ( a(i,k) ) ) then
        m = i
      end if
    end do

    ipvt(k) = m
    t = a(m,k)
    a(m,k) = a(k,k)
    a(k,k) = t

    if ( abs ( t) <= tolabs ) then
      return
    end if

    do i = k+1, n
      a(i,k) = a(i,k) / t
    end do

    do j = k+1, n

      t = a(m,j)
      a(m,j) = a(k,j)
      a(k,j) = t

      if ( t /= 0.0E+00 ) then

        do i = k+1, n
          a(i,j) = a(i,j) - a(i,k) * t
        end do

      end if

    end do

  end do

  if ( tolabs < abs ( a(n,n) ) ) then
    singlr = .false.
  end if

  return
end
subroutine lusol ( a, lda, n, ipvt, b )

!*****************************************************************************80
!
!! LUSOL solves a linear system with an LU factored matrix.
!
!  Purpose: 
!
!    Solve linear system A*X = B given LU factorization of A.
!    It is assumed that A is nonsingular.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(1:N,1:N), contains factors L, U output 
!    from routine LUFAC.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, integer ( kind = 4 ) N, the order of matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(1:N-1), the pivot indices from routine LUFAC.
!
!    Input/output, B(1:N).  On input, the right hand side vector.
!    On output, the solution vector X.
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) ipvt(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) t
!
!  Forward elimination
!
  do k = 1, n-1

    m = ipvt(k)
    t = b(m)
    b(m) = b(k)
    b(k) = t

    do i = k+1, n
      b(i) = b(i) - a(i,k) * t
    end do

  end do
!
!  Back substitution
!
  do k = n, 2, -1

    t = b(k) / a(k,k)
    b(k) = t

    do i = 1, k-1
      b(i) = b(i) - a(i,k) * t
    end do

  end do

  b(1) = b(1) / a(1,1)

  return
end
function mdf2 ( x, y, wsq, nev, ifv, listev, ivrt, edgval, vrtval, vcl )

!*****************************************************************************80
!
!! MDF2 evaluates the heuristic mesh distribution function at (X,Y).
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of point.
!
!    Input, real ( kind = 8 ) WSQ, the square of the width of the polygon 
!    containing (X,Y).
!
!    Input, integer ( kind = 4 ) NEV, IFV, LISTEV(1:NEV), output from routine PRMDF2.
!
!    Input, integer ( kind = 4 ) IVRT(1:*), output from DSMDF2.
!
!    Input, real ( kind = 8 ) EDGVAL(1:*), output from DSMDF2.
!
!    Input, real ( kind = 8 ) VRTVAL(1:*), output from DSMDF2.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Output, real ( kind = 8 ) MDF2, the reciprocal of square of length 
!    scale at (X,Y).
!
  integer ( kind = 4 ) nev

  real ( kind = 8 ) d
  real ( kind = 8 ) edgval(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifv
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) listev(nev)
  real ( kind = 8 ) mdf2
  real ( kind = 8 ) s
  real ( kind = 8 ) vcl(2,*)
  real ( kind = 8 ) vrtval(*)
  real ( kind = 8 ) wsq
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) y
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1

  s = wsq

  do i = 1, nev

    k = listev(i)

    if ( k < 0 ) then
      k = -k
      d = ( vcl(1,k) - x )**2 + ( vcl(2,k) - y )**2
      d = max ( 0.25D+00 * d, vrtval(k) )
      s = min ( s, d )
    else
      kp1 = k + 1
      if ( i == nev .and. 0 < ifv ) then
        kp1 = ifv
      end if
      j = ivrt(kp1)
      x0 = x - vcl(1,j)
      y0 = y - vcl(2,j)
      x1 = vcl(1,ivrt(k)) - vcl(1,j)
      y1 = vcl(2,ivrt(k)) - vcl(2,j)

      if ( x0 * x1 + y0 * y1 <= 0.0D+00 ) then
        d = x0**2 + y0**2
      else
        x0 = x0 - x1
        y0 = y0 - y1
        if ( 0.0D+00 <= x0 * x1 + y0 * y1 ) then
          d = x0**2 + y0**2
        else
          d = ( x1 * y0 - y1 * x0 )**2 / ( x1**2 + y1**2 )
        end if
      end if

      d = max ( 0.25D+00 * d, edgval(k) )
      s = min ( s, d )

    end if

  end do

  mdf2 = 1.0D+00 / s

  return
end
subroutine mfdec2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid,  &
  nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl,  &
  pvl, iang, ivrt, xivrt, widsq, edgval, vrtval, area, psi, iwk, wk, ierror )

!*****************************************************************************80
!
!! MFDEC2 subdivides polygons to decrease mesh distribution variation.
!
!  Purpose: 
!
!    Further subdivide convex polygons so that the variation
!    of heuristic or user-supplied mesh distribution function in
!    each polygon is limited.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical HFLAG, is .TRUE. if heuristic mdf, .FALSE. if 
!    user-supplied mdf.
!
!    Input, external UMDF, a user-supplied mdf, of the form
!
!      function umdf ( x, y )
!      double precision umdf
!      double precision x
!      double precision y
!
!    Input, real ( kind = 8 ) KAPPA, the mesh smoothness parameter in 
!    interval [0.0,1.0].
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in radians 
!    used to determine extra points as possible endpoints of separators.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter in 
!    radians used in accepting separators.
!
!    Input, real ( kind = 8 ) DMIN, a parameter used to determine if 
!    variation of mdf in polygon is 'sufficiently high'.
!
!    Input, integer ( kind = 4 ) NMIN, a parameter used to determine if 'sufficiently large'
!    number of triangles in polygon.
!
!    Input, integer ( kind = 4 ) NTRID, the desired number of triangles in mesh.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or 
!    positions used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions or 
!    positions used in HVL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of polygon vertices or positions 
!    used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM, AREA, 
!    PSI arrays.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should 
!    be about 3*NVRT + INT(2*PI/ANGSPC) where NVRT is maximum number of
!    vertices in a convex polygon of the (input) decomposition.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should 
!    be about NPOLG + 3*(NVRT + INT(2*PI/ANGSPC)) + 2.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, HVL(1:NPOLG), the head vertex list.
!
!    Input/output, PVL(1:4,1:NVERT), IANG(1:NVERT), the polygon vertex list and
!    interior angles.
!
!    Input, integer ( kind = 4 ) IVRT(1:NVERT), integer XIVRT(1:NPOLG+1), 
!    double precision WIDSQ(1:NPOLG), real ( kind = 8 ) EDGVAL(1:NVERT),
!    double precision VRTVAL(1:NVC), arrays output from routine DSMDF2;
!    if .NOT. HFLAG then only first two arrays exist.
!
!    Input/output, real ( kind = 8 ) AREA(1:NPOLG), the area of convex polygons 
!    in decomposition.
!
!    Output, real ( kind = 8 ) PSI(1:NPOLG), the mean mdf values in the
!    convex polygons
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 4, 5, 6, 7, 200, or 222.
!
  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert

  real ( kind = 8 ) alpha
  real ( kind = 8 ) angsp2
  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  real ( kind = 8 ) area(maxhv)
  real ( kind = 8 ) areapg
  real ( kind = 8 ) arearg
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) cosalp
  real ( kind = 8 ) ctrx
  real ( kind = 8 ) ctry
  real ( kind = 8 ) delta
  real ( kind = 8 ) dmin
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) edgval(nvert)
  logical hflag
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifv
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) indpvl
  real ( kind = 8 ) intreg
  integer ( kind = 4 ) ivrt(nvert)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) kappa
  integer ( kind = 4 ) l
  integer ( kind = 4 ) listev
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) maxn
  real ( kind = 8 ) mdfint
  integer ( kind = 4 ) mdftr
  real ( kind = 8 ) mean
  integer ( kind = 4 ) nev
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) ntrid
  real ( kind = 8 ) numer
  integer ( kind = 4 ) nvrt
  real ( kind = 8 ) nwarea
  integer ( kind = 4 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pi2
  real ( kind = 8 ) psi(maxhv)
  integer ( kind = 4 ) pvl(4,maxpv)
  real ( kind = 8 ) r
  integer ( kind = 4 ) regnum(maxhv)
  real ( kind = 8 ) sinalp
  real ( kind = 8 ) stdv
  integer ( kind = 4 ), parameter :: succ = 3
  real ( kind = 8 ) sumx
  real ( kind = 8 ) sumy
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) tol
  real ( kind = 8 ) umdf
  integer ( kind = 4 ) v
  real ( kind = 8 ) vcl(2,maxvc)
  real ( kind = 8 ) vrtval(nvc)
  integer ( kind = 4 ) w
  real ( kind = 8 ) widsq(npolg)
  real ( kind = 8 ) wk(maxwk)
  real ( kind = 8 ) wsq
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  integer ( kind = 4 ) xc
  integer ( kind = 4 ) xivrt(npolg+1)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  integer ( kind = 4 ) yc

  external umdf

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  WK(1:NPOLG) is used for MDF standard deviation in polygons.
!  Compute AREARG = area of region and INTREG = estimated integral
!  of MDF2(X,Y) or UMDF(X,Y).
!
  nvrt = 0
  do i = 1, npolg
    nvrt = max ( nvrt, xivrt(i+1) - xivrt(i) )
  end do

  if ( hflag .and. maxiw < 2 * nvrt ) then
    ierror = 6
    return
  end if

  if ( maxwk < npolg + 3 * nvrt + 2 ) then
    ierror = 7
    return
  end if

  listev = 1
  xc = npolg + 1
  yc = xc + nvrt + 1
  mdftr = yc + nvrt + 1
  arearg = 0.0D+00
  intreg = 0.0D+00
  nev = -1

  do i = 1, npolg

    if ( hflag ) then
      wsq = widsq(i)
      call prmdf2 ( i, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, iwk(listev) )
    end if

    if ( nev == 0 ) then
      psi(i) = 1.0D+00 / wsq
      wk(i) = 0.0D+00
      mdfint = psi(i) * area(i)
    else

      nvrt = xivrt(i+1) - xivrt(i)
      k = xivrt(i)

      sumx = 0.0D+00
      sumy = 0.0D+00
      do j = 0, nvrt-1
        l = ivrt(k)
        wk(xc+j) = vcl(1,l)
        wk(yc+j) = vcl(2,l)
        sumx = sumx + wk(xc+j)
        sumy = sumy + wk(yc+j)
        k = k + 1
      end do

      ctrx = sumx / real ( nvrt, kind = 8 )
      ctry = sumy / real ( nvrt, kind = 8 )

      do j = 0, nvrt-1
        wk(xc+j) = wk(xc+j) - ctrx
        wk(yc+j) = wk(yc+j) - ctry
      end do

      wk(xc+nvrt) = wk(xc)
      wk(yc+nvrt) = wk(yc)

      call intpg ( nvrt, wk(xc), wk(yc), ctrx, ctry, area(i), hflag, umdf, &
        wsq, nev, ifv, iwk(listev), ivrt, edgval, vrtval, vcl, mdfint, &
        psi(i), wk(i), wk(mdftr) )

    end if

    arearg = arearg + area(i)
    intreg = intreg + mdfint

  end do
!
!  If HFLAG, compute mean mdf values from KAPPA, etc. Scale PSI(I)'s
!  so that integral in region is 1. Determine which polygons need to
!  be further subdivided (indicated by negative PSI(I) value).
!
  if ( hflag ) then
    c1 = ( 1.0D+00 - kappa ) / intreg
    c2 = kappa / arearg
  else
    c1 = 1.0D+00 / intreg
    c2 = 0.0D+00
  end if

  do i = 1, npolg

    psi(i) = psi(i) * c1 + c2

    if ( psi(i) * dmin < c1 * wk(i) ) then
      if ( nmin < ntrid * psi(i) * area(i) ) then
        psi(i) = -psi(i)
      end if
    end if

  end do
!
!  Further subdivide polygons for which DMIN < STDV/MEAN and
!  NMIN < (estimated number of triangles).
!
  angsp2 = 2.0D+00 * angspc
  pi2 = 2.0D+00 * pi
  inc = int ( pi2 / angspc )
  nev = 0
  np = npolg
  xc = 1

  do i = 1, np

    if ( psi(i) < 0.0D+00 ) then

      if ( hflag ) then
        wsq = widsq(i)
        call prmdf2 ( i, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, &
          iwk(listev) )

      end if

      l = npolg + 1
      k = i

60    continue

      if ( npolg < k ) then
        go to 130
      end if

70       continue

         if ( 0.0D+00 <= psi(k) ) then
           go to 120
         end if

          nvrt = 0
          sumx = 0.0D+00
          sumy = 0.0D+00
          j = hvl(k)

          do

            nvrt = nvrt + 1
            m = pvl(loc,j)
            sumx = sumx + vcl(1,m)
            sumy = sumy + vcl(2,m)
            j = pvl(succ,j)

            if ( j == hvl(k) ) then
              exit
            end if

          end do

          ctrx = sumx / real ( nvrt, kind = 8 )
          ctry = sumy / real ( nvrt, kind = 8 )
          maxn = nvrt + inc

          if ( maxiw < nev + maxn + 1 ) then
            ierror = 6
            return
          else if ( maxwk < 3*maxn + 2 ) then
            ierror = 7
            return
          end if

          yc = xc + maxn + 1
          mdftr = yc + maxn + 1
          indpvl = listev + nev
          nvrt = 0
          m = pvl(loc,j)
          x1 = vcl(1,m) - ctrx
          y1 = vcl(2,m) - ctry
          wk(xc) = x1
          wk(yc) = y1
          theta1 = atan2(y1,x1)
          p = j
          iwk(indpvl) = j

90        continue

             j = pvl(succ,j)
             m = pvl(loc,j)
             x2 = vcl(1,m) - ctrx
             y2 = vcl(2,m) - ctry

             theta2 = atan2 ( y2, x2 )
             if ( theta2 < theta1 ) then
               theta2 = theta2 + pi2
             end if

             delta = theta2 - theta1

             if ( angsp2 <= delta ) then

               m = int ( delta / angspc )
               delta = delta / real ( m, kind = 8 )
               dx = x2 - x1
               dy = y2 - y1
               numer = x1 * dy - y1 * dx
               alpha = theta1

               do ii = 1, m-1
                 alpha = alpha + delta
                 cosalp = cos(alpha)
                 sinalp = sin(alpha)
                 r = numer / ( dy * cosalp - dx * sinalp )
                 nvrt = nvrt + 1
                 wk(xc+nvrt) = r * cosalp
                 wk(yc+nvrt) = r * sinalp
                 iwk(indpvl+nvrt) = -p
               end do

             end if

             nvrt = nvrt + 1
             wk(xc+nvrt) = x2
             wk(yc+nvrt) = y2
             x1 = x2
             y1 = y2
             theta1 = theta2
             p = j
             iwk(indpvl+nvrt) = j

          if ( j /= hvl(k) ) then
            go to 90
          end if


            call intpg ( nvrt, wk(xc), wk(yc), ctrx, ctry, area(k), hflag, &
              umdf, wsq, nev, ifv, iwk(listev), ivrt, edgval, vrtval, &
              vcl, mdfint, mean, stdv, wk(mdftr) )

            psi(k) = mean * c1 + c2

            if ( psi(k) * dmin < c1 * stdv ) then

               if ( nmin < ntrid * psi(k) * area(k) ) then

              call sepmdf ( angtol, nvrt, wk(xc), wk(yc), area(k), &
                          mean, wk(mdftr), iwk(indpvl), iang, i1, i2 )

              if ( i1 < 0 ) then

                if ( maxwk < yc + 3*nvrt ) then
                  ierror = 7
                  return
                end if

                call sepshp ( angtol, nvrt, wk(xc), wk(yc), &
                  iwk(indpvl), iang, i1, i2, wk(yc+nvrt+1), ierror )

                if ( ierror /= 0 ) then
                  return
                end if

              end if

              if ( i1 < 0 ) then
                ierror = 222
                return
              end if

              v = iwk(indpvl+i1)

              if ( v < 0 ) then
                 call insvr2 ( wk(xc+i1)+ctrx, wk(yc+i1)+ctry, -v, &
                   nvc, nvert, maxvc, maxpv, vcl, pvl, iang, v, ierror )
                 if ( ierror /= 0 ) then
                   return
                 end if
              end if

              w = iwk(indpvl+i2)

              if ( w < 0 ) then
                 call insvr2 ( wk(xc+i2)+ctrx, wk(yc+i2)+ctry, -w, &
                   nvc, nvert, maxvc, maxpv, vcl, pvl, iang, w, ierror )
                 if ( ierror /= 0 ) then
                   return
                 end if
              end if

          call insed2 ( v, w, npolg, nvert, maxhv, maxpv, vcl, &
            regnum, hvl, pvl, iang, ierror )

          if ( ierror /= 0 ) then
            return
          end if

          nvrt = 0
          j = hvl(k)

          do

            m = pvl(loc,j)
            wk(xc+nvrt) = vcl(1,m)
            wk(yc+nvrt) = vcl(2,m)
            nvrt = nvrt + 1
            j = pvl(succ,j)

            if ( j == hvl(k) ) then
              exit
            end if

          end do

          nwarea = areapg ( nvrt, wk(xc), wk(yc) ) * 0.5D+00
          area(npolg) = area(k) - nwarea
          area(k) = nwarea
          psi(k) = -psi(k)
          psi(npolg) = psi(k)

        end if

      end if

      go to 70

120   continue

      if ( k == i ) then
        k = l
      else
        k = k + 1
      end if

      go to 60

130   continue

    end if

  end do

  return
end
function minang ( xr, yr, xs, ys, ind, alpha, theta, vcl, pvl, iang )

!*****************************************************************************80
!
!! MINANG determines the minimum of the boundary angles for a separator.
!
!  Purpose: 
!
!    Determine the minimum of the 4 angles at the boundary
!    resulting from using edge joining vertices (XR,YR) and
!    (XS,YS) as a separator.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XR, YR, the coordinates of the reflex vertex.
!
!    Input, real ( kind = 8 ) XS, YS, the coordinates of other endpoint of 
!    possible separator.
!
!    Input, integer ( kind = 4 ) IND, if positive then (XS,YS) has index IND in PVL; else
!    (XS,YS) is on edge joining vertices with indices -IND
!    and SUCC(-IND) in PVL.
!
!    Input, real ( kind = 8 ) ALPHA, the polar angle of (XS,YS) with respect
!    to (XR,YR).
!
!    Input, real ( kind = 8 ) THETA, the interior angle at reflex vertex.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), real ( kind = 8 ) IANG(1:*), the polygon 
!    vertex list, interior angles.
!
!    Output, real ( kind = 8 ) MINANG, the minimum of the 4 angles in radians.
!
  real ( kind = 8 ) alpha
  real ( kind = 8 ) ang
  real ( kind = 8 ) angle
  real ( kind = 8 ) beta1
  real ( kind = 8 ) iang(*)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  real ( kind = 8 ) minang
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) pvl(4,*)
  integer ( kind = 4 ), parameter :: succ = 3
  real ( kind = 8 ) theta
  real ( kind = 8 ) vcl(2,*)
  real ( kind = 8 ) xr
  real ( kind = 8 ) xs
  real ( kind = 8 ) yr
  real ( kind = 8 ) ys

  if ( 0 < ind ) then
    j = pvl(succ,ind)
    ang = iang(ind)
  else
    j = pvl(succ,-ind)
    ang = pi
  end if

  l = pvl(loc,j)
  beta1 = angle ( xr, yr, xs, ys, vcl(1,l), vcl(2,l) )

  minang = min ( alpha, theta - alpha, ang - beta1, beta1 )

  return
end
subroutine mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )

!*****************************************************************************80
!
!! MMASEP chooses the best of four separators by the max-min angle criterion.
!
!  Purpose: 
!
!    Find best of four possible separators according to
!    max-min angle criterion.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter (in 
!    radians) for accepting separator.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the coordinates of 
!    polygon vertices in counter clockwise order where NVRT is number of 
!    vertices; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, integer ( kind = 4 ) INDPVL(0:NVRT), indices in PVL of vertices; INDPVL(I) = -K
!    if (XC(I),YC(I)) is extra vertex inserted on edge from K to PVL(SUCC,K).
!
!    Input, real ( kind = 8 ) IANG(1:*), the interior angle array.
!
!    Input, integer ( kind = 4 ) V(1:2), W(1:2), indices in XC, YC in range 0 to NVRT-1; 
!    four possible separators are V(I),W(J), I,J = 1,2.
!
!    Output, integer ( kind = 4 ) I1, I2, indices in range 0 to NVRT-1 of best separator
!    according to max-min angle criterion; I1 = -1
!    if no satisfactory separator is found.
!
  real ( kind = 8 ) alpha
  real ( kind = 8 ) angle
  real ( kind = 8 ) angmax
  real ( kind = 8 ) angmin
  real ( kind = 8 ) angtol
  real ( kind = 8 ) beta
  real ( kind = 8 ) delta
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real ( kind = 8 ) iang(*)
  integer ( kind = 4 ) indpvl(0:*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) tol
  integer ( kind = 4 ) v(2)
  integer ( kind = 4 ) w(2)
  real ( kind = 8 ) xc(0:*)
  real ( kind = 8 ) yc(0:*)

  tol = 100.0D+00 * epsilon ( tol )
  angmax = 0.0D+00

  do i = 1, 2

    l = v(i)
    k = indpvl(l)

    if ( 0 < k ) then
      alpha = iang(k)
    else
      alpha = pi
    end if

    do j = 1, 2

      m = w(j)

      if ( l /= m ) then

        k = indpvl(m)

        if ( 0 < k ) then
          beta = iang(k)
        else
          beta = pi
        end if

        gamma = angle ( xc(m), yc(m), xc(l), yc(l), xc(l+1), yc(l+1) )

        delta = angle ( xc(l), yc(l), xc(m), yc(m), xc(m+1), yc(m+1) )

        angmin = min ( gamma, alpha-gamma, delta, beta-delta )

        if ( angmax < angmin ) then
          angmax = angmin
          i1 = l
          i2 = m
        end if

      end if

    end do

  end do

  if ( angmax < angtol ) then
    i1 = -1
  end if

  return
end
subroutine mtredg ( utype, i1, i2, i3, ibndry, nt, til, tedg )

!*****************************************************************************80
!
!! MTREDG sets fields for a triangle as needed by routine TMERGE.
!
!
!  Purpose: 
!
!    Set fields for triangle as needed by routine TMERGE.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical UTYPE, is .TRUE. iff triangle contains two 'U' vertices.
!
!    Input, integer ( kind = 4 ) I1, I2, I3, the indices of 3 triangle vertices in VCL; 
!    the first two indices also belong to the next merge edge.
!
!    Input, integer ( kind = 4 ) IBNDRY, the index of boundary edge for TEDG.
!
!    Input/output, integer ( kind = 4 ) NT, the number of entries in TIL, TEDG so far.
!
!    Input/output, integer ( kind = 4 ) TIL(1:NT), the triangle incidence list.
!
!    Input/output, TEDG(1:NT), the triangle edge indices; see routine TMERGE.
!
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ibndry
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) tedg(3,*)
  integer ( kind = 4 ) til(3,*)
  logical utype

  nt = nt + 1
  til(1,nt) = i1
  til(2,nt) = i2
  til(3,nt) = i3
  tedg(1,nt) = nt

  if ( utype ) then
    tedg(2,nt) = nt - 1
    tedg(3,nt) = ibndry
  else
    tedg(2,nt) = ibndry
    tedg(3,nt) = nt - 1
  end if

  return
end
function prime ( k )

!*****************************************************************************80
!
!! PRIME returns a prime greater than a given integer K.
!
!  Purpose: 
!
!    Return a prime greater than or equal to K (if possible) from internal 
!    array of primes.  More primes can be added if desired.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, a positive integer.
!
!    Output, integer ( kind = 4 ) PRIME, the smallest prime greater than or equal to
!    K from internal array (or largest in array).
!
  integer ( kind = 4 ), parameter :: nprime = 150

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) prime
  integer ( kind = 4 ), save, dimension ( nprime ) :: primes = (/ &
    17,31,47,61,79,97,113,127,149,163,179,193,211,227,241, &
    257,271,293,307,331,353,379,401,431,457,479,503,541,563,587, &
    613,641,673,701,727,751,773,797,821,853,877,907,929,953,977, &
    1009,1049,1087,1123,1163,1201,1237,1277,1319,1361,1399,1433, &
    1471,1511,1543,1579,1613,1657,1699,1741,1783,1831,1873,1931, &
    1973,2017,2069,2129,2203,2267,2333,2389,2441,2503,2557,2609, &
    2663,2719,2789,2851,2917,2999,3061,3137,3209,3299,3371,3449, &
    3527,3613,3697,3779,3863,3947,4049,4211,4421,4621,4813,5011, &
    5227,5413,5623,5813,6011,6211,6421,6619,6823,7013,7211,7411, &
    7621,7817,8011,8219,8419,8623,8819,9011,9221,9413,9613,9811, &
    10037,10211,10427,10613,10831,11027,11213,11411,11617,11813, &
    12011,12211,12413,12611,12821,13033,13217,13411,13613,13829, &
    14011 /)
  integer ( kind = 4 ) u

  if ( k <= primes(1) ) then
    prime = primes(1)
    return
  else if ( primes(nprime) <= k ) then
    prime = primes(nprime)
    return
  end if
!
!  Use binary search to find prime greater than or equal to K.
!
  l = 1
  u = nprime

  do

    m = ( l + u ) / 2

    if ( k < primes(m) ) then
      u = m - 1
    else if ( primes(m) < k ) then
      l = m + 1
    else
      prime = primes(m)
      return
    end if

    if ( u < l ) then
      exit
    end if

  end do

  prime = primes(u+1)

  return
end
subroutine prmdf2 ( ipoly, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, listev )

!*****************************************************************************80
!
!! PRMDF2 preprocesses a mesh distribution function evaluation.
!
!  Purpose: 
!
!    Preprocessing step for evaluating mesh distribution
!    function in polygon IPOLY.  The edges and vertices for
!    which distances must be computed are determined.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IPOLY, the index of the polygon.
!
!    Input, real ( kind = 8 ) WSQ, the square of the width of polygon IPOLY.
!
!    Input, integer ( kind = 4 ) IVRT(1:*), the indices of polygon vertices in 
!    VCL, ordered by polygon.
!
!    Input, integer ( kind = 4 ) XIVRT(1:*), pointers to first vertex of each polygon 
!    in IVRT; vertices of polygon IPOLY are IVRT(I) for I from
!    XIVRT(IPOLY) to XIVRT(IPOLY+1)-1.
!
!    Input, real ( kind = 8 ) EDGVAL(1:*), a value associated with each edge 
!    of the decomposition.
!
!    Input, real ( kind = 8 ) VRTVAL(1:*), a value associated with each vertex 
!    of the decomposition.
!
!    Output, integer ( kind = 4 ) NEV, the number of edges and vertices for which distances
!    must be evaluated.
!
!    Output, integer ( kind = 4 ) IFV, the index of first vertex XIVRT(IPOLY) if LISTEV(NEV)
!    = XIVRT(IPOLY+1) - 1; 0 otherwise.
!
!    Output, integer ( kind = 4 ) LISTEV(1:*), an array of length 
!    <= [XIVRT(IPOLY+1)-XIVRT(IPOLY)]*2, containing indices of edges and 
!    vertices mentioned above; indices of vertices are negated.
!
  real ( kind = 8 ) edgval(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifv
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) listev(*)
  integer ( kind = 4 ) nev
  real ( kind = 8 ) vrtval(*)
  real ( kind = 8 ) wsq
  integer ( kind = 4 ) xivrt(*)

  ifv = 0
  nev = 0
  im1 = xivrt(ipoly+1) - 1
  l = im1

  do i = xivrt(ipoly), l

    j = ivrt(i)

    if ( vrtval(j) < min ( edgval(i), edgval(im1) ) ) then
      nev = nev + 1
      listev(nev) = -j
    end if

    if ( edgval(i) < wsq ) then
      nev = nev + 1
      listev(nev) = i
    end if

    im1 = i

  end do

  if ( 0 < nev ) then
    if ( listev(nev) == l ) then
      ifv = xivrt(ipoly)
    end if
  end if

  return
end
subroutine ptpolg ( dim, ldv, nv, inc, pgind, vcl, pt, nrml, dtol, inout )

!*****************************************************************************80
!
!! PTPOLG determines if a point is in, on or outside a polygon.
!
!  Purpose: 
!
!    Determine whether a point lies inside, outside, or on
!    boundary of a planar polygon in 2 or 3 dimensional space.
!    It is assumed that point lies in plane of polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM, the dimension of the polygon (2 or 3).
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL array in calling routine.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices in polygon.
!
!    Input, integer ( kind = 4 ) INC, the increment for PGIND indicating indices of polygon.
!
!    Input, integer ( kind = 4 ) PGIND(0:NV*INC), indices in VCL of polygon vertices are in
!    PGIND(0), PGIND(INC), ..., PGIND(NV*INC) with first and
!    last vertices identical.
!
!    Input, real ( kind = 8 ) VCL(1:DIM,1:*), the vertex coordinate list.
!
!    Input, real ( kind = 8 ) PT(1:DIM), the point for which in/out test is
!    applied.
!
!    Input, real ( kind = 8 ) NRML(1:3), the unit normal vector of plane 
!    containing polygon, with vertices oriented counter clockwise with
!    respect to the normal (used iff DIM = 3);
!    The normal is assumed to be (0,0,1) if DIM = 2.
!
!    Input, real ( kind = 8 ) DTOL, an absolute tolerance to determine 
!    whether a point is on a line or plane.
!
!    Output, integer ( kind = 4 ) INOUT, point PT is:
!    +1, inside the polygon, 
!     0, on boundary of polygon, 
!    -1, outside polygon;
!    -2 if error in input parameters.
!
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) ldv

  real ( kind = 8 ) cp(3)
  real ( kind = 8 ) de(3)
  real ( kind = 8 ) dir(3)
  real ( kind = 8 ) dist
  real ( kind = 8 ) dotp
  real ( kind = 8 ) dtol
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) inout
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  real ( kind = 8 ) len1
  real ( kind = 8 ) len2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) nr(4)
  real ( kind = 8 ) nrml(3)
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) pgind(0:*)
  real ( kind = 8 ) pt(dim)
  real ( kind = 8 ) rhs(3)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) sa
  integer ( kind = 4 ) sb
  real ( kind = 8 ) t
  real ( kind = 8 ) ta
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(ldv,*)

  tol = 100.0D+00 * epsilon ( tol )

  inout = -2

  if ( dim < 2 .or. 3 < dim ) then
    return
  end if
!
!  Direction of ray is from PT through midpoint of first edge
!  such that PT is not collinear with edge. NR is normal of plane
!  containing ray, which is also orthogonal to NRML.
!
  i = 0
  lb = pgind(0)

10 continue

    i = i + 1

    if ( nv <= i ) then
      return
    end if

    la = lb
    lb = pgind(i*inc)

    do j = 1, dim
      de(j) = vcl(j,lb) - vcl(j,la)
      dir(j) = pt(j) - vcl(j,la)
    end do

    if ( dim == 2 ) then
      len1 = de(1)**2 + de(2)**2
      len2 = dir(1)**2 + dir(2)**2
    else
      len1 = de(1)**2 + de(2)**2 + de(3)**2
      len2 = dir(1)**2 + dir(2)**2 + dir(3)**2
    end if

    if ( len1 == 0.0D+00 ) then
      go to 10
    else if ( len2 == 0.0D+00 ) then
      inout = 0
      return
    end if

    if ( dim == 2 ) then
      dotp = abs ( de(1) * dir(1) + de(2) * dir(2)) / sqrt(len1*len2)
    else if ( dim == 3 ) then
      dotp = abs ( de(1) * dir(1) + de(2) * dir(2) + de(3) * dir(3) ) &
        / sqrt(len1*len2)
    end if

    if ( 1.0D+00 - tol <= dotp ) then
      go to 10
    end if

    if ( dim == 2 ) then
      dir(1) = 0.5D+00 * ( vcl(1,la) + vcl(1,lb) ) - pt(1)
      dir(2) = 0.5D+00 * ( vcl(2,la) + vcl(2,lb) ) - pt(2)
      dist = sqrt ( dir(1)**2 + dir(2)**2 )
      dir(1) = dir(1) / dist
      dir(2) = dir(2) / dist
      dir(3) = 0.0D+00
      nr(1) = -dir(2)
      nr(2) = dir(1)
      nr(3) = 0.0D+00
      nr(4) = nr(1) * pt(1) + nr(2) * pt(2)
      dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) - nr(4)
    else if ( dim == 3 ) then
      dir(1) = 0.5D+00 * ( vcl(1,la) + vcl(1,lb) ) - pt(1)
      dir(2) = 0.5D+00 * ( vcl(2,la) + vcl(2,lb) ) - pt(2)
      dir(3) = 0.5D+00 * ( vcl(3,la) + vcl(3,lb) ) - pt(3)
      dist = sqrt ( dir(1)**2 + dir(2)**2 + dir(3)**2 )

      dir(1) = dir(1) / dist
      dir(2) = dir(2) / dist
      dir(3) = dir(3) / dist

      nr(1) = nrml(2)*dir(3) - nrml(3)*dir(2)
      nr(2) = nrml(3)*dir(1) - nrml(1)*dir(3)
      nr(3) = nrml(1)*dir(2) - nrml(2)*dir(1)
      nr(4) = nr(1)*pt(1) + nr(2)*pt(2) + nr(3)*pt(3)
      dist = nr(1)*vcl(1,lb)+nr(2)*vcl(2,lb)+nr(3)*vcl(3,lb) - nr(4)
    end if

    if ( 0.0D+00 < dist ) then
      sb = 1
    else
      sb = -1
    end if

    m = 1
    if ( abs ( dir(1) ) < abs ( dir(2) ) ) then
      m = 2
    end if

    if ( abs ( dir(m) ) < abs ( dir(3) ) ) then
      m = 3
    end if

    k = 1
!
!  For remaining edges of polygon, check whether ray intersects edge.
!  Vertices or edges lying on ray are handled by looking at preceding
!  and succeeding edges not lying on ray.
!
  n = i
  i = i + 1

30 continue

   la = lb
   lb = pgind(i*inc)
   sa = sb

   if ( dim == 2 ) then
     dist = nr(1) * vcl(1,lb) + nr(2)*vcl(2,lb) - nr(4)
   else if ( dim == 3 ) then
     dist = nr(1) * vcl(1,lb) + nr(2)*vcl(2,lb) + nr(3)*vcl(3,lb)- nr(4)
   end if

   if ( abs ( dist) <= dtol ) then
     sb = 0
   else if ( 0.0D+00 < dist ) then
     sb = 1
   else
     sb = -1
   end if

   s = sa * sb

   if ( s < 0 ) then

      if ( dim == 2 ) then
         de(1) = vcl(1,la) - vcl(1,lb)
         de(2) = vcl(2,la) - vcl(2,lb)
         rhs(1) = vcl(1,la) - pt(1)
         rhs(2) = vcl(2,la) - pt(2)
         t = ( rhs(1) * de(2) - rhs(2) * de(1) ) &
           / ( dir(1) * de(2) - dir(2) * de(1) )
      else if ( dim == 3 ) then
         de(1) = vcl(1,la) - vcl(1,lb)
         de(2) = vcl(2,la) - vcl(2,lb)
         de(3) = vcl(3,la) - vcl(3,lb)
         rhs(1) = vcl(1,la) - pt(1)
         rhs(2) = vcl(2,la) - pt(2)
         rhs(3) = vcl(3,la) - pt(3)
         cp(1) = dir(2) * de(3) - dir(3) * de(2)
         cp(2) = dir(3) * de(1) - dir(1) * de(3)
         cp(3) = dir(1) * de(2) - dir(2) * de(1)

         l = 1
         if ( abs ( cp(1) ) < abs ( cp(2) ) ) then
           l = 2
         end if
         if ( abs ( cp(l) ) < abs ( cp(3) ) ) then
           l = 3
         end if

         if ( l == 1 ) then
           t = ( rhs(2) * de(3) - rhs(3) * de(2) ) / cp(1)
         else if ( l == 2 ) then
           t = ( rhs(3) * de(1) - rhs(1) * de(3) ) / cp(2)
         else
           t = ( rhs(1) * de(2) - rhs(2) * de(1) ) / cp(3)
         end if

      end if

      if ( dtol < t ) then
         k = k + 1
      else if ( -dtol <= t ) then
         inout = 0
         return
      end if

   else if ( s == 0 ) then

      l = lb
40    continue
      i = i + 1

      if ( nv < i ) then
        i = 1
      end if

      if ( i == n) then
        return
      end if

      la = lb
      lb = pgind(i*inc)

    if ( dim == 2 ) then
      dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) - nr(4)
    else
      dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) + nr(3) * vcl(3,lb) - nr(4)
    end if

    if ( abs ( dist) <= dtol ) then
      go to 40
    else if ( 0.0D+00 < dist ) then
      sb = 1
    else
      sb = -1
    end if

    t = ( vcl(m,l) - pt(m) ) / dir(m)

    if ( abs ( t) <= dtol ) then
      inout = 0
      return
    end if

    if ( la /= l ) then
      ta = ( vcl(m,la) - pt(m) ) / dir(m)
      if ( abs ( ta ) <= dtol .or. t * ta < 0.0D+00 ) then
        inout = 0
        return
      end if
    end if

    if ( sa * sb < 0 .and. 0.0D+00 < t ) then
      k = k + 1
    end if

  end if

  i = i + 1

  if ( nv < i ) then
    i = 1
  end if

  if ( i /= n ) then
    go to 30
  end if
!
!  Point lies inside polygon if number of intersections K is odd.
!
  if ( mod ( k, 2 ) == 1 ) then
    inout = 1
  else
    inout = -1
  end if

  return
end
subroutine r8mat_print ( lda, m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Modified:
!
!    08 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) m
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 5
    jhi = min ( jlo + 4, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,5(i7,7x))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,5g14.6)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
function radians_to_degrees ( angle )

!*****************************************************************************80
!
!! RADIANS_TO_DEGREES converts an angle from radians to degrees.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, an angle in radians.
!
!    Output, real ( kind = 8 ) RADIANS_TO_DEGREES, the equivalent angle
!    in degrees.
!
  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radians_to_degrees

  radians_to_degrees = ( angle / pi ) * 180.0D+00

  return
end
subroutine randpt ( k, n, seed, axis, nptav, scale, trans, lda, a )

!*****************************************************************************80
!
!! RANDPT generates N random K-dimensional points from the uniform distribution.
!
!  Purpose: 
!
!    Generate N random K-dimensional points from the uniform distribution.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) N, the number of random points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for pseudo random number generator.
!
!    Input, integer ( kind = 4 ) AXIS, integer NPTAV; if AXIS < 1 or K < AXIS, then uniform 
!    random points are generated; if 1 <= AXIS <= K, then an average of NPTAV
!    uniform random points are generated with the same AXIS
!    coordinate on about N/NPTAV random parallel hyperplanes.
!
!    Input, real ( kind = 8 ) SCALE(K), TRANS(K), the scale and 
!    translation factors for coordinates 1 to K; the I-th coordinate of 
!    random point is R*SCALE(I) + TRANS(I) where 0 < R < 1.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling
!    routine; should be at least K.
!
!    Output, real ( kind = 8 ) A(LDA,N), an array of N uniform random 
!    K-dimensional points.
!
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) axis
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nptav
  real ( kind = 8 ) r
  real ( kind = 8 ) scale(k)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) trans(k)
  real ( kind = 8 ) urand

  if ( axis < 1 .or. k < axis ) then

    do j = 1, n
      do i = 1, k
        a(i,j) = trans(i) + scale(i) * urand ( seed )
      end do
    end do

  else

    m = int ( urand ( seed ) * 2.0D+00 * nptav + 0.5D+00 )
    r = urand ( seed ) * scale(axis) + trans(axis)

    do j = 1, n

      do i = 1, k
        if ( i == axis ) then
          a(i,j) = r
        else
          a(i,j) = urand ( seed ) * scale(i) + trans(i)
        end if
      end do

      m = m - 1

      if ( m <= 0 ) then
        m = int ( urand ( seed ) * 2.0D+00 * nptav + 0.5D+00 )
        r = urand ( seed ) * scale(axis) + trans(axis)
      end if

    end do

  end if

  return
end
subroutine resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw,  &
  maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

!*****************************************************************************80
!
!! RESVRT resolves a reflex vertex of a simply connected polygon.
!
!  Purpose: 
!
!    Resolve a reflex vertex of a simply connected polygon with
!    one or two separators. 
!
!    The reflex vertex must be a 'simple' vertex of the polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) VR, the index in PVL of reflex vertex.
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter used in
!    controlling the vertices to be considered as an endpoint of a separator.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter used in 
!    accepting separator(s).
!
!    Input/output, integer ( kind = 4 ) NVC, the number of positions used 
!    in VCL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should 
!    be about 3 times number of vertices in polygon.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should 
!    be about 5 times number of vertices in polygon.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT),
!    the polygon vertex list and interior angles.
!
!    Output, integer ( kind = 4 ) W1, the index in PVL of vertex which is the endpoint 
!    of separator in inner cone or right cone with respect to the reflex vertex.
!
!    Output, integer ( kind = 4 ) W2, is 0 if there is only one separator; else index 
!    in PVL of vertex which is endpoint of second separator in left cone.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 5, 6, 7, 206, 207, 208, 209, 210, or 212
!
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real ( kind = 8 ) angsep
  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ivis
  integer ( kind = 4 ) ivor
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) maxn
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) nvor
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) nvsvrt
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) theta
  integer ( kind = 4 ) v
  real ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vr
  integer ( kind = 4 ) w1
  integer ( kind = 4 ) w2
  real ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) wkang
  integer ( kind = 4 ) xc
  real ( kind = 8 ) xr
  integer ( kind = 4 ) xvor
  integer ( kind = 4 ) yc
  real ( kind = 8 ) yr
  integer ( kind = 4 ) yvor
!
!  Determine number of vertices in polygon containing reflex vertex.
!
  ierror = 0
  nvrt = 0
  v = vr

  do

    v = pvl(succ,v)

    if ( v == vr ) then
      exit
    end if

    nvrt = nvrt + 1

  end do

  maxn = nvrt + int ( iang(vr) / angspc )
  l = pvl(loc,vr)
  xr = vcl(1,l)
  yr = vcl(2,l)
!
!  Set up work arrays for routine VISPOL, and determine whether there
!  is enough workspace. XC, YC are d.p. arrays of length NVRT in WK,
!  used for the coordinates of the polygon containing the reflex
!  vertex. MAXN positions are reserved for XC, YC since this is the
!  maximum space required by routine VISVRT. IVIS is an integer array
!  of length MAXN in IWK. IVRT is an integer array of length NVRT in
!  IWK used temporarily for storing indices of vertices in PVL.
!
  if ( maxiw < maxn + nvrt ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESVRT - Fatal error!'
    write ( *, '(a)' ) '  MAXIW < MAXN + NVRT.'
    ierror = 6
    return
  end if

  if ( maxwk < maxn + maxn ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESVRT - Fatal error!'
    write ( *, '(a)' ) '  MAXWK < MAXN + MAXN.'
    ierror = 7
    return
  end if

  ivis = 1
  ivrt = ivis + maxn
  xc = 1
  yc = xc + maxn
  v = pvl(succ,vr)

  do i = 0, nvrt-1
    l = pvl(loc,v)
    wk(xc+i) = vcl(1,l)
    wk(yc+i) = vcl(2,l)
    iwk(ivrt+i) = v
    v = pvl(succ,v)
  end do

  call vispol ( xr, yr, nvrt-1, wk(xc), wk(yc), nvis, iwk(ivis), ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  XC, YC now contain visibility polygon coordinates. Update MAXN
!  and set up d.p. array THETA of length MAXN in WK for routine
!  VISVRT. Elements of IVIS are changed to indices of PVL after call.
!
  maxn = maxn - nvrt + nvis + 1
  theta = yc + maxn

  if ( maxwk < theta + maxn - 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESVRT - Fatal error!'
    write ( *, '(a)' ) '  MAXWK < THETA + MAXN - 1.'
    ierror = 7
    return
  end if

  call visvrt ( angspc, xr, yr, nvis, wk(xc), wk(yc), iwk(ivis), maxn-1, &
    nvsvrt, wk(theta) )

  wk(theta+nvsvrt) = iang(vr)

  do i = ivis, ivis+nvsvrt
    l = iwk(i)
    if ( 0 <= l ) then
      iwk(i) = iwk(ivrt+l)
    else
      iwk(i) = -iwk(ivrt-l-1)
    end if
  end do
!
!  XC, YC now contain coordinates of visible vertices to be considered
!  as an endpoint of a separator. Set up work arrays for routine
!  VORNBR. Integer array IVOR and d.p. arrays XVOR, YVOR, each of
!  length NVSVRT+1, are added at the end of IWK and WK arrays.
!
  ivor = ivis + nvsvrt + 1
  xvor = theta + nvsvrt + 1
  yvor = xvor + nvsvrt + 1

  if ( maxiw < ivor + nvsvrt ) then
    ierror = 6
    return
  end if

  if ( maxwk < yvor + nvsvrt ) then
    ierror = 7
    return
  end if

  call vornbr ( xr, yr, nvsvrt, wk(xc), wk(yc), nvor, iwk(ivor), wk(xvor), &
    wk(yvor), ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Set up the array WKANG of length NVOR+1 <= NVSVRT+1 in WK for
!  routine FNDSEP. Only Voronoi neighbors are considered as an
!  endpoint of a separator in the first call to FNDSEP. If the
!  minimum angle created at the boundary by the separator(s) is too
!  small, then a second call is made to FNDSEP in which all visible
!  vertices are considered as an endpoint of a separator.
!
  wkang = xvor
  if ( iwk(ivor+nvor) == nvsvrt ) then
    nvor = nvor - 1
  end if

  if ( iwk(ivor) == 0 ) then
    ivor = ivor + 1
    nvor = nvor - 1
  end if

  call fndsep ( angtol+angtol, xr, yr, nvsvrt, wk(xc), wk(yc), iwk(ivis), &
    wk(theta), nvor, iwk(ivor), vcl, pvl, iang, angsep, i1, i2, wk(wkang) )

  if ( angsep < angtol ) then

    ivrt = ivis + nvsvrt + 1

    do i = 1, nvsvrt-1
      iwk(ivrt+i-1) = i
    end do

    call fndsep ( angtol+angtol, xr, yr, nvsvrt, wk(xc), wk(yc), iwk(ivis), &
      wk(theta), nvsvrt-2, iwk(ivrt), vcl, pvl, iang, angsep, i1, i2, &
      wk(wkang) )

  end if
!
!  Insert endpoint(s) of separator(s) in vertex coordinate list and
!  polygon vertex list data structures, if they are not yet there.
!
  if ( i2 == -1 ) then
    w2 = 0
  else if ( iwk(ivis+i2) < 0 ) then
    call insvr2 ( wk(xc+i2), wk(yc+i2), -iwk(ivis+i2), nvc, nvert, maxvc, &
      maxpv, vcl, pvl, iang, w2, ierror )
    if ( ierror /= 0 ) then
      return
    end if
  else
    w2 = iwk(ivis+i2)
  end if

  if ( iwk(ivis+i1) < 0 ) then
    call insvr2 ( wk(xc+i1), wk(yc+i1), -iwk(ivis+i1), nvc, nvert, maxvc, &
      maxpv, vcl, pvl, iang, w1, ierror )
    if ( ierror /= 0 ) then
      return
    end if
  else
    w1 = iwk(ivis+i1)
  end if

  return
end
subroutine rotiar ( n, arr, shift )

!*****************************************************************************80
!
!! ROTIAR rotates the elements of an integer array.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of array.
!
!    Input/output, integer ( kind = 4 ) ARR(0:N-1), the array to be shifted.
!
!    Input, integer ( kind = 4 ) SHIFT, the amount of (left) shift or rotation; 
!    ARR(SHIFT) on input becomes ARR(0) on output.
!
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) arr(0:n-1)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  integer ( kind = 4 ) sh
  integer ( kind = 4 ) shift
  integer ( kind = 4 ) t

  sh = mod ( shift, n )

  if ( sh < 0 ) then
    sh = sh + n
  end if

  if ( sh == 0 ) then
    return
  end if

  a = n
  b = sh

  do
    r = mod ( a, b )
    a = b
    b = r
    if ( r <= 0 ) then
      exit
    end if
  end do

  m = n / a - 1

  do i = 0, a-1

    t = arr(i)
    k = i

    do j = 1, m
      l = k + sh
      if ( n <= l ) then
        l = l - n
      end if
      arr(k) = arr(l)
      k = l
    end do

    arr(k) = t

  end do

  return
end
subroutine rotipg ( xeye, yeye, nvrt, xc, yc, ierror )

!*****************************************************************************80
!
!! ROTIPG rotates the vertex indices of a simple polygon.
!
!  Purpose: 
!
!    Rotate the indices of the vertices of a simple polygon
!    and possibly insert one vertex so that (XC(0),YC(0)) is the
!    point on the horizontal line through (XEYE,YEYE) and on the
!    boundary of the polygon which is closest to and to the right
!    of (XEYE,YEYE). (XEYE,YEYE) is an eyepoint in the interior or
!    blocked exterior of the polygon. In the former (latter) case,
!    the vertices must be in counter clockwise (CW) order.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XEYE, YEYE, the coordinates of eyepoint.
!
!    Input/output, integer ( kind = 4 ) NVRT, the number of vertices on boundary of simple 
!    polygon.  On output, NVRT is increased by 1 if the closest vertex
!    is a new vertex.
!
!    Input/output, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertices 
!    of polygon in counter clockwise (or clockwise) order if eyepoint is 
!    interior (or blocked exterior);
!      (XC(0),YC(0)) = (XC(NVRT),YC(NVRT))
!    On output, the polygon vertices in same orientation
!    as input but with indices rotated and possibly
!      (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)) has been added.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 205.
!
  integer ( kind = 4 ) nvrt

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) irgt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) r
  real ( kind = 8 ) tol
  real ( kind = 8 ) xc(0:nvrt+1)
  real ( kind = 8 ) xeye
  real ( kind = 8 ) xint
  real ( kind = 8 ) xrgt
  real ( kind = 8 ) xt
  real ( kind = 8 ) yc(0:nvrt+1)
  real ( kind = 8 ) yeye
  real ( kind = 8 ) yeyemt
  real ( kind = 8 ) yeyept
  real ( kind = 8 ) yt

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  dy = 0.0D+00
  do i = 0, nvrt-1
    dy = max ( dy, abs ( yc(i+1) - yc(i) ) )
  end do

  yeyemt = yeye - tol * dy
  yeyept = yeye + tol * dy
  n = nvrt + 1
  irgt = n
  xrgt = 0.0D+00
!
!  Determine closest point on boundary which is to the right of
!  (XEYE,YEYE) and on the horizontal line through (XEYE,YEYE).
!  The closest point must be on an edge which intersects the
!  horizontal line and has (XEYE,YEYE) to its left.
!
  do i = 0, nvrt-1

    if ( yeyept < yc(i) .or. yc(i+1) < yeyemt ) then
      cycle
    end if

    if ( yc(i) < yeyemt .and. yeyept < yc(i+1) ) then
      xint = ( yeye - yc(i) ) * ( xc(i+1) - xc(i) ) &
        / ( yc(i+1) - yc(i) ) + xc(i)
      if ( xeye < xint ) then
         if ( xint < xrgt .or. irgt == n ) then
            irgt = -(i + 1)
            xrgt = xint
         end if
      end if
    else if ( yeyemt <= yc(i) .and. yeyept < yc(i+1) ) then
      if ( xeye < xc(i) ) then
         if ( xc(i) < xrgt .or. irgt == n ) then
            irgt = i
            xrgt = xc(i)
         end if
      end if
    else if ( yc(i) < yeyemt .and. yc(i+1) <= yeyept ) then
      if ( xeye < xc(i+1) ) then
         if ( xc(i+1) < xrgt .or. irgt == n ) then
            irgt = i + 1
            xrgt = xc(i+1)
         end if
      end if
    end if

  end do

  if ( irgt == n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROTIPG - Fatal error!'
    write ( *, '(a)' ) '  IRGT == N.'
    ierror = 205
    return
  end if

  if ( irgt == 0 .or. irgt == nvrt ) then
    return
  end if

  if ( irgt < 0 ) then
    irgt = -irgt
    do i = nvrt, irgt, -1
      xc(i+1) = xc(i)
      yc(i+1) = yc(i)
    end do
    xc(irgt) = xrgt
    yc(irgt) = yeye
    nvrt = nvrt + 1
  end if
!
!  Rotate the indices of the vertices so that (XC(IRGT),YC(IRGT))
!  becomes (XC(0),YC(0)). Compute A = GCD(NVRT,IRGT).
!
  a = nvrt
  b = irgt

  do

    r = mod ( a, b )
    a = b
    b = r

    if ( r <= 0 ) then
      exit
    end if

  end do

  m = nvrt / a - 1

  do i = 0, a-1

    xt = xc(i)
    yt = yc(i)
    k = i

    do j = 1, m
      l = k + irgt
      if ( nvrt <= l ) then
        l = l - nvrt
      end if
      xc(k) = xc(l)
      yc(k) = yc(l)
      k = l
    end do

    xc(k) = xt
    yc(k) = yt

  end do

  xc(nvrt) = xc(0)
  yc(nvrt) = yc(0)

  return
end
subroutine rotpg ( nvrt, xc, yc, i1, i2, ibot, costh, sinth )

!*****************************************************************************80
!
!! ROTPG rotates a convex polygon.
!
!  Purpose: 
!
!    Rotate convex polygon so that a line segment joining two
!    of its vertices is parallel to y-axis.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of
!    the convex polygon.
!
!    Input/output, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT).  The vertex 
!    coordinates in counter clockwise order;
!      (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!    On output, the rotated vertex coordinates; indices are
!    also rotated so that (XC(0),YC(0)) = (XC(NVRT),YC(NVRT))
!    is top vertex and (XC(IBOT),YC(IBOT)) is bottom vertex.
!
!    Input, integer ( kind = 4 ) I1, I2, the index of vertices of line segment; 
!    I1, I2 must be positive.
!
!    Output, integer ( kind = 4 ) IBOT, the index of bottom vertex.
!
!    Output, real ( kind = 8 ) COSTH, SINTH, the values COS(THETA) and 
!    SIN(THETA) where THETA in [-PI,PI] is the rotation angle.
!
  integer ( kind = 4 ) nvrt

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) costh
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) sinth
  real ( kind = 8 ) theta
  real ( kind = 8 ) tol
  real ( kind = 8 ) x0
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) y0
  real ( kind = 8 ) yc(0:nvrt)

  tol = 100.0D+00 * epsilon ( tol )

  itop = i1
  ibot = i2

  if ( yc(i1) == yc(i2) ) then
    if ( xc(i1) < xc(i2) ) then
      theta = -pi / 2.0D+00
    else
      theta = pi / 2.0D+00
    end if
  else
    if ( yc(i1) < yc(i2) ) then
      itop = i2
      ibot = i1
    end if
    theta = pi / 2.0D+00 &
      - atan2 ( yc(itop) - yc(ibot), xc(itop) - xc(ibot) )
  end if

  costh = cos(theta)
  sinth = sin(theta)

  do i = 1, nvrt
    x0 = xc(i)
    xc(i) = costh * x0 - sinth * yc(i)
    yc(i) = sinth * x0 + costh * yc(i)
  end do
!
!  Rotate indices.
!
  if ( itop /= nvrt ) then

    a = nvrt
    b = itop

    do

      r = mod ( a, b )
      a = b
      b = r

      if ( r <= 0 ) then
        exit
      end if

    end do

    m = nvrt / a - 1

    do i = 1, a

      x0 = xc(i)
      y0 = yc(i)
      k = i

      do j = 1, m
        l = k + itop
        if ( nvrt < l ) then
          l = l - nvrt
        end if
        xc(k) = xc(l)
        yc(k) = yc(l)
        k = l
      end do

      xc(k) = x0
      yc(k) = y0

    end do

    ibot = ibot - itop
    if ( ibot < 0 ) then
      ibot = ibot + nvrt
    end if

  end if

  xc(0) = xc(nvrt)
  yc(0) = yc(nvrt)

  return
end
subroutine sepmdf ( angtol, nvrt, xc, yc, arpoly, mean, mdftr, indpvl, &
  iang, i1, i2 )

!*****************************************************************************80
!
!! SEPMDF splits a polygon according to the mesh distribution function.
!
!  Purpose: 
!
!    Determine separator to split convex polygon into two
!    parts based on mesh distribution function.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter
!    (in radians).
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices in polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT),YC(0:NVRT), the coordinates of polygon
!    vertices in counter clockwise order, translated so that centroid is at 
!    origin; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, real ( kind = 8 ) ARPOLY, the area of polygon.
!
!    Input, real ( kind = 8 ) MEAN, the mean mdf value in polygon.
!
!    Input, real ( kind = 8 ) MDFTR(0:NVRT-1), the mean mdf value in each 
!    triangle of polygon; triangles are determined by polygon vertices 
!    and centroid.
!
!    Input, integer ( kind = 4 ) INDPVL(0:NVRT), the indices in PVL of vertices; 
!    INDPVL(I) = -K if (XC(I),YC(I)) is extra vertex inserted on edge from
!    K to PVL(SUCC,K).
!
!    Input, real ( kind = 8 ) IANG(1:*), the interior angle array.
!
!    Output, integer ( kind = 4 ) I1, I2, indices in range 0 to NVRT-1 of best separator
!    according to MDF and max-min angle criterion; I1 = -1
!    if no satisfactory separator is found.
!
  integer ( kind = 4 ) nvrt

  real ( kind = 8 ) angle
  real ( kind = 8 ) angtol
  real ( kind = 8 ) areatr
  real ( kind = 8 ) arpoly
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real ( kind = 8 ) iang(*)
  integer ( kind = 4 ) indpvl(0:nvrt)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real ( kind = 8 ) mdftr(0:nvrt-1)
  real ( kind = 8 ) mean
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sum2
  real ( kind = 8 ) tol
  integer ( kind = 4 ) v(2)
  integer ( kind = 4 ) w(2)
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) yc(0:nvrt)

  tol = 100.0D+00 * epsilon ( tol )
!
!  Determine triangle with highest mean mesh density; then determine
!  triangles adjacent to this triangle with mesh density at least MEAN
!  such that the area of these triangles is <= ARPOLY/2.
!  Note that twice the triangle area is computed.
!
  hi = 0
  do i = 1, nvrt-1
    if ( mdftr(hi) < mdftr(i) ) then
      hi = i
    end if
  end do

  sum2 = xc(hi) * yc(hi+1) - xc(hi+1) * yc(hi)

  l = hi - 1
  if ( l < 0 ) then
    l = nvrt - 1
  end if

  m = hi + 1
  if ( nvrt <= m ) then
    m = 0
  end if

  do

    if ( mdftr(m) <= mdftr(l) ) then
      i = l
    else
      i = m
    end if

    if ( mdftr(i) < mean ) then
      exit
    end if

    areatr = xc(i) * yc(i+1) - xc(i+1) * yc(i)
    sum2 = sum2 + areatr

    if ( arpoly < sum2 ) then
      exit
    end if

    if ( i == l ) then
      l = l - 1
      if ( l < 0 ) then
        l = nvrt - 1
      end if
    else
      m = m + 1
      if ( nvrt <= m ) then
        m = 0
      end if
    end if

  end do

  l = l + 1

  if ( nvrt <= l ) then
    l = 0
  end if
!
!  Interchange role of L and M depending on angle determined by
!  (XC(M),YC(M)), (0,0), and (XC(L),YC(L)).
!  Possible separators are L,M; L,M+1; L+1,M; L+1,M+1.
!
  if ( pi < angle ( xc(m), yc(m), 0.0D+00, 0.0D+00, xc(l), yc(l) ) ) then
    i = l
    l = m
    m = i
  end if

  v(1) = l
  v(2) = l - 1
  if ( v(2) < 0 ) then 
    v(2) = nvrt - 1
  end if

  w(1) = m
  w(2) = m + 1
  if ( nvrt <= w(2) ) then
    w(2) = 0
  end if

  call mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )

  return
end
subroutine sepshp ( angtol, nvrt, xc, yc, indpvl, iang, i1, i2, wk, ierror )

!*****************************************************************************80
!
!! SEPSHP splits a convex polygon according to shape.
!
!  Purpose: 
!
!    Determine separator to split convex polygon into two
!    parts based on shape (diameter) of polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter
!    (in radians).
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices in polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the coordinates of 
!    polygon vertices in counter clockwise order, translated so that 
!    centroid is at origin;
!    (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, integer ( kind = 4 ) INDPVL(0:NVRT), the indices in PVL of vertices;
!    INDPVL(I) = -K if (XC(I),YC(I)) is extra vertex inserted on edge from
!    K to PVL(SUCC,K).
!
!    Input, real ( kind = 8 ) IANG(1:*), the interior angle array.
!
!    Output, integer ( kind = 4 ) I1, I2, the indices in range 0 to NVRT-1 of best separator
!    according to shape and max-min angle criterion; I1 = -1
!    if no satisfactory separator is found.
!
!    Workspace, real ( kind = 8 ) WK(1:2*NVRT).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  For abnormal return,
!    IERROR is set to 200.
!
  integer ( kind = 4 ) nvrt

  real ( kind = 8 ) angtol
  real ( kind = 8 ) dist
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real ( kind = 8 ) iang(*)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indpvl(0:nvrt)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pimtol
  real ( kind = 8 ) tol
  integer ( kind = 4 ) v(2)
  integer ( kind = 4 ) w(2)
  real ( kind = 8 ) wk(2*nvrt)
  real ( kind = 8 ) xa
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) ya
  real ( kind = 8 ) yc(0:nvrt)

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Determine diameter of polygon. Possible separators endpoints (two
!  on each side of polygon) are nearest to perpendicular bisector of
!  diameter. (XA,YA) and (XA+DX,YA+DY) are on bisector. Distance to
!  bisector is proportional to two times triangle area.
!
  pimtol = pi - tol
  n = 0
  do i = 0, nvrt-1
    k = indpvl(i)
    if ( 0 < k ) then
      if ( iang(k) < pimtol ) then
         n = n + 1
         wk(n) = xc(i)
         wk(n+nvrt) = yc(i)
      end if
    end if
  end do

  call diam2 ( n, wk, wk(nvrt+1), i1, i2, dist, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  if ( i2 < i1 ) then
    i = i1
    i1 = i2
    i2 = i
  end if

  dx = wk(i2+nvrt) - wk(i1+nvrt)
  dy = wk(i1) - wk(i2)
  xa = 0.5D+00 * ( wk(i1) + wk(i2) - dx )
  ya = 0.5D+00 * ( wk(i1+nvrt) + wk(i2+nvrt) - dy )

  i = i1 - 1

20 continue

  if ( xc(i) == wk(i1) .and. yc(i) == wk(i1+nvrt) ) then
    i1 = i
  else
    i = i + 1
    go to 20
  end if

  i = max ( i2-1, i1+1 )

30 continue

  if ( xc(i) == wk(i2) .and. yc(i)  ==  wk(i2+nvrt) ) then
    i2 = i
  else
    i = i + 1
    go to 30
  end if

  i = i1 + 1

40 continue

  dist = dx * ( yc(i) - ya ) - dy * ( xc(i) - xa )

  if ( 0.0D+00 <= dist ) then
    v(1) = i - 1
    v(2) = i
  else
    i = i + 1
    go to 40
  end if

  i = i2 + 1

50 continue

  if ( nvrt <= i ) then
    i = 0
  end if

  dist = dx * ( yc(i) - ya ) - dy * ( xc(i) - xa )

  if ( dist <= 0.0D+00 ) then
    w(1) = i - 1
    w(2) = i
    if ( i <= 0 ) then
      w(1) = nvrt - 1
    end if
  else
    i = i + 1
    go to 50
  end if

  call mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )

  return
end
subroutine sfdwmf ( l, r, psi, indp, loch )

!*****************************************************************************80
!
!! SFDWMF sifts PSI(INDP(L)) down a heap.
!
!  Purpose: 
!
!    Sift PSI(INDP(L)) down heap which has maximum PSI value
!    at root of heap and is maintained by pointers in INDP.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the element of heap to be sifted down.
!
!    Input, integer ( kind = 4 ) R, the upper bound of heap.
!
!    Input, real ( kind = 8 ) PSI(1:*), the key values for heap.
!
!    Input/output, integer ( kind = 4 ) INDP(1:R), the indices of PSI which are 
!    maintained in heap.
!
!    Input/output, integer ( kind = 4 ) LOCH(1:*), the location of indices in heap 
!    (inverse of INDP).
!
  integer ( kind = 4 ) r

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indp(r)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) loch(*)
  real ( kind = 8 ) psi(*)
  real ( kind = 8 ) t

  i = l
  j = 2 * i
  k = indp(i)
  t = psi(k)

  do while ( j <= r )

    if ( j < r ) then
      if ( psi(indp(j)) < psi(indp(j+1)) ) then
        j = j + 1
      end if
    end if

    if ( psi(indp(j)) <= t ) then
      exit
    end if

    indp(i) = indp(j)
    loch(indp(i)) = i
    i = j
    j = 2 * i

  end do

  indp(i) = k
  loch(k) = i

  return
end
subroutine sfupmf ( r, psi, indp, loch )

!*****************************************************************************80
!
!! SFUPMF sifts PSI(INDP(R)) up a heap.
!
!  Purpose: 
!
!    Sift PSI(INDP(R)) up heap which has maximum PSI value
!    at root of heap and is maintained by pointers in INDP.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) R, the element of heap to be sifted up.
!
!    Input, real ( kind = 8 ) PSI(1:*), the key values for heap.
!
!    Input/output, integer ( kind = 4 ) INDP(1:R), the indices of PSI which are 
!    maintained in heap.
!
!    Input/output, integer ( kind = 4 ) LOCH(1:*), the location of indices in heap 
!    (inverse of INDP).
!
  integer ( kind = 4 ) r

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indp(r)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) loch(*)
  real ( kind = 8 ) psi(*)
  real ( kind = 8 ) t

  i = r
  j = int ( i / 2 )
  k = indp(i)
  t = psi(k)

  do

    if ( i <= 1 ) then
      exit
    end if

    if ( t <= psi(indp(j)) ) then
      exit
    end if

    indp(i) = indp(j)
    loch(indp(i)) = i
    i = j
    j = int ( i / 2 )

  end do

  indp(i) = k
  loch(k) = i

  return
end
subroutine shrnk2 ( nvrt, xc, yc, sdist, nshr, xs, ys, iedge, ierror )

!*****************************************************************************80
!
!! SHRNK2 shrinks a convex polygon.
!
!  Purpose: 
!
!    Shrink a convex polygon, with vertices given in counter clockwise
!    order and with all interior angles < PI, by distance SDIST(I)
!    for Ith edge, I = 0,...,NVRT-1.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of convex 
!    polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex coordinates 
!    in counter clockwise order;
!    (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, real ( kind = 8 ) SDIST(0:NVRT-1), the nonnegative shrink 
!    distances for edges.
!
!    Output, integer ( kind = 4 ) NSHR, the number of vertices on boundary of shrunken 
!    polygon; 0 if shrunken polygon is empty else 3 <= NSHR <= NVRT.
!
!    Output, real ( kind = 8 ) XS(0:NSHR), YS(0:NSHR), the coordinates of
!    shrunken polygon in counter clockwise order if NSHR is greater than 0; 
!    (XS(0),YS(0)) = (XS(NSHR),YS(NSHR)).
!
!    Output, integer ( kind = 4 ) IEDGE(0:NVRT), the indices of edges of shrunken polygon in
!    range from 0 to NVRT-1.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 202.
!
  integer ( kind = 4 ) nvrt

  real ( kind = 8 ) alpha
  logical first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge(0:nvrt)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nshr
  logical parall
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pi2
  real ( kind = 8 ) sdist(0:nvrt-1)
  real ( kind = 8 ) theta
  real ( kind = 8 ) tol
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) xs(0:nvrt)
  real ( kind = 8 ) yc(0:nvrt)
  real ( kind = 8 ) ys(0:nvrt)

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  pi2 = 2.0D+00 * pi
  alpha = atan2 ( yc(1)-yc(0), xc(1)-xc(0) )

  call xline ( xc(0), yc(0), xc(1), yc(1), xc(1), yc(1), xc(2), yc(2), &
    sdist(0), sdist(1), xs(1), ys(1), parall )

  if ( parall ) then
    ierror = 202
    go to 90
  end if

  iedge(0) = 0
  iedge(1) = 1
  i = 2
  j = 0
  nshr = 1
  first = .true.
!
!  First while loop processes edges subtending angle <= PI
!  with respect to first edge.
!
10 continue

  theta = atan2 ( yc(i+1)-yc(i), xc(i+1)-xc(i) ) - alpha
  if ( theta < 0.0D+00 ) then
    theta = theta + pi2
  end if

  if ( pi + tol < theta ) then
    go to 40
  end if

20 continue

  lr = lrline ( xs(nshr), ys(nshr), xc(i), yc(i), xc(i+1), yc(i+1), &
    sdist(i) )

  if ( lr < 0 ) then
    go to 30
  end if

  nshr = nshr - 1
  if ( 1 <= nshr ) then
    go to 20
  end if

30 continue

  if ( nshr < 1 .and. abs ( theta - pi ) <= tol ) then
    go to 90
  end if

  k = iedge(nshr)
  nshr = nshr + 1

  call xline ( xc(k), yc(k), xc(k+1), yc(k+1), xc(i), yc(i), xc(i+1), &
    yc(i+1), sdist(k), sdist(i), xs(nshr), ys(nshr), parall )

  if ( parall ) then
    ierror = 202
    go to 90
  end if

  iedge(nshr) = i
  i = i + 1
  go to 10
!
!  Second while loop processes remaining edges.
!
40 continue

  if ( first ) then
    first = .false.
    go to 50
  end if

  lr = lrline ( xs(j), ys(j), xc(i), yc(i), xc(i+1), yc(i+1), sdist(i) )

  if ( lr <= 0 ) then
    go to 70
  end if

50 continue

  if ( nshr <= j ) then
    go to 90
  end if

  lr = lrline ( xs(nshr), ys(nshr), xc(i), yc(i), xc(i+1), yc(i+1), sdist(i) )

  if ( 0 <= lr ) then
    nshr = nshr - 1
    go to 50
  end if

  k = iedge(nshr)
  nshr = nshr + 1

  call xline ( xc(k), yc(k), xc(k+1), yc(k+1), xc(i), yc(i), xc(i+1), &
    yc(i+1), sdist(k), sdist(i), xs(nshr), ys(nshr), parall )

  if ( parall ) then
    ierror = 202
    go to 90
  end if

  iedge(nshr) = i

60 continue

  lr = lrline ( xs(j+1), ys(j+1), xc(i), yc(i), xc(i+1), yc(i+1), sdist(i) )

  if ( 0 <= lr ) then
    j = j + 1
    go to 60
  end if

  k = iedge(j)
  call xline ( xc(k), yc(k), xc(k+1), yc(k+1), xc(i), yc(i), xc(i+1), &
    yc(i+1), sdist(k), sdist(i), xs(j), ys(j), parall )

  if ( parall ) then
    ierror = 202
    go to 90
  end if

  xs(nshr+1) = xs(j)
  ys(nshr+1) = ys(j)
  iedge(nshr+1) = iedge(j)

70 continue

  i = i + 1

  if ( i < nvrt ) then
    go to 40
  end if

  if ( 0 < j ) then
    do i = 0, nshr+1-j
      xs(i) = xs(i+j)
      ys(i) = ys(i+j)
      iedge(i) = iedge(i+j)
    end do
  end if

  nshr = nshr + 1 - j
  return

90 continue

  nshr = 0

  return
end
subroutine spdec2 ( angspc, angtol, nvc, npolg, nvert, nhole, nhola, maxvc,  &
  maxhv, maxpv, maxiw, maxwk, holv, vcl, regnum, hvl, pvl, iang, iwk, &
  wk, ierror )

!*****************************************************************************80
!
!! SPDEC2 decomposes a polygonal region with holes into simple polygons.
!
!  Discussion: 
!
!    Decompose general polygonal region with interfaces and
!    holes into simple polygons using vertex coordinate list,
!    head vertex list, and polygon vertex list data structures.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in radians 
!    used in controlling vertices to be considered as an endpoint of 
!    a separator.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter in radians 
!    used in accepting separator(s).
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates 
!    or positions used in the VCL array.
!
!    Input/output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions or 
!    positions used in HVL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of polygon vertices or positions 
!    used in PVL array.
!
!    Input, integer ( kind = 4 ) NHOLE, the number of holes and hole interfaces.
!
!    Input, integer ( kind = 4 ) NHOLA, the number of 'attached' holes; these holes are 
!    attached to the outer boundary of a subregion through vertices
!    or cut interfaces and have their edges in consecutive order on 
!    the boundary.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array, should 
!    be at least the number of vertex coordinates required for decomposition.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM arrays, 
!    should be at least the number of polygons required for decomposition.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays; 
!    should be at least the number of polygon vertices required for
!    decomposition.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should be 
!    about 3 times maximum number of vertices in any polygon.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be 
!    about 5 times maximum number of vertices in any polygon.
!
!    Input, integer ( kind = 4 ) HOLV(1:NHOLE*2+NHOLA), the indices in PVL of bottom or top 
!    vertex of holes; first (next) NHOLE entries are for top (bottom)
!    vertices of holes and hole interfaces, with top (bottom)
!    vertices sorted in decreasing (increasing) lexicographic
!    (y,x) order of coord; last NHOLA entries are for attached
!    holes; if bottom vertex of attached hole is a simple
!    vertex of boundary curve containing the hole then entry
!    contains index of bottom vertex otherwise entry contains
!    index of top vertex (which is simple).
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT), 
!    the polygon vertex list and interior angles; see routine DSPGDC for more 
!    details.  Note that the data structures should be as output from routines
!    DSMCPR or DSPGDC.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 4, 5, 6, 7, 206 to 210, 212, 218, or 219.
!
  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  logical ci
  logical cj
  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ) holv(*)
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) piptol
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  integer ( kind = 4 ), parameter :: succ = 3
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vr
  integer ( kind = 4 ) w1
  integer ( kind = 4 ) w2
  real ( kind = 8 ) wk(maxwk)

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  For each simple hole, find cut edge from top vertex of hole to
!  a point on the outer boundary above top vertex, and update
!  VCL, HVL, PVL, IANG.
!
  piptol = pi + tol

  do i = 1, nhole

    call jnhole ( holv(i), angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
      maxwk, vcl, hvl, pvl, iang, iwk, wk, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPDEC2 - Fatal error!'
      write ( *, '(a)' ) '  JNHOLE returned an error condition.'
      return
    end if

  end do
!
!  Resolve remaining vertices in HOLV array if they are reflex
!  vertices. These vertices may no longer be reflex if they are the
!  endpoint of a cut edge from the top vertex of another hole or
!  of a previous separator.
!
  do i = nhole+1, nhole+nhole+nhola

    vr = holv(i)

    if ( piptol < iang(vr) ) then

      call resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
        maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      call insed2 ( vr, w1, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
        pvl, iang, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      if ( 0 < w2 ) then
        call insed2 ( vr, w2, npolg, nvert, maxhv, maxpv, &
          vcl, regnum, hvl, pvl, iang, ierror )
        if ( ierror /= 0 ) then
          return
        end if
      end if

    end if

  end do

  if ( nhola == 0 ) then
    return
  end if
!
!  Check that polygons are simple. If polygon is simply-connected and
!  not simple then find a simple reflex vertex in polygon to resolve.
!
  p = 1

30 continue

  if ( npolg < p ) then
    return
  end if

  i = hvl(p)

  do

    if ( pvl(polg,pvl(edgv,i)) == p ) then
      go to 50
    end if

    i = pvl(succ,i)

    if ( i == hvl(p) ) then
      exit
    end if

  end do

  p = p + 1
  go to 30

50 continue

  ci = .true.

  do

    j = pvl(succ,i)
    cj = ( pvl(polg,pvl(edgv,j)) == p )

    if ( .not. ci .and. .not. cj .and. piptol < iang(j) ) then
      exit
    end if

    i = j
    ci = cj

  end do

  vr = j
  call resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
    maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call insed2 ( vr, w1, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
    pvl, iang, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  if ( 0 < w2 ) then

    call insed2 ( vr, w2, npolg, nvert, maxhv, maxpv, &
      vcl, regnum, hvl, pvl, iang, ierror )

    if ( ierror /= 0 ) then
      return
    end if

  end if

  go to 30

end
subroutine swapec ( i, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierror )

!*****************************************************************************80
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to triangulation.
!
!  Modified:
!
!    10 August 2003
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index in VCL of the new vertex.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input, integer ( kind = 4 ) MAXST, the maximum size available for the STACK array.
!
!    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are the
!    triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, real ( kind = 8 ) VCL(2,*), the coordinates of the vertices.
!
!    Input/output, integer ( kind = 4 ) TIL(3,*), the triangle incidence list.  
!    May be updated on output because of swaps.
!
!    Input/output, integer ( kind = 4 ) TNBR(3,*), the triangle neighbor list; negative 
!    values are used for links of the counter-clockwise linked list of boundary 
!    edges; May be updated on output because of swaps.
!
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(1:MAXST); on input, entries 1 through TOP
!    contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer ( kind = 4 ) IERROR is set to 8 for abnormal return.
!
  integer ( kind = 4 ) maxst

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  integer ( kind = 4 ) c
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fm1
  integer ( kind = 4 ) fp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(maxst)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real ( kind = 8 ) vcl(2,*)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Determine whether the triangles in the stack are Delaunay.
!  Ifnot, swap the diagonal edge of the convex quadrilateral.
!
  ierror = 0
  x = vcl(1,i)
  y = vcl(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( til(1,t) == i ) then
      e = 2
      b = til(3,t)
    else if ( til(2,t) == i ) then
      e = 3
      b = til(1,t)
    else
      e = 1
      b = til(2,t)
    end if

    a = til(e,t)
    u = tnbr(e,t)

    if ( tnbr(1,u) == t ) then
      f = 1
      c = til(3,u)
    else if ( tnbr(2,u) == t ) then
      f = 2
      c = til(1,u)
    else
      f = 3
      c = til(2,u)
    end if

    swap = diaedg ( x, y, vcl(1,a), vcl(2,a), vcl(1,c), vcl(2,c), &
      vcl(1,b), vcl(2,b) )

    if ( swap == 1 ) then

      em1 = i4_wrap ( e - 1, 1, 3 )
      ep1 = i4_wrap ( e + 1, 1, 3 )
      fm1 = i4_wrap ( f - 1, 1, 3 )
      fp1 = i4_wrap ( f + 1, 1, 3 )
      
      til(ep1,t) = c
      til(fp1,u) = i
      r = tnbr(ep1,t)
      s = tnbr(fp1,u)
      tnbr(ep1,t) = u
      tnbr(fp1,u) = t
      tnbr(e,t) = s
      tnbr(f,u) = r

      if ( 0 < tnbr(fm1,u) ) then
        top = top + 1
        stack(top) = u
      end if

      if ( 0 < s ) then

        if ( tnbr(1,s) == u ) then
          tnbr(1,s) = t
        else if ( tnbr(2,s) == u ) then
          tnbr(2,s) = t
        else
          tnbr(3,s) = t
        end if

        top = top + 1

        if ( maxst < top ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPEC - Fatal error!'
          write ( *, '(a)' ) '  Ran out of stack storage.'
          write ( *, '(a,i6)' ) '  MAXST = ', maxst
          write ( *, '(a,i6)' ) '  TOP =   ', top
          ierror = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( 0 < tnbr(ee,tt) )

          tt = tnbr(ee,tt)

          if ( til(1,tt) == a ) then
            ee = 3
          else if ( til(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tnbr(ee,tt) = l

      end if

      if ( 0 < r ) then

        if ( tnbr(1,r) == t ) then
          tnbr(1,r) = u
        else if ( tnbr(2,r) == t ) then
          tnbr(2,r) = u
        else
          tnbr(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( 0 < tnbr(ee,tt) )

          tt = tnbr(ee,tt)

          if ( til(1,tt) == b ) then
            ee = 3
          else if ( til(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tnbr(ee,tt) = l

      end if

    end if

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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
subroutine tmerge ( inter, nbl, ncr, chbl, chcr, ldv, vcl, til, tedg, &
  ierror )

!*****************************************************************************80
!
!! TMERGE forms triangles near the boundary by merging vertex chains.
!
!  Purpose: 
!
!    Form triangles in strip near boundary of polygon or
!    inside polygon by merging two chains of vertices.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical INTER, is .TRUE. iff at least one interior mesh vertex.
!
!    Input, integer ( kind = 4 ) NBL, the number of vertices on boundary cycle if INTER,
!    otherwise on left boundary chain.
!
!    Input, integer ( kind = 4 ) NCR, the number of vertices on closed walk if INTER,
!    otherwise on right boundary chain.
!
!    Input, integer ( kind = 4 ) CHBL(0:NBL), the indices in VCL of vertices on boundary
!    cycle or left boundary chain; if INTER, CHBL(NBL) = CHBL(0).
!
!    Input, integer ( kind = 4 ) CHCR(0:NCR), the indices in VCL of vertices on closed walk
!    or right boundary chain; if INTER, CHCR(NCR) = CHCR(0),
!    otherwise CHCR(0) is not referenced.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Output, integer ( kind = 4 ) TIL(1:3,1:NT), the triangle incidence list, where NT =
!    NBL + NCR - K where K = 0 if INTER, else K = 2.
!
!    Output, integer ( kind = 4 ) TEDG(1:3,1:NT), the TEDG(J,I) refers to edge with vertices
!    TIL(J:J+1,I) and contains index of merge edge or NBL+NCR+1 for edge of 
!    chains.  Note: It is assumed there is enough space in 2 arrays.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 230.
!
  implicit none

  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) nbl
  integer ( kind = 4 ) ncr

  integer ( kind = 4 ) chbl(0:nbl)
  integer ( kind = 4 ) chcr(0:ncr)
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibndry
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) in
  logical inter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lri
  integer ( kind = 4 ) lrip1
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) tedg(3,nbl+ncr)
  integer ( kind = 4 ) til(3,nbl+ncr)
  real ( kind = 8 ) vcl(ldv,*)
  real ( kind = 8 ) xi
  real ( kind = 8 ) xip1
  real ( kind = 8 ) xj
  real ( kind = 8 ) xjp1
  real ( kind = 8 ) yi
  real ( kind = 8 ) yip1
  real ( kind = 8 ) yj
  real ( kind = 8 ) yjp1

  ierror = 0
  ibndry = nbl + ncr + 1
  nt = 0

  if ( inter ) then

    nl = nbl
    nr = ncr
    i = 0
    j = 0

  else

    call mtredg ( .true., chbl(1), chcr(1), chbl(0), ibndry, nt, til, tedg )

    tedg(2,1) = ibndry
    if ( nbl + ncr <= 3 ) then
      return
    end if

    nl = nbl - 1
    nr = ncr - 1
    i = 1
    j = 1
    lri = 1
    lrip1 = 1

  end if
!
!  Main while loop for determining next triangle and edge.
!
10 continue

  if ( nl <= i .or. nr <= j ) then
    go to 20
  end if

   xi = vcl(1,chbl(i))
   yi = vcl(2,chbl(i))
   xip1 = vcl(1,chbl(i+1))
   yip1 = vcl(2,chbl(i+1))
   xj = vcl(1,chcr(j))
   yj = vcl(2,chcr(j))
   xjp1 = vcl(1,chcr(j+1))
   yjp1 = vcl(2,chcr(j+1))
   in = diaedg ( xjp1, yjp1, xj, yj, xi, yi, xip1, yip1 )

   if ( inter ) then
     lri = lrline ( xi, yi, xj, yj, xjp1, yjp1, 0.0D+00 )
     lrip1 = lrline ( xip1, yip1, xj, yj, xjp1, yjp1, 0.0D+00 )
   end if

   if ( in <= 0 .or. lri <= 0 .and. lrip1 <= 0 ) then

     call mtredg ( .true., chbl(i+1), chcr(j), chbl(i), ibndry, nt, til, tedg )
     i = i + 1

   else

     call mtredg ( .false., chbl(i), chcr(j+1), chcr(j), ibndry, nt, til, tedg )
     j = j + 1

   end if

  go to 10
!
!  Add remaining triangles at end of strip or bottom of polygon.
!
20 continue

  if ( i < nl ) then

    if ( .not. inter .and. j == nr ) then
      nl = nl + 1
    end if

    do

      call mtredg ( .true., chbl(i+1), chcr(j), chbl(i), ibndry, nt, til, tedg )
      i = i + 1

      if ( nl <= i ) then
        exit
      end if

    end do

  else
!
!  J < NR .OR. I = NL = J = NR = 1
!
    if ( .not. inter .and. i == nl ) then
      nr = nr + 1
    end if

40  continue

    call mtredg ( .false., chbl(i), chcr(j+1), chcr(j), ibndry, nt, til, tedg )

    if ( inter ) then

      lri = lrline ( vcl(1,chbl(i)), vcl(2,chbl(i)), &
        vcl(1,chcr(j+1)), vcl(2,chcr(j+1)), vcl(1,chcr(j)), &
        vcl(2,chcr(j)), 0.0D+00 )

      if ( 0 <= lri ) then
        ierror = 230
        return
      end if

    end if

    j = j + 1

    if ( j < nr ) then
      go to 40
    end if

  end if

  if ( inter ) then
    if ( tedg(2,1) == 0 ) then
      tedg(2,1) = nbl + ncr
    else
      tedg(3,1) = nbl + ncr
    end if
  end if

  return
end
subroutine triangulation_plot_eps ( file_name, g_num, g_xy, tri_num, nod_tri )

!*****************************************************************************80
!
!! TRIANGULATION_PLOT_EPS plots a triangulation of a pointset.
!
!  Discussion:
!
!    The triangulation is most usually a Delaunay triangulation,
!    but this is not necessary.
!
!    The data can be generated by calling RTRIS2, but this is not
!    necessary.
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) G_NUM, the number of points.
!
!    Input, real ( kind = 8 ) G_XY(2,G_NUM), the coordinates of the points.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) NOD_TRI(3,TRI_NUM), lists, for each triangle,
!    the indices of the points that form the vertices of the triangle.
!
  implicit none

  integer ( kind = 4 ) g_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) e
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) g
  real ( kind = 8 ) g_xy(2,g_num)
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nod_tri(3,tri_num)
  integer ( kind = 4 ) t
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) x_ps
  integer ( kind = 4 ) :: x_ps_max = 576
  integer ( kind = 4 ) :: x_ps_max_clip = 594
  integer ( kind = 4 ) :: x_ps_min = 36
  integer ( kind = 4 ) :: x_ps_min_clip = 18
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer ( kind = 4 ) y_ps
  integer ( kind = 4 ) :: y_ps_max = 666
  integer ( kind = 4 ) :: y_ps_max_clip = 684
  integer ( kind = 4 ) :: y_ps_min = 126
  integer ( kind = 4 ) :: y_ps_min_clip = 108

  x_max = maxval ( g_xy(1,1:g_num) )
  x_min = minval ( g_xy(1,1:g_num) )
  y_max = maxval ( g_xy(2,1:g_num) )
  y_min = minval ( g_xy(2,1:g_num) )
!
!  Plot the Delaunay triangulation.
!
! call get_unit ( file_unit )
  file_unit = 1

  open ( unit = file_unit, file = file_name, status = 'replace' )

  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: triangulation_plot_eps.f90'
  write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( file_unit, '(a)' ) '%%Pages: 1'
  write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%Bounding Box: ', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( file_unit, '(a)' ) '%%EndComments'
  write ( file_unit, '(a)' ) '%%BeginProlog'
  write ( file_unit, '(a)' ) '/inch {72 mul} def'
  write ( file_unit, '(a)' ) '%%EndProlog'
  write ( file_unit, '(a)' ) '%%Page: 1 1'
  write ( file_unit, '(a)' ) 'save'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
  write ( file_unit, '(a)' ) 'stroke'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to black.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the font and its size.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '/Times-Roman findfont'
  write ( file_unit, '(a)' ) '0.50 inch scalefont'
  write ( file_unit, '(a)' ) 'setfont'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Print a title.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '210  702  moveto'
  write ( file_unit, '(a)' ) '(Triangulation)  show'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a)' ) 'clip newpath'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to green.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.750  0.150 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw the nodes.'
  write ( file_unit, '(a)' ) '%'

  do g = 1, g_num
    x_ps = int ( &
      ( ( x_max - g_xy(1,g) ) * real ( x_ps_min, kind = 8 ) &
      + ( g_xy(1,g) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
      / ( x_max - x_min ) )
    y_ps = int ( &
      ( ( y_max - g_xy(2,g) ) * real ( y_ps_min, kind = 8 ) &
      + ( g_xy(2,g) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
      / ( y_max - y_min ) )
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 'newpath ', x_ps, y_ps, &
      ' 5 0 360 arc closepath fill'
  end do

  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to red.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw the triangles.'
  write ( file_unit, '(a)' ) '%'

  do t = 1, tri_num

    write ( file_unit, '(a)' ) 'newpath'

    do j = 1, 4

      e = i4_wrap ( j, 1, 3 )

      k = nod_tri(e,t)

      x_ps = int ( &
        ( ( x_max - g_xy(1,k) ) * real ( x_ps_min, kind = 8 ) &
        + ( g_xy(1,k) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max - x_min ) )

      y_ps = int ( &
        ( ( y_max - g_xy(2,k) ) * real ( y_ps_min, kind = 8 ) &
        + ( g_xy(2,k) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max - y_min ) )

      if ( j == 1 ) then
        write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' moveto'
      else
        write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'
      end if

    end do

    write ( file_unit, '(a)' ) 'stroke'

  end do

  write ( file_unit, '(a)' ) 'restore  showpage'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  End of page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%%Trailer'
  write ( file_unit, '(a)' ) '%%EOF'
  close ( unit = file_unit )

  return
end
subroutine trinbr ( nvc, ntri, til, tnbr, htsiz, maxedg, ht, edge, ierror )

!*****************************************************************************80
!
!! TRINBR determines the neighboring triangles of every triangle.
!
!  Purpose: 
!
!    Determine the neighboring triangle, if any, along each edge
!    of every triangle of triangulation of polygonal region.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVC, the number of vertices in triangulation.
!
!    Input, integer ( kind = 4 ) NTRI, the number of triangles in triangulation.
!
!    Input, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list; TIL(1:3,I) 
!    contains indices in VCL of 3 vertices of Ith triangle in counter 
!    clockwise order.
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT; should be a prime number
!    which is about NB where NB is number of boundary edges.
!
!    Input, integer ( kind = 4 ) MAXEDG, the maximum size available for EDGE array; should
!    be about 2*NB.
!
!    Output, integer ( kind = 4 ) TNBR(1:3,1:NTRI), the triangle neighbor list; positive 
!    elements are indices of TIL; zero elements indicate boundary edges.
!
!    Workspace, integer HT(0:HTSIZ-1), EDGE(1:4,1:MAXEDG), the hash table 
!    and edge records used to determine matching occurrences of triangle edges
!    by calling routine EDGHT.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 1.
!
  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) maxedg
  integer ( kind = 4 ) ntri

  integer ( kind = 4 ) e
  integer ( kind = 4 ) edge(4,maxedg)
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save, dimension ( 3 ) :: jp1 = (/ 2, 3, 1 /)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,ntri)
  integer ( kind = 4 ) tnbr(3,ntri)
  integer ( kind = 4 ) w

  ierror = 0
  hdfree = 0
  last = 0
  ht(0:htsiz-1) = 0

  do i = 1, ntri

    i3 = i * 3

    do j = 1, 3

      call edght ( til(j,i), til(jp1(j),i), i3+j-1, nvc, htsiz, maxedg, &
        hdfree, last, ht, edge, w, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      if ( 0 < w ) then
        t = w / 3
        e = mod ( w, 3 ) + 1
        tnbr(e,t) = i
        tnbr(j,i) = t
      end if

    end do

  end do

  do while ( hdfree /= 0 )
    edge(1,hdfree) = 0
    hdfree = edge(4,hdfree)
  end do

  do i = 1, last

    if ( edge(1,i) /= 0 ) then
      t = edge(3,i)/3
      e = mod ( edge(3,i), 3 ) + 1
      tnbr(e,t) = 0
    end if

  end do

  return
end
subroutine tripr2 ( nvc, npolg, nvert, maxvc, maxti, maxiw, maxwk, h, vcl,  &
  hvl, pvl, iang, ntri, til, vstart, vnum, tstart, iwk, wk, ierror )

!*****************************************************************************80
!
!! TRIPR2 generates triangles inside each convex polygon of a decomposition.
!
!  Purpose: 
!
!    Generate mesh vertices and triangles inside each convex
!    polygon of decomposition according to mesh spacings in H array
!    to get a triangulation of a polygonal region.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or positions
!    used in VCL array.
!
!    Input, integer ( kind = 4 ) NPOLG, the number of polygonal subregions or positions 
!    used in HVL array.
!
!    Input, integer ( kind = 4 ) NVERT, the number of polygon vertices or positions used 
!    in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array, should 
!    be greater than or equal to the number of mesh vertices in the 
!    triangulation of the region.
!
!    Input, integer ( kind = 4 ) MAXTI, the maximum size available for TIL array, should 
!    be greater than or equal to the number of triangles in the triangulation 
!    of region.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array, should 
!    be greater than or equal to 5*(NBC+NCW)+2, where NBC is maximum number 
!    of mesh edges on boundary of a polygon, NCW is maximum number of edges
!    on boundary of interior triangulation.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array, should 
!    be at least 5*NVRT+4 where NVRT is max no. of vertices in a polygon.
!
!    Input, real ( kind = 8 ) H(1:NPOLG), the mesh spacings for the polygons 
!    of the decomposition.
!
!    Input/output, real ( kind = 8 ) VCL(2,MAXVC), the vertex coordinates.
!
!    Input, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT), the 
!    polygon vertex list and interior angles; see routine DSPGDC for 
!    more details.
!
!    Output, integer ( kind = 4 ) NTRI, the number of triangles in triangulation of region.
!
!    Output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list; TIL(1:3,I) 
!    contains indices in VCL of 3 vertices of Ith triangle in counter
!    clockwise order.
!
!    Output, integer ( kind = 4 ) VSTART(1:NVERT), the start location in VCL for mesh 
!    vertices on each edge in PVL if there are any, else 0.
!
!    Output, integer ( kind = 4 ) VNUM(1:NVERT), the number of mesh vertices on interior
!    of each edge in PVL; entry is negated if mesh vertices are
!    listed in backward order in VCL.
!
!    Output, integer ( kind = 4 ) TSTART(1:NPOLG), the start location in TIL of triangles in
!    each polygon; TIL(1:3,I) for I=TSTRT(K),...,TSTRT(K+1)-1
!    are the triangles in the K-th polygon.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 6, 7, 9, 10, 200, 202, 230, or 231
!
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxti
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvert

  integer ( kind = 4 ) bndcyc
  real ( kind = 8 ) h(npolg)
  integer ( kind = 4 ) hvl(npolg)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(nvert)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) nbc
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvrt
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pimtol
  integer ( kind = 4 ) pvl(4,nvert)
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) til(3,maxti)
  real ( kind = 8 ) tol
  integer ( kind = 4 ) tstart(npolg)
  real ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vnum(nvert)
  integer ( kind = 4 ) vstart(nvert)
  real ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xc
  integer ( kind = 4 ) yc

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  ntri = 0
  pimtol = pi - tol


  call bedgmv ( nvc, npolg, nvert, maxvc, h, vcl, hvl, pvl, vstart, vnum, &
    ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPR2 - Fatal error!'
    write ( *, '(a)' ) '  Error return from BEDGMV.'
    return
  end if

  do k = 1, npolg

    nvrt = 0
    nbc = 0
    i = hvl(k)

    do

      if ( iang(i) < pimtol ) then
        nvrt = nvrt + 1
      endif

      nbc = nbc + 1 + abs ( vnum(i))
      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    if ( maxiw < nbc + 1 ) then
      ierror = 6
      return
    end if

    if ( maxwk < 2*nvrt + 2  ) then
      ierror = 7
      return
    end if

    xc = 1
    yc = xc + nvrt + 1
    bndcyc = 1

    do

      j = pvl(loc,i)

      if ( iang(i) < pimtol ) then
        wk(xc) = vcl(1,j)
        wk(yc) = vcl(2,j)
        xc = xc + 1
        yc = yc + 1
      end if

      iwk(bndcyc) = j
      bndcyc = bndcyc + 1

      if ( 0 <= vnum(i) ) then
        do j = vstart(i), vstart(i)+vnum(i)-1
          iwk(bndcyc) = j
          bndcyc = bndcyc + 1
        end do
      else
        do j = vstart(i)-vnum(i)-1, vstart(i),-1
          iwk(bndcyc) = j
          bndcyc = bndcyc + 1
        end do
      end if

      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    wk(xc) = wk(1)
    wk(yc) = wk(nvrt+2)
    iwk(bndcyc) = iwk(1)
    xc = 1
    yc = xc + nvrt + 1
    bndcyc = 1
    tstart(k) = ntri + 1
!
!  Generate a Delaunay triangulation inside a convex polygon.
!
    call trpolg ( nvrt, wk(xc), wk(yc), h(k), nbc, iwk(bndcyc), 2, nvc, ntri, &
      maxvc, maxti, maxiw-nbc-1, maxwk-2*nvrt-2, vcl, til, iwk(nbc+2), &
      wk(2*nvrt+3), ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIPR2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from TRPOLG.'
      write ( *, '(a,i6)' ) '  IERROR = ', ierror
      return
    end if

  end do

  return
end
subroutine trisiz ( ntrid, npolg, hvl, pvl, area, psi, h, indp, loch )

!*****************************************************************************80
!
!! TRISIZ smooths the mean mesh distribution function.
!
!  Purpose: 
!
!    Smooth PSI (mean mesh distribution function) values using
!    heap so that they differ by a factor of at most 4 in adjacent
!    polygons and then compute triangle sizes for each polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTRID, the desired number of triangles in mesh.
!
!    Input, integer ( kind = 4 ) NPOLG, the number of polygons or positions used in
!    HVL array.
!
!    Input, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), the polygon vertex list.
!
!    Input, real ( kind = 8 ) AREA(1:NPOLG), the area of convex polygons 
!    in decomposition.
!
!    Input/output, real ( kind = 8 ) PSI(1:NPOLG), the mean mdf values in 
!    the convex polygons.
!
!    Output, real ( kind = 8 ) H(1:NPOLG), the triangle size for
!    convex polygons.
!
!    Workspace, integer INDP(1:NPOLG), the indices of polygon or PSI which 
!    are maintained in heap according to PSI values.
!
!    Workspace, integer LOCH(1:NPOLG), the location of polygon indices in heap.
!
  integer ( kind = 4 ) npolg

  real ( kind = 8 ) area(npolg)
  integer ( kind = 4 ), parameter :: edgv = 4
  real ( kind = 8 ) factor
  real ( kind = 8 ) h(npolg)
  integer ( kind = 4 ) hvl(npolg)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indp(npolg)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) loch(npolg)
  integer ( kind = 4 ) ntrid
  integer ( kind = 4 ), parameter :: polg = 2
  real ( kind = 8 ) psi(npolg)
  integer ( kind = 4 ) pvl(4,*)
  integer ( kind = 4 ) r
  integer ( kind = 4 ), parameter :: succ = 3

  factor = 0.25D+00

  call i4vec_indicator ( npolg, indp )
  call i4vec_indicator ( npolg, loch )

  k = int ( npolg / 2 )

  do l = k, 1, -1
    call sfdwmf ( l, npolg, psi, indp, loch )
  end do

  do r = npolg, 2, -1

    j = indp(1)
    indp(1) = indp(r)
    loch(indp(1)) = 1
    call sfdwmf ( 1, r-1, psi, indp, loch )
    i = hvl(j)

    do

      k = pvl(edgv,i)

      if ( 0 < k ) then

        k = pvl(polg,k)

        if ( psi(k) < psi(j) * factor ) then
          psi(k) = psi(j) * factor
          call sfupmf ( loch(k), psi, indp, loch )
        end if

      end if

      i = pvl(succ,i)

      if ( i == hvl(j) ) then
        exit
      end if

    end do

  end do

  psi(1:npolg) = psi(1:npolg) / dot_product ( psi(1:npolg), area(1:npolg) )

  h(1:npolg) = sqrt ( 2.0D+00  / ( real ( ntrid, kind = 8 ) * psi(1:npolg) ) )

  return
end
subroutine trpolg ( nvrt, xc, yc, h, nbc, bndcyc, ldv, nvc, ntri, maxvc,  &
  maxti, maxiw, maxwk, vcl, til, iwk, wk, ierror )

!*****************************************************************************80
!
!! TRPOLG generates a Delaunay triangular mesh inside a convex polygon.
!
!  Discussion:
!
!    A quasi-uniform grid of spacing H is used.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of
!    convex polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex coordinates 
!    in counter clockwise order; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)); it is 
!    assumed that all interior angles are < PI.
!
!    Input, real ( kind = 8 ) H, the spacing of mesh vertices in polygon.
!
!    Input, integer ( kind = 4 ) NBC, the size of BNDCYC.
!
!    Input/output, integer ( kind = 4 ) BNDCYC(0:NBC), the indices in VCL of mesh 
!    vertices of boundary cycle; BNDCYC(0) = BNDCYC(NBC); 
!    contains (XC(I),YC(I)).
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of coordinates or positions used 
!    in VCL array.
!
!    Input/output, integer ( kind = 4 ) NTRI, the number of triangles or positions used 
!    in TIL.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXTI, the maximum size available for TIL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array, should 
!    be at least 6*(1 + INT(DIAM/H)) + 4*(NBC + NCW) where DIAM is
!    diameter of polygon, NCW is number of edges on boundary
!    of interior triangulation.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array, should 
!    be at least 3*NVRT+2.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 6, 7, 9, 10, 200, 202, 230, or 231.
!
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxti
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nbc
  integer ( kind = 4 ) nvrt

  integer ( kind = 4 ) bndcyc(0:nbc)
  real ( kind = 8 ) costh
  integer ( kind = 4 ) cwalk
  real ( kind = 8 ) dist
  real ( kind = 8 ) h
  real ( kind = 8 ) hs
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind
  logical inter
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) maxcw
  integer ( kind = 4 ) mbc
  integer ( kind = 4 ) ncw
  integer ( kind = 4 ) nshr
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) sdist
  real ( kind = 8 ) sinth
  real ( kind = 8 ) smdist
  integer ( kind = 4 ) sptr
  integer ( kind = 4 ) tedg
  integer ( kind = 4 ) til(3,maxti)
  real ( kind = 8 ) vcl(ldv,maxvc)
  real ( kind = 8 ) wk(maxwk)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) xi
  integer ( kind = 4 ) xs
  real ( kind = 8 ) y0
  real ( kind = 8 ) yc(0:nvrt)
  real ( kind = 8 ) yi
  real ( kind = 8 ) yr
  integer ( kind = 4 ) ys

  ierror = 0

  if ( maxiw < nvrt + 1 ) then
    ierror = 6
    return
  end if

  if ( maxwk < 3*nvrt + 2 ) then
    ierror = 7
    return
  end if

  xs = 1
  ys = xs + nvrt + 1
  sdist = ys + nvrt + 1
  iedge = 1
  hs = h / sqrt ( 2.0D+00 )
  wk(sdist:sdist+nvrt-1) = hs

  call shrnk2 ( nvrt, xc, yc, wk(sdist), nshr, wk(xs), wk(ys), iwk(iedge), &
    ierror )

  if ( ierror /= 0 ) then
    return
  end if

  inter = ( 0 < nshr )

  if ( inter ) then

    call diam2 ( nshr, wk(xs+1), wk(ys+1), i1, i2, dist, ierror )

    if ( ierror /= 0 ) then
      return
    end if

    call rotpg ( nshr, wk(xs), wk(ys), i1, i2, ibot, costh, sinth )

    maxcw = 6 * ( 1 + int ( ( wk(ys) - wk(ys+ibot) ) / h ) )

    if ( maxiw < maxcw + 1 ) then
      ierror = 6
      return
    end if

    cwalk = 1

    call inttri ( nshr, wk(xs), wk(ys), h, ibot, costh, sinth, ldv, nvc, ntri, &
      maxvc, maxti, maxcw, vcl, til, ncw, iwk(cwalk), ierror )

    if ( ierror /= 0 ) then
      return
    end if
!
!  Determine the mesh vertex which should be moved to front of
!  BNDCYC - closest to CWALK(0) and also with y-coordinate
!  greater than that of CWALK(0) when rotated if 0 < NCW.
!
    x0 = vcl(1,iwk(cwalk))
    y0 = vcl(2,iwk(cwalk))

    if ( 0 < ncw ) then
      yr = sinth * x0 + costh * y0
    end if

    smdist = 100000.0D+00 * h**2

    do i = 0, nbc-1

      xi = vcl(1,bndcyc(i))
      yi = vcl(2,bndcyc(i))

      if ( 0 < ncw ) then
        if ( sinth * xi + costh * yi <= yr ) then
          cycle
        end if
      end if

      dist = ( xi - x0 )**2 + ( yi - y0 )**2

      if ( dist < smdist ) then
        smdist = dist
        ind = i
      end if

    end do

    call rotiar ( nbc, bndcyc, ind )
    bndcyc(nbc) = bndcyc(0)
    nt = nbc + ncw
    tedg = cwalk + ncw + 1

  else

    call diam2 ( nvrt, xc(1), yc(1), i1, i2, dist, ierror )

    if ( ierror /= 0 ) then
      return
    end if

    ind = 0

    do

      if ( nbc <= ind ) then
        exit
      end if

      if ( xc(i1) == vcl(1,bndcyc(ind)) .and. &
           yc(i1) == vcl(2,bndcyc(ind)) ) then
        exit
      end if

      ind = ind + 1

    end do

    call rotiar ( nbc, bndcyc, ind )
    bndcyc(nbc) = bndcyc(0)
    mbc = 1

    do

      if ( nbc <= mbc ) then
        exit
      end if

      if ( xc(i2) == vcl(1,bndcyc(mbc)) .and. &
           yc(i2) == vcl(2,bndcyc(mbc)) ) then
        exit
      end if

      mbc = mbc + 1

    end do

    ind = nbc

    do i = mbc+1, mbc+(nbc-mbc-1)/2
      ind = ind - 1
      i1 = bndcyc(i)
      bndcyc(i) = bndcyc(ind)
      bndcyc(ind) = i1
    end do

    bndcyc(nbc) = bndcyc(mbc)
    nt = nbc - 2
    tedg = 1
!
!  Left boundary chain contains mesh vertices BNDCYC(0:MBC)
!  and right chain contains BNDCYC(0,MBC+1:NBC); MBC < NBC.
!
  end if

  if ( maxti < ntri + nt ) then
    ierror = 9
    return
  else if ( maxiw < tedg + 4*nt - 1 ) then
    ierror = 6
    return
  end if

  if ( inter ) then
    call tmerge ( inter, nbc, ncw, bndcyc, iwk(cwalk), ldv, vcl, &
      til(1,ntri+1), iwk(tedg), ierror )
  else
    call tmerge ( inter, mbc, nbc-mbc, bndcyc, bndcyc(mbc), ldv, vcl, &
      til(1,ntri+1), iwk(tedg), ierror )
  end if

  if ( ierror /= 0 ) then
    return
  end if

  sptr = tedg + 3 * nt

  call cvdtri ( inter, ldv, nt, vcl, til(1,ntri+1), iwk(tedg), iwk(sptr), &
    ierror )

  ntri = ntri + nt

  return
end
function umdf2 ( x, y )

!*****************************************************************************80
!
!! UMDF2 is a dummy mesh distribution function.
!
!  Purpose: 
!
!    Dummy user-supplied mesh distribution function which
!    is provided if heuristic mesh distribution function is used.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point.
!
!    Output, real ( kind = 8 ) UMDF2, the mesh distribution function value 
!    at (X,Y)
!
  real ( kind = 8 ) umdf2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  umdf2 = 1.0D+00

  return
end
function urand ( iy )

!*****************************************************************************80
!
!! URAND is a uniform random number generator.
!
!  Discussion:
!
!    URAND is a uniform random number generator based on theory and
!    suggestions given in D. E. Knuth (1969), Vol. 2. The integer IY
!    should be initialized to an arbitrary integer prior to the first
!    call to URAND. The calling program should not alter the value of
!    IY between subsequent calls to URAND. Values of URAND will be
!    returned in the interval (0,1).
!
!  Reference:
!
!    Forsythe, Malcolm, Moler, 
!    page 246.
!
!  Modified:
!
!    12 July 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IY, the seed value.
!
!    Output, real ( kind = 8 ) URAND, the random value.
!
  real ( kind = 8 ) halfm
  integer ( kind = 4 ), save :: ia = 0
  integer ( kind = 4 ), save :: ic = 0
  integer ( kind = 4 ), parameter :: itwo = 2
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: m2 = 0
  integer ( kind = 4 ), save :: mic = 0
  real ( kind = 8 ), save :: s = 0.0D+00
  real ( kind = 8 ) urand
!
!  If first entry, compute machine integer word length.
!
  if ( m2 == 0 ) then

    m = 1

    do

      m2 = m
      m = itwo * m2

      if ( m <= m2 ) then
        exit
      end if

    end do

    halfm = m2
!
!  Compute multiplier and increment for linear congruential method.
!
    ia = 8 * int ( halfm * atan ( 1.0D+00 ) / 8.0D+00 ) + 5
    ic = 2 * int ( halfm * ( 0.5D+00 - sqrt ( 3.0D+00 ) / 6.0D+00 ) ) + 1
    mic = ( m2 - ic ) + m2
!
!  S is the scale factor for converting to floating point.
!
    s = 0.5D+00 / halfm

  end if
!
!  Compute next random number.
!
  iy = iy * ia
!
!  The following statement is for computers which do not allow
!  integer overflow on addition.
!
  if ( mic < iy ) then
    iy = ( iy - m2 ) - m2
  end if

  iy = iy + ic
!
!  The following statement is for computers where the word
!  length for addition is greater than for multiplication.
!
  if ( m2 < iy / 2 ) then
    iy = (iy - m2) - m2
  end if
!
!  The following statement is for computers where integer
!  overflow affects the sign bit.
!
  if ( iy < 0 ) then
    iy = ( iy + m2 ) + m2
  end if

  urand = real ( iy, kind = 8 ) * s

  return
end
subroutine vbedg ( x, y, vcl, til, tnbr, ltri, ledg, rtri, redg )

!*****************************************************************************80
!
!! VBEDG determines visible boundary edges of a 2D triangulation.
!
!  Purpose: 
!
!    Determine boundary edges of 2D triangulation which are
!    visible from point (X,Y) outside convex hull.
!
!  Modified:
!
!    17 August 2003
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a 2D point outside
!    the convex hull.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list; negative 
!    values are used for links of counter clockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer ( kind = 4 ) LTRI, LEDG.  On input, if LTRI /= 0 then they 
!    are assumed to be as defined below and are not changed, else they are 
!    updated.  On output, LTRI is the index of the boundary triangle to the
!    left of leftmost boundary triangle visible from (X,Y), and LEDG is the
!    boundary edge of triangle LTRI to left of leftmost
!    boundary edge visible from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer ( kind = 4 ) RTRI, on input, the index of boundary triangle 
!    to begin search at.  On output, the index of rightmost boundary triangle 
!    visible from (X,Y).
!
!    Input/output, integer ( kind = 4 ) REDG.  On input, the edge of triangle RTRI that 
!    is visible from (X,Y).  On output, REDG has been updated so that this
!    is still true. 1 <= REDG <= 3.
!
  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  logical ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  real ( kind = 8 ) vcl(2,*)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Find rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -tnbr(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = til(e,t)

    if ( e <= 2 ) then
      b = til(e+1,t)
    else
      b = til(1,t)
    end if

    lr = lrline ( x, y, vcl(1,a), vcl(2,a), vcl(1,b), vcl(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e
 
  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = til(e,t)
    e = i4_wrap ( e-1, 1, 3 )

    do while ( 0 < tnbr(e,t) )

      t = tnbr(e,t)

      if ( til(1,t) == b ) then
        e = 3
      else if ( til(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = til(e,t)
    lr = lrline ( x, y, vcl(1,a), vcl(2,a), vcl(1,b), vcl(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end
subroutine vispol ( xeye, yeye, nvrt, xc, yc, nvis, ivis, ierror )

!*****************************************************************************80
!
!! VISPOL computes the visibility polygon.
!
!  Purpose: 
!
!    Compute the visibility polygon VP from an eyepoint in
!    the interior or blocked exterior of a simple polygon P or
!    on the boundary of a simply connected polygonal region P.
!    In the latter case, the interior angles at all vertices must
!    be strictly between 0 and 2*PI.
!
!  Discussion:
!
!    On input, XC and YC contain vertex coordinates of P. During
!    the algorithm, part of XC, YC is used as a stack, which, on
!    output, contains the vertex coordinates of VP. The stack
!    vertices overwrite the input vertices as the input vertices
!    are scanned. Elements of IVIS are set when vertices are added
!    to the stack; these values may have +NV or -NV added to them
!    to indicate that stack point has same angle as previous one.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe, R. B. Simpson, 
!    BIT 27 (1987), pages 458-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XEYE, YEYE, the coordinates of eyepoint; must 
!    be a simple vertex if it lies on the boundary (i.e. occurs only once).
!
!    Input, integer ( kind = 4 ) NVRT, the upper subscript of XC, YC (approximate
!    number of vertices).
!
!    Input/output, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT).  On input, if 
!    eyepoint is interior or blocked exterior then arrays contain coordinates 
!    in counter clockwise or clockwise order, respectively, with 
!    (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)); (XC(0),YC(0)) is a vertex visible from
!    (XEYE,YEYE), e.g. as computed by routine ROTIPG.  If eyepoint is a vertex
!    of P then arrays contain coordinates in counter clockwise order; 
!    (XC(0),YC(0)) is successor vertex of (XEYE,YEYE); (XC(NVRT),YC(NVRT)) is
!    predecessor vertex of (XEYE,YEYE).
!    On output, XC and YC contain the vertices of VP in counter clockwise order;
!    if eyepoint is interior or blocked exterior then
!    (XC(0),YC(0)) = (XC(NVIS),YC(NVIS)), else (XC(0),YC(0))
!    and (XC(NVIS),YC(NVIS)) are the successor and
!    predecessor vertices of (XEYE,YEYE) in VP.
!
!    Output, integer ( kind = 4 ) NVIS, the upper subscript of XC, YC on output (approximate
!    number of vertices of VP); NVIS <= NVRT.
!
!    Output, integer ( kind = 4 ) IVIS(0:NVIS), contains information about the vertices 
!    of VP with respect to the vertices of P; IVIS(I) = K if (XC(I),YC(I))
!    is the vertex of index K in the input polygon; IVIS(I) = -K if 
!    (XC(I),YC(I)) is on the interior of the edge joining vertices of index 
!    K-1 and K in input polygon
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 206, 207, 208, 209, or 210
!
  integer ( kind = 4 ) nvrt

  logical beye
  integer ( kind = 4 ) cur
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ivis(0:nvrt)
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) xe
  real ( kind = 8 ) xeye
  real ( kind = 8 ) xw
  real ( kind = 8 ) yc(0:nvrt)
  real ( kind = 8 ) ye
  real ( kind = 8 ) yeye
  real ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!     Variables in common block GVPVAR:
!        NV - NVRT
!        OPER - operation code 1 to 7 for LEFT, RIGHT, SCANA, SCANB,
!              SCANC, SCAND, FINISH
!        CUR - index of current vertex of P in XC, YC arrays
!        TOP - index of top vertex of stack in XC, YC arrays
!              (TOP <= CUR is always satisfied)
!        XE,YE - XEYE,YEYE
!        XW,YW - coordinates of point on last or second-last edge
!              processed (needed for routines VPSCNB, VPSCNC, VPSCND)
!        BEYE - .TRUE. iff eyepoint is on boundary
!
  ierror = 0
  beye = xc(0) /= xc(nvrt) .or. yc(0) /= yc(nvrt)
  nv = nvrt
  xe = xeye
  ye = yeye
  ivis(0) = 0
  cur = 1

  if ( beye ) then

    do

      lr = lrline ( xc(nv-1), yc(nv-1), xe, ye, xc(nv), yc(nv), 0.0D+00 )

      if ( lr /= 0 ) then
        exit
      end if
      nv = nv - 1

    end do

  end if

  do

    lr = lrline ( xc(cur), yc(cur), xe, ye, xc(0), yc(0), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    cur = cur + 1

  end do

  if ( lr == -1 ) then
    oper = 1
    if ( cur == 1 ) then
      top = 1
      ivis(1) = cur
    else if ( beye ) then
      top = 1
      xc(0) = xc(cur-1)
      yc(0) = yc(cur-1)
      ivis(0) = cur - 1
      xc(1) = xc(cur)
      yc(1) = yc(cur)
      ivis(1) = cur
    else
      top = 2
      xc(1) = xc(cur-1)
      yc(1) = yc(cur-1)
      ivis(1) = cur - 1 + nv
      xc(2) = xc(cur)
      yc(2) = yc(cur)
      ivis(2) = cur
    end if
  else
    oper = 3
    top = 0
    if ( beye .and. 1 < cur ) then
      xc(0) = xc(cur-1)
      yc(0) = yc(cur-1)
      ivis(0) = cur - 1
    end if
  end if
!
!  Angular displacement of stack points are in nondecreasing order,
!  with at most two consecutive points having the same displacement.
!
  do

    if ( oper == 1 ) then
      call vpleft ( xc, yc, ivis )
    else if ( oper == 2 ) then
      call vprght ( xc, yc, ivis, ierror )
    else if ( oper == 3 ) then
      call vpscna ( xc, yc, ivis, ierror )
    else if ( oper == 4 ) then
      call vpscnb ( xc, yc, ivis, ierror )
    else if ( oper == 5 ) then
      call vpscnc ( xc, yc, ivis, ierror )
    else
      call vpscnd ( xc, yc, ivis, ierror )
    end if

    if ( ierror /= 0 ) then
      nvis = top
      return
    end if

    if ( 6 < oper ) then
      exit
    end if

  end do
!
!  Add or subtract NV from those IVIS values which are used to
!  indicate that stack point has same angle as previous one.
!
  do i = 1, top

    if ( nv < ivis(i) ) then
      ivis(i) = ivis(i) - nv
    else if ( ivis(i) < -nv ) then
      ivis(i) = ivis(i) + nv
    end if

  end do

  nvis = top

  return
end
subroutine visvrt ( angspc, xeye, yeye, nvis, xc, yc, ivis, maxn, nvsvrt, &
  theta )

!*****************************************************************************80
!
!! VISVRT determines a list of visible vertices.
!
!  Purpose: 
!
!    Determine a list of visible vertices, ordered by
!    increasing "polar angle", on the boundary of the visibilty
!    polygon from boundary eyepoint (XEYE,YEYE).  This list
!    includes the vertices of visibility polygon such that a
!    line segment from (XEYE,YEYE) to vertex lies in interior
!    of polygon, as well as extra points on edges which subtend
!    an angle greater than or equal to 2*ANGSPC at (XEYE,YEYE).  
!
!    These extra points are at an equal angular spacing of at 
!    least ANGSPC and less than 2*ANGSPC.  The successor and predecessor 
!    of eyepoint are included in list.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in radians 
!    which controls how many extra points become visible vertices.
!
!    Input, real ( kind = 8 ) XEYE, YEYE, the coordinates of boundary eyepoint.
!
!    Input, integer ( kind = 4 ) NVIS, (number of vertices of visibility polygon) - 2.
!
!    Input/output, real ( kind = 8 ) XC(0:NVIS), YC(0:NVIS), on input, the 
!    coordinates of the vertices of visibility polygon in counter clockwise 
!    order; (XC(0),YC(0)) and (XC(NVIS),YC(NVIS)) are the successor and 
!    predecessor vertices of eyepoint in visibility polygon; at most 2
!    consecutive vertices have same polar angle with respect to eyepoint.
!    On output, coordinates of visible vertices which overwrite the 
!    input coordinates.
!
!    Input/output, IVIS(0:NVIS), on input, contains information about the 
!    vertices of XC, YC arrays with respect to the original polygon from
!    which visibility polygon is computed; if IVIS(I) is nonnegative
!    then (XC(I),YC(I)) has index I in original polygon;
!    if IVIS(I) < 0 then (XC(I),YC(I)) is on the edge
!    ending at vertex of index -IVIS(I) in original polygon;
!    indexing starts at 0 from successor of eyepoint.
!    On output, coordinates of visible vertices
!    which overwrite the input coordinates.
!
!    Input, integer ( kind = 4 ) MAXN, the upper bound on NVSVRT; should be at least
!    NVIS + INT(PHI/ANGSPC) where PHI is the interior angle at (XEYE,YEYE).
!
!    Output, integer ( kind = 4 ) NVSVRT, (number of visible vertices) - 1.
!
!    Output, real ( kind = 8 ) THETA(0:NVSVRT), the polar angles of visible 
!    vertices with respect to (XEYE,YEYE) at origin and (XC(0),YC(0))
!    on positive x-axis.
!
  integer ( kind = 4 ) maxn

  real ( kind = 8 ) alpha
  real ( kind = 8 ) ang
  real ( kind = 8 ) ang1
  real ( kind = 8 ) ang2
  real ( kind = 8 ) angdif
  real ( kind = 8 ) angle
  real ( kind = 8 ) angsp2
  real ( kind = 8 ) angspc
  real ( kind = 8 ) cosang
  integer ( kind = 4 ) cur
  real ( kind = 8 ) diff
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ivis(0:maxn)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) n
  real ( kind = 8 ) numer
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) nvsvrt
  real ( kind = 8 ) r
  real ( kind = 8 ) sinang
  real ( kind = 8 ) theta(0:maxn)
  real ( kind = 8 ) tol
  integer ( kind = 4 ) top
  real ( kind = 8 ) xc(0:maxn)
  real ( kind = 8 ) xeye
  real ( kind = 8 ) yc(0:maxn)
  real ( kind = 8 ) yeye

  tol = 100.0D+00 * epsilon ( tol )
!
!  Shift input vertices right, and possibly remove first and last
!  vertices due to collinearity with eyepoint.
!
  angsp2 = 2.0D+00 * angspc
  cur = maxn + 1
  n = maxn

  do i = nvis, 0, -1
    cur = cur - 1
    xc(cur) = xc(i)
    yc(cur) = yc(i)
    ivis(cur) = ivis(i)
  end do

  lr = lrline ( xc(cur+1), yc(cur+1), xeye, yeye, xc(cur), yc(cur), 0.0D+00 )

  if ( 0 <= lr ) then
    cur = cur + 1
    xc(0) = xc(cur)
    yc(0) = yc(cur)
    ivis(0) = ivis(cur)
  end if

  lr = lrline ( xc(n-1), yc(n-1), xeye, yeye, xc(n), yc(n), 0.0D+00 )

  if ( lr <= 0 ) then
    n = n - 1
  end if

  alpha = atan2 ( yc(0)-yeye, xc(0)-xeye )
  ang2 = 0.0D+00
  theta(0) = 0.0D+00
  top = 0
  cur = cur + 1
!
!  Process edge from vertices of indices CUR-1, CUR.
!
  do

    ang1 = ang2
    ang2 = angle ( xc(cur), yc(cur), xeye, yeye, xc(0), yc(0) )
    angdif = ang2 - ang1

    if ( angdif <= tol ) then
 
      diff = ( ( xc(cur) - xeye )**2 + ( yc(cur) - yeye)**2 ) - &
             ( ( xc(cur-1) - xeye )**2 + ( yc(cur-1) - yeye )**2 )

      if ( diff < 0.0D+00 ) then
        xc(top) = xc(cur)
        yc(top) = yc(cur)
        ivis(top) = ivis(cur)
        theta(top) = ang2
      end if

    else

      if ( angsp2 <= angdif ) then

        k = int ( angdif / angspc )
        ind = -abs ( ivis(cur))
        angdif = angdif / real ( k, kind = 8 )
        dx = xc(cur) - xc(cur-1)
        dy = yc(cur) - yc(cur-1)
        numer = ( xc(cur) - xeye ) * dy - ( yc(cur) - yeye ) * dx

        do i = 1, k-1
          top = top + 1
          theta(top) = ang1 + real ( i, kind = 8 ) * angdif
          ang = theta(top) + alpha
          cosang = cos(ang)
          sinang = sin(ang)
          r = numer / ( dy * cosang - dx * sinang )
          xc(top) = r * cosang + xeye
          yc(top) = r * sinang + yeye
          ivis(top) = ind
        end do
 
      end if

      top = top + 1
      xc(top) = xc(cur)
      yc(top) = yc(cur)
      ivis(top) = ivis(cur)
      theta(top) = ang2

    end if

    cur = cur + 1

    if ( n < cur ) then
      exit
    end if

  end do

  nvsvrt = top

  return
end
subroutine vornbr ( xeye, yeye, nvrt, xc, yc, nvor, ivor, xvor, yvor, ierror )

!*****************************************************************************80
!
!! VORNBR determines the Voronoi neighbors of an eyepoint.
!
!  Purpose: 
!
!    Determine the Voronoi neighbors of (XEYE,YEYE) from a list of vertices 
!    which are in increasing "polar angle" order.
!
!    The Voronoi neighbors are a sublist of this list.  The
!    Voronoi polygon is restricted to the sector formed from 
!    the edges joining (XEYE,YEYE) to the first and last vertices
!    of this list.  Each Voronoi neighbor corresponds to an edge
!    of the Voronoi polygon.
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XEYE, YEYE, the coordinates of the eyepoint.
!
!    Input, integer ( kind = 4 ) NVRT, (number of vertices in list) minus 1.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex 
!    coordinates from which Voronoi neighbors are determined; (XC(0),YC(0)),...,
!    (XC(NVRT),YC(NVRT)) are in increasing angular
!    displacement order with respect to (XEYE,YEYE).
!
!    Output, integer ( kind = 4 ) NVOR, (number of Voronoi neighbors) minus 1 [<= NVRT].
!
!    Output, integer ( kind = 4 ) IVOR(0:NVOR), the indices of Voronoi neighbors in XC, YC
!    arrays; 0 <= IVOR(0) < ... < IVOR(NVOR) <= NVRT.
!
!    Workspace, real ( kind = 8 ) XVOR(0:NVRT), YVOR(0:NVRT), arrays for
!    storing the vertex coordinates of the Voronoi polygon.
!
!    Output, integer ( kind = 4 ) IERROR, set to 212 if an error occurred.
!
  integer ( kind = 4 ) nvrt

  real ( kind = 8 ) a11
  real ( kind = 8 ) a12
  real ( kind = 8 ) a21
  real ( kind = 8 ) a22
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) det
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) im
  integer ( kind = 4 ) ivor(0:nvrt)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nvor
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs
  real ( kind = 8 ) xc(0:nvrt)
  real ( kind = 8 ) xeye
  real ( kind = 8 ) xi
  real ( kind = 8 ) xvor(0:nvrt)
  real ( kind = 8 ) yc(0:nvrt)
  real ( kind = 8 ) yeye
  real ( kind = 8 ) yi
  real ( kind = 8 ) yvor(0:nvrt)

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  k = 1
  m = 0
  ivor(0) = 0
  xvor(0) = ( xeye + xc(0) ) * 0.5D+00
  yvor(0) = ( yeye + yc(0) ) * 0.5D+00
!
!  Beginning of main loop
!
  do while ( k <= nvrt )
!
!  Determine the intersection of the perpendicular bisectors of edges 
!  from (XEYE,YEYE) to (XC(K),YC(K)) and 
!  from (XEYE,YEYE) to (XC(IM),YC(IM)).
!
     im = ivor(m)

     a11 = xc(k) - xeye
     a12 = yc(k) - yeye
     a21 = xc(im) - xeye
     a22 = yc(im) - yeye

     tolabs = tol * max ( abs ( a11 ), abs ( a12 ), abs ( a21 ), abs ( a22 ) )
     det = a11 * a22 - a21 * a12

     if ( abs ( det ) <= tolabs ) then
       ierror = 212
       return
     end if

     b1 = ( a11**2 + a12**2 ) * 0.5D+00
     b2 = ( a21**2 + a22**2 ) * 0.5D+00

     xi = ( b1 * a22 - b2 * a12 ) / det
     yi = ( b2 * a11 - b1 * a21 ) / det
!
!  Determine whether (XVOR(M+1),YVOR(M+1)) is to the left of or
!  on the directed line from (XEYE,YEYE) to (XVOR(M),YVOR(M)).
!
     xvor(m+1) = xi + xeye
     yvor(m+1) = yi + yeye
     lr = lrline ( xvor(m+1), yvor(m+1), xeye, yeye, xvor(m), yvor(m), 0.0D+00 )

     if ( lr <= 0 ) then
       m = m + 1
       ivor(m) = k
       k = k + 1
     else if ( 0 < m ) then
       m = m - 1
     else
!
!  Determine the intersection of edge from (XEYE,YEYE) to
!  (XC(0),YC(0)) and the perpendicular bisector of the edge
!  from (XEYE,YEYE) to (XC(K),YC(K)).
!
      a11 = xc(k) - xeye
      a12 = yc(k) - yeye
      a21 = yc(0) - yeye
      a22 = xeye - xc(0)
      tolabs = tol * max ( abs ( a11 ), abs ( a12 ), abs ( a21 ), abs ( a22 ) )
      det = a11 * a22 - a21 * a12

      if ( abs ( det) <= tolabs ) then
        ierror = 212
        return
      end if

      b1 = ( a11**2 + a12**2 ) * 0.5D+00
      b2 = 0.0D+00
      xi = ( b1 * a22 - b2 * a12 ) / det
      yi = ( b2 * a11 - b1 * a21 ) / det
      xvor(m) = xi + xeye
      yvor(m) = yi + yeye
      ivor(m) = k
      k = k + 1

    end if

  end do
!
!  The following short loop determines which vertices at the end
!  of list are not Voronoi neighbors.
!
  do

    lr = lrline ( xvor(m), yvor(m), xeye, yeye, xc(nvrt), yc(nvrt), 0.0D+00 )

    if ( 0 <= lr ) then
      exit
    end if

    m = m - 1
    if ( m < 0 ) then
      exit
    end if

  end do

  nvor = m

  return
end
subroutine vpleft ( xc, yc, ivis )

!*****************************************************************************80
!
!! VPLEFT is called by routine VISPOL for the LEFT operation (OPER = 1).
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!     Input and updated parameters:
!
!        XC,YC,IVIS - see comments in routine VISPOL
!
  logical beye
  integer ( kind = 4 ) cur
  logical intsct
  integer ( kind = 4 ) ivis(0:*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lr1
  integer ( kind = 4 ) lr2
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real ( kind = 8 ) xc(0:*)
  real ( kind = 8 ) xe
  real ( kind = 8 ) xu
  real ( kind = 8 ) xw
  real ( kind = 8 ) yc(0:*)
  real ( kind = 8 ) ye
  real ( kind = 8 ) yu
  real ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a left turn, S(TOP) = V(CUR), TOP <= CUR,
!  S(TOP-1) = V(CUR-1) or on interior of edge V(CUR-1)-V(CUR).
!
10 continue

  if ( cur == nv ) then
    oper = 7
    return
  end if

  if ( .not. beye .and. top <= 2 ) then
    go to 20
  end if
!
!  Check if angular displacement of stack chain is greater than or equal
!  to 2*PI or interior angle at boundary viewpoint.
!
  call xedge ( 1, xe, ye, xc(nv), yc(nv), xc(top-1), yc(top-1), xc(top), &
    yc(top), xu, yu, intsct )

  if ( intsct ) then

    oper = 4
    xw = xc(cur)
    yw = yc(cur)
    lr = lrline ( xc(top), yc(top), xe, ye, xc(nv), yc(nv), 0.0D+00 )

    if ( lr == -1 ) then
      xc(top) = xu
      yc(top) = yu
      ivis(top) = -cur
    end if

    return

  end if
!
!  Process next edge.
!
20 continue

  lr = lrline ( xc(cur+1), yc(cur+1), xe, ye, xc(cur), yc(cur), 0.0D+00 )

  if ( lr == -1 ) then

    cur = cur + 1
    top = top + 1
    xc(top) = xc(cur)
    yc(top) = yc(cur)
    ivis(top) = cur

  else

    j = cur + 1
    lr1 = lrline ( xc(j), yc(j), xc(top-1), yc(top-1), xc(cur), yc(cur), &
      0.0D+00 )

    if ( lr1 == 1 ) then

      oper = 3
      cur = j

    else

      if ( lr == 1 ) then
        lr2 = 1
        go to 40
      end if

      do

        j = j + 1
        lr2 = lrline ( xc(j), yc(j), xe, ye, xc(cur), yc(cur), 0.0D+00 )

        if ( lr2 /= 0 ) then
          exit
        end if

      end do

40    continue

      if ( lr2 == -1 ) then
        top = top + 1
        xc(top) = xc(j-1)
        yc(top) = yc(j-1)
        ivis(top) = j - 1 + nv
        top = top + 1
        xc(top) = xc(j)
        yc(top) = yc(j)
        ivis(top) = j
      else
        oper = 2
      end if

      cur = j

    end if

  end if
!
!  This test avoids extra subroutine calls.
!
  if ( oper == 1 ) then
    go to 10
  end if

  return
end
subroutine vprght ( xc, yc, ivis, ierror )

!*****************************************************************************80
!
!! VPRGHT is called by routine VISPOL for the RIGHT operation (OPER = 2).
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!     Input and updated parameters:
!        XC,YC,IVIS - see comments in routine VISPOL
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 206.
!
  logical beye
  integer ( kind = 4 ) case
  integer ( kind = 4 ) cur
  integer ( kind = 4 ) ierror
  logical intsct
  integer ( kind = 4 ) ivis(0:*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lr1
  integer ( kind = 4 ) lr2
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real ( kind = 8 ) xc(0:*)
  real ( kind = 8 ) xe
  real ( kind = 8 ) xu
  real ( kind = 8 ) xw
  real ( kind = 8 ) yc(0:*)
  real ( kind = 8 ) ye
  real ( kind = 8 ) yu
  real ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a right turn, EYE-S(TOP)-V(CUR) is a right
!  turn, EYE-S(TOP-1)-S(TOP) is a left turn, TOP < CUR, S(TOP) =
!  V(CUR-1) and S(TOP-1)-S(TOP)-V(CUR) is a left turn or S(TOP) is
!  not on edge V(CUR-1)-V(CUR) and V(CUR-1)-V(CUR) intersects
!  EYE-S(TOP).
!  Pop points from stack. If BEYE, it is not possible for
!  (XC(CUR),YC(CUR)) to be identical to any stack points.
!
  ierror = 0

10 continue

  case = 0
  j = top

20 continue

  if ( abs ( ivis(j)) <= nv ) then

    lr = lrline ( xc(cur), yc(cur), xe, ye, xc(j-1), yc(j-1), 0.0D+00 )

    if ( lr == -1 ) then

      case = 1

    else if ( lr == 0 ) then

      if ( abs ( ivis(j-1)) <= nv ) then
        j = j - 1
        case = 2
      else if ( (xc(j-1) - xe)**2 + (yc(j-1) - ye)**2 <= &
                (xc(j-2) - xe)**2 + (yc(j-2) - ye)**2 ) then
        j = j - 2
        case = 2
      else
        case = -1
      end if

    end if

  else if ( case == -1 ) then

    if ( ( xc(cur) - xe )**2 + ( yc(cur) - ye )**2 <= &
         ( xc(j-1) - xe )**2 + ( yc(j-1) - ye )**2 ) then
      j = j - 1
      case = 2
    else
      xw = xc(cur)
      yw = yc(cur)
      case = 3
    end if

  else

    call xedge ( 0, xc(cur-1), yc(cur-1), xc(cur), yc(cur), &
      xc(j-1), yc(j-1), xc(j), yc(j), xw, yw, intsct )

    if ( intsct ) then
      case = 3
    end if

  end if

  if ( 0 < case ) then
    go to 30
  end if

  j = j - 1
  if ( 1 <= j ) then
    go to 20
  end if
!
!  Error from no more edges in stack.
!
  ierror = 206
  return
!
!  Process next edge.
!
30 continue

  if ( case == 3 ) then

    oper = 6
    top = j - 1

  else

    top = j
    xw = xc(cur-1)
    yw = yc(cur-1)

    if ( case == 1 ) then
      call xedge ( 1, xe, ye, xc(cur), yc(cur), xc(top-1), yc(top-1), &
            xc(top), yc(top), xu, yu, intsct )
      xc(top) = xu
      yc(top) = yu
      ivis(top) = -abs ( ivis(top))
    end if

    lr = lrline ( xc(cur+1), yc(cur+1), xe, ye, xc(cur), yc(cur), 0.0D+00 )

    if ( lr == 1 ) then

      cur = cur + 1

    else

      j = cur + 1
      lr1 = lrline ( xc(j), yc(j), xw, yw, xc(cur), yc(cur), 0.0D+00 )

      if ( lr1 == -1 ) then

        oper = 5
        cur = j

      else

        if ( lr == -1 ) then
          lr2 = -1
          go to 50
        end if

        do

          j = j + 1
          lr2 = lrline ( xc(j), yc(j), xe, ye, xc(cur), yc(cur), 0.0D+00 )

          if ( lr2 /= 0 ) then
            exit
          end if

        end do

50      continue

        if ( lr2 == -1 ) then
          oper = 1
          top = top + 1
          xc(top) = xc(j-1)
          yc(top) = yc(j-1)
          ivis(top) = j - 1 + nv
          top = top + 1
          xc(top) = xc(j)
          yc(top) = yc(j)
          ivis(top) = j
        end if

        cur = j

      end if

    end if

  end if
!
!  This test avoids extra subroutine calls.
!
  if ( oper == 2 ) then
    go to 10
  end if

  return
end
subroutine vpscna ( xc, yc, ivis, ierror )

!*****************************************************************************80
!
!! VPSCNA is called by routine VISPOL for the SCANA operation (OPER = 3).
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input and updated parameters:
!    XC,YC,IVIS - see comments in routine VISPOL
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 207.
!
  logical beye
  integer ( kind = 4 ) case
  integer ( kind = 4 ) cur
  integer ( kind = 4 ) ierror
  logical intsct
  integer ( kind = 4 ) ivis(0:*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lr1
  integer ( kind = 4 ) lr2
  integer ( kind = 4 ) lr3
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real ( kind = 8 ) xc(0:*)
  real ( kind = 8 ) xe
  real ( kind = 8 ) xw
  real ( kind = 8 ) yc(0:*)
  real ( kind = 8 ) ye
  real ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a right turn or forward move, S(TOP) =
!  V(CUR-1) or EYE-S(TOP)-V(CUR-1) is a forward move and TOP = 0,
!  TOP < CUR; S(TOP-1)-S(TOP)-V(CUR) is a right turn if TOP is at least 1
!  or EYE-S(TOP)-V(CUR) is a right turn if TOP = 0.
!  If BEYE, it is possible that (XC(TOP),YC(TOP)) is a non-simple
!  vertex but any edge incident on this vertex encountered during
!  scan must be invisible from (XE,YE).
!
  ierror = 0
  k = cur

10 continue

  if ( xc(k+1) == xc(top) .and. yc(k+1)  ==  yc(top) ) then

    k = k + 2

  else

    call xedge ( 1, xe, ye, xc(top), yc(top), xc(k), yc(k), xc(k+1), &
      yc(k+1), xw, yw, intsct )

    if ( intsct ) then

      lr = lrline ( xc(k+1), yc(k+1), xe, ye, xc(k), yc(k), 0.0D+00 )

      if ( lr == 1 ) then

        if ( ( xw      - xe )**2 + ( yw      - ye )**2 <= &
             ( xc(top) - xe )**2 + ( yc(top) - ye )**2 ) then

          if ( 0 < top ) then
            case = 1
            go to 20
          end if

        else

          lr1 = lrline ( xc(k), yc(k), xe, ye, xc(top), yc(top), 0.0D+00 )

          if ( lr1 == -1 ) then
            case = 2
            go to 20
          end if

        end if

      else

        lr1 = lrline ( xc(k+1), yc(k+1), xe, ye, xc(top), yc(top), 0.0D+00 )

        if ( lr1 == -1 ) then
          case = 3
          go to 20
        end if

      end if

    end if

    k = k + 1

  end if

  if ( k < nv ) then
    go to 10
  end if
!
!  Error from unsuccessful scan.
!
  ierror = 207
  return
!
!  Process current edge.
!
20 continue

  if ( case == 3 ) then

    oper = 1
    cur = k + 1
    lr = lrline ( xc(k), yc(k), xe, ye, xc(top), yc(top), 0.0D+00 )
    top = top + 1

    if ( lr == 0 ) then
      xc(top) = xc(k)
      yc(top) = yc(k)
      ivis(top) = k + nv
    else
      xc(top) = xw
      yc(top) = yw
      ivis(top) = -(k + 1 + nv)
    end if

    top = top + 1
    xc(top) = xc(cur)
    yc(top) = yc(cur)
    ivis(top) = cur

  else if ( case == 1 ) then

    cur = k + 1
    lr = lrline ( xc(cur), yc(cur), xe, ye, xc(top), yc(top), 0.0D+00 )

    if ( lr == 1 ) then

      oper = 2

    else

      j = cur + 1
      lr1 = lrline ( xc(j), yc(j), xe, ye, xc(cur), yc(cur), 0.0D+00 )
      lr2 = lrline ( xc(j), yc(j), xc(k), yc(k), xc(cur), yc(cur), 0.0D+00 )

      if ( lr1 <= 0 .and. lr2 == -1 ) then

        oper = 5
        xw = xc(k)
        yw = yc(k)
        cur = j

      else

        if ( lr1 /= 0 ) then
          lr3 = lr1
          go to 40
        end if

        do

          j = j + 1
          lr3 = lrline ( xc(j), yc(j), xe, ye, xc(cur), yc(cur), 0.0D+00 )
 
          if ( lr3 /= 0 ) then
            exit
          end if

        end do

40      continue

        if ( lr3 == 1 ) then
          oper = 2
        else
          oper = 1
          top = top + 1
          xc(top) = xc(j-1)
          yc(top) = yc(j-1)
          ivis(top) = j - 1 + nv
          top = top + 1
          xc(top) = xc(j)
          yc(top) = yc(j)
          ivis(top) = j
        end if

        cur = j

      end if

    end if

  else

    oper = 6
    cur = k + 1
    lr = lrline ( xc(cur), yc(cur), xe, ye, xc(top), yc(top), 0.0D+00 )

    if ( lr == 0 ) then
      xw = xc(cur)
      yw = yc(cur)
    end if

  end if

  return
end
subroutine vpscnb ( xc, yc, ivis, ierror )

!*****************************************************************************80
!
!! VPSCNB is called by routine VISPOL for the SCANB operation (OPER = 4).
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!     Input and updated parameters:
!        XC,YC,IVIS - see comments in routine VISPOL
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 208.
!
  logical beye
  integer ( kind = 4 ) cur
  integer ( kind = 4 ) ierror
  logical intsct
  integer ( kind = 4 ) ivis(0:*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lr1
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs
  integer ( kind = 4 ) top
  real ( kind = 8 ) xc(0:*)
  real ( kind = 8 ) xe
  real ( kind = 8 ) xu
  real ( kind = 8 ) xw
  real ( kind = 8 ) yc(0:*)
  real ( kind = 8 ) ye
  real ( kind = 8 ) yu
  real ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  EYE-V(CUR-1)-V(CUR) is a left turn, S(TOP) = V(CUR) or S(TOP) is
!  on interior of edge V(CUR-1)-V(CUR), TOP <= CUR, S(TOP) has
!  angular displacement of 2*PI or interior angle at boundary eye.
!  (XW,YW) is the input version of (XC(CUR),YC(CUR)).
!  If BEYE, it is possible that (XC(TOP),YC(TOP)) is a non-simple
!  point but any edge containing this point encountered during scan
!  must be invisible from (XE,YE), except for 1 case where K = CUR.
!
  tolabs = tol * ( ( xc(nv) - xc(top) )**2 + ( yc(nv) - yc(top))**2 )
  k = cur

  if ( ivis(top) < 0 .or. k + 1 == nv ) then
    go to 10
  end if

  lr = lrline ( xc(k+1), yc(k+1), xe, ye, xc(top), yc(top), 0.0D+00 )

  lr1 = lrline ( xc(k+1), yc(k+1), xc(top-1), yc(top-1), xc(top), yc(top), &
    0.0D+00 )

  if ( lr == 1 .and. lr1  ==  -1 ) then
    oper = 2
    cur = k + 1
    return
  else
    k = k + 1
  end if

10 continue

  if ( k + 1 == nv ) then

    oper = 7
    cur = nv
    top = top + 1
    xc(top) = xc(nv)
    yc(top) = yc(nv)
    ivis(top) = nv
    return

  else

      if ( k == cur ) then
         call xedge ( 0, xc(nv), yc(nv), xc(top), yc(top), xw, yw, &
           xc(k+1), yc(k+1), xu, yu, intsct )
      else
         call xedge ( 0, xc(nv), yc(nv), xc(top), yc(top), xc(k), yc(k), &
           xc(k+1), yc(k+1), xu, yu, intsct )
      end if

      if ( intsct ) then
         if ( ( xc(top) - xu )**2 + ( yc(top) - yu )**2 <= tolabs ) then
           go to 20
         end if
         lr = lrline ( xc(k+1), yc(k+1), xe, ye, xc(nv), yc(nv), 0.0D+00 )
         if ( lr == 1 ) then
           oper = 2
           cur = k + 1
           return
         end if
      end if

20    continue

      k = k + 1

   end if

  if ( k < nv ) then
    go to 10
  end if
!
!  Error from unsuccessful scan.
!
  ierror = 208

  return
end
subroutine vpscnc ( xc, yc, ivis, ierror )

!*****************************************************************************80
!
!! VPSCNC is called by routine VISPOL for the SCANC operation (OPER = 5).
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input/output, XC, YC, IVIS - see comments in routine VISPOL
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 209.
!
  logical beye
  integer ( kind = 4 ) cur
  integer ( kind = 4 ) ierror
  logical intsct
  integer ( kind = 4 ) ivis(0:*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lr1
  integer ( kind = 4 ) lr2
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real ( kind = 8 ) xc(0:*)
  real ( kind = 8 ) xe
  real ( kind = 8 ) xp
  real ( kind = 8 ) xu
  real ( kind = 8 ) xw
  real ( kind = 8 ) yc(0:*)
  real ( kind = 8 ) ye
  real ( kind = 8 ) yp
  real ( kind = 8 ) yu
  real ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a left turn or forward move, EYE-V(CUR-2)-
!  V(CUR-1) is a right turn, V(CUR-2)-V(CUR-1)-V(CUR) is a left turn,
!  TOP < CUR-1, W = V(CUR-2), S(TOP) is not on V(CUR-1)-V(CUR), EYE-
!  S(TOP)-V(CUR-1) is a backward move, EYE-S(TOP-1)-S(TOP) is a left
!  turn. If BEYE, it is possible that V(CUR-1) is a non-simple point,
!  but intersection at (XC(TOP),YC(TOP)) cannot occur.
!
  ierror = 0
  xp = xc(cur-1)
  yp = yc(cur-1)
  k = cur

10 continue

  if ( xc(k+1) == xp .and. yc(k+1) == yp ) then

    go to 40

  else if ( xc(k) == xp .and. yc(k) == yp ) then

      j = k + 1
      lr = lrline ( xc(j), yc(j), xe, ye, xp, yp, 0.0D+00 )
      lr1 = lrline ( xc(j), yc(j), xw, yw, xp, yp, 0.0D+00 )

      if ( lr <= 0 .and. lr1 == -1 ) then
        go to 40
      end if

      if ( lr /= 0 ) then

        lr2 = lr

      else

        do

          j = j + 1
          lr2 = lrline ( xc(j), yc(j), xe, ye, xp, yp, 0.0D+00 )
 
          if ( lr2 /= 0 ) then
            exit
          end if

        end do

      end if

      if ( lr2 == 1 ) then
        oper = 2
      else
        oper = 1
        top = top + 1
        xc(top) = xc(j-1)
        yc(top) = yc(j-1)
        ivis(top) = j - 1 + nv
        top = top + 1
        xc(top) = xc(j)
        yc(top) = yc(j)
        ivis(top) = j
      end if

      cur = j
      return

   else

      call xedge ( 0, xp, yp, xc(top), yc(top), xc(k), yc(k), xc(k+1), &
        yc(k+1), xu, yu, intsct )

      if ( intsct ) then
        lr = lrline ( xc(k+1), yc(k+1), xe, ye, xp, yp, 0.0D+00 )
        if ( lr == 1 ) then
          oper = 2
          cur = k + 1
          return
        end if
      end if

   end if

40    continue

   k = k + 1

  if ( k < nv ) then
    go to 10
  end if
!
!  Error from unsuccessful scan.
!
  ierror = 209

  return
end
subroutine vpscnd ( xc, yc, ivis, ierror )

!*****************************************************************************80
!
!! VPSCND is called by routine VISPOL for the SCAND operation (OPER = 6).
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input/output, XC,YC,IVIS - see comments in routine VISPOL
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 210
!
  logical beye
  integer ( kind = 4 ) cur
  integer ( kind = 4 ) ierror
  logical intsct
  integer ( kind = 4 ) ivis(0:*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lr1
  integer ( kind = 4 ) lr2
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real ( kind = 8 ) xc(0:*)
  real ( kind = 8 ) xe
  real ( kind = 8 ) xp
  real ( kind = 8 ) xu
  real ( kind = 8 ) xw
  real ( kind = 8 ) yc(0:*)
  real ( kind = 8 ) ye
  real ( kind = 8 ) yp
  real ( kind = 8 ) yu
  real ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a right turn, S(TOP) is a V vertex not on
!  V(CUR-1)-V(CUR), TOP < CUR, W is intersection of V(CUR-1)-V(CUR)
!  and ray EYE-S(TOP), EYE-S(TOP)-W is a forward move, and
!  EYE-S(TOP-1)-S(TOP) is a left turn if TOP is at least 1.
!  If BEYE, it is possible that (XW,YW) is a non-simple point,
!  but intersection at (XC(TOP),YC(TOP)) cannot occur.
!
  ierror = 0
  xp = xc(cur-1)
  yp = yc(cur-1)
  k = cur

10 continue

  call xedge ( 0, xw, yw, xc(top), yc(top), xc(k), yc(k), xc(k+1), yc(k+1), &
    xu, yu, intsct )

  if ( intsct ) then

    lr = lrline ( xc(k+1), yc(k+1), xe, ye, xc(k), yc(k), 0.0D+00 )
    lr1 = lrline ( xc(k+1), yc(k+1), xe, ye, xc(top), yc(top), 0.0D+00 )

    if ( lr == -1 .and. lr1  ==  -1 ) then

      if ( xc(k) /= xw .or. yc(k) /= yw ) then
        go to 20
      end if

      lr2 = lrline ( xc(k+1), yc(k+1), xp, yp, xw, yw, 0.0D+00 )

      if ( lr2 == -1 ) then
        go to 30
      end if

20    continue

         oper = 1
         cur = k + 1
         lr2 = lrline ( xc(k), yc(k), xe, ye, xc(top), yc(top), 0.0D+00 )
         top = top + 1
         if ( lr2 == 0 ) then
          xc(top) = xc(k)
          yc(top) = yc(k)
          ivis(top) = k + nv
         else
          xc(top) = xu
          yc(top) = yu
          ivis(top) = -(k + 1 + nv)
         end if
         top = top + 1
         xc(top) = xc(cur)
         yc(top) = yc(cur)
         ivis(top) = cur
         return
      end if

  end if

30 continue

  k = k + 1

  if ( k < nv ) then
    go to 10
  end if
!
!  Error from unsuccessful scan.
!
  ierror = 210

  return
end
subroutine walkt2 ( x, y, ntri, vcl, til, tnbr, itri, iedg, ierror )

!*****************************************************************************80
!
!! WALKT2 searches for a triangle containing a point.
!
!  Purpose: 
!
!    Walk through neighboring triangles of a 2D Delaunay
!    triangulation until a triangle is found containing point (X,Y)
!    or (X,Y) is found to be outside the convex hull.  Search is
!    guaranteed to terminate for a Delaunay triangulation, else a
!    cycle may occur.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a 2D point.
!
!    Input, integer ( kind = 4 ) NTRI, the number of triangles in the triangulation; used 
!    to detect cycle.
!
!    Input, real ( kind = 8 ) VCL(2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) TIL(3,NTRI), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TNBR(3,NTRI), the triangle neighbor list.
!
!    Input/output, integer ( kind = 4 ) ITRI.  On input, the index of triangle to begin 
!    search at.  On output, the index of triangle that search ends at.
!
!    Output, integer ( kind = 4 ) IEDG, indicates the position of the point (X,Y) in
!    triangle ITRI.  A small tolerance is allowed in positions:
!    0, the interior of the triangle; 
!    1, interior of edge 1;
!    2, interior of edge 2;
!    3, interior or edge 3;
!    4, vertex 1;
!    5, vertex 2;
!    6, vertex 3;
!    -1, outside convex hull, past edge 1;
!    -2, outside convex hull, past edge 2;
!    -3, outside convex hull, past edge 3.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  On abnormal return,
!    IERROR is set to 226.
!
  integer ( kind = 4 ) ntri

  integer ( kind = 4 ) a
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cnt
  real ( kind = 8 ) det
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxa
  real ( kind = 8 ) dxb
  real ( kind = 8 ) dy
  real ( kind = 8 ) dya
  real ( kind = 8 ) dyb
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedg
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itri
  integer ( kind = 4 ) til(3,ntri)
  integer ( kind = 4 ) tnbr(3,ntri)
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(2,*)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  cnt = 0
  iedg = 0
  ierror = 0

  do

    cnt = cnt + 1

    if ( ntri < cnt ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WALKT2 - Fatal error!'
      write ( *, '(a)' ) '  All triangles have been searched.'
      ierror = 226
      return
    end if
!
!  Get the vertices of triangle ITRI.
!
    a = til(1,itri)
    b = til(2,itri)
    c = til(3,itri)
!
!  Using vertex C as a base, compute the distances to vertices A and B,
!  and the point (X,Y).
!
    dxa = vcl(1,a) - vcl(1,c)
    dya = vcl(2,a) - vcl(2,c)

    dxb = vcl(1,b) - vcl(1,c)
    dyb = vcl(2,b) - vcl(2,c)

    dx = x - vcl(1,c)
    dy = y - vcl(2,c)

    det = dxa * dyb - dya * dxb
!
!  Compute the barycentric coordinates of the point (X,Y) with respect
!  to this triangle.
!
    alpha = ( dx * dyb - dy * dxb ) / det
    beta = ( dxa * dy - dya * dx ) / det
    gamma = 1.0D+00 - alpha - beta
!
!  If the barycentric coordinates are all positive, then the point
!  is inside the triangle.
!
    if ( tol < alpha .and. tol < beta .and. tol < gamma ) then
      exit
    end if
!
!  If any barycentric coordinate is (strongly) negative with respect to
!  a side, and if that side is on the convex hull, the point is outside
!  the triangles, and we are done.
!
    if ( alpha < -tol ) then
      i = tnbr(2,itri)
      if ( i <= 0 ) then
        iedg = -2
        exit
      end if
    else if ( beta < -tol ) then
      i = tnbr(3,itri)
      if ( i <= 0 ) then
        iedg = -3
        exit
      end if
    else if ( gamma < -tol ) then
      i = tnbr(1,itri)
      if ( i <= 0 ) then
        iedg = -1
        exit
      end if
!
!  At least one barycentric coordinate is between -TOL and TOL,
!  and no barycentric coordinate is less than -TOL.  We are going
!  to assign the position to an edge or vertex.
!
    else if ( alpha <= tol ) then
      if ( beta <= tol ) then
        iedg = 6
      else if ( gamma <= tol ) then
        iedg = 5
      else
        iedg = 2
      end if
      exit
    else if ( beta <= tol ) then
      if ( gamma <= tol ) then
        iedg = 4
      else
        iedg = 3
      end if
      exit
    else
      iedg = 1
      exit
    end if
!
!  If we fell through, then at least one barycentric coordinate was negative
!  for a side of the current triangle, and that side has a neighboring
!  triangle I.  Let's go there.
!
    itri = i

  end do

  return
end
subroutine width2 ( nvrt, xc, yc, i1, i2, widsq, ierror )

!*****************************************************************************80
!
!! WIDTH2 finds the minimum breadth of a convex polygon.
!
!  Discussion:
!
!    WIDTH2 finds the width (minimum breadth) of a convex polygon with
!    vertices given in counter-clockwise order and with all interior 
!    angles < PI.
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices.
!
!    Input, real ( kind = 8 ) XC(1:NVRT), YC(1:NVRT), the vertex coordinates, in
!    counter-clockwise order.
!
!    Output, integer ( kind = 4 ) I1, I2, indices in XC, YC such that the width is 
!    the distance from vertex (XC(I1),YC(I1)) to the line joining 
!    (XC(I2),YC(I2)) and (XC(I2+1),YC(I2+1)), where index NVRT+1 
!    is same as 1.
!
!    Output, real ( kind = 8 ) WIDSQ, the square of the width of the polygon.
!
!    Output, integer ( kind = 4 ) IERROR, the error flag.
!    0, no error was detected.
!    201, an error was detected.
!
  integer ( kind = 4 ) nvrt

  integer ( kind = 4 ) a
  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) areatr
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) c1mtol
  real ( kind = 8 ) c1ptol
  real ( kind = 8 ) dist
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) m
  real ( kind = 8 ) tol
  real ( kind = 8 ) widsq
  real ( kind = 8 ) xc(nvrt)
  real ( kind = 8 ) yc(nvrt)

  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Find the first vertex which is farthest from the edge connecting
!  vertices NVRT and 1.
!
  c1mtol = 1.0D+00 - tol
  c1ptol = 1.0D+00 + tol
  j = nvrt
  jp1 = 1
  k = 2
  area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

  do

    area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k+1), yc(k+1) )

    if ( area2 <= area1 * c1ptol ) then
      exit
    end if

    area1 = area2
    k = k + 1

  end do

  m = k
  widsq = 0.0D+00
!
!  Find width = minimum distance of antipodal edge-vertex pairs.
!
  area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

  do

    kp1 = i4_wrap ( k+1, 1, nvrt )

    area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(kp1), yc(kp1) )

    if ( area1 * c1ptol < area2 ) then

      a = j
      b = k
      k = k + 1
      c = k
      if ( nvrt < c ) then
        c = 1
      end if
      area1 = area2

    else if ( area2 < area1 * c1mtol ) then

      a = k
      b = j
      c = jp1
      j = jp1
      jp1 = j + 1
      area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

    else

      a = k
      b = j
      c = jp1
      k = k + 1
      j = jp1
      jp1 = j + 1
      area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

    end if

    if ( m < j .or. nvrt < k ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WIDTH2 - Fatal error!'
      write ( *, '(a)' ) '  M < J or NVRT < K'
      write ( *, '(a,i6)' ) '  J = ', j
      write ( *, '(a,i6)' ) '  M = ', m
      write ( *, '(a,i6)' ) '  K = ', k
      write ( *, '(a,i6)' ) '  NVRT = ', nvrt
      ierror = 201
      return
    end if

    dx = xc(c) - xc(b)
    dy = yc(c) - yc(b)
    dist = ( ( yc(a) - yc(b) ) * dx - ( xc(a) - xc(b) ) * dy )**2 &
      / ( dx**2 + dy**2 )

    if ( dist < widsq .or. widsq <= 0.0D+00 ) then
      widsq = dist
      i1 = a
      i2 = b
    end if

    if ( j == m .and. k == nvrt ) then
      exit
    end if

  end do

  return
end
subroutine xedge ( mode, xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, xu, yu, &
  intsct )

!*****************************************************************************80
!
!! XEDGE determines if an edge intersects another edge or ray.
!
!  Discussion:
!
!    An edge is a finite line segment.  A ray is a semi-infinite line
!    segment.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MODE, is 0 for two edges, 1 (or nonzero) for a ray 
!    and an edge.
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, XW1, YW1, XW2, YW2, the
!    vertex coordinates;  an edge (ray) is from (XV1,YV1) to (thru) (XV2,YV2);
!    an edge joins vertices (XW1,YW1) and (XW2,YW2).
!
!    Output, real ( kind = 8 ) XU, YU, the coordinates of the point of 
!    intersection iff INTSCT is .TRUE.
!
!    Output, logical INTSCT, .TRUE. if the edges/ray are nondegenerate, not
!    parallel, and intersect, .FALSE. otherwise.
!
  real ( kind = 8 ) denom
  real ( kind = 8 ) dxv
  real ( kind = 8 ) dxw
  real ( kind = 8 ) dyv
  real ( kind = 8 ) dyw
  logical intsct
  integer ( kind = 4 ) mode
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) xw1
  real ( kind = 8 ) xw2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2
  real ( kind = 8 ) yw1
  real ( kind = 8 ) yw2

  tol = 100.0D+00 * epsilon ( tol )

  intsct = .false.
  dxv = xv2 - xv1
  dyv = yv2 - yv1
  dxw = xw2 - xw1
  dyw = yw2 - yw1
  tolabs = tol * max ( abs ( dxv ), abs ( dyv ), abs ( dxw ), abs ( dyw ) )
  denom = dyv * dxw - dxv * dyw

  if ( abs ( denom ) <= tolabs ) then
    return
  end if

  t = ( dyv * ( xv1 - xw1 ) - dxv * ( yv1 - yw1 ) ) / denom

  if ( t < -tol .or. 1.0D+00 + tol < t ) then
    return
  end if

  xu = xw1 + t * dxw
  yu = yw1 + t * dyw

  if ( abs ( dyv ) <= abs ( dxv ) ) then
    t = ( xu - xv1 ) / dxv
  else
    t = ( yu - yv1 ) / dyv
  end if

  if ( mode == 0 ) then
    if ( -tol <= t .and. t <= 1.0D+00 + tol ) then
      intsct = .true.
    end if
  else
    if ( -tol <= t ) then
      intsct = .true.
    end if
  end if

  return
end
subroutine xline ( xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, dv, dw, xu, &
  yu, parall )

!*****************************************************************************80
!
!! XLINE finds the intersection of lines parallel to two other lines.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, XW1, YW1, XW2, YW2, the 
!    vertex coordinates; the first line is parallel to and at signed distance 
!    DV to the left of directed line from (XV1,YV1) to (XV2,YV2);
!    second line is parallel to and at signed distance DW to
!    left of directed line from (XW1,YW1) to (XW2,YW2)
!
!    Input, real ( kind = 8 ) DV, DW, the signed distances (positive for left).
!
!    Output, real ( kind = 8 ) XU, YU, the coordinates of the point of
!    intersection, if PARALL is .FALSE.
!
!    Output, logical PARALL, is .TRUE. if the lines are parallel, or two 
!    points for a line are identical, .FALSE. otherwise.
!
  real ( kind = 8 ) a11
  real ( kind = 8 ) a12
  real ( kind = 8 ) a21
  real ( kind = 8 ) a22
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) det
  real ( kind = 8 ) dv
  real ( kind = 8 ) dw
  logical parall
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) xw1
  real ( kind = 8 ) xw2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2
  real ( kind = 8 ) yw1
  real ( kind = 8 ) yw2

  tol = 100.0D+00 * epsilon ( tol )

  parall = .true.

  a11 = yv2 - yv1
  a12 = xv1 - xv2
  a21 = yw2 - yw1
  a22 = xw1 - xw2
  tolabs = tol * max ( abs ( a11), abs ( a12), abs ( a21), abs ( a22) )
  det = a11 * a22 - a21 * a12

  if ( abs ( det ) <= tolabs ) then
    return
  end if

  b1 = xv1 * a11 + yv1 * a12 - dv * sqrt ( a11**2 + a12**2 )
  b2 = xw1 * a21 + yw1 * a22 - dw * sqrt ( a21**2 + a22**2 )

  xu = ( b1 * a22 - b2 * a12 ) / det
  yu = ( b2 * a11 - b1 * a21 ) / det

  parall = .false.

  return
end
