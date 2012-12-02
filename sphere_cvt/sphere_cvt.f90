program main

!*****************************************************************************80
!
!! MAIN is the main program for SPHERE_CVT.
!
!  Discussion:
!
!    SPHERE_CVT runs the spherical centroidal Voronoi tessellation code.
!
!  Diary:
!
!    01 May 2010: Program renamed from SCVT to SPHERE_CVT.
!
!    27 October 2004: DOUBLE PRECISION replaced by REAL ( KIND = 8 ).
!
!    19 June 2002: Adding Halton points as an initialization and sampling.
!    Something goes wrong and I get NaN's after a while.
!
!    11 June 2002: I added DISECTION_TYPE and DISECTION_NUM, to 
!    disect the soccerball center mesh using Delaunay or Voronoi
!    refinement.  I only tried DISECTION_NUM = 1, and I only
!    plotted the initial points, I didn't try to run a test case.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Qiang Du, Vance Faber, Max Gunzburger,
!    Centroidal Voronoi Tesselations: Applications and Algorithms,
!    SIAM Review, Volume 41, 1999, pages 637-676.
!
!    Douglas Hardin, Edward Saff,
!    Discretizing Manifolds via Minimum Energy Points,
!    Notices of the American Mathematical Society,
!
!    Edward Saff, Arno Kuijlaars,
!    Distributing Many Points on a Sphere,
!    The Mathematical Intelligencer,
!    Volume 19, Number 1, 1997, pages 5-11.
!
!  Parameters:
!
!    integer ( kind = 4 ) N, the number of Voronoi cells to generate.
!    A typical value is 256.
!
!    integer ( kind = 4 ) IT_MAX, the maximum number of correction iterations 
!    used in the Voronoi calculation.  A typical value is 10.
!
!    integer ( kind = 4 ) SAMPLE_NUM, the total number of sampling points 
!    tested.  A typical value is 5000 * N.
!
!    integer ( kind = 4 ) SAMPLE_BATCH, the maximum number of sample points to
!    generate at one time.  For problems where N is large, and so SAMPLE_NUM
!    is large, setting SAMPLE_BATCH to 1,000,000 or less avoids
!    memory problems.
!
!    Input, integer ( kind = 4 ) RANDOM_GENERATOR, indicates how the initial 
!    Voronoi cell generators are chosen.
!    0, random points.
!    1, spiral points.
!    2, soccer centers (N must be 32).
!    3, soccer vertices (N must be 60).
!    4, Halton points
!
!    integer ( kind = 4 ) SAMPLE_TYPE, how sampling is done.
!    1, random points.
!    2, Halton points.
!
!    integer ( kind = 4 ) DISECTION_TYPE, if the soccerball centers are used, 
!    and if 0 < DISECTION_NUM, then this variable specifies how the disection
!    is to be done.
!    0, add a point on the bisector between contiguous centers, adding
!       one center for every edge of the original grid;
!    1, add every Voronoi vertex of the original grid.
!
!    integer ( kind = 4 ) DISECTION_NUM, the number of disections to apply if 
!    the soccerball centers are used as initial generators.
!
!    real ( kind = 8 ) GENERATOR_INIT(DIM_NUM,N), the initial Voronoi cell
!    generators.
!
!    real ( kind = 8 ) GENERATOR(DIM_NUM,N), the Voronoi cell generators 
!    of the Voronoi tessellation, as approximated by the SPHERE_CVT algorithm.  
!    This is the output quantity of most interest.
!
!    integer ( kind = 4 ) SEED, determines how to initialize the RANDOM_NUMBER 
!    routine.  If SEED is zero on input, then RANDOM_INITIALIZE will make up a 
!    seed from the current double precision time clock reading.
!    If SEED is nonzero, then a reproducible sequence of random numbers
!    defined by SEED will be chosen.
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: b
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: c
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: centroid
  logical debug
  integer ( kind = 4 ) disection_num
  integer ( kind = 4 ), parameter :: disection_type = 0
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: generator
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: generator_init
  character ( len = 255 ) :: file_name = 'generators.xyz'
  integer ( kind = 4 ), parameter :: file_unit = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 2000
  integer ( kind = 4 ) :: n = 32
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) next
  integer ( kind = 4 ), parameter :: random_generator = 0
  integer ( kind = 4 ), parameter :: sample_batch = 1000000
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ), parameter :: sample_type = 1
  integer ( kind = 4 ) :: seed = 1

  debug = .false.
  disection_num = 0
!
!  Print introduction and options.
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_CVT'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A sample problem for the probabilistic'
  write ( *, '(a)' ) '  Spherical Centroidal Voronoi Tessellation algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Given a unit sphere in 3D, the problem is to determine'
  write ( *, '(a)' ) '  a set of GENERATORS, that is, a set of points which '
  write ( *, '(a)' ) '  lie on the sphere, and which implicitly define a'
  write ( *, '(a)' ) '  division of the surface into Voronoi cells.'
  write ( *, '(a)' ) '  It is also desired that each generator point actually'
  write ( *, '(a)' ) '  be the centroid of its cell.'
!
!  Initialize the random number generator.
!
  call random_initialize ( seed )
!
!  Initialize the Voronoi cell generators.
!
  allocate ( generator_init(1:3,1:n) )

  call generator_initialize ( random_generator, n, generator_init )
!
!  Carry out disection if requested.
!
!  I would REALLY like to do the disection in a subroutine.
!  But I don't know allocatable arrays well enough to be sure I can do this!
!
  if ( random_generator == 2 ) then

    if ( 0 < disection_num ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Disection option = ', disection_type
      write ( *, '(a,i8)' ) '  Number of disections = ', disection_num
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Initial N = ', n
      write ( *, '(a)' ) ' '

      n3 = n
      allocate ( c(1:3,1:n3) )
      c(1:3,1:n3) = generator_init(1:3,1:n)
      deallocate ( generator_init )

      do i = 1, disection_num

        n1 = n3

        allocate ( a(1:3,1:n1) )
        a(1:3,1:n1) = c(1:3,1:n1)

        deallocate ( c )

       if ( disection_type == 0 ) then

          n2 = 3 * n1 - 6
          allocate ( b(1:3,1:n2) )

          call delaunay_midpoints ( n1, a, n2, b )

        else if ( disection_type == 1 ) then

          n2 = 2 * n1 - 4
          allocate ( b(1:3,1:n2) )
          call voronoi_vertices ( n1, a, n2, b )

        end if
 
        n3 = n1 + n2
        allocate ( c(1:3,1:n3) )
        c(1:3,1:n1) = a(1:3,1:n1)
        c(1:3,n1+1:n1+n2) = b(1:3,1:n2)

 
        deallocate ( a )
        deallocate ( b )
 
        write ( *, '(a,i8,a,i8)' ) '  Disection ', i, ' results in N = ', n3

      end do

      n = n3
      allocate ( generator_init(1:3,1:n) )
      generator_init(1:3,1:n) = c(1:3,1:n3)
      deallocate ( c )

    end if

  end if
!
!  Now we are sure we know what N is, so we can allocate other things.
!
  allocate ( centroid(1:3,1:n) )
  allocate ( generator(1:3,1:n) )
  sample_num = 1000 * n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_CVT Algorithm parameters:'
  write ( *, '(a)' ) '-------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of Voronoi cells to generate: ', n
  write ( *, '(a,i8)' ) '  Number of iterations to determine CVT: ', it_max
  write ( *, '(a,i12)' ) '  Number of sampling points:           ', sample_num
  write ( *, '(a,i12)' ) '  Sampling is done in batches of size  ', sample_batch
  if ( sample_type == 1 ) then
    write ( *, '(a)' ) '  Sample is done by RANDOM points.'
  else if ( sample_type == 2 ) then
    write ( *, '(a)' ) '  Sample is done by HALTON points.'
  end if
!
!  Write initial generators to a file.
!
  call r8mat_write ( 'initial.xyz', 3, n, generator )

  generator(1:3,1:n) = generator_init(1:3,1:n)
!
!  Carry out the SPHERE_CVT iteration, which drives the Voronoi generators 
!  and Voronoi centroids closer and closer.
!
  next = 1

  do it = 1, it_max

    call sphere_cvt_centroid ( n, generator, sample_num, sample_batch, &
      sample_type, centroid )

    if ( it == next ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 'STEP = ', it
      write ( *, '(a)' ) '  Discrepancy between generator and centroid'

      call motion ( n, generator, centroid )

      if ( debug ) then
        next = next + 1
      else if ( it < 5 ) then
        next = next + 1
      else if ( it_max < 11 ) then
        next = next + 1
      else
        next = ( ( ( 10 * it ) / it_max ) + 1 ) * ( it_max / 10 )
      end if
      
      next = min ( next, it_max )

    end if

    generator(1:3,1:n) = centroid(1:3,1:n)

  end do
!
!  Write generators to files.
!
  call r8mat_write ( file_name, 3, n, generator )
!
!  Determine motion of generators from initial to final.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Total motion of generators from initial to final.'
  write ( *, '(a)' ) ' '

  call motion ( n, generator, generator_init )
!
!  Free memory.
!
  deallocate ( centroid )
  deallocate ( generator )
  deallocate ( generator_init )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_CVT_MAIN'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine delaunay_midpoints ( nc, centers, nm, midpoints )

!*****************************************************************************80
!
!! DELAUNAY_MIDPOINTS returns the midpoints of a Delaunay triangulation.
!
!  Discussion:
!
!    A set of NC points on the unit sphere is given.
!
!    The appropriate STRIPACK routines are called to generate the
!    Voronoi diagram of the points.
!
!    For each pair of centers that are joined by a line of a Delaunay
!    trianglation, we want to compute the midpoint.
!
!    If we start with NC center points, then the Delaunay triangulation
!    will have:
!
!    * NC vertices
!    * the number of edges will be E = ( 3 * F ) / 2 because each
!      triangle adds 3 sides, but each side is added twice.
!    * Euler's formula will read F - ( 3 * F ) / 2 + NC = 2
!
!    Therefore, the number of triangles will be:
!
!      F = 2 * ( NC - 2 )
!
!    and the number of edges will be 
!
!      E = 3 * ( NC - 2 )
!
!    hence the number of midpoints will be
!
!      NM = E = 3 * ( NC - 2 ).
!
!    If we add the NM midpoints to our NC = NC(0) center points to create a
!    refined mesh of NC(1) points, then:
!
!      NC(1) = NC(0) + 3 * ( NC(0) - 2 ) = 4 * ( NC(0) - 2 ) + 2
!
!    which can be rewritten as:
!
!      NC(1) - 2 = 4 * ( NC(0) - 2 )
!
!    and if this process is repeated K times,
!
!      NC(K) - 2 = 4**K * ( NC(0) - 2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NC, the number of "center" points.
!
!    Input, real ( kind = 8 ) CENTERS(3,NC), the coordinates of the center
!    points.
!
!    Output, integer ( kind = 4 ) NM, the number of midpoints.
!
!    Output, real ( kind = 8 ) MIDPOINTS(3,NM), the coordinates 
!    of the midpoints.
!
  implicit none

  integer ( kind = 4 ) nc

  real ( kind = 8 ) centers(3,nc)
  real ( kind = 8 ) ds(nc)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) iwk(2*nc)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lend(nc)
  integer ( kind = 4 ) list(6*nc)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lptr(6*nc)
  integer ( kind = 4 ) ltri(6,2*nc-4)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nm
  integer ( kind = 4 ) nt
  real ( kind = 8 ) midpoints(3,3*nc-6)
  real ( kind = 8 ) x(nc)
  real ( kind = 8 ) y(nc)
  real ( kind = 8 ) z(nc)
!
!  Copy the points out.
!
  x(1:nc) = centers(1,1:nc)
  y(1:nc) = centers(2,1:nc)
  z(1:nc) = centers(3,1:nc)
!
!  Create the triangulation.
!
  call trmesh ( nc, x, y, z, list, lptr, lend, lnew, iwk, iwk(nc+1), ds, &
    ierror )

  if ( ierror == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DELAUNAY_MIDPOINTS - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first 3 nodes are collinear.'
    stop
  end if

  if ( 0 < ierror ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DELAUNAY_MIDPOINTS - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Create a triangle list.
!
  call trlist ( nc, list, lptr, lend, 6, nt, ltri, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DELAUNAY_MIDPOINTS - Fatal error!'
    write ( *, '(a)' ) '  Error in TRLIST.'
    stop
  end if
!
!  We expect that every pair of nodes (N1,N2) connected by a Delaunay
!  side will be listed also as (N2,N1).  To avoid repetition, we select
!  the listing in which N1 < N2.
!
  nm = 0

  do j = 1, nt

    do i = 1, 3

      n1 = ltri(i,j)
      ip1 = i4_wrap ( i+1, 1, 3 )
      n2 = ltri(ip1,j)

      if ( n1 < n2 ) then
        nm = nm + 1
        midpoints(1:3,nm) = 0.5D+00 * ( centers(1:3,n1) + centers(1:3,n2) )
      end if

    end do

  end do
!
!  Normalize the points.
!
  do j = 1, nm
    midpoints(1:3,j) = midpoints(1:3,j) / sqrt ( sum ( midpoints(1:3,j)**2 ) )
  end do

  return
end
subroutine find_closest ( dim_num, n, sample_num, s, r, nearest )

!*****************************************************************************80
!
!! FIND_CLOSEST finds the nearest R point to each S point.
!
!  Discussion:
!
!    This routine finds the closest Voronoi cell generator by checking every
!    one.  For problems with many cells, this process can take the bulk
!    of the CPU time.  Other approaches, which group the cell generators into
!    bins, can run faster by a large factor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of cell generators.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sample points.
!
!    Input, real ( kind = 8 ) S(DIM_NUM,SAMPLE_NUM), the points to be checked.
!
!    Input, real ( kind = 8 ) R(DIM_NUM,N), the cell generators.
!
!    Output, integer ( kind = 4 ) NEAREST(SAMPLE_NUM), the index of the nearest
!    cell generators.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) sample_num

  real ( kind = 8 ) dist_sq_min
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) jr
  integer ( kind = 4 ) js
  integer ( kind = 4 ) nearest(sample_num)
  real ( kind = 8 ) r(dim_num,n)
  real ( kind = 8 ) s(dim_num,sample_num)

  do js = 1, sample_num

    dist_sq_min = huge ( dist_sq_min )
    nearest(js) = -1

    do jr = 1, n

      dist_sq = sum ( ( r(1:dim_num,jr) - s(1:dim_num,js) )**2 )

      if ( dist_sq < dist_sq_min ) then
        dist_sq_min = dist_sq
        nearest(js) = jr
      end if

    end do

  end do

  return
end
subroutine generator_initialize ( random_generator, n, generator )

!*****************************************************************************80
!
!! GENERATOR_INITIALIZE sets initial values for the generators.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANDOM_GENERATOR, indicates the initialization.
!    0, random points.
!    1, spiral points.
!    2, soccer centers (N must be 32).
!    3, soccer vertices (N must be 60).
!    4, Halton points.
!
!    Input, integer ( kind = 4 ) N, the number of cell generatorrs.
!
!    Output, real ( kind = 8 ) GENERATOR(3,N), the initial values for 
!    the cell generators.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) generator(3,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) random_generator

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initialize the Voronoi cell generators by...'
  if ( random_generator == 0 ) then
    write ( *, '(a)' ) '  ...RANDOM POINTS.'
  else if ( random_generator == 1 ) then
    write ( *, '(a)' ) '  ...SPIRAL POINTS'
  else if ( random_generator == 2 ) then
    write ( *, '(a)' ) '  ...SOCCER BALL CENTERS (Requires N = 32)'
  else if ( random_generator == 3 ) then
    write ( *, '(a)' ) '  ...SOCCER BALL VERTICES (Requires N = 60)'
  else if ( random_generator == 4 ) then
    write ( *, '(a)' ) '  ...HALTON points'
  end if

  if ( random_generator == 0 ) then

    call sphere_unit_samples_3d ( n, generator )

  else if ( random_generator == 1 ) then

    call sphere_unit_spiralpoints_3d ( n, generator )

  else if ( random_generator == 2 ) then

    if ( n == 32 ) then

      call soccer_centers ( generator )

      do j = 1, n
        generator(1:3,j) = generator(1:3,j) &
          / sqrt ( sum ( generator(1:3,j)**2 ) )
      end do

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPHERE_CVT - Fatal error!'
      write ( *, '(a)' ) '  This option requires N = 32.'
      stop

    end if

  else if ( random_generator == 3 ) then

    if ( n == 60 ) then
      call soccer_vertices ( generator )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPHERE_CVT - Fatal error!'
      write ( *, '(a)' ) '  This option requires N = 60.'
      stop
    end if

  else if ( random_generator == 4 ) then

    call sphere_unit_haltons_3d ( n, generator )

  end if

  call r8mat_transpose_print ( 3, n, generator, '  Initial generators:' )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2008
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
subroutine halton_memory ( action, name, dim_num, value )

!*****************************************************************************80
!
!! HALTON_MEMORY sets or returns quantities associated with the Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action.
!    'GET' means get the value of a particular quantity.
!    'SET' means set the value of a particular quantity.
!    'INC' means increment the value of a particular quantity.
!          (Only the SEED can be incremented.)
!
!    Input, character ( len = * ) NAME, the name of the quantity.
!    'BASE' means the Halton base or bases.
!    'DIM_NUM' means the spatial dimension.
!    'SEED' means the current Halton seed.
!
!    Input/output, integer ( kind = 4 ) DIM_NUM, the dimension of the quantity.
!    If ACTION is 'SET' and NAME is 'BASE', then DIM_NUM is input, and
!    is the number of entries in VALUE to be put into BASE.
!
!    Input/output, integer ( kind = 4 ) VALUE(DIM_NUM), contains a value.
!    If ACTION is 'SET', then on input, VALUE contains values to be assigned
!    to the internal variable.
!    If ACTION is 'GET', then on output, VALUE contains the values of
!    the specified internal variable.
!    If ACTION is 'INC', then on input, VALUE contains the increment to
!    be added to the specified internal variable.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ), allocatable, save :: base(:)
  logical, save :: first_call = .true.
  integer ( kind = 4 ) i
  character ( len = * ) name
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), save :: dim_num_save = 0
  integer ( kind = 4 ) prime
  integer ( kind = 4 ), save :: seed = 1
  integer ( kind = 4 ) value(*)

  if ( first_call ) then
    dim_num_save = 1
    allocate ( base(dim_num_save) )
    base(1) = 2
    first_call = .false.
  end if
!
!  Set
!
  if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name(1:1) == 'B' .or. name(1:1) == 'b' ) then

      if ( dim_num_save /= dim_num ) then
        deallocate ( base )
        dim_num_save = dim_num
        allocate ( base(dim_num_save) )
      end if

      base(1:dim_num) = value(1:dim_num)

    else if ( name(1:1) == 'N' .or. name(1:1) == 'n' ) then

      if ( dim_num_save /= value(1) ) then
        deallocate ( base )
        dim_num_save = value(1)
        allocate ( base(dim_num_save) )
        do i = 1, dim_num_save
          base(i) = prime ( i )
        end do
      else
        dim_num_save = value(1)
      end if

    else if ( name(1:1) == 'S' .or. name(1:1) == 's' ) then

      seed = value(1)

    end if
!
!  Get
!
  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name(1:1) == 'B' .or. name(1:1) == 'b' ) then
      value(1:dim_num_save) = base(1:dim_num_save)
    else if ( name(1:1) == 'N' .or. name(1:1) == 'n' ) then
      value(1) = dim_num_save
    else if ( name(1:1) == 'S' .or. name(1:1) == 's' ) then
      value(1) = seed
    end if
!
!  Increment
!
  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name(1:1) == 'S' .or. name(1:1) == 's' ) then
      seed = seed + value(1)
    end if

  end if

  return
end
subroutine halton_vector_sequence ( dim_num, n, r )

!*****************************************************************************80
!
!! HALTON_VECTOR_SEQUENCE: next N elements in the vector Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the element.
!
!    Input, integer ( kind = 4 ) N, the number of elements desired.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the next N elements of the 
!    current vector Halton sequence.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base(dim_num)
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value(1)

  call halton_memory ( 'GET', 'SEED', 1, value )
  seed = value(1)

  call halton_memory ( 'GET', 'BASE', dim_num, base )

  call i4_to_halton_vector_sequence ( seed, base, dim_num, n, r )

  value(1) = n
  call halton_memory ( 'INC', 'SEED', 1, value )

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_to_halton_vector_sequence ( seed, base, dim_num, n, r )

!*****************************************************************************80
!
!! I4_TO_HALTON_VECTOR_SEQUENCE computes N elements of a vector Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    Numerische Mathematik,
!    Volume 2, pages 84-90.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, the index of the desired element.
!    Only the absolute value of SEED is considered.  SEED = 0 is allowed,
!    and returns R = 0.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the Halton bases, which should 
!    be distinct prime numbers.  This routine only checks that each base
!    is greater than 1.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the sequence.
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the SEED-th through (SEED+N-1)-th
!    elements of the Halton sequence for the given bases.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  real ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed2(n)
!
!  We assume that N is large compared to DIM_NUM, so our implicit inner
!  loop is on N.
!
  do i = 1, dim_num

    call i4vec_indicator ( n, seed2 )

    seed2(1:n) = seed2(1:n) - 1 + abs ( seed )

    r(i,1:n) = 0.0D+00

    if ( base(i) <= 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I_TO_HALTON_VECTOR_SEQUENCE - Fatal error!'
      write ( *, '(a)' ) '  An input base BASE(I) is <= 1!'
      write ( *, '(a,i8)' ) '  For index I = ', i
      write ( *, '(a,i8)' ) '  BASE(I) = ', base(i)
      stop
    end if

    base_inv = 1.0D+00 / real ( base(i), kind = 8 )

    do while ( any ( seed2(1:n) /= 0 ) )
      digit(1:n) = mod ( seed2(1:n), base(i) )
      r(i,1:n) = r(i,1:n) + real ( digit(1:n), kind = 8 ) * base_inv
      base_inv = base_inv / real ( base(i), kind = 8 )
      seed2(1:n) = seed2(1:n) / base(i)
    end do

  end do

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

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
    i4_wrap = ilo + i4_modp ( ival - ilo, wide )
  end if

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector A(I)=I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine motion ( n, x1, x2 )

!*****************************************************************************80
!
!! MOTION computes the "motion" between two sets of points on the sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X1(3,N), X2(3,N), the sets of points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ) arc_cosine
  integer ( kind = 4 ) i
  real ( kind = 8 ) motion_average
  real ( kind = 8 ) motion_max
  real ( kind = 8 ) motion_min
  real ( kind = 8 ) motion_total
  real ( kind = 8 ) x1(3,n)
  real ( kind = 8 ) x2(3,n)

  motion_max = - huge ( motion_max )
  motion_min =   huge ( motion_min )
  motion_total = 0.0D+00

  do i = 1, n

    angle = arc_cosine ( dot_product ( x1(1:3,i), x2(1:3,i) ) )
!
!  This call to ACOS results in occasional NaN's...
!
!   angle = acos ( dot_product ( x1(1:3,i), x2(1:3,i) ) )

    motion_max = max ( motion_max, angle )
    motion_min = min ( motion_min, angle )
    motion_total = motion_total + angle

  end do

  motion_average = motion_total / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Generator motions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Total =   ', motion_total
  write ( *, '(a,g14.6)' ) '  Average = ', motion_average
  write ( *, '(a,g14.6)' ) '  Minimum = ', motion_min
  write ( *, '(a,g14.6)' ) '  Maximum = ', motion_max

  return
end
function prime ( n )

!*****************************************************************************80
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Discussion:
!
!    PRIME_MAX is 1600, and the largest prime stored is 13499.
!
!    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, pages 870-873.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired prime number.
!    In general, is should be true that 0 <= N <= PRIME_MAX.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!
!    Output, integer ( kind = 4 ) PRIME, the N-th prime.  If N is out of range,
!    PRIME is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: prime_max = 1600

  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( prime_max ) :: npvec
  integer ( kind = 4 ) prime

  if ( icall == 0 ) then

    icall = 1

    npvec(1:100) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

    npvec(101:200) = (/ &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

    npvec(201:300) = (/ &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

    npvec(301:400) = (/ &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

    npvec(401:500) = (/ &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

    npvec(501:600) = (/ &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

    npvec(601:700) = (/ &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

    npvec(701:800) = (/ &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

    npvec(801:900) = (/ &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

    npvec(901:1000) = (/ &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

    npvec(1001:1100) = (/ &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

    npvec(1101:1200) = (/ &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

    npvec(1201:1300) = (/ &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

    npvec(1301:1400) = (/ &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

    npvec(1401:1500) = (/ &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

   npvec(1501:1600) = (/ &
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, &
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, &
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, &
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, &
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, &
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, &
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, &
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, &
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, &
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /)

  end if

  if ( n == -1 ) then
    prime = prime_max
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= prime_max ) then
    prime = npvec(n)
  else
    prime = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
  end if

  return
end
subroutine r83vec_unit_l2 ( n, x )

!*****************************************************************************80
!
!! R83VEC_UNIT_L2 makes each R83 vector in an R83VEC have unit L2 norm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vectors.
!
!    Input/output, real ( kind = 8 ) X(3,N), the coordinates of N vectors.
!    On output, the nonzero vectors have been scaled to have unit L2 norm.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) norm
  real ( kind = 8 ) x(3,n)

  do i = 1, n

    norm = sqrt ( sum ( x(1:3,i)**2 ) )

    if ( norm /= 0.0D+00 ) then
      x(1:3,i) = x(1:3,i) / norm
    end if

  end do

  return
end
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
  character ( len = * )  title

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
!    10 September 2009
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

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 8 ) t
!
!  Initialize the random number seed.
!
  call random_seed ( )
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )

  if ( seed /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) '  Initialize RANDOM_NUMBER with user SEED = ', seed

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) &
      '  Initialize RANDOM_NUMBER with arbitrary SEED = ', seed

  end if
!
!  Now set the seed.
!
  seed_vector(1:seed_size) = seed

  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do

  return
end
subroutine sphere_cvt_centroid ( n, generator, sample_num, sample_batch, &
  sample_type, centroid )

!*****************************************************************************80
!
!! SPHERE_CVT_CENTROID computes the centroids of the regions.
!
!  Discussion:
!
!    The routine is given a set of points, called "generators", which
!    define a tessellation of the region into Voronoi cells.  Each point
!    defines a cell.  Each cell, in turn, has a centroid, but it is
!    unlikely that the centroid and the generator coincide.
!
!    Each time this CVT iteration is carried out, an attempt is made
!    to modify the generators in such a way that they are closer and
!    closer to being the centroids of the Voronoi cells they generate.
!
!    A large number of sample points are generated, and the nearest generator
!    is determined.  A count is kept of how many points were nearest to each
!    generator.  Once the sampling is completed, the location of all the
!    generators is adjusted.  This step should decrease the discrepancy
!    between the generators and the centroids.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi cells.
!
!    Input, real ( kind = 8 ) GENERATOR(3,N), the Voronoi cell generators. 
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the total number of sample points.
!
!    Input, integer ( kind = 4 ) SAMPLE_BATCH, the maximum size of a single
!    batch of sample points.
!
!    Output, real ( kind = 8 ) CENTROID(3,N), the Voronoi cell centroids.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) sample_batch

  real ( kind = 8 ) centroid(3,n)
  real ( kind = 8 ) generator(3,n)
  integer ( kind = 4 ) count(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest(sample_batch)
  real ( kind = 8 ) norm
  integer ( kind = 4 ), parameter :: random_sampler = 0
  real ( kind = 8 ) sample(3,sample_batch)
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) sample_size
  integer ( kind = 4 ) sample_sofar
  integer ( kind = 4 ) sample_type
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value(1)

  centroid(1:3,1:n) = 0.0D+00
  count(1:n) = 0

  sample_sofar = 0

  if ( sample_type == 2 ) then
    seed = 1
    value(1) = seed
    call halton_memory ( 'SET', 'SEED', 1, value )
  end if

  do 
    
    sample_size = sample_batch
    sample_size = min ( sample_size, sample_num - sample_sofar )

    if ( sample_type == 1 ) then
      call sphere_unit_samples_3d ( sample_size, sample )
    else if ( sample_type == 2 ) then
      call sphere_unit_haltons_3d ( sample_size, sample )
    end if

    call find_closest ( 3, n, sample_size, sample, generator, nearest )

    do j = 1, sample_size

      centroid(1:3,nearest(j)) = centroid(1:3,nearest(j)) + sample(1:3,j)

      count(nearest(j)) = count(nearest(j)) + 1

    end do

    sample_sofar = sample_sofar + sample_size

    if ( sample_num <= sample_sofar ) then
      exit
    end if

  end do
!
!  Average.
!
  do j = 1, n
    if ( count(j) /= 0 ) then
      centroid(1:3,j) = centroid(1:3,j) / real ( count(j), kind = 8 )
    end if
  end do
!
!  Normalize.
!
  do j = 1, n
    if ( count(j) /= 0 ) then
      norm = sqrt ( sum ( centroid(1:3,j)**2 ) )
      if ( norm /= 0.0D+00 ) then
        centroid(1:3,j) = centroid(1:3,j) / norm
      end if
    end if
  end do

  return
end
subroutine soccer_centers ( x )

!*****************************************************************************80
!
!! SOCCER_CENTERS returns the centers of the truncated icosahedron in 3D.
!
!  Discussion:
!
!    The shape is a truncated icosahedron, which is the design used
!    on a soccer ball.  There are 12 pentagons and 20 hexagons.
!
!    The centers are computed by averaging the vertices.  Note that
!    the vertices of the shape lie on the unit sphere, but the face
!    centers do not.  To force this to happen, simply normalize each point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X(3,32), the coordinates of the 32 face 
!    centers of a truncated icosahedron.
!
  implicit none

  real ( kind = 8 ) x(3,32)

  x(1:3,1:32) = reshape ( &
  (/  0.176084,  0.333333,  0.854729, &
     -0.374587, -0.042441,  0.854729, &
      0.185551, -0.271913,  0.900685, &
      0.681916, -0.042441,  0.637077, &
     -0.371103,  0.543825,  0.697234, &
     -0.209088, -0.650456,  0.637077, &
      0.616466,  0.543825,  0.493784, &
      0.443867, -0.650456,  0.502561, &
     -0.804597,  0.222222,  0.419426, &
      0.086406,  0.830237,  0.419426, &
     -0.731143, -0.375774,  0.493784, &
      0.904859,  0.222222,  0.067258, &
      0.866776, -0.375774,  0.164595, &
     -0.519688,  0.761567,  0.150394, &
      0.033908, -0.944118,  0.164595, &
     -0.536815, -0.761567,  0.067258, &
      0.536815,  0.761567, -0.067258, &
      0.519688, -0.761567, -0.150394, &
     -0.866776,  0.375774, -0.164595, &
     -0.904859, -0.222222, -0.067258, &
     -0.033908,  0.944118, -0.164595, &
      0.731143,  0.375774, -0.493784, &
      0.804597, -0.222222, -0.419426, &
     -0.443866,  0.650456, -0.502561, &
     -0.086406, -0.830237, -0.419426, &
     -0.616466, -0.543825, -0.493784, &
      0.209088,  0.650456, -0.637077, &
      0.371103, -0.543825, -0.697234, &
     -0.681916,  0.042441, -0.637077, &
      0.374587,  0.042441, -0.854729, &
     -0.185551,  0.271913, -0.900685, &
     -0.176084, -0.333333, -0.854729 /), (/ 3, 32 /) )

  return
end
subroutine soccer_vertices ( x )

!*****************************************************************************80
!
!! SOCCER_VERTICES returns the vertices of the truncated icosahedron in 3D.
!
!  Discussion:
!
!    The shape is a truncated icosahedron, which is the design used
!    on a soccer ball.  There are 12 pentagons and 20 hexagons.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X(3,60), the coordinates of the 60 vertices of
!    a truncated icosahedron.
!
  implicit none

  real ( kind = 8 ) x(3,60)

  x(1:3,1:60) = reshape ( &
    (/ -1.00714,    0.153552,   0.067258, &
       -0.960284,   0.0848813, -0.33629,  &
       -0.95172,   -0.153552,   0.33629,  &
       -0.860021,   0.529326,   0.150394, &
       -0.858,     -0.290893,  -0.470806, &
       -0.849436,  -0.529326,   0.201774, &
       -0.802576,  -0.597996,  -0.201774, &
       -0.7842,     0.418215,  -0.502561, &
       -0.749174,  -0.0848813,  0.688458, &
       -0.722234,   0.692896,  -0.201774, &
       -0.657475,   0.597996,   0.502561, &
       -0.602051,   0.290893,   0.771593, &
       -0.583675,  -0.692896,   0.470806, &
       -0.579632,  -0.333333,  -0.771593, &
       -0.52171,   -0.418215,   0.771593, &
       -0.505832,   0.375774,  -0.803348, &
       -0.489955,  -0.830237,  -0.33629,  &
       -0.403548,   0.,        -0.937864, &
       -0.381901,   0.925138,  -0.201774, &
       -0.352168,  -0.666667,  -0.688458, &
       -0.317142,   0.830237,   0.502561, &
       -0.271054,  -0.925138,   0.33629,  &
       -0.227464,   0.333333,   0.937864, &
       -0.224193,  -0.993808,  -0.067258, &
       -0.179355,   0.993808,   0.150394, &
       -0.165499,   0.608015,  -0.803348, &
       -0.147123,  -0.375774,   0.937864, &
       -0.103533,   0.882697,  -0.502561, &
       -0.0513806,  0.666667,   0.771593, &
        0.0000000,  0.,         1.021,    &
        0.0000000,  0.,        -1.021,    &
        0.0513806, -0.666667,  -0.771593, &
        0.103533,  -0.882697,   0.502561, &
        0.147123,   0.375774,  -0.937864, &
        0.165499,  -0.608015,   0.803348, &
        0.179355,  -0.993808,  -0.150394, &
        0.224193,   0.993808,   0.067258, &
        0.227464,  -0.333333,  -0.937864, &
        0.271054,   0.925138,  -0.33629,  &
        0.317142,  -0.830237,  -0.502561, &
        0.352168,   0.666667,   0.688458, &
        0.381901,  -0.925138,   0.201774, &
        0.403548,   0.,         0.937864, &
        0.489955,   0.830237,   0.33629,  &
        0.505832,  -0.375774,   0.803348, &
        0.521710,   0.418215,  -0.771593, &
        0.579632,   0.333333,   0.771593, &
        0.583675,   0.692896,  -0.470806, &
        0.602051,  -0.290893,  -0.771593, &
        0.657475,  -0.597996,  -0.502561, &
        0.722234,  -0.692896,   0.201774, &
        0.749174,   0.0848813, -0.688458, &
        0.784200,  -0.418215,   0.502561, &
        0.802576,   0.597996,   0.201774, &
        0.849436,   0.529326,  -0.201774, &
        0.858000,   0.290893,   0.470806, &
        0.860021,  -0.529326,  -0.150394, &
        0.951720,   0.153552,  -0.33629,  &
        0.960284,  -0.0848813,  0.33629,  &
        1.007140,  -0.153552,  -0.067258 /), (/ 3, 60 /) )

  return
end
subroutine sphere_unit_haltons_3d ( n, x )

!*****************************************************************************80
!
!! SPHERE_UNIT_HALTONS_3D picks a Halton point on the unit sphere in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X(3), the sample point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) phi(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) u(2,n)
  integer ( kind = 4 ) value(1)
  real ( kind = 8 ) vdot(n)
  real ( kind = 8 ) x(3,n)

  value(1) = 0
  call halton_memory ( 'SET', 'SEED', 1, value )

  value(1) = 3
  call halton_memory ( 'SET', 'DIM_NUM', 1, value )

  call halton_vector_sequence ( 2, n, u )
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  vdot(1:n) = 2.0D+00 * u(1,1:n) - 1.0E+00

  phi(1:n) = acos ( vdot(1:n) )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  theta(1:n) = 2.0D+00 * pi * u(2,1:n)

  x(1,1:n) = cos ( theta(1:n) ) * sin ( phi(1:n) )
  x(2,1:n) = sin ( theta(1:n) ) * sin ( phi(1:n) )
  x(3,1:n) = cos ( phi(1:n) )

  return
end
subroutine sphere_unit_samples_3d ( n, x )

!*****************************************************************************80
!
!! SPHERE_UNIT_SAMPLES_3D picks a random point on the unit sphere in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Output, real ( kind = 8 ) X(3,N), the random points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) phi(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) x(3,n)
!
!  Pick a uniformly random value between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  call random_number ( harvest = phi(1:n) )
  phi(1:n) = 2.0D+00 * phi(1:n) - 1.0D+00
  phi(1:n) = acos ( phi(1:n) )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  call random_number ( harvest = theta(1:n) )
  theta(1:n) = 2.0D+00 * pi * theta(1:n)

  x(1,1:n) = cos ( theta(1:n) ) * sin ( phi(1:n) )
  x(2,1:n) = sin ( theta(1:n) ) * sin ( phi(1:n) )
  x(3,1:n) = cos ( phi(1:n) )

  return
end
subroutine sphere_unit_spiralpoints_3d ( n, x )

!*****************************************************************************80
!
!! SPHERE_UNIT_SPIRALPOINTS_3D produces spiral points on the unit sphere in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Saff, Arno Kuijlaars,
!    Distributing Many Points on a Sphere,
!    The Mathematical Intelligencer,
!    Volume 19, Number 1, 1997, pages 5-11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points to create.
!
!    Output, real ( kind = 8 ) X(3,N), the coordinates of the points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cosphi
  integer ( kind = 4 ) i
  real ( kind = 8 ) sinphi
  real ( kind = 8 ) theta
  real ( kind = 8 ), parameter :: two_pi = 2.0D+00 * 3.141592653589793D+00
  real ( kind = 8 ) x(3,n)

  do i = 1, n

    cosphi = ( real ( n - i,     kind = 8 ) * ( -1.0D+00 )   &
             + real (     i - 1, kind = 8 ) * ( +1.0D+00 ) ) &
             / real ( n     - 1, kind = 8 )

    sinphi = sqrt ( 1.0D+00 - cosphi**2 )

    if ( i == 1 .or. i == n ) then
      theta = 0.0D+00
    else
      theta = theta + 3.6D+00 / ( sinphi * sqrt ( real ( n, kind = 8 ) ) )
      theta = mod ( theta, two_pi )
    end if

    x(1,i) = sinphi * cos ( theta )
    x(2,i) = sinphi * sin ( theta )
    x(3,i) = cosphi

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

  character ( len = 8  ) ampm
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
subroutine voronoi_vertices ( nc, centers, nv, vertices )

!*****************************************************************************80
!
!! VORONOI_VERTICES returns the vertices of a Voronoi diagram.
!
!  Discussion:
!
!    A set of NC points on the unit sphere is given.
!
!    The appropriate STRIPACK routines are called to generate the
!    Voronoi diagram of the points.
!
!    The vertices of this diagram are returned.
!
!    If we start with NC center points, then the Delaunay triangulation
!    will have:
!
!    * NC vertices
!    * the number of edges will be E = ( 3 * F ) / 2 because each
!      triangle adds 3 sides, but each side is added twice.
!    * Euler's formula will read F - ( 3 * F ) / 2 + NC = 2
!
!    Therefore, the number of triangles will be:
!
!      F = 2 * ( NC - 2 )
!
!    But each Delaunay triangle corresponds to a Voronoi vertex.  
!
!    If we add the NV vertices to our NC = NC(0) center points to create a
!    refined mesh of NC(1) points, then:
!
!      NC(1) = NC(0) + 2 * ( NC(0) - 2 ) = 3 * ( NC(0) - 2 ) + 2
!
!    which can be rewritten as:
!
!      NC(1) - 2 = 3 * ( NC(0) - 2 )
!
!    and if this process is repeated K times,
!
!      NC(K) - 2 = 3**K * ( NC(0) - 2 )
! 
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NC, the number of "center" points.
!
!    Input, real ( kind = 8 ) CENTERS(3,NC), the coordinates of the
!    center points.
!
!    Output, integer ( kind = 4 ) NV, the number of vertices.  NV = 2 * NC - 4.
!
!    Output, real ( kind = 8 ) VERTICES(3,NV), the coordinates of the vertices.
!
  implicit none

  integer ( kind = 4 ) nc

  real ( kind = 8 ) centers(3,nc)
  real ( kind = 8 ) ds(nc)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iwk(2*nc)
  integer ( kind = 4 ) lbtri(6,nc)
  integer ( kind = 4 ) lend(nc)
  integer ( kind = 4 ) list(6*nc)
  integer ( kind = 4 ) listc(6*nc)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lptr(6*nc)
  integer ( kind = 4 ) ltri(9,2*nc-4)
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nv
  real ( kind = 8 ) rc(2*nc-4)
  real ( kind = 8 ) vertices(3,2*nc-4)
  real ( kind = 8 ) x(nc)
  real ( kind = 8 ) xc(2*nc-4)
  real ( kind = 8 ) y(nc)
  real ( kind = 8 ) yc(2*nc-4)
  real ( kind = 8 ) z(nc)
  real ( kind = 8 ) zc(2*nc-4)
!
!  Copy the points out.
!
  x(1:nc) = centers(1,1:nc)
  y(1:nc) = centers(2,1:nc)
  z(1:nc) = centers(3,1:nc)
!
!  Create the triangulation.
!
  call trmesh ( nc, x, y, z, list, lptr, lend, lnew, iwk, iwk(nc+1), ds, &
    ierror )

  if ( ierror == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_VERTICES - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first 3 nodes are collinear.'
    stop
  end if

  if ( 0 < ierror ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_VERTICES - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Create a triangle list.
!
  call trlist ( nc, list, lptr, lend, 9, nt, ltri, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_VERTICES - Fatal error!'
    write ( *, '(a)' ) '  Error in TRLIST.'
    stop
  end if
!
!  Construct the Voronoi diagram.
!
!  Note that the triangulation data structure is altered if 0 < NB.
!
  call crlist ( nc, nc, x, y, z, list, lend, lptr, lnew, &
    lbtri, listc, nb, xc, yc, zc, rc, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_VERTICES - Fatal error!'
    write ( *, '(a)' ) '  Error in CRLIST.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Pack up the vertices.
!
  nv = 2 * nc - 4

  vertices(1,1:nv) = xc(1:nv)
  vertices(2,1:nv) = yc(1:nv)
  vertices(3,1:nv) = zc(1:nv)

  return
end
