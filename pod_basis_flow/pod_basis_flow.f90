program main

!*****************************************************************************80
!
!! MAIN is the main program for POD_BASIS_FLOW.
!
!  Discussion:
!
!    POD_BASIS_FLOW extracts a "good" basis from a set of data.
!
!    What we really want to do is take a thousand "points" and put
!    them into 8 or 9 clusters.  Each cluster will be represented
!    by a special generator point (which need not be one of the original
!    data points).  The cost of a particular clustering is the
!    sum of the squares of the distances of each data point to its
!    generator point.
!
!    The goal of this particular application is then to use the
!    generator points as a good basis for representing the dynamics
!    of the system that generated the original data.
!
!    Or, in other words, we're going to solve a PDE, save a bunch
!    of solutions, pick out a few representative ones, and assume
!    that the most interesting behavior of solutions to the PDE
!    is incorporated in the representatives.
!
!
!    The data to be examined is assumed to be stored in a file.
!
!    The file is assumed to contain a number of records, with each
!    record stored on its own line.
!
!    Each record, in turn, contains a fixed number of data values,
!    and can be regarded as a vector, or a point in M dimensional space.
!
!    The program will try to cluster the data, that is, to organize
!    the data by defining a number of cluster centers, which are
!    also points in M dimensional space, and assigning each record
!    to the cluster associated with a particular center.
!
!    The method of assigning data aims to minimize the cluster energy,
!    which is taken to be the sum of the squares of the distances of
!    each data point from its cluster center.
!
!    In some contexts, it makes sense to use the usual Euclidean sort
!    of distance.  In others, it may make more sense to replace each
!    data record by a normalized version, and to assign distance
!    by computing angles between the unit vectors.
!
!  Diary:
!
!    08 October 2004: Janet requests the ability to input multiple
!    sets of input files, as was done earlier in CVT_BASIS_FLOW.
!    In order to do this, I decided to copy as much of the internal
!    coding of CVT_BASIS_FLOW as possible.  This meant, in particular,
!    that the range and definition of the values of RUN, for instance,
!    changed entirely.  In the long run, it's better for the two codes
!    to agree about such matters.
!
!    13 September 2004: Janet requests a "no comment" option on the
!    output files.
!
!    17 July 2004: Janet requests a value of RUN_TYPE which indicates
!    that the user has already done the preprocessing, so that no
!    steady state solution has to be subtracted by the program.
!    This will be RUN_TYPE 3.
!
!    07 July 2004: Max and Lee want CPU timings for the 4, 8 and 16
!    vector cases, with the T-Cell data.  Timing only the "interesting"
!    part of the calculation, and running on the Alpha, I get:
!
!       4: 108.057 seconds
!       8: 108.354 seconds
!      16: 108.399 seconds
!
!    18 July 2003: Max wants to look at runs in which half the
!    data vectors are dropped.  I called this RUN_TYPE 2.
!
!    10 July 2003: We agreed that, if the mass matrix M is used
!    to precondition the SVD system ( X -> L' * X ), then the
!    vectors extracted from the SVD system must be premultiplied
!    by the inverse of L' before output!  Lee first noticed this
!    when he pointed out that the data vectors were zero at boundary
!    nodes, but the output vectors were not (but should be).
!
!    08 July 2003: Converted to double precision at Lee's request.
!
!    05 June 2003: I wasn't working so hard on this code for a while,
!    partly because of a surge of work on quasirandom datasets, and
!    partly because I did not have any confidence in the mass matrix
!    calculations.  Today, I corrected a missing factor of 2 in some
!    of the basis functions, and compared the mass matrix computed by
!    FEMPACK for the same sample 9 node problem, and they are the same,
!    so I can get back on the job again.
!
!    25 April 2003: I finally inserted the last bit of code that
!    multiplies (optionally) the snapshot matrix by the Cholesky
!    factor of the mass matrix.  Assuming I computed everything
!    correctly, I'm hoping it makes little actual difference.
!    Opportunities for missteps abound, starting with the fact
!    that Lee and I differ on the convention for ordering the nodes
!    of a 6 node quadratic triangle, that he sent me element data
!    for a different resolution of the TCELL problem than I had
!    solution data for; proceeding through my attempts to reconstruct
!    what I imagine his quadratic basis function evaluation looks like,
!    sailing on through the tricks of storing, factoring, transposing
!    and premultiplying by the mass matrix and its Choleksy factor...
!    All for...probably minimal change...if I'm lucky.
!
!    30 March 2003: Lee sent a memo.  From what I can make of it
!    so far, he wants to consider the problem in terms of a finite
!    element stiffness matrix; instead of analyzing snapshots of
!    solution values, we want to analyze snapshots of solution
!    coefficient vectors.  For our finite elements, aren't these
!    the same thing?
!
!    22 March 2003: I'm able to get vectors.  Today, I pulled out
!    16 POD vectors from the TCELL data, and I'm telling H C Lee
!    about it.
!
!    15 March 2003: Max wants a POD_BASIS_FLOW code similar to CVT_BASIS.
!    The idea is to take as input the same set of PDE solutions,
!    CAVITY, INOUT, and TCELL, compute the singular value decomposition,
!    and output (some of) the singular values and associated left
!    singular vectors.
!
!    Also, Max says do NOT normalize the data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford,
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum,
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Third Edition,
!    SIAM, 1999,
!    LC: QA76.73.F25L36
!
!    John Burkardt, Max Gunzburger, Hyung-Chun Lee,
!    Centroidal Voronoi Tessellation-Based Reduced-Order
!    Modelling of Complex Systems,
!    SIAM Journal on Scientific Computing,
!    Volume 28, Number 2, 2006, pages 459-484.
!
!    Gal Berkooz, Philip Holmes, John Lumley,
!    The proper orthogonal decomposition in the analysis
!    of turbulent flows,
!    Annual Review of Fluid Mechanics,
!    Volume 25, 1993, pages 539-575.
!
!    Lawrence Sirovitch,
!    Turbulence and the dynamics of coherent structures, Parts I-III,
!    Quarterly of Applied Mathematics,
!    Volume XLV, Number 3, 1987, pages 561-590.
!
  implicit none

  integer ( kind = 4 ), parameter :: uv0_file_max = 10

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  integer ( kind = 4 ) bandwidth
  character ( len = 80 ) basis_file_name
  integer ( kind = 4 ) basis_num
  logical comment
  character comment_char
  integer ( kind = 4 ) dim_num
  character ( len = 100 ) element_file_name
  integer ( kind = 4 ) element_num
  logical file_exist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ) norm
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrhs
  real ( kind = 8 ), allocatable, dimension ( :, :) :: point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) run_type
  logical s_eqi
  character ( len = 11 ) s_of_i4
  character ( len = 100 ) steady_file_name
  real ( kind = 8 ) steady_max
  real ( kind = 8 ) steady_norm
  real ( kind = 8 ), allocatable, dimension ( : ) :: sval
  real ( kind = 4 ) t1
  real ( kind = 4 ) t2
  integer ( kind = 4 ) thin
  real ( kind = 4 ) time_mass
  real ( kind = 4 ) time_svd
  real ( kind = 8 ), allocatable, dimension ( : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: u_steady
  character ( len = 100 ) uv_file_name
  integer ( kind = 4 ) uv_file_num
  character ( len = 100 ) uv0_file_name(uv0_file_max)
  integer ( kind = 4 ) uv0_file_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: v
  real ( kind = 8 ), allocatable, dimension ( : ) :: v_steady
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  character ( len = 100 ) xy_file_name
  integer ( kind = 4 ) xy_lines
  integer ( kind = 4 ) xy_values
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min

  call timestamp ( )

  time_mass = 0.0D+00
  time_svd = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POD_BASIS_FLOW'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Given a PDE for which:'
  write ( *, '(a)' ) '    M is the dimension of each solution vector,'
  write ( *, '(a)' ) '    N is the number of solution vectors,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Set up A, the M by N matrix of solution vectors,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get A = U * S * V'', the singular value decomposition.'
  write ( *, '(a)' ) ' '
!
!  Get the run type
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The variable RUN_TYPE determines preprocessing:'
  write ( *, '(a)' ) '* 1, NO steady state file, no preprocessing;'
  write ( *, '(a)' ) '* 2, NO steady state file, no preprocessing;'
  write ( *, '(a)' ) '* 3,    steady state file;'
  write ( *, '(a)' ) '    subtract 1/3 SS from solution  1'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 2 through 201'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 202 through 401.'
  write ( *, '(a)' ) '* 4,    steady state file;'
  write ( *, '(a)' ) '    subtract 1/3 SS from solution  1'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 2 through 201'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 202 through 401.'
  write ( *, '(a)' ) '    discard the EVEN numbered data files.'
  write ( *, '(a)' ) '* 5,    steady state file;'
  write ( *, '(a)' ) '    subtract 1/3 SS from solution  1'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 2 through 201'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 202 through 401.'
  write ( *, '(a)' ) '    normalize the data.'
  write ( *, '(a)' ) '* 6,    steady state file;'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 1 through 250'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 251 through 500.'
  write ( *, '(a)' ) '    do NOT normalize the result.'
  write ( *, '(a)' ) '* 7,    steady state file;'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 1 through 250'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 251 through 500.'
  write ( *, '(a)' ) '    NORMALIZE the result.'
  write ( *, '(a)' ) '    discard the ODD numbered data files.'
  write ( *, '(a)' ) '* 8,    steady state file;'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 1 through 250'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 251 through 500.'
  write ( *, '(a)' ) '    do NOT normalize the result.'

  call i4_input ( 'What is the run type (1-8)?', run_type, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the run type.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  RUN_TYPE = ', run_type
!
!  What is the basis size?
!
  call i4_input ( 'What is the requested basis size?', basis_num, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the basis size.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  BASIS_NUM = ', basis_num
!
!  Get the XY data file.
!
  call s_input ( '  What is the XY data file name?', xy_file_name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the XY data file name.'
    stop
  end if

  call data_size ( xy_file_name, xy_lines, xy_values, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  Error reading the XY data file.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The file "' // trim ( xy_file_name ) // '" contains ' &
    // trim ( s_of_i4 ( xy_lines ) ) // ' lines.'
!
!  Allocate space for some arrays.
!
  node_num = xy_lines

  allocate ( u(1:node_num) )
  allocate ( u_steady(1:node_num) )
  allocate ( v(1:node_num) )
  allocate ( v_steady(1:node_num) )
  allocate ( x(1:node_num) )
  allocate ( y(1:node_num) )
!
!  Read in X and Y.
!
  call data_d2_read ( xy_file_name, node_num, x, y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  Error reading from the XY data file.'
    stop
  end if
!
!  Extract some interesting data.
!
  x_min = minval ( x(1:node_num) )
  x_max = maxval ( x(1:node_num) )
  y_min = minval ( y(1:node_num) )
  y_max = maxval ( y(1:node_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X minimum : ', x_min
  write ( *, '(a,g14.6)' ) '  X maximum : ', x_max
  write ( *, '(a,g14.6)' ) '  Y minimum : ', y_min
  write ( *, '(a,g14.6)' ) '  Y maximum : ', y_max
!
!  Get the steady state file name.
!
  call s_input ( 'What is the name of the steady state file or "none"?', &
    steady_file_name, ierror )

  if ( s_eqi ( steady_file_name, 'none' ) ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  No steady state file.'

  else
!
!  Read in the steady state solution.
!
    call data_d2_read ( steady_file_name, node_num, u_steady, v_steady, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
      write ( *, '(a)' ) '  Error reading from the steady state file.'
      stop
    end if

    steady_max = &
      maxval ( sqrt ( u_steady(1:node_num)**2 + v_steady(1:node_num)**2 ) )
    steady_norm = &
      sqrt ( sum ( u_steady(1:node_num)**2 + v_steady(1:node_num)**2 ) )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Steady state information was read from'
    write ( *, '(a)' ) '  the file "' // trim ( steady_file_name ) // '".'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Steady max norm = ', steady_max
    write ( *, '(a,g14.6)' ) '  Steady l2 norm =  ', steady_norm

  end if

  uv_file_num = 0
  uv0_file_num = 0
  i = 0
!
!  Get the UV0 file name.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It is time to read sets of solution files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We assume that a series of solution files exists,'
  write ( *, '(a)' ) '  with "consecutive" names, like'
  write ( *, '(a)' ) '    fred001.txt, fred002,txt, ...'
  write ( *, '(a)' ) '  Just specify the FIRST name in the series, and'
  write ( *, '(a)' ) '  the program will read them all.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program can read another series of data if'
  write ( *, '(a)' ) '  you specify yet another first name, or you can'
  write ( *, '(a)' ) '  type "none" when there are no more file series.'
  write ( *, '(a)' ) ' '

  do while ( uv0_file_num < uv0_file_max )

    i = i + 1

    if ( i == 1 ) then
      call s_input ( 'What is the first solution file (in the first series)?', &
        uv0_file_name(i), ierror )
    else
      call s_input ( &
        'What is the first solution file (in the NEXT series) or "none"?', &
        uv0_file_name(i), ierror )
    end if

    if ( s_eqi ( uv0_file_name(i), 'none' ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The user indicated the end of the series.'
      exit
    end if

    uv0_file_num = uv0_file_num + 1
!
!  Presumably, all the solution files have the same name as the first
!  solution file, but with a numerical increment.  To begin with, simply count
!  the number of files.
!
    uv_file_name = uv0_file_name(i)

    do

      if ( .not. file_exist ( uv_file_name ) ) then
        exit
      end if

      uv_file_num = uv_file_num + 1

      call file_name_inc ( uv_file_name )

    end do

  end do

  if ( uv_file_num == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  There do not seem to be any solution files.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  The number of initial solution files is ', uv0_file_num
  write ( *, '(a,i8)' ) &
    '  The total number of solution files is ', uv_file_num
!
!  Now we have enough information to set up a data structure.
!
!  Determine the spatial dimension (columns) and number of points (rows).
!
  dim_num = 2 * node_num

  point_num = uv_file_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data is stored in an M by N matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The "spatial" dimension M is   ', dim_num
  write ( *, '(a,i8)' ) '  The number of data points N is ', point_num
!
!  Allocate space for the POINT array.
!
  allocate ( point(1:dim_num,1:point_num) )
!
!  Now read the data from the individual files, process it if necessary,
!  and gather it into a single array called POINT.
!
  write ( *, '(a)' ) ' '
  j = 0

  do i = 1, uv0_file_num

    uv_file_name = uv0_file_name(i)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Processing files starting with "' &
      // trim ( uv0_file_name(i) ) // '".'

    do

      if ( .not. file_exist ( uv_file_name ) ) then
        exit
      end if

      j = j + 1
      if ( j == 1 ) then
        write ( *, '(2x,a,i8,2x,a)' ) 'Reading file ', j, trim ( uv_file_name )
      else if ( j == uv_file_num ) then
        write ( *, '(2x,a,i8,2x,a)' ) 'Reading file ', j, trim ( uv_file_name )
      end if

      call data_d2_read ( uv_file_name, node_num, u, v, ierror )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
        write ( *, '(a)' ) '  Error reading from a UV file.'
        stop
      end if
!
!  For RUN_TYPE = 3 (or 4 or 5), H C Lee wants to
!  subtract 1/3 of SS from solution  1.
!  subtract 5/3 of SS from solutions 2-201,
!  subtract 1/3 of SS from solutions 202-401.
!
      if ( run_type == 3 .or. run_type == 4 .or. run_type == 5 ) then

        if ( 1 <= j .and. j <= 1 ) then
          u(1:node_num) = u(1:node_num) - u_steady(1:node_num) / 3.0D+00
          v(1:node_num) = v(1:node_num) - v_steady(1:node_num) / 3.0D+00
        else if ( 2 <= j .and. j <= 201 ) then
          u(1:node_num) = u(1:node_num) &
            - 5.0D+00 * u_steady(1:node_num) / 3.0D+00
          v(1:node_num) = v(1:node_num) &
            - 5.0D+00 * v_steady(1:node_num) / 3.0D+00
        else if ( 202 <= j .and. j <= 401 ) then
          u(1:node_num) = u(1:node_num) - u_steady(1:node_num) / 3.0D+00
          v(1:node_num) = v(1:node_num) - v_steady(1:node_num) / 3.0D+00
        end if
!
!  For RUN_TYPE = 6 or 7 or 8, H C Lee wants to
!  subtract 5/3 of SS from solutions 1-250,
!  subtract 1/3 of SS from solutions 251-500.
!
      else if ( run_type == 6 .or. run_type == 7 .or. run_type == 8 ) then

        if ( 1 <= j .and. j <= 250 ) then
          u(1:node_num) = u(1:node_num) &
            - 5.0D+00 * u_steady(1:node_num) / 3.0D+00
          v(1:node_num) = v(1:node_num) &
            - 5.0D+00 * v_steady(1:node_num) / 3.0D+00
        else if ( 251 <= j .and. j <= 500 ) then
          u(1:node_num) = u(1:node_num) - u_steady(1:node_num) / 3.0D+00
          v(1:node_num) = v(1:node_num) - v_steady(1:node_num) / 3.0D+00
        end if

      end if

      point(1:2*node_num-1:2,j) = u(1:node_num)
      point(2:2*node_num  :2,j) = v(1:node_num)
!
!  19 June 2002: At Max's request, normalize each vector after it's
!  been processed.
!
      if ( run_type == 5 .or. run_type == 7 ) then

        norm = sqrt ( sum ( point(1:2*node_num,j)**2 ) )

        if ( 0.0D+00 < norm ) then
          point(1:2*node_num,j) = point(1:2*node_num,j) / norm
        end if

      end if

      call file_name_inc ( uv_file_name )

    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  All the data has been read into POINT.'
!
!  For RUN_TYPE = 4, thin the data by dropping EVEN values.
!
  if ( run_type == 4 ) then

    thin = 2

    point_num2 = 1 + ( point_num - 1 ) / thin

    point(1:2*node_num,1:point_num2) = point(1:2*node_num,1:point_num:thin)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  RUN_TYPE = 4:'
    write ( *, '(a)' ) '  Thin out the input data points.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Thinning increment is       ', thin
    write ( *, '(a,i8)' ) '  Original input data size is ', point_num
    write ( *, '(a,i8)' ) '  Thinned data size is        ', point_num2

    point_num = point_num2
!
!  For RUN_TYPE = 8, thin the data by dropping ODD values.
!
  else if ( run_type == 8 ) then

    thin = 2

    point_num2 = point_num / thin

    point(1:2*node_num,1:point_num2) = point(1:2*node_num,2:point_num:thin)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  RUN_TYPE = 8:'
    write ( *, '(a)' ) '  Thin out the input data points.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Thinning increment is       ', thin
    write ( *, '(a,i8)' ) '  Original input data size is ', point_num
    write ( *, '(a,i8)' ) '  Thinned data size is        ', point_num2

    point_num = point_num2

  end if
!
!----------------------------------------------------------------------------
!
!  Precondition with the mass matrix M?
!
!  If so, instead of the system X' * X, we wish to study the
!  system X' * M * X.
!
!  However, it's easy to apply either eigenanalsis to X'*X or
!  SVD to X, however, if we are dealing with the system X' * M * X,
!  it's easy to apply eigenanalysis, but not obvious what system
!  we should apply SVD to.
!
!  In order to make this look like the simpler problem, we need
!  to split M up into appropriate factors and regroup.
!
!  Thus, we will:
!  * compute the (symmetric positive definite) mass matrix M,
!  * determine the Cholesky factorization M = L * L',
!  * premultiply the data matrix X by L'.
!
!  Thus, we regard the system
!    X' * M * X
!  as the system
!    X' * L * L' * X
!  or
!    (L'*X)' * (L'*X)
!
!  Now it's obvious that we will apply SVD to L'*X.
!
!  THEN (10 July) we must postprocess the vectors we extract from
!  the columns of U by essentially multiplying them by inverse (L').
!
!  So if we overwrite X by L'X, we can proceed as though nothing
!  is really different.
!
!----------------------------------------------------------------------------
!
  call s_input ( &
    '  Enter element file for mass matrix preconditioning or "None".', &
    element_file_name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the element file name.'
    stop
  end if

  if ( s_eqi ( element_file_name, 'None' ) ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The user indicated there is no element file.'

  else

    call data_size ( element_file_name, element_num, npe, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
      write ( *, '(a)' ) '  Input error reading the element file.'
      stop
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of elements = ', element_num
    write ( *, '(a,i8)' ) '  Number of nodes per element = ', npe

    allocate ( node(1:npe,1:element_num) )

    call data_ivec_read ( element_file_name, npe, element_num, node, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
      write ( *, '(a)' ) '  Error reading the element file.'
      stop
    end if

    call bandwidth_determine ( npe, element_num, node, bandwidth )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The bandwidth of the matrix is ', bandwidth
!
!  Allocate storage for A.
!
    allocate ( a(bandwidth+1,node_num) )
!
!  Compute A.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Compute the mass matrix.'

    call cpu_time ( t1 )

    call mass_matrix ( node_num, npe, element_num, node, bandwidth, x, y, a )

    call cpu_time ( t2 )
    time_mass = time_mass + t2 - t1
!
!  Print the mass matrix, if small.
!
    if ( node_num < 10 ) then
      call dpbl_print ( node_num, bandwidth, a, '  The mass matrix.' )
    end if

    deallocate ( node )
!
!  Factor A into A = L * L'.
!
!  L is a band lower triangular matrix, with the same
!  bandwidth as A.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Cholesky-factor the mass matrix.'

    call cpu_time ( t1 )

    call dpbtrf ( 'L', node_num, bandwidth, a, bandwidth+1, info )

    call cpu_time ( t2 )
    time_mass = time_mass + t2 - t1

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POD_BASIS_FLOW - Fatal error!'
      write ( *, '(a)' ) '  DPBTRF reports the mass matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      stop
    end if

    if ( node_num < 10 ) then
      call dblt_print ( node_num, bandwidth, a, '  The L factor:' )
    end if
!
!  Premultiply POINT by L'.
!  Each column of POINT contains U and V coefficients.
!  We need to delicately extract U, compute L'*U, and put it back;
!  then repeat for V.
!
    do j = 1, point_num

      u(1:node_num) = point(1:2*node_num-1:2,j)

      if ( node_num < 10 ) then
        call r8vec_print ( node_num, u, '  U' )
      end if

      call cpu_time ( t1 )

      call dtbmv ( 'L', 'T', 'N', node_num, bandwidth, a, bandwidth+1, u, 1 )

      call cpu_time ( t2 )
      time_mass = time_mass + t2 - t1

      if ( node_num < 10 ) then
        call r8vec_print ( node_num, u, '  L''*U' )
      end if

      point(1:2*node_num-1:2,j) = u(1:node_num)

      v(1:node_num) = point(2:2*node_num  :2,j)

      if ( node_num < 10 ) then
        call r8vec_print ( node_num, v, '  V' )
      end if

      call cpu_time ( t1 )

      call dtbmv ( 'L', 'T', 'N', node_num, bandwidth, a, bandwidth+1, v, 1 )

      call cpu_time ( t2 )
      time_mass = time_mass + t2 - t1

      point(2:2*node_num  :2,j) = v(1:node_num)

      if ( node_num < 10 ) then
        call r8vec_print ( node_num, v, '  L''*V' )
      end if

    end do

  end if
!
!----------------------------------------------------------------------------
!
!  Compute the SVD of the data.
!
!----------------------------------------------------------------------------
!
  allocate ( sval(1:basis_num) )

  call cpu_time ( t1 )

  call singular_vectors ( dim_num, point_num, basis_num, point, sval )

  call cpu_time ( t2 )
  time_svd = time_svd + t2 - t1
!
!----------------------------------------------------------------------------
!
!  If mass matrix preconditioning was used, we have to
!  "unprecondition" now.  Each left singular vector is decomposable
!  into a "U" and "V" portion.  Both the "U" and "V" portion were
!  multiplied by the Cholesky factor L' before the SVD was applied.
!  Therefore, to get back to the original data world, we now have
!  to reverse that transformation.
!
!----------------------------------------------------------------------------
!
  if ( s_eqi ( element_file_name, 'None' ) ) then

  else
!
!  Premultiply POINT by inverse ( L' ).
!  Each column of POINT contains U and V coefficients.
!  We need to delicately extract U, solve L'*U(new)=U(old), and put it back;
!  then repeat for V.
!
    nrhs = 1

    call cpu_time ( t1 )

    do j = 1, basis_num

      u(1:node_num) = point(1:2*node_num-1:2,j)

      call dtbtrs ( 'L', 'T', 'N', node_num, bandwidth, nrhs, a, &
        bandwidth+1, u, node_num, info )

      point(1:2*node_num-1:2,j) = u(1:node_num)

      v(1:node_num) = point(2:2*node_num  :2,j)

      call dtbtrs ( 'L', 'T', 'N', node_num, bandwidth, nrhs, a, &
        bandwidth+1, v, node_num, info )

      point(2:2*node_num  :2,j) = v(1:node_num)

    end do

    call cpu_time ( t2 )
    time_mass = time_mass + t2 - t1

    deallocate ( a )

  end if
!
!----------------------------------------------------------------------------
!
!  Write vectors to files.
!
!----------------------------------------------------------------------------
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POD_BASIS_FLOW:'
  write ( *, '(a)' ) '  Ready to write the left singular vectors to files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Do you want comments in the header of the file?'
  write ( *, '(a)' ) '  (These begin with the "#" character.) (Y/N)'

  call s_input ( '  Enter Y or N:', comment_char, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS_FLOW - Warning!'
    write ( *, '(a)' ) '  Input error reading the comment option.'
    write ( *, '(a)' ) '  We will assume comments are acceptable.'
    comment_char = 'Y'
  end if

  if ( comment_char == 'Y' .or. comment_char == 'y' ) then
    comment = .true.
  else
    comment = .false.
  end if

  basis_file_name = 'pod_000.txt'

  do j = 1, basis_num

    call file_name_inc ( basis_file_name )

    if ( j == 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Writing first file ' // trim ( basis_file_name )
    end if

    if ( j == basis_num ) then
      write ( *, '(a)' ) '  Writing last file  ' // trim ( basis_file_name )
    end if

    call basis_write ( basis_file_name, dim_num, sval(j), point(1:dim_num,j), &
      comment )

  end do

  deallocate ( point )
  deallocate ( sval )
  deallocate ( u )
  deallocate ( u_steady )
  deallocate ( v )
  deallocate ( v_steady )
  deallocate ( x )
  deallocate ( y )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CPU time in seconds:'
  write ( *, '(a,g14.6)' ) '    For Mass matrix: ', time_mass
  write ( *, '(a,g14.6)' ) '    For SVD:         ', time_svd
  write ( *, '(a,g14.6)' ) '    Total:           ', time_mass + time_svd
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POD_BASIS_FLOW'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine bandwidth_determine ( npe, element_num, node, bandwidth )

!*****************************************************************************80
!
!! BANDWIDTH_DETERMINE computes the lower bandwidth of a finite element matrix.
!
!  Discussion:
!
!    The finite element matrix is assumed to be structured in such a
!    way that row I represents an equation associated with the unknown
!    at node I, and column I stores coefficients associated with
!    the unknown at node I as it appears in various equations.
!
!    Further, it is assumed that the I-th equation, associated with
!    node and unknown I, involves only those nodes and unknowns J
!    with the property that there is an element K that includes both
!    nodes.
!
!    Thus, the (half) bandwidth calculation simply involves finding the
!    greatest difference between two nodes in the same element and adding 1.
!
!    A diagonal matrix will have bandwidth 1.  A tridiagonal matrix
!    will have bandwidth 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPE, the number of nodes per element.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) NODE(NPE,ELEMENT_NUM), the nodes that make up each element.
!
!    Output, integer ( kind = 4 ) BANDWIDTH, the (half) bandwidth of the matrix.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) npe

  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) elem
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) node(npe,element_num)

  bandwidth = 1

  do elem = 1, element_num
    do i1 = 1, npe
      n1 = node(i1,elem)
      do i2 = 1, npe
        n2 = node(i2,elem)
        bandwidth = max ( bandwidth, n2 + 1 - n1 )
      end do
    end do
  end do

  return
end
subroutine basis_write ( file_out_name, m, s, x, comment )

!*****************************************************************************80
!
!! BASIS_WRITE writes a basis vector to a file.
!
!  Discussion:
!
!    The basis vector X(M) is assumed to represent (M/2) pairs of
!    horizontal and vertical velocities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file to write.
!
!    Input, integer ( kind = 4 ) M, the number of data items.
!
!    Input, real S, the associated singular value.
!
!    Input, real X(M), the data values.
!
!    Input, logical COMMENT, is TRUE if comments may be inserted into the file.
!
  implicit none

  integer ( kind = 4 ) m

  logical comment
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) i
  real ( kind = 8 ) s
  character ( len = 40 ) string
  real ( kind = 8 ) x(m)

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace' )

  call timestring ( string )

  if ( comment ) then
    write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
    write ( file_out_unit, '(a)'       ) '#  created by BASIS_WRITE.F90,'
    write ( file_out_unit, '(a)'       ) '#  part of program POD_BASIS_FLOW.F90.'
    write ( file_out_unit, '(a)'       ) '#'
    write ( file_out_unit, '(a)'       ) '#  Created on ' // trim ( string )
    write ( file_out_unit, '(a)'       ) '#'
    write ( file_out_unit, '(a)'       ) '#  Velocities ( U, V ).'
    write ( file_out_unit, '(a)'       ) '#'
    write ( file_out_unit, '(a,i8)'    ) '#  Number of records = ', m / 2
    write ( file_out_unit, '(a,g14.6)' ) '#  Singular value S =  ', s
    write ( file_out_unit, '(a,g14.6)' ) '#  Epsilon =           ', epsilon ( s )
    write ( file_out_unit, '(a)'       ) '#'
  end if

  do i = 1, m, 2

    write ( file_out_unit, '(2e25.15)' ) x(i), x(i+1)

  end do

  close ( unit = file_out_unit )

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT is TRUE if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine data_ivec_read ( file_in_name, m, n, a, ierror )

!*****************************************************************************80
!
!! DATA_IVEC_READ reads an dataset of integer vectors stored in a file.
!
!  Discussion:
!
!    The data set can be thought of as an integer M by N array.
!
!    Each column of the array corresponds to one data "item".
!
!    The data is stored in a file, one column at a time.
!
!    Each data item begins on a new line.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the file to read.
!
!    Input, integer ( kind = 4 ) M, the size of each data item.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Output, integer ( kind = 4 ) A(M,N), the data.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) file_in_line
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  character ( len = 80 ) line
  integer ( kind = 4 ) n2

  ierror = 0

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, iostat = ios, &
    status = 'old' )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_IVEC_READ - Fatal error!'
    write ( *, '(a)' ) '  Error opening the data file.'
    return
  end if

  n2 = 0
  file_in_line = 0

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_IVEC_READ - Fatal error!'
      write ( *, '(a)' ) '  Error reading the data file.'
      return
    end if

    file_in_line = file_in_line + 1

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      n2 = n2 + 1

      last = 0

      do i = 1, m

        call s_to_i4 ( line(last+1:), a(i,n2), ierror, length )

        if ( ierror /= 0 ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DATA_IVEC_READ - Fatal error!'
          write ( *, '(a)' ) '  Error reading data.'
          write ( *, '(a)' ) '  Reading line number ', file_in_line
          write ( *, '(a)' ) '  Data row number ', i
          write ( *, '(a)' ) '  Data column ', n2
          return
        end if

        last = last + length

      end do

      if ( n <= n2 ) then
        exit
      end if

    end if

  end do

  close ( unit = file_in_unit )

  return
end
subroutine data_d2_read ( file_name, n, x, y, ierror )

!*****************************************************************************80
!
!! DATA_D2_READ reads pairs of double precision numbers stored in a file.
!
!  Discussion:
!
!    The data set can be thought of as a real M by 2 array.
!
!    Each row of the array corresponds to one data "item".
!
!    The data is stored in a file, one row (pair of values) at a time.
!
!    Each row (pair of values) begins on a new line.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to read.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the data values.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  character ( len = 80 ) line
  integer ( kind = 4 ) n2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  ierror = 0

  call get_unit ( input )

  open ( unit = input, file = file_name, iostat = ios, status = 'old' )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_D2_READ - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  x(1:n) = huge ( x(1) )
  y(1:n) = huge ( y(1) )

  n2 = 0

  do

    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_D2_READ - Fatal error!'
      write ( *, '(a)' ) '  Error reading the file.'
      return
    end if

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      n2 = n2 + 1

      last = 0
      call s_to_r8 ( line(last+1:), x(n2), ierror, length )

      if ( ierror /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_D2_READ - Fatal error!'
        write ( *, '(a)' ) '  Error decoding information in the file.'
        return
      end if

      last = last + length

      call s_to_r8 ( line(last+1:), y(n2), ierror, length )

      if ( ierror /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_D2_READ - Fatal error!'
        write ( *, '(a)' ) '  Error decoding information in the file.'
        return
      end if

      if ( n2 == n ) then
        exit
      end if

    end if

  end do

  close ( unit = input )

  return
end
subroutine data_size ( file_name, m, n, ierror )

!*****************************************************************************80
!
!! DATA_SIZE counts the size of a data set stored in a file.
!
!  Discussion:
!
!    Blank lines and comment lines (which begin with '#') are ignored).
!
!    All other lines are assumed to represent data items.
!
!    This routine assumes that each data line contains exactly the
!    same number of values, which are separated by spaces.
!
!    (This means this routine cannot handle cases where a data item
!    extends over more than one line, or cases where data is squeezed
!    together with no spaces, or where commas are used as separators,
!    but with space on both sides.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to read.
!
!    Output, integer ( kind = 4 ) M, the number of nonblank, noncomment lines.
!
!    Output, integer ( kind = 4 ) N, the number of values per line.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max
  integer ( kind = 4 ) n_min
  integer ( kind = 4 ) n_word

  ierror = 0
  m = 0
  n_max = - huge ( n_max )
  n_min = huge ( n_min )

  call get_unit ( input )

  open ( unit = input, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '  "' // trim ( file_name ) // '".'
    stop
  end if

  do

    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      m = m + 1

      call s_word_count ( line, n_word )

      n_max = max ( n_max, n_word )
      n_min = min ( n_min, n_word )

    end if

  end do

  if ( n_max /= n_min ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Number of words per line varies.'
    write ( *, '(a,i8)' ) '  Minimum is ', n_min
    write ( *, '(a,i8)' ) '  Maximum is ', n_max
    n = 0
  else
    n = n_min
  end if

  close ( unit = input )

  return
end
subroutine dblt_check ( n, ml, ierror )

!*****************************************************************************80
!
!! DBLT_CHECK checks the dimensions of a banded lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth of the matrix.
!    ML must be at least 0, and no greater than N-1.
!
!    Output, integer ( kind = 4 ) IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if ML is illegal;
!    IERROR = IERROR + 2 if N is illegal.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  ierror = 0

  if ( ml < 0 .or. n - 1 < ml ) then
    ierror = ierror + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DBLT_CHECK - Illegal ML = ', ml
  end if

  if ( n <= 0 ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DBLT_CHECK - Illegal N = ', n
    return
  end if

  return
end
subroutine dblt_print ( n, ml, a, title )

!*****************************************************************************80
!
!! DBLT_PRINT prints a band lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!    ML must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the N by N band matrix, stored
!    in band lower triangle mode.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call dblt_print_some ( n, ml, a, 1, 1, n, n )

  return
end
subroutine dblt_print_some ( n, ml, a, ilo, jlo, ihi, jhi )

!*****************************************************************************80
!
!! DBLT_PRINT_SOME prints some of a band lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, the upper (and lower) bandwidth.
!    ML must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the N by N band lower triangular
!    matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
!
!  Check the dimensions.
!
  call dblt_check ( n, ml, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DBLT_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) 'Columns', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j <= i .and. i <= j + ml ) then
          aij = a(i-j+1,j)
        else
          aij = 0.0D+00
        end if

        if ( ml < i-j .or. 0 < j-i ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine dpbl_check ( n, mu, ierror )

!*****************************************************************************80
!
!! DPBL_CHECK checks the dimensions of a positive definite symmetric band matrix.
!
!  Discussion:
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth of the matrix.
!    MU must be at least 0, and SHOULD BE no greater than N-1.
!
!    Output, integer ( kind = 4 ) IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if MU is illegal;
!    IERROR = IERROR + 2 if N is illegal.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  ierror = 0

  if ( mu < 0 ) then
    ierror = ierror + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DPBL_CHECK - Illegal MU < 0 = ', mu
  end if

  if ( n - 1 < mu ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DPBL_CHECK - Warning!'
    write ( *, '(a,i8)' ) '  Not advisable to have N - 1 < MU!'
    write ( *, '(a,i8)' ) '  MU = ', mu
    write ( *, '(a,i8)' ) '  N =  ', n
  end if

  if ( n <= 0 ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DPBL_CHECK - Illegal N = ', n
    return
  end if

  return
end
subroutine dpbl_print ( n, mu, a, title )

!*****************************************************************************80
!
!! DPBL_PRINT prints a symmetric banded matrix.
!
!  Discussion:
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the N by N band matrix, stored
!    in positive definite symmetric band storage lower triangle mode.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call dpbl_print_some ( n, mu, a, 1, 1, n, n )

  return
end
subroutine dpbl_print_some ( n, mu, a, ilo, jlo, ihi, jhi )

!*****************************************************************************80
!
!! DPBL_PRINT_SOME prints some of a symmetric banded matrix.
!
!  Discussion:
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the N by N band matrix, stored
!    in positive definite symmetric band storage lower triangle mode.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
!
!  Check the dimensions.
!
  call dpbl_check ( n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DPBL_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) 'Columns', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + mu )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j .and. j <= i + mu ) then
          aij = a(j-i+1,i)
        else if ( j <= i .and. i <= j + mu ) then
          aij = a(i-j+1,j)
        else
          aij = 0.0D+00
        end if

        if ( mu < i-j .or. mu < j-i ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
function file_exist ( file_name )

!*****************************************************************************80
!
!! FILE_EXIST reports whether a file exists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_EXIST, is TRUE if the file exists.
!
  implicit none

  character ( len = * ) file_name
  logical file_exist

  inquire ( file = file_name, exist = file_exist )

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
subroutine i4_input ( string, value, ierror )

!*****************************************************************************80
!
!! I4_INPUT prints a prompt string and reads an integer from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#') or is
!    blank, the routine ignores that line, and tries to read the next one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, integer ( kind = 4 ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  character ( len = 80 ) line
  character ( len = * ) string
  integer ( kind = 4 ) value

  ierror = 0
  value = huge ( value )
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      return
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Extract integer information from the string.
!
    call s_to_i4 ( line, value, ierror, last )

    if ( ierror /= 0 ) then
      value = huge ( value )
      return
    end if

    exit

  end do

  return
end
subroutine i4_range_input ( string, value1, value2, ierror )

!*****************************************************************************80
!
!! I4_RANGE_INPUT reads a pair of integers from the user, representing a range.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#') or is blank,
!    the routine ignores that line, and tries to read the next one.
!
!    The pair of integers may be separated by spaces or a comma or both.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, integer ( kind = 4 ) VALUE1, VALUE2, the values entered by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  character, parameter :: comma = ','
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) last2
  character ( len = 80 ) line
  character, parameter :: space = ' '
  character ( len = * ) string
  integer ( kind = 4 ) value1
  integer ( kind = 4 ) value2

  ierror = 0
  value1 = huge ( value1 )
  value2 = huge ( value2 )
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      return
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Remove commas.
!
    call s_rep_ch ( line, comma, space )
!
!  Extract integer information from the string.
!
    call s_to_i4 ( line, value1, ierror, last )

    if ( ierror /= 0 ) then
      value1 = huge ( value1 )
      return
    end if

    call s_to_i4 ( line(last+1:), value2, ierror, last2 )

    if ( ierror /= 0 ) then
      value2 = huge ( value2 )
      return
    end if

    exit

  end do

  return
end
subroutine mass_matrix ( node_num, npe, element_num, node, bandwidth, x, y, a )

!*****************************************************************************80
!
!! MASS_MATRIX computes the mass matrix.
!
!  Discussion:
!
!    I want to compute the mass matrix associated with velocity.
!
!      A(I,J) = integral ( PHI(I)(X,Y) * PHI(J)(X,Y) ) d Region
!
!    where PHI(I) and PHI(J) are the shape functions associated with
!    the I-th and J-th variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPE, the number of nodes per element.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) NODE(NPE,ELEMENT_NUM), the nodes that make up each element.
!
!    Input, integer ( kind = 4 ) BANDWIDTH, the half bandwidth of the matrix.
!
!    Input, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the coordinates
!    of the nodes.
!
!    Output, real ( kind = 8 ) A(BANDWIDTH+1,NODE_NUM), the mass matrix.
!
  implicit none

  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) npe
  integer ( kind = 4 ), parameter :: nquad = 13

  real ( kind = 8 ) a(bandwidth+1,node_num)
  real ( kind = 8 ) area
  integer ( kind = 4 ) element
  real ( kind = 8 ) eta
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jq
  integer ( kind = 4 ) node(npe,element_num)
  integer ( kind = 4 ) norder
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) p2
  integer ( kind = 4 ) p3
  real ( kind = 8 ) refqbf
  integer ( kind = 4 ) rule
  real ( kind = 8 ) weight(nquad)
  real ( kind = 8 ) wi
  real ( kind = 8 ) wj
  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xtab(nquad)
  real ( kind = 8 ) y(node_num)
  real ( kind = 8 ) ytab(nquad)
!
!  Zero out the matrix.
!
  a(1:bandwidth+1,1:node_num) = 0.0D+00
!
!  Get the weights and abscissas for a unit triangle.
!
  rule = 12
  call triangle_unit_set ( rule, norder, xtab, ytab, weight )
!
!  For each element.
!
  do element = 1, element_num

    p1 = node(1,element)
    p2 = node(2,element)
    p3 = node(3,element)

    area = 0.5D+00 * abs ( &
        x(p1) * ( y(p2) - y(p3) ) &
      + x(p2) * ( y(p3) - y(p1) ) &
      + x(p3) * ( y(p1) - y(p2) ) )

    if ( area == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MASS_MATRIX - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero area for element ', element
      stop
    end if
!
!  For each quadrature point in the element...
!
    do iquad = 1, nquad

      xsi = xtab(iquad)
      eta = ytab(iquad)
!
!  For each basis function PHI(I) associated with a node in the element,
!
      do iq = 1, 6

        ip = node(iq,element)
        wi = refqbf ( iq, xsi, eta )
!
!  For each "neighbor" basis function PHI(J) associated with a node in
!  the element.
!
        do jq = 1, 6

          jp = node(jq,element)

          if ( jp <= ip ) then
            wj = refqbf ( jq, xsi, eta )
            a(ip+1-jp,jp) = a(ip+1-jp,jp) + area * weight(iquad) * wi * wj
          end if

        end do
      end do
    end do
  end do

  return
end
subroutine node_t6 ( r, s )

!*****************************************************************************80
!
!! NODE_T6 returns the basis nodes for a 6 node triangle.
!
!  Diagram:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(6), S(6), the coordinates of the basis nodes.
!
  implicit none

  real ( kind = 8 ) r(6)
  real ( kind = 8 ) s(6)

  r(1:6) = (/ 0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.0D+00 /)
  s(1:6) = (/ 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00 /)

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints a double precision vector.
!
!  Discussion:
!
!    If all the entries are integers, the data if printed
!    in integer format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
    do i = 1, n
      write ( *, '(i8,i8)' ) i, int ( a(i) )
    end do
  else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
    do i = 1, n
      write ( *, '(i8,f14.6)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i8,g14.6)' ) i, a(i)
    end do
  end if

  return
end
function refqbf ( iq, xsi, eta )

!*****************************************************************************80
!
!! REFQBF evaluates a reference element quadratic basis function.
!
!  Discussion:
!
!    There are six possible quadratic basis functions.  This routine
!    evaluates just one of them, and its X and Y derivatives, at a
!    particular point in a particular element, by referring to the
!    reference triangle.
!
!    The point we are interested in is referred to by its coordinates
!    in the reference triangle.  That is, we are given coordinates
!    (XSI, ETA), even though, physically, we are interested
!    in points in (X, Y) space.
!
!    Here is a graph of the (XSI, ETA) reference triangle.
!
!          ^
!          |
!      1.0 +    3
!          |    |\
!      0.5 |    6 5
!          |    |  \
!      0.0 +    1-4-2
!          |
!          +----+---+--->
!               0 0 1
!               . . .
!               0 5 0
!
!                XSI
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IQ, the local node number, between 1 and
!    6, whose basis function is being evaluated.
!
!    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
!    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!    Output, real ( kind = 8 ) REFQBF, the value of the basis function
!    PSI(IQ)(XSI,ETA).
!
  implicit none

  real ( kind = 8 ) eta
  integer ( kind = 4 ) iq
  real ( kind = 8 ) refqbf
  real ( kind = 8 ) w
  real ( kind = 8 ) xsi
!
!  W(1)(XSI,ETA) = 0 if XSI + ETA = 0.5 or XSI + ETA = 1.
!  W(1)(0.0,0.0) = 1.
!
  if ( iq == 1 ) then

    w = 2.0D+00 * ( 0.5D+00 - xsi - eta ) * ( 1.0D+00 - xsi - eta )
!
!  W(2)(XSI,ETA) = 0 if XSI=0 or XSI=0.5.
!  W(2)(1.0,0.0) = 1.
!
  else if ( iq == 2 ) then

    w = 2.0D+00 * xsi * ( xsi - 0.5D+00 )
!
!  W(3)(XSI,ETA) = 0 if ETA = 0, or ETA = 0.5.
!  W(3)(0.0,1.0) = 1.
!
  else if ( iq == 3 ) then

    w = 2.0D+00 * eta * ( eta - 0.5D+00 )
!
!  W(4)(XSI,ETA) = 0 if XSI = 0 or XSI + ETA = 1.
!  W(4)(0.5,0.0) = 1.
!
  else if ( iq == 4 ) then

    w = 4.0D+00 * xsi * ( 1.0D+00 - xsi - eta )
!
!  W(5)(XSI,ETA) = 0 if ETA = 0 or XSI = 0.
!  W(5)(0.5,0.5) = 1.
!
  else if ( iq == 5 ) then

    w = 4.0D+00 * eta * xsi
!
!  W(6)(XSI,ETA) = 0 if ETA = 0 or XSI + ETA = 1.
!  W(6)(0.0,0.5) = 1.
!
  else if ( iq == 6 ) then

    w = 4.0D+00 * eta * ( 1.0D+00 - xsi - eta )

  end if

  refqbf = w

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_input ( string, value, ierror )

!*****************************************************************************80
!
!! S_INPUT prints a prompt string and reads a string from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#'), or is blank,
!    the routine ignores that line, and tries to read the next one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, character ( len = * ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) string
  character ( len = * ) value

  ierror = 0
  value = ' '
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) value

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S_INPUT: Fatal error!'
      write ( *, '(a)' ) '  Input error!'
      stop
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( value(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( value ) == 0 ) then
      cycle
    end if

    exit

  end do

  return
end
function s_of_i4 ( i )

!*****************************************************************************80
!
!! S_OF_I4 converts an integer to a left-justified string.
!
!  Example:
!
!         I  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, an integer to be converted.
!
!    Output, character ( len = 11 ) S_OF_I4, the representation of the
!    integer ( kind = 4 ).  The integer will be left-justified.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  character ( len = 11 ) s
  character ( len = 11 ) s_of_i4

  s = ' '

  ilo = 1
  ihi = 11
!
!  Make a copy of the integer.
!
  ival = i
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do j = 1, ihi
        s(j:j) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

  s_of_i4 = s

  return
end
subroutine s_rep_ch ( s, c1, c2 )

!*****************************************************************************80
!
!! S_REP_CH replaces all occurrences of one character by another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string.
!
!    Input, character C1, C2, the character to be replaced, and the
!    replacement character.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  character ( len = * ) s

  do i = 1, len ( s )
    if ( s(i:i) == c1 ) then
      s(i:i) = c2
    end if
  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2D-1'           0.2
!    '23.45'           23.45
!    '-4.2D+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal number.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( ihave > 1 ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( ihave >= 6 .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + dble ( ndig )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + dble ( ndig )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. lchar+1 >= nchar ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

  return
end
subroutine singular_vectors ( m, n, basis_num, a, sval )

!*****************************************************************************80
!
!! SINGULAR_VECTORS computes the desired singular values.
!
!  Discussion:
!
!    The LAPACK SVD routine DGESVD is used to compute the singular
!    value decomposition:
!
!      A = U * S * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford,
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum,
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Third Edition,
!    SIAM, 1999,
!    LC: QA76.73.F25L36
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) BASIS_NUM, the number of basis vectors to be extracted.
!
!    Input, real POINT(DIM_NUM,POINT_NUM), the data points.
!
  implicit none

  integer ( kind = 4 ) basis_num
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldvt
  character jobu
  character jobvt
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) s(min(m,n))
  real ( kind = 8 ) s_partial
  real ( kind = 8 ) s_total
  real ( kind = 8 ) sval(basis_num)
  real ( kind = 8 ) u(1,1)
  real ( kind = 8 ) vt(1,1)
  real ( kind = 8 ) work(3*min(m,n)+max(max(m,n),2*min(m,n)))

  lwork = 3 * min ( m, n ) + max ( max ( m, n ), 2 * min ( m, n ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINGULAR_VECTORS:'
  write ( *, '(a)' ) '  For an MxN matrix A in general storage,'
  write ( *, '(a)' ) '  we call the LAPACK routine'
  write ( *, '(a)' ) '    DGESVD'
  write ( *, '(a)' ) '  which computes the singular value decomposition:'
  write ( *, '(a)' ) '    A = U * S * V'''
!
!  Compute the eigenvalues and eigenvectors.
!
  jobu = 'O'
  jobvt = 'N'
  lda = m
  ldu = m
  ldvt = n

  s(1:max(m,n)) = 0.0D+00

  call dgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &
    lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SINGULAR_VECTORS - Warning:'
    write ( *, '(a,i8)' ) '  DGESVD returned nonzero INFO = ', info
    return
  end if

  s_total = sum ( s(1:max(m,n)) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Singular value total = ', s_total

  sval(1:basis_num) = s(1:basis_num)
  s_partial = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Singular values, and Normalized, and Running normalized sum:'
  write ( *, '(a)' ) ' '

  do i = 1, basis_num
    if ( s_total /= 0.0D+00 ) then
      s_partial = s_partial + sval(i) / s_total
      write ( *, '(i8,g14.6,2x,f10.6,2x,f10.6)' ) &
        i, sval(i), sval(i) / s_total, s_partial
    else
      write ( *, '(i8,g14.6,2x,f10.6,2x,f10.6)' ) &
        i, sval(i), 0.0D+00, 0.0D+00
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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
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
  character ( len = * ) string
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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine triangle_unit_set ( rule, norder, xtab, ytab, weight )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SET sets weights and abscissas for quadrature in a unit triangle.
!
!  Integration region:
!
!    Points (X,Y) such that
!
!      0 <= X,
!      0 <= Y, and
!      X + Y <= 1.
!
!  Graph:
!
!      ^
!    1 | *
!      | |\
!    Y | | \
!      | |  \
!    0 | *---*
!      +------->
!        0 X 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2000
!
!  Author:
!
!    John Burkardt
!
!  References:
!
!    H R Schwarz,
!    Methode der Finiten Elemente,
!    Teubner Studienbuecher, 1980.
!
!    Strang and Fix,
!    An Analysis of the Finite Element Method,
!    Prentice Hall, 1973, page 184.
!
!    Arthur H Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971.
!
!    O C Zienkiewicz,
!    The Finite Element Method,
!    McGraw Hill, Third Edition, 1977, page 201.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!     1, NORDER =  1, precision 1, Zienkiewicz #1.
!     2, NORDER =  3, precision 2, Strang and Fix formula #1.
!     3, NORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
!     4, NORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
!     5, NORDER =  6, precision 3, Strang and Fix formula #4.
!     6, NORDER =  6, precision 3, Stroud formula T2:3-1.
!     7, NORDER =  6, precision 4, Strang and Fix formula #5.
!     8, NORDER =  7, precision 4, Strang and Fix formula #6.
!     9, NORDER =  7, precision 5, Strang and Fix formula #7,
!        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
!    10, NORDER =  9, precision 6, Strang and Fix formula #8.
!    11, NORDER = 12, precision 6, Strang and Fix formula #9.
!
!    Output, integer ( kind = 4 ) NORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), YTAB(NORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) norder
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) weight(*)
  real ( kind = 8 ) xtab(*)
  real ( kind = 8 ) ytab(*)
  real ( kind = 8 ) z
!
!  1 point, precision 1.
!
  if ( rule == 1 ) then

    a = 1.0D+00 / 3.0D+00
    w = 1.0D+00

    norder = 1
    xtab(1) = a
    ytab(1) = a
    weight(1) = w
!
!  3 points, precision 2, Strang and Fix formula #1.
!
  else if ( rule == 2 ) then

    a = 1.0D+00
    b = 3.0D+00
    c = 4.0D+00
    d = 6.0D+00

    norder = 3
    xtab(1:3) =   (/ c, a, a /) / d
    ytab(1:3) =   (/ a, c, a /) / d
    weight(1:3) = (/ a, a, a /) / b
!
!  3 points, precision 2, Strang and Fix formula #2.
!
  else if ( rule == 3 ) then

    a = 0.5D+00
    b = 1.0D+00
    c = 1.0D+00 / 3.0D+00
    z = 0.0D+00

    norder = 3
    xtab(1:3) =   (/ z, a, a /)
    ytab(1:3) =   (/ a, z, a /)
    weight(1:3) = (/ c, c, c /)
!
!  4 points, precision 3, Strang and Fix formula #3.
!
  else if ( rule == 4 ) then

    a =   6.0D+00
    b =  10.0D+00
    c =  18.0D+00
    d =  25.0D+00
    e = -27.0D+00
    f =  30.0D+00
    g =  48.0D+00

    norder = 4
    xtab(1:4) =   (/ b, c, a, a /) / f
    ytab(1:4) =   (/ b, a, c, a /) / f
    weight(1:4) = (/ e, d, d, d /) / g
!
!  6 points, precision 3, Strang and Fix formula #4.
!
  else if ( rule == 5 ) then

    a = 0.659027622374092D+00
    b = 0.231933368553031D+00
    c = 0.109039009072877D+00
    w = 1.0D+00 / 6.0D+00

    norder = 6
    xtab(1:6) =   (/ a, a, b, b, c, c /)
    ytab(1:6) =   (/ b, c, a, c, a, b /)
    weight(1:6) = (/ w, w, w, w, w, w /)
!
!  6 points, precision 3, Stroud T2:3-1.
!
  else if ( rule == 6 ) then

    a = 0.0D+00
    b = 0.5D+00
    c = 2.0D+00 /  3.0D+00
    d = 1.0D+00 /  6.0D+00
    v = 1.0D+00 / 30.0D+00
    w = 3.0D+00 / 10.0D+00

    norder = 6
    xtab(1:6) =   (/ a, b, b, c, d, d /)
    ytab(1:6) =   (/ b, a, b, d, c, d /)
    weight(1:6) = (/ v, v, v, w, w, w /)
!
!  6 points, precision 4, Strang and Fix, formula #5.
!
  else if ( rule == 7 ) then

    a = 0.816847572980459D+00
    b = 0.091576213509771D+00
    c = 0.108103018168070D+00
    d = 0.445948490915965D+00
    v = 0.109951743655322D+00
    w = 0.223381589678011D+00

    norder = 6
    xtab(1:6) =   (/ a, b, b, c, d, d /)
    ytab(1:6) =   (/ b, a, b, d, c, d /)
    weight(1:6) = (/ v, v, v, w, w, w /)
!
!  7 points, precision 4, Strang and Fix formula #6.
!
  else if ( rule == 8 ) then

    a = 1.0D+00 / 3.0D+00
    c = 0.736712498968435D+00
    d = 0.237932366472434D+00
    e = 0.025355134551932D+00
    v = 0.375D+00
    w = 0.104166666666667D+00

    norder = 7
    xtab(1:7) =   (/ a, c, c, d, d, e, e /)
    ytab(1:7) =   (/ a, d, e, c, e, c, d /)
    weight(1:7) = (/ v, w, w, w, w, w, w /)
!
!  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
!
  else if ( rule == 9 ) then

    a = 1.0D+00 / 3.0D+00
    b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    c = ( 6.0D+00 -           sqrt ( 15.0D+00 ) ) / 21.0D+00
    d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    e = ( 6.0D+00 +           sqrt ( 15.0D+00 ) ) / 21.0D+00
    u = 0.225D+00
    v = ( 155.0D+00 - sqrt ( 15.0D+00 ) ) / 1200.0D+00
    w = ( 155.0D+00 + sqrt ( 15.0D+00 ) ) / 1200.0D+00

    norder = 7
    xtab(1:7) =   (/ a, b, c, c, d, e, e /)
    ytab(1:7) =   (/ a, c, b, c, e, d, e /)
    weight(1:7) = (/ u, v, v, v, w, w, w /)
!
!  9 points, precision 6, Strang and Fix formula #8.
!
  else if ( rule == 10 ) then

    a = 0.124949503233232D+00
    b = 0.437525248383384D+00
    c = 0.797112651860071D+00
    d = 0.165409927389841D+00
    e = 0.037477420750088D+00

    u = 0.205950504760887D+00
    v = 0.063691414286223D+00

    norder = 9
    xtab(1:9) =   (/ a, b, b, c, c, d, d, e, e /)
    ytab(1:9) =   (/ b, a, b, d, e, c, e, c, d /)
    weight(1:9) = (/ u, u, u, v, v, v, v, v, v /)
!
!  12 points, precision 6, Strang and Fix, formula #9.
!
  else if ( rule == 11 ) then

    a = 0.873821971016996D+00
    b = 0.063089014491502D+00
    c = 0.501426509658179D+00
    d = 0.249286745170910D+00
    e = 0.636502499121399D+00
    f = 0.310352451033785D+00
    g = 0.053145049844816D+00

    u = 0.050844906370207D+00
    v = 0.116786275726379D+00
    w = 0.082851075618374D+00

    norder = 12
    xtab(1:12) =   (/ a, b, b, d, c, d, e, e, f, f, g, g /)
    ytab(1:12) =   (/ b, a, b, c, d, d, f, g, e, g, e, f /)
    weight(1:12) = (/ u, u, u, v, v, v, w, w, w, w, w, w /)
!
!  13 points, precision 7, Strang and Fix, formula #10.
!
  else if ( rule == 12 ) then

    a = 0.479308067841923D+00
    b = 0.260345966079038D+00
    c = 0.869739794195568D+00
    d = 0.065130102902216D+00
    e = 0.638444188569809D+00
    f = 0.312865496004875D+00
    g = 0.048690315425316D+00
    h = 1.0D+00 / 3.0D+00
    t = 0.175615257433204D+00
    u = 0.053347235608839D+00
    v = 0.077113760890257D+00
    w = -0.149570044467670D+00

    norder = 13
    xtab(1:13) =   (/ a, b, b, c, d, d, e, e, f, f, g, g, h /)
    ytab(1:13) =   (/ b, a, b, d, c, d, f, g, e, g, e, f, h /)
    weight(1:13) = (/ t, t, t, u, u, u, v, v, v, v, v, v, w /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_UNIT_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
    stop

  end if

  return
end
