program main

!*****************************************************************************80
!
!! MAIN is the main routine for the CVT_BASIS_FLOW program.
!
!  Discussion:
!
!    What we really want to do is take a thousand points and put
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
!    Each record, in turn, contains a fixed number of data values
!    that describe a particular gene expression experiment.
!
!    Each record will be regarded as a point in N dimensional space.
!
!    The program will try to cluster the data, that is, to organize
!    the data by defining a number of cluster centers, which are
!    also points in N dimensional space, and assigning each record
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
!    19 July 2005: Lee requests that the program be revised, so that
!    the input velocity data can be either in two column format, or
!    one column format.  The easiest way to do this is to allow
!    "any" column format...  I set "COMP_NUM" to the number of components
!    of the solution associated with a node.  I think this is working now.
!
!    26 June 2005: Lee has made a new version of INOUT data, called
!    INOUT_FLOW2.  The grid is 61 by 61 instead of 41 by 41.
!    There are 800 time steps instead of 500.  I think that the
!    treatment of ALPHA is the same.
!
!    I discovered that all the old data doesn't know about the use
!    of more than one set of solution files.  This is a headache I
!    will not try to treat today!  I'll just fix RUN_92 and proceed.
!
!    21 September 2004: Janet and I both forgot that the program now
!    allows more than one set of solution files.  I added more explanation
!    to the printout.
!
!    17 July 2004: Janet has two sets of files.  She'd like to be able
!    to process either set, or both sets, without having to rename.
!    I am trying a modification of the code that allows up to 10
!    such sets of data to be handled.  This new code takes slightly
!    different input (after the first UV file family, you must now
!    say "none" if that's all, or enter the name of the next file family.)
!    Therefore, none of the old datasets will work properly with this
!    code until they've had an extra "none" statement inserted!
!
!    14 July 2004: Janet has data from a new problem.  She needs to
!    preprocess the data herself, rather than having the code subtract
!    off a multiple of the steady state.  I modified the code so that
!    run_type 1 or 2 should work, along with a specification of "none"
!    for the name of the steady state file.
!
!    07 July 2004: Max and Lee want CPU timing data for 4, 8 and 16 vectors
!    for the T-Cell.  The numbers I came up with, on the Alpha,
!    timing only the "interesting" code, are:
!
!       4:  338 seconds
!       8:  891 seconds
!      16: 2386 seconds
!
!    15 July 2003: This code still has vestigial code in it to handle
!    the case where the cluster size is to vary from CLUSTER_LO to
!    CLUSTER_HI.  To accommodate this feature, the code does other
!    things that discommode other more useful routines.  In particular,
!    the CLUSTER vector, which assigns points to cluster centers, is
!    not returned from ANALYSIS_RAW, because you might have been looking
!    at a range of cluster sizes...(but of course I haven't done that
!    in a long time!).  So rip that feature out, and get right with
!    the code lord!  This is making your preconditioning work a little
!    iffy, because now you assign clusters in the preconditioned world,
!    throw those assigments away, unprecondition, and then try to
!    recover the assignments...
!
!    14 July 2003: H C Lee requests that Mass Matrix preconditioning
!    used in POD_BASIS be copied into this code as well.  This
!    requires one extra input, the name of the element file, or the
!    word NONE, which means I have to update ALL THE OLD INPUT FILES
!    to say "NONE".
!
!    22 June 2002: Lee says drop double cell problem, add INOUT
!    (same form as very first problem, but now more correct(??), and
!    CAVITY.
!
!    21 June 2002: Corrected allocation of CLUSTER from CLUSTER_HI
!    to POINT_NUM.
!
!    Added run type 9 for Lee's 600 solutions of "double" cell problem.
!
!    20 June 2002: Lee requested real ( kind = 8 ) calculations.
!    Strangely enough, printing out the generators with the D25.15
!    format seemed to cause some data to print as 0.xxxx-313 instead of
!    0.xxxxE-313.  I tried to fix this.
!
!    20 June 2002: Lee has 500 solutions for a new region, the T-Cell.
!    Preprocessing is to subtract 5/3 of steady from 1-250, and
!    1/3 of steady from 251 to 500.  Also, normalize the result.
!    This will be Run type 6.
!
!    19 June 2002: Max requests a run in which data, after preprocessing,
!    is normalized, that is, simply divided by its length.  OK.
!
!    19 June 2002: Lee requests no "#" comment records in output files.  OK.
!
!    24 April 2002: Lee says I should have subtracted STEADY.DAT, not
!    UP000.DAT.  OK, redo.
!
!    23 April 2002: Got new data from Lee, reran Run #3 and sent results.
!    There are no singleton clusters anymore, and the clusters actually
!    comprise time segments.
!
!    19 April 2002: Lee complains that some of the generators do not
!    have zero boundary conditions as they should.  I believe that if
!    all points have zero boundary conditions, then their generators
!    will as well.  So now I have to go back and figure out:
!    * are these generators meaningless or discountable for some reason?
!    * did he give me the wrong ranges or coefficients, so that the
!      preprocessing I did did not fix the BC?
!    * or is there something else going on.
!
!    ...I discovered that the bad generators are #4 and #6, and that
!    these both correspond to clusters of 1 point, and that these
!    are data points #100 and #1 respectively.  Highly nonrandom values.
!    I asked Lee to ponder this.
!
!    ...I discovered I was accidentally NOT preprocessing data point #1.
!    I'll rerun, but why #100 wasn't getting fixed I still don't know.
!
!    18 April 2002: I still can't locate the problem.  I found a version
!    of the code that ran on March 30.  I recompiled it and ran it,
!    and it runs fine.  I note that the problem seems to be occurring
!    in HMEANS, and in particular in the first step, identification
!    of the nearest cluster.  Let me see if I can fix that.
!
!    ...Looks like I forgot to allocate CLUSTER, which I just moved
!    into the main program for some good reason.  Let's see.
!
!    Another problem may be the treatment of null clusters.  I've
!    added the variable NULL_CLUSTER_POLICY, currently set to 0,
!    which means just ignore null clusters.
!
!    ...FINALLY, things got better.  The code was running OVERNIGHT
!    and not getting results.  Now it runs in 5 minutes.
!
!    I had to be extra careful about null clusters - stupid junk things!
!    Now, under policy 0, I do not adjust the centers of null clusters
!    after detecting them, and I just make sure I divide the coordinates
!    by max ( cluster_pop, 1 ) to avoid division by 0 problems.  Also,
!    I am starting to think about policy 1 (replace center by a data point)
!    and 2 (replace center by a random point in convex hull of data),
!    but after this SECOND miserable experience with trying to avoid
!    zero clusters, I say just run more damn initial random configurations
!    and lump it.
!
!
!    17 April 2002: Trying to reason out the normalized data code.
!    The KMEANS algorithm has to change a lot.
!
!
!    16 April 2002: For Lee's Run #3, the code is running forever.
!    I'm worried that the code to avoid empty clusters is the cause
!    of this.  I just killed a run that had taken more than 6 hours
!    of time without result yet.  I will turn off zero-cluster avoidance,
!    and increase from 15 to 30 random starts, and I bet it runs
!    in less than half an hour.  The zero cluster business may have
!    been exacerbated the simple normalization that Lee prescribed,
!    subtracting a multiple of the steady state, which probably brings
!    the solution data closer together.
!
!
!    15 April 2002: Modified code by inserting DATA_SIZE and DATA_D2_READ,
!    which allow blanks and # comments in data file.  This is because
!    I'm starting to have multiple data files, and I would like internal
!    comments.
!
!    22 March 2002: Returning to code after having worked on KMEANS
!    algorithms.  HMEANS is the generic name for the cluster/average
!    method and KMEANS for the approach which considers moving each
!    point to every cluster.  I'm sticking with the RAW approach for now.
!
!    06 March 2002:  With the current CVT algorithm, Mr Lee's data
!    finds a local minimum most of the time, with an energy of about
!    180, and a much better minimum occasionally, with an energy of 63.
!    I would estimate 90% of the solutions are the local minimum.
!    Whether or not this is wise, I have spent a lot of time trying
!    to think of improvements to the algorithm.  I have hunted up
!    two Applied Statistics algorithms, and looked at HMEANS/KMEANS
!    from the Melendez book, and tried to think of better initial
!    guesses for the cluster centers.  Progress has been uncertain and
!    slow.
!
!    28 February 2002: Max says the CVT iteration should converge.
!    Rather than computing 50 random initial configurations to
!    30 iterations, try running one configuration to convergence.
!    Interesting things:
!    *  about how many steps does this take?
!    *  Are we sure that pretty much any initial configuration will drive
!       down to the same energy?
!    Also, send Hyung-Chun Lee the first set of generators, so he can get
!    an idea of what's coming.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: v0_file_max = 10

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) a_dot_p
  integer bandwidth
  integer, allocatable, dimension ( : ) :: cluster
  real ( kind = 8 ), allocatable, dimension (:,:) :: cluster_center
  integer, parameter :: cluster_max = 10
  integer cluster_hi
  integer cluster_it_max
  integer cluster_lo
  integer cluster_num
  integer cluster_range
  real ( kind = 8 ), allocatable, dimension ( : ) :: cluster_size
  real ( kind = 8 ), allocatable, dimension ( : ) :: cluster_size_inv
  logical comment
  character comment_char
  integer comp_num
  integer dim_num
  real ( kind = 8 ) r8vec_norm2
  character ( len = 100 ) element_file_name
  integer element_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: energy
  integer energy_it_max
  logical file_exist
  character ( len = 100 ) gen_file_name
  integer gen_unit
  integer i
  integer ierror
  integer info
  integer j
  integer, allocatable, dimension ( :, : ) :: node
  integer node_num
  real ( kind = 8 ) norm
  integer normal
  integer npe
  integer nrhs
  integer, parameter :: null_cluster_policy = 1
  real ( kind = 8 ), allocatable, dimension ( :, :) :: point
  integer point_num
  integer point_num2
  integer run_type
  logical s_eqi
  character ( len = 11 ) s_of_i4
  integer seed
  character ( len = 100 ) steady_file_name
  real ( kind = 8 ) steady_max
  real ( kind = 8 ) steady_norm
  integer swap_num
  real t1
  real t2
  integer thin
  real time_cvt
  real time_mass
  real ( kind = 8 ), allocatable, dimension ( : ) :: v
  real ( kind = 8 ), allocatable, dimension ( : ) :: v_steady
  character ( len = 100 ) v_file_name
  integer v_file_num
  character ( len = 100 ) v0_file_name(v0_file_max)
  integer v0_file_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  real ( kind = 8 ), allocatable, dimension ( :,: ) :: xy
  integer xy_dim
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  character ( len = 100 ) xy_file_name
  integer xy_values
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min

  call timestamp ( )

  time_cvt = 0.0D+00
  time_mass = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_BASIS_FLOW'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Arrange a set of PDE solution data into clusters.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This NEW version of the code handles up to'
  write ( *, '(i6,a)' ) v0_file_max, ' families of data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Null cluster policy:'
  write ( *, '(a)' ) '  0, do nothing, accept null clusters;'
  write ( *, '(a)' ) '  1, reset center to a random data point;'
  write ( *, '(a)' ) '  2, reset center to random point in hull;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NULL_CLUSTER_POLICY = ', null_cluster_policy

  seed = 123456789

  call random_initialize ( seed )
!
!  Get the run type.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The variable RUN_TYPE determines preprocessing:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* 1, NO steady state file, no preprocessing;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* 2, NO steady state file, no preprocessing;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* 3,    steady state file;'
  write ( *, '(a)' ) '    subtract 1/3 SS from solution  1'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 2 through 201'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 202 through 401.'
  write ( *, '(a)' ) '    do NOT normalize the data.'
  write ( *, '(a)' ) '    use all the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* 4,    steady state file;'
  write ( *, '(a)' ) '    subtract 1/3 SS from solution  1'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 2 through 201'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 202 through 401.'
  write ( *, '(a)' ) '    do NOT normalize the data;'
  write ( *, '(a)' ) '    then discard the EVEN numbered data files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* 5,    steady state file;'
  write ( *, '(a)' ) '    subtract 1/3 SS from solution  1'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 2 through 201'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 202 through 401.'
  write ( *, '(a)' ) '    then NORMALIZE the data;'
  write ( *, '(a)' ) '    use all the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* 6,    steady state file;'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 1 through 250'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 251 through 500.'
  write ( *, '(a)' ) '    do NOT normalize the data.'
  write ( *, '(a)' ) '    use all the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* 7,    steady state file;'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 1 through 250'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 251 through 500.'
  write ( *, '(a)' ) '    then NORMALIZE the data;'
  write ( *, '(a)' ) '    use all the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* 8,    steady state file;'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 1 through 250'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 251 through 500.'
  write ( *, '(a)' ) '    do NOT normalize the data.'
  write ( *, '(a)' ) '    then discard the ODD numbered data files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '* 9,    steady state file;'
  write ( *, '(a)' ) '    subtract 5/3 SS from solutions 1 through 400'
  write ( *, '(a)' ) '    subtract 1/3 SS from solutions 401 through 800.'
  write ( *, '(a)' ) '    do NOT normalize the data.'
  write ( *, '(a)' ) '    use all the data.'

  call i4_input ( 'What is the run type (1 through 9)?', run_type, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the run type.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  RUN_TYPE = ', run_type
!
!  Get the XY data file.
!
  call s_input ( '  What is the XY data file name (or NONE)?', &
    xy_file_name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the XY data file name.'
    stop
  end if

  if ( s_eqi ( xy_file_name, 'NONE' ) ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  No XY file.'

  else

    call dtable_header_read ( xy_file_name, xy_dim, node_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The file "' // trim ( xy_file_name ) &
      // '" contains the coordinates of ' &
      // trim ( s_of_i4 ( node_num ) ) // ' with spatial dimension ' &
      // trim ( s_of_i4 ( xy_dim ) ) // '.'
!
!  Allocate space for some arrays.
!
    allocate ( xy(1:xy_dim,1:node_num) )
!
!  Read in the spatial coordinates.
!
    call dtable_data_read ( xy_file_name, xy_dim, node_num, xy )
!
!  Extract some interesting data.
!
    x_min = minval ( xy(1,1:node_num) )
    x_max = maxval ( xy(1,1:node_num) )
    y_min = minval ( xy(2,1:node_num) )
    y_max = maxval ( xy(2,1:node_num) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  X minimum : ', x_min
    write ( *, '(a,g14.6)' ) '  X maximum : ', x_max
    write ( *, '(a,g14.6)' ) '  Y minimum : ', y_min
    write ( *, '(a,g14.6)' ) '  Y maximum : ', y_max

  end if

  comp_num = 0
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
    call dtable_header_read ( steady_file_name, comp_num, node_num )

    allocate ( v(comp_num*node_num) )
    allocate ( v_steady(comp_num*node_num) )

    call dtable_data_read ( steady_file_name, comp_num, node_num, v_steady )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Steady state information was read from'
    write ( *, '(a)' ) '  the file "' // trim ( steady_file_name ) // '".'

  end if

  v_file_num = 0
  v0_file_num = 0
  i = 0
!
!  Get the V0 file name.
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

  do while ( v0_file_num < v0_file_max )

    i = i + 1

    if ( i == 1 ) then
      call s_input ( 'What is the first solution file (in the first series)?', &
        v0_file_name(i), ierror )
    else
      call s_input ( &
        'What is the first solution file (in the NEXT series) or "NONE"?', &
        v0_file_name(i), ierror )
    end if

    if ( s_eqi ( v0_file_name(i), 'NONE' ) ) then
      exit
    end if

    v0_file_num = v0_file_num + 1
!
!  Presumably, all the solution files have the same name as the first
!  solution file, but with a numerical increment.  To begin with, simply count
!  the number of files.
!
    v_file_name = v0_file_name(i)

    do

      if ( .not. file_exist ( v_file_name ) ) then
        exit
      end if

      v_file_num = v_file_num + 1

      call file_name_inc ( v_file_name )

    end do

  end do

  if ( v_file_num == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  There do not seem to be any solution files.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) &
    '  The number of initial solution files is ', v0_file_num
  write ( *, '(a,i6)' ) &
    '  The total number of solution files is ', v_file_num
!
!  If there was no steady state file, then you need to read the
!  first data file in order to determine the number of components.
!
  if ( comp_num == 0 ) then

    call dtable_header_read ( v0_file_name(1), comp_num, node_num )

    allocate ( v(comp_num*node_num) )
    allocate ( v_steady(comp_num*node_num) )

  end if
!
!  Now we have enough information to set up a data structure.
!
!  Determine the spatial dimension (columns) and number of points (rows).
!
  dim_num = comp_num * node_num

  point_num = v_file_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data is stored in an M by N matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The "spatial" dimension M is   ', dim_num
  write ( *, '(a,i6)' ) '  The number of data points N is ', point_num
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

  do i = 1, v0_file_num

    v_file_name = v0_file_name(i)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Processing files starting with "' &
      // trim ( v0_file_name(i) ) // '".'

    do

      if ( .not. file_exist ( v_file_name ) ) then
        exit
      end if

      j = j + 1

      if ( .false. ) then
        write ( *, '(2x,i6,2x,a)' ) j, trim ( v_file_name )
      end if

      call dtable_data_read ( v_file_name, comp_num, node_num, v )
!
!  For RUN_TYPE = 3 (or 4 or 5), H C Lee wants to
!  subtract 1/3 of SS from solution  1.
!  subtract 5/3 of SS from solutions 2-201,
!  subtract 1/3 of SS from solutions 202-401.
!
      if ( run_type == 3 .or. run_type == 4 .or. run_type == 5 ) then

        if ( 1 <= j .and. j <= 1 ) then
          v(1:dim_num) = v(1:dim_num) -  v_steady(1:dim_num) / 3.0D+00
        else if ( 2 <= j .and. j <= 201 ) then
          v(1:dim_num) = v(1:dim_num) &
            - 5.0D+00 * v_steady(1:dim_num) / 3.0D+00
        else if ( 202 <= j .and. j <= 401 ) then
          v(1:dim_num) = v(1:dim_num) - v_steady(1:dim_num) / 3.0D+00
        end if
!
!  For RUN_TYPE = 6 or 7 or 8, H C Lee wants to
!  subtract 5/3 of SS from solutions 1-250,
!  subtract 1/3 of SS from solutions 251-500.
!
      else if ( run_type == 6 .or. run_type == 7 .or. run_type == 8 ) then

        if ( 1 <= j .and. j <= 250 ) then
          v(1:dim_num) = v(1:dim_num) &
            - 5.0D+00 * v_steady(1:dim_num) / 3.0D+00
        else if ( 251 <= j .and. j <= 500 ) then
          v(1:dim_num) = v(1:dim_num) - v_steady(1:dim_num) / 3.0D+00
        end if
!
!  For RUN_TYPE = 9, H C Lee wants to
!  subtract 5/3 of SS from solutions 1-400,
!  subtract 1/3 of SS from solutions 401-800.
!
      else if ( run_type == 9 ) then

        if ( 1 <= j .and. j <= 400 ) then
          v(1:dim_num) = v(1:dim_num) &
            - 5.0D+00 * v_steady(1:dim_num) / 3.0D+00
        else if ( 401 <= j .and. j <= 800 ) then
          v(1:dim_num) = v(1:dim_num) - v_steady(1:dim_num) / 3.0D+00
        end if

      end if

      point(1:dim_num,j) = v(1:dim_num)
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

      call file_name_inc ( v_file_name )

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
    write ( *, '(a,i6)' ) '  Thinning increment is       ', thin
    write ( *, '(a,i6)' ) '  Original input data size is ', point_num
    write ( *, '(a,i6)' ) '  Thinned data size is        ', point_num2

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
    write ( *, '(a,i6)' ) '  Thinning increment is       ', thin
    write ( *, '(a,i6)' ) '  Original input data size is ', point_num
    write ( *, '(a,i6)' ) '  Thinned data size is        ', point_num2

    point_num = point_num2

  end if
!
!  Get the range of cluster sizes to check.
!
  call i4_range_input ( 'Enter lower and upper number of clusters', &
    cluster_lo, cluster_hi, ierror )

  write ( *, * ) 'DEBUG: CLUSTER_LO = ', cluster_lo, &
    ' cluster_hi = ', cluster_hi

  if ( point_num < cluster_hi ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS_FLOW - Warning!'
    write ( *, '(a)' ) '  CLUSTER_HI exceeds number of points.'
    write ( *, '(a,i12)' ) '  CLUSTER_HI = ', cluster_hi
    write ( *, '(a,i12)' ) '  POINT_NUM =  ', point_num
    cluster_hi = point_num
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program will try to determine the minimum energy'
  write ( *, '(a)' ) '  of a clustering, for cluster sizes in the range:'
  write ( *, '(2x,2i6)' ) cluster_lo, cluster_hi

  allocate ( cluster(1:point_num) )
  allocate ( cluster_center(1:dim_num,1:cluster_hi) )
  allocate ( energy(1:cluster_hi) )
!
!  Get the number of different random starting cluster configurations.
!
  call i4_input ( &
    'Enter the number of different random cluster configurations to check', &
    cluster_it_max, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each number of clusters, the number of'
  write ( *, '(a)' ) '  distinct initial random configurations to be checked'
  write ( *, '(a,i6)' ) '  will be  ', cluster_it_max
!
!  Get the number of energy iterations for a particular random start:
!
  call i4_input ( 'Enter the number of energy iterations', &
    energy_it_max, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each initial random configuration, the number of'
  write ( *, '(a)' ) '  times the program will recompute the cluster centers,'
  write ( *, '(a,i6)' ) '  cluster components, and energy is ', energy_it_max
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
  call s_input ( &
    'Enter element file for mass matrix preconditioning or "none".', &
    element_file_name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS_FLOW - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the element file name.'
    stop
  end if

  if ( s_eqi ( element_file_name, 'none' ) ) then

  else

    call itable_header_read ( element_file_name, npe, element_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of elements = ', element_num
    write ( *, '(a,i6)' ) '  Number of nodes per element = ', npe

    allocate ( node(1:npe,1:element_num) )

    call itable_data_read ( element_file_name, npe, element_num, node )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_BASIS_FLOW - Fatal error!'
      write ( *, '(a)' ) '  Error reading the element file.'
      stop
    end if

    call bandwidth_determine ( npe, element_num, node, bandwidth )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  The bandwidth of the matrix is ', bandwidth
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

    call mass_matrix ( node_num, npe, element_num, node, bandwidth, xy, a )

    call cpu_time ( t2 )
    time_mass = time_mass + t2 - t1

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
      write ( *, '(a)' ) 'CVT_BASIS_FLOW - Fatal error!'
      write ( *, '(a)' ) '  DPBTRF reports the mass matrix is singular!'
      write ( *, '(a,i6)' ) '  The value of INFO is ', info
      stop
    end if
!
!  Premultiply POINT by L'.
!  Each column of POINT contains U and V coefficients.
!  We need to delicately extract U, compute L'*U, and put it back;
!  then repeat for V.
!
    allocate ( v1(1:node_num) )

    do j = 1, point_num

      do i = 1, comp_num

        v1(1:node_num) = point(i:comp_num*(node_num-1)+i:comp_num,j)

        call cpu_time ( t1 )

        call dtbmv ( 'L', 'T', 'N', node_num, bandwidth, a, bandwidth+1, v1, 1 )

        call cpu_time ( t2 )
        time_mass = time_mass + t2 - t1

        point(i:comp_num*(node_num-1)+i:comp_num,j) = v1(1:node_num)

      end do

    end do

    deallocate ( v1 )

  end if
!
!  Get the normalization option:
!
  call i4_input ( 'Enter 0 to use raw data, 1 to use normalized data.', &
    normal, ierror )

  if ( normal == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NORMAL = 0'
    write ( *, '(a)' ) '  Data will NOT be normalized.'
  else if ( normal == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NORMAL = 1'
    write ( *, '(a)' ) '  Data WILL be normalized.'
  end if
!
!  EITHER:
!
!  A) Compute the CVT of the raw data
!     (Notice that I do not return the CLUSTER assignments!
!     That's why I have to recompute them lower down.)
!
!  OR
!
!  B) Compute the CVT of normalized data
!
!
  if ( normal == 0 ) then

    call cpu_time ( t1 )

    call analysis_raw ( dim_num, point_num, point, cluster_lo, &
      cluster_hi, cluster_it_max, energy_it_max, energy, cluster_center, &
      null_cluster_policy, seed )

    call cpu_time ( t2 )
    time_cvt = time_cvt + t2 - t1
!
!  Preprocess the data.
!  In particular, we plan to consider linear combinations of solutions,
!  with the steady state solution being our fundamental solution.
!
!
!  THOUGHTS:
!
!  Now that I think about it, shouldn't column 1 be completely 0?
!  Does the steady state solution have any special status?
!  It does if I'm subtracting it from all the others.
!  It may deserve that status, particularly from physical reasons.
!  Then I'm doing a Voronoi tessellation of some affine space.
!  But in that case, I don't want to keep the steady state solution
!  in column 1, but rather, all 0's.
!
!  So I think I've convinced myself that what I am computing are
!  generators for the affine space of perturbations from the steady state.
!
!  I save a copy of V0 separately.
!  I include 0 as a data point in column 1.
!
!  In fact, now I think I need to constrain the code to always include
!  0 as a generator.
!
!
  else if ( normal == 1 ) then
!
!  Normalize the steady state solution (column 1).
!
    call r8vec_unit_euclidean ( dim_num, point(1:dim_num,1) )
!
!  Subtract the projection of the normalized steady state solution
!  from the later solutions, and normalize them.
!
    do j = 2, point_num

      a_dot_p = dot_product ( point(1:dim_num,j), point(1:dim_num,1) )

      point(1:dim_num,j) = point(1:dim_num,j) - a_dot_p * point(1:dim_num,1)

      call r8vec_unit_euclidean ( dim_num, point(1:dim_num,j) )

    end do
!
!  Now zap column 1.
!
    point(1:dim_num,1) = 0.0D+00
!
!  Report any zero columns.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Zero columns:'
    write ( *, '(a)' ) '(Expecting only column 1)'
    write ( *, '(a)' ) ' '
    do j = 1, point_num
      if ( r8vec_norm2 ( dim_num, point(1:dim_num,j) ) == 0.0D+00 ) then
        write ( *, '(i6)' ) j
      end if
    end do
!
!  Analyze the data.
!
    call cpu_time ( t1 )

    call analysis_normal ( dim_num, point_num, point, cluster_lo, &
      cluster_hi, cluster_it_max, energy_it_max, energy, cluster_center, &
      null_cluster_policy )

    call cpu_time ( t2 )
    time_cvt = time_cvt + t2 - t1
!
!  Bad value of NORMAL.
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS_FLOW - Fatal error!'
    write ( *, '(a,i6)' ) '  Bad value of NORMAL = ', normal
    stop

  end if
!
!  If mass matrix preconditioning was used, we have to
!  "unprecondition" now.  Each left singular vector is decomposable
!  into a "U" and "V" portion.  Both the "U" and "V" portion were
!  multiplied by the Cholesky factor L' before the SVD was applied.
!  Therefore, to get back to the original data world, we now have
!  to reverse that transformation.
!
  if ( s_eqi ( element_file_name, 'None' ) ) then

  else
!
!  Premultiply CLUSTER_CENTER by inverse ( L' ).
!  Each column of CLUSTER_CENTER contains U and V coefficients.
!  We need to delicately extract U, solve L'*U(new)=U(old), and put it back;
!  then repeat for V.
!
    nrhs = 1
    allocate ( v1(1:node_num) )


    call cpu_time ( t1 )

    do j = 1, cluster_hi

      do i = 1, comp_num

        v1(1:node_num) = cluster_center(i:comp_num*(node_num-1)+i:comp_num,j)

        call dtbtrs ( 'L', 'T', 'N', node_num, bandwidth, nrhs, a, &
          bandwidth+1, v1, node_num, info )

        cluster_center(i:comp_num*(node_num-1)+i:comp_num,j) = v1(1:node_num)

      end do

    end do
!
!  Premultiply POINT by inverse ( L' ).
!  Each column of CLUSTER_CENTER contains U and V coefficients.
!  We need to delicately extract U, solve L'*U(new)=U(old), and put it back;
!  then repeat for V.
!
    do j = 1, point_num

      do i = 1, comp_num

        v1(1:node_num) = point(i:comp_num*(node_num-1)+i:comp_num,j)

        call dtbtrs ( 'L', 'T', 'N', node_num, bandwidth, nrhs, a, &
          bandwidth+1, v1, node_num, info )

        point(i:comp_num*(node_num-1)+i:comp_num,j) = v1(1:node_num)

      end do

    end do

    call cpu_time ( t2 )
    time_mass = time_mass + t2 - t1

    deallocate ( a )
    deallocate ( v1 )

  end if
!
!  Multiple cluster values, presumably we want an energy plot.
!
  if ( cluster_lo < cluster_hi ) then

    allocate ( cluster_size(1:cluster_hi) )
    allocate ( cluster_size_inv(1:cluster_hi) )

    do i = 1, cluster_hi
      cluster_size(i) = real ( i, kind = 8 )
    end do

    cluster_size_inv(1:cluster_hi) = 1.0D+00 / cluster_size(1:cluster_hi)
      cluster_range = cluster_hi + 1 - cluster_lo

    call data_to_gnuplot ( 'raw.txt', cluster_range, &
      cluster_size(cluster_lo:cluster_hi), energy(cluster_lo:cluster_hi) )

    call data_to_gnuplot ( 'raw2.txt', cluster_range, &
      cluster_size_inv(cluster_lo:cluster_hi), energy(cluster_lo:cluster_hi) )

    deallocate ( cluster_size )
    deallocate ( cluster_size_inv )
!
!  If there was only one cluster value, presumably we want to
!  * print out the population of each cluster;
!  * write out the generators.
!
  else
!
!  Why am I recomputing the cluster assignments here?
!
    cluster_num = cluster_hi

    if ( normal == 0 ) then

      cluster(1:point_num) = 1

      call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
        cluster_center, point, cluster, swap_num )

    else

!     call nearest_cluster_normal ( dim_num, point_num, cluster_num, &
!       cluster_center, point, cluster, swap_num )

    end if

    call cluster_census ( dim_num, point_num, cluster_num, cluster_center, &
      point, cluster )

    if ( .false. ) then
      call cluster_list ( point_num, cluster )
    end if
!
!  Write vectors to files.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS_FLOW:'
    write ( *, '(a)' ) '  Ready to write the cluster generators to files.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Do you want comments in the header of the file?'
    write ( *, '(a)' ) '  (These begin with the "#" character.) (Y/N)'

    call s_input ( '  Enter Y or N:', comment_char, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_BASIS_FLOW - Fatal error!'
      write ( *, '(a)' ) '  Input error reading the comment option.'
      write ( *, '(a)' ) '  We will assume comments are acceptable.'
      comment_char = 'Y'
    end if

    if ( comment_char == 'Y' .or. comment_char == 'y' ) then
      comment = .true.
      write ( *, * ) 'The output files will include a commented header.'
    else
      comment = .false.
      write ( *, * ) 'The output files will NOT include a commented header.'
    end if

    gen_file_name = 'gen_000.txt'
    call get_unit ( gen_unit )

    do j = 1, cluster_hi

      call file_name_inc ( gen_file_name )

      if ( j == 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Writing first file ' // trim ( gen_file_name )
      end if

      if ( j == cluster_hi ) then
        write ( *, '(a)' ) '  Writing last file  ' // trim ( gen_file_name )
      end if

      v(1:comp_num*node_num) = cluster_center(1:comp_num*node_num,j)

      call dtable_write ( gen_file_name, comp_num, node_num, v, comment )

    end do

  end if

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( energy )
  deallocate ( point )
  deallocate ( v )
  if ( allocated ( v_steady ) ) then
    deallocate ( v_steady )
  end if
  if ( allocated ( xy ) ) then
    deallocate ( xy )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CPU_TIME (seconds):'
  write ( *, '(a,g14.6)' ) '  Mass matrix: ', time_mass
  write ( *, '(a,g14.6)' ) '  CVT:         ', time_cvt
  write ( *, '(a,g14.6)' ) '  Total:       ', time_mass + time_cvt
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_BASIS_FLOW'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine analysis_normal ( dim_num, point_num, point, cluster_lo, &
  cluster_hi, cluster_it_max, it_max, energy, cluster_center, &
  null_cluster_policy )

!*****************************************************************************80
!
!! ANALYSIS_NORMAL computes the energy for a range of number of clusters.
!
!  Discussion:
!
!    This version of the analysis routine is for "normalized" data.  That is,
!    all the solution vectors are replaced by their normalized differences
!    with the steady state solution.  Thus, our space is essentially
!    a sphere, and the distance between two points is now the angle
!    between them (in the plane defined by the two points and the
!    steady state solution).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORMAL:
!    0, analyze the raw data.
!    1, analyze the normalized data.
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, integer CLUSTER_LO, CLUSTER_HI, the low and high
!    cluster numbers that define a range of clusters to check.
!
!    Input, integer CLUSTER_IT_MAX, the number of different random
!    startup configurations to try.
!
!    Input, integer ENERGY_IT_MAX, the maximum number of energy iterations.
!
!    Output, real ( kind = 8 ) ENERGY(CLUSTER_HI), contains in entries
!    CLUSTER_LO through CLUSTER_HI the estimated minimum energy over
!    all clusterings of this size.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_HI), contains the
!    generators for the minimum cluster energy calculation, but only
!    for the last cluster size considered, namely, CLUSTER_HI.
!
!    Input, integer NULL_CLUSTER_POLICY, specifies what to do if a
!    null cluster is encountered.
!    0, do nothing.
!    1, reset center of null cluster to a random data point;
!    2, reset center of null cluster to a random point in the data hull;
!
  implicit none

  integer cluster_hi
  integer dim_num
  integer point_num

  integer, dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension (dim_num,cluster_hi) :: cluster_center
  real ( kind = 8 ), dimension (dim_num,cluster_hi) :: cluster_center_test
  integer cluster_it
  integer cluster_it_max
  integer cluster_lo
  integer cluster_num
  real ( kind = 8 ), dimension (cluster_hi) :: energy
  real ( kind = 8 ) energy_max(cluster_hi)
  real ( kind = 8 ) energy_min(cluster_hi)
  real ( kind = 8 ) energy_test
  integer i
  integer it
  integer it_max
  integer null_cluster_policy
  real ( kind = 8 ), dimension (dim_num,point_num) :: point
  real ( kind = 8 ), dimension (dim_num) :: r_min
  real ( kind = 8 ), dimension (dim_num) :: r_max
  real ( kind = 8 ) t

  if ( .false. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ANALYSIS_NORMAL:'
    write ( *, '(a)' ) '  Point data:'
    write ( *, '(a)' ) ' '
    call point_print ( dim_num, point_num, point )
  end if
!
!  Compute the minimum and maximum component values.
!
  do i = 1, dim_num
    r_min(i) = minval ( point(i,1:point_num) )
    r_max(i) = maxval ( point(i,1:point_num) )
  end do

  if ( dim_num <= 20 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Minimum and maximum data values:'
    write ( *, '(a)' ) ' '
    do i = 1, dim_num
      write ( *, '(i6,2f10.4)' ) i, r_min(i), r_max(i)
    end do
  end if
!
!  Consider a range of clusters.
!
  do cluster_num = cluster_lo, cluster_hi
!
!  For each cluster size, try several random starting configurations.
!
    energy_min(cluster_num) = huge ( energy_min )
    energy_max(cluster_num) = 0.0D+00

    do cluster_it = 1, cluster_it_max

      call hmeans_normal ( dim_num, point_num, cluster_num, it_max, it, &
        point, cluster, cluster_center_test, energy_test, null_cluster_policy )

      call kmeans_normal ( dim_num, point_num, cluster_num, it_max, it, &
        point, cluster, cluster_center_test, energy_test, null_cluster_policy )

      if ( energy_test < energy_min(cluster_num) ) then
        energy_min(cluster_num) = energy_test
        cluster_center(1:dim_num,1:cluster_hi) = &
          cluster_center_test(1:dim_num,1:cluster_hi)
      end if

      energy_max(cluster_num) = max ( energy_max(cluster_num), energy_test )

    end do

    energy(cluster_num) = energy_min(cluster_num)

  end do
!
!  Report energy ranges for the various cluster sizes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ANALYSIS_NORMAL:'
  write ( *, '(a)' ) '  Computed energy range for given cluster size:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (The minimum and maximum should be close if'
  write ( *, '(a)' ) '  we''re taking enough iterations.)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cluster  Minimum    Maximum'
  write ( *, '(a)' ) '  Size     Energy     Energy'
  write ( *, '(a)' ) ' '

  do cluster_num = cluster_lo, cluster_hi
    write ( *, '(i9,2f12.4)' ) &
      cluster_num, energy_min(cluster_num), energy_max(cluster_num)
  end do
!
!  Report best energy for the various cluster sizes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Energy table:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Cluster              Energy'
  write ( *, '(a)' ) 'Size      Energy     /point    Sqrt(E/Pt)'
  write ( *, '(a)' ) ' '

  do cluster_num = cluster_lo, cluster_hi
    t = energy(cluster_num) / dble ( point_num )
    write ( *, '(i9,3f12.4)' ) cluster_num, energy(cluster_num), t, sqrt ( t )
  end do

  return
end
subroutine analysis_raw ( dim_num, point_num, point, cluster_lo, cluster_hi, &
  cluster_it_max, energy_it_max, energy, cluster_center, null_cluster_policy, &
  seed )

!*****************************************************************************80
!
!! ANALYSIS_RAW computes the energy for a range of number of clusters.
!
!  Discussion:
!
!    This version of the analysis routine is for "raw" data.  That is,
!    all the solution vectors are treated as points in Euclidean space
!    with the usual distance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORMAL:
!    0, analyze the raw data.
!    1, analyze the normalized data.
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, integer CLUSTER_LO, CLUSTER_HI, the low and high
!    cluster numbers that define a range of clusters to check.
!
!    Input, integer CLUSTER_IT_MAX, the number of different random
!    startup configurations to try.
!
!    Input, integer ENERGY_IT_MAX, the maximum number of energy iterations.
!
!    Output, real ( kind = 8 ) ENERGY(CLUSTER_HI), contains in entries
!    CLUSTER_LO through CLUSTER_HI the estimated minimum energy over
!    all clusterings of this size.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_HI), contains the
!    generators for the minimum cluster energy calculation, but only
!    for the last cluster size considered, namely, CLUSTER_HI.
!
!    Input, integer NULL_CLUSTER_POLICY, specifies what to do if a
!    null cluster is encountered.
!    0, do nothing.
!    1, reset center of null cluster to a random data point;
!    2, reset center of null cluster to a random point in the data hull;
!
!    Input/output, integer SEED, a seed for the random number generator.
!
  implicit none

  integer cluster_hi
  integer dim_num
  integer point_num

  integer, dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension (dim_num,cluster_hi) :: cluster_center
  real ( kind = 8 ), dimension (dim_num,cluster_hi) :: cluster_center_test
  integer cluster_it
  integer cluster_it_max
  integer cluster_lo
  integer cluster_num
  logical, parameter :: debug = .true.
  real ( kind = 8 ), dimension (cluster_hi) :: energy
  integer energy_it_max
  real ( kind = 8 ) energy_max(cluster_hi)
  real ( kind = 8 ) energy_min(cluster_hi)
  real ( kind = 8 ) energy_test
  integer i
  integer it
  integer null_cluster_policy
  real ( kind = 8 ), dimension (dim_num,point_num) :: point
  real ( kind = 8 ), dimension (dim_num) :: r_min
  real ( kind = 8 ), dimension (dim_num) :: r_max
  integer seed
  real ( kind = 8 ) t

  if ( .false. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ANALYSIS_RAW:'
    write ( *, '(a)' ) '  Point data:'
    write ( *, '(a)' ) ' '
    call point_print ( dim_num, point_num, point )
  end if
!
!  Compute the minimum and maximum component values.
!
  do i = 1, dim_num
    r_min(i) = minval ( point(i,1:point_num) )
    r_max(i) = maxval ( point(i,1:point_num) )
  end do

  if ( dim_num <= 20 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Minimum and maximum data values:'
    write ( *, '(a)' ) ' '
    do i = 1, dim_num
      write ( *, '(i6,2f10.4)' ) i, r_min(i), r_max(i)
    end do
  end if
!
!  Consider a range of clusters.
!
  do cluster_num = cluster_lo, cluster_hi

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Number of clusters allowed: ', cluster_num
!
!  For each cluster size, try several random starting configurations.
!
    energy_min(cluster_num) = huge ( energy_min )
    energy_max(cluster_num) = 0.0D+00

    do cluster_it = 1, cluster_it_max

      write ( *, '(a)' ) ' '
      write ( *, '(i6)' ) cluster_it

      call cluster_initialize_raw ( dim_num, point_num, cluster_num, point, &
        cluster, cluster_center_test, energy_test )

      it = 0
      write ( *, '(a,g14.6,i6)' ) 'Initial_RAW  ', energy_test, it

      call hmeans_raw ( dim_num, point_num, cluster_num, energy_it_max, &
        it, point, cluster, cluster_center_test, energy_test, &
        null_cluster_policy, seed )

      write ( *, '(a,g14.6,i6)' ) 'HMEANS_RAW   ', energy_test, it

      call kmeans_raw ( dim_num, point_num, cluster_num, energy_it_max, &
        it, point, cluster, cluster_center_test, energy_test, &
        null_cluster_policy, seed )

      write ( *, '(a,g14.6,i6)' ) 'KMEANS_RAW   ', energy_test, it

      if ( energy_test < energy_min(cluster_num) ) then
        energy_min(cluster_num) = energy_test
        cluster_center(1:dim_num,1:cluster_hi) = &
          cluster_center_test(1:dim_num,1:cluster_hi)
      end if

      energy_max(cluster_num) = max ( energy_max(cluster_num), energy_test )

    end do

    energy(cluster_num) = energy_min(cluster_num)

  end do
!
!  Report energy ranges for the various cluster sizes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ANALYSIS_RAW:'
  write ( *, '(a)' ) '  Computed energy range for given cluster size:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (The minimum and maximum should be close if'
  write ( *, '(a)' ) '  we''re taking enough iterations.)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number'
  write ( *, '(a)' ) '  of       Minimum      Maximum'
  write ( *, '(a)' ) '  Clusters Energy       Energy'
  write ( *, '(a)' ) ' '

  do cluster_num = cluster_lo, cluster_hi
    write ( *, '(i7,2f14.4)' ) &
      cluster_num, energy_min(cluster_num), energy_max(cluster_num)
  end do
!
!  Report best energy for the various cluster sizes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Energy table:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Number'
  write ( *, '(a)' ) 'of                   Energy'
  write ( *, '(a)' ) 'Clusters  Energy     /point      Sqrt(E/Pt)'
  write ( *, '(a)' ) ' '

  do cluster_num = cluster_lo, cluster_hi
    t = energy(cluster_num) / dble ( point_num )
    write ( *, '(i7,3f14.4)' ) cluster_num, energy(cluster_num), t, sqrt ( t )
  end do

  return
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
!    Input, integer NPE, the number of nodes per element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer NODE(NPE,ELEMENT_NUM), the nodes that make up each element.
!
!    Output, integer BANDWIDTH, the (half) bandwidth of the matrix.
!
  implicit none

  integer element_num
  integer npe

  integer bandwidth
  integer elem
  integer i1
  integer i2
  integer n1
  integer n2
  integer node(npe,element_num)

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
  integer itemp

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
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
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
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine cluster_census ( dim_num, point_num, cluster_num, cluster_center, &
  point, cluster )

!*****************************************************************************80
!
!! CLUSTER_CENSUS computes and prints the population of each cluster.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the cluster
!    generators.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, integer CLUSTER(POINT_NUM), the cluster to which each
!    point is assigned.
!
  implicit none

  integer cluster_num
  integer dim_num
  integer point_num

  integer cluster(point_num)
  real ( kind = 8 ), dimension (dim_num,cluster_num) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  integer, dimension ( cluster_num ) :: cluster_max
  integer, dimension ( cluster_num ) :: cluster_min
  integer, dimension ( cluster_num ) :: cluster_population
  integer i
  integer j
  integer percent1
  integer percent2
  real ( kind = 8 ), dimension (dim_num,point_num) :: point
  integer swap_num
  real ( kind = 8 ) total_energy

  call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
    cluster_center, point, cluster, swap_num )

  cluster_energy(1:cluster_num) = 0.0D+00
  cluster_population(1:cluster_num) = 0

  do i = 1, point_num

    j = cluster(i)

    cluster_population(j) = cluster_population(j) + 1

    cluster_energy(j) = cluster_energy(j) &
      + sum ( ( cluster_center(1:dim_num,j) - point(1:dim_num,i) )**2 )

  end do

  total_energy = sum ( cluster_energy(1:cluster_num) )

  cluster_min(1:cluster_num) = point_num + 1
  cluster_max(1:cluster_num) = 0

  do i = 1, point_num
    j = cluster(i)
    cluster_min(j) = min ( cluster_min(j), i )
    cluster_max(j) = max ( cluster_max(j), i )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLUSTER_CENSUS'
  write ( *, '(a)' ) '  Individual cluster population and energy'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Index    Population   Percentage   Energy  Percentage  Min  Max'
  write ( *, '(a)' ) ' '
  do i = 1, cluster_num
    percent1 = int ( dble ( 100 * cluster_population(i) ) / dble ( point_num ) )
    percent2 = int ( 100.0D+00 * cluster_energy(i) ) / total_energy
    write ( *, '(1x,i6,8x,i6,10x,i3,g14.6,4x,i3,2i5)' ) &
      i, cluster_population(i), &
      percent1, cluster_energy(i), percent2, cluster_min(i), cluster_max(i)
  end do

  percent1 = 100
  percent2 = 100

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               ------          ---  ------------    ---'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a5,8x,i6,10x,i3,g14.6,4x,i3,2i5)' ) 'Total', &
  sum ( cluster_population(1:cluster_num) ), &
    percent1, sum ( cluster_energy(1:cluster_num) ), percent2, 1, point_num

  return
end
subroutine cluster_initialize_raw ( dim_num, point_num, cluster_num, point, &
  cluster, cluster_center, energy )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_RAW initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each cluster center is a random convex combination
!    of the data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates of
!    the points.
!
!    Output, integer CLUSTER(POINT_NUM), the clusters to which points
!    are assigned.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
!    Output, real ( kind = 8 ) ENERGY, the energy of the clustering.
!
  implicit none

  integer cluster_num
  integer dim_num
  integer point_num

  integer cluster(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) column_sum
  real ( kind = 8 ) energy
  real ( kind = 8 ) factor(point_num,cluster_num)
  integer j
  real ( kind = 8 ) point(dim_num,point_num)
  integer swap_num
!
!  Get a PxC block of random factors.
!
  call random_number ( harvest = factor(1:point_num,1:cluster_num) )
!
!  Make each column of factors have unit sum.
!
  do j = 1, cluster_num
    column_sum = sum ( factor(1:point_num,j) )
    factor(1:point_num,j) = factor(1:point_num,j) / column_sum
  end do
!
!  Set centers = points * factors.
!
  cluster_center(1:dim_num,1:cluster_num) = &
    matmul ( point(1:dim_num,1:point_num), factor(1:point_num,1:cluster_num) )
!
!  Assign points to the nearest centers.
!
  call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
    cluster_center, point, cluster, swap_num )
!
!  Determine the energy of the clustering.
!
  call energy_raw ( dim_num, point_num, cluster_num, point, &
    cluster_center, cluster, energy )

  return
end
subroutine cluster_list ( point_num, cluster )

!*****************************************************************************80
!
!! CLUSTER_LIST prints out the assignments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER(POINT_NUM), the cluster to which each
!    point is assigned.
!
  implicit none

  integer point_num

  integer cluster(point_num)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I  Cluster(I)'
  write ( *, '(a)' ) ' '
  do i = 1, point_num
    write ( *, '(i4,8x,i4)' ) i, cluster(i)
  end do

  return
end
subroutine data_to_gnuplot ( file_name, n, x, y )

!*****************************************************************************80
!
!! DATA_TO_GNUPLOT writes data to a file suitable for processing by GNUPLOT.
!
!  Discussion:
!
!    Once the data file is written, the GNUPLOT program can be used
!    to display the data, using commands like:
!
!      set term post default
!      set grid
!      set yrnage [0:*]
!      set title "Number of Clusters vs Total Energy, Normalized Data"
!      plot 'file_name'
!      quit
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file
!    to which the data is written.
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the data.
!
  integer n

  character ( len = * ) file_name
  integer i
  integer iunit
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace' )

  do i = 1, n
    write ( iunit, '(g16.8,2x,g16.8)' ) x(i), y(i)
  end do

  close ( unit = iunit )

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
  integer digit

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
!    Input, integer DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
function distance_normal_sq ( n, x, y )

!*****************************************************************************80
!
!! DISTANCE_NORMAL_SQ computes the distance between normalized vectors.
!
!  Discussion:
!
!    Actually, it computes the SQUARE of the distance.
!
!    For "normalized" vectors, we assume that an appropriate multiple
!    of the steady state solution has been subtracted, and that what
!    remains can be regarded as perturbations from 0 whose direction
!    (but not magnitude) is of importance.
!
!    The distance between two such vectors, then, is the angle between
!    them, or, for our purposes, the absolute value of the sine of the
!    angle.  (This allows us to identify vectors V and -V, for instance).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the vectors whose distance is
!    to be computed.
!
!    Output, real ( kind = 8 ) DISTANCE_NORMAL_SQ, the square of the
!    "distance" between the vectors.
!
  integer n

  real ( kind = 8 ) cosine
  real ( kind = 8 ) distance_normal_sq
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_norm
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y_norm

  x_norm = sqrt ( dot_product ( x(1:n), x(1:n) ) )
  y_norm = sqrt ( dot_product ( y(1:n), y(1:n) ) )
  cosine = dot_product ( x(1:n), y(1:n) ) / ( x_norm * y_norm )

  distance_normal_sq = 1.0D+00 - cosine**2

  return
end
subroutine dtable_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! DTABLE_DATA_READ reads data from a double precision table file.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer m
  integer n

  integer ierror
  character ( len = * ) input_filename
  integer input_unit
  integer ios
  integer j
  character ( len = 255 ) line
  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTABLE_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i6)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i6)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i6)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine dtable_data_write ( output_unit, m, n, table )

!*****************************************************************************80
!
!! DTABLE_DATA_WRITE writes data to a double precision table file.
!
!  Discussion:
!
!    This routine writes a single line of output for each point,
!    containing its spatial coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the output unit.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer m
  integer n

  integer output_unit
  integer j
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Create the format string.
!
  write ( string, '(a1,i6,a1,i6,a1,i6,a1)' ) '(', m, 'g', 14, '.', 6, ')'
  call s_blank_delete ( string )

  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do

  return
end
subroutine dtable_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! DTABLE_HEADER_READ reads the header from a double precision table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer M, spatial dimension.
!
!    Output, integer N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer m
  integer n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine dtable_header_write ( output_filename, output_unit, m, n )

!*****************************************************************************80
!
!! DTABLE_HEADER_WRITE writes the header to a double precision table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output filename.
!
!    Input, integer OUTPUT_UNIT, the output unit.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
  implicit none

  integer m
  integer n
  character ( len = * ) output_filename
  integer output_unit
  character ( len = 40 ) string
  real ( kind = 8 ), parameter :: x = 1.0D+00

  call timestring ( string )

  write ( output_unit, '(a)'       ) '#  ' // trim ( output_filename )
  write ( output_unit, '(a)'       ) '#  created by TABLE_IO.F90'
  write ( output_unit, '(a)'       ) '#  at ' // trim ( string )
  write ( output_unit, '(a)'       ) '#'
  write ( output_unit, '(a,i3)'    ) '#  Spatial dimension M = ', m
  write ( output_unit, '(a,i6)'    ) '#  Number of points N = ', n
  write ( output_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) = ', &
    epsilon ( x )
  write ( output_unit, '(a)'       ) '#'

  return
end
subroutine dtable_write ( output_filename, m, n, table, comment )

!*****************************************************************************80
!
!! DTABLE_WRITE writes a double precision table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output filename.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
!    Input, logical COMMENT, is TRUE if the header comments are to be included.
!
  implicit none

  integer m
  integer n

  logical comment
  integer ios
  character ( len = * ) output_filename
  integer output_unit
  real ( kind = 8 ) table(m,n)

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  if ( comment ) then
    call dtable_header_write ( output_filename, output_unit, m, n )
  end if

  call dtable_data_write ( output_unit, m, n, table )

  close ( unit = output_unit )

  return
end
subroutine energy_normal ( dim_num, point_num, cluster_num, point, &
  cluster_center, cluster, energy )

!*****************************************************************************80
!
!! ENERGY_NORMAL computes the total energy of a given clustering.
!
!  Discussion:
!
!    For "normalized" vectors, we assume that an appropriate multiple
!    of the steady state solution has been subtracted, and that what
!    remains can be regarded as perturbations from 0 whose direction
!    (but not magnitude) is of importance.
!
!    The total energy function is the sum of the cluster energies.
!
!    The energy of a cluster is the sum of the "distances" of
!    each point in the cluster to the center point of the cluster.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    center points.
!
!    Input, integer CLUSTER(POINT_NUM), indicates the cluster to which
!    each data point belongs.
!
!    Output, real ( kind = 8 ) ENERGY, the total energy of the clustering.
!
  implicit none

  integer cluster_num
  integer dim_num
  integer point_num

  integer, dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  real ( kind = 8 ) distance_normal_sq
  real ( kind = 8 ) energy
  integer i
  integer j
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point

  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_energy(j) = cluster_energy(j) &
      + distance_normal_sq ( dim_num, point(1,i), cluster_center(1,j) )
  end do

  energy = sum ( cluster_energy(1:cluster_num) )

  return
end
subroutine energy_raw ( dim_num, point_num, cluster_num, point, &
  cluster_center, cluster, energy )

!*****************************************************************************80
!
!! ENERGY_RAW computes the total energy of a given clustering.
!
!  Discussion:
!
!    This routine is used with the raw data.  No normalization is
!    done to the data.
!
!    The total energy function is the sum of the cluster energies.
!
!    The energy of a cluster is the sum of the squares of the distances of
!    each point in the cluster to the center point of the cluster.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    center points.
!
!    Input, integer CLUSTER(POINT_NUM), indicates the cluster to which
!    each data point belongs.
!
!    Output, real ( kind = 8 ) ENERGY, the total energy of the clustering.
!
  implicit none

  integer cluster_num
  integer dim_num
  integer point_num

  integer, dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  real ( kind = 8 ) energy
  integer i
  integer j
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point

  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_energy(j) = cluster_energy(j) + sum ( &
      ( cluster_center(1:dim_num,j) - point(1:dim_num,i) )**2 )
  end do

  energy = sum ( cluster_energy(1:cluster_num) )

  return
end
subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer column_num
  logical got_one
  character ( len = * ) input_filename
  integer input_unit
  integer ios
  character ( len = 256 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( input_filename )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

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
  integer i
  integer lens

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
subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ROW_NUM, the number of rows found.
!
  implicit none

  integer bad_num
  integer comment_num
  integer ierror
  character ( len = * ) input_filename
  integer input_unit
  integer ios
  character ( len = 100 ) line
  integer record_num
  integer row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

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
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
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
subroutine hmeans_normal ( dim_num, point_num, cluster_num, it_max, it, &
  point, cluster, cluster_center, energy, null_cluster_policy )

!*****************************************************************************80
!
!! HMEANS_NORMAL seeks the minimal energy of a cluster of a given size.
!
!  Discussion:
!
!    This routine works with the normalized data.
!
!    The data for the H-Means problem is a set of N points X in
!    M-dimensions, and a desired number of clusters K.
!
!    The goal is to determine K points Z, called cluster centers, so that
!    if we associate each point X with its nearest Z value, we minimize
!    the standard deviation or cluster energy.  Writing CLUSTER(I) to
!    indicate the index of the nearest cluster center to point X(I), the
!    energy can be written as:
!
!      Energy = Sum ( 1 <= I <= N ) ( distance ( X(I), Z(CLUSTER(I)) ) )**2
!
!    where
!
!      distance ( X, Z ) = | sine ( angle between X and Z ) |
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez and Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, integer IT_MAX, the maximum number of iterations allowed.
!
!    Output, integer IT, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Output, integer CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the centers associated with the minimal energy clustering.
!
!    Output, real ( kind = 8 ) ENERGY, the total energy associated with
!    the minimal energy clustering.
!
!    Input, integer NULL_CLUSTER_POLICY, specifies what to do if a
!    null cluster is encountered.
!    0, do nothing.
!    1, reset center of null cluster to a random data point;
!    2, reset center of null cluster to a random point in the data hull;
!
  implicit none

  integer cluster_num
  integer dim_num
  integer point_num

  integer, dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  integer, dimension ( cluster_num ) :: cluster_population
  real ( kind = 8 ) column_sum
  logical, parameter :: debug = .false.
  real ( kind = 8 ) energy
  real ( kind = 8 ) factor(point_num)
  integer i
  integer it
  integer it_max
  integer j
  integer null_cluster_num
  integer null_cluster_policy
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer swap

  do it = 1, it_max
!
!  #1: Assign each point to the cluster of its nearest center.
!
    swap = 0

    do i = 1, point_num

      call nearest_cluster_normal ( dim_num, cluster_num, cluster_center, &
        point(1:dim_num,i), j )

      if ( j /= cluster(i) ) then
        cluster(i) = j
        swap = swap + 1
      end if

    end do
!
!  #2: Determine the energy of the new clustering with the current centers.
!
    call energy_normal ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, energy )

    if ( .false. ) then
      write ( *, '(i6,3x,14x,3x,g14.6)' ) it, energy
    end if
!
!  #3: Determine the populations of the new clusters.
!  If null clusters are not OK, we need to handle any empty clusters.
!  If any cluster is empty, set its center to a random point within
!  the convex hull of the data, reassign points, and repeat the entire
!  process.
!
    do

      cluster_population(1:cluster_num) = 0

      do i = 1, point_num
        j = cluster(i)
        cluster_population(j) = cluster_population(j) + 1
      end do

      null_cluster_num = 0

      do j = 1, cluster_num

        if ( cluster_population(j) == 0 ) then
          null_cluster_num = null_cluster_num + 1
          call random_number ( harvest = factor(1:point_num) )
          column_sum = sum ( factor(1:point_num) )
          factor(1:point_num) = factor(1:point_num) / column_sum
          cluster_center(1:dim_num,j) = &
            matmul ( point(1:dim_num,1:point_num), factor(1:point_num) )
        end if

      end do

      if ( null_cluster_num == 0 ) then
        exit
      end if

      if ( null_cluster_policy == 0 ) then
        exit
      end if

      if ( .false. ) then
        write ( *, '(a,i6)' ) &
          'HMEANS_NORMAL, number of empty clusters = ', null_cluster_num
      end if
!
!  Resort the points.
!
      do i = 1, point_num

        call nearest_cluster_normal ( dim_num, cluster_num, cluster_center, &
          point(1:dim_num,i), j )

        cluster(i) = j

      end do

    end do
!
!  #4: Recompute the cluster centers as the centroids of the points
!  in the cluster.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
      + point(1:dim_num,i)
  end do

  do i = 1, cluster_num
    cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
      / dble ( max ( cluster_population(i), 1 ) )
  end do
!
!  #5: Determine the energy of the current clustering with the new centers.
!
    call energy_normal ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, energy )

    if ( .false. ) then
      write ( *, '(i6,3x,g14.6)' ) it, energy
    end if

    if ( swap == 0 .and. 1 < it ) then
      exit
    end if

  end do

  return
end
subroutine hmeans_raw ( dim_num, point_num, cluster_num, it_max, it, &
  point, cluster, cluster_center, energy, null_cluster_policy, seed )

!*****************************************************************************80
!
!! HMEANS_RAW seeks the minimal energy of a cluster of a given size.
!
!  Discussion:
!
!    This routine works with the raw data, and does not do any
!    normalization.
!
!    The data for the H-Means problem is a set of N points X in
!    M-dimensions, and a desired number of clusters K.
!
!    The goal is to determine K points Z, called cluster centers, so that
!    if we associate each point X with its nearest Z value, we minimize
!    the standard deviation or cluster energy.  Writing CLUSTER(I) to
!    indicate the index of the nearest cluster center to point X(I), the
!    energy can be written as:
!
!      Energy = Sum ( 1 <= I <= N ) distance ( X(I), Z(CLUSTER(I)) )**2
!
!    where
!
!      distance ( X - Z ) = Sqrt ( Sum ( 1 <= J <= M ) ( X(J) - Z(J) )**2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez and Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, integer IT_MAX, the maximum number of iterations allowed.
!
!    Output, integer IT, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Output, integer CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the centers associated with the minimal energy clustering.
!
!    Output, real ( kind = 8 ) ENERGY, the total energy associated with
!    the minimal energy clustering.
!
!    Input, integer NULL_CLUSTER_POLICY, specifies what to do if a
!    null cluster is encountered.
!    0, do nothing.
!    1, reset center of null cluster to a random data point;
!    2, reset center of null cluster to a random point in the data hull;
!
!    Input/output, integer SEED, a seed for the random number generator.
!
  implicit none

  integer cluster_num
  integer dim_num
  integer point_num

  integer, dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  integer, dimension ( cluster_num ) :: cluster_population
  real ( kind = 8 ) column_sum
  logical, parameter :: debug = .false.
  real ( kind = 8 ) energy
  real ( kind = 8 ) factor(point_num)
  integer i
  integer i4_uniform
  integer it
  integer it_max
  integer j
  integer k
  integer null_cluster_num
  integer null_cluster_policy
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer seed
  integer swap_num

  do it = 1, it_max
!
!  #1: Assign each point to the cluster of its nearest center.
!
    call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
      cluster_center, point, cluster, swap_num )
!
!  #2: Determine the energy of the new clustering with the current centers.
!
    call energy_raw ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, energy )

    if ( .false. ) then
      write ( *, '(i6,3x,14x,3x,g14.6)' ) it, energy
    end if
!
!  #3: Determine the populations of the new clusters.
!  If null clusters are not OK, we need to handle any empty clusters.
!  If any cluster is empty, set its center to a random point within
!  the convex hull of the data, reassign points, and repeat the entire
!  process.
!
    do

      cluster_population(1:cluster_num) = 0

      do i = 1, point_num
        j = cluster(i)
        cluster_population(j) = cluster_population(j) + 1
      end do

      null_cluster_num = 0

      do j = 1, cluster_num
        if ( cluster_population(j) == 0 ) then
          null_cluster_num = null_cluster_num + 1
        end if
      end do

      if ( null_cluster_num == 0 ) then
        exit
      end if

      if ( null_cluster_policy == 0 ) then
        exit
      end if

      if ( debug ) then
        write ( *, '(a,i6)' ) &
          'HMEANS_RAW, number of empty clusters = ', null_cluster_num
      end if

      if ( null_cluster_policy == 1 ) then

        do j = 1, cluster_num
          if ( cluster_population(j) == 0 ) then
            k = i4_uniform ( 1, point_num, seed )
            cluster_center(1:dim_num,j) = point(1:dim_num,k)
          end if
        end do

      else if ( null_cluster_policy == 2 ) then

        do j = 1, cluster_num
          if ( cluster_population(j) == 0 ) then
            call random_number ( harvest = factor(1:point_num) )
            column_sum = sum ( factor(1:point_num) )
            factor(1:point_num) = factor(1:point_num) / column_sum
            cluster_center(1:dim_num,j) = &
              matmul ( point(1:dim_num,1:point_num), factor(1:point_num) )
          end if
        end do

      end if
!
!  Resort the points.
!
      call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
        cluster_center, point, cluster, swap_num )

    end do
!
!  #4: Recompute the cluster centers as the centroids of the points
!  in the cluster.
!
    cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

    do i = 1, point_num
      j = cluster(i)
      cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
        + point(1:dim_num,i)
    end do

    do i = 1, cluster_num
      cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
        / dble ( max ( cluster_population(i), 1 ) )
    end do
!
!  #5: Determine the energy of the current clustering with the new centers.
!
    call energy_raw ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, energy )

    if ( .false. ) then
      write ( *, '(i6,3x,g14.6)' ) it, energy
    end if

    if ( swap_num == 0 .and. 1 < it ) then
      exit
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
!    Output, integer VALUE, the value input by the user.
!
!    Output, integer IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ierror
  integer last
  character ( len = 80 ) line
  character ( len = * ) string
  integer value

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
!    Output, integer VALUE1, VALUE2, the values entered by the user.
!
!    Output, integer IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  character, parameter :: comma = ','
  integer ierror
  integer last
  integer last2
  character ( len = 80 ) line
  character, parameter :: space = ' '
  character ( len = * ) string
  integer value1
  integer value2

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
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
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
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine itable_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! ITABLE_DATA_READ reads data from an integer table file.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine
!    will return after reading N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Output, integer TABLE(M,N), the table data.
!
  implicit none

  integer m
  integer n

  integer ierror
  character ( len = * ) input_filename
  integer input_unit
  integer ios
  integer j
  character ( len = 255 ) line
  integer table(m,n)
  integer x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ITABLE_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i6)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i6)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i6)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_i4vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine itable_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! ITABLE_HEADER_READ reads the header from an integer table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer M, spatial dimension.
!
!    Output, integer N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer m
  integer n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an integer vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer big
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(i6,1x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(i6,1x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i6,i11)' ) i, a(i)
    end do
  end if

  return
end
subroutine kmeans_normal ( dim_num, point_num, cluster_num, it_max, it, &
  point, cluster, cluster_center, energy, null_cluster_policy )

!*****************************************************************************80
!
!! KMEANS_NORMAL tries to improve a partition of points.
!
!  Discussion:
!
!    This routine works with the "normalized" data.
!
!    For "normalized" vectors, we assume that an appropriate multiple
!    of the steady state solution has been subtracted, and that what
!    remains can be regarded as perturbations from 0 whose direction
!    (but not magnitude) is of importance.
!
!    The data for the K-Means problem is a set of N points X in
!    M-dimensions, and a desired number of clusters K.
!
!    The goal is to determine K points Z, called cluster centers, so that
!    if we associate each point X with its nearest Z value, we minimize
!    the standard deviation or cluster energy.  Writing CLUSTER(I) to
!    indicate the index of the nearest cluster center to point X(I), the
!    energy can be written as:
!
!      Energy = Sum ( 1 <= I <= N ) ( distance ( X(I), Z(CLUSTER(I)) ) )**2
!
!    where
!
!      distance ( X, Z ) = | sine ( angle between X and Z ) |
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez and Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, integer IT_MAX, the maximum number of iterations allowed.
!
!    Output, integer IT, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input/output, integer CLUSTER(POINT_NUM), the cluster to which
!    each point belongs.  On output, these may have been altered.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the centers associated with the clustering.  On output, these may
!    have been altered.
!
!    Input, integer NULL_CLUSTER_POLICY, specifies what to do if a
!    null cluster is encountered.
!    0, do nothing.
!    1, reset center of null cluster to a random data point;
!    2, reset center of null cluster to a random point in the data hull;
!
  implicit none

  integer cluster_num
  integer dim_num
  integer point_num

  integer ci
  integer cj
  integer, dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  integer, dimension ( cluster_num ) :: cluster_population
  real ( kind = 8 ) column_sum
  logical, parameter :: debug = .false.
  real ( kind = 8 ), dimension ( cluster_num ) :: distsq
  real ( kind = 8 ) energy
  real ( kind = 8 ), dimension ( point_num ) :: factor
  integer i
  integer it
  integer it_max
  integer j
  integer list(1)
  integer null_cluster_num
  integer null_cluster_policy
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer swap
  integer swap_total
!
!  For each observation, calculate the distance from each cluster
!  center, and assign to the nearest.
!
  do

    do i = 1, point_num

      call nearest_cluster_normal ( dim_num, cluster_num, cluster_center, &
        point(1:dim_num,i), j )

      if ( j /= cluster(i) ) then
        cluster(i) = j
      end if

    end do
!
!  Determine the cluster populations.
!
    cluster_population(1:cluster_num) = 0
    do i = 1, point_num
      j = cluster(i)
      cluster_population(j) = cluster_population(j) + 1
    end do
!
!  If a cluster is empty, give it a new cluster center, and restart
!  the process.
!  WARNING: This can take a long time!
!
    null_cluster_num = 0

    do j = 1, cluster_num
      if ( cluster_population(j) == 0 ) then
        null_cluster_num = null_cluster_num + 1
        call random_number ( harvest = factor(1:point_num) )
        column_sum = sum ( factor(1:point_num) )
        factor(1:point_num) = factor(1:point_num) / column_sum
        cluster_center(1:dim_num,j) = &
          matmul ( point(1:dim_num,1:point_num), factor(1:point_num) )
      end if
    end do

    if ( null_cluster_num == 0 ) then
      exit
    end if

    if ( null_cluster_policy == 0 ) then
      exit
    end if

    if ( debug ) then
      write ( *, '(a,i6)' ) &
        'KMEANS_NORMAL - Number of null clusters = ', null_cluster_num
    end if

  end do
!
!  Recompute the cluster centers as the centroids of the points
!  in the cluster.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
      + point(1:dim_num,i)
  end do

  do i = 1, cluster_num
    cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
      / dble ( max ( cluster_population(i), 1 ) )
  end do
!
!  Carry out the iteration.
!
  it = 0
  swap_total = 0

  do

    call energy_normal ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, energy )

    if ( .false. ) then
      write ( *, '(i6,3x,14x,3x,g14.6)' ) it, energy
    end if

    if ( it_max <= it ) then
      exit
    end if

    it = it + 1

    swap = 0

    do i = 1, point_num

      ci = cluster(i)

      if ( cluster_population(ci) <= 1 ) then
        cycle
      end if

      do cj = 1, cluster_num

        if ( cj == ci ) then

          distsq(cj) = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
            * dble ( cluster_population(cj) ) &
            / dble ( cluster_population(cj) - 1 )

        else if ( cluster_population(cj) == 0 ) then

          cluster_center(1:dim_num,cj) = point(1:dim_num,i)
          distsq(cj) = 0.0D+00

        else

          distsq(cj) = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
            * dble ( cluster_population(cj) ) &
            / dble ( cluster_population(cj) + 1 )

        end if

      end do
!
!  Find the index of the minimum value of DISTSQ.
!
      list = minloc ( distsq(1:cluster_num) )
!
!  If that's not the cluster to which point I now belongs, move it there.
!
      if ( list(1) == ci ) then
        cycle
      end if

      cj = list(1)

      cluster_center(1:dim_num,ci) = &
        ( dble ( cluster_population(ci) ) * cluster_center(1:dim_num,ci) &
        - point(1:dim_num,i) ) / dble ( cluster_population(ci) - 1 )

      cluster_center(1:dim_num,cj) = &
        ( dble ( cluster_population(cj) ) * cluster_center(1:dim_num,cj) &
        + point(1:dim_num,i) ) / dble ( cluster_population(cj) + 1 )

      cluster_population(ci) = cluster_population(ci) - 1
      cluster_population(cj) = cluster_population(cj) + 1

      cluster(i) = cj

      swap = swap + 1
      swap_total = swap_total + 1

    end do

    if ( swap == 0 ) then
      exit
    end if

  end do

  return
end
subroutine kmeans_raw ( dim_num, point_num, cluster_num, it_max, it, &
  point, cluster, cluster_center, energy, null_cluster_policy, seed )

!*****************************************************************************80
!
!! KMEANS_RAW tries to improve a partition of points.
!
!  Discussion:
!
!    This routine works with the raw data, and does not do any
!    normalization.
!
!    The data for the K-Means problem is a set of N points X in
!    M-dimensions, and a desired number of clusters K.
!
!    The goal is to determine K points Z, called cluster centers, so that
!    if we associate each point X with its nearest Z value, we minimize
!    the standard deviation or cluster energy.  Writing CLUSTER(I) to
!    indicate the index of the nearest cluster center to point X(I), the
!    energy can be written as:
!
!      Energy = Sum ( 1 <= I <= N ) distance ( X(I), Z(CLUSTER(I)) )**2
!
!    where
!
!      distance ( X - Z ) = Sqrt ( Sum ( 1 <= J <= M ) ( X(J) - Z(J) )**2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez and Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, integer IT_MAX, the maximum number of iterations allowed.
!
!    Output, integer IT, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input/output, integer CLUSTER(POINT_NUM), the cluster to which
!    each point belongs.  On output, these may have been altered.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the centers associated with the clustering.  On output, these may
!    have been altered.
!
!    Input, integer NULL_CLUSTER_POLICY, specifies what to do if a
!    null cluster is encountered.
!    0, do nothing.
!    1, reset center of null cluster to a random data point;
!    2, reset center of null cluster to a random point in the data hull;
!
!    Input/output, integer SEED, a seed for the random number generator.
!
  implicit none

  integer cluster_num
  integer dim_num
  integer point_num

  integer ci
  integer cj
  integer, dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  integer, dimension ( cluster_num ) :: cluster_population
  real ( kind = 8 ) column_sum
  logical, parameter :: debug = .false.
  real ( kind = 8 ), dimension ( cluster_num ) :: distsq
  real ( kind = 8 ) energy
  real ( kind = 8 ), dimension ( point_num ) :: factor
  integer i
  integer i4_uniform
  integer it
  integer it_max
  integer j
  integer k
  integer list(1)
  integer null_cluster_num
  integer null_cluster_policy
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer seed
  integer swap
  integer swap_num
  integer swap_total
!
!  For each observation, calculate the distance from each cluster
!  center, and assign to the nearest.
!
  do

    call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
      cluster_center, point, cluster, swap_num )
!
!  Determine the cluster populations.
!
    cluster_population(1:cluster_num) = 0
    do i = 1, point_num
      j = cluster(i)
      cluster_population(j) = cluster_population(j) + 1
    end do

    if ( debug ) then
      call i4vec_print ( cluster_num, cluster_population, &
        '  KMEANS_RAW: Input Pops:' )
    end if
!
!  If a cluster is empty, give it a new cluster center, and restart
!  the process.
!  WARNING: This can take a long time!
!
    null_cluster_num = 0

    do j = 1, cluster_num
      if ( cluster_population(j) == 0 ) then
        null_cluster_num = null_cluster_num + 1
      end if
    end do

    if ( null_cluster_num == 0 ) then
      exit
    end if

    if ( null_cluster_policy == 0 ) then
      exit
    end if

    if ( debug ) then
      write ( *, '(a,i6)' ) &
        'KMEANS_RAW - Number of null clusters = ', null_cluster_num
    end if

    if ( null_cluster_policy == 1 ) then

      do j = 1, cluster_num
        if ( cluster_population(j) == 0 ) then
          k = i4_uniform ( 1, point_num, seed )
          cluster_center(1:dim_num,j) = point(1:dim_num,k)
        end if
      end do

    else if ( null_cluster_policy == 2 ) then

      do j = 1, cluster_num
        if ( cluster_population(j) == 0 ) then
          call random_number ( harvest = factor(1:point_num) )
          column_sum = sum ( factor(1:point_num) )
          factor(1:point_num) = factor(1:point_num) / column_sum
          cluster_center(1:dim_num,j) = &
            matmul ( point(1:dim_num,1:point_num), factor(1:point_num) )
        end if
      end do

    end if

  end do
!
!  Calculate the new cluster centers.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
      + point(1:dim_num,i)
  end do

  do i = 1, cluster_num
    cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
      / dble ( max ( cluster_population(i), 1 ) )
  end do
!
!  Carry out the iteration.
!
  it = 0
  swap_total = 0

  do
!
!  Determine the energy.
!
    call energy_raw ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, energy )

    if ( .false. ) then
      write ( *, '(i6,3x,14x,3x,g14.6)' ) it, energy
    end if

    if ( it_max <= it ) then
      exit
    end if

    it = it + 1

    swap = 0

    do i = 1, point_num

      ci = cluster(i)

      if ( cluster_population(ci) <= 1 ) then
        cycle
      end if

      do cj = 1, cluster_num

        if ( cj == ci ) then

          distsq(cj) = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
            * dble ( cluster_population(cj) ) &
            / dble ( cluster_population(cj) - 1 )

        else if ( cluster_population(cj) == 0 ) then

          cluster_center(1:dim_num,cj) = point(1:dim_num,i)
          distsq(cj) = 0.0D+00

        else

          distsq(cj) = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
            * dble ( cluster_population(cj) ) &
            / dble ( cluster_population(cj) + 1 )

        end if

      end do
!
!  Find the index of the minimum value of DISTSQ.
!
      list = minloc ( distsq(1:cluster_num) )
!
!  If that's not the cluster to which point I now belongs, move it there.
!
      if ( list(1) == ci ) then
        cycle
      end if

      cj = list(1)

      cluster_center(1:dim_num,ci) = &
        ( dble ( cluster_population(ci) ) * cluster_center(1:dim_num,ci) &
        - point(1:dim_num,i) ) / dble ( cluster_population(ci) - 1 )

      cluster_center(1:dim_num,cj) = &
        ( dble ( cluster_population(cj) ) * cluster_center(1:dim_num,cj) &
        + point(1:dim_num,i) ) / dble ( cluster_population(cj) + 1 )

      cluster_population(ci) = cluster_population(ci) - 1

      cluster_population(cj) = cluster_population(cj) + 1

      cluster(i) = cj

      swap = swap + 1
      swap_total = swap_total + 1

    end do

    if ( swap == 0 ) then
      exit
    end if

  end do
!
!  Determine the cluster populations.
!
  if ( debug ) then

    cluster_population(1:cluster_num) = 0
    do i = 1, point_num
      j = cluster(i)
      cluster_population(j) = cluster_population(j) + 1
    end do

    call i4vec_print ( cluster_num, cluster_population, &
      '  KMEANS_RAW: Output Pops:' )

  end if

  return
end
subroutine mass_matrix ( node_num, npe, element_num, node, bandwidth, xy, a )

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
!    19 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NPE, the number of nodes per element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer NODE(NPE,ELEMENT_NUM), the nodes that make up each element.
!
!    Input, integer BANDWIDTH, the half bandwidth of the matrix.
!
!    Input, real ( kind = 8 ) XY(2,NODE_NUM), the nodes.
!
!    Output, real ( kind = 8 ) A(BANDWIDTH+1,NODE_NUM), the mass matrix.
!
  implicit none

  integer bandwidth
  integer element_num
  integer node_num
  integer npe
  integer, parameter :: nquad = 13

  real ( kind = 8 ) a(bandwidth+1,node_num)
  real ( kind = 8 ) area
  integer element
  real ( kind = 8 ) eta
  integer ip
  integer iq
  integer iquad
  integer jp
  integer jq
  integer node(npe,element_num)
  integer norder
  integer p1
  integer p2
  integer p3
  real ( kind = 8 ) refqbf
  integer rule
  real ( kind = 8 ) weight(nquad)
  real ( kind = 8 ) wi
  real ( kind = 8 ) wj
  real ( kind = 8 ) xy(2,node_num)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xtab(nquad)
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
        xy(1,p1) * ( xy(2,p2) - xy(2,p3) ) &
      + xy(1,p2) * ( xy(2,p3) - xy(2,p1) ) &
      + xy(1,p3) * ( xy(2,p1) - xy(2,p2) ) )

    if ( area == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MASS_MATRIX - Fatal error!'
      write ( *, '(a,i6)' ) '  Zero area for element ', element
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
subroutine nearest_cluster_normal ( dim_num, cluster_num, cluster_center, &
  point, nearest )

!*****************************************************************************80
!
!! NEAREST_CLUSTER_NORMAL finds the cluster nearest to a data point.
!
!  Discussion:
!
!    This routine uses the "normalized" data.
!
!    An appropriate multiple of the steady state solution has been
!    subtracted from each data vector, so that it satisfies the
!    zero boundary condition of the problem.
!
!    The distance between two solutions is the absolute value of
!    the sine of the angle between them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer CLUSTER_NUM, the number of center points.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    coordinates of the center points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM), the data point to be checked.
!
!    Output, integer NEAREST, the index of the center point closest to
!    the data point.
!
  implicit none

  integer cluster_num
  integer dim_num

  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_new
  real ( kind = 8 ) distance_normal_sq
  integer j
  integer nearest
  real ( kind = 8 ), dimension ( dim_num ) :: point

  dist = huge ( dist )
  nearest = 0

  do j = 1, cluster_num

    dist_new = distance_normal_sq ( dim_num, point, cluster_center(1,j) )

    if ( dist_new < dist ) then
      dist = dist_new
      nearest = j
    end if

  end do

  return
end
subroutine nearest_cluster_raw ( dim_num, point_num, cluster_num, &
  cluster_center, point, cluster, swap_num )

!*****************************************************************************80
!
!! NEAREST_CLUSTER_RAW finds the cluster nearest to a data point.
!
!  Discussion:
!
!    This routine uses the "raw" data.  Data is not normalized.
!    Distance is the usual Euclidean distance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer CLUSTER_NUM, the number of center points.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    coordinates of the center points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points
!    to be checked.
!
!    Input/output, integer CLUSTER(POINT_NUM).  On input, the cluster to
!    which each point was assigned.  On output, the cluster to which
!    each point has been reassigned.
!
!    Output, integer SWAP_NUM, the number of times a point was moved
!    from its input cluster to a different cluster.
!
  implicit none

  integer cluster_num
  integer dim_num
  integer point_num

  integer, dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_new
  integer i
  integer j
  integer nearest
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer swap_num

  swap_num = 0

  do i = 1, point_num

    dist = huge ( dist )
    nearest = 0

    do j = 1, cluster_num

      dist_new = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

      if ( dist_new < dist ) then
        dist = dist_new
        nearest = j
      end if

    end do

    if ( nearest /= cluster(i) ) then
      swap_num = swap_num + 1
    end if

    cluster(i) = nearest

  end do

  return
end
subroutine nearest_point ( dim_num, cluster_num, center, point, nearest )

!*****************************************************************************80
!
!! NEAREST_POINT finds the center point nearest a data point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer CLUSTER_NUM, the number of center points.
!
!    Input, real ( kind = 8 ) CENTER(DIM_NUM,CLUSTER_NUM), the
!    center points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM), the data point to be checked.
!
!    Output, integer NEAREST, the index of the center point closest to
!    the data point.
!
  implicit none

  integer cluster_num
  integer dim_num

  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: center
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_new
  integer j
  integer nearest
  real ( kind = 8 ), dimension ( dim_num ) :: point

  dist = huge ( dist )
  nearest = 0

  do j = 1, cluster_num

    dist_new = sum ( ( point(1:dim_num) - center(1:dim_num,j) )**2 )

    if ( dist_new < dist ) then
      dist = dist_new
      nearest = j
    end if

  end do

  return
end
subroutine point_generate ( point_dist, dim_num, r_min, r_max, point_num, &
  point )

!*****************************************************************************80
!
!! POINT_GENERATE generates data points for the problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_DIST, the point distribution to use.
!    1, use a uniform random distribution.
!    2, use a uniform grid of points.  (This hasn't been set up properly
!       except for 1D!).
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, real ( kind = 8 ) R_MIN(DIM_NUM), R_MAX(DIM_NUM), the coordinates
!    of the minimum and maximum corners of the region.
!
!    Input, integer POINT_NUM, the number of points to generate.
!
!    Output, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates
!    of the points.
!
  implicit none

  integer dim_num
  integer point_num

  integer i
  real ( kind = 8 ) point(dim_num,point_num)
  integer point_dist
  real ( kind = 8 ) r
  real ( kind = 8 ) r_max(dim_num)
  real ( kind = 8 ) r_min(dim_num)

  if ( point_dist == 1 ) then

    call random_number ( harvest = point(1:dim_num,1:point_num) )

    do i = 1, dim_num
      point(i,1:point_num) = r_min(i) + &
        point(i,1:point_num) * ( r_max(i) - r_min(i) )
    end do

  else if ( point_dist == 2 ) then

    do i = 1, point_num

      if ( 1 < point_num ) then
        r = dble ( i - 1 ) / ( point_num - 1 )
      else
        r = 0.5D+00
      end if

      point(1:dim_num,i) = ( 1.0D+00 - r ) * r_min(1:dim_num) &
                                     + r   * r_max(1:dim_num)
    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINT_GENERATE - Fatal error!'
    write ( *, '(a)' ) '  Meaningless input value of point distribution,'
    write ( *, '(a,i6)' ) '  POINT_DIST = ', point_dist
    stop

  end if

  return
end
subroutine point_print ( dim_num, point_num, point )

!*****************************************************************************80
!
!! POINT_PRINT prints out the values of the data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates of
!    the points.
!
  implicit none

  integer dim_num
  integer point_num

  integer i
  real ( kind = 8 ) point(dim_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Data points:'
  write ( *, '(a)' ) ' '

  do i = 1, point_num
    write ( *, '(i6)' ) i
    write ( *, '(8f10.4)' ) point(1:dim_num,i)
  end do

  return
end
function r8vec_norm2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM2 returns the 2-norm of a vector.
!
!  Discussion:
!
!    The vector 2-norm is defined as:
!
!      R8VEC_NORM2 = Sqrt ( Sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose 2-norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM2, the 2-norm of A.
!
  implicit none

  integer n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm2

  r8vec_norm2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
subroutine r8vec_range_input ( string, dim_num, value1, value2, ierror )

!*****************************************************************************80
!
!! R8VEC_RANGE_INPUT reads two DP vectors from the user, representing a range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Input, integer DIM_NUM, the number of dimensions.
!
!    Output, real ( kind = 8 ) VALUE1(DIM_NUM), VALUE2(DIM_NUM), the values
!    entered by the user.
!
!    Output, integer IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer dim_num

  integer ierror
  character ( len = * ) string
  real ( kind = 8 ) value1(dim_num)
  real ( kind = 8 ) value2(dim_num)

  ierror = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )
  read ( *, *, iostat = ierror ) value1(1:dim_num), value2(1:dim_num)

  return
end
subroutine r8vec_unit_euclidean ( n, a )

!*****************************************************************************80
!
!! R8VEC_UNIT_EUCLIDEAN normalizes a N-vector in the Euclidean norm.
!
!  Discussion;
!
!    The euclidean norm is also sometimes called the l2 or
!    least squares norm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.
!
!    Input/output, A(N), the vector to be normalized.
!
  implicit none

  integer n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) norm

  norm = sqrt ( sum ( a(1:n)**2 ) )

  if ( norm /= 0.0D+00 ) then
    a(1:n) = a(1:n) / norm
  end if

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
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator.
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  integer i
  integer seed
  integer, allocatable :: seed_vector(:)
  integer seed_size
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
!    Input, integer IQ, the local node number, between 1 and
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
  integer iq
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
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer get
  integer put
  integer nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

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
  integer i
  integer len1
  integer len2
  integer lenc
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
!    Output, integer IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ierror
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
!    Input, integer I, an integer to be converted.
!
!    Output, character ( len = 11 ) S_OF_I4, the representation of the
!    integer.  The integer will be left-justified.
!
  implicit none

  character c
  integer i
  integer idig
  integer ihi
  integer ilo
  integer ipos
  integer ival
  integer j
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
  integer i
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
!    Output, integer IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
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
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an integer vector from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, integer IVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer n

  integer i
  integer ierror
  integer ilo
  integer ivec(n)
  integer length
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

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
!    part of the real ( kind = 8 ) number.
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
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
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
      else if ( 1 < ihave ) then
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
      else if ( 6 <= ihave .and. ihave <= 8 ) then
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
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
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
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer n

  integer i
  integer ierror
  integer ilo
  integer lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

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
!    Output, integer NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer i
  integer lens
  integer nword
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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = * ) string
  character ( len = 10 ) time
  integer values(8)
  integer y
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
!! TRIANGLE_UNIT_SET sets a quadrature rule in a unit triangle.
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
!    Input, integer RULE, the index of the rule.
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
!    Output, integer NORDER, the order of the rule.
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
  integer norder
  integer rule
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
    write ( *, '(a,i6)' ) '  Illegal value of RULE = ', rule
    stop

  end if

  return
end
