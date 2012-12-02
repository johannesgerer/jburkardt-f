program main

!*****************************************************************************80
!
!! MAIN is the main program for BUFFON_LAPLACE.
!
!  Discussion:
!
!    This program uses MPI to do a Buffon-Laplace simulation in parallel.
!
!    This is an example of an "embarassingly parallel" computation.  Each
!    processor does the same randomized computation, counting "successes",
!    and returning the number of successes to the master process.
!
!    The particular point made by this example is that each process must
!    use a different stream of random numbers, and that this can be done
!    if each process uses a different seed to initialize the random number
!    generator, and that this, in turn, can be done by simply adding the
!    process's rank to a common seed value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!    Sudarshan Raghunathan,
!    Making a Supercomputer Do What You Want: High Level Tools for
!    Parallel Programming,
!    Computing in Science and Engineering,
!    Volume 8, Number 5, September/October 2006, pages 70-80.
!
  implicit none

  include 'mpif.h'

  real ( kind = 8 ) :: a = 1.0D+00
  real ( kind = 8 ) :: b = 1.0D+00
  integer ( kind = 4 ) :: buffon_laplace_simulate
  integer ( kind = 4 ):: hit_num
  integer ( kind = 4 ):: hit_total
  integer ( kind = 4 ):: ierr
  real ( kind = 8 ) :: l = 1.0D+00
  integer, parameter :: master = 0
  real ( kind = 8 ) :: pdf_estimate
  real ( kind = 8 ), parameter :: pi = 3.141592653589793238462643D+00
  real ( kind = 8 ) :: pi_error
  real ( kind = 8 ) :: pi_estimate
  integer ( kind = 4 ):: process_num
  integer ( kind = 4 ):: process_rank
  integer ( kind = 4 ):: seed
  integer ( kind = 4 ):: trial_num = 100000
  integer ( kind = 4 ):: trial_total
  real ( kind = 8 ) :: random_value
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )

  if ( ierr /= MPI_SUCCESS ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BUFFON_LAPLACE: Warning!'
    write ( *, '(a,i8)' ) '  MPI_INIT returns IERR = ', ierr
    call MPI_Finalize ( ierr )
    stop
  end if
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, process_num, ierr )
!
!  Get the rank of this processor.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, process_rank, ierr )
!
!  The master process sets the value of the parameters,
!  and broadcasts them.
!
  if ( process_rank == master ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BUFFON_LAPLACE - Master process:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An MPI example program to estimate PI'
    write ( *, '(a)' ) '  using the Buffon-Laplace needle experiment.'
    write ( *, '(a)' ) '  On a grid of cells of  width A and height B,'
    write ( *, '(a)' ) '  a needle of length L is dropped at random.'
    write ( *, '(a)' ) '  We count the number of times it crosses'
    write ( *, '(a)' ) '  at least one grid line, and use this to estimate '
    write ( *, '(a)' ) '  the value of PI.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of processes is ', process_num
    write ( *, '(a)' ) ' '
    write ( *, '(a,f14.6)' ) '  Cell width A =    ', a
    write ( *, '(a,f14.6)' ) '  Cell height B =   ', b
    write ( *, '(a,f14.6)' ) '  Needle length L = ', l
  end if
!
!  Each process sets a random number seed.
!
  seed = 123456789 + process_rank * 100
!
!  The process of setting the random number seed in FORTRAN90 is, 
!  at least, uniformly specified by the language,
!  but lamentably Byzantine and cumbersome.  
!
!  I put all this rigmarole in a subroutine.
!
  call random_initialize ( seed )
!
!  Just to make sure that we're all doing different things, have each
!  process print out its rank, seed value, and a first test random value.
!
  call random_number ( harvest = random_value )

  write ( *, '(2x,i8,2x,i12,2x,g14.6)' ) process_rank, seed, random_value
!
!  Each process now carries out TRIAL_NUM trials, and then
!  sends the value back to the master process.
!
  hit_num = buffon_laplace_simulate ( a, b, l, trial_num ) 

  call MPI_Reduce ( hit_num, hit_total, 1, MPI_INTEGER, MPI_SUM, master, &
    MPI_COMM_WORLD, ierr )
!
!  The master process can now estimate PI.
!
  if ( process_rank == master ) then

    trial_total = trial_num * process_num

    pdf_estimate = real ( hit_total, kind = 8 ) / real ( trial_total, kind = 8 )

    if ( hit_total == 0 ) then

      pi_estimate = huge ( pi_estimate )

    else

      pi_estimate = l * ( 2.0D+00 * ( a + b ) - l ) / ( a * b * pdf_estimate )

    end if

    pi_error = abs ( pi - pi_estimate )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Trials      Hits    Estimated PDF  ' // &
      '     Estimated Pi         Error'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8,2x,i8,2x,g20.12,2x,g20.12,2x,g20.12)' ) &
      trial_total, hit_total, pdf_estimate, pi_estimate, pi_error

  end if
!
!  Terminate MPI.
!
  call MPI_Finalize ( ierr )
!
!  Terminate.
!
  if ( process_rank == master ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BUFFON_LAPLACE - Master process:'
    write ( *, '(a)' ) '  Normal end of execution.'
  end if

  stop
end
function buffon_laplace_simulate ( a, b, l, trial_num )

!*****************************************************************************80
!
!! BUFFON_LAPLACE_SIMULATE simulates a Buffon-Laplace needle experiment.
!
!  Discussion:
!
!    In the Buffon-Laplace needle experiment, we suppose that the plane has been
!    tiled into a grid of rectangles of width A and height B, and that a
!    needle of length L is dropped "at random" onto this grid.
!
!    We may assume that one end, the "eye" of the needle falls at the point
!    (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].
!
!    ANGLE, the angle that the needle makes is taken to be uniformly random.
!    The point of the needle, (X2,Y2), therefore lies at
!
!      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
!
!    The needle will have crossed at least one grid line if any of the
!    following are true:
!
!      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.
!
!    This routine simulates the tossing of the needle, and returns the number
!    of times that the needle crossed at least one grid line.
!
!    If L is larger than sqrt ( A*A + B*B ), then the needle will
!    cross every time, and the computation is uninteresting.  However, if
!    L is smaller than this limit, then the probability of a crossing on
!    a single trial is
!
!      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( PI * A * B )
!
!    and therefore, a record of the number of hits for a given number of
!    trials can be used as a very roundabout way of estimating PI.
!    (Particularly roundabout, since we actually will use a good value of
!    PI in order to pick the random angles!)
!
!    Note that this routine will try to generate 5 * TRIAL_NUM random
!    double precision values at one time, using automatic arrays.
!    When I tried this with TRIAL_NUM = 1,000,000, the program failed,
!    because of internal system limits on such arrays.
!
!    Such a problem could be avoided by using a DO loop running through
!    each trial individually, but this tend to run much more slowly than
!    necessary.
!
!    Since this routine invokes the FORTRAN90 random number generator,
!    the user should initialize the random number generator, particularly
!    if it is desired to control whether the sequence is to be varied
!    or repeated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Sudarshan Raghunathan,
!    Making a Supercomputer Do What You Want: High Level Tools for
!    Parallel Programming,
!    Computing in Science and Engineering,
!    Volume 8, Number 5, September/October 2006, pages 70-80.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the horizontal and vertical dimensions
!    of each cell of the grid.  0 <= A, 0 <= B.
!
!    Input, real ( kind = 8 ) L, the length of the needle.
!    0 <= L <= min ( A, B ).
!
!    Input, integer TRIAL_NUM, the number of times the needle is
!    to be dropped onto the grid.
!
!    Output, integer BUFFON_LAPLACE_SIMULATE, the number of times the needle
!    crossed at least one line of the grid of cells.
!
!  Local Parameters:
!
!    Local, integer BATCH_SIZE, specifies the number of trials to be done
!    in a single batch.  Setting BATCH_SIZE to 1 will be very slow.
!    Replacing it by TRIAL_NUM would be fine except that your system
!    may have a limit on the size of automatic arrays.  We have set a default
!    value of 10,000 here which should be large enough to be efficient
!    but small enough not to annoy the system.
!
  implicit none

  integer, parameter :: batch_size = 10000
  integer trial_num

  real ( kind = 8 ) a
  real ( kind = 8 ) angle(batch_size)
  real ( kind = 8 ) b
  integer batch
  integer buffon_laplace_simulate
  integer hits
  real ( kind = 8 ) l
  integer n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793238462643D+00
  real ( kind = 8 ) x1(batch_size)
  real ( kind = 8 ) x2(batch_size)
  real ( kind = 8 ) y1(batch_size)
  real ( kind = 8 ) y2(batch_size)

  hits = 0

  do batch = 1, trial_num, batch_size

    n = min ( batch_size, trial_num + 1 - batch )
!
!  Randomly choose the location of the eye of the needle in [0,0]x[A,B],
!  and the angle the needle makes.
!
    call random_number ( harvest = x1(1:n) )
    call random_number ( harvest = y1(1:n) )
    call random_number ( harvest = angle(1:n) )

    x1(1:n) = a * x1(1:n)
    y1(1:n) = b * y1(1:n)
    angle(1:n) = 2.0D+00 * pi * angle(1:n)
!
!  Compute the location of the point of the needle.
!
    x2(1:n) = x1(1:n) + l * cos ( angle(1:n) )
    y2(1:n) = y1(1:n) + l * sin ( angle(1:n) )
!
!  Count the end locations that lie outside the cell.
!
    hits = hits + count (      x2(1:n) <= 0.0 .or. &
                          a <= x2(1:n)        .or. &
                               y2(1:n) <= 0.0 .or. &
                          b <= y2(1:n) )

  end do

  buffon_laplace_simulate = hits

  return
end
subroutine random_initialize ( seed_input )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee.
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator
!    a bit, this routine goes through the tedious process of getting the
!    size of the random number seed, making up values based on the current
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer SEED_INPUT, the user's "suggestion" for a seed.
!    However, if the input value is 0, the routine will come up with
!    its own "suggestion", based on the system clock.
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  logical, parameter :: debug = .false.
  integer i
  integer seed
  integer seed_input
  integer, allocatable :: seed_vector(:)
  integer seed_size
  real ( kind = 8 ) t
  integer, parameter :: warm_up = 100

  seed = seed_input
!
!  Initialize the random seed routine.
!
  call random_seed ( )
!
!  Determine the size of the random number seed vector.
!
  call random_seed ( size = seed_size )
!
!  Allocate a vector of the right size to be used as a random seed.
!
  allocate ( seed_vector(seed_size) )
!
!  If the user supplied a SEED value, use that.
!
!  Otherwise, use the system clock value to make up a value that is
!  likely to change based on when this routine is called.
!
  if ( seed /= 0 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
        seed
    end if

  end if
!
!  Set the seed vector.  We don't know the significance of the
!  individual entries of the internal seed vector, so we'll just set
!  all entries to SEED.
!
  seed_vector(1:seed_size) = seed
!
!  Now call RANDOM_SEED, and tell it to use this seed vector.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times just to "warm it up".
!
  do i = 1, warm_up
    call random_number ( harvest = t )
  end do

  return
end
