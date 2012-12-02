program main

!*****************************************************************************80
!
!  Purpose:
!
!    MAIN is the main program for HEATED_PLATE_WORKSHARE.
!
!  Discussion:
!
!    This code solves the steady state heat equation on a rectangular region.
!
!    The sequential version of this program needs approximately
!    18/eps iterations to complete. 
!
!
!    The physical region, and the boundary conditions, are suggested
!    by this diagram;
!
!                   W = 0
!             +------------------+
!             |                  |
!    W = 100  |                  | W = 100
!             |                  |
!             +------------------+
!                   W = 100
!
!    The region is covered with a grid of M by N nodes, and an N by N
!    array W is used to record the temperature.  The correspondence between
!    array indices and locations in the region is suggested by giving the
!    indices of the four corners:
!
!                  I = 0
!          [0][0]-------------[0][N-1]
!             |                  |
!      J = 0  |                  |  J = N-1
!             |                  |
!        [M-1][0]-----------[M-1][N-1]
!                  I = M-1
!
!    The steady state solution to the discrete heat equation satisfies the
!    following condition at an interior grid point:
!
!      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )
!
!    where "Central" is the index of the grid point, "North" is the index
!    of its immediate neighbor to the "north", and so on.
!   
!    Given an approximate solution of the steady state heat equation, a
!    "better" solution is given by replacing each interior point by the
!    average of its 4 neighbors - in other words, by using the condition
!    as an ASSIGNMENT statement:
!
!      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )
!
!    If this process is repeated often enough, the difference between successive 
!    estimates of the solution will go to zero.
!
!    This program carries out such an iteration, using a tolerance specified by
!    the user, and writes the final estimate of the solution to a file that can
!    be used for graphic processing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2010
!
!  Author:
!
!    Original FORTRAN90 version by Michael Quinn.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Quinn,
!    Parallel Programming in C with MPI and OpenMP,
!    McGraw-Hill, 2004,
!    ISBN13: 978-0071232654,
!    LC: QA76.73.C15.Q55.
!
!  Local parameters:
!
!    Local, double precision DIFF, the norm of the change in the solution from 
!    one iteration to the next.
!
!    Local, double precision MEAN, the average of the boundary values, used 
!    to initialize the values of the solution in the interior.
!
!    Local, double precision U(M,N), the solution at the previous iteration.
!
!    Local, double precision W(M,N), the solution computed at the latest 
!    iteration.
!
  use omp_lib

  implicit none

  integer, parameter :: m = 500
  integer, parameter :: n = 500

  double precision diff
  double precision :: eps = 0.001D+00
  integer i
  integer iterations
  integer iterations_print
  integer j
  integer k
  double precision mean
  double precision u(m,n)
  double precision w(m,n)
  double precision wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEATED_PLATE_WORKSHARE'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
  write ( *, '(a)' ) '  A program to solve for the steady state temperature distribution'
  write ( *, '(a)' ) '  over a rectangular plate.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a)' ) '  Spatial grid of ', m, ' by ', n, ' points.'
  write ( *, '(a,g14.6)' ) &
    '  The iteration will repeat until the change is <= ', eps
  write ( *, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) &
    '  The number of threads available    = ', omp_get_max_threads ( )
!
!  Set the boundary values, which don't change. 
!
!$omp parallel shared ( mean, w )

!$omp workshare
  w(2:m-1,1) = 100.0
  w(2:m-1,n) = 100.0
  w(m,1:n) = 100.0
  w(1,1:n) =   0.0
!$omp end workshare
!
!  Average the boundary values, to come up with a reasonable
!  initial value for the interior.
!
!  OpenMP comment: the previous WORKSHARE must be complete
!  before we can start this one.  That's why they are in
!  separate WORKSHARE's.
!
!$omp workshare
  mean = ( &
      sum ( w(2:m-1,1) ) &
    + sum ( w(2:m-1,n) ) &
    + sum ( w(m,1:n)   ) &
    + sum ( w(1,1:n)   ) ) &
    / dble ( 2 * m + 2 * n - 4 )
!$omp end workshare
!
!  Initialize the interior solution to the mean value.
!
!$omp workshare
  w(2:m-1,2:n-1) = mean
!$omp end workshare

!$omp end parallel

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  MEAN = ', mean
!
!  Iterate until the new solution W differs from the old solution U
!  by no more than EPS.
!
  iterations = 0
  iterations_print = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Iteration  Change'
  write ( *, '(a)' ) ' '

  wtime = omp_get_wtime ( )

  diff = eps

  do while ( eps <= diff )

!$omp parallel shared ( diff, u, w )

!$omp workshare
  u = w
!$omp end workshare

!$omp workshare
  forall ( i = 2 : m - 1, j = 2 : n - 1 )
    w(i,j) = 0.25 * ( &
        u(i-1,j) &
      + u(i+1,j) &
      + u(i,j-1) &
      + u(i,j+1) )
  end forall
!$omp end workshare

!$omp workshare
    diff = maxval ( abs ( u - w ) )
!$omp end workshare

!$omp end parallel

    iterations = iterations + 1

    if ( iterations == iterations_print ) then
      write ( *, '(2x,i8,2x,g14.6)' ) iterations, diff
      iterations_print = 2 * iterations_print
    end if

  end do

  wtime = omp_get_wtime ( ) - wtime

  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8,2x,g14.6)' ) iterations, diff
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Error tolerance achieved.'
  write ( *, '(a,g14.6)' ) '  Wall clock time = ', wtime
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEATED_PLATE_WORKSHARE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
