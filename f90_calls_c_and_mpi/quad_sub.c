# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <mpi.h>

double quad_compute ( int n, int *fComm );
double f ( double x );

/******************************************************************************/

double quad_compute ( int n, int *fComm )

/******************************************************************************/
/*
  Purpose:

    QUAD_COMPUTE helps to estimate an integral.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2011

  Author:

    John Burkardt

  Reference:

    William Gropp, Ewing Lusk, Anthony Skjellum,
    Using MPI: Portable Parallel Programming with the
    Message-Passing Interface,
    Second Edition,
    MIT Press, 1999,
    ISBN: 0262571323.

  Parameters:

    Input, int N, the number of evaluation points to use.

    Input, int *FCOMM, the FORTRAN MPI Communicator ID.

    Output, double QUAD_COMPUTE, the integral estimate made by this process.
*/
{
  MPI_Comm comm;
  double h;
  int i;
  int id;
  int ierr;
  int n_part;
  int p;
  double q_part;
  double x;
/*
  Convert Fortran MPI Communicator ID to C.
*/
  comm = MPI_Comm_f2c ( *fComm );
/*
  Determine this processes's rank.
*/
  ierr = MPI_Comm_rank ( comm, &id );
/*
  Determine the number of processes.
*/
  ierr = MPI_Comm_size ( comm, &p );
/*
  The process integrates F(X) over a subinterval determined by its process ID.
*/
  h = 1.0 / ( double ) n;

  q_part = 0.0;
  n_part = 0;

  for ( i = id + 1; i <= n; i = i + p ) 
  {
    x = ( double ) ( 2 * i - 1 )
      / ( double ) ( 2 * n     );

    n_part = n_part + 1;
    q_part = q_part + f ( x );
  }
  q_part = q_part * h;

  printf ( "\n" );
  printf ( "QUAD_COMPUTE - Process %d:\n", id );
  printf ( "  Points used = %d\n", n_part );
  printf ( "  Estimate %f.\n", q_part );

  return q_part;
}
/******************************************************************************/

double f ( double x )

/******************************************************************************/
/*
  Purpose:

    F evaluates the function F(X) which we are integrating.

  Discussion:

    Integral ( 0 <= X <= 1 ) 4/(1+X*X) dX = PI

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which to evaluate F.

    Output, double F, the value of F(X).
*/
{
  double value;

  value = 4.0 / ( 1.0 + x * x );

  return value;
}
           
