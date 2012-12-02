# include <cstdlib>
# include <iostream>
# include <cmath>
# include <mpi.h>

extern "C" 
{
  double quad_compute ( int n );
  double f ( double x );
}
//***************************************************************************80

double quad_compute ( int n )

//***************************************************************************80
//
//  Purpose:
//
//    QUAD_COMPUTE helps to estimate an integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: Portable Parallel Programming with the
//    Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points to use.
//
//    Output, double QUAD_COMPUTE, the integral estimate made by this process.
//
{
  double h;
  int i;
  int id;
  int n_part;
  int p;
  double q_part;
  double x;
//
//  Determine this processes's rank.
//
  id = MPI::COMM_WORLD.Get_rank ( );
//
//  Determine the number of processes.
//
  p = MPI::COMM_WORLD.Get_size (  );
//
//  The process integrates F(X) over a subinterval determined by its process ID.
//
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
//****************************************************************************80

double f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F evaluates the function F(X) which we are integrating.
//
//  Discussion:
//
//    Integral ( 0 <= X <= 1 ) 4/(1+X*X) dX = PI
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which to evaluate F.
//
//    Output, double F, the value of F(X).
//
{
  double value;

  value = 4.0 / ( 1.0 + x * x );

  return value;
}
           
