/* 
  Discussion:

    This is a modification of code supplied by David Crawford of Sandia Labs.
    The routine DIATOM_SETUP() is intended to be called by a FORTRAN program,
    and to initialize and set up the DIATOM geometry structures.

    This routine functions as an interface between higher level FORTRAN
    routines, and the DIATOM library, which is written primarily in C.

*/

#include <stdlib.h>
#include <stdio.h>

#define INTEGER int
#define REAL double

/* Declare external functions. Usage is explained in main(). */

extern void diatom_mpp(INTEGER nprocs, INTEGER nid,
            void (*MPSIRY)(INTEGER *NVALS, INTEGER *IARRAY),
            void (*MPRMAX)(REAL *X1, REAL *X2));

extern void diatom_msg(
            void (*Adiom_log_message)(char* c_array),
            void (*Adiom_fatal_error)(char* c_array),
            void (*Adiom_bombed)(char* c_array));

extern void diatom_mesh(INTEGER Acoord_type, INTEGER Amesh_type,
                 REAL xmin, REAL xmax,
                 REAL ymin, REAL ymax,
                 REAL zmin, REAL zmax );

extern diatom_parse_input(char *, char *);

extern INTEGER diatom_point_test(REAL x, REAL y, REAL z, REAL dr, REAL *value);

void diatom_per_cycle_init(INTEGER init, REAL time, REAL dt, REAL *x1, 
  REAL *y1, REAL *z1, REAL *x2, REAL *y2, REAL *z2, INTEGER *err);

FILE *diatom_log_file;

void log_message(char *message)
{
  fprintf(diatom_log_file,"%s\n",message);
}

void fatal_message(char *message)
{
  printf("%s\n",message);
  exit(0);
}

void bomb_message(char *message)
{
  printf("%s\n",message);
  exit(1);
}

void dummy_broadcast_int(int *nvals, int *array)
{
}

void dummy_max_real(double *x1, double *x2)
{
  *x2 = *x1;
}

double my_random()
/*
  JVB: 

  The original code was
    "return (double) random() / RAND_MAX;"
  but this was not returning values between 0 and 1.
  I guess RAND_MAX, which is 32767 on my DEC ALPHA, is set for a signed 
  integer*2 value.
*/
{
  return (double) random() / ( 2 * RAND_MAX * RAND_MAX );
}

void diatom_setup()
{
  FILE *boundary;
  int boundary_n;
  double dr;
  int err;
  int exterior_n;
  int flag;
  int i;
  int igm;
  int interior_n;
  FILE *interior;
  int j;
  int mesh;
  double tmp;
  double value;
  double x;
  double xmax;
  double xmin;
  double y;
  double ymax;
  double ymin;
  double z;
  double zmax;
  double zmin;
/* 
  Establish callbacks for message passing...dummies for now 
*/

  diatom_mpp ( 1, 0, dummy_broadcast_int, dummy_max_real );

/* 
  Establish callbacks for printing diagnostic and error messages
  See the routines above and modify as needed. 
*/

  diatom_msg ( log_message, fatal_message, bomb_message );

/* 
  Establish the global coordinate system:
*/

/*
  Specify the geometry.
  10 = 1-d rectilinear
  11 = 1-d cylindrical
  12 = 1-d spherical
  20 = 2-d rectilinear
  21 = 2-d cylindrical
  30 = 3-d rectilinear
*/
  igm = 30;
/*
  Set the mesh type:

    0 (bricks only for now - others will be added later)
*/
  mesh = 0;
/*
  Set the limits of the bounding box.
*/
  xmin = 0.0;
  xmax = 100.0;
  ymin = 0.0;
  ymax = 100.0;
  zmin = 0.0;
  zmax = 20.0;

  printf ( "\n" );
  printf ( "DIATOM_TEST calling DIATOM_MESH.\n" );
  printf ( "    Requesting IGM = %d.\n", igm );
  printf ( "    Requesting MESH = %d.\n", mesh );
  printf ( "    [XMIN, XMAX] = [%f, %f]\n", xmin, xmax );
  printf ( "    [YMIN, YMAX] = [%f, %f]\n", ymin, ymax );
  printf ( "    [ZMIN, ZMAX] = [%f, %f]\n", zmin, zmax );

  diatom_mesh ( igm, mesh, xmin, xmax, ymin, ymax, zmin, zmax );
/* 
  Open a log file for diatom diagnostics only...not important
  but comforting to know that diatom is working properly. 
*/

  diatom_log_file = fopen ( "diatom_test.log", "w" );

/* 
  Parse the user input.

  DIATOM stores the user geometry in a linked list representation.
  Diagnostic information from the parsing is sent to the log file
  via the message-passing established in log_message routine  above.

  The first parameter below is the filename containing the user-defined
  geometry definitions. The second parameter is the keyword that delineates
  the beginning of user input within the file. In the example below, user
  input is contained within the file "diatom_test.in" starting with a
  line containing "diatom" and ending with a line containing "endd"
  ("end"+first character of "diatom"). 
*/

  diatom_parse_input ( "diatom_test.in", "diatom" );

/* 
  Perform time-dependent initialization. For setting up
  initial conditions, as we're doing here, we call this just once with
  TIME and DT both zero (the 2nd and 3rd parameters). 
*/

  diatom_per_cycle_init ( 1, 0, 0, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &err );

/* 
  Close the log file 
*/

  fclose(diatom_log_file);

  printf ( "\n" );
  printf ( "DIATOM_SETUP:\n" );
  printf ( "  The DIATOM geometry has been set up!\n" );

  return;
}
