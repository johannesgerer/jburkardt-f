#include <math.h>
#include "diom_globals.h"
#include "diom_procedures.h"

void diatom_point_test2 ( REAL *x, REAL *y, REAL *z, REAL *dr, REAL *value,
  INTEGER *ival )
{
/* 
  DIATOM_POINT_TEST2

  Discussion:

    This is a modified version of DIATOM_POINT_TEST, as supplied by David
    Crawford.  In this version, the arguments X, Y, Z and DR are pointers.
    This facilitates calling the routine from FORTRAN.

    This routine functions as an interface between higher level FORTRAN
    routines, and the DIATOM library, which is written primarily in C.

  Modified:

    12 April 2001

  Parameters:

    Input, REAL *X, *Y, *Z, pointers to the (X,Y,Z) coordinates of a point.

    Input, REAL *DR, pointer to a tolerance used when determining if a
    point is in, on, or outside the region.

    Input, REAL *VALUE, pointer to an output quantity, a material density 
    to be computed by this routine.

    Input, INTEGER *IVAL, pointer to an output quantity, an index to be 
    computed by this routine, which is 1 if the point (X,Y,Z) is in the 
    region, -1 if it is on the boundary of the region, and 0 if it is 
    outside the region.
*/
 static struct diom *curr_diom;
 struct geom_object *curr_geom;
 REAL dt = 0;
 INTEGER i;
 int init = 1;
 REAL io;
 INTEGER kbad;
 INTEGER mat_id;
 REAL *mmtmp;
 REAL time = 0;
 REAL tmp;
 REAL vf;
 REAL vf_all;
 REAL x1;
 REAL x2;
 REAL y1;
 REAL y2;
 REAL z1;
 REAL z2;

/* printf ( "DIATOM_POINT_TEST2: xyz = %f %f %f.\n", *x, *y, *z ); */
/*  printf ( "DIATOM_POINT_TEST2: dr = %f.\n", *dr ); */

 x1 = *x - *dr;
 y1 = *y - *dr;
 z1 = *z - *dr;

 x2 = *x + *dr;
 y2 = *y + *dr;
 z2 = *z + *dr;

 curr_diom = head_diom;

 if ( curr_diom != NULL ) {
  curr_geom = curr_diom->geom_list;
 } else {
  *ival = 0;
  printf ( "\n" );
  printf ( "DIATOM_POINT_TEST2: Bad News!\n" );
  printf ( "  The data structure is not set up.\n" );
  return;
 }

 vf_all=0;
 *value=0;

 while (curr_diom != NULL) {
  vf=0;
  if (((init == 1) && (curr_diom->off == time)) || 
      ((init == 1) && (curr_diom->off < time) && (curr_diom->speed == 0)) ||
      ((init == 1) && (curr_diom->anywhere == 2)) ||
      ((init == 0) && (time+dt >= curr_diom->on) && (time < curr_diom->off))) {
   if (diom_brick_VF != NULL) {
     if ((curr_diom->anywhere == 1) || ((curr_diom->anywhere == 2) && (init == 1)) ||
         (in_ureg(number_dimensions,x1, y1, z1, x2, y2, z2) != 0)) {
       curr_geom = curr_diom->geom_list;
       while (curr_geom != NULL) {
        mmtmp = curr_geom->minmax;
        if ((mmtmp[0]<x2) && (mmtmp[3]>x1) &&
          (((mmtmp[1]<y2) && (mmtmp[4]>y1)) || (number_dimensions<2)) &&
          (((mmtmp[2]<z2) && (mmtmp[5]>z1)) || (number_dimensions<3))) {
         /*curr_geom->GSUB(curr_geom->G,&curr_geom->NG,x,y,z,dr,dr,dr,&io);
         vf += curr_geom->VFADD*fabs(io);*/
         vf += curr_geom->VFADD*diom_brick_VF(curr_geom, x1, y1, z1, x2, y2, z2, 0);
         if (curr_diom->IP[1] < 2) {
          if (vf < 0.0) vf = 0.0;
          if (vf > 1.0) vf = 1.0;
         } else {
          if (vf < -1.0) vf = -1.0;
          if (vf > 0.0) vf = 0.0;
         }
        }
        curr_geom = curr_geom->next;
       }
     }
   }
  }
  if (vf>0) {
   if (vf_all<1) {
    tmp = 1-vf_all;
    if (tmp<vf) vf=tmp;
    *value += vf*curr_diom->P[0];
    vf_all += vf;
   }
  } else {
   if (vf_all>0) {
    tmp = vf_all;
    vf = abs(vf);
    if (tmp<vf) vf=tmp;
    *value -= vf*curr_diom->P[0];
    vf_all -= vf;
   }
  }
  curr_diom = curr_diom->next;
 }

/*
  Exterior.
*/
  *ival = 0;
/*
  Interior.
*/
  if ( vf_all > 0.999999 ) {
    *ival = 1;
/*
  Boundary
*/
  } else if ( vf_all > 0.0 ) {
    *value /= vf_all;
    *ival = -1;
  }

  return;
}
