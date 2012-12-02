/* $Id: diom_globals.h,v 1.2 2000/11/02 14:53:31 dacrawf Exp $ */

/* diom_globals.h */

#include <stddef.h>

/* The following ifdef's were added for compatibility with ALEGRA */

#ifdef LC_FLINK
#define UNDERSCORE 1
#endif

#ifdef UCFLINK
#define UPCASE 1
#endif

#ifdef LCFLINK
#define LOWCASE 1
#endif

#ifdef UPCASE
#include <fortran.h>
#endif

#define INTEGER int
#define REAL double

/************ Structures *************/


struct transform {
 REAL tvector[3];
 REAL tmatrix[9];
 INTEGER ncopies;
 INTEGER xtype;
 struct transform *next,*previous;
};

struct geom_object {
 char shape_name[20];  /* name of geometry...useful for diagnostics */
 INTEGER NG;           /* length of geometry array */
 REAL *G;              /* geometry array */
 INTEGER TYPE;         /* object type */
 INTEGER AUX_TYPE;     /* aux object type */
 INTEGER VFADD;        /* +1 = add, -1 = delete */
 INTEGER iter;         /* iteration level */
 REAL minmax[6];
 char *filename;       /* optional filename */ 
 INTEGER mat_id;       /* optional material id redirect */
 void (*GINIT)(REAL *G,INTEGER *NG, REAL *XM);
 void (*GSUB)(REAL *G,INTEGER *NG,REAL *X,REAL *Y,REAL *Z,
           REAL *DX,REAL *DY,REAL *DZ,REAL *XIO);
 void (*GTRANS)(REAL *G,INTEGER *NG,REAL *DX,REAL *DY,REAL *DZ);
 void (*GXFRM)(REAL *G,INTEGER *NG,REAL *TMATRX);
 struct geom_object *next;
};

#define GCONST = 6.67e-8; /* cgs units */

struct gstable {
 INTEGER initialized;
 INTEGER Type;
 REAL x0[3];
 REAL Gravity;
 INTEGER NL;
 REAL LMIN;
 REAL DL;
 INTEGER NMAT;
 REAL **GSD,*GST,*GSP;
 void (*EOSSUB)(INTEGER *,REAL *,REAL *,REAL *,REAL *, REAL *, REAL *, INTEGER *,
                REAL *, INTEGER *);
 struct gstable *next, *previous;
};

struct diom {
 char *diom_id;
 struct geom_object *geom_list,*curr_geom;
 struct transform *transformation;  /* transformation of this package */
 REAL *P;             /* Property array */
 INTEGER *TP;         /* Property table array (optional) */
 INTEGER *IP;         /* Integer property array 0 = mat i.d., 1 = package type*/
 INTEGER NP;          /* Length of real property array */
 INTEGER NIP;         /* length of integer property array */
 INTEGER iter;        /* iteration depth for VF calculations */
 INTEGER grading;     /* type of grading */
 REAL grading_P0[3];
 REAL grading_A0[3];
 REAL center[3];
 REAL avel[3];
 REAL mvel[3];
 struct gstable *gstable_ptr;
 INTEGER anywhere;    /* 1 = diom can be placed anywhere in mesh */
                      /*     regardless of update regions */
                      /* 2 = same as above except only on cycle 0 */
 INTEGER eostag;      /* 0 = don't use CTH-like EOS calculations in diom lib */
                      /* 1 = pressure dependent EOS */
                      /* 2 = temperature dependent EOS */
 REAL minmax[6];
 REAL R0[3];
 REAL V0[3];
 INTEGER speed, spin_flag, center_flag, mvel_flag;
 REAL displacement, rotation;
 REAL on, off, dt, t0;
 struct diom *next,*previous;
};

struct update_region {
 REAL UR[12];
 struct update_region *next;
};

struct diom_vector {
 REAL x,y,z;
};

/************Globals************/

extern INTEGER global_coords;       /* 1dr, 2dc, etc. */
extern INTEGER number_dimensions;   /* 1, 2 or 3 */
extern INTEGER mesh_type;           /* 0 = brick, ... */
extern REAL diom_mesh_min[3];
extern REAL diom_mesh_max[3];

#define NTPROP 300   /* This value is MAXISV+100 */
#define NTIPRP 300   /* This value is MAXISV+100 */
#define NTGEOM 1000

extern struct diom *top_diom, *stack_diom,
                   *head_diom, *curr_diom;

extern struct transform *top_transform, *stack_transform;

extern struct update_region *top_update, *curr_update;

extern struct gstable *top_grav, *stack_grav, *match_grav;

extern size_t total_num_bytes;
extern size_t max_num_bytes;
extern INTEGER memory_alloc_error;

/* malloc global */

extern void * (*diom_malloc)(size_t num_bytes);
extern void (*diom_free_mem)(void *);

/* Parallel processing globals */

extern INTEGER number_of_processors;
extern INTEGER this_processor_id;
extern INTEGER master;
extern void (*broadcast_int_array)(INTEGER *num_vals, INTEGER *int_array);
extern void (*broadcast_max)(REAL *src, REAL *dest);

/* Message processing globals */

extern void (*diom_log_message)(char *c_array);
extern void (*diom_fatal_error)(char *c_array);
extern void (*diom_bombed)(char *c_array);

/* EOS globals */

extern void (*diom_temp_dep_EOS)(INTEGER *,REAL *,REAL *,REAL *,REAL *, REAL *, REAL *, INTEGER *, REAL *, INTEGER *);
extern void (*diom_pres_dep_EOS)(INTEGER *,REAL *,REAL *,REAL *,REAL *, REAL *, REAL *, INTEGER *, REAL *, INTEGER *);
extern REAL *diom_ref_density;
extern REAL *diom_ref_temperature;
extern INTEGER diom_number_mats;
extern INTEGER *diom_mat_id;

/* Volume fraction routine globals. */

extern REAL (*diom_brick_VF)(struct geom_object *vf_geom, REAL x0, REAL y0,
                             REAL z0, REAL x1, REAL y1, REAL z1, INTEGER n);

/* Table subroutine global */

extern void (*diom_table_sub)(REAL *XOUT,INTEGER *NTAB,REAL *XIN);

/* Gravity global */

extern REAL diom_gravity;
