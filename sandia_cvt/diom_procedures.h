/* $Id: diom_procedures.h,v 1.1 2000/03/23 18:33:52 dacrawf Exp $ */

/* diom_procedures.h */

REAL bit_volume(REAL *t0, REAL *t1, REAL *t2, REAL *t3, REAL *bmin, REAL *bmax);
char * diom_alloc(size_t num_bytes);
char * diom_realloc(void * bptr, size_t old_num_bytes, size_t num_bytes);
void diom_memory_diagnostics(INTEGER *max_bytes, INTEGER *num_bytes, INTEGER *err);
struct geom_object *geom_add(struct geom_object *dest, char *shape, INTEGER obj_type,
               INTEGER aux_type, INTEGER VFadd, INTEGER iter,
               void (*GINIT)(REAL *G,INTEGER *NG, REAL *XM),
               void (*GSUB)(REAL *G,INTEGER *NG,REAL *X,REAL *Y,REAL *Z,
                      REAL *DX,REAL *DY,REAL *DZ,REAL *XIO),
               void (*GTRANS)(REAL *G,INTEGER *NG,REAL *DX,REAL *DY,REAL *DZ),
               void (*GXFRM)(REAL *G,INTEGER *NG,REAL *TMATRX) );
void diom_free(void * bptr, size_t num_bytes);
void diom_init(struct diom *dest, char *Adiom_id);
struct diom * diom_find(char *Adiom_id);
REAL * diom_get_Lagrangian_velocity();
void diom_object(char *shape,INTEGER obtyp,INTEGER aux_type,INTEGER vfadd,INTEGER *ngeom);
void diom_get_shape(char *shape_name, INTEGER *type, INTEGER *aux_type);
void diom_box_to_triangle(REAL *geometry);
void diom_box_to_tet(REAL *geometry);
INTEGER XCoord_LCoord(struct diom *tmp, REAL x);
REAL get_VF3R(struct geom_object *vf_geom, REAL x0, REAL y0,
              REAL z0, REAL x1, REAL y1, REAL z1, INTEGER n);
REAL get_VF2R(struct geom_object *vf_geom, REAL x0, REAL y0,
              REAL z0, REAL x1, REAL y1, REAL z1, INTEGER n);
REAL get_VF2C(struct geom_object *vf_geom, REAL x0, REAL y0,
              REAL z0, REAL x1, REAL y1, REAL z1, INTEGER n);
REAL get_VF1R(struct geom_object *vf_geom, REAL x0, REAL y0,
              REAL z0, REAL x1, REAL y1, REAL z1, INTEGER n);
REAL get_VF1C(struct geom_object *vf_geom, REAL x0, REAL y0,
              REAL z0, REAL x1, REAL y1, REAL z1, INTEGER n);
REAL get_VF1S(struct geom_object *vf_geom, REAL x0, REAL y0,
              REAL z0, REAL x1, REAL y1, REAL z1, INTEGER n);
void add_ureg(REAL *AUR);
INTEGER in_ureg(INTEGER IGM, REAL x0, REAL y0, REAL z0, REAL x1, REAL y1, REAL z1);
void gstable_init(struct gstable *astable, REAL time, REAL dt, INTEGER init, INTEGER *err);
INTEGER nest_push(INTEGER IVAL);
INTEGER nest_pop();
INTEGER nest_value();
REAL diom_min(REAL x, REAL y);
REAL diom_max(REAL x, REAL y);
void diom_matrix_mult(REAL *amatrix, REAL *bmatrix);
void diom_matrix_vector_mult(REAL *dvector, REAL *matrix, REAL *svector);
void diom_grading(INTEGER Agrading, REAL *AP1, REAL *AP2, INTEGER *err);
void diom_table_function(void (*ATABSUB)(REAL *XOUT,INTEGER *NTAB,REAL *XIN));
void diom_tet(REAL *p);
void diom_box_to_tet(REAL *p);
void diom_box_to_triangle(REAL *p);
void diom_duplicate(char *Adiom_id);
void diom_get_id(char *Adiom_id);
void diom_pop();
void diom_add();
void diom_remove();
void diom_first();
void diom_next();
void diom_set_property_arrays(REAL *AP, INTEGER *ATP, INTEGER ANP);
void diom_get_property_arrays(REAL *AP, INTEGER *ATP, INTEGER *ANP);
void diom_set_int_property_array(INTEGER *AIP, INTEGER ANIP);
void diom_get_int_property_array(INTEGER *AIP, INTEGER *ANIP);
void diom_set_iteration_depth(INTEGER Aiter);
void diom_on(REAL Aon);
void diom_off(REAL Aoff);
void diom_dt(REAL Adt);
void diom_grading(INTEGER Agrading, REAL AP1[3], REAL AP2[3], INTEGER *err);
void diom_get_grading(INTEGER *grading, REAL P0[3], REAL A0[3]);
void diom_anywhere(INTEGER Aanywhere);
void diom_Lagrangian_velocity(INTEGER i, REAL AV);
void diom_set_eos_tag(INTEGER Aeostag);
void diom_clean();
void diom_init(struct diom *dest, char *Adiom_id);
void diom_cp(struct diom *dest, struct diom *src, char *dest_id);
void diom_delete(struct diom *tmp);
void diom_set_geometry(REAL *AG, INTEGER ANG);
void diom_get_geometry(REAL *AG, INTEGER *ANG);
void diom_set_geometry_filename(char *filename);
void diom_get_geometry_filename(struct geom_object *curr_geom, char *filename);
void diom_set_geometry_mat_id(INTEGER Aid);
void diom_log_geometry(struct diom *adiom);
INTEGER diom_get_geometry_mat_id(struct geom_object *curr_geom);
void diom_set_center(REAL *Acg);
void diom_set_avel(REAL *Aavel);
void diom_set_mvel(REAL *Amvel);
void diom_get_shape(char *shape_name, INTEGER *type, INTEGER *aux_type);
void diom_free(void * bptr, size_t num_bytes);
void diom_object(char *shape,INTEGER obtyp,INTEGER aux_type,INTEGER vfadd,INTEGER *ngeom);
REAL diom_get_on();
REAL diom_get_off();
REAL diom_min(REAL x, REAL y);
REAL diom_max(REAL x, REAL y);
INTEGER diom_get_iteration_depth();
INTEGER diom_get_eos_tag();
REAL diatom_volume_fraction(INTEGER init, INTEGER start_flag, INTEGER *stop_flag,
           REAL x1, REAL y1, REAL z1, REAL x2, REAL y2, REAL z2,
           REAL time, REAL dt, REAL *rp, INTEGER *ip);
void diatom_per_cycle_init(INTEGER init, REAL time, REAL dt, REAL *x1, REAL *y1, REAL *z1,
                           REAL *x2, REAL *y2, REAL *z2, INTEGER *err);
void transform_init(struct transform *tmp, INTEGER count);
void transform_delete(struct transform *tmp);
void transform_pop();
void translate_push(REAL *Avector, INTEGER count);
void scale_push(REAL Ascalar, INTEGER count);
void rotate_push(INTEGER dir, REAL Angle, INTEGER count);
/*void vector_transform(struct transform *tmp1, REAL dvector[3], REAL svector[3]);*/
void vector_rotate_stretch(struct transform *tmp1, REAL dvector[3], REAL svector[3]);
void scalar_stretch(struct transform *tmp1, REAL *dscalar, REAL scalar);
struct transform * collapse_transform();
void construct_arbitrary_rotation(REAL vector[3], REAL *amatrix);
void rotate_arbitrary_push(REAL vector[3], INTEGER count);
struct geom_object *geom_add(struct geom_object *list, char *shape, INTEGER obj_type,
               INTEGER aux_type, INTEGER VFadd, INTEGER iter,
               void (*GINIT)(REAL *G,INTEGER *NG, REAL *XM),
               void (*GSUB)(REAL *G,INTEGER *NG,REAL *X,REAL *Y,REAL *Z,
                      REAL *DX,REAL *DY,REAL *DZ,REAL *XIO),
               void (*GTRANS)(REAL *G,INTEGER *NG,REAL *DX,REAL *DY,REAL *DZ),
               void (*GXFRM)(REAL *G,INTEGER *NG,REAL *TMATRX) );
void geom_list_delete(struct geom_object *tmp);
void gstable_push(INTEGER TYPE, REAL AGRAV, REAL X0[3],
            INTEGER NL, INTEGER NMAT, void (*EOS)());
void gstable_pop();
void gstable_init(struct gstable *astable, REAL time, REAL dt, INTEGER init, INTEGER *err);
void get_grav_DT(struct gstable *tmp, REAL xc, REAL yc, REAL zc,
                 INTEGER mat, REAL *density, REAL *temp);
void diom_hash_init(INTEGER nbins, INTEGER bin_length, INTEGER length);
void diom_hash_add(INTEGER id, REAL *vector);
void diom_hash_delete();
void diom_hash_get(INTEGER id, INTEGER *length, REAL *vector);
void diom_hash_count(INTEGER id, INTEGER *idc);
void diom_hash_total_count(INTEGER *count);
void diom_hash_set(INTEGER id, REAL *vector);
void diom_read_fnf_shell(char *fname, INTEGER mat, struct diom *curr_diom, 
  struct geom_object *curr_geom);
void diom_read_fnf_tet(char *fname, INTEGER mat, struct diom *curr_diom,
  struct geom_object *curr_geom);
void diom_read_topo(char *fname, struct geom_object *curr_geom);
void broadcast_integer(INTEGER *int_val);
void broadcast_real(REAL *real_val);
void broadcast_real_array(INTEGER *num_vals, REAL *real_val);
void broadcast_min(REAL *src, REAL *dest);
void diom_push_do(long int fpos, char *instr, REAL start, REAL end, REAL delta);
void diom_pop_do();
INTEGER diom_parse_do(char *instr);
INTEGER diom_end_do(long int *fpos);
void diom_set_var(char *var_name, REAL var_value);
INTEGER diom_get_var(char *var_name);
REAL diom_expression(char *tmpstr);
