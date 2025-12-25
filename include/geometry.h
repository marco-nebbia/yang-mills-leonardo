#ifndef GEOMETRY_H
#define GEOMETRY_H

#include"macro.h"
#include<stdbool.h>  // For boolean data type (bool, true, false)

typedef struct Geometry {
   int d_size[STDIM];

   long d_volume;           // total volume
   double d_inv_vol;        // 1 / tot. volume
   long d_space_vol;        // spatial component of the volume
   double d_inv_space_vol;  // 1 / spatial component of the volume

   long **d_nnp;      // d_nnp[r][i] = next neighbour (on the local lattice) in dir.  i of the site r
   long **d_nnm;      // d_nnm[r][i] = next neighbour (on the local lattice) in dir. -i of the site r
   int *d_timeslice;  // d_timeslice[r]  = time component of r
   long *d_spacecomp; // d_spacecomp[r]  = space component of r
   long **d_tsp;      // d_tsp[t][rsp] = r such that d_timeslice[r]=t and d_spacecomp[r]=rsp
} Geometry;

// these are the functions to be used in switching between different indices
extern long (*cart_to_si)(int const * const cartcoord, Geometry const * const geo); // cartesian coordinates -> single index
extern void (*si_to_cart)(int *cartcoord, long si, Geometry const * const geo);     // single index -> cartesian coordinates
extern long (*lex_to_si)(long lex, Geometry const * const geo);          // lexicographic -> single index
extern long (*si_to_lex)(long si, Geometry const * const geo);           // lexicographic -> single index
extern long (*sisp_and_t_to_si_compute)(long sisp, int t, Geometry const * const geo);            // single index spatial and time -> single index tot
extern void (*si_to_sisp_and_t_compute)(long *sisp, int *t, long si, Geometry const * const geo); // single index tot -> single index spatial and time

// general functions
void init_geometry(Geometry *geo, int insize[STDIM]);
void free_geometry(Geometry *geo);

// check if link is on border face of STDIM=1
bool check_link_on_border(Geometry const * const geo, long lexeo, long j); 

// next neighbour in + direction
inline long nnp(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnp[r][i];
  }

// next neighbour in - direction
inline long nnm(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnm[r][i];
  }

// single index spatial and time -> single index tot
inline long sisp_and_t_to_si(Geometry const * const geo, long sisp, int t)
  {
  return geo->d_tsp[t][sisp];
  }

// single index tot -> single index spatial and time
inline void si_to_sisp_and_t(long *sisp, int *t, Geometry const * const geo, long si)
  {
  *sisp=geo->d_spacecomp[si];
  *t=geo->d_timeslice[si];
  }

// for debug
void test_geometry(Geometry const * const geo);

//------------ these are not to be used outside geometry.c ----------------

long cart_to_lex(int const * const cartcoord, Geometry const * const geo);   // cartesian coordinates -> lexicographic index
void lex_to_cart(int *cartcoord, long lex, Geometry const * const geo);      // lexicographic index -> cartesian coordinates

long cart_to_lexeo(int const * const cartcoord, Geometry const * const geo); // cartesian coordinates -> lexicographic eo index
void lexeo_to_cart(int *cartcoord, long lexeo, Geometry const * const geo);  // lexicographic eo index -> cartesian coordinates

long lex_to_lexeo(long lex, Geometry const * const geo);                     //  lexicographic index -> lexicographic eo index
long lexeo_to_lex(long lexeo, Geometry const * const geo);                   //  lexicographic eo index -> lexicographic index

long cartsp_to_lexsp(int const * const ccsp, Geometry const * const geo); // spatial cartesian coordinates -> spatial lexicographic index
void lexsp_to_cartsp(int *ccsp, long lexsp, Geometry const * const geo);  // spatial lexicographic index -> spatial cartesian coordinates

long cartsp_to_lexeosp(int const * const ccsp, Geometry const * const geo);  // spatial cartesian coordinates -> spatial lexicographic eo index
void lexeosp_to_cartsp(int *ccsp, long lexeosp, Geometry const * const geo); // spatial lexicographic eo index -> spatial cartesian coordinates

long lexsp_to_lexeosp(long lexsp, Geometry const * const geo);     //  spatial lexicographic index -> spatial lexicographic eo index
long lexeosp_to_lexsp(long lexeosp, Geometry const * const geo);   //  spatial lexicographic eo index -> spatial lexicographic index

long lexeosp_and_t_to_lexeo(long lexeosp, int t, Geometry const * const geo);    // lexicographic eo spatial and time -> lexicographic eo index
void lexeo_to_lexeosp_and_t(long *lexeosp, int *t, long lexeo, Geometry const * const geo); // lex. eo index -> lex. eo spatial and t

long cartorth_to_lexorth(int const * const ccorth, Geometry const * const geo); // spatial cartesian coordinates -> spatial lexicographic index
void lexorth_to_cartorth(int *ccorth, long lexorth, Geometry const * const geo);  // spatial lexicographic index -> spatial cartesian coordinates

long cartorth_to_lexeoorth(int const * const ccorth, Geometry const * const geo);  // spatial cartesian coordinates -> spatial lexicographic eo index
void lexeoorth_to_cartorth(int *ccorth, long lexeoorth, Geometry const * const geo); // spatial lexicographic eo index -> spatial cartesian coordinates

long lexorth_to_lexeoorth(long lexorth, Geometry const * const geo);     //  spatial lexicographic index -> spatial lexicographic eo index
long lexeoorth_to_lexorth(long lexeoorth, Geometry const * const geo);   //  spatial lexicographic eo index -> spatial lexicographic index

long lexeoorth_and_dir_to_lexeosp(long lexeosp, int dir, Geometry const * const geo);
void lexeosp_to_lexeoorth_and_dir(long *lexeoorth, int *dir, long lexeosp, Geometry const * const geo);

#endif
