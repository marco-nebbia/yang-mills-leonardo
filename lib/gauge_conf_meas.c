#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

// 0 <-> no debug
// 1 <-> Polyakov loop debug
// 2 <-> Polyakov loop debug
// 3 <-> Polyakov correlator debug
#ifndef GAUGE_DEBUG
#define GAUGE_DEBUG 0
#endif

#ifndef DEBUG
#define DEBUG
#endif

#include"../include/macro.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/flavour_matrix.h"
#include"../include/gparam.h"
#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"
#include"../include/su2_monopoles.h"
#include"../include/sun_monopoles.h"


// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  long r,
                  int i,
                  int j)
   {
   GAUGE_GROUP matrix;

   #ifdef DEBUG
   if(r >= geo->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(&matrix, &(GC->lattice[r][i]));
   times_equal(&matrix, &(GC->lattice[r][j]));

   return retr(&matrix);
   }


// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double complex plaquettep_complex(Gauge_Conf const * const GC,
                                  Geometry const * const geo,
                                  long r,
                                  int i,
                                  int j)
   {
   GAUGE_GROUP matrix;

   #ifdef DEBUG
   if(r >= geo->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(&matrix, &(GC->lattice[r][i]));
   times_equal(&matrix, &(GC->lattice[r][j]));

   return retr(&matrix)+I*imtr(&matrix);
   }


// computation of the plaquette (matrix) in position r and positive directions i,j
void plaquettep_matrix(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       long r,
                       int i,
                       int j,
                       GAUGE_GROUP *matrix)
   {
   #ifdef DEBUG
   if(r >= geo->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

//
//       ^ j
//       |   (3)
//       +---<---+
//       |       |
//   (4) V       ^ (2)
//       |       |
//       +--->---+---> i
//       r   (1)
//

   equal(matrix, &(GC->lattice[r][i]));
   times_equal(matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(matrix, &(GC->lattice[r][j]));
   }


// compute the four-leaf clover in position r, in the plane i,j and save it in M
void clover(Gauge_Conf const * const GC,
            Geometry const * const geo,
            long r,
            int i,
            int j,
            GAUGE_GROUP *M)
   {
   GAUGE_GROUP aux;
   long k, p;

   #ifdef DEBUG
   if(r >= geo->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i >= STDIM || j >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

   zero(M);

//
//                   i ^
//                     |
//              (7)    |     (2)
//         +-----<-----++----->-----+
//         |           ||           |
//         |           ||           |
//    (6)  ^       (8) V^ (1)       V (3)
//         |           ||           |
//         |   (5)     || r   (4)   |
//       k +-----<-----++-----<-----+------>   j
//         +----->-----++----->-----+
//         |    (12)   ||   (13)    |
//         |           ||           |
//    (11) ^       (9) V^ (16)       V (14)
//         |           ||           |
//         |           ||           |
//         +------<----++-----<-----+
//              (10)   p      (15)
//
   // avanti-avanti
   equal(&aux, &(GC->lattice[r][i]) );                           // 1
   times_equal(&aux, &(GC->lattice[nnp(geo, r, i)][j]) );        // 2
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, r, j)][i]) );    // 3
   times_equal_dag(&aux, &(GC->lattice[r][j]) );                 // 4
   plus_equal(M, &aux);

   k=nnm(geo, r, j);

   // avanti-indietro
   equal_dag(&aux, &(GC->lattice[k][j]) );                       // 5
   times_equal(&aux, &(GC->lattice[k][i]) );                     // 6
   times_equal(&aux, &(GC->lattice[nnp(geo, k, i)][j]) );        // 7
   times_equal_dag(&aux, &(GC->lattice[r][i]) );                 // 8
   plus_equal(M, &aux);

   p=nnm(geo, r, i);

   // indietro-indietro
   equal_dag(&aux, &(GC->lattice[p][i]) );                       // 9
   times_equal_dag(&aux, &(GC->lattice[nnm(geo, k, i)][j]) );    // 10
   times_equal(&aux, &(GC->lattice[nnm(geo, k, i)][i]) );        // 11
   times_equal(&aux, &(GC->lattice[k][j]) );                     // 12
   plus_equal(M, &aux);

   // indietro-avanti
   equal(&aux, &(GC->lattice[r][j]) );                            // 13
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, p, j)][i]) );     // 14
   times_equal_dag(&aux, &(GC->lattice[p][j]) );                  // 15
   times_equal(&aux, &(GC->lattice[p][i]) );                      // 16
   plus_equal(M, &aux);
   }


// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               double *plaqs,
               double *plaqt)
   {
   long r;
   double ps, pt;

   ps=0.0;
   pt=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt) reduction(+ : ps)
   #endif
   for(r=0; r<(geo->d_volume); r++)
      {
      int i, j;
      i=0;
      for(j=1; j<STDIM; j++)
         {
         pt+=plaquettep(GC, geo, r, i, j);
         }
     
      for(i=1; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            ps+=plaquettep(GC, geo, r, i, j);
            }
         }
      }

   if(STDIM>2)
     {
     ps*=geo->d_inv_vol;
     ps/=((double) (STDIM-1)*(STDIM-2)/2);
     }
   else
     {
     ps=0.0;
     }

   pt*=geo->d_inv_vol;
   pt/=((double) STDIM-1);

   *plaqs=ps;
   *plaqt=pt;
   }

// compute the mean plaquettes with normal vector
// parallel or orthogonal to time direction for a
// specific slice of the lattice 
void plaquette_slice_time(Gauge_Conf const * const GC,
                          Geometry const * const geo,
                          int slice,
                          double *plaqs,
                          double *plaqt)
   {
   long r, rsp;
   double ps, pt;

   ps=0.;
   pt=0.;
   
   #ifdef DEBUG
   // check if site is inside the lattice
   if(slice>=(geo->d_size[0]))
      {
      fprintf(stderr, "Slice outside lattice: %d >= %d (%s, %d)\n", slice, geo->d_size[0], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   else if(slice<0)
      {
      fprintf(stderr, "Slice outside lattice: %d < 0 (%s, %d)\n", slice, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   #endif

   int i,j;

   for(rsp=0; rsp<(geo->d_space_vol); rsp++)
      {
      r=sisp_and_t_to_si(geo, rsp, slice);
      
      i=0;
      for(j=1; j<STDIM; j++)
         {
         pt+=plaquettep(GC, geo, r, i, j);
         }
     
      for(i=1; i<STDIM; i++)
         {
         j=i+1;

         while(j<STDIM)
            {
            ps+=plaquettep(GC, geo, r, i, j);
            j++;
            }
         }
      }
   
   if(STDIM>2)
     {
     ps*=geo->d_inv_space_vol;
     ps/=((double) (STDIM-1)*(STDIM-2)/2);
     }
   else
     {
     ps=0.0;
     }

   pt*=geo->d_inv_space_vol;
   pt/=((double) STDIM-1);

   *plaqs=ps;
   *plaqt=pt;
   }

// compute the clover discretization of
// sum_{\mu\nu}  Tr(F_{\mu\nu}F_{\mu\nu})/2
void clover_disc_energy(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        double *energy)
  {
  long r;
  double ris;

  ris=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
  #endif
  for(r=0; r<geo->d_volume; r++)
     {
     int i, j;
     GAUGE_GROUP aux1, aux2;

     for(i=0; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           clover(GC, geo, r, i, j, &aux1);

           ta(&aux1);
           equal(&aux2, &aux1);
           times_equal(&aux1, &aux2);
           ris+=-NCOLOR*retr(&aux1)/16.0;
           }
        }
     }

  *energy=ris*geo->d_inv_vol;
  }


// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              double *repoly,
              double *impoly)
   {
   long rsp;
   double rep, imp;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<geo->d_space_vol; rsp++)
      {
      long r;
      int i;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(i=0; i<geo->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

      rep+=retr(&matrix);
      imp+=imtr(&matrix);
      }

   *repoly=rep*geo->d_inv_space_vol;
   *impoly=imp*geo->d_inv_space_vol;
   }


// compute the mean Polyakov loop in the adjoint representation (the trace of)
void polyakov_adj(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  double *repoly,
                  double *impoly)
   {
   long rsp;
   double rep, imp;
   double complex tr;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<geo->d_space_vol; rsp++)
      {
      long r;
      int i;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(i=0; i<geo->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }
      tr=NCOLOR*retr(&matrix)+NCOLOR*imtr(&matrix)*I;

      #if NCOLOR==1
        (void) tr;
        rep+=0.0;
      #else
        rep+=(cabs(tr)*cabs(tr)-1)/(NCOLOR*NCOLOR-1);
      #endif

      imp+=0.0;
      }

   *repoly=rep*geo->d_inv_space_vol;
   *impoly=imp*geo->d_inv_space_vol;
   }

// compute the Polyakov loop (the trace of) starting from the site r (its time column)
void polyakov_fixed_site(Gauge_Conf const * const GC,
                         Geometry const * const geo,
                         long r,
                         double *repoly,
                         double *impoly)
   {
   GAUGE_GROUP matrix;

   one(&matrix);

   long rsp;
   int t;

   #if (GAUGE_DEBUG == 1 || GAUGE_DEBUG == 2)
   int cartcoord[STDIM];
   #endif
   
   // rsp is single index for spatial coordinates,
   // t is time coordinate
   si_to_sisp_and_t(&rsp, &t, geo, r);

   #if GAUGE_DEBUG == 1
   si_to_cart(cartcoord, r, geo);
   printf("Single index coordinate: %ld\n", r);
   printf("Single index spatial coordinate: %ld\n", rsp);
   for(int i=0; i<STDIM; i++)
      {
      printf("Cartesian coordinate %d: %d\n", i, cartcoord[i]);
      }
   #endif

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif

   // be sure to start from bottom of the lattice in time dimension
   r=sisp_and_t_to_si(geo, rsp, 0);


   #if GAUGE_DEBUG == 1
   si_to_cart(cartcoord, r, geo);
   printf("After changing time coordinate:\n");
   printf("Single index coordinate: %ld\n", r);
   printf("Single index spatial coordinate: %ld\n", rsp);
   for(int i=0; i<STDIM; i++)
      {
      printf("Cartesian coordinate %d: %d\n", i, cartcoord[i]);
      }
   #endif
   
   int x0=0;

   // last link is outside lattice if Dirichlet BC are imposed
   // on both bottom and top faces, so it doesn't take part in
   // Polyakov loop, so if DIRICHLET_MODE=1 then x0 stops at
   // geo->d_size[0]-1
   for(x0=0; x0<geo->d_size[0]; x0++)
      {
      #if DIRICHLET_MODE == 1
      if(x0==geo->d_size[0]-1) break;
      #endif
      times_equal(&matrix, &(GC->lattice[r][0]));
      
      #if GAUGE_DEBUG == 2
      si_to_cart(cartcoord, r, geo);
      printf("r value: %ld\n", r);
      for(int i=0; i<STDIM; i++)
         {
         printf("x_%d value: %d\n", i, cartcoord[i]);
         }
      printf("Gauge variable real value: %f\n", (GC->lattice[r][0]).comp[0]);
      printf("Gauge variable imaginary values: %f  %f  %f\n", (GC->lattice[r][0]).comp[1], 
                                                              (GC->lattice[r][0]).comp[2], 
                                                              (GC->lattice[r][0]).comp[3]);
      printf("Matrix real value: %f\n", matrix.comp[0]);
      printf("Matrix imaginary values: %f  %f  %f\n", matrix.comp[1], 
                                                      matrix.comp[2], 
                                                      matrix.comp[3]);
      printf("\n");
      #endif
     
      r=nnp(geo, r, 0);
      }

   if(x0!=geo->d_size[0]-DIRICHLET_MODE)
      {
      fprintf(stderr, "Links multiplied: %d (should be %d)\n", 
            x0, geo->d_size[0] - DIRICHLET_MODE);
      exit(EXIT_FAILURE);
      }
   
   *repoly=retr(&matrix);
   *impoly=imtr(&matrix);

   #if GAUGE_DEBUG == 2
   printf("Repoly value: %f\n", *repoly);
   printf("Impoly value: %f\n", *impoly);
   printf("\n");
   #endif
   }

void polyakov_horizontal_fixed_site(Gauge_Conf const * const GC, 
                                    Geometry const * const geo,
                                    long r,
                                    int dir,
                                    double *repoly,
                                    double *impoly)
   {
   #ifdef DEBUG
   if(dir>STDIM)
      {
      fprintf(stderr, "Direction of horizontal Polyakov loop is greater than space-time dimensions! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   else if(dir==0)
      {
      fprintf(stderr, "Direction of horizontal Polyakov loop is time direction! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }   
   #endif

   // initialize gauge matrix and coordinates
   GAUGE_GROUP matrix;
   int cartcoord[STDIM];
   
   one(&matrix);
   si_to_cart(cartcoord, r, geo);

   // start from the 0th link in direction dir
   cartcoord[dir]=0;

   r=cart_to_si(cartcoord, geo);

   for(int x_dir=0; x_dir<(geo->d_size[dir]); x_dir++)
      {
      times_equal(&matrix, &(GC->lattice[r][dir]));
      
      #if GAUGE_DEBUG == 2
      si_to_cart(cartcoord, r, geo);
      printf("r value: %ld\n", r);
      for(int i=0; i<STDIM; i++)
         {
         printf("x_%d value: %d\n", i, cartcoord[i]);
         }
      printf("Gauge variable real value: %f\n", (GC->lattice[r][dir]).comp[0]);
      printf("Gauge variable imaginary values: %f  %f  %f\n", (GC->lattice[r][0]).comp[1], 
                                                              (GC->lattice[r][0]).comp[2], 
                                                              (GC->lattice[r][0]).comp[3]);
      printf("Matrix real value: %f\n", matrix.comp[0]);
      printf("Matrix imaginary values: %f  %f  %f\n", matrix.comp[1], 
                                                      matrix.comp[2], 
                                                      matrix.comp[3]);
      printf("\n");
      #endif
     
      r=nnp(geo, r, dir);
      }
   
   *repoly=retr(&matrix);
   *impoly=imtr(&matrix);

   #if GAUGE_DEBUG == 2
   printf("Repoly value: %f\n", *repoly);
   printf("Impoly value: %f\n", *impoly);
   printf("\n");
   #endif
      
   }
                     
// compute the mean Polyakov loop and its powers (trace of) in the presence of trace deformation
void polyakov_for_tracedef(Gauge_Conf const * const GC,
                           Geometry const * const geo,
                           double *repoly,
                           double *impoly)
   {
   long rsp;
   double **rep, **imp;
   int j, err;
   long i;

   for(j=0;j<(int)floor(NCOLOR/2);j++)
      {
      repoly[j]=0.0;
      impoly[j]=0.0;
      }

   err=posix_memalign((void**)&rep, (size_t)DOUBLE_ALIGN, (size_t) geo->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   err=posix_memalign((void**)&imp, (size_t)DOUBLE_ALIGN, (size_t) geo->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   for(i=0; i<geo->d_space_vol; i++)
      {
      err=posix_memalign((void**)&(rep[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      err=posix_memalign((void**)&(imp[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   for(i=0; i<geo->d_space_vol; i++)
      {
      for(j=0; j<(int)floor(NCOLOR/2); j++)
         {
         rep[i][j] = 0.0;
         imp[i][j] = 0.0;
         }
      }

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp)
   #endif
   for(rsp=0; rsp<geo->d_space_vol; rsp++)
      {
      long r;
      int k;
      GAUGE_GROUP matrix, matrix2;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(k=0; k<geo->d_size[0]; k++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

       rep[rsp][0] = retr(&matrix);
       imp[rsp][0] = imtr(&matrix);

       equal(&matrix2, &matrix);

      for(k=1; k<(int)floor(NCOLOR/2.0); k++)
         {
         times_equal(&matrix2, &matrix);
         rep[rsp][k] = retr(&matrix2);
         imp[rsp][k] = imtr(&matrix2);
         }
      }

    for(j=0; j<(int)floor(NCOLOR/2); j++)
       {
       for(i=0; i<geo->d_space_vol; i++)
          {
          repoly[j] += rep[i][j];
          impoly[j] += imp[i][j];
          }
       }

   for(j=0; j<(int)floor(NCOLOR/2.0); j++)
      {
      repoly[j] *= geo->d_inv_space_vol;
      impoly[j] *= geo->d_inv_space_vol;
      }

   for(i=0; i<geo->d_space_vol; i++)
      {
      free(rep[i]);
      free(imp[i]);
      }
   free(rep);
   free(imp);
   }

// product of two Polyakov loops in spatial direction dir
// averaged on all spatial sites
void polyakov_correlator_dir(Gauge_Conf const * const GC,
                             Geometry const * const geo,
                             int d,
                             int dir,
                             double *re,
                             double *im)
   {
   double repoly1, repoly2;
   double impoly1, impoly2;
   double reprod, improd;

   int cartcoord[STDIM];
   long r, rsp;

   #ifdef DEBUG
   if(dir==0)
      {
      fprintf(stderr, "Direction of correlator of Polyakov loops is time, should be a spatial dimension! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(dir>=STDIM)
      {
      fprintf(stderr, "Direction of correlator of Polyakov loops is greater than space-time dimensions! %d>=%d (%s, %d)\n", dir, STDIM, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(d<=0)
      {
      fprintf(stderr, "Distance between Polyakov loops <=0! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(d>geo->d_size[dir])
      {
      fprintf(stderr, "Distance between Polyakov loops greater than lattice size! %d>%d (%s, %d)\n", d, geo->d_size[dir], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   #endif

   reprod=0.;
   improd=0.;

   for(rsp=0; rsp<(geo->d_space_vol); rsp++)
      {
      r=sisp_and_t_to_si(geo, rsp, 0);

      polyakov_fixed_site(GC, geo, r, &repoly1, &impoly1);

      #if GAUGE_DEBUG == 3
      fprintf(stderr, "TEST single Polyakov: re=%.12g im=%.12g\n", repoly1, impoly1);
      #endif

      si_to_cart(cartcoord, r, geo);
      cartcoord[dir]+=d;

      if(cartcoord[dir]>=(geo->d_size[dir])) cartcoord[dir]-=(geo->d_size[dir]);
      r=cart_to_si(cartcoord, geo);

      polyakov_fixed_site(GC, geo, r, &repoly2, &impoly2);

      #if GAUGE_DEBUG == 3
      fprintf(stderr, "TEST single Polyakov: re=%.12g im=%.12g\n", repoly2, impoly2);
      #endif

      #if NCOLOR == 2
      // if gauge group is SU(2), trace is real
      reprod+=repoly1*repoly2;
      improd=0.;
      #elif NCOLOR > 2
      // if gauge group is SU(N) with N!=2, trace is complex
      // so conjugate the 2nd loop to take it in downward sense
      // impoly2->-impoly2
      // and use complex multiplication rules
      reprod+=repoly1*repoly2+impoly1*impoly2;
      improd+=-repoly1*impoly2+repoly2*impoly1;
      #endif
      }

   *re=reprod*(geo->d_inv_space_vol);
   *im=improd*(geo->d_inv_space_vol);
   }

// product of two Polyakov loops parallel to
// the Dirichlet faces in the bulk of the lattice
/*void polyakov_product_in_bulk_dir2(Gauge_Conf const * const GC,
                              Geometry const * const geo,
                              int d,
                              double *re,
                              double *im)
   {
   double repoly1, repoly2;
   double impoly1, impoly2;
   double reprod, improd;
   
   #ifdef DEBUG
   if(d<=0)
      {
      fprintf(stderr, "Distance between Polyakov loops <=0\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   #endif

   // be sure to choose a number in the lattice
   int minsize=geo->d_size[2];

   #if (STDIM > 3)
   for(int i=3; i<STDIM; i++)
      {
      if(geo->d_size[i]<minsize) minsize=geo->d_size[i];
      }
   #endif

   int cartcoord[STDIM];
   long r;

   cartcoord[0]=0; // bottom of lattice
   // in the bulk, halfway between the
   // two Dirichlet faces
   cartcoord[1]=((geo->d_size[1])/2)-1;
   cartcoord[2]=0;

   for(int i=3; i<STDIM && i!=2; i++)
      {
      cartcoord[i]=0; // x_i=0
      }
   
   reprod=0.;
   improd=0.;

   // x_0 is direction of Polyakov loop
   // x_1 has Dirichlet boundary condition 
   // so no translation invariance
   #if STDIM == 3  
   for(int x2=0; x2<geo->d_size[2]; x2++)
      {
      cartcoord[2]=x2;

      r=cart_to_si(cartcoord, geo);
      
      // compute Polyakov loop from site r that has x_2
      polyakov_fixed_site(GC, geo, r, &repoly1, &impoly1);

      // moving to x_2+d
      cartcoord[2]+=d;
      if(cartcoord[2]>=geo->d_size[2]) cartcoord[2]-=geo->d_size[2];
      r=cart_to_si(cartcoord, geo);

      // compute Polyakov loop from site r' that has x_2=l+d
      polyakov_fixed_site(GC, geo, r, &repoly2, &impoly2);
      
      #if NCOLOR == 2
      // if gauge group is SU(2), trace is real
      reprod+=repoly1*repoly2;
      #elif
      // if gauge group is SU(N) with N!=2, trace is complex
      // so conjugate the 2nd loop to take it in downward sense
      // impoly2->-impoly2
      // and use complex multiplication rules
      reprod+=repoly1*repoly2+impoly1*impoly2;
      improd+=-repoly1*impoly2+repoly2*impoly1;
      #endif

      #if GAUGE_DEBUG == 3
      printf("Repoly1: %f\n", repoly1);
      printf("Repoly2: %f\n", repoly2);
      printf("Correlator: %f\n\n", reprod);
      #endif
      }
   
   reprod/=((double) geo->d_size[2]);

   #elif STDIM == 4
   for(int x2=0; x2<(geo->d_size[2]); x2++)
      {
      cartcoord[2]=x2;

      for(int x3=0; x3<(geo->d_size[3]); x3++)
         {
         cartcoord[3]=x3;
         r=cart_to_si(cartcoord, geo);

         // compute Polyakov loop from site r that has x_2=l
         polyakov_fixed_site(GC, geo, r, &repoly1, &impoly1);

         // moving to x_2=l+d
         cartcoord[2]+=d;
         if(cartcoord[2]>=geo->d_size[2]) cartcoord[2]-=geo->d_size[2];
         r=cart_to_si(cartcoord, geo);

         // compute Polyakov loop from site r' that has x_2=l+d
         polyakov_fixed_site(GC, geo, r, &repoly2, &impoly2);
         
         #if NCOLOR == 2
         // if gauge group is SU(2), trace is real
         reprod+=repoly1*repoly2;
         #elif
         // if gauge group is SU(N) with N!=2, trace is complex
         // so conjugate the 2nd loop to take it in downward sense
         // impoly2->-impoly2
         // and use complex multiplication rules
         reprod+=repoly1*repoly2+impoly1*impoly2;
         improd+=-repoly1*impoly2+repoly2*impoly1;
         #endif

         #if GAUGE_DEBUG == 3
         printf("Repoly1: %f\n", repoly1);
         printf("Repoly2: %f\n", repoly2);
         printf("Correlator: %f\n\n", reprod);
         #endif
         }
      }
   
   // average for number of sites
   reprod/=(((double) geo->d_size[2])*((double) geo->d_size[3]));
   improd/=(((double) geo->d_size[3])*((double) geo->d_size[3]));  
   #endif
   
   *re=reprod;
   *im=improd;
   }*/

// correlator with distance d in direction dircorr of two horizontal Polyakov loops
// in direction dirpoly averaged on all orthogonal directions to dirpoly
// on a timeslice slice
/*void polyakov_horizontal_correlator_timeslice_dir(Gauge_Conf const * const GC,
                                                  Geometry const * const geo,
                                                  int d,
                                                  int dirpoly,
                                                  int dircorr,
                                                  int slice,
                                                  double *re,
                                                  double *im)
   {
   double repoly1, repoly2;
   double impoly1, impoly2;
   double reprod, improd;
   
   int cartcoord[STDIM];
   long r, rsp, rsport, space_vol_ort;

   #ifdef DEBUG
   if(slice>=geo->d_size[0])
      {
      fprintf(stderr, "Timeslice to compute spatial Polyakov loops on is greater than time lattice size! %d>=%d (%s, %d)\n", slice, geo->d_size[0], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(dirpoly==0)
      {
      fprintf(stderr, "Direction of Polyakov spatial loop is chosen to be actually time dimension! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   else if(dirpoly==dircorr)
      {
      fprintf(stderr, "Direction of correlator distance is chosen to be the same direction of Polyakov spatial loops! %d=%d (%s, %d)\n", dirpoly, dircorr, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(d<=0)
      {
      fprintf(stderr, "Distance between spatial Polyakov loops <=0! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   else if(d>geo->d_size[dircorr])
      {
      fprintf(stderr, "Distance between spatial Polyakov loops greater than lattice size! %d>%d (%s, %d)\n", d, geo->d_size[dircorr], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(dircorr==0 && slice+d>geo->d_size[dircorr])
      {
      fprintf(stderr, "When computing spatial Polyakov loops correlators in time direction, it crosses the vertical face! %d+%d>%d (%s, %d)\n", d, slice, geo->d_size[dircorr], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   #endif

   space_vol_ort=(long) (geo->d_space_vol)/(geo->d_size[dirpoly]);

   reprod=0.;
   improd=0.;

   for(rsport=0; rsport<space_vol_ort; rsport++)
      {
      rsp=sisport_and_par_to_sisp(geo, rsport, 0, dirpoly);
      r=sisp_and_t_to_si(geo, rsp, slice);

      polyakov_horizontal_fixed_site(GC, geo, r, dirpoly, &repoly1, &impoly1);

      si_to_cart(cartcoord, r, geo);
      cartcoord[dircorr]+=d;
      if(cartcoord[dircorr]>=(geo->d_size[dir])) cartcoord[dircorr]-=(geo->d_size[dir]);

      r=cart_to_si(cartcoord, geo);

      polyakov_horizontal_fixed_site(GC, geo, r, dirpoly, &repoly2, &impoly2);

      #if NCOLOR == 2
      // if gauge group is SU(2), trace is real
      reprod+=repoly1*repoly2;
      improd=0.;
      #elif
      // if gauge group is SU(N) with N!=2, trace is complex
      // so conjugate the 2nd loop to take it in downward sense
      // impoly2->-impoly2
      // and use complex multiplication rules
      reprod+=repoly1*repoly2+impoly1*impoly2;
      improd+=-repoly1*impoly2+repoly2*impoly1;
      #endif
      }

   *re=reprod*(geo->d_size[dirpoly])*(geo->d_inv_space_vol);
   *im=improd*(geo->d_size[dirpoly])*(geo->d_inv_space_vol);
   }

// correlator with distance d in evey spatial direction of two horizontal Polyakov loops
// in direction dirpoly averaged on all spatial orthogonal directions to dirpoly
// on a timeslice slice
void polyakov_horizontal_correlator_timeslice(Gauge_Conf const * const GC,
                                              Geometry const * const geo,
                                              int d,
                                              int slice,
                                              double *re,
                                              double *im)
   {
   double repoly, impoly;
   double repoly_memo, impoly_memo;

   #ifdef DEBUG
   if(slice>=geo->d_size[0])
      {
      fprintf(stderr, "Timeslice to compute spatial Polyakov loops on is greater than time lattice size! %d>=%d (%s, %d)\n", slice, geo->d_size[0], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(d<=0)
      {
      fprintf(stderr, "Distance between spatial Polyakov loops <=0! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   #endif

   repoly_memo=0.;
   impoly_memo=0.;

   for(int i=1; i<STDIM; i++)
      {
      for(int j=1; j<STDIM; j++)
         {
         if(j==i) continue;
         
         polyakov_horizontal_correlator_timeslice_dir(GC, geo, d, i, j, slice, &repoly, &impoly);

         repoly_memo+=repoly;
         impoly_memo+=impoly;
         }
      }

   *re=(double) repoly/((STDIM-1)*(STDIM-2));
   *im=(double) impoly/((STDIM-1)*(STDIM-2));
   }*/


void polyakov_horizontal_corr_bulk_dir1(Gauge_Conf const * const GC,
                                        Geometry const * const geo,
                                        int d,
                                        double *re,
                                        double *im)
   {
   int cartcoord[STDIM], l;
   long r;
   double repoly1, impoly1;
   double repoly2, impoly2;
   double reprod, improd;

   #ifdef DEBUG
   if(d<=0)
      {
      fprintf(stderr, "Distance between spatial Polyakov loops <=0! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   else if(d>=(geo->d_size[0]))
      {
      fprintf(stderr, "Distance between spatial Polyakov loops greater than lattice time size! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE); 
      }
   #endif

   l=(int) (geo->d_size[0]-d)/2;

   cartcoord[0]=l-1;
   cartcoord[1]=0;
   cartcoord[2]=0;

   reprod=0.;
   improd=0.;

   for(int x2=0; x2<(geo->d_size[2]); x2++)
      {
      cartcoord[2]=x2;

      r=cart_to_si(cartcoord, geo);

      polyakov_horizontal_fixed_site(GC, geo, r, 1, &repoly1, &impoly1);

      cartcoord[2]+=d;
      if(cartcoord[2]>=(geo->d_size[2])) cartcoord[2]-=(geo->d_size[2]);

      r=cart_to_si(cartcoord, geo);

      polyakov_horizontal_fixed_site(GC, geo, r, 1, &repoly2, &impoly2);

      #if NCOLOR == 2
      // if gauge group is SU(2), trace is real
      reprod+=repoly1*repoly2;
      improd=0.;
      #elif
      // if gauge group is SU(N) with N!=2, trace is complex
      // so conjugate the 2nd loop to take it in downward sense
      // impoly2->-impoly2
      // and use complex multiplication rules
      reprod+=repoly1*repoly2+impoly1*impoly2;
      improd+=-repoly1*impoly2+repoly2*impoly1;
      #endif
      }

   *re=(double) reprod/(geo->d_size[2]);
   *im=(double) improd/(geo->d_size[2]);
   }


// compute the Wilson loop of dimensions d_i*d_j on directions (i,j)
// starting from lattice site r
void wilson_fixed_site(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       long r,
                       int d_i,
                       int d_j,
                       int i,
                       int j,
                       double *re,
                       double *im)
   {
   GAUGE_GROUP matrix;
   int x_i, x_j;

   #ifdef DEBUG
   if(i>=STDIM)
      {
      fprintf(stderr, "Spatial direction chosen for the Wilson loop greater than STDIM! %d>=%d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(j>=STDIM)
      {
      fprintf(stderr, "Spatial direction chosen for the Wilson loop greater than STDIM! %d>=%d (%s, %d)\n", j, STDIM, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(d_i>geo->d_size[i])
      {
      fprintf(stderr, "Wilson loop size exceeds corresponding lattice size! %d>%d (%s, %d)\n", d_i, geo->d_size[i], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if(d_j>geo->d_size[j])
      {
      fprintf(stderr, "Wilson loop size exceeds corresponding lattice size! %d>%d (%s, %d)\n", d_j, geo->d_size[j], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   #endif

   one(&matrix);

   // first horizontal product
   for(x_i=0; x_i<d_i; x_i++)
            {
            times_equal(&matrix, &(GC->lattice[r][i]));
            r=nnp(geo, r, i);
            }
        
   // first vertical product
   for(x_j=0; x_j<d_j; x_j++)
            {
            times_equal(&matrix, &(GC->lattice[r][j]));
            r=nnp(geo, r, j);
            }
         
   // second horizontal product, reverse order so multiply for the adjoint
   for(x_i=d_i-1; x_i>=0; x_i--)
            {
            // we move first this time because the SU(2) gauge variable
            // U[r][-i] corresponds to the adjoint of the previous U^{dag}[r-1][i]
            r=nnm(geo, r, i);
            times_equal_dag(&matrix, &(GC->lattice[r][i]));            
            }
        
   // second vertical product, reverse order so multiply for the adjoint
   for(x_j=d_j-1; x_j>=0; x_j--)
            {
            // we move first this time because the SU(2) gauge variable
            // U[r][-j] corresponds to the adjoint of the previous U^{dag}[r-1][j]
            r=nnm(geo, r, j);
            times_equal_dag(&matrix, &(GC->lattice[r][j]));
            }

   *re=retr(&matrix);
   *im=imtr(&matrix);
   }

// compute the mean Wilson loop of size R x T 
// where R is a distance in direction dir
// and T is a distance in time direction
// on a whole slice in time dimension
void wilson_slice_time_dir(Gauge_Conf const * const GC,
                           Geometry const * const geo,
                           int R,
                           int T,
                           int slice,
                           int dir,
                           double *re,
                           double *im)
   {
   long r, rsp;
   double rewilson, imwilson, rewil_memo, imwil_memo;

   #ifdef DEBUG
   // check if direction is not greater than STDIM or null
   if(dir>=STDIM)
      {
      fprintf(stderr, "Spatial direction chosen for the Wilson loop greater than STDIM! %d>=%d (%s, %d)\n", dir, STDIM, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   else if(dir==0)
      {
      fprintf(stderr, "Spatial direction chosen for the Wilson loop is actually time dimension! dir=%d (%s, %d)\n", dir, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   // check if slice to compute Wilson loop on
   // is inside the lattice
   if(slice>=(geo->d_size[0]))
      {
      fprintf(stderr, "Time slice to compute Wilson loop on greater than time size! %d>=%d (%s, %d)\n", slice, (geo->d_size[0]), __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   // check if sizes of Wilson loop are smaller
   // of the corresponding lattice sizes
   if(R>(geo->d_size[dir]))
      {
      fprintf(stderr, "Spatial size R of Wilson loop greater than corresponding spatial lattice size! %d>%d (%s, %d)\n", R, (geo->d_size[dir]), __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   if((slice+T)>(geo->d_size[0]))
      {
      fprintf(stderr, "Temporal size T of Wilson loop crossing the corresponding lattice border! %d+%d>%d (%s, %d)\n", T, slice, geo->d_size[0], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   #endif

   rewilson=0.;
   imwilson=0.;

   rewil_memo=0.;
   imwil_memo=0.;

   for(rsp=0; rsp<(geo->d_space_vol); rsp++)
      {
      r=sisp_and_t_to_si(geo, rsp, slice);
      
      wilson_fixed_site(GC, geo, r, R, T, dir, 0, &rewilson, &imwilson);

      rewil_memo+=rewilson;
      imwil_memo+=imwilson;
      }

   *re=rewil_memo*(geo->d_inv_space_vol);
   *im=imwil_memo*(geo->d_inv_space_vol);
   }

// compute the mean Wilson loop of size R x T
// on every spatial direction
// where R is a distance in space
// and T is a distance in time 
// on a whole slice in time dimension
void wilson_slice_time(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       int R,
                       int T,
                       int slice,
                       double *re,
                       double *im)
   {
   double rewilson, imwilson, rewil_memo, imwil_memo;

   #ifdef DEBUG
   // check if slice to compute Wilson loop on is inside the lattice
   if(slice>=(geo->d_size[0]))
      {
      fprintf(stderr, "Time slice to compute Wilson loop on greater than time size! %d>=%d (%s, %d)\n", slice, (geo->d_size[0]), __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   // check if sizes of Wilson loop are smaller
   // of the corresponding lattice sizes
   for(int dir=1; dir<STDIM; dir++)
      {
      if(R>(geo->d_size[dir]))
      {
         fprintf(stderr, "Spatial size R of Wilson loop greater than corresponding spatial lattice size! %d>%d (%s, %d)\n", R, (geo->d_size[dir]), __FILE__, __LINE__);
         exit(EXIT_FAILURE);
      }
      }
   if((slice+T)>(geo->d_size[0]))
      {
      fprintf(stderr, "Temporal size T of Wilson loop crossing the corresponding lattice border! %d+%d>%d (%s, %d)\n", T, slice, geo->d_size[0], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   #endif

   rewilson=0.;
   imwilson=0.;

   rewil_memo=0.;
   imwil_memo=0.;

   for(int i=1; i<STDIM; i++)
      {
      wilson_slice_time_dir(GC, geo, R, T, slice, i, &rewilson, &imwilson);

      rewil_memo+=rewilson;
      imwil_memo+=imwilson;
      }
   
   *re=(double) (rewil_memo/(STDIM-1));
   *im=(double) (imwil_memo/(STDIM-1));
   }

// compute the local topological charge at point r
// see readme for more details
double loc_topcharge(Gauge_Conf const * const GC,
                     Geometry const * const geo,
                     GParam const * const param,
                     long r)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     (void) GC;
     (void) geo;
     (void) param;
     (void) r;
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double ris=0.0; // initialized just to avoid compiler warnings

   #if (STDIM==4)
     GAUGE_GROUP aux1, aux2, aux3;
     double real1, real2, loc_charge;
     const double chnorm=1.0/(128.0*PI*PI);
     int i, dir[4][3], sign;

     dir[0][0] = 0;
     dir[0][1] = 0;
     dir[0][2] = 0;

     dir[1][0] = 1;
     dir[1][1] = 2;
     dir[1][2] = 3;

     dir[2][0] = 2;
     dir[2][1] = 1;
     dir[2][2] = 1;

     dir[3][0] = 3;
     dir[3][1] = 3;
     dir[3][2] = 2;

     sign=-1;
     loc_charge=0.0;

     for(i=0; i<3; i++)
        {
        clover(GC, geo, r, dir[0][i], dir[1][i], &aux1);
        clover(GC, geo, r, dir[2][i], dir[3][i], &aux2);

        times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
        real1=retr(&aux3)*NCOLOR;

        times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
        real2=retr(&aux3)*NCOLOR;

        loc_charge+=((double) sign*(real1-real2));
        sign=-sign;
        }
     ris=(loc_charge*chnorm);
   #endif

   return ris;
   }


// compute the topological charge
// see readme for more details
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double ris;
   long r;

   ris=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
   #endif
   for(r=0; r<(geo->d_volume); r++)
      {
      ris+=loc_topcharge(GC, geo, param, r);
      }

   return ris;
   }


// compute GParam::d_nummeas values of the topological charge after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topcharge_cooling(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
                       double *meanplaq)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   if(param->d_coolsteps>0)  // if using cooling
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, geo);
     // helperconf is a copy of the configuration
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        cooling(&helperconf, geo, param->d_coolsteps);

        ris=topcharge(&helperconf, geo, param);
        charge[iter]=ris;

        plaquette(&helperconf, geo, &plaqs, &plaqt);
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }

     free_gauge_conf(&helperconf, geo);
     }
   else   // no cooling
     {
     double ris, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, geo, param);
     plaquette(GC, geo, &plaqs, &plaqt);
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        charge[iter]=ris;
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }
     } 
   }


void perform_measures_localobs(Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep,
                               FILE *monofilep)
   {
   double plaqs, plaqt;
   //double plaqpar, plaqort;
   double polyre, polyim;
   double wilre, wilim;
   //double prod_polyre, prod_polyim;
   int dist_max=param->d_dist_poly;

   //plaquette(GC, geo, &plaqs, &plaqt);
   //polyakov(GC, geo, &polyre, &polyim);

   fprintf(datafilep, "%ld ", GC->update_index);
   //fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);

   // quality check
   /*for(int dist=1; dist<=dist_max; dist++)
      {
      polyakov_horizontal_corr_bulk_dir1(GC, geo, dist, &polyre, &polyim);
      fprintf(datafilep, "%d  %.12g %.12g ", dist,  polyre, polyim);
      }*/

   // to measure correlators between spins on bottom face
   for(int dist=1; dist<=dist_max; dist++)
      {
      wilre=0.;
      wilim=0.;
      
      wilson_slice_time(GC, geo, dist, 1, 0, &wilre, &wilim);
      fprintf(datafilep, "%d %.12g %.12g ", dist, wilre, wilim);
      }

   // to measure correlators between spins of PCM
   /*for(int dist=1; dist<=dist_max; dist++)
      {
      wilre=0.;
      wilim=0.;
      
      wilson_slice_time(GC, geo, dist, 1, (geo->d_size[0])-1, &wilre, &wilim);
      fprintf(datafilep, "%d %.12g %.12g ", dist, wilre, wilim);
      }*/
      
   // to measure correlators between Polyakov loops
   for(int dist=1; dist<=dist_max; dist++)
      {   
      for(int dir=1; dir<STDIM; dir++)
         {
         polyre=0.;
         polyim=0.;

         polyakov_correlator_dir(GC, geo, dist, dir, &polyre, &polyim);
         fprintf(datafilep, "%d %d %.12g %.12g ", dist, dir,  polyre, polyim);
         }
      }

   /*for(int dist=1; dist<=dist_max; dist++)
      {
      for(int slice=1; slice<geo->d_size[0]-1; slice++)
         {
         polyre=0.;
         polyim=0.;

         polyakov_horizontal_correlator_timeslice(GC, geo, dist, slice, &polyre, &polyim);
         fprintf(datafilep, "%d %d %.12g %.12g ", dist, slice,  polyre, polyim);
         }
      }*/

   for(int slice=0; slice<(geo->d_size[0]); slice++)
      {
      plaqs=0.;
      plaqt=0.;

      plaquette_slice_time(GC, geo, slice, &plaqs, &plaqt);
      fprintf(datafilep, "%d %.12g %.12g ", slice,  plaqs, plaqt);
      }
   
   // topological observables
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
     int i, err;
     double *charge, *meanplaq, charge_nocooling;

     charge_nocooling=topcharge(GC, geo, param);
     fprintf(datafilep, " %.12g ", charge_nocooling);

     err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     topcharge_cooling(GC, geo, param, charge, meanplaq);
     for(i=0; i<param->d_coolrepeat; i++)
        {
        fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
        }
     fprintf(datafilep, "\n");

     free(charge);
     free(meanplaq);
   #else
     fprintf(datafilep, "\n");
   #endif
   fflush(datafilep);

   // monopole observables
   if(param->d_mon_meas == 1)
     {
     #if STDIM==4
     int subg, subgnum;
     Gauge_Conf helperconf;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, geo);
     alloc_diag_proj_stuff(&helperconf, geo);

     // MAG gauge fixing
     max_abelian_gauge_fix(&helperconf, geo);
 
     //diagonal projection
     diag_projection(&helperconf, geo);
   
     //loop on all the U(1) subgroups
     if(NCOLOR>1)
       {
       subgnum=NCOLOR-1;
       }
     else
       {
       subgnum=1;
       }
     for(subg=0; subg<subgnum; subg++)
        {
        // extract the abelian component subg and save it to GC->u1_subg
        U1_extract(&helperconf, geo, subg);

        // compute monopole observables
        monopoles_obs(&helperconf, geo, param, subg, monofilep);
        }

     free_diag_proj_stuff(&helperconf, geo);
     free_gauge_conf(&helperconf, geo);

     fflush(monofilep);
     #else
     (void) monofilep;
     #endif
     }
   }


// perform local observables in the case of trace deformation, it computes all the order parameters
void perform_measures_localobs_with_tracedef(Gauge_Conf const * const GC,
                                             Geometry const * const geo,
                                             GParam const * const param,
                                             FILE *datafilep,
                                             FILE *monofilep)
   {
   int i;
   double plaqs, plaqt, polyre[NCOLOR/2+1], polyim[NCOLOR/2+1]; // +1 just to avoid warning if NCOLOR=1

   plaquette(GC, geo, &plaqs, &plaqt);
   polyakov_for_tracedef(GC, geo, polyre, polyim);

   fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);

   for(i=0; i<(int)floor(NCOLOR/2); i++)
      {
      fprintf(datafilep, "%.12g %.12g ", polyre[i], polyim[i]);
      }

   // topological observables
   #if(STDIM==4)
     int err;
     double *charge, *meanplaq, charge_nocooling;

     charge_nocooling=topcharge(GC, geo, param);

     fprintf(datafilep, "%.12g ", charge_nocooling);

     err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     topcharge_cooling(GC, geo, param, charge, meanplaq);
     for(i=0; i<param->d_coolrepeat; i++)
        {
        fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
        }
     fprintf(datafilep, "\n");

     free(charge);
     free(meanplaq);
   #else
     fprintf(datafilep, "\n");
   #endif

   fflush(datafilep);

   // monopole observables
   if(param->d_mon_meas == 1)
     {
     #if(STDIM==4)
     Gauge_Conf helperconf;
     int subg, subgnum;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, geo);
     alloc_diag_proj_stuff(&helperconf, geo);

     // MAG gauge fixing
     max_abelian_gauge_fix(&helperconf, geo);

     //diagonal projection
     diag_projection(&helperconf, geo);

     //loop on all the U(1) subgroups
     if(NCOLOR>1)
       {
       subgnum=NCOLOR-1;
       }
     else
       {
       subgnum=1;
       }
     for(subg=0; subg<subgnum; subg++)
        {
        // extract the abelian component subg and save it to GC->u1_subg
        U1_extract(&helperconf, geo, subg);

        // compute monopole observables
        monopoles_obs(&helperconf, geo, param, subg, monofilep);
        }

     free_diag_proj_stuff(&helperconf, geo);
     free_gauge_conf(&helperconf, geo);

     fflush(monofilep);
     #else
     (void) monofilep;
     #endif
     }
   }


// compute the average value of \sum_{flavours} Re(H_x U_{x,mu} H_{x+mu})
void higgs_interaction(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       double *he)
  {
  long r;
  double ris=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
  #endif
  for(r=0; r<(geo->d_volume); r++)
     {
     int i;
     double aux=0.0;
     GAUGE_VECS v1;
     GAUGE_GROUP matrix;

     for(i=0; i<STDIM; i++)
        {
        equal(&matrix, &(GC->lattice[r][i]));

        matrix_times_vector_all_vecs(&v1, &matrix, &(GC->higgs[nnp(geo, r, i)]));
        aux+=re_scal_prod_vecs(&(GC->higgs[r]), &v1);
        }

     ris+=aux;
     }

  ris/=(double) STDIM;
  ris*=geo->d_inv_vol;

  *he=ris;
  }


// compute flavour related observables
//
// flavour matrices Qh and Dh HAVE TO BE INITIALIZED before calling this function
//
// tildeG0=ReTr[(\sum_x Q_x)(\sum_y Q_y)]/volume/NHIGGS
// tildeGminp=ReTr[(\sum_x Q_xe^{ipx})(\sum_y Q_ye^{-ipy)]/volume/NHIGGS
//
// tildeG0 is susceptibility/NHIGGS, tildeGminp is used to compute the 2nd momentum correlation function
//
// tildeD0=conj(\sum_x D_x) (\sum_y D_y) / volume
// tildeDminp=(\sum_x D_x e^{ipx}) conj(\sum_y D_y e^{ipy}) /volume
//
// tildeD0 is a U1 susceptibility, tildeDminp is used to compute the 2nd momentum correlation function
void compute_flavour_observables(Gauge_Conf const * const GC,
                                 Geometry const * const geo,
                                 double *tildeG0,
                                 double *tildeGminp,
                                 double *tildeD0,
                                 double *tildeDminp)
  {
  int coord[STDIM];
  long r;
  const double p = 2.0*PI/(double)geo->d_size[1];
  double complex D, Dp;
  FMatrix Q, Qp, Qmp, tmp1, tmp2;

  // Q =sum_x Q_x
  // Qp=sum_x e^{ipx}Q_x
  // Qmp=sum_x e^{-ipx}Q_x
  //
  // D, Dp and are the analogous of Q and Qp for D

  D=0.0+0.0*I;
  Dp=0.0+0.0*I;

  zero_FMatrix(&Q);
  zero_FMatrix(&Qp);
  zero_FMatrix(&Qmp);
  for(r=0; r<(geo->d_volume); r++)
     {
     equal_FMatrix(&tmp1, &(GC->Qh[r]));
     equal_FMatrix(&tmp2, &tmp1);

     plus_equal_FMatrix(&Q, &tmp1);
     D+=(GC->Dh[r]);

     si_to_cart(coord, r, geo);

     times_equal_complex_FMatrix(&tmp1, cexp(I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qp, &tmp1);
     Dp+=((GC->Dh[r]) * cexp(I*((double)coord[1])*p) );

     times_equal_complex_FMatrix(&tmp2, cexp(-I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qmp, &tmp2);
     }

  equal_FMatrix(&tmp1, &Q);
  times_equal_FMatrix(&tmp1, &Q);

  *tildeG0=retr_FMatrix(&tmp1)*geo->d_inv_vol;
  *tildeD0=creal(conj(D)*D)*geo->d_inv_vol;

  equal_FMatrix(&tmp1, &Qp);
  times_equal_FMatrix(&tmp1, &Qmp);
  *tildeGminp=retr_FMatrix(&tmp1)*geo->d_inv_vol;
  *tildeDminp=creal(Dp*conj(Dp))*geo->d_inv_vol;
  }


// compute correlators of flavour observables
//
// flavour matrices Qh and Dh HAVE TO BE INITIALIZED before calling this function
//
// corrQQ is the correlato ReTr[Q_x Q_{x+d}]/N_higgs
// corr0string0 is the correlator \sum_f Re[hf^{dag} U_{x,1}U_{x+1,1}....Q_{x+d-1,1} hf], where hf is the f-th flavour
// corr0string1 is the correlator Re[h0^{dag} U_{x,1}U_{x+1,1}....U_{x+d-1,1} h1], where h1 is the second flavour
void compute_flavour_observables_corr(Gauge_Conf const * const GC,
                                      Geometry const * const geo,
                                      double *corrQQ,
                                      double *corr0string0,
                                      double *corr0string1)
  {
  int dist;
  long r;
  double accumulator1, accumulator2;

  for(dist=0; dist<geo->d_size[1]; dist++)
     {
     accumulator1=0.0;

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : accumulator1)
     #endif
     for(r=0; r<(geo->d_volume); r++)
        {
        int i;
        long r1;
        FMatrix tmp1;

        equal_FMatrix(&tmp1, &(GC->Qh[r]));
        r1=r;
        for(i=0; i<dist; i++)
           {
           r1=nnp(geo, r1, 1);
           }
        times_equal_FMatrix(&tmp1, &(GC->Qh[r1]));
        accumulator1+=retr_FMatrix(&tmp1);
        }
     accumulator1*=geo->d_inv_vol;
     corrQQ[dist]=accumulator1;

     accumulator1=0.0;
     accumulator2=0.0;

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : accumulator1) reduction(+ : accumulator2)
     #endif
     for(r=0; r<(geo->d_volume); r++)
        {
        int i;
        long r1;
        GAUGE_VECS phi1, phi2;
        GAUGE_GROUP U;

        equal_vecs(&phi1, &(GC->higgs[r]));
        r1=r;
        one(&U);
        for(i=0; i<dist; i++)
           {
           times_equal(&U, &(GC->lattice[r1][1]));
           r1=nnp(geo, r1, 1);
           }
        matrix_times_vector_all_vecs(&phi2, &U, &(GC->higgs[r1]));
        accumulator1+=re_scal_prod_vecs(&phi1, &phi2);
        #if NHIGGS >1
         accumulator2+=re_scal_prod_single_vecs(&phi1, &phi2, 0, 1);
        #else
         accumulator2+=0.0;
        #endif
        }
     accumulator1*=geo->d_inv_vol;
     accumulator2*=geo->d_inv_vol;

     corr0string0[dist]=accumulator1;
     corr0string1[dist]=accumulator2;
     }
  }


void perform_measures_higgs(Gauge_Conf *GC,
                            Geometry const * const geo,
                            FILE *datafilep)
   {
   double plaqs, plaqt, polyre, polyim, he, tildeG0, tildeGminp, tildeD0, tildeDminp;
   long r;

   plaquette(GC, geo, &plaqs, &plaqt);
   polyakov(GC, geo,  &polyre, &polyim);
   higgs_interaction(GC, geo, &he);

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(geo->d_volume); r++)
      {
      init_FMatrix_vecs(&(GC->Qh[r]), &(GC->higgs[r]));
      GC->Dh[r] = HiggsU1Obs_vecs(&(GC->higgs[r]));
      }

   compute_flavour_observables(GC,
                               geo,
                               &tildeG0,
                               &tildeGminp,
                               &tildeD0,
                               &tildeDminp);

   fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
   fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
   fprintf(datafilep, "%.12g ", he);
   fprintf(datafilep, "%.12g %.12g ", tildeG0, tildeGminp);
   fprintf(datafilep, "%.12g %.12g ", tildeD0, tildeDminp);

   /*
   // for correlators

   int err, i;
   double *corrQQ, *corr0string0, *corr0string1;
   err=posix_memalign((void**) &(corrQQ), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[1] * sizeof(double));
   err+=posix_memalign((void**) &(corr0string0), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[1] * sizeof(double));
   err+=posix_memalign((void**) &(corr0string1), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[1] * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating the correlators! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   compute_flavour_observables_corr(GC,
                                    geo,
                                    param,
                                    corrQQ,
                                    corr0string0,
                                    corr0string1);
   for(i=0; i<param->d_size[1]; i++)
      {
      fprintf(datafilep, "%.12g ", corrQQ[i]);
      }
   for(i=0; i<param->d_size[1]; i++)
      {
      fprintf(datafilep, "%.12g ", corr0string0[i]);
      }
    for(i=0; i<param->d_size[1]; i++)
      {
      fprintf(datafilep, "%.12g ", corr0string1[i]);
      }

   free(corrQQ);
   free(corr0string0);
   free(corr0string1);
   */

   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


// this is a function to be used just to test some fine points
// most notably TrP^2=TrQ^2+1/NHIGGS
void perform_measures_higgs_for_testing(Gauge_Conf *GC,
                                        Geometry const * const geo,
                                        FILE *datafilep)
   {
   double plaqs, plaqt, polyre, polyim, he, p2, tildeG0, tildeGminp, tildeD0, tildeDminp;
   long r;

   plaquette(GC, geo, &plaqs, &plaqt);
   polyakov(GC, geo, &polyre, &polyim);
   higgs_interaction(GC, geo, &he);

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(geo->d_volume); r++)
      {
      init_FMatrix_vecs(&(GC->Qh[r]), &(GC->higgs[r]));
      GC->Dh[r] = HiggsU1Obs_vecs(&(GC->higgs[r]));
      }

   fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
   fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
   fprintf(datafilep, "%.12g ", he);

   compute_flavour_observables(GC,
                               geo,
                               &tildeG0,
                               &tildeGminp,
                               &tildeD0,
                               &tildeDminp);

   fprintf(datafilep, "%.12g %.12g ", tildeG0, tildeGminp);
   fprintf(datafilep, "%.12g %.12g ", tildeD0, tildeDminp);

   p2=0.0;
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+: p2)
   #endif
   for(r=0; r<(geo->d_volume); r++)
      {
      FMatrix tmp1, tmp2;
      equal_FMatrix(&tmp1, &(GC->Qh[r]));
      equal_FMatrix(&tmp2, &tmp1);
      times_equal_FMatrix(&tmp1, &tmp2);
      p2+=retr_FMatrix(&tmp1)*NHIGGS+1./NHIGGS;
      }
   p2*=geo->d_inv_vol;

   fprintf(datafilep, "%.12g ", p2);

   fprintf(datafilep, "\n");

   fflush(datafilep);
   }



// fix maximal abelian gauge
// following the procedure described in
// C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ]
void max_abelian_gauge_fix(Gauge_Conf *GC,
                           Geometry const * const geo)
   {
   int i, dir;
   long r;
   double lambda[NCOLOR];
   const double overrelaxparam=1.85; // 1.0 means no overrelaxation
   const double target=1.0e-8;
   double nondiag, nondiagaux;

   // inizialize the matrix lambda = diag((N-1)/2, (N-1)/2-1, ..., -(N-1)/2)
   for(i=0; i<NCOLOR; i++)
      {
      lambda[i] = ( (double) NCOLOR -1.)/2. - (double) i;
      }

   nondiag=1;
   while(nondiag > target)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
        #endif
        for(r=0; r<geo->d_volume/2; r++)
           {
           GAUGE_GROUP G_mag, help, X_links[2*STDIM];   // X_links contains the 2*STDIM links used in the computation of X(n)

           // initialize X_links[2*STDIM] with the 2*STDIM links surrounding the point r
           // links 0 to (STDIM-1) are forward, while links STDIM to (2*STDIM-1) are backwards.
           for(dir=0; dir<STDIM; dir++)
              {
              equal(&(X_links[dir]), &(GC->lattice[r][dir]));
              equal(&(X_links[dir+STDIM]), &(GC->lattice[nnm(geo, r, dir)][dir]));
              }

           comp_MAG_gauge_transformation(X_links, lambda, overrelaxparam, &G_mag);
 
           // apply the gauge transformation
           for(dir=0; dir<STDIM; dir++)
              {
              times(&help, &G_mag, &(GC->lattice[r][dir]));
              equal(&(GC->lattice[r][dir]), &help);

              times_equal_dag(&(GC->lattice[nnm(geo, r, dir)][dir]), &G_mag);
              }
           }

        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
        #endif
        for(r=geo->d_volume/2; r<geo->d_volume; r++)
           {
           GAUGE_GROUP G_mag, help, X_links[2*STDIM];   // X_links contains the 2*STDIM links used in the computation of X(n)

           // initialize X_links[2*STDIM] with the 2*STDIM links surrounding the point r
           // links 0 to (STDIM-1) are forward, while links STDIM to (2*STDIM-1) are backwards.
           for(dir=0; dir<STDIM; dir++)
              {
              equal(&(X_links[dir]), &(GC->lattice[r][dir]));
              equal(&(X_links[dir+STDIM]), &(GC->lattice[nnm(geo, r, dir)][dir]));
              }

           comp_MAG_gauge_transformation(X_links, lambda, overrelaxparam, &G_mag);

           // apply the gauge transformation
           for(dir=0; dir<STDIM; dir++)
              {
              times(&help, &G_mag, &(GC->lattice[r][dir]));
              equal(&(GC->lattice[r][dir]), &help);

              times_equal_dag(&(GC->lattice[nnm(geo, r, dir)][dir]), &G_mag);
              }
           }

        // nondiagaux is the sum of the squares of the out-diagonal terms
        nondiagaux=0;

        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r, dir)  reduction(+ : nondiagaux)
        #endif
        for(r=0; r<geo->d_volume; r++)
           {
           GAUGE_GROUP X_links[2*STDIM];   // X_links contains the 2*STDIM links used in the computation of X(n)
           double counter;

           for(dir=0; dir<STDIM; dir++)
              {
              equal(&(X_links[dir]), &(GC->lattice[r][dir]));
              equal(&(X_links[dir+STDIM]), &(GC->lattice[nnm(geo, r, dir)][dir]));
              }
           comp_outdiagnorm_of_X(X_links, lambda, &counter);
           nondiagaux += counter;
           }
     
        nondiag = nondiagaux * geo->d_inv_vol / (double)NCOLOR / (double) NCOLOR;

        // printf("%g  %g\n", nondiag, nondiag/target);
        // fflush(stdout);
        }

   // unitarize all the links
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif
   for(r=0; r<(geo->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         }
      }
   } 


// extract the diagonal part of the links after gauge fixing.
// the phases are saved in GC->diag_proj but these are NOT the monopole phases (see U1_extract)
void diag_projection(Gauge_Conf *GC,
                     Geometry const * const geo)
   {
   int dir; 
   long r;

   for(r=0;r<geo->d_volume;r++)
      {
      for(dir=0;dir<STDIM;dir++)
         {
         diag_projection_single_site(GC, &(GC->lattice[r][dir]), r, dir);
         }
      }
   }


// extract the abelian components of the link
// following the procedure described in
// Bonati, D'Elia https://arxiv.org/abs/1308.0302
// and save them in GC->u1_subg
//
// also intialize GC->uflag to zero
void U1_extract(Gauge_Conf *GC, 
                Geometry const * const geo,
                int subg)
   { 
   int dir, i;
   long r;    

   for(r=0;r<geo->d_volume;r++)
      {
      for(dir=0;dir<STDIM;dir++)
         {
         GC->u1_subg[r][dir] = 0.0;
         for(i=0;i<=subg;i++)
            {
            GC->u1_subg[r][dir] += GC->diag_proj[r][dir][i];
            }

         GC->uflag[r][dir] = 0;
         }
      }
   }


// Compute the forward derivative of the abelian part of the plaquette Fjk in direction i.
// the angle is chosen in between -pi and pi.
void Di_Fjk(Gauge_Conf *GC,
            Geometry const * const geo,
            long r,
            int idir,
            int jdir,
            int kdir,
            double *DiFjk)

   {
   double prpi, pr; // pr -> plaquette at site r, prpi plaquette at site r+idir

//
//       ^ k
//       |
//       +---<---+
//       |       |
//       V       ^         pr
//       |       |
//       +--->---+---> j
//       r
//

   pr  = GC->u1_subg[r][jdir] - GC->u1_subg[r][kdir];
   pr += GC->u1_subg[nnp(geo, r, jdir)][kdir] - GC->u1_subg[nnp(geo, r, kdir)][jdir];


//
//       ^ k
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)        prpi
//       |       |
//       +--->---+---> j
//      r+i    (4)
//

   r=nnp(geo, r, idir);

   prpi  = GC->u1_subg[r][jdir] - GC->u1_subg[r][kdir];
   prpi += GC->u1_subg[nnp(geo, r, jdir)][kdir] - GC->u1_subg[nnp(geo, r, kdir)][jdir];
 
   *DiFjk = 2.0*(atan(tan(prpi/2.0)) - atan(tan(pr/2.0)));
   }


// compute the DeGrand-DeTar currents
int DeGrand_current(Gauge_Conf *GC,
                    Geometry const * const geo,
                    long r,
                    int dir)
   {
   if(STDIM!=4)
     {
     fprintf(stderr, "Wrong number of dimensions! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double der1, der2, der3;
   int ris;

   if(dir == 0)
     {
     Di_Fjk(GC, geo, r, 1, 2, 3, &der1);
     Di_Fjk(GC, geo, r, 3, 1, 2, &der2);
     Di_Fjk(GC, geo, r, 2, 1, 3, &der3);
   
     ris = (int) round( ((der1 + der2 - der3)/PI2) );
     }
   else if(dir ==1)
          {
          Di_Fjk(GC, geo, r, 3, 2, 0, &der1);
          Di_Fjk(GC, geo, r, 0, 3, 2, &der2);
          Di_Fjk(GC, geo, r, 2, 3, 0, &der3);

          ris = (int) round( ((der1 + der2 - der3)/PI2) );
          }
   else if(dir == 2)
          {
          Di_Fjk(GC, geo, r, 3, 0, 1, &der1);
          Di_Fjk(GC, geo, r, 0, 1, 3, &der2);
          Di_Fjk(GC, geo, r, 1, 0, 3, &der3);

          ris = (int) round( ((der1 + der2 - der3)/PI2) );
          }
   else
     {
     Di_Fjk(GC, geo, r, 0,2,1, &der1);
     Di_Fjk(GC, geo, r, 2,1,0, &der2);
     Di_Fjk(GC, geo, r, 1,2,0, &der3);

     ris = (int) round( ((der1 + der2 - der3)/PI2) );
     }

   return ris;
   } 


// search for monopole wrappings passing from r_tback
// this function can be invoked in two different ways
//
// or r=nnp(geo, r_tback, 0) and DeGrand_current(GC, geo, r_tback, 0)!=0 (forward case)
// or r=nnm(geo, r_tback, 0) and DeGrand_current(GC, geo, r, 0)!=0  (backward case)
//
// GC->uflag[][] is initialized in monopole_obs
//
// nonzero DeGrand_current(GC, geo, nnp(geo, r, dir), dir ) are associated to uflag[r][dir]
//
// num_wrap = number of wrappings
void wrap_search(Gauge_Conf *GC,
                 Geometry const * const geo,
                 GParam const * const param,
                 long r,
                 long r_tback,
                 int *num_wrap)
   {
   #if STDIM!=4
     fprintf(stderr, "Wrong number of dimensions! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
   #endif

   int dir, n_mu;

   if(r == r_tback)
     {
     return;
     }
   else
     {
     // forward case
     for(dir=0; dir<STDIM; dir++)
        {
        n_mu=DeGrand_current(GC, geo, nnp(geo, r, dir), dir);

        // if not all the monopole currents have been followed
        if(n_mu > GC->uflag[r][dir])
          {
          GC->uflag[r][dir] += 1;

          if( (geo->d_timeslice[r] == geo->d_size[0]-1) && (dir == 0) )
            {
            *num_wrap += 1;
            }

          wrap_search(GC, geo, param, nnp(geo, r, dir), r_tback, num_wrap);

          return;
          }
        }

     //backward case
     for(dir=0;dir<STDIM;dir++)
        {
        n_mu=DeGrand_current(GC, geo, r, dir);

        if(n_mu < GC->uflag[nnm(geo, r, dir)][dir])
          {
          GC->uflag[nnm(geo, r, dir)][dir] -= 1;

          if( (geo->d_timeslice[r] == 0) && (dir == 0) )
            {
            *num_wrap -= 1;
            }

          wrap_search(GC, geo, param, nnm(geo, r, dir), r_tback, num_wrap);

          return;
          }
        }
     }
   }


// GC->uflag[][] has to be initialized to zero before calling this function
// (when GC->uflag is allocated it is also initialized to zero)
void monopoles_obs(Gauge_Conf *GC, 
                   Geometry const * const geo,
                   GParam const * const param, 
                   int subg, 
                   FILE* monofilep)
   {
   double mean_wrap;
   long r, rsp, r_tback, r_tbackback;
   int n_mu, num_wrap, mono_charge;
   int cartcoord[4];

   mean_wrap = 0.0;     // mean value of monopole wraps for unit volume

   for(rsp=0; rsp<geo->d_space_vol; rsp++)
      {
      r = sisp_and_t_to_si(geo, rsp, 1);                             // t=1 slice
      r_tback = sisp_and_t_to_si(geo, rsp, 0);                       // t=0 slice
      r_tbackback = sisp_and_t_to_si(geo, rsp, geo->d_size[0]-1);  // t=T-1 slice

      // check the t=1 temporal slice to find monopoles currents
      n_mu=DeGrand_current(GC, geo, r, 0);

      // start following monopole charge in forward direction. Maximum lattice charge is +2 so we try twice
      for(mono_charge = 0; mono_charge<2; mono_charge++)
         {
         // nonzero DeGrand_current(GC, geo, nnp(geo, r, dir), dir ) are associated to uflag[r][dir]
         if(n_mu > GC->uflag[r_tback][0])
           {
           GC->uflag[r_tback][0] += 1;

           num_wrap = 0;
           wrap_search(GC, geo, param, r, r_tback, &num_wrap);

           mean_wrap += abs(num_wrap);

           lexeo_to_cart(cartcoord, r_tback, geo);
           if(n_mu == 1)
             {
             fprintf(monofilep, "%ld ", GC->update_index);

             for(int k = 0; k< 4; k++)
                {
                fprintf(monofilep, "%d ", cartcoord[k]);
                }
             fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
             }
           else if(GC->uflag[r_tback][0] == 1) // this is to print only once monopole of charge +2
                  {
                  fprintf(monofilep, "%ld ", GC->update_index);

                  for(int k = 0; k<4; k++)
                     {
                     fprintf(monofilep, "%d ", cartcoord[k]);
                     }
                  fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
                  }
           }
         }

      n_mu=DeGrand_current(GC, geo, r_tback, 0);

      // start following monopole charge in backward direction. Maximum lattice charge is +2 so we try twice
      for(mono_charge = 0; mono_charge<2; mono_charge++)
         {
         // nonzero DeGrand_current(GC, geo, nnp(geo, r, dir), dir ) are associated to uflag[r][dir]
         if(n_mu < GC->uflag[r_tbackback][0])
           {
           GC->uflag[r_tbackback][0] -= 1;

           num_wrap = -1;
           wrap_search(GC, geo, param, r_tbackback, r_tback, &num_wrap);

           lexeo_to_cart(cartcoord, r_tback, geo);
           if(n_mu == -1)
             {
             fprintf(monofilep, "%ld ", GC->update_index);

             for(int k = 0; k<4; k++)
                {
                fprintf(monofilep, "%d ", cartcoord[k]);
                }
             fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
             }
           else if(GC->uflag[r][0] == -1)  // this is to print only once monopole of charge +2
                  {
                  fprintf(monofilep, "%ld ", GC->update_index);

                  for(int k=0; k<4; k++)
                     {
                     fprintf(monofilep, "%d ", cartcoord[k]);
                     }
                  fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
                  }
           }
         }
      }
   }


#endif












