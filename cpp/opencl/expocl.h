#ifndef __EXPOCL_H__
#define __EXPOCL_H__

#include <stdio.h>

#include "../shared/expo.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include "util.h"

#ifdef DOUBLE
typedef double clReal;
#else
typedef float clReal;
#endif

typedef struct _expo_cl {
  cl_int           error;          // Used to handle error codes
  cl_platform_id   platform;
  cl_context       context;
  cl_command_queue queue;
  cl_device_id     device;
  cl_program       program;
  cl_event         last_event;

  cl_kernel        inf_norm_kernel;
  cl_kernel        two_norm_kernel;
  cl_kernel        tridiag_kernel;
  cl_kernel        vector_add_kernel;
  cl_kernel        vector_copy_kernel;
  cl_kernel        vector_smult_kernel;
  cl_kernel        vector_one_norm_tridiag_mat_kernel;
  cl_kernel        expo_trans_kernel;
  cl_kernel        expo_samp_kernel;
  cl_kernel        neg_vals_kernel;

  cl_kernel        tridiag_sir_kernel;
  cl_kernel        expo_sir_trans_kernel;
  cl_kernel        expo_sir_samp_kernel;

  cl_event (*calc)(struct _expo_cl* ecl, expo_type* expo, char trans, clReal alpha);
  void (*trans)(struct _expo_cl* ecl, const expo_type* expo);
  void (*samp)(struct _expo_cl* ecl, const expo_type* expo);

  size_t global_ws;
  size_t local_ws;

  int size;
  int mem_size;

  /* Memory items */

  cl_mem ld;      /* lambda vector */
  cl_mem md;      /* mu vector */
  cl_mem pd;      /* psi vector */

  cl_mem Ad;      /* matrix */
  cl_mem xd;      /* work vector 1 */
  cl_mem yd;      /* work vector 2 */
  cl_mem bd;      /* result vector */

  cl_mem id;      /* index vector */
} expo_cl;

void clreal_to_double(int n, clReal* b, double* x);
void double_to_clreal(int n, double* x, clReal* b);

void cl_afun_power(char trans, double alpha, int m, clReal* x, 
                   expo_type* expo, expo_cl* ecl);
int expocl_normam(double alpha, int m, double* c, int* mv,
                  double* wrk, int* iwrk, expo_type* expo, expo_cl* ecl);

/* OpenCL functions *********************************************************/

cl_int expocl_finish(expo_cl* ecl);
double expocl_profile_last_event(expo_cl* ecl);

/****************************************************************************/

expo_cl* expocl_init(const expo_type* expo, cl_device_type dtype, int local_ws);
void expocl_free(expo_cl* ecl);

void expocl_init_kernel(expo_cl* ecl, const expo_type* expo);
void expocl_free_kernel(expo_cl* ecl);

void expocl_copy_mat(expo_cl* ecl, const expo_type* expo);
void expocl_copy_index(expo_cl* ecl, const expo_type* expo);
void expocl_copy_vectors(expo_cl* ecl, const expo_type* expo);

void expocl_create_buffers(expo_cl* ecl);
void expocl_release_buffers(expo_cl* ecl);

void expocl_create_vector_buffers(expo_cl* ecl);
void expocl_release_vector_buffers(expo_cl* ecl);

void expocl_create_index_buffer(expo_cl* ecl, const expo_type* expo);
void expocl_release_index_buffer(expo_cl* ecl);

void expocl_copy_b(expo_cl* ecl, clReal* x);
void expocl_copy_x(expo_cl* ecl, clReal* x);

void expocl_read_b(expo_cl* ecl, clReal* b);
void expocl_read_x(expo_cl* ecl, clReal* x);
void expocl_read_y(expo_cl* ecl, clReal* y);

void expocl_swap_xy(expo_cl* ecl);

cl_event expocl_sis_calc(expo_cl* ecl, expo_type* expo, char trans, clReal alpha);

clReal expocl_inf_norm_b(expo_cl* ecl);
clReal expocl_inf_norm_x(expo_cl* ecl);

clReal expocl_two_norm_b(expo_cl* ecl);
clReal expocl_two_norm_x(expo_cl* ecl);

clReal expocl_one_norm_mat(expo_cl* ecl);

void expocl_add_bx(expo_cl* ecl);
void expocl_add_by(expo_cl* ecl);

void expocl_copy_bx(expo_cl* ecl);
void expocl_copy_xb(expo_cl* ecl);
void expocl_copy_xy(expo_cl* ecl);

void expocl_smult_b(expo_cl* ecl, clReal alpha);

void expocl_trans_b(expo_cl* ecl, const expo_type* expo);
void expocl_samp_b(expo_cl* ecl, const expo_type* expo);

int expocl_neg_vals(expo_cl* ecl);

/* SIR SPECIFIC FUNCTIONS ***************************************************/

cl_event expocl_sir_calc(expo_cl* ecl, expo_type* expo, char trans, clReal alpha);
void expocl_sir_samp_b(expo_cl* ecl, const expo_type* expo);
void expocl_sir_trans_b(expo_cl* ecl, const expo_type* expo);

#endif
