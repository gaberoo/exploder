#include "expocl.h"

void clreal_to_double(int n, clReal* b, double* x) {
  int i;
  for (i = 0; i < n; ++i) x[i] = b[i];
}

void double_to_clreal(int n, double* x, clReal* b) {
  int i;
  for (i = 0; i < n; ++i) b[i] = x[i];
}

cl_int expocl_finish(expo_cl* ecl) { return clFinish(ecl->queue); }

/****************************************************************************/

double expocl_profile_last_event(expo_cl* ecl) {
  cl_ulong time_start;
  cl_ulong time_end;
  clGetEventProfilingInfo(ecl->last_event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
  clGetEventProfilingInfo(ecl->last_event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
  return time_end - time_start;
}

/****************************************************************************/

expo_cl* expocl_init(const expo_type* expo, cl_device_type dtype, int local_ws) 
{
  expo_cl* ecl = (expo_cl*) malloc(sizeof(expo_cl));
  ecl->calc = &expocl_sis_calc;
  ecl->samp = &expocl_samp_b;
  ecl->trans = &expocl_trans_b;

  // Platform
  ecl->error = oclGetPlatformID(&ecl->platform);
  checkErr(ecl->error,"Platform ID");

  // Device
  ecl->error = clGetDeviceIDs(ecl->platform, dtype, 1, &ecl->device, NULL);
  checkErr(ecl->error,"Device ID");

  // Context
  ecl->context = clCreateContext(0, 1, &ecl->device, NULL, NULL, &ecl->error);
  checkErr(ecl->error,"Create context");

  // Command-queue
  // ecl->queue = clCreateCommandQueue(ecl->context, ecl->device, 0, &ecl->error);
  ecl->queue = clCreateCommandQueue(ecl->context, ecl->device, CL_QUEUE_PROFILING_ENABLE, &ecl->error);
  checkErr(ecl->error,"Create command queue");

  size_t max_gsize = 0;
  ecl->error = clGetDeviceInfo(ecl->device,CL_DEVICE_MAX_WORK_GROUP_SIZE,
                               sizeof(size_t), &max_gsize, NULL);

  ecl->local_ws = (local_ws < max_gsize) ? local_ws : max_gsize;
  ecl->global_ws = shrRoundUp(ecl->local_ws, expo->dim);

  ecl->size = expo->dim;
  ecl->mem_size = ecl->size*sizeof(clReal);

  return ecl;
}

/****************************************************************************/

void expocl_free(expo_cl* ecl) {
  clReleaseCommandQueue(ecl->queue);
  clReleaseContext(ecl->context);
}

/****************************************************************************/

void expocl_init_kernel(expo_cl* ecl, const expo_type* expo) {
  char* kernel_str = NULL;

  char type = 's';
#ifdef DOUBLE
  type = 'd';
#endif

  // Currently loading kernels from file. Should be hard-coded in the future.
  load_kernel_source("expmv_kernel.cl",&kernel_str,type);

  const char* kstr = kernel_str;
  // const char* kstr = program_str;
  // const char** kstr = (const char**) kernel_str;

  ecl->program = clCreateProgramWithSource(ecl->context, 1, 
                                           &kstr, NULL, &ecl->error);
  checkErr(ecl->error,"Create program");

  // Builds the program
  ecl->error = clBuildProgram(ecl->program, 1, &ecl->device, NULL, NULL, NULL);

  if (ecl->error != CL_SUCCESS) {
    printf("%s\n",kernel_str);

    cl_build_status status;
    clGetProgramBuildInfo(ecl->program, ecl->device, 
                          CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), 
                          &status, NULL);

    if (status == CL_BUILD_NONE) {
      printf("Not built.\n");
    } else if (status == CL_BUILD_ERROR) {
      printf("Build error.\n");
    }

    size_t log_size;
    clGetProgramBuildInfo(ecl->program, ecl->device, CL_PROGRAM_BUILD_LOG, 
                          0, NULL, &log_size);
    printf("Build log length = %lu\n",log_size);

    // Second call to get the log
    char* build_log = (char*) malloc((log_size+1)*sizeof(char));
    clGetProgramBuildInfo(ecl->program, ecl->device, CL_PROGRAM_BUILD_LOG, 
                          log_size, build_log, NULL);
    build_log[log_size] = '\0';

    printf("%s\n",build_log);

    print_info();

    free(build_log);
    checkErr(ecl->error,"Build program");
  }

  free(kernel_str);

  /* extracting the kernels */

  printf("SIS kernels "); fflush(stdout);
  ecl->tridiag_kernel = clCreateKernel(ecl->program, "tridiag_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'tridiag'");

  ecl->inf_norm_kernel = clCreateKernel(ecl->program, "inf_norm_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'inf norm'");

  ecl->two_norm_kernel = clCreateKernel(ecl->program, "two_norm_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'two norm'");

  ecl->vector_add_kernel = clCreateKernel(ecl->program, "vector_add_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 3");

  ecl->vector_copy_kernel = clCreateKernel(ecl->program, "vector_copy_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 4");

  ecl->vector_smult_kernel = clCreateKernel(ecl->program, "vector_smult_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 5");

  ecl->vector_one_norm_tridiag_mat_kernel = clCreateKernel(ecl->program, "one_norm_tridiag_mat_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'one_norm_mat'");

  ecl->expo_trans_kernel = clCreateKernel(ecl->program, "expo_trans_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'expo_trans_kernel'");

  ecl->expo_samp_kernel = clCreateKernel(ecl->program, "expo_samp_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'expo_samp_kernel'");

  ecl->neg_vals_kernel = clCreateKernel(ecl->program, "neg_vals_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'neg_vals_kernel'");

  printf("SIR kernels "); fflush(stdout);
  ecl->tridiag_sir_kernel = clCreateKernel(ecl->program, "tridiag_sir_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'tridiag index'");

  ecl->expo_sir_trans_kernel = clCreateKernel(ecl->program, "expo_sir_trans_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'expo_sir_trans_kernel'");

  ecl->expo_sir_samp_kernel = clCreateKernel(ecl->program, "expo_sir_samp_kernel", &ecl->error);
  checkErr(ecl->error,"Create kernel 'expo_sir_samp_kernel'");
}

/****************************************************************************/

void expocl_free_kernel(expo_cl* ecl) {
  clReleaseKernel(ecl->vector_smult_kernel);
  clReleaseKernel(ecl->vector_add_kernel);
  clReleaseKernel(ecl->vector_copy_kernel);
  clReleaseKernel(ecl->inf_norm_kernel);
  clReleaseKernel(ecl->two_norm_kernel);
  clReleaseKernel(ecl->tridiag_sir_kernel);
  clReleaseKernel(ecl->tridiag_kernel);
  clReleaseKernel(ecl->vector_one_norm_tridiag_mat_kernel);
  clReleaseKernel(ecl->expo_trans_kernel);
  clReleaseKernel(ecl->expo_samp_kernel);
  clReleaseKernel(ecl->neg_vals_kernel);
  clReleaseKernel(ecl->expo_sir_trans_kernel);
  clReleaseKernel(ecl->expo_sir_samp_kernel);
  clReleaseProgram(ecl->program);
}

/****************************************************************************/

void expocl_copy_mat(expo_cl* ecl, const expo_type* expo) {
  clReal* A = NULL;
#ifdef DOUBLE
  A = expo->mat;
#else
  A = (clReal*) malloc(3*ecl->mem_size);
  size_t i;
  for (i = 0; i < 3*ecl->size; ++i) A[i] = (clReal) expo->mat[i];
#endif

  ecl->error = clEnqueueWriteBuffer(ecl->queue, ecl->Ad, CL_TRUE, 0, 
                                    3*(ecl->mem_size), A, 0, NULL, NULL);
  checkErr(ecl->error,"Write matrix to device");

#ifndef DOUBLE
  free(A);
#endif
}

/****************************************************************************/

void expocl_copy_index(expo_cl* ecl, const expo_type* expo) {
  if (expo->lookup != NULL) {
    // printf("%d %d %lu\n",expo->num_states,ecl->size,sizeof(int));
    // printf("%d %d\n",expo->lookup[1],expo->lookup[ecl->size*expo->num_states-1]);
    ecl->error = clEnqueueWriteBuffer(ecl->queue, ecl->id, CL_TRUE, 0, 
                                      expo->num_states*ecl->size*sizeof(int), 
                                      expo->lookup, 0, NULL, NULL);
    checkErr(ecl->error,"Write index to device");
  } else {
    fprintf(stderr,"Lookup table not initialized!\n");
  }
}

/****************************************************************************/

void expocl_copy_vectors(expo_cl* ecl, const expo_type* expo) {
  clReal* l = NULL;
  clReal* m = NULL;
  clReal* p = NULL;

#ifdef DOUBLE
  l = expo->lambdaVec;
  m = expo->muFun;
  p = expo->psiFun;
#else
  size_t i;
  if (expo->lambdaVec != NULL) {
    l = (clReal*) malloc(ecl->mem_size);
    for (i = 0; i < ecl->size; ++i) l[i] = (clReal) expo->lambdaVec[i];
  } else {
    printf("Lambda is NULL!\n");
  }

  if (expo->muFun != NULL) {
    m = (clReal*) malloc(ecl->mem_size);
    for (i = 0; i < ecl->size; ++i) m[i] = (clReal) expo->muFun[i];
  } else {
    printf("Mu is NULL!\n");
  }

  if (expo->psiFun != NULL) {
    p = (clReal*) malloc(ecl->mem_size);
    for (i = 0; i < ecl->size; ++i) p[i] = (clReal) expo->psiFun[i];
  } else {
    printf("Psi is NULL!\n");
  }
#endif

  if (l != NULL) {
    ecl->error = clEnqueueWriteBuffer(ecl->queue, ecl->ld, CL_TRUE, 0, 
                                      ecl->mem_size, l, 0, NULL, NULL);
    checkErr(ecl->error,"Write vectors to device");
  }

  if (m != NULL) {
    ecl->error = clEnqueueWriteBuffer(ecl->queue, ecl->md, CL_TRUE, 0, 
                                      ecl->mem_size, m, 0, NULL, NULL);
    checkErr(ecl->error,"Write vectors to device");
  }

  if (p != NULL) {
    ecl->error = clEnqueueWriteBuffer(ecl->queue, ecl->pd, CL_TRUE, 0, 
                                      ecl->mem_size, p, 0, NULL, NULL);
    checkErr(ecl->error,"Write vectors to device");
  }

#ifndef DOUBLE
  if (l != NULL) free(l);
  if (m != NULL) free(m);
  if (p != NULL) free(p);
#endif
}

/****************************************************************************/

void expocl_create_buffers(expo_cl* ecl) {
  cl_mem_flags flags = CL_MEM_READ_WRITE;
  ecl->xd = clCreateBuffer(ecl->context, flags, ecl->mem_size, NULL, &ecl->error);
  checkErr(ecl->error,"x buffer");
  ecl->yd = clCreateBuffer(ecl->context, flags, ecl->mem_size, NULL, &ecl->error);
  checkErr(ecl->error,"y buffer");
  ecl->bd = clCreateBuffer(ecl->context, flags, ecl->mem_size, NULL, &ecl->error);
  checkErr(ecl->error,"b buffer");

  flags = CL_MEM_READ_ONLY;
  ecl->Ad = clCreateBuffer(ecl->context, flags, 3*ecl->mem_size, NULL, &ecl->error);
  checkErr(ecl->error,"A buffer");
}

/****************************************************************************/

void expocl_create_index_buffer(expo_cl* ecl, const expo_type* expo) {
  cl_mem_flags flags = CL_MEM_READ_ONLY;
  ecl->id = clCreateBuffer(ecl->context, flags, 
                           expo->num_states*ecl->size*sizeof(int), 
                           NULL, &ecl->error);
  checkErr(ecl->error,"index buffer");
}

/****************************************************************************/

void expocl_create_vector_buffers(expo_cl* ecl) {
  cl_mem_flags flags = CL_MEM_READ_ONLY;
  ecl->ld = clCreateBuffer(ecl->context, flags, ecl->mem_size, NULL, &ecl->error);
  checkErr(ecl->error,"l buffer");
  ecl->md = clCreateBuffer(ecl->context, flags, ecl->mem_size, NULL, &ecl->error);
  checkErr(ecl->error,"m buffer");
  ecl->pd = clCreateBuffer(ecl->context, flags, ecl->mem_size, NULL, &ecl->error);
  checkErr(ecl->error,"p buffer");
}

/****************************************************************************/

void expocl_release_buffers(expo_cl* ecl) {
  clReleaseMemObject(ecl->xd);
  clReleaseMemObject(ecl->yd);
  clReleaseMemObject(ecl->bd);
  clReleaseMemObject(ecl->Ad);
}

/****************************************************************************/

void expocl_release_vector_buffers(expo_cl* ecl) {
  clReleaseMemObject(ecl->ld);
  clReleaseMemObject(ecl->md);
  clReleaseMemObject(ecl->pd);
}

/****************************************************************************/

void expocl_release_index_buffer(expo_cl* ecl) {
  clReleaseMemObject(ecl->id);
}

/****************************************************************************/

void expocl_copy_b(expo_cl* ecl, clReal* x) {
  ecl->error = clEnqueueWriteBuffer(ecl->queue, ecl->bd, CL_TRUE, 0, 
                                    ecl->mem_size, x, 0, NULL, NULL);
  checkErr(ecl->error,"Write array to device");
}

/****************************************************************************/

void expocl_copy_x(expo_cl* ecl, clReal* x) {
  ecl->error = clEnqueueWriteBuffer(ecl->queue, ecl->xd, CL_TRUE, 0, 
                                    ecl->mem_size, x, 0, NULL, NULL);
  checkErr(ecl->error,"Write array to device");
}

/****************************************************************************/

void expocl_read_b(expo_cl* ecl, clReal* x) {
  clEnqueueReadBuffer(ecl->queue, ecl->xd, CL_TRUE, 0, ecl->mem_size,
                      x, 0, NULL, NULL);
}

/****************************************************************************/

void expocl_read_x(expo_cl* ecl, clReal* x) {
  clEnqueueReadBuffer(ecl->queue, ecl->xd, CL_TRUE, 0, ecl->mem_size,
                      x, 0, NULL, NULL);
}

/****************************************************************************/

void expocl_read_y(expo_cl* ecl, clReal* y) {
  clEnqueueReadBuffer(ecl->queue, ecl->yd, CL_TRUE, 0, ecl->mem_size,
                      y, 0, NULL, NULL);
}

/****************************************************************************/

void expocl_swap_xy(expo_cl* ecl) {
  cl_mem tmp = ecl->xd;
  ecl->xd = ecl->yd;
  ecl->yd = tmp;
}

/****************************************************************************/

cl_event expocl_sis_calc(expo_cl* ecl, expo_type* expo, char trans, clReal alpha) {
  cl_event event;
  // Enqueuing parameters
  // Note that we inform the size of the cl_mem object, 
  // not the size of the memory pointed by it
  ecl->error  = clSetKernelArg(ecl->tridiag_kernel, 0, sizeof(cl_mem), &ecl->Ad);
  ecl->error |= clSetKernelArg(ecl->tridiag_kernel, 1, sizeof(cl_mem), &ecl->xd);
  ecl->error |= clSetKernelArg(ecl->tridiag_kernel, 2, sizeof(cl_mem), &ecl->yd);
  ecl->error |= clSetKernelArg(ecl->tridiag_kernel, 3, sizeof(char), &trans);
  ecl->error |= clSetKernelArg(ecl->tridiag_kernel, 4, sizeof(clReal), &alpha);
  ecl->error |= clSetKernelArg(ecl->tridiag_kernel, 5, sizeof(int), &ecl->size);
  checkErr(ecl->error,"Enque parameters");

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->tridiag_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, &event);
  checkErr(ecl->error,"Enqueue");
  // clFinish(ecl->queue);
  return event;
}

/****************************************************************************/

clReal expocl_inf_norm(expo_cl* ecl, cl_mem* buffer) {
  size_t local_mem_size = ecl->local_ws*sizeof(clReal);
  size_t group_size = ecl->global_ws/ecl->local_ws;

  ecl->error  = clSetKernelArg(ecl->inf_norm_kernel, 0, sizeof(cl_mem), buffer);
  ecl->error |= clSetKernelArg(ecl->inf_norm_kernel, 1, local_mem_size, NULL);
  ecl->error |= clSetKernelArg(ecl->inf_norm_kernel, 2, sizeof(int), &ecl->size);
  ecl->error |= clSetKernelArg(ecl->inf_norm_kernel, 3, sizeof(cl_mem), &ecl->yd);
  checkErr(ecl->error,"Failed to set kernel arguments");

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->inf_norm_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  if (ecl->error != CL_SUCCESS) {
    fprintf(stderr,"Global WS = %lu, Local WS = %lu.\n",ecl->global_ws,ecl->local_ws);
    checkErr(ecl->error,"Enqueue");
  }
  clFinish(ecl->queue);

  // Reading back
  clReal* res = (clReal*) malloc(group_size*sizeof(clReal));
  clEnqueueReadBuffer(ecl->queue, ecl->yd, CL_TRUE, 0, group_size*sizeof(clReal),
                      res, 0, NULL, NULL);

  int i;
  clReal infNorm = res[0];
  for (i = 1; i < group_size; ++i) {
    if (res[i] > infNorm) infNorm = res[i];
  }
  free(res);

  return infNorm;
}

clReal expocl_inf_norm_b(expo_cl* ecl) { return expocl_inf_norm(ecl,&ecl->bd); }
clReal expocl_inf_norm_x(expo_cl* ecl) { return expocl_inf_norm(ecl,&ecl->xd); }

/****************************************************************************/

clReal expocl_two_norm(expo_cl* ecl, cl_mem* buffer) {
  size_t local_mem_size = ecl->local_ws*sizeof(clReal);
  size_t group_size = ecl->global_ws/ecl->local_ws;

  ecl->error  = clSetKernelArg(ecl->two_norm_kernel, 0, sizeof(cl_mem), buffer);
  ecl->error |= clSetKernelArg(ecl->two_norm_kernel, 1, local_mem_size, NULL);
  ecl->error |= clSetKernelArg(ecl->two_norm_kernel, 2, sizeof(int), &ecl->size);
  ecl->error |= clSetKernelArg(ecl->two_norm_kernel, 3, sizeof(cl_mem), &ecl->yd);
  checkErr(ecl->error,"Failed to set kernel arguments");

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->two_norm_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  if (ecl->error != CL_SUCCESS) {
    fprintf(stderr,"Global WS = %lu, Local WS = %lu.\n",ecl->global_ws,ecl->local_ws);
    checkErr(ecl->error,"Enqueue");
  }
  clFinish(ecl->queue);

  // Reading back
  clReal norm = 0.0;
  clReal* res = (clReal*) malloc(group_size*sizeof(clReal));
  clEnqueueReadBuffer(ecl->queue, ecl->yd, CL_TRUE, 0, group_size*sizeof(clReal),
                      res, 0, NULL, NULL);
  int i;
  for (i = 1; i < group_size; ++i) res[0] += res[i];
  norm = sqrt(res[0]);
  free(res);

  return norm;
}

clReal expocl_two_norm_b(expo_cl* ecl) { return expocl_two_norm(ecl,&ecl->bd); }
clReal expocl_two_norm_x(expo_cl* ecl) { return expocl_two_norm(ecl,&ecl->xd); }

/****************************************************************************/

clReal expocl_one_norm_mat(expo_cl* ecl) {
  size_t local_mem_size = ecl->local_ws*sizeof(clReal);

  ecl->error  = clSetKernelArg(ecl->inf_norm_kernel, 0, sizeof(cl_mem), &ecl->Ad);
  ecl->error |= clSetKernelArg(ecl->inf_norm_kernel, 1, local_mem_size, NULL);
  ecl->error |= clSetKernelArg(ecl->inf_norm_kernel, 2, sizeof(int), &ecl->size);
  ecl->error |= clSetKernelArg(ecl->inf_norm_kernel, 3, sizeof(cl_mem), &ecl->yd);
  checkErr(ecl->error,"Failed to set kernel arguments");

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->inf_norm_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  if (ecl->error != CL_SUCCESS) {
    fprintf(stderr,"Global WS = %lu, Local WS = %lu.\n",ecl->global_ws,ecl->local_ws);
    checkErr(ecl->error,"Enqueue");
  }
  clFinish(ecl->queue);

  // Reading back
  clReal res = 0.0;
  clEnqueueReadBuffer(ecl->queue, ecl->yd, CL_TRUE, 0, sizeof(clReal),
                      &res, 0, NULL, NULL);
  return res;
}

/****************************************************************************/

void expocl_copy(expo_cl* ecl, cl_mem* b1, cl_mem* b2) {
  ecl->error  = clSetKernelArg(ecl->vector_copy_kernel, 0, sizeof(cl_mem), b1);
  ecl->error |= clSetKernelArg(ecl->vector_copy_kernel, 1, sizeof(cl_mem), b2);
  ecl->error |= clSetKernelArg(ecl->vector_copy_kernel, 2, sizeof(int), &ecl->size);
  checkErr(ecl->error,"Failed to set kernel arguments");

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->vector_copy_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  checkErr(ecl->error,"Enqueue");
}

void expocl_copy_bx(expo_cl* ecl) { expocl_copy(ecl,&ecl->bd,&ecl->xd); }
void expocl_copy_xb(expo_cl* ecl) { expocl_copy(ecl,&ecl->xd,&ecl->bd); }
void expocl_copy_xy(expo_cl* ecl) { expocl_copy(ecl,&ecl->xd,&ecl->yd); }

/****************************************************************************/

void expocl_smult_b(expo_cl* ecl, clReal alpha) {
  ecl->error  = clSetKernelArg(ecl->vector_smult_kernel, 0, sizeof(cl_mem), &ecl->bd);
  ecl->error |= clSetKernelArg(ecl->vector_smult_kernel, 1, sizeof(clReal), &alpha);
  ecl->error |= clSetKernelArg(ecl->vector_smult_kernel, 2, sizeof(int), &ecl->size);

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->vector_smult_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  checkErr(ecl->error,"Enqueue");
}

/****************************************************************************/

void expocl_add(expo_cl* ecl, cl_mem* u, cl_mem* v) {
  ecl->error  = clSetKernelArg(ecl->vector_add_kernel, 0, sizeof(cl_mem), u);
  ecl->error |= clSetKernelArg(ecl->vector_add_kernel, 1, sizeof(cl_mem), v);
  ecl->error |= clSetKernelArg(ecl->vector_add_kernel, 2, sizeof(int), &ecl->size);

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->vector_add_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  checkErr(ecl->error,"Enqueue");
}

void expocl_add_bx(expo_cl* ecl) { expocl_add(ecl,&ecl->bd,&ecl->xd); }
void expocl_add_by(expo_cl* ecl) { expocl_add(ecl,&ecl->bd,&ecl->yd); }

/****************************************************************************/

void expocl_trans_b(expo_cl* ecl, const expo_type* expo) {
  ecl->error  = clSetKernelArg(ecl->expo_trans_kernel, 0, sizeof(cl_mem), &ecl->bd);
  ecl->error |= clSetKernelArg(ecl->expo_trans_kernel, 1, sizeof(cl_mem), &ecl->xd);
  ecl->error |= clSetKernelArg(ecl->expo_trans_kernel, 2, sizeof(cl_mem), &ecl->ld);
  ecl->error |= clSetKernelArg(ecl->expo_trans_kernel, 3, sizeof(int), &ecl->size);

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->expo_trans_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  checkErr(ecl->error,"Enqueue");
}

/****************************************************************************/

void expocl_samp_b(expo_cl* ecl, const expo_type* expo) {
  ecl->error  = clSetKernelArg(ecl->expo_samp_kernel, 0, sizeof(cl_mem), &ecl->bd);
  ecl->error |= clSetKernelArg(ecl->expo_samp_kernel, 1, sizeof(cl_mem), &ecl->xd);
  ecl->error |= clSetKernelArg(ecl->expo_samp_kernel, 2, sizeof(cl_mem), &ecl->pd);
  ecl->error |= clSetKernelArg(ecl->expo_samp_kernel, 3, sizeof(int), &ecl->size);
  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->expo_samp_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  checkErr(ecl->error,"Enqueue");
}

/****************************************************************************/

int expocl_neg_vals(expo_cl* ecl) {
  size_t local_mem_size = ecl->local_ws*sizeof(int);
  size_t group_size = ecl->global_ws/ecl->local_ws;

  ecl->error  = clSetKernelArg(ecl->neg_vals_kernel, 0, sizeof(cl_mem), &ecl->bd);
  ecl->error  = clSetKernelArg(ecl->neg_vals_kernel, 1, sizeof(cl_mem), &ecl->xd);
  ecl->error |= clSetKernelArg(ecl->neg_vals_kernel, 2, local_mem_size, NULL);
  ecl->error |= clSetKernelArg(ecl->neg_vals_kernel, 3, sizeof(int), &ecl->size);
  ecl->error |= clSetKernelArg(ecl->neg_vals_kernel, 4, sizeof(cl_mem), &ecl->yd);
  checkErr(ecl->error,"Failed to set kernel arguments");

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->neg_vals_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  if (ecl->error != CL_SUCCESS) {
    fprintf(stderr,"Global WS = %lu, Local WS = %lu.\n",ecl->global_ws,ecl->local_ws);
    checkErr(ecl->error,"Enqueue");
  }
  clFinish(ecl->queue);

  // Reading back
  int neg_vals = 0;
  clReal* res = (clReal*) malloc(group_size*sizeof(clReal));
  clEnqueueReadBuffer(ecl->queue, ecl->yd, CL_TRUE, 0, group_size*sizeof(clReal),
                      res, 0, NULL, NULL);
  int i;
  for (i = 1; i < group_size; ++i) res[0] += res[i];
  neg_vals = res[0];
  free(res);

  return neg_vals;
}

/****************************************************************************/

cl_event expocl_sir_calc(expo_cl* ecl, expo_type* expo, char trans, clReal alpha) {
  int N = expo->N_max;
  cl_event event;

  ecl->error  = clSetKernelArg(ecl->tridiag_sir_kernel, 0, sizeof(cl_mem), &ecl->Ad);
  ecl->error |= clSetKernelArg(ecl->tridiag_sir_kernel, 1, sizeof(cl_mem), &ecl->xd);
  ecl->error |= clSetKernelArg(ecl->tridiag_sir_kernel, 2, sizeof(cl_mem), &ecl->id);
  ecl->error |= clSetKernelArg(ecl->tridiag_sir_kernel, 3, sizeof(cl_mem), &ecl->yd);
  ecl->error |= clSetKernelArg(ecl->tridiag_sir_kernel, 4, sizeof(char), &trans);
  ecl->error |= clSetKernelArg(ecl->tridiag_sir_kernel, 5, sizeof(clReal), &alpha);
  ecl->error |= clSetKernelArg(ecl->tridiag_sir_kernel, 6, sizeof(int), &N);
  ecl->error |= clSetKernelArg(ecl->tridiag_sir_kernel, 7, sizeof(int), &ecl->size);
  checkErr(ecl->error,"Enque parameters");

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->tridiag_sir_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, &event);
  checkErr(ecl->error,"Enqueue");
  return event;
}

/****************************************************************************/

void expocl_sir_trans_b(expo_cl* ecl, const expo_type* expo) {
  ecl->error  = clSetKernelArg(ecl->expo_sir_trans_kernel, 0, sizeof(cl_mem), &ecl->bd);
  ecl->error |= clSetKernelArg(ecl->expo_sir_trans_kernel, 1, sizeof(cl_mem), &ecl->xd);
  ecl->error |= clSetKernelArg(ecl->expo_sir_trans_kernel, 2, sizeof(cl_mem), &ecl->id);
  ecl->error |= clSetKernelArg(ecl->expo_sir_trans_kernel, 3, sizeof(cl_mem), &ecl->ld);
  ecl->error |= clSetKernelArg(ecl->expo_sir_trans_kernel, 4, sizeof(int), &expo->N_max);
  ecl->error |= clSetKernelArg(ecl->expo_sir_trans_kernel, 5, sizeof(int), &ecl->size);

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->expo_sir_trans_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  checkErr(ecl->error,"Enqueue");
}

/****************************************************************************/

void expocl_sir_samp_b(expo_cl* ecl, const expo_type* expo) {
  ecl->error  = clSetKernelArg(ecl->expo_sir_samp_kernel, 0, sizeof(cl_mem), &ecl->bd);
  ecl->error |= clSetKernelArg(ecl->expo_sir_samp_kernel, 1, sizeof(cl_mem), &ecl->xd);
  ecl->error |= clSetKernelArg(ecl->expo_sir_samp_kernel, 2, sizeof(cl_mem), &ecl->id);
  ecl->error |= clSetKernelArg(ecl->expo_sir_samp_kernel, 3, sizeof(cl_mem), &ecl->pd);
  ecl->error |= clSetKernelArg(ecl->expo_sir_samp_kernel, 4, sizeof(int), &expo->N_max);
  ecl->error |= clSetKernelArg(ecl->expo_sir_samp_kernel, 5, sizeof(int), &ecl->size);

  ecl->error = clEnqueueNDRangeKernel(ecl->queue, ecl->expo_sir_samp_kernel, 1, NULL, 
                                      &ecl->global_ws, &ecl->local_ws, 
                                      0, NULL, NULL);
  checkErr(ecl->error,"Enqueue");
}


