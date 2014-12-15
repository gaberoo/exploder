#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <OpenCL/opencl.h>
#include "util.h"

#include "tridiag_kernel_str.h"

typedef float clReal;

void checkErr(cl_int err, const char * name) {
  if (err != CL_SUCCESS) {
    fprintf(stderr,"ERROR: %s (%d)\n",name,err);
    exit(EXIT_FAILURE);
  }
}

size_t shrRoundUp(size_t localWorkSize, size_t numItems) {
  size_t result = localWorkSize;
  while (result < numItems) result += localWorkSize;
  checkErr((result >= numItems) && ((result % localWorkSize) == 0) ? CL_SUCCESS : -1,
           "invalid post-condition");
  return result;
}

int main(void) {
  cl_platform_id   platform;
  cl_command_queue queue;
  cl_device_id     device;
  cl_program       program;
  cl_kernel        kernel;

  cl_device_type   dtype = CL_DEVICE_TYPE_GPU;

  cl_int error = oclGetPlatformID(&platform);
  checkErr(error,"Platform ID");

  // Device
  error = clGetDeviceIDs(platform, dtype, 1, &device, NULL);
  checkErr(error,"Device ID");

  // Context
  cl_context context = clCreateContext(0, 1, &device, NULL, NULL, &error);
  checkErr(error,"Create context");

  // Command-queue
  queue = clCreateCommandQueue(context, device, 0, &error);
  checkErr(error,"Create command queue");

  size_t N = 10000000;

  size_t local_ws = 512;
  size_t global_ws = shrRoundUp(local_ws, N);

  int size = N;
  int mem_size = size*sizeof(clReal);

  const char filename[] = "expmv.cl";
  char *kernel_str;
  long len;
  FILE *kernel_file = fopen(filename, "rb");
  fseek(kernel_file,0,SEEK_END);
  len = ftell(kernel_file);
  rewind(kernel_file);
  kernel_str = malloc(len*sizeof(char));
  fread(kernel_str,sizeof(char),len,kernel_file);
  fclose(kernel_file);

  const char* kstr = kernel_str;
  program = clCreateProgramWithSource(context, 1, &kstr, NULL, &error);
  checkErr(error,"Create program");

  // Builds the program
  error = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
  checkErr(error,"Build program");

  // Shows the log
  char* build_log;
  size_t log_size;

  // First call to know the proper size
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

  // Second call to get the log
  build_log = (char*) malloc((log_size+1)*sizeof(char));
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, build_log, NULL);
  build_log[log_size] = '\0';
  free(build_log);

  printf("Build successful.\n");

  // Extracting the kernel
  kernel = clCreateKernel(program, "inf_norm_kernel", &error);
  checkErr(error,"Create kernel");

  clReal* x = (clReal*) malloc(N*sizeof(clReal));
  clReal* y = (clReal*) malloc(N*sizeof(clReal));

  int i;
  for (i = 0; i < N; ++i) { 
    x[i] = 100*exp(-(i-20.0)*(i-20.0)); 
    y[i] = 0.0; 
  }

  cl_mem_flags flags = CL_MEM_READ_WRITE;
  cl_mem xd = clCreateBuffer(context, flags, mem_size, NULL, &error);
  checkErr(error,"Create buffer");

  error = clEnqueueWriteBuffer(queue, xd, CL_TRUE, 0, 
                               mem_size, x, 0, NULL, NULL);
  checkErr(error,"Write array to device");

  flags = CL_MEM_READ_WRITE;
  cl_mem yd = clCreateBuffer(context, flags, mem_size, NULL, NULL);

  // Enqueuing parameters
  // Note that we inform the size of the cl_mem object, 
  // not the size of the memory pointed by it
  error  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &xd);
  error |= clSetKernelArg(kernel, 1, local_ws*sizeof(clReal), NULL);
  error |= clSetKernelArg(kernel, 2, sizeof(size_t), &N);
  error |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &yd);
  checkErr(error,"Failed to set kernel arguments");

  error = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, 
                                 &global_ws, &local_ws, 
                                 0, NULL, NULL);
  checkErr(error,"Enqueue");

  // Reading back
  clEnqueueReadBuffer(queue, yd, CL_TRUE, 0, mem_size, y, 0, NULL, NULL);

  i = 0;
  // for (i = 0; i < N; ++i) 
  printf("%g %g\n",x[i],y[i]);

  clReleaseMemObject(xd);
  clReleaseMemObject(yd);

  free(x);
  free(y);

  clReleaseCommandQueue(queue);
  clReleaseContext(context);

  return 0;
}

