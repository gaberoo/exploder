#include "expocl.h"

int main() {
  expo_type* expo = expo_type_alloc();

  printf("Initializing device...");
  expo_cl* ecl = expocl_init(expo,CL_DEVICE_TYPE_GPU,128);
  printf("done.\n");

  char type = 's';
#ifdef DOUBLE
  type = 'd';
#endif

  char* kernel_str = NULL;
  load_kernel_source("expmv_kernel.cl",&kernel_str,type);

  const char* kstr = kernel_str;
  // const char* kstr = program_str;
  // const char** kstr = (const char**) kernel_str;

  ecl->program = clCreateProgramWithSource(ecl->context, 1, 
                                           &kstr, NULL, &ecl->error);
  checkErr(ecl->error,"Create program");

  // Builds the program
  ecl->error = clBuildProgram(ecl->program, 1, &ecl->device, NULL, NULL, NULL);

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

  expocl_free(ecl);

  return 0;
}
