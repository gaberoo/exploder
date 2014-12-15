#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

void checkErr(cl_int err, const char * name);
size_t shrRoundUp(size_t localWorkSize, size_t numItems);
void print_info();
int load_kernel_source(const char* filename, char** str, char type);

cl_int oclGetPlatformID(cl_platform_id* clSelectedPlatformID);

