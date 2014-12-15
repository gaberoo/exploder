#include "util.h"

/****************************************************************************/

void checkErr(cl_int err, const char * name) {
  if (err != CL_SUCCESS) {
    fprintf(stderr,"ERROR: %s (%d)\n",name,err);
    exit(EXIT_FAILURE);
  }
}

/****************************************************************************/

size_t shrRoundUp(size_t localWorkSize, size_t numItems) {
  size_t result = localWorkSize;
  while (result < numItems) result += localWorkSize;
  checkErr((result >= numItems) && ((result % localWorkSize) == 0) ? CL_SUCCESS : -1,
           "invalid post-condition");
  return result;
}

//NVIDIA's code follows
//license issues probably prevent you from using this, but shouldn't take long
//to reimplement
//////////////////////////////////////////////////////////////////////////////
//! Gets the platform ID for NVIDIA if available, otherwise default to platform 0
//!
//! @return the id 
//! @param clSelectedPlatformID         OpenCL platform ID
//////////////////////////////////////////////////////////////////////////////
cl_int oclGetPlatformID(cl_platform_id* clSelectedPlatformID) {
  char chBuffer[1024];
  cl_uint num_platforms;
  cl_platform_id* clPlatformIDs;
  cl_int ciErrNum;
  *clSelectedPlatformID = NULL;
  cl_uint i = 0;

  // Get OpenCL platform count
  ciErrNum = clGetPlatformIDs (0, NULL, &num_platforms);
  if (ciErrNum != CL_SUCCESS) {
    printf(" Error %i in clGetPlatformIDs Call !!!\n\n", ciErrNum);
    return -1000;
  } else {
    if (num_platforms == 0) {
      printf("No OpenCL platform found!\n\n");
      return -2000;
    } else {
      // if there's a platform or more, make space for ID's
      clPlatformIDs = (cl_platform_id*) malloc(num_platforms * sizeof(cl_platform_id));
      if (clPlatformIDs == NULL) {
        printf("Failed to allocate memory for cl_platform ID's!\n\n");
        return -3000;
      }

      // get platform info for each platform and trap the NVIDIA platform if found
      ciErrNum = clGetPlatformIDs (num_platforms, clPlatformIDs, NULL);
      for (i = 0; i < num_platforms; ++i) {
        ciErrNum = clGetPlatformInfo (clPlatformIDs[i], CL_PLATFORM_NAME, 1024, &chBuffer, NULL);
        if (ciErrNum == CL_SUCCESS) {
          if (strstr(chBuffer, "NVIDIA") != NULL) {
            *clSelectedPlatformID = clPlatformIDs[i];
            break;
          }
        }
      }

      // default to zeroeth platform if NVIDIA not found
      if (*clSelectedPlatformID == NULL) {
        // printf("WARNING: NVIDIA OpenCL platform not found - defaulting to first platform!\n\n");
        // printf("selected platform: %d\n", 0);
        *clSelectedPlatformID = clPlatformIDs[0];
      }

      free(clPlatformIDs);
    }
  }

  return CL_SUCCESS;
}

/****************************************************************************/

int load_kernel_source(const char* filename, char** str, char type) {
  long len;
  FILE *kernel_file = fopen(filename, "rb");
  if (kernel_file == NULL) {
    // printf("Couldn't find kernel file: %s!\n",filename);
    return -1001;
  }

  const char doubleType[] = "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\ntypedef double clReal;\n\n";
  const char  floatType[] = "typedef float  clReal;\n";
  size_t typeLen;

  fseek(kernel_file,0,SEEK_END);
  len = ftell(kernel_file);
  rewind(kernel_file);

  switch (type) {
    case 'd':
      typeLen = strlen(doubleType);
      *str = malloc((len+1+typeLen)*sizeof(char));
      strcpy(*str,doubleType);
      break;
    case 's':
    default:
      typeLen = strlen(floatType);
      *str = malloc((len+1+typeLen)*sizeof(char));
      strcpy(*str,floatType);
      break;
  }

  fread(*str+typeLen,sizeof(char),len,kernel_file);
  (*str)[len+typeLen] = '\0';

  fclose(kernel_file);
  return 1;
}

/****************************************************************************/

void print_info() {
	cl_platform_id platforms[32];
	cl_uint num_platforms;
	char vendor[1024];

	cl_device_id devices[32];
	cl_uint num_devices;

	char deviceName[1024];
	cl_uint numberOfCores;
	cl_long amountOfMemory;
	cl_uint clockFreq;
	cl_ulong maxAlocatableMem;
	cl_ulong localMem;
	cl_bool	available;

  printf("OpenCL Info\n");
  printf("-----------\n\n");

	clGetPlatformIDs (32, platforms, &num_platforms);
	printf("Number of platforms: %u\n", num_platforms);

  unsigned int i;
	for (i = 0; i < num_platforms; ++i) {
		printf("  Platform %u:\n", i);

		clGetPlatformInfo (platforms[i], CL_PLATFORM_VENDOR, sizeof(vendor), vendor, NULL);
		clGetDeviceIDs (platforms[i], CL_DEVICE_TYPE_ALL, sizeof(devices), devices, &num_devices);

		printf("    Platform Vendor:   %s\n", vendor);
		printf("    Number of devices: %u\n\n", num_devices);

    unsigned int j;
		for (j = 0; j < num_devices; j++) {
			clGetDeviceInfo(devices[j], CL_DEVICE_NAME, sizeof(deviceName), deviceName, NULL);
			clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, sizeof(vendor), vendor, NULL);
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(numberOfCores), &numberOfCores, NULL);
			clGetDeviceInfo(devices[j], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(amountOfMemory), &amountOfMemory, NULL);
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clockFreq), &clockFreq, NULL);
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(maxAlocatableMem), &maxAlocatableMem, NULL);
			clGetDeviceInfo(devices[j], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(localMem), &localMem, NULL);
			clGetDeviceInfo(devices[j], CL_DEVICE_AVAILABLE, sizeof(available), &available, NULL);

			printf("    Device %u:\n", j);
			printf("      Name:                    %s\n", deviceName);
			printf("      Vendor:                  %s\n", vendor);
			printf("      Available:               %s\n", available ? "Yes" : "No");
			printf("      Compute Units:           %u\n", numberOfCores);
			printf("      Clock Frequency:         %u MHz\n", clockFreq);
			printf("      Global Memory:           %0.00f MB\n", (double) amountOfMemory/1048576);
			printf("      Max Allocateable Memory: %0.00f MB\n", (double) maxAlocatableMem/1048576);
			printf("      Local Memory:            %u kB\n\n", (unsigned int) localMem);
		}
	}
}
