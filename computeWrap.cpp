
#include "cl.hpp"

 
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
 
#undef max
#undef min

#define real double
int maxWorkGroupSize = 256;
size_t maxThreadsSize = 256;
cl::CommandQueue queue;

std::string load_file(const std::string &file_name, size_t max_size = 0x100000)
{
	FILE *fp = fopen(file_name.c_str(), "rb");
	if (!fp)
	{
		// print some error or throw exception here
		return std::string();
	}
	char *source = new char[max_size];
	memset(source, 0, max_size);
	size_t source_size = fread(source, 1, max_size, fp);
	fclose(fp);
	if (!source_size)
	{
		delete[] source;
		// print some error or throw exception here
		return std::string();
	}
	std::string result(source);
	delete[] source;
	return result;
}

cl::Context getContext()
{
	std::vector<cl::Platform> platforms;
	auto err = cl::Platform::get(&platforms);
	printf("number of platforms: %d\n", platforms.size());
	if (platforms.size() == 0) {
		printf("Platform size 0\n");
	}

	cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0 };
	return   cl::Context(CL_DEVICE_TYPE_GPU, properties);
}

cl::Program build_Program(cl::Context context )
{ 

 
	auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
 
	printf("number of devices %d\n", devices.size());
	  devices[0];

	  //devices[0].getInfo(CL_DEVICE_VENDOR, &maxWorkGroupSize);
	  maxWorkGroupSize = 256;
	  printf("number of max work items size  %d\n", maxWorkGroupSize);

	  maxThreadsSize = 1;
	  devices[0].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE , &maxThreadsSize);
	  printf("number of max processors  %d\n", maxThreadsSize);

	  std::string  extensions;
	  devices[0].getInfo(CL_DEVICE_EXTENSIONS, &extensions);
	  std::cout << extensions << std::endl;


	const char * pcode = R"V(
  kernel void vectorAdd(global const int *inputA, global const int *inputB, global int *output ,    const int x)
  {
     output[get_global_id(0)] = inputA[get_global_id(0)]+ inputB[get_global_id(0)] + x  ;
  }		;

 )V";

	std::string fullSource = load_file("likeHood.cl");
	
	cl::Program vectorWrapper(context , cl::STRING_CLASS(fullSource.c_str()), false);

 
	std::vector<cl::Program> programs;

	cl_int  err_compile = vectorWrapper.compile();

	if (err_compile != CL_SUCCESS)
	{

		printf("done building program\n");
		std::cout << "Build Status: " << vectorWrapper.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(devices[0]) << std::endl;
		std::cout << "Build Options:\t" << vectorWrapper.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(devices[0]) << std::endl;
		std::cout << "Build Log:\t " << vectorWrapper.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << std::endl;
		return 	cl::Program();
	}


	programs.push_back(vectorWrapper);
	cl::Program gpuProgram = cl::linkProgram(programs);

	printf("done building program\n");
	cl_build_status link_status = gpuProgram.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(devices[0]);
	if (link_status != CL_SUCCESS)
	{
		std::cout << "Build Status: " << gpuProgram.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(devices[0]) << std::endl;
		std::cout << "Build Options:\t" << gpuProgram.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(devices[0]) << std::endl;
		std::cout << "Build Log:\t " << gpuProgram.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << std::endl;
		 return 	cl::Program();

	}

	printf("Build Program\n");
	cl_int err;
	 // queue = cl::CommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &err);
	queue = cl::CommandQueue(context, devices[0], 0 , &err);

	return gpuProgram;

	 
}


 

struct LikelihoodParameter
{
	real alpha, T, ap, tp, tau1, tau2;
	real result;
	LikelihoodParameter(real _alpha, real _T, real _ap, real _tp, real _tau1, real _tau2);
};


 

 void computeParams( std::vector<LikelihoodParameter> &params)
{
	 static cl::Context context = getContext();
	 static cl::Program gpuProgram = build_Program(context);

	 auto LikeHoodKernel =
		 cl::make_kernel<
		 cl::Buffer&,
		 cl::Buffer&,
		 int
		 >(gpuProgram, "LikelihoodList");




 
	 int numParams  = params.size();
 
	 if (numParams > maxWorkGroupSize)
	 {
		 while (params.size() % maxWorkGroupSize != 0)
		 {
			 params.emplace_back(0, 0, 0, 0, 0, 0);
		 }
	 }
	 int actualWorkGroupSize = std::min(numParams, maxWorkGroupSize);

	 int arraySize = params.size();
	 //std::vector<double> output(numParams, 0xdeadbeef);
	 std::vector<real> output(arraySize, 1e-99);
	 cl::Buffer inputABuffer(queue , begin(params), end(params), true);
	 cl::Buffer outputBuffer(queue, begin(output), end(output), false);

	 cl::Kernel kernel(gpuProgram, "LikelihoodList");

	 cl::Event event;

	  
	 auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
	// printf("number of devices %d\n", devices.size());
	 
	 

	 //LikeHoodKernel( cl::EnqueueArgs(event, cl::NDRange(arraySize), cl::NDRange(actualWorkGroupSize)), inputABuffer, outputBuffer, numParams );

	 kernel.setArg(0, inputABuffer);
	 kernel.setArg(1, outputBuffer);
	 kernel.setArg(2, numParams);


	 cl::Event evt;
	cl_int  err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(arraySize), cl::NDRange(actualWorkGroupSize) );
	 
	queue.finish();

	 //auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
	// cl:: CommandQueue queue =cl::CommandQueue(context, devices[0]);
	 //cl::Kernel kernel(gpuProgram, "LikelihoodList");
	// queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(arraySize),  cl::NDRange(actualWorkGroupSize));

	cl::copy(queue,outputBuffer, begin(output), end(output));

	//queue.finish();

	 
	 for (int i = 0; i < numParams; ++i) 
	 {
		 double rr = output[i];
		 //printf(" %15.12g \n", rr );

		 params[i].result = output[i];
	 }
 
	 return ;

}