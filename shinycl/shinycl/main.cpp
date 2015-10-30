//
//  main.cpp
//  shinycl
//
//  Created by Cianciosa, Mark R. on 10/28/15.
//  Copyright © 2015 Cianciosa, Mark R. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

#if defined (__APPLE__) || defined (MACOSX)
#include <OpenCL/OpenCL.h>
#else
#include <CL/opencl.h>
#endif

int main(int argc, const char * argv[]) {
    
    cl_int error = 0;
    
    const size_t nx = 40;
    const size_t ny = 40;
    
    const cl_float xmin = -0.5;
    const cl_float xmax =  0.5;
    const cl_float ymin = -0.5;
    const cl_float ymax =  0.5;
    
    const cl_float dx = (xmax - xmin)/(nx - 1);
    const cl_float dy = (ymax - ymin)/(ny - 1);
    
    const cl_float kper = 1.0;
    const cl_float kpar = 1.0E3;
    
    const size_t end_time = 1000;
    
    cl_uint num_platforms;
    cl_platform_id *platform_ids;
    error = clGetPlatformIDs(0, NULL, &num_platforms);
    if (num_platforms != 0) {
        platform_ids = new cl_platform_id[num_platforms];
        error |= clGetPlatformIDs(num_platforms, platform_ids, NULL);
    } else {
        platform_ids = new cl_platform_id[1];
        platform_ids[0] = NULL;
    }
    
    for (cl_uint i = 0; i < num_platforms; i++) {
        size_t platform_info_size;
        error |= clGetPlatformInfo(platform_ids[i], CL_PLATFORM_NAME, 0, NULL, &platform_info_size);
        
        char *info = new char[platform_info_size];
        error |= clGetPlatformInfo(platform_ids[i], CL_PLATFORM_NAME, platform_info_size, info, NULL);
       
        std::cout << info << std::endl;
        delete [] info;
    }
    
    cl_uint num_devices;
    cl_device_id *device_ids;
    error |= clGetDeviceIDs(platform_ids[0], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
    if (num_devices == 0) {
        error |= clGetDeviceIDs(NULL, CL_DEVICE_TYPE_CPU, 0, NULL, &num_devices);
        if (num_devices == 0) {
            std::cout << "Error: Could not find any valid OpenCL devices." << std::endl;
            exit(1);
        }
        
        device_ids = new cl_device_id[num_devices];
        error |= clGetDeviceIDs(NULL, CL_DEVICE_TYPE_CPU, num_devices, device_ids, NULL);
    } else {
        device_ids = new cl_device_id[num_devices];
        error |= clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, num_devices, device_ids, NULL);
    }
    
    for (size_t i = 0; i < num_devices; i++) {
        size_t device_info_size;
        error |= clGetDeviceInfo(device_ids[0], CL_DEVICE_NAME, 0, NULL, &device_info_size);
        
        char *info = new char[device_info_size];
        error |= clGetDeviceInfo(device_ids[0], CL_DEVICE_NAME, device_info_size, info, NULL);
        
        std::cout << info << std::endl;
        delete [] info;
    }
    
    if (error) {
        std::cout << "Failed to get platforms and devices" << std::endl;
        exit(1);
    }
    
    // Apparently you can set up a contect with multiple devices. This may be usefull for running all
    // three methods at once.
    cl_context context = clCreateContext(NULL, num_devices, device_ids, NULL, NULL, &error);
    cl_command_queue queue0 = clCreateCommandQueue(context, device_ids[0], 0, &error);
    
    std::ifstream cl_source("shiny.cl", std::fstream::in);
    std::filebuf *fileBuffer = cl_source.rdbuf();
    const size_t cl_source_size = fileBuffer->pubseekoff(0, cl_source.end, cl_source.in);
    fileBuffer->pubseekpos(0, cl_source.in);
    
    char *source_buffer =  new char[cl_source_size];
    
    fileBuffer->sgetn(source_buffer, cl_source_size);
    cl_source.close();
    
    cl_program program = clCreateProgramWithSource(context, 1, const_cast<const char **> (&source_buffer), &cl_source_size, &error);
    clBuildProgram(program, num_devices, device_ids, NULL, NULL, NULL);
    delete [] source_buffer;
    
    size_t build_info_size;
    error |= clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &build_info_size);
    if (build_info_size > 0) {
        char *info = new char[build_info_size];
        error |= clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, build_info_size, info, NULL);
        
        std::cout << info << std::endl;
        
        delete [] info;
    }
    
    cl_kernel init_kernel = clCreateKernel(program, "init", &error);
    cl_kernel stepTimeGunterSym_kernel = clCreateKernel(program, "stepTimeGunterSym", &error);
    cl_kernel update_kernel = clCreateKernel(program, "update", &error);
    cl_kernel analytic_kernel = clCreateKernel(program, "analytic", &error);
    
    if (error) {
        std::cout << "Failed to get create kernels" << std::endl;
        exit(1);
    }
    
    const size_t buffersize = sizeof(cl_float)*nx*ny;
    
    cl_mem t_mem_object = clCreateBuffer(context, CL_MEM_WRITE_ONLY, buffersize, NULL, &error);
    cl_mem t_next_mem_object = clCreateBuffer(context, CL_MEM_WRITE_ONLY, buffersize, NULL, &error);
    cl_mem t_analytic_mem_object = clCreateBuffer(context, CL_MEM_WRITE_ONLY, buffersize, NULL, &error);

    if (error) {
        std::cout << "Failed to create memory objects" << std::endl;
        exit(1);
    }
    
    size_t work_group_size;
    error = clGetKernelWorkGroupInfo(init_kernel, device_ids[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(work_group_size), &work_group_size, NULL);
    
    size_t global_x_size = work_group_size;
    size_t global_y_size = work_group_size;
    
    while (nx > global_x_size) {
        global_x_size += work_group_size;
    }
    while (ny > global_y_size) {
        global_y_size += work_group_size;
    }
    
    const size_t global_work_size[2] = {global_x_size, global_y_size};
    const size_t local_work_size[2] = {1, 1};
    
    error  = clSetKernelArg(init_kernel, 0, sizeof(cl_mem), &t_mem_object);
    error |= clSetKernelArg(init_kernel, 1, sizeof(cl_mem), &t_next_mem_object);
    error |= clSetKernelArg(init_kernel, 2, sizeof(size_t), &nx);
    error |= clSetKernelArg(init_kernel, 3, sizeof(size_t), &ny);
    error |= clEnqueueNDRangeKernel(queue0, init_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL);

    if (error) {
        std::cout << CL_INVALID_PROGRAM_EXECUTABLE << std::endl;
        std::cout << CL_INVALID_COMMAND_QUEUE << std::endl;
        std::cout << CL_INVALID_KERNEL << std::endl;
        std::cout << CL_INVALID_CONTEXT << std::endl;
        std::cout << CL_INVALID_KERNEL_ARGS << std::endl;
        std::cout << CL_INVALID_WORK_DIMENSION << std::endl;
        std::cout << CL_INVALID_WORK_GROUP_SIZE << std::endl;
        
        std::cout << "Failed to enqueue init_kernel" << std::endl;
        exit(1);
    }
    
    const cl_float dt = 0.9/8.0*(dx*dx + dy*dy)/std::max(fabsf(kper), fabsf(kpar));
    
    std::cout << std::endl << "dt = " << dt << std::endl << std::endl;
    
    error  = clSetKernelArg(stepTimeGunterSym_kernel, 0, sizeof(cl_mem), &t_mem_object);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 1, sizeof(cl_mem), &t_next_mem_object);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 2, sizeof(size_t), &nx);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 3, sizeof(size_t), &ny);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 4, sizeof(cl_float), &xmin);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 5, sizeof(cl_float), &ymin);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 6, sizeof(cl_float), &dt);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 7, sizeof(cl_float), &dx);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 8, sizeof(cl_float), &dy);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 9, sizeof(cl_float), &kper);
    error |= clSetKernelArg(stepTimeGunterSym_kernel, 10, sizeof(cl_float), &kpar);
    
    error |= clSetKernelArg(update_kernel, 0, sizeof(cl_mem), &t_mem_object);
    error |= clSetKernelArg(update_kernel, 1, sizeof(cl_mem), &t_next_mem_object);
    error |= clSetKernelArg(update_kernel, 2, sizeof(size_t), &nx);
    error |= clSetKernelArg(update_kernel, 3, sizeof(size_t), &ny);
    
    for (size_t time = 0; time <= end_time; time++) {
        error |= clEnqueueNDRangeKernel(queue0, stepTimeGunterSym_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL);
        error |= clEnqueueNDRangeKernel(queue0, update_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL);
    }
    
    if (error) {
        std::cout << "Failed to enqueue solve symmetric" << std::endl;
        exit(1);
    }
    
    const cl_float final_time = dt*end_time;
    
    error  = clSetKernelArg(analytic_kernel, 0, sizeof(cl_mem), &t_analytic_mem_object);
    error |= clSetKernelArg(analytic_kernel, 1, sizeof(size_t), &nx);
    error |= clSetKernelArg(analytic_kernel, 2, sizeof(size_t), &ny);
    error |= clSetKernelArg(analytic_kernel, 3, sizeof(cl_float), &xmin);
    error |= clSetKernelArg(analytic_kernel, 4, sizeof(cl_float), &ymin);
    error |= clSetKernelArg(analytic_kernel, 5, sizeof(cl_float), &final_time);
    error |= clSetKernelArg(analytic_kernel, 6, sizeof(cl_float), &dx);
    error |= clSetKernelArg(analytic_kernel, 7, sizeof(cl_float), &dy);
    error |= clSetKernelArg(analytic_kernel, 8, sizeof(cl_float), &kper);
    
    error |= clEnqueueNDRangeKernel(queue0, analytic_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL);
    
    if (error) {
        std::cout << "Failed to enqueue analytic_kernel" << std::endl;
        exit(1);
    }
    
    cl_float *t_host_mem = new cl_float[nx*ny];
    cl_float *t_analytic_host_mem = new cl_float[nx*ny];
    
    error = clFinish(queue0);
    
    error |= clEnqueueReadBuffer(queue0, t_mem_object, CL_TRUE, 0, buffersize, t_host_mem, 0, NULL, NULL);
    error |= clEnqueueReadBuffer(queue0, t_analytic_mem_object, CL_TRUE, 0, buffersize, t_analytic_host_mem, 0, NULL, NULL);
    
    if (error) {
        std::cout << "Failed to read data buffers" << std::endl;
        exit(1);
    }

    for (size_t i = 0, e = nx*ny; i < e; i++) {
        std::cout << t_host_mem[i] << " " << t_analytic_host_mem[i] << std::endl;
    }
     
    std::cout << std::endl;
    
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            size_t k = i*ny + j;
            std::cout << t_host_mem[k] << " ";
        }
        std::cout << std::endl;
    }
    
    delete [] platform_ids;
    delete [] device_ids;
    
    delete [] t_host_mem;
    clReleaseMemObject(t_mem_object);
    delete [] t_analytic_host_mem;
    clReleaseMemObject(t_analytic_mem_object);
    clReleaseMemObject(t_next_mem_object);
    
    clRetainContext(context);
    
    return 0;
}
