//
//  main.cpp
//  shinycl
//
//  Created by Cianciosa, Mark R. on 10/28/15.
//  Copyright © 2015 Cianciosa, Mark R. All rights reserved.
//

#include <iostream>
#include <OpenCL/OpenCL.h>
#include <algorithm>
#include <cmath>

#include "shiny.cl.h"

int main(int argc, const char * argv[]) {
    
    dispatch_queue_t gpu_queue = gcl_create_dispatch_queue(CL_DEVICE_TYPE_GPU, NULL);
    
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
    
    void *t_mem_object = gcl_malloc(sizeof(cl_float)*nx*ny, NULL, CL_MEM_WRITE_ONLY);
    void *t_next_mem_object = gcl_malloc(sizeof(cl_float)*nx*ny, NULL, CL_MEM_WRITE_ONLY);
    void *t_analytic_mem_object = gcl_malloc(sizeof(cl_float)*nx*ny, NULL, CL_MEM_WRITE_ONLY);
    
    __block size_t work_group_size;
    dispatch_sync(gpu_queue, ^{
        gcl_get_kernel_block_workgroup_info(init_kernel, CL_KERNEL_WORK_GROUP_SIZE, sizeof(work_group_size), &work_group_size, NULL);
    });
    
    size_t global_x_size = work_group_size;
    size_t global_y_size = work_group_size;
    
    while (nx > global_x_size) {
        global_x_size += work_group_size;
    }
    while (ny > global_y_size) {
        global_y_size += work_group_size;
    }
    
    cl_ndrange range = {
        2,
        {0, 0, 0},
        {global_x_size, global_y_size, 0},
        {global_x_size/work_group_size, global_y_size/work_group_size, 0}
    };
    
    dispatch_sync(gpu_queue, ^{
        init_kernel(&range, (cl_float *)t_mem_object, (cl_float *)t_next_mem_object, nx, ny);
    });

    const cl_float dt = 0.9/8.0*(dx*dx + dy*dy)/std::max(fabsf(kper), fabsf(kpar));
    
    std::cout << std::endl << "dt = " << dt << std::endl << std::endl;
    
    for (size_t time = 0; time <= end_time; time++) {
        dispatch_async(gpu_queue, ^{
            stepTimeGunterSym_kernel(&range, (cl_float *)t_mem_object, (cl_float *)t_next_mem_object,
                                     nx, ny, xmin, ymin, dt, dx, dy, kper, kpar);
            update_kernel(&range, (cl_float *)t_mem_object, (cl_float *)t_next_mem_object, nx, ny);
        });
    }

    dispatch_async(gpu_queue, ^{
        const cl_float final_time = dt*end_time;
        analytic_kernel(&range, (cl_float *)t_analytic_mem_object,
                        nx, ny, xmin, ymin, final_time, dx, dy, kper);
    });
    
    cl_float *t_host_mem = new cl_float[nx*ny];
    cl_float *t_analytic_host_mem = new cl_float[nx*ny];
    dispatch_sync(gpu_queue, ^{
        gcl_memcpy(t_analytic_host_mem, t_analytic_mem_object, sizeof(cl_float)*nx*ny);
        gcl_memcpy(t_host_mem, t_mem_object, sizeof(cl_float)*nx*ny);
    });
    
    for (size_t i = 0, e = nx*ny; i < e; i++) {
        std::cout << t_host_mem[i] << " " << t_analytic_host_mem[i] << std::endl;
    }
    
//    for (size_t i = 0; i < nx; i++) {
//        for (size_t j = 0; j < ny; j++) {
//            size_t k = i*ny + j;
//            std::cout << t_host_mem[k] << " ";
//        }
//        std::cout << std::endl;
//    }
    
    gcl_free(t_mem_object);
    gcl_free(t_next_mem_object);
    
    dispatch_release(gpu_queue);
    
    return 0;
}
