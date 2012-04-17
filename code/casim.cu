/****************************************************************************
** Copyright (c) 2012 Hugues Bruant <hugues@cmu.edu>
** All rights reserved.
**
** This file may be used under the terms of the GNU General Public License
** version 3 as published by the Free Software Foundation.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
****************************************************************************/

#include "casim.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/scan.h>

struct CAParams {
    unsigned short type;
    uint8_t B, S;
    
    unsigned int width;
    unsigned int height;
    uint8_t *cell0;
    uint8_t *cell1;
};

__constant__ CAParams caParams;


CASim::CASim(const char *rule) {
    // TODO: parse rule
    //rule.type = ;
    // rule.B = ;
    // rule.S = ;
    
    generation = 0;
    width = height = 0;
    cell0 = cell1 = 0;
    
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);

    printf("Initializing CUDA for CudaRenderer\n");
    printf("Found %d CUDA devices\n", deviceCount);

    for (int i=0; i<deviceCount; i++) {
        cudaDeviceProp deviceProps;
        cudaGetDeviceProperties(&deviceProps, i);
        printf("Device %d: %s\n", i, deviceProps.name);
        printf("   SMs:        %d\n", deviceProps.multiProcessorCount);
        printf("   Global mem: %.0f MB\n", static_cast<float>(deviceProps.totalGlobalMem) / (1024 * 1024));
        printf("   CUDA Cap:   %d.%d\n", deviceProps.major, deviceProps.minor);
    }
}

CASim::~CASim() {
    if (cell0) cudaFree(cell0);
    if (cell1) cudaFree(cell1);
}

void CASim::step(int n) {
    for (int i = 0; i < n; ++i) {
        // TODO : kernel launch
        ++generation;
    }
}

bool CASim::setCells(unsigned int width, unsigned int height, uint8_t max,
                     const uint8_t *cells) {
    // check that the input respects the maximum number of states of the CA
    if ((rule.type == Rule::WireWorld && max > 3) ||
        (rule.type == Rule::LifeLike && max > 1)) {
        return false;
    }
    
    if (cell0) cudaFree(cell0);
    if (cell1) cudaFree(cell1);
    
    this->width = width;
    this->height = height;
    
    generation = 0;
    
    // Alloc double buffers
    size_t sz = sizeof(uint8_t) * width * height;
    cudaMalloc(&cell0, sz);
    cudaMalloc(&cell1, sz);
    
    // set global CA params
    CAParams params;
    params.width = width;
    params.height = height;
    params.type = rule.type;
    params.B = rule.B;
    params.S = rule.S;
    params.cell0 = cell0;
    params.cell1 = cell1;
    cudaMemcpyToSymbol(caParams, &params, sizeof(CAParams));
    
    // copy initial state to first buffer
    // for now simple 1:1 mapping, each byte being one cell but may change in
    // the future to allow further optimization
    cudaMemcpy(cell0, cells, sz, cudaMemcpyHostToDevice);
    
    return true;
}

void CASim::getCells(uint8_t *cells) {
    // copy current state from appropriate buffer
    // for now simple 1:1 mapping, each byte being one cell but may change in
    // the future to allow further optimization
    size_t sz = sizeof(uint8_t) * width * height;
    cudaMemcpy(cells, (generation & 1) ? cell1 : cell0, sz, cudaMemcpyDeviceToHost);
}
