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

////////////////////////////////////////////////////////////////////////////////

// naive, embarassingly parallel CA update
__global__ void kernelCAUpdateNaive(uint8_t *in, uint8_t *out) {
    const uint8_t B = caParams.B;
    const uint8_t S = caParams.S;
    
    
}

// more sophisticated worklist-based approach
// may or may not be faster depending on the state
__global__ void kernelCAUpdateWorkList() {

}

////////////////////////////////////////////////////////////////////////////////

void CASim::Rule::parse(const char *s) {
    bool ok = false;
    if (!strcmp(s, "wire")) {
        type = WireWorld;
        ok = true;
    } else if (s[0] == 'B') {
        // B[1-8]+S[0-8]+
        type = LifeLike;
        B = 0;
        S = 0;
        int i = 1;
        while (s[i] >= '1' && s[i] <= '8')
            B |= (1 << (s[i++] - '0'));
        if (s[i] == 'S') {
            ++i;
            while (s[i] >= '1' && s[i] <= '8')
                S |= (1 << (s[i++] - '0'));
            ok = s[i] == '\0';
        }
    }
    
    if (!ok) {
        fprintf(stderr, "Invalid rule, falling back to Life\n");
        type = Rule::LifeLike;
        B = 1 << 3;
        S = (1 << 2) | (1 << 3);
    }
}

const char* CASim::Rule::toString() const {
    static char buffer[20];
    if (type == WireWorld)
        return "wire";
    char *p = buffer;
    *p++ = 'B';
    for (int i = 1; i <= 8; ++i)
        if (B & (1 << i))
            *p++ = '0' + i;
    *p++ = 'S';
    for (int i = 0; i <= 8; ++i)
        if (S & (1 << i))
            *p++ = '0' + i;
    *p = '\0';
    return buffer;
}

////////////////////////////////////////////////////////////////////////////////

CASim::CASim(const char *rule) {
    this->rule.parse(rule);
    
    generation = 0;
    width = height = 0;
    cell0 = cell1 = 0;
    
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);

    printf("Initializing CUDA for CASim\n");
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
    fprintf(stderr, "Running %i steps of %s\n", n, rule.toString());
    
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
