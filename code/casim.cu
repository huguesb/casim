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
#include "cycleTimer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/scan.h>

struct CARule {
    enum Type {
        LifeLike,
        WireWorld
    };
    
    CARule(const char *s);
    const char* toString() const;
    
    unsigned short type;
    uint16_t B, S; // parameters for Life-like
};

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

enum {
    kMacroCellWidth = 16,
    kMacroCellHeight = 16,
    
};

// naive, embarassingly parallel CA update
template <uint16_t T>
__global__ void kernelCAUpdateNaive(uint8_t *in, size_t ipitch,
                                    uint8_t *out, size_t opitch) {
    __shared__ uint8_t ocells[(kMacroCellHeight+2)*(kMacroCellWidth+8)];
    
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    
    // cooperative load of relevant old cell values to shared memory
    uint8_t *ip = in + row * ipitch + col * 4;
    
    uint8_t *sp = ocells + (threadIdx.y+1) * (kMacroCellWidth+8) + (threadIdx.x+1) * 4;
    
    // load inner cells cooperatively
    if (row < caParams.height && col < caParams.width) {
        // TODO: what if width is not a multiple of block width?
        *reinterpret_cast<uint32_t*>(sp) = *reinterpret_cast<uint32_t*>(ip);
        
        // load left and right boundary cells
        if (threadIdx.x == 0) {
            sp[-1] = blockIdx.x == 0 ? 0 : ip[-1];
        } else if (threadIdx.x == blockDim.x-1) {
            sp[4] = blockIdx.x == gridDim.x-1 ? 0 : ip[4];
        }
        
        // load top and bottom boundary cells
        if (threadIdx.y == 0) {
            sp -= (kMacroCellWidth+8);
            if (blockIdx.y != 0) {
                ip -= ipitch;
                *reinterpret_cast<uint32_t*>(sp) = *reinterpret_cast<uint32_t*>(ip);
                if (threadIdx.x == 0) {
                    sp[-1] = blockIdx.x == 0 ? 0 : ip[-1];
                } else if (threadIdx.x == blockDim.x-1) {
                    sp[4] = blockIdx.x == gridDim.x-1 ? 0 : ip[4];
                }
            } else {
                *reinterpret_cast<uint32_t*>(sp) = 0;
                if (threadIdx.x == 0) {
                    sp[-1] = 0;
                } else if (threadIdx.x == blockDim.x-1) {
                    sp[4] = 0;
                }
            }
            sp += (kMacroCellWidth+8);
        } else if (threadIdx.y == blockDim.y-1) {
            sp += (kMacroCellWidth+8);
            if (blockIdx.y != gridDim.y-1) {
                ip += ipitch;
                *reinterpret_cast<uint32_t*>(sp) = *reinterpret_cast<uint32_t*>(ip);
                if (threadIdx.x == 0) {
                    sp[-1] = blockIdx.x == 0 ? 0 : ip[-1];
                } else if (threadIdx.x == blockDim.x-1) {
                    sp[4] = blockIdx.x == gridDim.x-1 ? 0 : ip[4];
                }
            } else {
                *reinterpret_cast<uint32_t*>(sp) = 0;
                if (threadIdx.x == 0) {
                    sp[-1] = 0;
                } else if (threadIdx.x == blockDim.x-1) {
                    sp[4] = 0;
                }
            }
            sp -= (kMacroCellWidth+8);
        }
    } else {
        *reinterpret_cast<uint32_t*>(sp) = 0;
        
        // load left and right boundary cells
        if (threadIdx.x == 0) {
            sp[-1] = 0;
        } else if (threadIdx.x == blockDim.x-1) {
            sp[4] = 0;
        }
    }
    
    __syncthreads();
    
    if (row >= caParams.height || col >= caParams.width)
        return;
    
    uint8_t *op = out + row * opitch + col * 4;
    
    if (T == CARule::WireWorld) {
        // WireWorld
        // TODO:
        
    } else {
        // Life-like
        for (int i = 0; i < 4; ++i) {
            int stride = kMacroCellWidth+8;
            // compute number of live neighbours
            int sum =
            (sp[i-stride-1] & 1) + (sp[i-stride] & 1) + (sp[i-stride+1] & 1) +
            (sp[i       -1] & 1)                      + (sp[i       +1] & 1) +
            (sp[i+stride-1] & 1) + (sp[i+stride] & 1) + (sp[i+stride+1] & 1);
            
            op[i] = ((1 << sum) & (sp[i] ? caParams.S : caParams.B)) ? 1 : 0;
        }
    }
}

// more sophisticated worklist-based approach
// may or may not be faster depending on the state
__global__ void kernelCAUpdateWorkList() {

}

////////////////////////////////////////////////////////////////////////////////

CARule::CARule(const char *s) {
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
        type = LifeLike;
        B = 1 << 3;
        S = (1 << 2) | (1 << 3);
    }
}

const char* CARule::toString() const {
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
    this->rule = new CARule(rule);
    
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
    fprintf(stderr, "Running %i steps of %s\n", n, rule->toString());
    double ref = CycleTimer::currentSeconds();
    
    for (int i = 0; i < n; ++i) {
        // TODO : kernel launch
        dim3 updateBlockDim(kMacroCellWidth / 4, kMacroCellHeight, 1);
        dim3 updateGridDim((width / 4 + updateBlockDim.x - 1) / updateBlockDim.x,
                           (height + updateBlockDim.y - 1) / updateBlockDim.y);
        
        if (rule->type == CARule::WireWorld)
            kernelCAUpdateNaive<CARule::WireWorld><<<updateGridDim, updateBlockDim>>>(
                (generation & 1) ? cell1 : cell0, (generation & 1) ? pitch1 : pitch0,
                (generation & 1) ? cell0 : cell1, (generation & 1) ? pitch0 : pitch1);
        else
            kernelCAUpdateNaive<CARule::LifeLike><<<updateGridDim, updateBlockDim>>>(
                (generation & 1) ? cell1 : cell0, (generation & 1) ? pitch1 : pitch0,
                (generation & 1) ? cell0 : cell1, (generation & 1) ? pitch0 : pitch1);
        
        cudaDeviceSynchronize();
        ++generation;
    }
    
    double end = CycleTimer::currentSeconds();
    fprintf(stderr, "Elapsed : %lf ms (%lf / step)\n",
            (end - ref) * 1000.0, ((end - ref) * 1000.0) / (double)n);
}

bool CASim::setCells(unsigned int width, unsigned int height, uint8_t max,
                     const uint8_t *cells) {
    // check that the input respects the maximum number of states of the CA
    uint8_t maxState = rule->type == CARule::WireWorld ? 3 : 1;
    if (max > maxState) {
        fprintf(stderr, "Input max value outside of CA bounds (%u > %u).\n",
                (unsigned int)max, (unsigned int)maxState);
        return false;
    }
    
    if (cell0) cudaFree(cell0);
    if (cell1) cudaFree(cell1);
    
    this->width = width;
    this->height = height;
    
    generation = 0;
    
    // Alloc double buffers
    // TODO: pad buffers to simplify kernels?
    cudaMallocPitch(&cell0, &pitch0, width, height);
    cudaMallocPitch(&cell1, &pitch1, width, height);
    
    // set global CA params
    CAParams params;
    params.width = width;
    params.height = height;
    params.type = rule->type;
    params.B = rule->B;
    params.S = rule->S;
    params.cell0 = cell0;
    params.cell1 = cell1;
    cudaMemcpyToSymbol(caParams, &params, sizeof(CAParams));
    
    // copy initial state to first buffer
    // for now simple 1:1 mapping, each byte being one cell but may change in
    // the future to allow further optimization
    cudaMemcpy2D(cell0, pitch0, cells, width, width, height, cudaMemcpyHostToDevice);
    
    return true;
}

void CASim::getCells(uint8_t *cells) {
    // copy current state from appropriate buffer
    // for now simple 1:1 mapping, each byte being one cell but may change in
    // the future to allow further optimization
    cudaMemcpy2D(cells, width,
                 (generation & 1) ? cell1 : cell0,
                 (generation & 1) ? pitch1 : pitch0,
                 width, height, cudaMemcpyDeviceToHost);
}
