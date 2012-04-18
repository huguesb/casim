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
    
    uint32_t mid;
    
    // load inner cells cooperatively
    if (row < caParams.height && col < caParams.width) {
        mid = *reinterpret_cast<uint32_t*>(ip);
        *reinterpret_cast<uint32_t*>(sp) = mid;
        
        // load left and right boundary cells
        if (threadIdx.x == 0) {
            sp[-1] = blockIdx.x == 0 ? 0 : ip[-1];
        } else if (threadIdx.x == blockDim.x-1) {
            sp[4] = blockIdx.x == gridDim.x-1 ? 0 : ip[4];
        }
        
        // load top and bottom boundary cells
        if (threadIdx.y == 0) {
            sp -= (kMacroCellWidth+8);
            ip -= ipitch;
            *reinterpret_cast<uint32_t*>(sp) = *reinterpret_cast<uint32_t*>(ip);
            if (threadIdx.x == 0) {
                sp[-1] = blockIdx.x == 0 ? 0 : ip[-1];
            } else if (threadIdx.x == blockDim.x-1) {
                sp[4] = blockIdx.x == gridDim.x-1 ? 0 : ip[4];
            }
            sp += (kMacroCellWidth+8);
        } else if (threadIdx.y == blockDim.y-1) {
            sp += (kMacroCellWidth+8);
            ip += ipitch;
            *reinterpret_cast<uint32_t*>(sp) = *reinterpret_cast<uint32_t*>(ip);
            if (threadIdx.x == 0) {
                sp[-1] = blockIdx.x == 0 ? 0 : ip[-1];
            } else if (threadIdx.x == blockDim.x-1) {
                sp[4] = blockIdx.x == gridDim.x-1 ? 0 : ip[4];
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
        uint32_t nval = 0;
        int stride = kMacroCellWidth+8;
#if 0
        for (int i = 0; i < 4; ++i) {
            // compute number of live neighbours
            int sum =
            (sp[i-stride-1] & 1) + (sp[i-stride] & 1) + (sp[i-stride+1] & 1) +
            (sp[i       -1] & 1)                      + (sp[i       +1] & 1) +
            (sp[i+stride-1] & 1) + (sp[i+stride] & 1) + (sp[i+stride+1] & 1);
            
            int ne = ((1 << sum) & (sp[i] ? caParams.S : caParams.B)) ? 1 : 0;
            nval |= ne << (8*i);
        }
#else
        // TODO: check that little-endian is certain on all GPUs
        uint32_t top = *reinterpret_cast<uint32_t*>(sp-stride);
        //uint32_t mid = *reinterpret_cast<uint32_t*>(sp);
        uint32_t bot = *reinterpret_cast<uint32_t*>(sp+stride);
        
        int atop = (top & 0xff), btop = ((top >> 8) & 0xff), ctop = ((top >> 16) & 0xff), dtop  = ((top >> 24) & 0xff);
        int amid = (mid & 0xff), bmid = ((mid >> 8) & 0xff), cmid = ((mid >> 16) & 0xff), dmid  = ((mid >> 24) & 0xff);
        int abot = (bot & 0xff), bbot = ((bot >> 8) & 0xff), cbot = ((bot >> 16) & 0xff), dbot  = ((bot >> 24) & 0xff);
        
        int sum =
            sp[-stride-1] + atop + btop +
            sp[       -1]        + bmid +
            sp[ stride-1] + abot + bbot;
        
        //int na = ((1 << sum) & (amid ? caParams.S : caParams.B)) ? 1 : 0;
        nval |= ((1 << sum) & (amid ? caParams.S : caParams.B)) ? 1 : 0;
        
        sum =
            atop + btop + ctop +
            amid        + cmid +
            abot + bbot + cbot;
        
//         int nb = ((1 << sum) & (bmid ? caParams.S : caParams.B)) ? 1 : 0;
        nval |= ((1 << sum) & (bmid ? caParams.S : caParams.B)) ? (1 << 8) : 0;
        
        sum =
            btop + ctop + dtop +
            bmid        + dmid +
            bbot + cbot + dbot;
        
//         int nc = ((1 << sum) & (cmid ? caParams.S : caParams.B)) ? 1 : 0;
        nval |= ((1 << sum) & (cmid ? caParams.S : caParams.B)) ? (1 << 16) : 0;
        
        sum =
            ctop + dtop + sp[-stride+4] +
            cmid        + sp[        4] +
            cbot + dbot + sp[ stride+4];
        
//         int nd = ((1 << sum) & (dmid ? caParams.S : caParams.B)) ? 1 : 0;
        nval |= ((1 << sum) & (dmid ? caParams.S : caParams.B)) ? (1 << 24) : 0;
        
//         nval = (nd << 24) | (nc << 16) | (nb << 8) | na;
#endif
        *reinterpret_cast<uint32_t*>(op) = nval;
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
        dim3 updateBlockDim(kMacroCellWidth / 4, kMacroCellHeight, 1);
        dim3 updateGridDim((width / 4 + updateBlockDim.x - 1) / updateBlockDim.x,
                           (height + updateBlockDim.y - 1) / updateBlockDim.y);
        
        uint8_t *src = (generation & 1) ? cell1 : cell0;
        uint8_t *dst = (generation & 1) ? cell0 : cell1;
        size_t spitch = (generation & 1) ? pitch1 : pitch0;
        size_t dpitch = (generation & 1) ? pitch0 : pitch1;
        
        if (rule->type == CARule::WireWorld)
            kernelCAUpdateNaive<CARule::WireWorld><<<updateGridDim, updateBlockDim>>>(
                src + spitch, spitch, dst + dpitch, dpitch);
        else
            kernelCAUpdateNaive<CARule::LifeLike><<<updateGridDim, updateBlockDim>>>(
                src + spitch, spitch, dst + dpitch, dpitch);
        
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
    // buffers are padded as follows to simplify kernels :
    // * add a top row (always set to 0)
    // * add a bottom row (always set to 0)
    // * ensure the width of the allocated array is a multiple of 16
    int wpad = (width & 15);
    if (wpad)
        wpad = 16 - wpad;
    cudaMallocPitch(&cell0, &pitch0, width + wpad, height+2);
    cudaMallocPitch(&cell1, &pitch1, width + wpad, height+2);
    
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
    cudaMemcpy2D(cell0+pitch0, pitch0, cells, width, width, height, cudaMemcpyHostToDevice);
    
    // reset padding rows
    cudaMemset(cell0, 0, width);
    cudaMemset(cell0+(height+1)*pitch0, 0, width);
    cudaMemset(cell1, 0, width);
    cudaMemset(cell1+(height+1)*pitch1, 0, width);
    
    // reset padding columns
    if (wpad) {
//         for (int i = 0; i < height+2; ++i) {
//             cudaMemset(cell0+i*pitch0+width, 0, wpad);
//             cudaMemset(cell1+i*pitch1+width, 0, wpad);
//         }
        cudaMemset2D(cell0+width, pitch0, 0, wpad, height+2);
    }
    
    return true;
}

void CASim::getCells(uint8_t *cells) {
    // copy current state from appropriate buffer
    // for now simple 1:1 mapping, each byte being one cell but may change in
    // the future to allow further optimization
    uint8_t *src = (generation & 1) ? cell1 : cell0;
    size_t spitch = (generation & 1) ? pitch1 : pitch0;
    
    cudaMemcpy2D(cells, width, src+spitch, spitch,
                 width, height, cudaMemcpyDeviceToHost);
}
