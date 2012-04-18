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
#include <thrust/unique.h>

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
    uint16_t B, S;
    
    unsigned int width;
    unsigned int height;
    unsigned int *workOffset;
};

__constant__ CAParams caParams;

////////////////////////////////////////////////////////////////////////////////

enum {
    kMacroCellWidth = 32,
    kMacroCellHeight = 16
};

// core CA update
// responsible for updating 4 adjacent cells in a rows
// * mid is the old values of the 4 cells to update (little endian)
// * sp is the shared memory array containing the old values of all cells in the
// macrocell being updated plus all boundary cells (Moore neighbourhood). Aligned
// on 4-byte boundary
// * return the new values of the cells (same encoding)
template <uint16_t T>
inline __device__ uint32_t CAUpdateCore(uint32_t mid, uint8_t *sp) {
    uint32_t nval = 0;
    if (T == CARule::WireWorld) {
        // WireWorld
        // TODO:
        
    } else {
        // Life-like
        const int S = caParams.S, B = caParams.B;
        
        const int stride = kMacroCellWidth+8;
        // NOTE: assume GPU in little-endian mode
        // TODO: get driver to ensure little-endian mode during setup
        const uint32_t top = *reinterpret_cast<uint32_t*>(sp-stride);
        const uint32_t bot = *reinterpret_cast<uint32_t*>(sp+stride);
        
        const int atop = (top & 0xff), btop = ((top >> 8) & 0xff), ctop = ((top >> 16) & 0xff), dtop  = ((top >> 24) & 0xff);
        const int amid = (mid & 0xff), bmid = ((mid >> 8) & 0xff), cmid = ((mid >> 16) & 0xff), dmid  = ((mid >> 24) & 0xff);
        const int abot = (bot & 0xff), bbot = ((bot >> 8) & 0xff), cbot = ((bot >> 16) & 0xff), dbot  = ((bot >> 24) & 0xff);
        
        int sum =
            sp[-stride-1] + atop + btop +
            sp[       -1]        + bmid +
            sp[ stride-1] + abot + bbot;
        
        nval |= ((1 << sum) & (amid ? S : B)) ? 1 : 0;
        
        sum =
            atop + btop + ctop +
            amid        + cmid +
            abot + bbot + cbot;
        
        nval |= ((1 << sum) & (bmid ? S : B)) ? (1 << 8) : 0;
        
        sum =
            btop + ctop + dtop +
            bmid        + dmid +
            bbot + cbot + dbot;
        
        nval |= ((1 << sum) & (cmid ? S : B)) ? (1 << 16) : 0;
        
        sum =
            ctop + dtop + sp[-stride+4] +
            cmid        + sp[        4] +
            cbot + dbot + sp[ stride+4];
        
        nval |= ((1 << sum) & (dmid ? S : B)) ? (1 << 24) : 0;
    }
    return nval;
}

//
//
inline __device__ uint32_t loadCells(int row, int col, uint8_t *sp, uint8_t *ip, size_t pitch) {
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
            ip -= pitch;
            *reinterpret_cast<uint32_t*>(sp) = *reinterpret_cast<uint32_t*>(ip);
            if (threadIdx.x == 0) {
                sp[-1] = blockIdx.x == 0 ? 0 : ip[-1];
            } else if (threadIdx.x == blockDim.x-1) {
                sp[4] = blockIdx.x == gridDim.x-1 ? 0 : ip[4];
            }
            sp += (kMacroCellWidth+8);
        } else if (threadIdx.y == blockDim.y-1) {
            sp += (kMacroCellWidth+8);
            ip += pitch;
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
    return mid;
}

// naive, embarassingly parallel CA update
template <uint16_t T>
__global__ void kernelCAUpdateNaive(uint8_t *in,  uint8_t *out, size_t pitch) {
    __shared__ uint8_t ocells[(kMacroCellHeight+2)*(kMacroCellWidth+8)];
    
    unsigned int row = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int col = blockIdx.x * blockDim.x + threadIdx.x;
    
    uint8_t *ip = in + row * pitch + col * 4;
    
    // cooperative load of relevant old cell values to shared memory
    uint8_t *sp = ocells + (threadIdx.y+1) * (kMacroCellWidth+8) + (threadIdx.x+1) * 4;
    
    uint32_t mid = loadCells(row, col, sp, ip, pitch);
    
    // wait for shared array to be fully initialized
    __syncthreads();
    
    // only update relevant cells
    if (row >= caParams.height || col >= caParams.width)
        return;
    
    uint32_t nval = CAUpdateCore<T>(mid, sp);
    uint8_t *op = out + row * pitch + col * 4;
    *reinterpret_cast<uint32_t*>(op) = nval;
}

// fill initial worklist (all cells to be updated)
__global__ void kernelInitWorklist(uint32_t *work, uint32_t nCellX, uint32_t nCellY) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nCellX*nCellY)
        return;
    work[2*idx+0] = (idx / nCellX) * kMacroCellHeight;
    work[2*idx+1] = (idx % nCellX) * (kMacroCellWidth / 4);
}

// more sophisticated worklist-based approach
// may or may not be faster depending on the state
template <uint16_t T>
__global__ void kernelCAUpdateWorkList(uint32_t *iwork, uint32_t *owork,
                                       uint8_t *in, uint8_t *out, size_t pitch) {
    __shared__ uint32_t base[2];
    __shared__ uint32_t update[9];
    __shared__ uint8_t ocells[(kMacroCellHeight+2)*(kMacroCellWidth+8)];
    
    unsigned int idx = blockIdx.x;
    unsigned int tidx = threadIdx.y * blockDim.x + threadIdx.x;
    
    // cooperatively load macrocell position from worklist
    if (tidx < 2)
        base[tidx] = iwork[2*idx+tidx];
    __syncthreads();
    
    // derive row/col from macrocell position
    unsigned int row = base[0] + threadIdx.y;
    unsigned int col = base[1] + threadIdx.x;
    
    uint8_t *ip = in + row * pitch + col * 4;
    
    // cooperative load of relevant old cell values to shared memory
    uint8_t *sp = ocells + (threadIdx.y+1) * (kMacroCellWidth+8) + (threadIdx.x+1) * 4;
    
    uint32_t mid = loadCells(row, col, sp, ip, pitch);
    
    // wait for shared array to be fully initialized
    __syncthreads();
    
    if (row >= caParams.height || col >= caParams.width)
        return;
    
    // update cells
    uint32_t nval = CAUpdateCore<T>(mid, sp);
    uint8_t *op = out + row * pitch + col * 4;
    *reinterpret_cast<uint32_t*>(op) = nval;
    
    // determine modification
    uint32_t mod = mid ^ nval;
    
    // reset mod flags
    if (tidx < 9)
        update[tidx] = 0;
    
    __syncthreads();
    
    // compute mod flags
    // This code relies on the fact that if multiple threads try to write to the
    // same shared memory address one will always succeed. Which does write is
    // undefined but since we're always writing the same value it doesn't matter.
    if (mod) {
        update[4] = 1;
        
        if (threadIdx.x == 0 && (mod & 0xff)) {
            update[3] = 1;
            if (threadIdx.y == 0)
                update[0] = 1;
            else if (threadIdx.y == blockDim.y-1)
                update[6] = 1;
        } else if (threadIdx.x == blockDim.x-1 && (mod & 0xff000000)) {
            update[5] = 1;
            if (threadIdx.y == 0)
                update[2] = 1;
            else if (threadIdx.y == blockDim.y-1)
                update[8] = 1;
        }
        
        if (threadIdx.y == 0)
            update[1] = 1;
        else if (threadIdx.y == blockDim.y-1)
            update[7] = 1;
        
    }
    
    __syncthreads();
    
    unsigned int cmod = update[4];
    
    // fill worklist based on mod flags
    if (tidx == 0 && cmod) {
        unsigned int count = update[0] + update[1] + update[2] +
                             update[3] + cmod      + update[5] + 
                             update[6] + update[7] + update[8];
        
        unsigned int offset = atomicAdd(caParams.workOffset, count);
        uint32_t *ow = owork + 2*offset;
        bool lastRow = row == (caParams.height-kMacroCellHeight);
        bool lastCol = col == (caParams.width-kMacroCellWidth/4);
        if (row) {
            if (update[0] && col)       { ow[0] = row-kMacroCellHeight; ow[1] = col-kMacroCellWidth/4; ow += 2; }
            if (update[1])              { ow[0] = row-kMacroCellHeight; ow[1] = col;                   ow += 2; }
            if (update[2] && !lastCol)  { ow[0] = row-kMacroCellHeight; ow[1] = col+kMacroCellWidth/4; ow += 2; }
        }
        if (update[3] && col)           { ow[0] = row;                  ow[1] = col-kMacroCellWidth/4; ow += 2; }
        if (cmod)                       { ow[0] = row;                  ow[1] = col;                   ow += 2; }
        if (update[5] && !lastCol)      { ow[0] = row;                  ow[1] = col+kMacroCellWidth/4; ow += 2; }
        if (!lastRow) {
            if (update[6] && col)       { ow[0] = row+kMacroCellHeight; ow[1] = col-kMacroCellWidth/4; ow += 2; }
            if (update[7])              { ow[0] = row+kMacroCellHeight; ow[1] = col;                   ow += 2; }
            if (update[8] && !lastCol)  { ow[0] = row+kMacroCellHeight; ow[1] = col+kMacroCellWidth/4; }
        }
    }
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
    work0 = work1 = 0;
    workOffset = 0;
    
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
    if (work0) cudaFree(work0);
    if (work1) cudaFree(work1);
    if (workOffset) cudaFree(workOffset);
}

static int cmp(const void *a, const void *b) {
    int diff = *((const int*)a) - *((const int*)b);
    return diff ? diff : ((const int*)a)[1] - ((const int*)b)[1];
}

void CASim::step(int n) {
    fprintf(stderr, "Running %i steps of %s\n", n, rule->toString());
    double ref = CycleTimer::currentSeconds();
    
    uint64_t *mwork = 0;
    uint32_t wAmount, wAlloc = 0;
    cudaMemcpy(&wAmount, workOffset, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    thrust::device_ptr<uint64_t> wit1(reinterpret_cast<uint64_t*>(work1));
    thrust::device_ptr<uint64_t> wit0(reinterpret_cast<uint64_t*>(work0));
    
    for (int i = 0; i < n; ++i) {
        uint8_t *scell = (generation & 1) ? cell1 : cell0;
        uint8_t *dcell = (generation & 1) ? cell0 : cell1;
        dim3 updateBlockDim(kMacroCellWidth / 4, kMacroCellHeight, 1);
        
#if 0
        dim3 updateGridDim((width / 4 + updateBlockDim.x - 1) / updateBlockDim.x,
                           (height + updateBlockDim.y - 1) / updateBlockDim.y);
        
        if (rule->type == CARule::WireWorld)
            kernelCAUpdateNaive<CARule::WireWorld><<<updateGridDim, updateBlockDim>>>(
                scell + pitch, dcell + pitch, pitch);
        else
            kernelCAUpdateNaive<CARule::LifeLike><<<updateGridDim, updateBlockDim>>>(
                scell + pitch, dcell + pitch, pitch);
        
        cudaMemcpy(&wAmount, workOffset, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
#else
        // when worklist contains all the macrocells, the update phase alone
        // is more than twice as slow as the naive approach as it needs to read
        // from the worklist
        // TODO: switch back and forth between naive and worklist depending on
        // the size of the worklist?
        
        //fprintf(stderr, "updating %u macrocells\n", wAmount);
        
        cudaMemset(workOffset, 0, sizeof(uint32_t));
        cudaDeviceSynchronize();
        
        int wLim = 1024;
        
        for (int i = 0; i < (wAmount + wLim-1)/wLim; ++i) {
            dim3 updateGridDim(min(wLim, wAmount - i*wLim), 1);
            if (rule->type == CARule::WireWorld)
                kernelCAUpdateWorkList<CARule::WireWorld><<<updateGridDim, updateBlockDim>>>(
                    work0 + i*2*wLim, work1, scell + pitch, dcell + pitch, pitch);
            else
                kernelCAUpdateWorkList<CARule::LifeLike><<<updateGridDim, updateBlockDim>>>(
                    work0 + i*2*wLim, work1, scell + pitch, dcell + pitch, pitch);
            
        }
        cudaDeviceSynchronize();
        
        uint32_t wNext;
        cudaMemcpy(&wNext, workOffset, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        if (!wNext) {
            fprintf(stderr, "Reached still life\n");
            break;
        }
        
        // TODO: find best empirical threshold
        const uint32_t threshold = 1024;
        
        if (wNext < threshold) {
            // too small : no point doing the sort and duplicate removal on GPU
            
            // make sure the CPU-side worklist buffer is large enough
            if (wNext > wAlloc) {
                wAlloc = wNext;
                mwork = (uint64_t*)realloc(mwork, wAlloc * sizeof(uint64_t));
            }
            
            // copy worklist from GPU
            cudaMemcpy(mwork, work1, wNext * sizeof(uint64_t), cudaMemcpyDeviceToHost);
            
            // sort worklist
            qsort(mwork, wNext, sizeof(uint64_t), cmp);
            
            // remove duplicates
            uint64_t last = *mwork;
            uint64_t *h = mwork, *q = mwork+1, *e = mwork + wNext;
            while (++h != e) {
                uint64_t cur = *h;
                if (cur != last) {
                    last = cur;
                    if (q != h)
                        *q = cur;
                    ++q;
                }
            }
            wNext = q - mwork;
            
            // copy filtered worklist to GPU for next step
            cudaMemcpy(work0, mwork, wNext * sizeof(uint64_t), cudaMemcpyHostToDevice);
            wAmount = wNext;
        } else {
            // sort worklist
            thrust::sort(wit1, wit1 + wNext);
            // remove duplicates
            wAmount = thrust::unique_copy(wit1, wit1 + wNext, wit0) - wit0;
        }
#endif
        ++generation;
        cudaDeviceSynchronize();
    }
    
    free(mwork);
    
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
    if (work0) cudaFree(work0);
    if (work1) cudaFree(work1);
    
    double tmp, ref = CycleTimer::currentSeconds();
    
    this->width = width;
    this->height = height;
    
    generation = 0;
    
    cudaError_t err;
    
    // Alloc double buffers
    // buffers are padded as follows to simplify kernels :
    // * add a top row (always set to 0)
    // * add a bottom row (always set to 0)
    // * ensure the width of the allocated array is a multiple of 16
    int wpad = (width & (kMacroCellWidth - 1));
    if (wpad)
        wpad = kMacroCellWidth - wpad;
    size_t pitch0, pitch1;
    err = cudaMallocPitch(&cell0, &pitch0, width + wpad, height+2);
    if (err) { fprintf(stderr, "Unable to allocate buffer\n"); return false; }
    err = cudaMallocPitch(&cell1, &pitch1, width + wpad, height+2);
    if (err) { fprintf(stderr, "Unable to allocate buffer\n"); return false; }
    
    // alloc worklists (macrocell granularity)
    int nCellX = (width + (kMacroCellWidth - 1)) / kMacroCellWidth;
    int nCellY = (height + (kMacroCellHeight - 1)) / kMacroCellHeight;
    uint32_t nCell = nCellX * nCellY;
    err = cudaMalloc(&work0, nCell * 2 * sizeof(uint32_t));
    if (err) { fprintf(stderr, "Unable to allocate worklist\n"); return false; }
    err = cudaMalloc(&work1, nCell * 9 * 2 * sizeof(uint32_t));
    if (err) { fprintf(stderr, "Unable to allocate worklist\n"); return false; }
    if (!workOffset)
        cudaMalloc(&workOffset, sizeof(uint32_t));
    
    tmp = CycleTimer::currentSeconds();
    fprintf(stderr, "Alloc buffers : %lf ms\n", (tmp - ref)*1000.0);
    
    if (pitch0 != pitch1) {
        fprintf(stderr, "Pitch mismatch\n");
        return false;
    }
    pitch = pitch0;
    
    // set global CA params
    CAParams params;
    params.type = rule->type;
    params.B = rule->B;
    params.S = rule->S;
    params.width = width;
    params.height = height;
    params.workOffset = workOffset;
    cudaMemcpyToSymbol(caParams, &params, sizeof(CAParams));
    
    // copy initial state to first buffer
    // for now simple 1:1 mapping, each byte being one cell but may change in
    // the future to allow further optimization
    cudaMemcpy2D(cell0+pitch, pitch, cells, width, width, height, cudaMemcpyHostToDevice);
    
    // reset padding rows
    cudaMemset(cell0, 0, width);
    cudaMemset(cell0+(height+1)*pitch, 0, width);
    cudaMemset(cell1, 0, width);
    cudaMemset(cell1+(height+1)*pitch, 0, width);
    
    // reset padding columns
    if (wpad) {
//         for (int i = 0; i < height+2; ++i) {
//             cudaMemset(cell0+i*pitch+width, 0, wpad);
//             cudaMemset(cell1+i*pitch+width, 0, wpad);
//         }
        cudaMemset2D(cell0+width, pitch, 0, wpad, height+2);
    }
    
    
    // prepare initial worklist
    dim3 blockDim(128, 1);
    dim3 gridDim((nCell + 127) / 128, 1);
    kernelInitWorklist<<<gridDim, blockDim>>>(work0, nCellX, nCellY);
    cudaMemcpy(workOffset, &nCell, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    
    double end = CycleTimer::currentSeconds();
    fprintf(stderr, "Initialized CA in %lf ms\n", (end - ref)*1000.0);
    
    return true;
}

void CASim::getCells(uint8_t *cells) {
    // copy current state from appropriate buffer
    // for now simple 1:1 mapping, each byte being one cell but may change in
    // the future to allow further optimization
    uint8_t *src = (generation & 1) ? cell1 : cell0;
    
    cudaMemcpy2D(cells, width, src+pitch, pitch,
                 width, height, cudaMemcpyDeviceToHost);
}
