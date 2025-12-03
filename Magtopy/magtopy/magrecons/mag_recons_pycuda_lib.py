#!/usr/bin/env python3
# Author: Marisel Di Pietro Martinez
# Any work making use of these codes should refer to, and cite, the
# following paper:
# Di Pietro Martínez et al., Phys. Rev. B 107, 094425 (2023)
# [https://doi.org/10.1103/PhysRevB.107.094425]

import numpy as np
import pycuda.autoinit
import pycuda.gpuarray as gpuarray
from pycuda.compiler import SourceModule

#block_size = 10
block_3d = (16,8,8) # = 1024 total
block_2d = (32,32,1) # = 1024 total

def grid_3d(L):
    grid_x = int((L + block_3d[0] - 1)/block_3d[0])
    grid_y = int((L + block_3d[1] - 1)/block_3d[1])
    grid_z = int((L + block_3d[2] - 1)/block_3d[2])
    return (grid_x, grid_y, grid_z)

def grid_2d(L):
    grid_x = int((L + block_2d[0] - 1)/block_2d[0])
    grid_y = int((L + block_2d[1] - 1)/block_2d[1])
    return (grid_x, grid_y, 1)


mod = SourceModule("""
__global__ void rotate_mag(float* MZ_rot, const float* MX, const float* MY, const float* MZ, const float* RHO, const unsigned int imageWidth, const unsigned int imageHeight, const unsigned int imageDepth, const float R3_1, const float R3_2, const float R3_3)
{
    // compute thread dimension
    const unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
    const unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

    //// compute target address
    const unsigned int idx = z + y * imageDepth + x * imageHeight * imageDepth;
    
    if (idx < imageWidth*imageWidth*imageWidth){
        MZ_rot[idx] = MX[idx]*R3_1 + MY[idx]*R3_2 + MZ[idx]*R3_3 + RHO[idx];
    }

}

__global__ void rotate_z(float* trg, const float* src, const unsigned int L, const float angle)
{    
    // compute thread dimension
    const unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
    const unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

    //// compute target address
    const unsigned int idx = z + y * L + x * L * L;
    
    if (idx < L*L*L){
        const float center =  L/2;
        
        const float x_original = (x - center)  * cos(angle) + (y - center) * sin(angle) + center;
        const float y_original = - (x - center)  * sin(angle) + (y - center) * cos(angle) + center;

        const int p = (int) (floor(x_original));
        const int q = (int) (floor(y_original));

        const float a = x_original - (float) p;
        const float b = y_original - (float) q;
        
        if ((p >= 0) && (p < L) && (q >= 0) && (q < L)){
            const unsigned int new_idx = z + q * L + p * L * L;
            trg[idx] = trg[idx] + (1-a)*(1-b)*src[new_idx];
            }
        if ((p >= 0) && (p < L) && ((q+1) >= 0) && ((q+1) < L)){
            const unsigned int new_idx = z + (q+1) * L + p * L * L;
            trg[idx] = trg[idx] + (1-a)*b*src[new_idx];
            }
        if (((p+1) >= 0) && ((p+1) < L) && (q >= 0) && (q < L)){
            const unsigned int new_idx = z + q * L + (p+1) * L * L;
            trg[idx] = trg[idx] + a*(1-b)*src[new_idx];
            }

        if (((p+1) >= 0) && ((p+1) < L) && ((q+1) >= 0) && ((q+1) < L)){
            const unsigned int new_idx = z + (q+1) * L + (p+1) * L * L;
            trg[idx] = trg[idx] + a*b*src[new_idx];
            }
    }
}

__global__ void rotate_y(float* trg, const float* src, const unsigned int L, const float angle)
{    
    // compute thread dimension
    const unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
    const unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

    //// compute target address
    const unsigned int idx = z + y * L + x * L * L;
    
    if (idx < L*L*L){    
        const float center =  L/2;
        
        const float x_original = (x - center)  * cos(angle) + (z - center) * sin(angle) + center;
        const float z_original = - (x - center)  * sin(angle) + (z - center) * cos(angle) + center;

        const int p = (int) (floor(x_original));
        const int q = (int) (floor(z_original));

        const float a = x_original - (float) p;
        const float b = z_original - (float) q;
        
        if ((p >= 0) && (p < L) && (q >= 0) && (q < L)){
            const unsigned int new_idx = q + y * L + p * L * L;
            trg[idx] = trg[idx] + (1-a)*(1-b)*src[new_idx];
            }
        if ((p >= 0) && (p < L) && ((q+1) >= 0) && ((q+1) < L)){
            const unsigned int new_idx = q + 1 + y * L + p * L * L;
            trg[idx] = trg[idx] + (1-a)*b*src[new_idx];
            }
        if (((p+1) >= 0) && ((p+1) < L) && (q >= 0) && (q < L)){
            const unsigned int new_idx = q + y * L + (p+1) * L * L;
            trg[idx] = trg[idx] + a*(1-b)*src[new_idx];
            }

        if (((p+1) >= 0) && ((p+1) < L) && ((q+1) >= 0) && ((q+1) < L)){
            const unsigned int new_idx = q + 1 + y * L + (p+1) * L * L;
            trg[idx] = trg[idx] + a*b*src[new_idx];
            }
    }
}

__global__ void sum_z(float* proj, const float* MZ_rot, const unsigned int imageWidth, const unsigned int imageHeight, const unsigned int imageDepth)
{
    // compute thread dimension
    const unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
    
    //// compute target address
    const unsigned int idx = y + x * imageHeight;
    
    if (idx < imageWidth*imageWidth){
        proj[idx] = 0;
        for(int z = 0; z < imageDepth; z++)
            proj[idx] += MZ_rot[idx * imageDepth + z];
    }
}

__global__ void back_sum_z(float* MZ, const float* proj, const unsigned int imageWidth, const unsigned int imageHeight, const unsigned int imageDepth)
{
    // compute thread dimension
    const unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
    const unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

    //// compute target address
    const unsigned int idx = y + x * imageHeight;
    
    const unsigned int idx3 = idx * imageDepth + z;
    if (idx3 < imageWidth*imageWidth*imageWidth){
        MZ[idx3] = proj[idx];
    }
}
""")

rotate_z = mod.get_function("rotate_z")
rotate_y = mod.get_function("rotate_y")
rotate_mag = mod.get_function("rotate_mag")
sum_z = mod.get_function("sum_z")
back_sum_z = mod.get_function("back_sum_z")

def rot_z(MZ, phi):
    #prepare data to pass to CUDA function
    L = np.int32(MZ.shape[0])
    phi = np.float32(phi)
    dest_gpu = gpuarray.zeros_like(MZ)
    
    #rotate along z axis (counter-clockwise)
    rotate_z(dest_gpu.gpudata, MZ.gpudata, L, phi, block=block_3d, grid=grid_3d(L))
    
    return dest_gpu

def rot_y(MZ, theta):
    #prepare data to pass to CUDA function
    L = np.int32(MZ.shape[0])
    theta = np.float32(theta)
    dest_gpu = gpuarray.zeros_like(MZ)

    #rotate along y axis (counter-clockwise)
    rotate_y(dest_gpu.gpudata, MZ.gpudata, L, theta, block=block_3d, grid=grid_3d(L))
    
    return dest_gpu

def mag_rotate(MX, MY, MZ, RHO, R3_1, R3_2, R3_3):
    L = np.int32(MZ.shape[0])
    MZ_rot = gpuarray.zeros_like(MZ)
    
    R3_1 = np.float32(R3_1)
    R3_2 = np.float32(R3_2)
    R3_3 = np.float32(R3_3)
    
    rotate_mag(MZ_rot.gpudata, MX.gpudata, MY.gpudata, MZ.gpudata, RHO.gpudata, L, L, L, R3_1, R3_2, R3_3, block=block_3d, grid=grid_3d(L))
    
    return MZ_rot

def project(MZ_rot):
    L = np.int32(MZ_rot.shape[0])
    proj = gpuarray.zeros((L,L), dtype=np.float32)
    sum_z(proj.gpudata, MZ_rot.gpudata, L, L, L, block=block_2d, grid=grid_2d(L))
    
    return proj

def back_project(proj):
    L = np.int32(proj.shape[0])
    MZ = gpuarray.zeros((L,L,L), dtype=np.float32)
    back_sum_z(MZ.gpudata, proj.gpudata, L, L, L, block=block_3d, grid=grid_3d(L))
    
    return MZ
