#pragma once
#include <assert.h>
#include "GP_Attribute.h"
#include "GP_Context.h"
#include "Kernels/GP_Query_Dispatcher.cuh"
#include "Util/GP_Vector.h"
#include <Eigen/Core>

#define OFFSET(row, col, ld) ((row) * (ld) + (col))
#define FLOAT4(pointer) (reinterpret_cast<float4*>(&(pointer))[0])

//__global__ void naiveSgemm(
//    float * __restrict__ a, float * __restrict__ b, float * __restrict__ c,
//    const int M, const int N, const int K) {
//
//    int n = blockIdx.x * blockDim.x + threadIdx.x;
//    int m = blockIdx.y * blockDim.y + threadIdx.y;
//    if (m < M && n < N) {
//        float psum = 0.0;
//#pragma unroll
//        for (int k = 0; k < K; k++) {
//            psum += a[OFFSET(m, k, K)] * b[OFFSET(k, n, N)];
//        }
//        c[OFFSET(m, n, N)] = psum;
//    }
//}
//
//__global__ void mySgemmV1Aligned(
//    float * __restrict__ a, float * __restrict__ b, float * __restrict__ c,
//    const int M, const int N, const int K) {
//
//    const int BM = 128;
//    const int BN = 128;
//    const int BK = 8;
//    const int TM = 8;
//    const int TN = 8;
//
//    const int bx = blockIdx.x;
//    const int by = blockIdx.y;
//    const int tx = threadIdx.x;
//    const int ty = threadIdx.y;
//    const int tid = ty * blockDim.x + tx;
//
//    __shared__ float s_a[BM][BK];
//    __shared__ float s_b[BK][BN];
//
//    float r_c[TM][TN] = {0.0};
//
//    int load_a_smem_m = tid >> 1;
//    int load_a_smem_k = (tid & 1) << 2;
//    int load_b_smem_k = tid >> 5;
//    int load_b_smem_n = (tid & 31) << 2;
//
//    int load_a_gmem_m = by * BM + load_a_smem_m;
//    int load_b_gmem_n = bx * BN + load_b_smem_n;
//
//    for (int bk = 0; bk < (K + BK - 1) / BK; bk++) {
//        int load_a_gmem_k = bk * BK + load_a_smem_k;
//        int load_a_gmem_addr = OFFSET(load_a_gmem_m, load_a_gmem_k, K);
//        FLOAT4(s_a[load_a_smem_m][load_a_smem_k]) = FLOAT4(a[load_a_gmem_addr]);
//        int load_b_gmem_k = bk * BK + load_b_smem_k;
//        int load_b_gmem_addr = OFFSET(load_b_gmem_k, load_b_gmem_n, N);
//        FLOAT4(s_b[load_b_smem_k][load_b_smem_n]) = FLOAT4(b[load_b_gmem_addr]);
//
//        __syncthreads();
//
//#pragma unroll
//        for (int k = 0; k < BK; k++) {
//#pragma unroll
//            for (int m = 0; m < TM; m++) {
//#pragma unroll
//                for (int n = 0; n < TN; n++) {
//                    int comp_a_smem_m = ty * TM + m;
//                    int comp_b_smem_n = tx * TN + n;
//                    r_c[m][n] += s_a[comp_a_smem_m][k] * s_b[k][comp_b_smem_n];
//                }
//            }
//        }
//
//        __syncthreads();
//    }
//
//#pragma unroll
//    for (int i = 0; i < TM; i++) {
//        int store_c_gmem_m = by * BM + ty * TM + i;
//#pragma unroll
//        for (int j = 0; j < TN; j += 4) {
//            int store_c_gmem_n = bx * BN + tx * TN + j;
//            int store_c_gmem_addr = OFFSET(store_c_gmem_m, store_c_gmem_n, N);
//            FLOAT4(c[store_c_gmem_addr]) = FLOAT4(r_c[i][j]);
//        }
//    }
//}
//
//__global__ void mySgemmV2Aligned(
//    float * __restrict__ a, float * __restrict__ b, float * __restrict__ c,
//    const int M, const int N, const int K) {
//
//    const int BM = 128;
//    const int BN = 128;
//    const int BK = 8;
//    const int TM = 8;
//    const int TN = 8;
//
//    const int bx = blockIdx.x;
//    const int by = blockIdx.y;
//    const int tx = threadIdx.x;
//    const int ty = threadIdx.y;
//    const int tid = ty * blockDim.x + tx;
//
//    __shared__ float s_a[BK][BM];
//    __shared__ float s_b[BK][BN];
//
//    float r_load_a[4];
//    float r_load_b[4];
//    float r_comp_a[TM];
//    float r_comp_b[TN];
//    float r_c[TM][TN] = {0.0};
//
//    int load_a_smem_m = tid >> 1;
//    int load_a_smem_k = (tid & 1) << 2;
//    int load_b_smem_k = tid >> 5;
//    int load_b_smem_n = (tid & 31) << 2;
//
//    int load_a_gmem_m = by * BM + load_a_smem_m;
//    int load_b_gmem_n = bx * BN + load_b_smem_n;
//
//    for (int bk = 0; bk < (K + BK - 1) / BK; bk++) {
//
//        int load_a_gmem_k = bk * BK + load_a_smem_k;
//        int load_a_gmem_addr = OFFSET(load_a_gmem_m, load_a_gmem_k, K);
//        int load_b_gmem_k = bk * BK + load_b_smem_k;
//        int load_b_gmem_addr = OFFSET(load_b_gmem_k, load_b_gmem_n, N);
//        FLOAT4(r_load_a[0]) = FLOAT4(a[load_a_gmem_addr]);
//        FLOAT4(r_load_b[0]) = FLOAT4(b[load_b_gmem_addr]);
//
//        s_a[load_a_smem_k    ][load_a_smem_m] = r_load_a[0];
//        s_a[load_a_smem_k + 1][load_a_smem_m] = r_load_a[1];
//        s_a[load_a_smem_k + 2][load_a_smem_m] = r_load_a[2];
//        s_a[load_a_smem_k + 3][load_a_smem_m] = r_load_a[3];
//        FLOAT4(s_b[load_b_smem_k][load_b_smem_n]) = FLOAT4(r_load_b[0]);
//
//        __syncthreads();
//
//#pragma unroll
//        for (int tk = 0; tk < BK; tk++) {
//            FLOAT4(r_comp_a[0]) = FLOAT4(s_a[tk][ty * TM / 2         ]);
//            FLOAT4(r_comp_a[4]) = FLOAT4(s_a[tk][ty * TM / 2 + BM / 2]);
//            FLOAT4(r_comp_b[0]) = FLOAT4(s_b[tk][tx * TN / 2         ]);
//            FLOAT4(r_comp_b[4]) = FLOAT4(s_b[tk][tx * TN / 2 + BN / 2]);
//
//#pragma unroll
//            for (int tm = 0; tm < TM; tm++) {
//#pragma unroll
//                for (int tn = 0; tn < TN; tn++) {
//                    r_c[tm][tn] += r_comp_a[tm] * r_comp_b[tn];
//                }
//            }
//        }
//
//        __syncthreads();
//    }
//
//#pragma unroll
//    for (int i = 0; i < TM / 2; i++) {
//        int store_c_gmem_m = by * BM + ty * TM / 2 + i;
//        int store_c_gmem_n = bx * BN + tx * TN / 2;
//        int store_c_gmem_addr = OFFSET(store_c_gmem_m, store_c_gmem_n, N);
//        FLOAT4(c[store_c_gmem_addr]) = FLOAT4(r_c[i][0]);
//        FLOAT4(c[store_c_gmem_addr + BN / 2]) = FLOAT4(r_c[i][4]);
//    }
//#pragma unroll
//    for (int i = 0; i < TM / 2; i++) {
//        int store_c_gmem_m = by * BM + BM / 2 + ty * TM / 2 + i;
//        int store_c_gmem_n = bx * BN + tx * TN / 2;
//        int store_c_gmem_addr = OFFSET(store_c_gmem_m, store_c_gmem_n, N);
//        FLOAT4(c[store_c_gmem_addr]) = FLOAT4(r_c[i + TM / 2][0]);
//        FLOAT4(c[store_c_gmem_addr + BN / 2]) = FLOAT4(r_c[i + TM / 2][4]);
//    }
//}

#define MAT_OFFSET()

//<<<(GeometryNum, RowNum, 1), (ColNum, 1, 1)>>>
//
template <typename T, uint16_t ColNum>
__global__ static void IGA_matvec_A_multiply_v0(const gpmesh::SplineContext context,
                                                   const gpmesh::AssemblyContext assembly_context,
                                                   const gpmesh::SplineAttributeBase<T> in,
                                                   const gpmesh::SplineAttributeBase<T> out
)
{
    using namespace gpmesh;

    //load all information in shrd_mem
    uint32_t geometry_id = blockIdx.x;
    uint32_t row_id = blockIdx.y;

    uint16_t col_id = threadIdx.x;

    if (geometry_id >= context.get_num_geometry() || col_id >= ColNum)
    {
        return;
    }

    AssemblyInfo assembly_info = assembly_context.get_assembly_info()[geometry_id];
    double* localMats = assembly_info.local_mats;

    assert(assembly_info.geometry_id == geometry_id);

    __shared__ T assembly_data[ColNum];
    __shared__ T in_data[ColNum];
    __shared__ T res;

    res = 0;

    assembly_data[col_id] = localMats[row_id + 192 * col_id];

    uint16_t ci = col_id % 64;
    uint16_t cd = col_id / 64;

    in_data[col_id] = in(geometry_id, ci, cd);

    __syncthreads();

    ci = row_id % 64;
    cd = row_id / 64;
}

//<<<(GeometryNum, 1, 1), (RowNum, 1, 1)>>>
//合并访存，输入数据 shared memory优化
template <typename T, uint16_t n>
__global__ static void IGA_matvec_A_multiply_v1(const gpmesh::SplineContext context,
                                                const gpmesh::AssemblyContext assembly_context,
                                                const gpmesh::SplineAttributeBase<T> in,
                                                const gpmesh::SplineAttributeBase<T> out
)
{
    using namespace gpmesh;
    uint32_t geometry_id = blockIdx.x;
    uint16_t row_id = threadIdx.x;

    if (geometry_id >= context.get_num_geometry() || row_id >= 192)
    {
        return;
    }

    AssemblyInfo assembly_info = assembly_context.get_assembly_info()[geometry_id];
    T* localMats = assembly_info.local_mats;

    assert(assembly_info.geometry_id == geometry_id);

    __shared__ T s_in[n];

    s_in[row_id] = in.RawOperator(geometry_id, row_id);

    __syncthreads();

    T res = 0;

    #pragma unroll
    for(uint16_t col_id = 0; col_id < 192; ++col_id)
    {
        res += localMats[row_id + 192 * col_id] * s_in[col_id];
    }

    out.RawOperator(geometry_id, row_id) = res;
}

template <typename T, uint16_t n>
__global__ static void IGA_matvec_A_multiply_v1(const gpmesh::SplineContext context,
                                                const gpmesh::AssemblyContext assembly_context,
                                                const gpmesh::SplineAttributeGlobal<T> in,
                                                const gpmesh::SplineAttributeGlobal<T> out
)
{
    using namespace gpmesh;
    uint32_t geometry_id = blockIdx.x;
    uint16_t row_id = threadIdx.x;

    if (geometry_id >= context.get_num_geometry() || row_id >= 192)
    {
        return;
    }

    AssemblyInfo assembly_info = assembly_context.get_assembly_info()[geometry_id];
    T* localMats = assembly_info.local_mats;

    assert(assembly_info.geometry_id == geometry_id);

    __shared__ T s_in[n];

    s_in[row_id] = in.RawOperator(geometry_id, row_id);

    __syncthreads();

    T res = 0;

#pragma unroll
    for(uint16_t col_id = 0; col_id < 192; ++col_id)
    {
        res += localMats[row_id + 192 * col_id] * s_in[col_id];
    }

    out.RawOperator(geometry_id, row_id) = res;
}


//<<<(GeometryNum, 1, 1), (RowNum, 1, 1)>>>
//合并访存, 取消共享内存优化
template <typename T, uint16_t n>
__global__ static void IGA_matvec_A_multiply_v2(const gpmesh::SplineContext context,
                                                const gpmesh::AssemblyContext assembly_context,
                                                const gpmesh::SplineAttributeBase<T> in,
                                                const gpmesh::SplineAttributeBase<T> out
)
{
    using namespace gpmesh;
    uint32_t geometry_id = blockIdx.x;
    uint16_t row_id = threadIdx.x;

    if (geometry_id >= context.get_num_geometry() || row_id >= 192)
    {
        return;
    }

    AssemblyInfo assembly_info = assembly_context.get_assembly_info()[geometry_id];
    T* localMats = assembly_info.local_mats;

    assert(assembly_info.geometry_id == geometry_id);

    T res = 0;

#pragma unroll
    for(uint16_t col_id = 0; col_id < 192; ++col_id)
    {
        res += localMats[row_id + 192 * col_id] * in.RawOperator(geometry_id, col_id);
    }

    out.RawOperator(geometry_id, row_id) = res;
}

//初始化 PR
template<typename T>
__global__ static void init_PR(const gpmesh::SplineAttributeGlobal<T> B,
                               const gpmesh::SplineAttributeGlobal<T> S,
                               const gpmesh::SplineAttributeGlobal<T> P,
                               const gpmesh::SplineAttributeGlobal<T> R,
                               uint32_t num_cpoints, uint32_t num_attributes)
{
    const size_t tid = blockDim.x * blockIdx.x + threadIdx.x;

    if(tid < 4)
    {
        P.m_d_global_attr[tid] = 0.0;
    }
}

template <typename T>
__global__ static void init_PR(T* B, T* S, T* P, T* R, uint32_t num_elements)
{
    const size_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid < num_elements)
    {
        //printf("B[%d] = %f in gpu\n", tid, B[tid]);
        P[tid] = B[tid] - S[tid];
        R[tid] = P[tid];
    }
}

//Y = alpha * X + beta * Y
template <typename T>
__global__ void axpy(const gpmesh::SplineAttributeGlobal<T> Y,
          const gpmesh::SplineAttributeGlobal<T> X,
          const T                           alpha,
          const T                           beta,
          const uint32_t                    num_elements,
          cudaStream_t                      stream = NULL)
{
    assert(X.m_num_cpoints == Y.m_num_cpoints);
    size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    T thread_val = 0;

    while(tid < num_elements)
    {

        Y.m_d_global_attr[tid] = alpha * X.m_d_global_attr[tid] + beta * Y.m_d_global_attr[tid];

        tid += blockDim.x * gridDim.x;
    }
}




