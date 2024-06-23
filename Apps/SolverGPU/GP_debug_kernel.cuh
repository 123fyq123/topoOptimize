#pragma once

#include "GP_Attribute.h"
#include "GP_Context.h"
#include "Kernels/GP_Query_Dispatcher.cuh"
#include "Util/GP_Vector.h"

template <typename T, uint32_t blockThreads>
__global__ static void  debug_kernel(const gpmesh::Context  context)
{
    using namespace gpmesh;
    auto test = [&](gpmesh::VertexHandle vertex_handle, VertexIterator& vv_iter) {

    };

    query_block_dispatcher<Op::V_SURROUND_V, blockThreads>(context, test, true);
}

template <typename T, uint32_t blockThreads>
__global__ static void  debug_surround_kernel(const gpmesh::Context  context,
                                                 gpmesh::VertexAttribute<T> surround_size)
{
    using namespace gpmesh;
    auto test = [&](gpmesh::VertexHandle vertex_handle, VertexIterator& vv_iter) {
        surround_size(vertex_handle) = vv_iter.size();
    };

    query_block_dispatcher<Op::V_SURROUND_V, blockThreads>(context, test, true);
}

template <typename T, uint32_t blockThreads>
__global__ static void  debug_surround_kerne_v2(const gpmesh::Context  context,
                                               gpmesh::VertexAttribute<T> surround_size,
                                               gpmesh::VertexAttribute<uint32_t> surround_info)
{
    using namespace gpmesh;
    auto test = [&](gpmesh::VertexHandle vertex_handle, VertexIterator& vv_iter) {
        int n = vv_iter.size();

        surround_size(vertex_handle) = n;

        for (int i = 0; i < n; ++i)
        {
            std::pair<uint32_t, uint16_t> p_l_pair = vv_iter[i].unpack();
            atomicAdd(&surround_info(vertex_handle, i), p_l_pair.first * 100000 + p_l_pair.second);
            //atomicAdd(&surround_info(vertex_handle, i), p_l_pair.second);
        }
    };

    query_block_dispatcher<Op::V_SURROUND_V, blockThreads>(context, test, true);
}
