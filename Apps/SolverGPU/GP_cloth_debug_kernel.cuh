#pragma once

#include "GP_Attribute.h"
#include "GP_Context.h"
#include "Kernels/GP_Query_Dispatcher.cuh"
#include "Util/GP_Vector.h"

#include <Eigen/Core>

template <typename T, uint32_t blockThreads>
__global__ static void  cloth_debug_kernel(const gpmesh::Context  context,
                                          gpmesh::VertexAttribute<T> surround_size,
                                          gpmesh::VertexAttribute<uint32_t> surround_info)
{
    using namespace gpmesh;
    auto test = [&](gpmesh::VertexHandle vertex_handle, VertexIterator& vv_iter) {
        int n = vv_iter.size();

        surround_size(vertex_handle) = n;
        std::pair<uint32_t, uint16_t> center_pair = vertex_handle.unpack();

        //printf("center vid-%u\n", center_pair.second);

        for (int i = 0; i < n; ++i)
        {
            std::pair<uint32_t, uint16_t> p_l_pair = vv_iter[i].unpack();
            if (p_l_pair.second == INVALID16)
            {
                surround_size(vertex_handle)--;
                continue;
            }

            atomicAdd(&surround_info(vertex_handle, i), p_l_pair.first * 100000 + p_l_pair.second);
            //atomicAdd(&surround_info(vertex_handle, i), p_l_pair.second);
        }
    };

    query_block_dispatcher<Op::V_V_SURROUND_UNCLOSED, blockThreads>(context, test, true);
}

template <typename T, uint32_t blockThreads>
__global__ static void  cloth_explicit_integration_kernel(const gpmesh::Context  context,
                                          gpmesh::VertexAttribute<T> x,
                                          gpmesh::VertexAttribute<T> ox,
                                          gpmesh::VertexAttribute<T> v,
                                          gpmesh::VertexAttribute<T> f,
                                          float dt,
                                          float damping,
                                          float stiffness
                                                         )
{
    using namespace gpmesh;
    auto explicit_lambda = [&](gpmesh::VertexHandle vertex_handle, VertexIterator& vv_iter) {
        int n = vv_iter.size();

        v(vertex_handle, 0) *= expf(-dt * damping);

        Eigen::Vector3f force(0.0f, 0.0f, -9.8f);

        for (uint32_t i = 0; i < n; ++i)
        {
            std::pair<uint32_t, uint16_t> p_l_pair = vv_iter[i].unpack();
            if (p_l_pair.second == INVALID16) continue;

            Eigen::Vector3f disp(0.0f, 0.0f, 0.0f);
            Eigen::Vector3f rest_disp(0.0f, 0.0f, 0.0f);
            for (uint32_t dim = 0; dim < 3; ++dim)
            {
                disp[dim] = x(vertex_handle, dim) - x(vv_iter[i], dim);
                rest_disp[dim] = ox(vertex_handle, dim) - ox(vv_iter[i], dim);
            }
            force += -stiffness * (disp.norm() - rest_disp.norm()) * disp.normalized();
        }

        for (uint32_t dim = 0; dim < 3; ++dim)
        {
            v(vertex_handle, dim) += dt * force[dim];
        }
    };

    query_block_dispatcher<Op::V_V_SURROUND_UNCLOSED, blockThreads>(context, explicit_lambda, true);
}