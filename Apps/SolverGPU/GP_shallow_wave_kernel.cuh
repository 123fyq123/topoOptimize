#pragma once

#include "GP_Attribute.h"
#include "GP_Context.h"
#include "Kernels/GP_Query_Dispatcher.cuh"
#include "Util/GP_Vector.h"
#include "GP_Mesh_Writer.h"

inline float rand_FloatRange(float a, float b)
{
    return ((b - a) * ((float)rand() / RAND_MAX)) + a;
}

template <typename T, uint32_t blockThreads>
__global__ static void  shallow_wave(const gpmesh::Context      context,
                                             gpmesh::VertexAttribute<T> old_h,
                                             gpmesh::VertexAttribute<T> h,
                                             gpmesh::VertexAttribute<T> new_h,
                                             float damping = 0.996f,
                                             float rate = 0.005f
                                    )
{
    using namespace gpmesh;
    auto shallow_wave_lambda = [&](gpmesh::VertexHandle vertex_id, VertexIterator& fv) {
        // get the face's three vertices coordinates

        uint16_t thread_iter_time = 1;
        for (uint16_t k = 0; k < thread_iter_time; ++k)
        {
            new_h(vertex_id) = h(vertex_id) + (h(vertex_id) - old_h(vertex_id)) * damping;

            uint16_t neighbor_v_size = fv.size();

            for (uint16_t i = 0; i < neighbor_v_size; ++i)
            {
                new_h(vertex_id) += (h(fv[i]) - h(vertex_id)) * rate;
            }

//            old_h(vertex_id) = h(vertex_id);
//
//            __syncthreads();
//            h(vertex_id) = new_h(vertex_id);
        }
    };

    query_block_dispatcher<Op::VV, blockThreads>(context, shallow_wave_lambda);
}

