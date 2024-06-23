#pragma once

#include "GP_Attribute.h"
#include "GP_Context.h"
#include "Kernels/GP_Query_Dispatcher.cuh"
#include "Util/GP_Vector.h"
#include "GP_Mesh_Writer.h"

/**
 * partial_voronoi_area()
 */
template <typename T>
__host__ __device__ __forceinline__ T
partial_voronoi_area(const gpmesh::Vector<3, T>& p,  // center
                     const gpmesh::Vector<3, T>& q,  // before center
                     const gpmesh::Vector<3, T>& r)  // after center

{
    // compute partial Voronoi area of the center vertex that is associated with
    // the triangle p->q->r (oriented ccw)
    using namespace gpmesh;

    // Edge vector p->q
    const Vector<3, T> pq = q - p;

    // Edge vector q->r
    const Vector<3, T> qr = r - q;

    // Edge vector p->r
    const Vector<3, T> pr = r - p;

    // compute and check triangle area
    T triangle_area = cross(pq, pr).norm();

    if (triangle_area <= std::numeric_limits<T>::min()) {
        return -1;
    }


    // dot products for each corner (of its two emanating edge vectors)
    T dotp = dot(pq, pr);
    T dotq = -dot(qr, pq);
    T dotr = dot(qr, pr);
    if (dotp < 0.0) {
        return 0.25 * triangle_area;
    }

    // angle at q or r obtuse
    else if (dotq < 0.0 || dotr < 0.0) {
        return 0.125 * triangle_area;
    }

    // no obtuse angles
    else {
        // cot(angle) = cos(angle)/sin(angle) = dot(A,B)/norm(cross(A,B))
        T cotq = dotq / triangle_area;
        T cotr = dotr / triangle_area;

        // clamp cot(angle) by clamping angle to [1,179]
        clamp_cot(cotq);
        clamp_cot(cotr);

        return 0.125 * (pr.norm2() * cotq + pq.norm2() * cotr);
    }

    return -1;
}

template <typename T>
__device__ __forceinline__ T
partial_voronoi_area(const gpmesh::VertexHandle&       p_id,  // center
                     const gpmesh::VertexHandle&       q_id,  // before center
                     const gpmesh::VertexHandle&       r_id,  // after center
                     const gpmesh::VertexAttribute<T>& X)
{
    // compute partial Voronoi area of the center vertex that is associated with
    // the triangle p->q->r (oriented ccw)
    using namespace gpmesh;
    
    const Vector<3, T> p(X(p_id, 0), X(p_id, 1), X(p_id, 2));
    const Vector<3, T> q(X(q_id, 0), X(q_id, 1), X(q_id, 2));
    const Vector<3, T> r(X(r_id, 0), X(r_id, 1), X(r_id, 2));

    return partial_voronoi_area(p, q, r);
}


template <typename T, uint32_t blockThreads>
__global__  static void init_B(const gpmesh::Context      context,
                              gpmesh::VertexAttribute<T> X,
                              gpmesh::VertexAttribute<T> B,
                              const bool use_uniform_laplace)
{
    using namespace gpmesh;

    auto init_lambda = [&](VertexHandle& p_id, const VertexIterator& iter) {
        if (use_uniform_laplace) {
            const T valence = static_cast<T>(iter.size());
            B(p_id, 0)      = X(p_id, 0) * valence;
            B(p_id, 1)      = X(p_id, 1) * valence;
            B(p_id, 2)      = X(p_id, 2) * valence;
        } else {

            // using Laplace weights
            T v_weight = 0;

            // this is the last vertex in the one-ring (before r_id)
            VertexHandle q_id = iter.back();

            for (uint32_t v = 0; v < iter.size(); ++v) {
                // the current one ring vertex
                VertexHandle r_id = iter[v];

                T tri_area = partial_voronoi_area(p_id, q_id, r_id, X);

                v_weight += (tri_area > 0) ? tri_area : 0.0;

                q_id = r_id;
            }
            v_weight = 0.5 / v_weight;

            B(p_id, 0) = X(p_id, 0) / v_weight;
            B(p_id, 1) = X(p_id, 1) / v_weight;
            B(p_id, 2) = X(p_id, 2) / v_weight;
        }
    };

    //if not use_uniform_laplace
    //need iter neighboring vertex in order


    query_block_dispatcher<Op::VV, blockThreads>(
        context, init_lambda, !use_uniform_laplace);
}

template <typename T, uint32_t blockThreads>
__global__  static void test_V_Surround_V(const gpmesh::Context      context,
                              gpmesh::VertexAttribute<T> X,
                              gpmesh::VertexAttribute<T> B,
                              const bool use_uniform_laplace)
{
    using namespace gpmesh;

    auto init_lambda = [&](VertexHandle& p_id, const VertexIterator& iter) {
        if (use_uniform_laplace) {
            const T valence = static_cast<T>(iter.size());
            B(p_id, 0)      = X(p_id, 0) * valence;
            B(p_id, 1)      = X(p_id, 1) * valence;
            B(p_id, 2)      = X(p_id, 2) * valence;
        } else {

            // using Laplace weights
            T v_weight = 0;

            // this is the last vertex in the one-ring (before r_id)
            VertexHandle q_id = iter.back();

            for (uint32_t v = 0; v < iter.size(); ++v) {
                // the current one ring vertex
                VertexHandle r_id = iter[v];

                T tri_area = partial_voronoi_area(p_id, q_id, r_id, X);

                v_weight += (tri_area > 0) ? tri_area : 0.0;

                q_id = r_id;
            }
            v_weight = 0.5 / v_weight;

            B(p_id, 0) = X(p_id, 0) / v_weight;
            B(p_id, 1) = X(p_id, 1) / v_weight;
            B(p_id, 2) = X(p_id, 2) / v_weight;
        }
    };

    //if not use_uniform_laplace
    //need iter neighboring vertex in order


    query_block_dispatcher<Op::V_SURROUND_V, blockThreads>(
        context, init_lambda, !use_uniform_laplace);
}