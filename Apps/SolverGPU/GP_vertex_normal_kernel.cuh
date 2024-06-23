#pragma once

#include "GP_Attribute.h"
#include "GP_Context.h"
#include "Kernels/GP_Query_Dispatcher.cuh"
#include "Util/GP_Vector.h"

template <typename T, uint32_t blockThreads>
__global__ static void  compute_vertex_normal(const gpmesh::Context      context,
                                             gpmesh::VertexAttribute<T> coords,
                                             gpmesh::VertexAttribute<T> normals)
{
    using namespace gpmesh;
    auto vn_lambda = [&](gpmesh::FaceHandle face_id, gpmesh::VertexIterator& fv) {
        // get the face's three vertices coordinates
        Vector<3, T> c0(coords(fv[0], 0), coords(fv[0], 1), coords(fv[0], 2));
        Vector<3, T> c1(coords(fv[1], 0), coords(fv[1], 1), coords(fv[1], 2));
        Vector<3, T> c2(coords(fv[2], 0), coords(fv[2], 1), coords(fv[2], 2));
        Vector<3, T> c3(coords(fv[3], 0), coords(fv[3], 1), coords(fv[3], 2));

        // the four edges length
        Vector<4, T> l(dist2(c0, c1), dist2(c1, c2), dist2(c2, c3), dist2(c3, c0));

        Vector<3, T> n = cross(c1 - c0, c3 - c0);

        for(uint32_t i = 0; i < 3; ++i)
        {
            atomicAdd(&normals(fv[0], i), n[i] / (l[0] + l[3]));
        }

        n = cross(c2 - c1, c0 - c1);

        for(uint32_t i = 0; i < 3; ++i)
        {
            atomicAdd(&normals(fv[1], i), n[i] / (l[1] + l[0]));
        }

        n = cross(c3 - c2, c1 - c2);

        for(uint32_t i = 0; i < 3; ++i)
        {
            atomicAdd(&normals(fv[2], i), n[i] / (l[2] + l[1]));
        }

        n = cross(c2 - c3, c0 - c3);

        for(uint32_t i = 0; i < 3; ++i)
        {
            atomicAdd(&normals(fv[3], i), n[i] / (l[3] + l[2]));
        }
    };

    query_block_dispatcher<Op::FV, blockThreads>(context, vn_lambda);
}

template <typename T, uint32_t blockThreads>
__global__ static void  compute_vertex_normal_v2(const gpmesh::Context      context,
                                             gpmesh::VertexAttribute<T> coords,
                                             gpmesh::VertexAttribute<T> normals)
{
    using namespace gpmesh;
    auto vn_lambda = [&](gpmesh::FaceHandle face_id, VertexIterator& fv) {
        // get the face's three vertices coordinates
        Vector<3, T> c0(coords(fv[0], 0), coords(fv[0], 1), coords(fv[0], 2));
        Vector<3, T> c1(coords(fv[1], 0), coords(fv[1], 1), coords(fv[1], 2));
        Vector<3, T> c2(coords(fv[2], 0), coords(fv[2], 1), coords(fv[2], 2));
        Vector<3, T> c3(coords(fv[3], 0), coords(fv[3], 1), coords(fv[3], 2));

        Vector<3, T> n0 = cross(c1 - c0, c3 - c0);
        Vector<3, T> n1 = cross(c2 - c1, c0 - c1);
        Vector<3, T> n2 = cross(c3 - c2, c1 - c2);
        Vector<3, T> n3 = cross(c0 - c3, c2 - c3);

        // the four edges length
        Vector<4, T> l(dist2(c0, c1), dist2(c1, c2), dist2(c2, c3), dist2(c3, c0));

        for(uint32_t i = 0; i < 3; ++i)
        {
            atomicAdd(&normals(fv[0], i), n0[i]);
            atomicAdd(&normals(fv[1], i), n1[i]);
            atomicAdd(&normals(fv[2], i), n2[i]);
            atomicAdd(&normals(fv[3], i), n3[i]);
        }
    };

    query_block_dispatcher<Op::FV, blockThreads>(context, vn_lambda);
}

template <typename T, uint32_t blockThreads>
__global__ static void  compute_vertex_neighbor_size(const gpmesh::Context      context,
                                             gpmesh::VertexAttribute<T> coords,
                                             gpmesh::VertexAttribute<T> neighbor_sizes)
{
    using namespace gpmesh;
    auto neighbor_length_lambda = [&](gpmesh::FaceHandle face_id, VertexIterator& fv) {
        // get the face's three vertices coordinates
        Vector<3, T> c0(coords(fv[0], 0), coords(fv[0], 1), coords(fv[0], 2));
        Vector<3, T> c1(coords(fv[1], 0), coords(fv[1], 1), coords(fv[1], 2));
        Vector<3, T> c2(coords(fv[2], 0), coords(fv[2], 1), coords(fv[2], 2));
        Vector<3, T> c3(coords(fv[3], 0), coords(fv[3], 1), coords(fv[3], 2));

        Vector<3, T> origin(0, 0, 0);
        // the four edges length
        Vector<4, T> l(dist2(c0, origin),
                       dist2(c1, origin),
                       dist2(c2, origin),
                       dist2(c3, origin));

//        Vector<4, T> l(c0[0] + c0[1] + c0[2],
//                       c1[0] + c1[1] + c1[2],
//                       c2[0] + c2[1] + c2[2],
//                       c3[0] + c3[1] + c3[2]);


        atomicAdd(&neighbor_sizes(fv[0], 0), (l[0] + l[1] + l[2] + l[3]));
        atomicAdd(&neighbor_sizes(fv[1], 0), (l[0] + l[1] + l[2] + l[3]));
        atomicAdd(&neighbor_sizes(fv[2], 0), (l[0] + l[1] + l[2] + l[3]));
        atomicAdd(&neighbor_sizes(fv[3], 0), (l[0] + l[1] + l[2] + l[3]));


    };

    query_block_dispatcher<Op::FV, blockThreads>(context, neighbor_length_lambda);
}

template <typename T, uint32_t blockThreads>
__global__ static void  compute_vertex_topology(const gpmesh::Context      context,
                                                    gpmesh::VertexAttribute<T> coords,
                                                    gpmesh::VertexAttribute<T> topology_info)
{
    using namespace gpmesh;
    auto neighbor_length_lambda = [&](gpmesh::FaceHandle face_id, VertexIterator& fv) {
        // get the face's three vertices coordinates
        Vector<3, T> c0(coords(fv[0], 0), coords(fv[0], 1), coords(fv[0], 2));
        Vector<3, T> c1(coords(fv[1], 0), coords(fv[1], 1), coords(fv[1], 2));
        Vector<3, T> c2(coords(fv[2], 0), coords(fv[2], 1), coords(fv[2], 2));
        Vector<3, T> c3(coords(fv[3], 0), coords(fv[3], 1), coords(fv[3], 2));

        Vector<3, T> origin(0, 0, 0);
        // the four edges length
        //        Vector<4, T> l(dist2(c0, origin),
        //                       dist2(c1, origin),
        //                       dist2(c2, origin),
        //                       dist2(c3, origin));

        Vector<4, T> l(c0[0] + c0[1] + c0[2],
                       c1[0] + c1[1] + c1[2],
                       c2[0] + c2[1] + c2[2],
                       c3[0] + c3[1] + c3[2]);

        float t = fv.size() + fv[0].unpack().first + fv[0].unpack().second
                  + fv[1].unpack().first + fv[1].unpack().second
                  + fv[2].unpack().first + fv[2].unpack().second
                  + fv[3].unpack().first + fv[3].unpack().second;

        atomicAdd(&topology_info(fv[0], 0), t);
        atomicAdd(&topology_info(fv[1], 0), t);
        atomicAdd(&topology_info(fv[2], 0), t);
        atomicAdd(&topology_info(fv[3], 0), t);
    };

    query_block_dispatcher<Op::FV, blockThreads>(context, neighbor_length_lambda);
}

template <typename T, uint32_t blockThreads>
__global__ static void  compute_face_vertex_order(const gpmesh::Context      context,
                                               gpmesh::VertexAttribute<T> face_vids)
{
    using namespace gpmesh;
    auto face_vids_lambda = [&](gpmesh::FaceHandle face_id, VertexIterator& fv) {

        for (int i = 0; i < 4; ++i)
        {
            std::pair<uint32_t, uint16_t> p_l_pair = fv[i].unpack();
            atomicAdd(&face_vids(fv[i], 0), p_l_pair.first * 100000 + p_l_pair.second);
        }
    };

    query_block_dispatcher<Op::FV, blockThreads>(context, face_vids_lambda);
}

template <typename T, uint32_t blockThreads>
__global__ static void  compute_fvids(const gpmesh::Context      context,
                                                 gpmesh::FaceAttribute<T> fvids)
{
    using namespace gpmesh;
    auto fvids_lambda = [&](gpmesh::FaceHandle face_id, VertexIterator& fv) {

        for (int i = 0; i < 4; ++i)
        {
            std::pair<uint32_t, uint16_t> p_l_pair = fv[i].unpack();
            atomicAdd(&fvids(face_id, i), p_l_pair.first * 100000 + p_l_pair.second);
        }
    };

    query_block_dispatcher<Op::FV, blockThreads>(context, fvids_lambda);
}