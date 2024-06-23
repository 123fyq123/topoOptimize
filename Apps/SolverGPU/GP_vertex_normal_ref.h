#pragma once
#include <vector>
#include "Util/GP_Report.h"

template <typename T>
inline void vertex_normal_ref(const std::vector<std::vector<uint32_t>>& Faces,
                              const std::vector<std::vector<T>>&        Verts,
                              std::vector<T>& vertex_normal)
{
    uint32_t num_vertices = Verts.size();
    uint32_t num_faces    = Faces.size();

    memset((void*)vertex_normal.data(), 0, vertex_normal.size() * sizeof(T));

    T        edge_len[4];
    uint32_t v[4];
    T        fn[4][3];

    auto l2_norm_sq = [](const T ax0,
                         const T ax1,
                         const T ax2,
                         const T bx0,
                         const T bx1,
                         const T bx2) {
        // compute (xa0-xb0)*(xa0-xb0) + (xa1-xb1)*(xa1-xb1) +
        // (xa2-xb2)*(xa2-xb2)
        T x0 = ax0 - bx0;
        T x1 = ax1 - bx1;
        T x2 = ax2 - bx2;
        return x0 * x0 + x1 * x1 + x2 * x2;
    };

    auto cross_product =
        [](T xv1, T yv1, T zv1, T xv2, T yv2, T zv2, T& xx, T& yy, T& zz) {
            xx = yv1 * zv2 - zv1 * yv2;
            yy = zv1 * xv2 - xv1 * zv2;
            zz = xv1 * yv2 - yv1 * xv2;
        };

    for (uint32_t f = 0; f < num_faces; ++f) {

        v[0] = Faces[f][0];
        v[1] = Faces[f][1];
        v[2] = Faces[f][2];
        v[3] = Faces[f][3];

        cross_product(Verts[v[1]][0] - Verts[v[0]][0],
                      Verts[v[1]][1] - Verts[v[0]][1],
                      Verts[v[1]][2] - Verts[v[0]][2],
                      Verts[v[3]][0] - Verts[v[0]][0],
                      Verts[v[3]][1] - Verts[v[0]][1],
                      Verts[v[3]][2] - Verts[v[0]][2],
                      fn[0][0],
                      fn[0][1],
                      fn[0][2]);

        cross_product(Verts[v[2]][0] - Verts[v[1]][0],
                      Verts[v[2]][1] - Verts[v[1]][1],
                      Verts[v[2]][2] - Verts[v[1]][2],
                      Verts[v[0]][0] - Verts[v[1]][0],
                      Verts[v[0]][1] - Verts[v[1]][1],
                      Verts[v[0]][2] - Verts[v[1]][2],
                      fn[1][0],
                      fn[1][1],
                      fn[1][2]);

        cross_product(Verts[v[3]][0] - Verts[v[2]][0],
                      Verts[v[3]][1] - Verts[v[2]][1],
                      Verts[v[3]][2] - Verts[v[2]][2],
                      Verts[v[1]][0] - Verts[v[2]][0],
                      Verts[v[1]][1] - Verts[v[2]][1],
                      Verts[v[1]][2] - Verts[v[2]][2],
                      fn[2][0],
                      fn[2][1],
                      fn[2][2]);

        cross_product(Verts[v[2]][0] - Verts[v[3]][0],
                      Verts[v[2]][1] - Verts[v[3]][1],
                      Verts[v[2]][2] - Verts[v[3]][2],
                      Verts[v[0]][0] - Verts[v[3]][0],
                      Verts[v[0]][1] - Verts[v[3]][1],
                      Verts[v[0]][2] - Verts[v[3]][2],
                      fn[3][0],
                      fn[3][1],
                      fn[3][2]);

        edge_len[0] = l2_norm_sq(Verts[v[0]][0],
                                 Verts[v[0]][1],
                                 Verts[v[0]][2],
                                 Verts[v[1]][0],
                                 Verts[v[1]][1],
                                 Verts[v[1]][2]);  // v0-v1

        edge_len[1] = l2_norm_sq(Verts[v[1]][0],
                                 Verts[v[1]][1],
                                 Verts[v[1]][2],
                                 Verts[v[2]][0],
                                 Verts[v[2]][1],
                                 Verts[v[2]][2]);  // v1-v2

        edge_len[2] = l2_norm_sq(Verts[v[2]][0],
                                 Verts[v[2]][1],
                                 Verts[v[2]][2],
                                 Verts[v[3]][0],
                                 Verts[v[3]][1],
                                 Verts[v[3]][2]);  // v2-v3


        edge_len[3] = l2_norm_sq(Verts[v[3]][0],
                                 Verts[v[3]][1],
                                 Verts[v[3]][2],
                                 Verts[v[0]][0],
                                 Verts[v[0]][1],
                                 Verts[v[0]][2]);  // v3-v0


        for (uint32_t i = 0; i < 4; ++i) {
            uint32_t k    = (i + 3) % 3;
            uint32_t base = 3 * v[i];

            for (uint32_t l = 0; l < 3; ++l) {
                vertex_normal[base + l] += fn[i][l] / (edge_len[i] + edge_len[k]);
            }
        }
    }
}


template <typename T>
inline void vertex_neighbor_sizes_ref(
    const std::vector<std::vector<uint32_t>>& Faces,
    const std::vector<std::vector<T>>&        Verts,
    std::vector<T>& vertex_neighbor_sizes
    )
{
    uint32_t num_vertices = Verts.size();
    uint32_t num_faces    = Faces.size();
    memset(vertex_neighbor_sizes.data(), 0, vertex_neighbor_sizes.size() * sizeof(T));

    auto l2_norm_sq = [](const T ax0,
                         const T ax1,
                         const T ax2,
                         const T bx0,
                         const T bx1,
                         const T bx2) {
        // compute (xa0-xb0)*(xa0-xb0) + (xa1-xb1)*(xa1-xb1) +
        // (xa2-xb2)*(xa2-xb2)
        T x0 = ax0 - bx0;
        T x1 = ax1 - bx1;
        T x2 = ax2 - bx2;
        return x0 * x0 + x1 * x1 + x2 * x2;
    };

    for (int i = 0; i < num_faces; ++i)
    {
        auto vids = Faces[i];

        std::vector<float> f_length(4, 0);

        for (int j = 0; j < vids.size(); ++j)
        {
            f_length[j] = l2_norm_sq(
                Verts[vids[j]][0],
                Verts[vids[j]][1],
                Verts[vids[j]][2],
                0,
                0,
                0);

//            f_length[j] =
//                Verts[vids[j]][0] +
//                Verts[vids[j]][1] +
//                Verts[vids[j]][2];
        }

        for (int j = 0; j < vids.size(); ++j)
        {
            vertex_neighbor_sizes[vids[j]] +=
                f_length[0] + f_length[1] + f_length[2] + f_length[3];
        }
    }

}


template <typename T>
inline void vertex_normal_ref_v2(const std::vector<std::vector<uint32_t>>& Faces,
                              const std::vector<std::vector<T>>&        Verts,
                              std::vector<T>& vertex_normal)
{
    uint32_t num_vertices = Verts.size();
    uint32_t num_faces    = Faces.size();

    memset((void*)vertex_normal.data(), 0, vertex_normal.size() * sizeof(T));

    T        edge_len[4];
    uint32_t v[4];
    T        fn[4][3];

    auto l2_norm_sq = [](const T ax0,
                         const T ax1,
                         const T ax2,
                         const T bx0,
                         const T bx1,
                         const T bx2) {
        // compute (xa0-xb0)*(xa0-xb0) + (xa1-xb1)*(xa1-xb1) +
        // (xa2-xb2)*(xa2-xb2)
        T x0 = ax0 - bx0;
        T x1 = ax1 - bx1;
        T x2 = ax2 - bx2;
        return x0 * x0 + x1 * x1 + x2 * x2;
    };

    auto cross_product =
        [](T xv1, T yv1, T zv1, T xv2, T yv2, T zv2, T& xx, T& yy, T& zz) {
            xx = yv1 * zv2 - zv1 * yv2;
            yy = zv1 * xv2 - xv1 * zv2;
            zz = xv1 * yv2 - yv1 * xv2;
        };

    for (uint32_t f = 0; f < num_faces; ++f) {

        v[0] = Faces[f][0];
        v[1] = Faces[f][1];
        v[2] = Faces[f][2];
        v[3] = Faces[f][3];

        cross_product(Verts[v[1]][0] - Verts[v[0]][0],
                      Verts[v[1]][1] - Verts[v[0]][1],
                      Verts[v[1]][2] - Verts[v[0]][2],
                      Verts[v[3]][0] - Verts[v[0]][0],
                      Verts[v[3]][1] - Verts[v[0]][1],
                      Verts[v[3]][2] - Verts[v[0]][2],
                      fn[0][0],
                      fn[0][1],
                      fn[0][2]);

        cross_product(Verts[v[2]][0] - Verts[v[1]][0],
                      Verts[v[2]][1] - Verts[v[1]][1],
                      Verts[v[2]][2] - Verts[v[1]][2],
                      Verts[v[0]][0] - Verts[v[1]][0],
                      Verts[v[0]][1] - Verts[v[1]][1],
                      Verts[v[0]][2] - Verts[v[1]][2],
                      fn[1][0],
                      fn[1][1],
                      fn[1][2]);

        cross_product(Verts[v[3]][0] - Verts[v[2]][0],
                      Verts[v[3]][1] - Verts[v[2]][1],
                      Verts[v[3]][2] - Verts[v[2]][2],
                      Verts[v[1]][0] - Verts[v[2]][0],
                      Verts[v[1]][1] - Verts[v[2]][1],
                      Verts[v[1]][2] - Verts[v[2]][2],
                      fn[2][0],
                      fn[2][1],
                      fn[2][2]);

        cross_product(Verts[v[0]][0] - Verts[v[3]][0],
                      Verts[v[0]][1] - Verts[v[3]][1],
                      Verts[v[0]][2] - Verts[v[3]][2],
                      Verts[v[2]][0] - Verts[v[3]][0],
                      Verts[v[2]][1] - Verts[v[3]][1],
                      Verts[v[2]][2] - Verts[v[3]][2],
                      fn[3][0],
                      fn[3][1],
                      fn[3][2]);

        edge_len[0] = l2_norm_sq(Verts[v[0]][0],
                                 Verts[v[0]][1],
                                 Verts[v[0]][2],
                                 Verts[v[1]][0],
                                 Verts[v[1]][1],
                                 Verts[v[1]][2]);  // v0-v1

        edge_len[1] = l2_norm_sq(Verts[v[1]][0],
                                 Verts[v[1]][1],
                                 Verts[v[1]][2],
                                 Verts[v[2]][0],
                                 Verts[v[2]][1],
                                 Verts[v[2]][2]);  // v1-v2

        edge_len[2] = l2_norm_sq(Verts[v[2]][0],
                                 Verts[v[2]][1],
                                 Verts[v[2]][2],
                                 Verts[v[3]][0],
                                 Verts[v[3]][1],
                                 Verts[v[3]][2]);  // v2-v3


        edge_len[3] = l2_norm_sq(Verts[v[3]][0],
                                 Verts[v[3]][1],
                                 Verts[v[3]][2],
                                 Verts[v[0]][0],
                                 Verts[v[0]][1],
                                 Verts[v[0]][2]);  // v3-v0


        for (uint32_t i = 0; i < 4; ++i) {
            //uint32_t k    = (i + 3) % 3;
            uint32_t base = 3 * v[i];

            for (uint32_t l = 0; l < 3; ++l) {
                vertex_normal[base + l] += fn[i][l];
            }
        }
    }
}