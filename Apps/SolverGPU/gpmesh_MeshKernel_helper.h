#pragma once

#include <MeshKernel/Mesh.h>
#include <vector>

inline void convert_to_mesh_kernel(std::vector<std::vector<float>>& vertices, std::vector<std::vector<uint32_t>>& fv, MeshKernel::SurfaceMesh& mesh)
{
    std::vector<Vex> _vertices;
    _vertices.reserve(vertices.size());

    for (int i = 0 ; i < vertices.size(); ++i)
    {
        _vertices.emplace_back(vertices[i][0], vertices[i][1], vertices[i][2]);
    }

    std::vector<std::vector<VH>> _faces;
    _faces.resize(fv.size());
    for (int i = 0; i < fv.size(); ++i) {
        for (int j = 0; j < fv[i].size(); ++j)
        {
            _faces[i].push_back(VH(fv[i][j]));
        }
    }

    mesh = MeshKernel::SurfaceMesh(_vertices, _faces);
}