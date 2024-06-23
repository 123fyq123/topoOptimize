#pragma once

#include <string>
#include <vector>
#include <fstream>
#include "Util/GP_Log.h"
#include "cuda_runtime.h"

template <typename T, int N>
void outputVerticesInfo2File(std::string filename, std::vector<T> vertices_Info)
{
    if (vertices_Info.size() % N != 0)
    {
        GPMESH_ERROR("vertices info layout error");
        return;
    }
    int v_size = vertices_Info.size() / N;

    std::ofstream output(filename);

    if(!output.is_open())
    {
        GPMESH_ERROR("file is not open");
        return;
    }

    for (int i = 0; i < v_size; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            output << vertices_Info[i * N + j] << " ";
        }

        output << std::endl;
    }

    output.close();
}
