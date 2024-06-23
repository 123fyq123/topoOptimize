#pragma once
#include <iostream>
#include <fstream>
#include <vector>

inline void write_mesh(std::vector<std::vector<float>> vertices, std::vector<std::vector<uint32_t>> fv, std::string filepath)
{
    std::ofstream output(filepath);

    if (!output.is_open())
    {
        throw std::runtime_error("cannot open file: " + filepath);
    }

    for (int i = 0; i < vertices.size(); ++i)
    {
        output << "v";

        for (int j = 0; j < vertices[i].size(); ++j)
        {
            output << " " << vertices[i][j];
        }

        output << std::endl;
    }

    for (int i = 0; i < fv.size(); ++i)
    {
        output << "f";

        for (int j = 0; j < fv[i].size(); ++j)
        {
            output << " " << fv[i][j] + 1;
        }

        output << std::endl;
    }

    output.close();
}