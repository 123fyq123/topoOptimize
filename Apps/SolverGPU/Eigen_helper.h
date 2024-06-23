#pragma once

#include <vector>
#include <Eigen/Core>

template<typename T>
void convert_eigen_matrix_to_vector(Eigen::MatrixX<T>& mat, std::vector<T>& vec)
{
    vec.resize(mat.rows() * mat.cols());

    for (int i = 0; i < mat.rows(); ++i)
    {
        for (int j = 0; j < mat.cols(); ++j)
        {
            vec[i * mat.cols() + j] = mat.operator()(i, j);
        }
    }
}

template<typename T>
void convert_displacement_coefs_to_tripleVec(Eigen::MatrixX<T>& mat,
                                             std::vector<std::vector<std::vector<T>>>& vec)
{
    int n = mat.rows();
    int num_attribute = mat.cols();

    const int local_size = 64;

    assert(n % local_size == 0);
    int geometry_size = n / local_size;

    vec = std::vector<std::vector<std::vector<T>>>(geometry_size,
                                                   std::vector<std::vector<T>>(local_size, std::vector<T>(num_attribute)));

    for (int g = 0; g < geometry_size; ++g)
    {
        for (int local_idx = 0; local_idx < local_size; ++local_idx)
        {
            for (int a = 0; a < num_attribute; ++a)
            {
                vec[g][local_idx][a] = mat(local_idx + local_size * g, a);
            }
        }
    }

    return;
}
