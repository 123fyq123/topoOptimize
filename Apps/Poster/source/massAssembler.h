#pragma once
#include"SparseSystem.h"
#include"ToolFunc.h"
namespace Poster_Needed {
    class massAssembler {
        double density;
        std::vector<Poster_Needed::CBSplineSolid> geometry;
        sparseSystem m_system;
        Eigen::SparseMatrix<double> m_matrix;
    public:
        massAssembler(std::vector<Poster_Needed::CBSplineSolid> _geometry, double density, sparseSystem system);

        void assemble();

        const Eigen::SparseMatrix<double> matrix() const;
    };
}