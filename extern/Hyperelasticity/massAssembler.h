#pragma once
#include"SparseSystem.h"
#include"ToolFunc.h"
class massAssembler{
	double density;
	std::vector<CBSplineSolid> geometry;
	SparseSystem m_system;
	Eigen::SparseMatrix<double>  m_matrix;
public:
	massAssembler(std::vector<CBSplineSolid> _geometry, double density, SparseSystem system);
	void assemble();
    void assemble(std::vector<Eigen::MatrixXd>& mass_mats);
	const Eigen::SparseMatrix<double> matrix() const;
};