#pragma once
#include"stiffAssembler.h"
class Hyperelastic {

public:

    bool bUsePoint = true;

public:
	Hyperelastic(stiffAssembler* assembleInstance);

	Hyperelastic(stiffAssembler* assembleInstance, const Eigen::MatrixXd initSolution);

	std::string status();

	void reset();

	bool computeUpdate();

	void solve();

	void setSolverType(int solveType);

	solveInfo defaultInfo();

	void setMaterial(const pair<int, vector<double>> material);

	Eigen::MatrixXd solution();

	void setInfoNone();

protected:
	stiffAssembler* m_assembleInstance = nullptr;
	Eigen::MatrixXd m_solVector;

	size_t m_numIterations;

	newton_status m_status;

	std::pair<int, vector<double>> m_material;
	
	const double tol = 1e-10;

	solveInfo info;

	std::ofstream myfile;

public:
    double m_residualNorm;

    double m_initResidualNorm;

    double m_updateNorm;

    double m_initUpdateNorm;
};