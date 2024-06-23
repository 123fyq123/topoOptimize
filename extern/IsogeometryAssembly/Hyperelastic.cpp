#include"Hyperelastic.h"

Hyperelastic::Hyperelastic(stiffAssembler* assembleInstance):m_assembleInstance(assembleInstance),info(defaultInfo()) {
	m_solVector.setZero(assembleInstance->numDofs(), 1);
	reset();
}

Hyperelastic::Hyperelastic(stiffAssembler* assembleInstance, const Eigen::MatrixXd initSolution): m_assembleInstance(assembleInstance),
info(defaultInfo()) {
	m_solVector = initSolution;
	reset();
}

std::string Hyperelastic::status()
{
	std::string statusString;
	if (m_status == newton_status::converged)
		statusString = "牛顿迭代法在 " +
		std::to_string(m_numIterations) + "次迭代后收敛.";
	else if (m_status == newton_status::interrupted)
		statusString = "牛顿迭代法在 " +
		std::to_string(m_numIterations) + "次迭代后中断.";
	else if (m_status == newton_status::bad_solution)
		statusString = "牛顿迭代法在 " +
		std::to_string(m_numIterations) + " 次迭代后中断，原因是当前求得解无效";
	else if (m_status == newton_status::working)
		statusString = "迭代次数: " + std::to_string(m_numIterations) +
		", 更新向量范数: " + std::to_string(m_updateNorm) +
		", 与初始向量比: " + std::to_string(m_updateNorm / m_initUpdateNorm) +
		", 残差向量范数: " + std::to_string(m_residualNorm) +
		", 与初始残差比: " + std::to_string(m_residualNorm / m_initResidualNorm);
	return statusString;
}

void Hyperelastic::reset()
{
	m_status = newton_status::working;
	m_numIterations = 0;
	m_residualNorm = 0;
	m_initResidualNorm = 1;
	m_updateNorm = 0;
	m_initUpdateNorm = 1;
	m_material = { material_law::saint_venant_kirchhoff,{0.} };
	m_assembleInstance->setMaterialLaw(m_material);
}

bool Hyperelastic::computeUpdate()
{
	m_assembleInstance->assemble(m_solVector);
	Eigen::VectorXd updateVector;
	if (info.SolveType == linear_solver::LU) {
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(m_assembleInstance->getMatrix());
		updateVector = solver.solve(m_assembleInstance->getRhs());
	}
	if (info.SolveType == linear_solver::LDLT) {
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(m_assembleInstance->getMatrix());
		updateVector = solver.solve(m_assembleInstance->getRhs());
	}
	m_solVector += updateVector;
	m_residualNorm = m_assembleInstance->getRhs().norm();
	m_updateNorm = updateVector.norm();

	if (m_numIterations == 0)
	{
		m_initUpdateNorm = m_updateNorm;
		m_initResidualNorm = m_residualNorm;
	}

	m_numIterations++;

	return true;
}

void Hyperelastic::solve()
{
	reset();
	int maxIter = info.MaxIters;
	while (m_status == newton_status::working && m_numIterations < maxIter) {
		if (!computeUpdate()) {       // 求解失败
			m_status = newton_status::bad_solution;
			m_numIterations = maxIter + 1;
			break;
		}
		if (info.Vervosity == newton_verbosity::all)
			cout << status() << endl;
		if (m_residualNorm < info.AbsTol || m_updateNorm < info.AbsTol ||
			m_residualNorm / m_initResidualNorm < info.RelTol ||m_updateNorm / m_initUpdateNorm < info.RelTol)
			m_status = newton_status::converged;
	}

	if (info.Vervosity != newton_verbosity::none) cout << status() << endl;
}

void Hyperelastic::setSolverType(int solveType)
{
	info.SolveType = solveType;
}

solveInfo Hyperelastic::defaultInfo()
{
	solveInfo info;
	info.AbsTol = 1e-12;
	info.RelTol = 1e-9;
	info.MaxIters = 50;
	info.SolveType = Solver::LDLT;
	info.Vervosity = newton_verbosity::all;
	return info;
}

void Hyperelastic::setMaterial(const pair<int, vector<double>> material)
{
	m_material = material;
	m_assembleInstance->setMaterialLaw(m_material);
}

Eigen::MatrixXd Hyperelastic::solution()
{
	return m_solVector;
}
