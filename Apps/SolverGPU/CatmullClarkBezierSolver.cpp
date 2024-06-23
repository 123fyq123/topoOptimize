#include "CatmullClarkBezierSolver.h"
#include <assert.h>
#include <iostream>
#include "Timer.h"
#include <omp.h>
#include <set>


void transformGradients(const Eigen::MatrixXd& jac, int q, const Eigen::MatrixXd& allGrads, Eigen::MatrixXd& trfGradsK, int dim)
{
	const int numGrads = allGrads.rows() / dim;
	const Eigen::Map<const Eigen::MatrixXd> grads_k(allGrads.col(q).data(), dim, numGrads);
	trfGradsK.noalias() = jac.inverse() * grads_k;
}

CatmullClarkBezierLinerElasticlSolver::CatmullClarkBezierLinerElasticlSolver(CatmullClarkBezier _model, IGA::CCBezierBdy _bdy) : m_model(_model), m_bdy(_bdy){
	getC();
	std::cout << "getC end\n";
	initAll();
	std::cout << "initAll end\n";
};
void CatmullClarkBezierLinerElasticlSolver::pushToMatrix(int patch, const std::vector<int>& localIndex, const Eigen::MatrixXd& localMats) {
	int N_D = localIndex.size();
	for (int m = 0; m < N_D; ++m) {
		for (int n = 0; n < N_D; ++n) {
			for (int di = 0; di < 3; di++) {
				for (int dj = 0; dj < 3; dj++) {

					int r = m_model.global_id[patch][m] + di * m_model.control_vertices.size();
					int c = m_model.global_id[patch][n] + dj * m_model.control_vertices.size();

					m_bigmatrix[r][m_colIDs[r][c]] += localMats(N_D * di + m, N_D * dj + n);
				}
			}
		}
	}
}
void CatmullClarkBezierLinerElasticlSolver::initAll() {
	m_colIDs.resize(m_model.control_vertices.size() * 3);
	m_bigmatrix.resize(m_model.control_vertices.size() * 3);
	m_bigrhs.resize(m_model.control_vertices.size() * 3, 0);
	std::cout << "m_model.control_vertices.size() = " << m_model.control_vertices.size() <<"\n";

	for (int i = 0; i < m_model.global_id.size(); ++i) {

		int N_D = 64;

		for (int m = 0; m < N_D; ++m) {
			for (int n = 0; n < N_D; ++n) {
				for (int di = 0; di < 3; di++) {
					for (int dj = 0; dj < 3; dj++) {

						int r = m_model.global_id[i][m] + di * m_model.control_vertices.size();
						int c = m_model.global_id[i][n] + dj * m_model.control_vertices.size();
						
						if (!m_colIDs[r].count(c)) {
							m_colIDs[r][c] = m_bigmatrix[r].size();
							m_bigmatrix[r].push_back(0.0);
						}
					}
				}
			}
		}
	}
}
void CatmullClarkBezierLinerElasticlSolver::Solving() {
	
	std::vector<double> p_u{ 0.0694318442,0.3300094782,0.6699905218,0.9305681558 }; 
	std::vector<double> p_v{ 0.0694318442,0.3300094782,0.6699905218,0.9305681558 };
	std::vector<double> p_w{ 0.0694318442,0.3300094782,0.6699905218,0.9305681558 };

	std::vector<double> w_u{ 0.3478548451 / 2 ,0.6521451549 / 2, 0.6521451549 / 2, 0.3478548451 / 2};
	std::vector<double> w_v{ 0.3478548451 / 2 ,0.6521451549 / 2, 0.6521451549 / 2, 0.3478548451 / 2};
	std::vector<double> w_w{ 0.3478548451 / 2 ,0.6521451549 / 2, 0.6521451549 / 2, 0.3478548451 / 2};

	std::vector<double> uknots = { 0, 1 };
	std::vector<double> vknots = { 0, 1 };
	std::vector<double> wknots = { 0, 1 };

	int gaussNum_u = 4;
	int gaussNum_v = 4;
	int gaussNum_w = 4;
	int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;

	std::vector<Eigen::MatrixXd> basisValue;
	PreComputeBasisValue(p_u, p_v, p_w, basisValue);
	// ?????????
	std::cout << "bigmatrix constructing\n";
//#pragma omp parallel for
	for (int patch_id = 0; patch_id < m_model.global_id.size(); ++patch_id) {
		int local_basis_cnt = 64;
		auto& local_index = m_model.global_id[patch_id];
		Eigen::MatrixXd localMats;
		localMats.setZero(3 * local_basis_cnt, 3 * local_basis_cnt);

		for (int kk = 0; kk < gaussNum_w; ++kk) {
			for (int jj = 0; jj < gaussNum_v; ++jj) {
				for (int ii = 0; ii < gaussNum_u; ++ii) {
					int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
					Eigen::MatrixXd jac, P;
					getJacobi(basisValue[1].col(q), local_index, jac);
					
					double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jac.determinant());
					Eigen::MatrixXd physGrad;
					transformGradients(jac, q, basisValue[1], physGrad, 3);
					// ?????????? ???
					double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;
					for (int j = 0; j < local_basis_cnt; j++) {
						m_Rx = physGrad.col(j)(0);
						m_Ry = physGrad.col(j)(1);
						m_Rz = physGrad.col(j)(2);

						
						for (int k = 0; k < local_basis_cnt; k++) {
							m_Qx = physGrad.col(k)(0);
							m_Qy = physGrad.col(k)(1);
							m_Qz = physGrad.col(k)(2);

							// localMats(j, k) += weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
							localMats(j, k) += weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
							//localMats(j, k + N_D) += weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
							localMats(j, k + local_basis_cnt) += weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
							//localMats(j, k + 2 * N_D) += weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qy);
							localMats(j, k + 2 * local_basis_cnt) += weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qx);
							//localMats(j + N_D, k) += weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
							localMats(j + local_basis_cnt, k) += weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
							//localMats(j + N_D, k + N_D) += weight * (C(0, 0) * m_Ry * m_Qy + C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
							localMats(j + local_basis_cnt, k + local_basis_cnt) += weight * (C(0, 0) * m_Ry * m_Qy + C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
							//localMats(j + N_D, k + 2 * N_D) += weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
							localMats(j + local_basis_cnt, k + 2 * local_basis_cnt) += weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
							//localMats(j + 2 * N_D, k) += weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
							localMats(j + 2 * local_basis_cnt, k) += weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
							//localMats(j + 2 * N_D, k + N_D) += weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
							localMats(j + 2 * local_basis_cnt, k + local_basis_cnt) += weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
							//localMats(j + 2 * N_D, k + 2 * N_D) += weight * (C(0, 0) * m_Rz * m_Qz + C(3, 3) * (m_Ry * m_Qy + m_Rx * m_Qx));
							localMats(j + 2 * local_basis_cnt, k + 2 * local_basis_cnt) += weight * (C(0, 0) * m_Rz * m_Qz + C(3, 3) * (m_Ry * m_Qy + m_Rx * m_Qx));

						}
					}

				}
			}
		}

//#pragma omp critical
		pushToMatrix(patch_id, local_index, localMats);
	}


	std::cout << "rhs constructing\n";
	int global_dof_cnt = 3 * m_model.control_vertices.size();
	// ?????????
	for (auto& bdy : m_bdy) {
		// ??????
		if (bdy.bdy_type == IGA::BdyBase::Neumann) {
			int cur_global_id = bdy.global_base_id;
			for (int i = 0; i < 3; ++i) {
				m_bigrhs[cur_global_id + i * m_model.control_vertices.size()] = bdy.bdy_value[i];
			}
		}
	}

	std::cout << "set Dirichlet\n";
	// ????????????
	std::set<int> dof_dirichlet;
	for (auto& bdy : m_bdy) {
		// ??????
		if (bdy.bdy_type == IGA::BdyBase::Dirichlet) {
			int cur_global_id = bdy.global_base_id;
			for (int i = 0; i < 3; ++i) {
				int cur_dof_id = cur_global_id + i * m_model.control_vertices.size();
				m_bigrhs[cur_dof_id] = bdy.bdy_value[i];

				dof_dirichlet.insert(cur_dof_id);

				// ???????0??1
				for (auto mp : m_colIDs[cur_dof_id]) {
					int o_dof_id = mp.first;
					int k = mp.second;
					if (cur_dof_id == o_dof_id) m_bigmatrix[cur_dof_id][k] = 1;
					else m_bigmatrix[cur_dof_id][k] = 0;

				}
			}
		}
	}
	// ???????0
	for (int row_dof_id = 0; row_dof_id < 3 * m_model.control_vertices.size(); ++row_dof_id) {
		
		for (auto mp : m_colIDs[row_dof_id]) {
			int o_dof_id = mp.first;
			int k = mp.second;
			if (dof_dirichlet.count(o_dof_id) && row_dof_id != o_dof_id) {
				m_bigmatrix[row_dof_id][k] = 0;
			}

		}
	}

	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	std::vector<T> coeffs;
	for (int i = 0; i < global_dof_cnt; ++i) {
		for (auto mp : m_colIDs[i]) {
			int o_dof_id = mp.first;
			int k = mp.second;
			coeffs.push_back(T(i, o_dof_id, m_bigmatrix[i][k]));
		}
	}
	{
		Timer timer;
		std::cout << "prepare solving...\n";
		Eigen::setNbThreads(8);  // ?????С
		Eigen::initParallel();   // ??????????


		std::cout << "Current threads: " << Eigen::nbThreads() << std::endl;
		SpMat A(global_dof_cnt, global_dof_cnt);
		A.setFromTriplets(coeffs.begin(), coeffs.end());
		Eigen::VectorXd b = Eigen::VectorXd::Map(m_bigrhs.data(), m_bigrhs.size());

		Eigen::PardisoLU<SpMat> solver;
		std::cout << "start solving...\n";
		solver.compute(A);
		//solver.compute(A.selfadjointView<Eigen::Lower>());
		/*solver.analyzePattern(A);
		solver.factorize(A)*/;

		Eigen::VectorXd eigen_u = solver.solve(b);

		m_u = std::vector<double>(eigen_u.data(), eigen_u.data() + eigen_u.size());

		std::cout << "end solving...\n";
	}
}
void CatmullClarkBezierLinerElasticlSolver::SolvingWithDensity(const std::vector<double>& density, std::vector<std::vector<std::vector<double>>>& local_matrixs) {

	std::vector<double> p_u{ 0.0694318442,0.3300094782,0.6699905218,0.9305681558 };
	std::vector<double> p_v{ 0.0694318442,0.3300094782,0.6699905218,0.9305681558 };
	std::vector<double> p_w{ 0.0694318442,0.3300094782,0.6699905218,0.9305681558 };

	std::vector<double> w_u{ 0.3478548451 / 2 ,0.6521451549 / 2, 0.6521451549 / 2, 0.3478548451 / 2 };
	std::vector<double> w_v{ 0.3478548451 / 2 ,0.6521451549 / 2, 0.6521451549 / 2, 0.3478548451 / 2 };
	std::vector<double> w_w{ 0.3478548451 / 2 ,0.6521451549 / 2, 0.6521451549 / 2, 0.3478548451 / 2 };

	std::vector<double> uknots = { 0, 1 };
	std::vector<double> vknots = { 0, 1 };
	std::vector<double> wknots = { 0, 1 };

	int gaussNum_u = 4;
	int gaussNum_v = 4;
	int gaussNum_w = 4;
	int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;

	std::vector<Eigen::MatrixXd> basisValue;
	PreComputeBasisValue(p_u, p_v, p_w, basisValue);
	// ?????????
	std::cout << "bigmatrix constructing\n";
	local_matrixs.resize(m_model.global_id.size());
	//#pragma omp parallel for
	for (int patch_id = 0; patch_id < m_model.global_id.size(); ++patch_id) {
		int local_basis_cnt = 64;
		auto& local_index = m_model.global_id[patch_id];
		Eigen::MatrixXd localMats;
		localMats.setZero(3 * local_basis_cnt, 3 * local_basis_cnt);

		double rho = (0.0001 + (1 - 0.0001) * pow(density[patch_id], 3));
		for (int kk = 0; kk < gaussNum_w; ++kk) {
			for (int jj = 0; jj < gaussNum_v; ++jj) {
				for (int ii = 0; ii < gaussNum_u; ++ii) {
					int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
					Eigen::MatrixXd jac, P;
					getJacobi(basisValue[1].col(q), local_index, jac);

					double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jac.determinant());
					weight *= rho;
					Eigen::MatrixXd physGrad;
					transformGradients(jac, q, basisValue[1], physGrad, 3);
					// ?????????? ???
					double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;
					for (int j = 0; j < local_basis_cnt; j++) {
						m_Rx = physGrad.col(j)(0);
						m_Ry = physGrad.col(j)(1);
						m_Rz = physGrad.col(j)(2);


						for (int k = 0; k < local_basis_cnt; k++) {
							m_Qx = physGrad.col(k)(0);
							m_Qy = physGrad.col(k)(1);
							m_Qz = physGrad.col(k)(2);

							//localMats(j, k) += weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
							localMats(j, k) += weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
							//localMats(j, k + N_D) += weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
							localMats(j, k + local_basis_cnt) += weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
							//localMats(j, k + 2 * N_D) += weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qy);
							localMats(j, k + 2 * local_basis_cnt) += weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qx);
							//localMats(j + N_D, k) += weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
							localMats(j + local_basis_cnt, k) += weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
							//localMats(j + N_D, k + N_D) += weight * (C(0, 0) * m_Ry * m_Qy + C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
							localMats(j + local_basis_cnt, k + local_basis_cnt) += weight * (C(0, 0) * m_Ry * m_Qy + C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
							//localMats(j + N_D, k + 2 * N_D) += weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
							localMats(j + local_basis_cnt, k + 2 * local_basis_cnt) += weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
							//localMats(j + 2 * N_D, k) += weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
							localMats(j + 2 * local_basis_cnt, k) += weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
							//localMats(j + 2 * N_D, k + N_D) += weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
							localMats(j + 2 * local_basis_cnt, k + local_basis_cnt) += weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
							//localMats(j + 2 * N_D, k + 2 * N_D) += weight * (C(0, 0) * m_Rz * m_Qz + C(3, 3) * (m_Ry * m_Qy + m_Rx * m_Qx));
							localMats(j + 2 * local_basis_cnt, k + 2 * local_basis_cnt) += weight * (C(0, 0) * m_Rz * m_Qz + C(3, 3) * (m_Ry * m_Qy + m_Rx * m_Qx));

						}
					}

				}
			}
		}
		local_matrixs[patch_id].resize(64 * 3, std::vector<double>(64 * 3, 0));

		for (int i = 0; i < 64 * 3; ++i) {
			for (int j = 0; j < 64 * 3; ++j) {
				local_matrixs[patch_id][i][j] = localMats(i, j);
			}
		}
		
		//#pragma omp critical
		pushToMatrix(patch_id, local_index, localMats);
	}


	std::cout << "rhs constructing\n";
	int global_dof_cnt = 3 * m_model.control_vertices.size();
	// ?????????
	for (auto& bdy : m_bdy) {
		// ??????
		if (bdy.bdy_type == IGA::BdyBase::Neumann) {
			int cur_global_id = bdy.global_base_id;
			for (int i = 0; i < 3; ++i) {
				m_bigrhs[cur_global_id + i * m_model.control_vertices.size()] = bdy.bdy_value[i];
			}
		}
	}

	std::cout << "set Dirichlet\n";
	// ????????????
	std::set<int> dof_dirichlet;
	for (auto& bdy : m_bdy) {
		// ??????
		if (bdy.bdy_type == IGA::BdyBase::Dirichlet) {
			int cur_global_id = bdy.global_base_id;
			for (int i = 0; i < 3; ++i) {
				int cur_dof_id = cur_global_id + i * m_model.control_vertices.size();
				m_bigrhs[cur_dof_id] = bdy.bdy_value[i];

				dof_dirichlet.insert(cur_dof_id);

				// ???????0??1
				for (auto mp : m_colIDs[cur_dof_id]) {
					int o_dof_id = mp.first;
					int k = mp.second;
					if (cur_dof_id == o_dof_id) m_bigmatrix[cur_dof_id][k] = 1;
					else m_bigmatrix[cur_dof_id][k] = 0;

				}
			}
		}
	}
	// ???????0
	for (int row_dof_id = 0; row_dof_id < 3 * m_model.control_vertices.size(); ++row_dof_id) {

		for (auto mp : m_colIDs[row_dof_id]) {
			int o_dof_id = mp.first;
			int k = mp.second;
			if (dof_dirichlet.count(o_dof_id) && row_dof_id != o_dof_id) {
				m_bigmatrix[row_dof_id][k] = 0;
			}

		}
	}

	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	std::vector<T> coeffs;
	for (int i = 0; i < global_dof_cnt; ++i) {
		for (auto mp : m_colIDs[i]) {
			int o_dof_id = mp.first;
			int k = mp.second;
			coeffs.push_back(T(i, o_dof_id, m_bigmatrix[i][k]));
		}
	}
	{
		Timer timer;
		std::cout << "prepare solving...\n";
		Eigen::setNbThreads(8);  // ?????С
		Eigen::initParallel();   // ??????????


		std::cout << "Current threads: " << Eigen::nbThreads() << std::endl;
		SpMat A(global_dof_cnt, global_dof_cnt);
		A.setFromTriplets(coeffs.begin(), coeffs.end());
		Eigen::VectorXd b = Eigen::VectorXd::Map(m_bigrhs.data(), m_bigrhs.size());

		Eigen::PardisoLU<SpMat> solver;
		std::cout << "start solving...\n";
		solver.compute(A);
		//solver.compute(A.selfadjointView<Eigen::Lower>());
		/*solver.analyzePattern(A);
		solver.factorize(A)*/;

		Eigen::VectorXd eigen_u = solver.solve(b);

		m_u = std::vector<double>(eigen_u.data(), eigen_u.data() + eigen_u.size());

		std::cout << "end solving...\n";
	}
}
void CatmullClarkBezierLinerElasticlSolver::getJacobi(
	const Eigen::MatrixXd& basisValue,
	const std::vector<int>& localIndex,
	Eigen::MatrixXd& jac) {
	int len = localIndex.size();

	const Eigen::Map<const Eigen::MatrixXd> Nder(basisValue.data(), 3, len);
	Eigen::MatrixXd X;

	
	X.resize(localIndex.size(), 3);
	for (size_t i = 0; i < localIndex.size(); ++i) {
		X.row(i) = m_model.control_vertices[localIndex[i]];
	}
	jac = Nder * X;
}


void CatmullClarkBezierLinerElasticlSolver::PreComputeBasisValue(  
							std::vector<double>& u,
							std::vector<double>& v, 
							std::vector<double>& w, 
							std::vector<Eigen::MatrixXd>& basisValue) {
	basisValue.resize(2);
	basisValue[0].resize(64, 64);
	basisValue[1].resize(3 * 64, 64);
	std::vector<Eigen::VectorXd> singleValue;

	for (int k = 0; k < w.size(); ++k) {
		for (int j = 0; j < v.size(); ++j) {
			for (int i = 0; i < u.size(); ++i) {
				EvalCatmullClarkBezier2Ders(u[i], v[j], w[k], singleValue);
				basisValue[0].col(k * u.size() * v.size() + j * u.size() + i) = singleValue[0];
				basisValue[1].col(k * u.size() * v.size() + j * u.size() + i) = singleValue[1];
			}
		}
	}
}

double** init2DArray(int x, int y)
{
	/*
	* Allocate space for two-dimensional array
	*/

	double** array = (double**)malloc(x * sizeof(double*));

	int c;
	for (c = 0; c < x; c++)
	{
		array[c] = (double*)malloc(y * sizeof(double));
	}

	return array;
}

void free2Darray(double** array, int x)
{
	/*
	* Free space for two-dimensional array
	*/

	int c;
	for (c = 0; c < x; c++)
	{
		free(array[c]);
	}
	free(array);
}

void nrbAllBsplineBasisFuns(int i, double u, int p, double* knot, double* N)
{
	/*
	* Calculate all non zero basis functions at point u
	*/

	int j, r;
	double saved, temp;
	double* lef = (double*)malloc(sizeof(double) * (p + 1));
	double* rig = (double*)malloc(sizeof(double) * (p + 1));

	N[0] = 1.0;
	for (j = 1; j <= p; j++)
	{
		lef[j] = u - knot[i + 1 - j];
		rig[j] = knot[i + j] - u;
		saved = 0.0;
		for (r = 0; r < j; r++)
		{
			temp = N[r] / (rig[r + 1] + lef[j - r]);
			N[r] = saved + rig[r + 1] * temp;
			saved = lef[j - r] * temp;
		}
		N[j] = saved;
	}

	free(lef);
	free(rig);
}

void nrbAllBsplineBasisFunsDers(int i, double u, int p, double* knot, int nd, double** ders)
{
	/*
	* Calculate all non-zero derivatives of the b-spline functions
	*/

	int j, r, k, s1, s2, rk, pk, j1, j2;
	double saved, temp, d;
	double* lef = (double*)malloc(sizeof(double) * (p + 1));
	double* rig = (double*)malloc(sizeof(double) * (p + 1));
	double** ndu = init2DArray(p + 1, p + 1);
	double** a = init2DArray(p + 1, p + 1);

	ndu[0][0] = 1.0;
	for (j = 1; j <= p; j++)
	{
		lef[j] = u - knot[i + 1 - j];
		rig[j] = knot[i + j] - u;
		saved = 0.0;
		for (r = 0; r < j; r++)
		{
			ndu[j][r] = rig[r + 1] + lef[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + rig[r + 1] * temp;
			saved = lef[j - r] * temp;
		}
		ndu[j][j] = saved;
	}

	for (j = 0; j <= p; j++)
		ders[0][j] = ndu[j][p];

	for (r = 0; r <= p; r++)
	{
		s1 = 0; s2 = 1;
		a[0][0] = 1.0;

		for (k = 1; k <= nd; k++)
		{
			d = 0.0;
			rk = r - k;
			pk = p - k;
			if (r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
			if (rk >= -1)  j1 = 1;
			else          j1 = -rk;
			if (r - 1 <= pk) j2 = k - 1;
			else          j2 = p - r;

			for (j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}

			if (r <= pk)
			{
				a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			ders[k][r] = d;
			j = s1; s1 = s2; s2 = j;
		}
	}

	r = p;
	for (k = 1; k <= nd; k++)
	{
		for (j = 0; j <= p; j++)
			ders[k][j] *= r;
		r *= (p - k);
	}

	free(lef);
	free(rig);

	free2Darray(ndu, p + 1);
	free2Darray(a, p + 1);
}

void CatmullClarkBezierLinerElasticlSolver::EvalCatmullClarkBezier2Ders(double& u,
	double& v,
	double& w,
	std::vector<Eigen::VectorXd>& singleValue) 
{
	int tol = 1e-15;
	int u_ctrlPtsNum = 3;
	int v_ctrlPtsNum = 3;
	int w_ctrlPtsNum = 3;
	singleValue.resize(2);

	int uspan = 3;
	int vspan = 3;
	int wspan = 3;
	
	double* Nu = (double*)malloc(sizeof(double) * (3 + 1));
	double* Nv = (double*)malloc(sizeof(double) * (3 + 1));
	double* Nw = (double*)malloc(sizeof(double) * (3 + 1));
	double** Nu_ders = init2DArray(u_ctrlPtsNum + 1, 3 + 1);
	double** Nv_ders = init2DArray(v_ctrlPtsNum + 1, 3 + 1);
	double** Nw_ders = init2DArray(w_ctrlPtsNum + 1, 3 + 1);

	double u_knots[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
	double v_knots[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
	double w_knots[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };

	nrbAllBsplineBasisFuns(uspan, u, 3, u_knots, Nu);
	nrbAllBsplineBasisFuns(vspan, v, 3, v_knots, Nv);
	nrbAllBsplineBasisFuns(wspan, w, 3, w_knots, Nw);
	nrbAllBsplineBasisFunsDers(uspan, u, 3, u_knots, u_ctrlPtsNum, Nu_ders);
	nrbAllBsplineBasisFunsDers(vspan, v, 3, v_knots, v_ctrlPtsNum, Nv_ders);
	nrbAllBsplineBasisFunsDers(wspan, w, 3, w_knots, w_ctrlPtsNum, Nw_ders);

	int vind, uind, wind;
	double wgt = 1, we = 0.0, fac, dwdxi = 0.0, dwdeta = 0.0, dwdgemma = 0.0;
	std::vector<double> res;
	for (int k = 0; k <= 3; k++) {
		wind = wspan - 3 + k;
		for (int j = 0; j <= 3; j++)
		{
			vind = vspan - 3 + j;
			for (int i = 0; i <= 3; i++)
			{
				uind = uspan - 3 + i;
				we += Nu[i] * Nv[j] * Nw[k] * wgt;
				dwdxi += Nu_ders[1][i] * Nv[j] * Nw[k] * wgt;
				dwdeta += Nu[i] * Nv_ders[1][j] * Nw[k] * wgt;
				dwdgemma += Nu[i] * Nv[j] * Nw_ders[1][k] * wgt;
			}
		}
	}
	singleValue[0].resize((3 + 1) * (3 + 1) * (3 + 1));
	singleValue[1].resize(3 * (3 + 1) * (3 + 1) * (3 + 1));
	int k = 0;
	for (int t = 0; t <= 3; t++) {
		wind = wspan - 3 + t;
		for (int j = 0; j <= 3; j++)
		{
			vind = vspan - 3 + j;
			for (int i = 0; i <= 3; i++)
			{
				uind = uspan - 3 + i;
				fac = wgt / (we * we);
				singleValue[0][k] = Nu[i] * Nv[j] * Nw[t] * fac * we;
				singleValue[1][3 * k] = (Nu_ders[1][i] * Nv[j] * Nw[t] * we - Nu[i] * Nv[j] * Nw[t] * dwdxi) * fac;
				singleValue[1][3 * k + 1] = (Nu[i] * Nv_ders[1][j] * Nw[t] * we - Nu[i] * Nv[j] * Nw[t] * dwdeta) * fac;
				singleValue[1][3 * k + 2] = (Nu[i] * Nv[j] * Nw_ders[1][t] * we - Nu[i] * Nv[j] * Nw[t] * dwdgemma) * fac;
				++k;
			}
		}
	}
	free(Nu);
	free(Nv);
	free(Nw);
	free2Darray(Nu_ders, (u_ctrlPtsNum + 1));
	free2Darray(Nv_ders, (v_ctrlPtsNum + 1));
	free2Darray(Nw_ders, (w_ctrlPtsNum + 1));
}

Eigen::Vector3d CatmullClarkBezierLinerElasticlSolver::getU(int patch_id, double u, double v, double w) {
	Eigen::Vector3d xx;
	return xx;
}
void CatmullClarkBezierLinerElasticlSolver::getC() {
	Eigen::MatrixXd R = Eigen::MatrixXd::Identity(3, 3);
	double a = lambda, b = mu;
	short a_dim = R.cols();
	short dimTensor = a_dim * (a_dim + 1) / 2;
	C.resize(dimTensor, dimTensor);
	// voigt indices
	Eigen::MatrixX2i v(dimTensor, 2);
	if (a_dim == 2)
		v << 0, 0,
		1, 1,
		0, 1;
	if (a_dim == 3)
		v << 0, 0,
		1, 1,
		2, 2,
		0, 1,
		1, 2,
		0, 2;
	// componentwise definition
	for (short i = 0; i < dimTensor; ++i)
		for (short j = 0; j < dimTensor; ++j)
			C(i, j) = a * R(v(i, 0), v(i, 1)) * R(v(j, 0), v(j, 1)) +
			b * (R(v(i, 0), v(j, 0)) * R(v(i, 1), v(j, 1)) +
				R(v(i, 0), v(j, 1)) * R(v(i, 1), v(j, 0)));
}