#include "stiffAssembler.h"

stiffAssembler::stiffAssembler(std::vector<CBSplineSolid> m_patch, std::vector<boundary_condition> m_dirichletSides, std::vector<boundary_condition> m_neumannSides, std::vector<interFace> m_inf)
    :geometry(m_patch),dirichletSides(m_dirichletSides), neumanSides(m_neumannSides)
{
	globalIndex.resize(geometry.size());
	computerIndex();
	gl_coefs_index.resize(controlpoints);
	int count = 0;
	for (int p = 0; p < globalIndex.size(); ++p) {
		const vector<CPoint>& coefs = geometry[p].coefs();
		coefsNum += coefs.size();
		for (size_t m = 0; m < coefs.size(); ++m)
		{
			gl_coefs_index[globalIndex[p][m]].push_back(count);   // 全局 到 控制点序号 的索引
			count++;
		}
	}
	lambda = E * pr / ((1. + pr) * (1. - 2. * pr));    //拉梅常数
	mu = E / (2. * (1. + pr));
	sparseValue.resize(controlpoints * dim);      
	matrix.resize(controlpoints * dim, controlpoints * dim);
	rhs.resize(controlpoints * dim);
	setDirichlet();
}

//void stiffAssembler::computerIndex()
//{
//	std::vector<interFace> inf;
//	computerTopogy(inf);
//	std::vector<unordered_set<int>> interface_index(geometry.size());            // 存储每片索引点
//	std::vector<std::set<pair<int, int>>> connect_vertex(inf.size());  // 存储接口的重合索引点
//	// 提取连接边
//	for (int i = 0; i < inf.size(); ++i) { // 遍历所有接口
//		auto& it = inf[i];
//		std::vector<int>& inf1 = it.block1;
//		std::vector<int>& inf2 = it.block2;
//		if (inf1.size() != inf2.size()) {
//			cout << "接口错误" << "\n";
//			return;
//		}
//		for (int ii = 0, jj = 0; ii < inf1.size(); ++ii, ++jj) {
//			interface_index[it.p1].insert(inf1[ii]);
//			interface_index[it.p2].insert(inf2[jj]);
//			connect_vertex[i].insert({ inf1[ii],inf2[jj] });
//		}
//	}
//
//	// 构造全局和局部索引关系
//	for (int i = 0; i < globalIndex.size(); ++i) {
//		int m = geometry[i].getSizeUControlPoints(), n = geometry[i].getSizeVControlPoints(),
//			r = geometry[i].getSizeWControlPoints();
//		patchInfo.push_back(m * n * r);
//		globalIndex[i].resize(m * n * r,-1);
//		coefsNum += m * n * r;
//		for (int j = 0; j < globalIndex[i].size(); ++j) {
//			if (!interface_index[i].count(j))
//				globalIndex[i][j] = controlpoints++;
//		}
//	}
//
//	// 加入接口处的索引
//
//	std::map<LocInd, int> connect_point;
//	int ttt = 0;
//	for (int i = 0; i < inf.size(); ++i) {
//		auto& it = inf[i];
//		std::vector<int>& inf1 = it.block1;
//		std::vector<int>& inf2 = it.block2;
//		if (inf1.size() != inf2.size()) {
//			cout << "接口错误" << "\n";
//			return;
//		}
//		for (int ii = 0, jj = 0; ii < inf1.size(); ++ii, ++jj) {
//			LocInd ppp;
//			ppp.x = it.loc[ii].getX();
//			ppp.y = it.loc[ii].getY();
//			ppp.z = it.loc[ii].getZ();
//			
//			if (connect_point.count(ppp)) {
//				globalIndex[it.p2][inf2[jj]] = globalIndex[it.p1][inf1[ii]] = connect_point[ppp];
//			}
//			else {
//				connect_point[ppp] = controlpoints;
//				globalIndex[it.p2][inf2[jj]] = globalIndex[it.p1][inf1[ii]] = controlpoints++;
//				ttt++;
//			}
//		}
//	}
//}

void stiffAssembler::computerIndex() {
	
	std::vector<LocInd> indexSelect;
	for (size_t i = 0; i < geometry.size();++i)
	{
		auto coefs = geometry[i].coefs();
		globalIndex[i].resize((int)coefs.size());
		int index = 0;
		for (CPoint& it : coefs) {
			LocInd temp;
			temp.p = i; temp.x = it.getX();
			temp.y = it.getY(); temp.z = it.getZ();
			temp.index = index++;
			indexSelect.emplace_back(temp);
		}
	}
	std::sort(indexSelect.begin(), indexSelect.end());
	int p = indexSelect[0].p, loc = indexSelect[0].index;
	globalIndex[p][loc] = controlpoints++;
	for (size_t i = 1; i < indexSelect.size();++i) {
		p = indexSelect[i].p; loc = indexSelect[i].index;
		if (indexSelect[i] == indexSelect[i - 1])
			globalIndex[p][loc] = controlpoints - 1;
		else
			globalIndex[p][loc] = controlpoints++;
	}
}

int stiffAssembler::numDofs()
{
	return controlpoints * dim;
}

//bool stiffAssembler::computerTopogy(std::vector<interFace>& inf)
//{
//	std::vector<CFace> blockFace;
//
//	for (int i = 0; i < geometry.size(); ++i) {
//		findBoundaryLoc(i, blockFace);
//	}
//	std::sort(blockFace.begin(), blockFace.end());
//
//	for (int i = 0; i < blockFace.size()-1; ++i)
//	{
//		if (blockFace[i] == blockFace[i + 1]) {
//			interFace temp;
//			temp.p1 = blockFace[i].p; temp.dir1 = blockFace[i].dir;
//			temp.p2 = blockFace[i + 1].p; temp.dir2 = blockFace[i + 1].dir;
//			auto& loc1 = blockFace[i].locInfo, & loc2 = blockFace[i + 1].locInfo;
//			for (int t = 0; t < loc1.size(); ++t) {
//				temp.block1.emplace_back(loc1[t].index);
//				temp.block2.emplace_back(loc2[t].index);
//				temp.loc.emplace_back(CPoint(loc1[t].x, loc1[t].y, loc1[t].z));
//			}
//			inf.emplace_back(temp);
//		}
//	}
//	return true;
//}

//void stiffAssembler::findBoundaryLoc(int blockIndex, std::vector<CFace>& blockFace)
//{
//	std::vector<int> res;
//	auto& geo = geometry[blockIndex];
//	int um = geo.getSizeUControlPoints(), vm = geo.getSizeVControlPoints(),
//		wm = geo.getSizeWControlPoints();
//	auto points = geo.getControlPoints();
//	int allCon = um * vm * wm;
//	//east
//	std::vector<CFace> face(6);
//	for (int i = 0; i < 6; ++i) {
//		face[i].p = blockIndex;
//		face[i].dir = i+1;
//	}
//	int count = 0;
//	for (int k = 0; k < wm; ++k)
//	{
//		for (int j = 0; j < vm; ++j)
//		{
//			for (int i = 0; i < um; ++i)
//			{
//				LocInd temp;
//				temp.x = points[i][j][k].getX();
//				temp.y = points[i][j][k].getY();
//				temp.z = points[i][j][k].getZ();
//				temp.index = count++;
//				if (i == 0) face[0].locInfo.push_back(temp);
//				if (i == um - 1) face[1].locInfo.push_back(temp);
//				if (j == 0) face[2].locInfo.push_back(temp);
//				if (j == vm - 1) face[3].locInfo.push_back(temp);
//				if (k == 0) face[4].locInfo.push_back(temp);
//				if (k == wm - 1) face[5].locInfo.push_back(temp);
//			}
//		}
//	}
//	for (int i = 0; i < 6; ++i) {
//		std::sort(face[i].locInfo.begin(), face[i].locInfo.end());
//		blockFace.emplace_back(face[i]);
//	}
//	
//}

void stiffAssembler::stiffMatrix()
{
	std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
	std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
	int N_D;
	//遍历所有面片
	Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dim, dim), C;
	ToolFunc::setC(C, identity, lambda, mu);
	for (int i = 0; i < geometry.size(); ++i) {
		// 某一片模型
		CBSplineSolid* geo = &geometry[i];
		std::vector<double> uknots = geo->getUKnots(),
			vknots = geo->getVKnots(), wknots = geo->getWKnots();
		int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
			gaussNum_w = geo->getWDegree() + 1;
		int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
		p_u.resize(gaussNum_u);
		p_v.resize(gaussNum_v);
		p_v.resize(gaussNum_w);
		std::vector<Eigen::MatrixXd> basisValue(2);
		std::vector<int> localIndex(gaussNum);
		for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // 单元u方向循环
		{
			ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
			for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // 单元v方向循环
			{
				//高斯点选择
				ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
				for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval)
				{
					ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
					//局部索引号
					geo->active(u_interval, v_interval, w_interval, localIndex);
					N_D = localIndex.size();
					//基函数
					geo->eval_AllDers(p_u, p_v, p_w, basisValue);
					Eigen::MatrixXd localMats;
					localMats.setZero(dim * N_D, dim * N_D);
					for (int kk = 0; kk < gaussNum_w; ++kk) {
						for (int jj = 0; jj < gaussNum_v; ++jj) {
							for (int ii = 0; ii < gaussNum_u; ++ii) {
								int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
								Eigen::MatrixXd jac, P;
								geo->getJacobi(p_u[ii], p_v[jj], p_w[kk], localIndex, jac);
								double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jac.determinant());
								Eigen::MatrixXd physGrad, B_i, B_j, tempK, K;
								ToolFunc::transformGradients(jac, q, basisValue[1], physGrad, dim);
								// 局部矩阵计算 装配
								for (int j = 0; j < N_D; j++) {
									ToolFunc::setB(B_i, identity, physGrad.col(j));
									tempK = B_i.transpose() * C;
									for (int k = 0; k < N_D; k++) {
										ToolFunc::setB(B_j, identity, physGrad.col(k));
										K = tempK * B_j;
										for (int di = 0; di < dim; di++) {
											for (int dj = 0; dj < dim; dj++) {
												localMats(di * N_D + j, dj * N_D + k) += weight * K(di, dj);
											}
										}
									}
								}

							}
						}
					}
					// 单元刚度矩阵组装总刚度矩阵
					for (int m = 0; m < N_D; ++m) {
						for (int n = 0; n < N_D; ++n) {
							for (int di = 0; di < dim; di++) {
								for (int dj = 0; dj < dim; dj++) {
									//if (fabs(localMats(N_D * di + m, N_D * dj + n)) < tol) continue;
									int r = globalIndex[i][localIndex[m]] + di * controlpoints,
										c = globalIndex[i][localIndex[n]] + dj * controlpoints;
									if (globalDirichlet[r] == 1 || globalDirichlet[c] == 1) continue;

									if (!sparseValue[r].count(c)) sparseValue[r][c] = 0;
									sparseValue[r][c] += localMats(N_D * di + m, N_D * dj + n);
								}
							}
						}
					}
				}
			}
		}
	}
	std::vector< Eigen::Triplet<double> > tripletList;
	for (int begin = 0; begin < sparseValue.size(); ++begin) {
		for (auto& col : sparseValue[begin]) {
			tripletList.emplace_back(Eigen::Triplet<double>( begin,col.first,col.second ));
		}
	}
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}



void stiffAssembler::stiffNonlinearMatrix(const Eigen::MatrixXd& displacement, bool isAssembleMatrix)
{
	std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
	std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
	int N_D,begin = 0;
	//遍历所有面片
	Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dim, dim), C;
	if(m_materLaw.first == material_law::saint_venant_kirchhoff && isAssembleMatrix)
		ToolFunc::setC(C, identity, lambda, mu);
	for (int i = 0; i < geometry.size(); ++i) {
		// 某一片模型
		CBSplineSolid* geo = &geometry[i];
		
		int num = geometry[i].getSizeControlPoints();
		Eigen::MatrixXd blockU = displacement.block(begin, 0, num, 3);
		begin += num;

		std::vector<double> uknots = geo->getUKnots(),
			vknots = geo->getVKnots(), wknots = geo->getWKnots();
		int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
			gaussNum_w = geo->getWDegree() + 1;
		int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
		p_u.resize(gaussNum_u);
		p_v.resize(gaussNum_v);
		p_v.resize(gaussNum_w);
		std::vector<Eigen::MatrixXd> basisValue(2);
		std::vector<int> localIndex(gaussNum);
		Eigen::MatrixXd forceValues;
		for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // 单元u方向循环
		{
			ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
			for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // 单元v方向循环
			{
				//高斯点选择
				ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
				for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval)
				{
					ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
					//局部索引号
					geo->active(u_interval, v_interval, w_interval, localIndex);
					N_D = localIndex.size();
					//基函数
					geo->eval_AllDers(p_u, p_v, p_w, basisValue);

					//TODO
					bodyforce.eval(p_u, p_v, p_w, forceValues);

					Eigen::MatrixXd localMats, localRhs, jacGeo, jacDis, F, RCG, E, S, materialTangent;
					Eigen::VectorXd Svec, geometricTangentTemp, localResidual;
					Eigen::MatrixXd physGrad, phyDispJac, materialTangentTemp,B_i ,B_j;
					if(isAssembleMatrix)
						localMats.setZero(dim * N_D, dim * N_D);
					localRhs.setZero(dim * N_D, 1);
					for (int kk = 0; kk < gaussNum_w; ++kk) {
						for (int jj = 0; jj < gaussNum_v; ++jj) {
							for (int ii = 0; ii < gaussNum_u; ++ii) {
								int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
								
								//geo->getJacobi(p_u[ii], p_v[jj], p_w[kk], localIndex, jacGeo);  // 几何雅可比X/ξ  位移雅可比 U/ξ
								geo->getJacobiU(p_u[ii], p_v[jj], p_w[kk], localIndex, blockU, jacGeo, jacDis);
								double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jacGeo.determinant()); //权值

								ToolFunc::transformGradients(jacGeo, q, basisValue[1], physGrad, dim);
								phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());


								F = identity + phyDispJac;   // F = I + du/dx

								double J = F.determinant();

								RCG = F.transpose() * F;
								E = 0.5 * (RCG - identity);

								if (m_materLaw.first == material_law::saint_venant_kirchhoff) {
									S = lambda * E.trace() * identity + 2 * mu * E;       // 应力
									ToolFunc::voigtStressToK(S, Svec);                              // 应力张量
								}
								else {    // 其它材料模型

								}

								// 局部矩阵计算 装配
								for (size_t j = 0; j < N_D; j++) {
									 //B = BL + BN 即小位移矩阵+大位移矩阵
									ToolFunc::setB(B_i, F, physGrad.col(j));
									if (isAssembleMatrix) {
										materialTangentTemp = B_i.transpose() * C;
										geometricTangentTemp = S * physGrad.col(j);
										//loop over active basis functions (v_j)
										for (size_t k = 0; k < N_D; ++k)
										{
											ToolFunc::setB(B_j, F, physGrad.col(k));
											materialTangent.setZero(materialTangentTemp.rows(), B_j.cols());
											for (int mm = 0; mm < materialTangentTemp.rows(); mm++) {
												for (int aa = 0; aa < B_j.rows(); aa++) {
													for (int nn = 0; nn < B_j.cols(); nn++) {
														materialTangent(mm, nn) += materialTangentTemp(mm, aa) * B_j(aa, nn);
													}
												}
											} // 切线刚度矩阵第一项 位移项

											double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);
											// K_tg = K_tg_mat + I*K_tg_geo;
											for (short d = 0; d < dim; ++d)
												materialTangent(d, d) += geometricTangent;

											for (short di = 0; di < dim; ++di)
												for (short dj = 0; dj < dim; ++dj)
													localMats(di * N_D + j, dj * N_D + k) += weight * materialTangent(di, dj);
										}
									}

									localResidual = B_i.transpose() * Svec;
									for (short d = 0; d < dim; d++)
										localRhs(d * N_D + j) -= weight * localResidual(d);
								}
								for (short d = 0; d < dim; ++d)
									localRhs.middleRows(d * N_D, N_D).noalias() += weight * forceScaling * forceValues(d, q) * basisValue[0].col(q);

							}
						}
					}
					// 右端项
					for (int m = 0; m < N_D; ++m) {
						for (int di = 0; di < dim; di++) {
							int r = globalIndex[i][localIndex[m]] + di * controlpoints;
							if (globalDirichlet[r] == 1) continue;
							rhs(globalIndex[i][localIndex[m]] + di * controlpoints) += localRhs(N_D * di + m);
						}
					}


					// 单元刚度矩阵组装总刚度矩阵
					if (isAssembleMatrix) {
						for (int m = 0; m < N_D; ++m) {
							for (int n = 0; n < N_D; ++n) {
								for (int di = 0; di < dim; di++) {
									for (int dj = 0; dj < dim; dj++) {
										//if (fabs(localMats(N_D * di + m, N_D * dj + n)) < tol) continue;
										int r = globalIndex[i][localIndex[m]] + di * controlpoints,
											c = globalIndex[i][localIndex[n]] + dj * controlpoints;
										if (globalDirichlet[r] == 1 || globalDirichlet[c] == 1) continue;

										if (!sparseValue[r].count(c)) sparseValue[r][c] = 0;
										sparseValue[r][c] += localMats(N_D * di + m, N_D * dj + n);
									}
								}
							}
						}
					}


				}
			}
		}
	}
	std::vector< Eigen::Triplet<double> > tripletList;
	for (int begin = 0; begin < sparseValue.size(); ++begin) {
		for (auto& col : sparseValue[begin]) {
			tripletList.emplace_back(Eigen::Triplet<double>(begin, col.first, col.second));
		}
	}
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

void stiffAssembler::rhsVector()
{
	std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
	std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
	int N_D;
	auto it = neumanSides.begin();
	for (; it != neumanSides.end(); ++it) {
		int patchIndex = it->patchSide();
		boxSide side = it->side();
		CBSplineSolid* geo = &geometry[patchIndex];
		std::vector<double> uknots = geo->getUMKnots(),
			vknots = geo->getVMKnots(), wknots = geo->getWMKnots();
		int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
			gaussNum_w = geo->getWDegree() + 1;

		int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
		p_u.resize(gaussNum_u);
		p_v.resize(gaussNum_v);
		p_w.resize(gaussNum_w);
		//确定单元范围
		int u_interval_min = 1, u_interval_max = uknots.size() - 2;
		int v_interval_min = 1, v_interval_max = vknots.size() - 2;
		int w_interval_min = 1, w_interval_max = wknots.size() - 2;
		int u_add = 0, v_add = 0, w_add = 0;
		if (side.index() == boundary::west) {
			u_interval_min = 0;
			u_interval_max = 1;
			u_add = 1;
		}
		else if (side.index() == boundary::east) {
			u_interval_min = u_interval_max;
			u_interval_max += 1;
			u_add = -1;
		}
		else if (side.index() == boundary::south) {
			v_interval_min = 0;
			v_interval_max = 1;
			v_add = 1;
		}
		else if(side.index() == boundary::north) {
			v_interval_min = v_interval_max;
			v_interval_max += 1;
			v_add = -1;
		}
		else if (side.index() == boundary::front) {
			w_interval_min = 0;
			w_interval_max = 1;
			w_add = 1;
		}
		else {
			w_interval_min = w_interval_max;
			w_interval_max += 1;
			w_add = -1;
		}

		std::vector<int> localIndex(gaussNum);

		for (int w_interval = w_interval_min; w_interval < w_interval_max; ++w_interval) {
			ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
			for (int v_interval = v_interval_min; v_interval < v_interval_max; ++v_interval){
				ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
				for (int u_interval = u_interval_min; u_interval < u_interval_max; ++u_interval) {
					ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
					//高斯点选择
					geo->active(u_interval - 1 + u_add, v_interval - 1 + v_add,w_interval-1+w_add, localIndex);
					N_D = localIndex.size();
					Eigen::MatrixXd neumanValues, basisValue, localRhs;
					it->eval(p_u, p_v,p_w, neumanValues);
					geo->eval_into(p_u, p_v,p_w, basisValue);
					localRhs.setZero(dim * N_D, 1);
					for (int kk = 0; kk < gaussNum_w; ++kk) {
						for (int jj = 0; jj < gaussNum_v; ++jj) {
							for (int ii = 0; ii < gaussNum_u; ++ii) {
								int q = kk * gaussNum_v * gaussNum_u + jj * gaussNum_u + ii;
								Eigen::MatrixXd jac;
								geo->getJacobi(p_u[ii], p_v[jj], p_w[kk], localIndex, jac);
								Eigen::VectorXd unormal;
								ToolFunc::outerNormal(jac, side, unormal);
								double weight = w_u[ii] * w_v[jj] *w_w[kk] * unormal.norm();
								for (short d = 0; d < dim; ++d)
									localRhs.middleRows(d * N_D, N_D).noalias() += weight * neumanValues(d, q) * basisValue.col(q);
							}
						}
					}
					//cout << localRhs << "\n\n";
					// 单元右端项组装总右端项
					for (int m = 0; m < N_D; ++m) {
						for (int di = 0; di < dim; di++) {
							int r = globalIndex[patchIndex][localIndex[m]] + di * controlpoints;
							if (globalDirichlet[r] == 1) continue;
							rhs(globalIndex[patchIndex][localIndex[m]] + di * controlpoints) += localRhs(N_D * di + m);
						}
					}

				}
			}
		}
	}
}

void stiffAssembler::setDirichlet()
{
	globalDirichlet.resize(dim*controlpoints,0);
	auto it = dirichletSides.begin();
	for (; it != dirichletSides.end(); ++it) {
		int patchIndex = it->patchSide();
		boxSide side = it->side();
		CBSplineSolid* geo = &geometry[patchIndex];
		int m = geo->getSizeUControlPoints(), n = geo->getSizeVControlPoints(),
			o = geo->getSizeWControlPoints();
		if (side.index() == boundary::east || side.index() == boundary::west) {
			int u = side.index() == boundary::east ? 1 : 0;
			for (int k = 0; k < o; ++k) {
				for (int j = 0; j < n; ++j) {
					int localIndex = k*m*n + j * m + u * (m - 1);
					dir_index[globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = it->getValue();
					globalDirichlet[globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = 1;
				}
			}
		}
		else if(side.index() == boundary::south || side.index() == boundary::north) {
			int v = side.index() == boundary::north ? 1 : 0;
			for (int k = 0; k < o; ++k) {
				for (int i = 0; i < m; ++i) {
					int localIndex = k*m*n + i + v * ((n - 1) * m);
					dir_index[globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = it->getValue();
					globalDirichlet[globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = 1;
				}
			}
		}
		else {
			int w = side.index() == boundary::back ? 1 : 0;
			for (int j = 0; j < n; ++j) {
				for (int i = 0; i < m; ++i) {
					int localIndex = j * m + i + w * ((o - 1) * m * n);
					dir_index[globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = it->getValue();
					globalDirichlet[globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = 1;
				}
			}
		}
	}
}

void stiffAssembler::assemble()
{
	matrix.setZero();
	rhs.setZero();
	vector<unordered_map<int, double>>(controlpoints*dim).swap(sparseValue);

	stiffMatrix();
	rhsVector();

	for (auto& it : dir_index) {
		rhs(it.first) = it.second;
		matrix.coeffRef(it.first, it.first) = 1;
	}

	matrix.makeCompressed();

	if (matrix.rows() != rhs.rows()) {
		cout << "ERROR matrix size != rhs size" << endl;
		return;
	}
}

void stiffAssembler::assemble(const Eigen::VectorXd& solution, bool isAssemblerMatrix)
{
	Eigen::MatrixXd displacement;
	constructCoefs(solution, displacement);  // 构造初始解

	if (isAssemblerMatrix)
		matrix.setZero();
	rhs.setZero();
	vector<unordered_map<int, double>>(controlpoints * dim).swap(sparseValue);
	stiffNonlinearMatrix(displacement, isAssemblerMatrix);
	rhsVector();
	for (auto& it : dir_index) {
		rhs(it.first) = it.second;
		matrix.coeffRef(it.first, it.first) = 1;
	}
	auto mat = matrix.toDense();
	
	//std::ofstream myfile;
	//myfile.open("D:/Resource/FileDataInput/threeD/occ.csv");

	//for (int r = 0; r < mat.rows(); ++r) {
	//	for (int c = 0; c < mat.cols(); ++c)
	//		myfile << mat(r, c) << ",";
	//	myfile << "\n";
	//}
	//myfile.close();

	matrix.makeCompressed();
	if (matrix.rows() != rhs.rows()) {
		cout << "ERROR matrix size != rhs size" << endl;
		return;
	}
}

void stiffAssembler::setMaterialLaw(const pair<int, vector<double>> material)
{
	m_materLaw = material;
}

void stiffAssembler::updateCoefs(Eigen::MatrixXd& coefs)
{
	int begin = 0;
	for (size_t i = 0; i < geometry.size(); ++i) {	
		const vector<CPoint>& coef = geometry[i].coefs();
		int num = geometry[i].getSizeControlPoints();
		geometry[i].updateCoefs(coefs.block(begin,0,num,3));
		begin += num;
	}
}

void stiffAssembler::constructCoefs(const Eigen::VectorXd& solution, Eigen::MatrixXd& displacement)
{
	if (solution.rows() != dim * controlpoints) {
		cout << "solution length is not Equal to the size of contropoints" << endl;
		return;
	}
	displacement.resize(coefsNum, dim);
	for (int i = 0; i < controlpoints; ++i) {
		for (auto index : gl_coefs_index[i]) {
			for (int di = 0; di < dim; ++di) {
				displacement(index, di) = solution(i + di * controlpoints);
			}
		}
	}
}

const Eigen::SparseMatrix<double>& stiffAssembler::getMatrix() const
{
	return matrix;
}

const Eigen::VectorXd& stiffAssembler::getRhs() const
{
	return rhs;
}
