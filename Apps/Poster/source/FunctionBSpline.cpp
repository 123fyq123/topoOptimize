#include"FunctionBSpline.h"
#include"FileHandler.h"
#include<omp.h>
#include<fstream>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>
#include "Timer.h"
#include<numeric>
#include <algorithm>

#ifndef PI
#define PI 3.1415926
#endif // !PI

void FunctionBSpline::BdyFunLine(std::vector< std::vector<double> >& m_coef)
{
	int m_nptnum;
	std::vector<double> m_scoef;
	std::vector<double> m_sscoef;
	for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		int m_uptnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		int m_vptnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		int m_wptnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		int sameuvw = m_wptnum;
		int m_ncptnum = m_uptnum * m_vptnum * m_wptnum;
		V3f m_pt;
		double m_value;
		m_scoef.clear();
		m_scoef.resize(m_ncptnum, 0);
		m_pt = m_solids.Local_BSplineSolid[i].vertices_[m_solids.Control_Point_id[i][0][0][0]];//(u,v,w)=(0,0,0)
		m_value = ExactSol(m_pt.x, m_pt.y, m_pt.z);
		m_scoef[m_solids.Control_Point_id[i][0][0][0]] = m_value;

		m_pt = m_solids.Local_BSplineSolid[i].vertices_[m_solids.Control_Point_id[i][0][0][sameuvw - 1]];//(u,v,w)=(1,0,0)
		m_value = ExactSol(m_pt.x, m_pt.y, m_pt.z);
		m_scoef[m_solids.Control_Point_id[i][0][0][sameuvw - 1]] = m_value;

		m_pt = m_solids.Local_BSplineSolid[i].vertices_[m_solids.Control_Point_id[i][0][sameuvw - 1][0]];//(u,v,w)=(0,1,0)
		m_value = ExactSol(m_pt.x, m_pt.y, m_pt.z);
		m_scoef[m_solids.Control_Point_id[i][0][sameuvw - 1][0]] = m_value;

		m_pt = m_solids.Local_BSplineSolid[i].vertices_[m_solids.Control_Point_id[i][0][sameuvw - 1][sameuvw - 1]];//(u,v,w)=(1,1,0)
		m_value = ExactSol(m_pt.x, m_pt.y, m_pt.z);
		m_scoef[m_solids.Control_Point_id[i][0][sameuvw - 1][sameuvw - 1]] = m_value;

		m_pt = m_solids.Local_BSplineSolid[i].vertices_[m_solids.Control_Point_id[i][sameuvw - 1][0][0]];//(u,v,w)=(0,0,1)
		m_value = ExactSol(m_pt.x, m_pt.y, m_pt.z);
		m_scoef[m_solids.Control_Point_id[i][sameuvw - 1][0][0]] = m_value;

		m_pt = m_solids.Local_BSplineSolid[i].vertices_[m_solids.Control_Point_id[i][sameuvw - 1][0][sameuvw - 1]];//(u,v,w)=(1,0,1)
		m_value = ExactSol(m_pt.x, m_pt.y, m_pt.z);
		m_scoef[m_solids.Control_Point_id[i][sameuvw - 1][0][sameuvw - 1]] = m_value;

		m_pt = m_solids.Local_BSplineSolid[i].vertices_[m_solids.Control_Point_id[i][sameuvw - 1][sameuvw - 1][0]];//(u,v,w)=(0,1,1)
		m_value = ExactSol(m_pt.x, m_pt.y, m_pt.z);
		m_scoef[m_solids.Control_Point_id[i][sameuvw - 1][sameuvw - 1][0]] = m_value;

		m_pt = m_solids.Local_BSplineSolid[i].vertices_[m_solids.Control_Point_id[i][sameuvw - 1][sameuvw - 1][sameuvw - 1]];//(u,v,w)=(1,1,1)
		m_value = ExactSol(m_pt.x, m_pt.y, m_pt.z);
		m_scoef[m_solids.Control_Point_id[i][sameuvw - 1][sameuvw - 1][sameuvw - 1]] = m_value;

		m_sscoef.clear();
		std::vector<int> m_para(3, 0); std::vector<int> m_bdy(3, 0);
		//近似第1条边界上的数值解
		m_bdy[0] = 1; m_para[1] = 0; m_para[2] = 0;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_uptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][k][0][0]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第2条边界上的数值解
		m_bdy[1] = 1; m_para[0] = 0; m_para[2] = 0;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][0][k][0]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第3条边界上的数值解
		m_bdy[0] = 1; m_para[1] = 1; m_para[2] = 0;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_uptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][k][sameuvw - 1][0]] = m_sscoef[k - 1];


		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第4条边界上的数值解
		m_bdy[1] = 1; m_para[0] = 1; m_para[2] = 0;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][sameuvw - 1][k][0]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第5条边界上的数值解
		m_bdy[2] = 1; m_para[0] = 1; m_para[1] = 0;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][sameuvw - 1][0][k]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第6条边界上的数值解
		m_bdy[2] = 1; m_para[0] = 0; m_para[1] = 0;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][0][0][k]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第7条边界上的数值解
		m_bdy[2] = 1; m_para[0] = 0; m_para[1] = 1;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][0][sameuvw - 1][k]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第8条边界上的数值解
		m_bdy[2] = 1; m_para[0] = 1; m_para[1] = 1;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][sameuvw - 1][sameuvw - 1][k]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第9条边界上的数值解
		m_bdy[0] = 1; m_para[1] = 0; m_para[2] = 1;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][k][0][sameuvw - 1]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第10条边界上的数值解
		m_bdy[1] = 1; m_para[0] = 0; m_para[2] = 1;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][0][k][sameuvw - 1]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第11条边界上的数值解
		m_bdy[0] = 1; m_para[1] = 1; m_para[2] = 1;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][k][sameuvw - 1][sameuvw - 1]] = m_sscoef[k - 1];

		m_sscoef.clear(); m_para.resize(3, 0); m_bdy.resize(3, 0);
		//近似第12条边界上的数值解
		m_bdy[1] = 1; m_para[0] = 1; m_para[2] = 1;
		m_solids.SolveCoef(i, m_scoef, m_para, m_bdy, m_sscoef);
		for (int k = 1; k < m_vptnum - 1; k++)
			m_scoef[m_solids.Control_Point_id[i][sameuvw - 1][k][sameuvw - 1]] = m_sscoef[k - 1];

		m_coef[i] = m_scoef;

		/*int num = 0;
		for (int aa = 0; aa < m_scoef.size();aa++)
		{
			if (m_scoef[aa] != 0)
				num++;
		}*/
	}
	////内边界的基对应的值赋为0
	//for (int i = 0; i<m_faces.size(); i++)
	//{
	//	int m_uCpt = m_faces[i].m_U.size() - m_faces[i].m_udegree - 1;
	//	int m_vCpt = m_faces[i].m_V.size() - m_faces[i].m_vdegree - 1;
	//	if (m_faces[i].m_edgeNbour[0].m_edgeid != -1)
	//	{
	//		for (int j = 1; j<m_vCpt - 1; j++)
	//			m_coef[i][j] = 0;
	//		if (m_faces[i].m_edgeNbour[1].m_edgeid != -1)
	//			m_coef[i][0] = 0;
	//		if (m_faces[i].m_edgeNbour[3].m_edgeid != -1)
	//			m_coef[i][m_vCpt - 1] = 0;
	//	}
	//	if (m_faces[i].m_edgeNbour[1].m_edgeid != -1)
	//	{
	//		for (int j = 1; j<m_uCpt - 1; j++)
	//			m_coef[i][j*m_vCpt] = 0;
	//		if (m_faces[i].m_edgeNbour[2].m_edgeid != -1)
	//			m_coef[i][(m_uCpt - 1)*m_vCpt] = 0;
	//		if (m_faces[i].m_edgeNbour[0].m_edgeid != -1)
	//			m_coef[i][0] = 0;
	//	}
	//	if (m_faces[i].m_edgeNbour[2].m_edgeid != -1)
	//	{
	//		for (int j = 1; j<m_vCpt - 1; j++)
	//			m_coef[i][(m_uCpt - 1)*m_vCpt + j] = 0;
	//		if (m_faces[i].m_edgeNbour[1].m_edgeid != -1)
	//			m_coef[i][(m_uCpt - 1)*m_vCpt] = 0;
	//		if (m_faces[i].m_edgeNbour[3].m_edgeid != -1)
	//			m_coef[i][m_uCpt*m_vCpt - 1] = 0;
	//	}
	//	if (m_faces[i].m_edgeNbour[3].m_edgeid != -1)
	//	{
	//		for (int j = 1; j<m_uCpt - 1; j++)
	//			m_coef[i][j*m_vCpt + m_vCpt - 1] = 0;
	//		if (m_faces[i].m_edgeNbour[2].m_edgeid != -1)
	//			m_coef[i][m_uCpt*m_vCpt - 1] = 0;
	//		if (m_faces[i].m_edgeNbour[0].m_edgeid != -1)
	//			m_coef[i][m_vCpt - 1] = 0;
	//	}
	//}
}

void FunctionBSpline::Global_MatrixAndRhs(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs)
{
	int m_GlobalSize = m_solids.Global_BSplineSolid.vertices_.size();
	/*for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		m_GlobalSize += m_solids.Local_BSplineSolid[i].vertices_.size();
	}*/
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.resize(m_GlobalSize, 0);

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//m_u[0] = 0; m_u[m_GaussNum - 1] = 1;
			//m_u[1] = 0.3;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//m_v[0] = 0; m_v[m_GaussNum - 1] = 1;
				//m_v[1] = 0.3;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//m_w[0] = 0; m_w[m_GaussNum - 1] = 1;
					//m_w[1] = 0.3;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								//m_BaseIndex = ww * m_vPtnum*m_uPtnum + j * m_uPtnum + ii;
								m_BaseIndex = m_solids.Control_Point_id[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								m_bigrhs[m_solids.CP_Local_to_Global[i][m_BaseIndex]] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								/*int check1 = m_solids.CP_Local_to_Global[i][m_baseindex[ii]];
								int check2 = m_solids.CP_Local_to_Global[i][m_baseindex[j]];*/
								m_bigmatrix[m_solids.CP_Local_to_Global[i][m_baseindex[ii]]][m_solids.CP_Local_to_Global[i][m_baseindex[j]]] += (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
							}
						}
					}

				}
			}
		}
	}
}

void FunctionBSpline::Global_MatrixAndRhs(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix
	, std::vector< std::vector<double>>& local_bigrhs)
{
	int m_GlobalSize = m_solids.Global_BSplineSolid.vertices_.size();
	local_matrix.resize(m_solids.Local_BSplineSolid.size());
	local_bigrhs.resize(m_solids.Local_BSplineSolid.size());
	for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		int size = m_solids.Local_BSplineSolid[i].vertices_.size();
		local_matrix[i].resize(size, std::vector<double>(size, 0));
		local_bigrhs[i].resize(size, 0);
	}
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.resize(m_GlobalSize, 0);

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								//m_BaseIndex = ww * m_vPtnum*m_uPtnum + j * m_uPtnum + ii;
								m_BaseIndex = m_solids.Control_Point_id[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								local_bigrhs[i][m_BaseIndex] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								m_bigrhs[m_solids.CP_Local_to_Global[i][m_BaseIndex]] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8 * (umin + (1 - umin) * pow(material_distribution[i], p));
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double cal_temp = (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii]][m_baseindex[j]] += cal_temp;
								m_bigmatrix[m_solids.CP_Local_to_Global[i][m_baseindex[ii]]][m_solids.CP_Local_to_Global[i][m_baseindex[j]]] +=
									cal_temp * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
							/*m_bigmatrix[m_solids.CP_Local_to_Global[i][m_baseindex[ii]]][m_solids.CP_Local_to_Global[i][m_baseindex[j]]] =
								local_matrix[i][m_baseindex[ii]][m_baseindex[j]] * (umin + (1 - umin)*pow(material_distribution[i], p));*/
						}
					}
				}
			}
		}
	}
}

void FunctionBSpline::Global_MatrixAndRhs1(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix
	, std::vector< std::vector<double>>& local_bigrhs, std::map<int, int>& base_check)
{
	int m_GlobalSize = 0;
	int global_ver_num = m_solids.Global_BSplineSolid.vertices_.size();
	for (int i = 0; i < global_ver_num; i++)
	{
		if (!m_solids.Global_BSplineSolid.is_bdy(VertexHandle(i)))
		{
			m_GlobalSize++;
		}
	}
	local_matrix.resize(m_solids.Local_BSplineSolid.size());
	local_bigrhs.resize(m_solids.Local_BSplineSolid.size());
	for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		int size = m_solids.Local_BSplineSolid[i].vertices_.size();
		local_matrix[i].resize(size, std::vector<double>(size, 0));
		local_bigrhs[i].resize(size, 0);
	}
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.resize(m_GlobalSize, 0);


	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	int base_num = 0;
	for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_BSplineSolid.size() << std::endl;
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵


		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								bool base_bdy_check = false;
								m_BaseIndex = m_solids.Control_Point_id[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								if (base_check.find(m_solids.CP_Local_to_Global[i][m_BaseIndex]) == base_check.end()
									&& !m_solids.Global_BSplineSolid.is_bdy(m_solids.CP_Local_to_Global[i][m_BaseIndex]))
								{
									base_check.insert(std::pair<int, int>(m_solids.CP_Local_to_Global[i][m_BaseIndex], base_num));
									base_num++;
									base_bdy_check = true;
								}
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								local_bigrhs[i][m_BaseIndex] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								if (base_bdy_check == true)
									m_bigrhs[base_check[m_solids.CP_Local_to_Global[i][m_BaseIndex]]] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double cal_temp = (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii]][m_baseindex[j]] += cal_temp;
								if (base_check.find(m_solids.CP_Local_to_Global[i][m_baseindex[ii]]) != base_check.end()
									&& base_check.find(m_solids.CP_Local_to_Global[i][m_baseindex[j]]) != base_check.end())
									m_bigmatrix[base_check[m_solids.CP_Local_to_Global[i][m_baseindex[ii]]]][base_check[m_solids.CP_Local_to_Global[i][m_baseindex[j]]]] +=
									cal_temp * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
						}
					}

				}
			}
		}
	}
}

void FunctionBSpline::GaussIntegalPoint(int n, std::vector<double>& w, std::vector<double>& gausspoint)
{
	int iter, m, el;
	int mm;
	double nn;
	for (int i = 0; i < n; i++)
	{
		w.push_back(0);	gausspoint.push_back(0);
	}
	iter = 2;
	m = (n + 1) / 2; el = n * (n + 1);
	mm = 4 * m - 1;
	nn = (1 - (1 - 1.0 / n) / (8 * n * n));
	std::vector<double> t, x0;
	double tt;
	for (int i = 3; i <= mm; i += 4)
	{
		tt = (PI / (4 * n + 2)) * i;
		t.push_back(tt);
		x0.push_back(nn * cos(tt));
	}
	std::vector<double> pkm1, pkp1, t1, pk;
	int x0_size = x0.size();
	std::vector<double> den, d1, dpn, d2pn, d3pn, d4pn, u, v, h, p, dp, fx;
	t1.resize(x0_size);  pkm1.resize(x0_size);
	pk.resize(x0_size);  t1.resize(x0_size);
	pkp1.resize(x0_size); den.resize(x0_size);
	dpn.resize(x0_size); d2pn.resize(x0_size);
	d3pn.resize(x0_size); d4pn.resize(x0_size);
	d1.resize(x0_size);  u.resize(x0_size);
	v.resize(x0_size);   h.resize(x0_size);
	p.resize(x0_size);   fx.resize(x0_size);
	dp.resize(x0_size);

	for (int i = 0; i < iter; i++)
	{
		for (int j = 0; j < x0_size; j++)
		{
			pkm1[j] = 1;
			pk[j] = x0[j];
		}
		for (int k = 2; k <= n; k++)
			for (int r = 0; r < x0_size; r++)
			{
				t1[r] = x0[r] * pk[r];
				pkp1[r] = t1[r] - pkm1[r] - 1.0 * (t1[r] - pkm1[r]) / k + t1[r];
				pkm1[r] = pk[r]; pk[r] = pkp1[r];
			}
		for (int j = 0; j < x0_size; j++)
		{
			den[j] = 1.0 - x0[j] * x0[j];
			d1[j] = n * (pkm1[j] - x0[j] * pk[j]);
			dpn[j] = d1[j] / den[j];
			d2pn[j] = (2.0 * x0[j] * dpn[j] - el * pk[j]) / den[j];
			d3pn[j] = (4 * x0[j] * d2pn[j] + (2 - el) * dpn[j]) / den[j];
			d4pn[j] = (6 * x0[j] * d3pn[j] + (6 - el) * d2pn[j]) / den[j];
			u[j] = pk[j] / dpn[j]; v[j] = d2pn[j] / dpn[j];
			h[j] = -u[j] * (1 + 0.5 * u[j] * (v[j] + u[j] * (v[j] * v[j] - u[j] * d3pn[j] / (3 * dpn[j]))));
			p[j] = pk[j] + h[j] * (dpn[j] + 0.5 * h[j] * (d2pn[j] + 1.0 * h[j] / 3 * (d3pn[j] + 0.25 * h[j] * d4pn[j])));
			dp[j] = dpn[j] + h[j] * (d2pn[j] + 0.5 * h[j] * (d3pn[j] + h[j] * d4pn[j] / 3));
			h[j] = h[j] - p[j] / dp[j]; x0[j] = x0[j] + h[j];
		}
	}

	for (int j = 0; j < x0_size; j++)
	{
		gausspoint[j] = -x0[j] - h[j];
		fx[j] = d1[j] - h[j] * el * (pk[j] + h[j] / 2 * (dpn[j] + h[j] / 3 * (d2pn[j] + h[j] / 4 * (d3pn[j] + 0.2 * h[j] * d4pn[j]))));
		w[j] = 2 * (1.0 - gausspoint[j] * gausspoint[j]) / (fx[j] * fx[j]);
	}
	if (m + m > n) gausspoint[m - 1] = 0;
	if ((m + m) != n) m = m - 1;
	std::vector<int> jj, n1j;
	jj.resize(m); n1j.resize(m);
	for (int j = 0; j < m; j++)
	{
		jj[j] = j; n1j[j] = n - jj[j] - 1;
		gausspoint[n1j[j]] = -gausspoint[jj[j]];
		w[n1j[j]] = w[jj[j]];
	}
}

void FunctionBSpline::SolveLinearEquation(int mtype, const std::vector<std::vector<double>>& coeff_matrix,
	const std::vector<double>& load_vect, std::vector<double>& m_solute)
{
	int matrix_size = coeff_matrix.size();

	int* ia = new int[matrix_size + 1];   //next row start index
	std::vector<int> jav;
	std::vector<double> val;

	ia[0] = 1;
	if (11 == mtype)//real unsymmetric
	{
		for (int i = 0; i < matrix_size; ++i)
		{
			for (int j = 0; j < matrix_size; ++j)
			{
				if (fabs(coeff_matrix[i][j]) > 10e-8)
				{
					val.push_back(coeff_matrix[i][j]);
					jav.push_back(j + 1);

				}
			}
			ia[i + 1] = jav.size() + 1;
		}
	}
	if (-2 == mtype)
	{
		for (int i = 0; i < matrix_size; ++i)
		{
			for (int j = i; j < matrix_size; ++j)//for (int j=0;j < matrix_size-i;++j)
			{
				if (fabs(coeff_matrix[i][j]) > 10e-8)
				{
					val.push_back(coeff_matrix[i][j]);
					jav.push_back(j + 1);//jav.push_back(j+1+i);

				}
			}
			ia[i + 1] = jav.size() + 1;
		}
		std::cout << "symmetric matrix!!!:";
	}


	int n = val.size();      //non zero elements

	int* ja = new int[n];
	double* a = new double[n];
	double* b = new double[matrix_size];
	double* x = new double[matrix_size];   /* RHS and solution vectors. */

	for (int i = 0; i < matrix_size; ++i)
		b[i] = load_vect[i];

	for (int i = 0; i < n; ++i)
	{
		ja[i] = jav[i];
		a[i] = val[i];
	}

	int nrhs = 1;                       /* Number of right hand sides. */
	//int mtype = 11;                     /* Real unsymmetric matrix */						

	void* pt[64];                           /* Internal solver memory pointer pt, */
	int iparm[64];                          /* Pardiso control parameters. */
	int maxfct, mnum, phase, error, msglvl;
	int i;
	double ddum; /* Double dummy *//* Auxiliary variables. */
	int idum; /* Integer dummy. */

	for (i = 0; i < 64; i++)
		iparm[i] = 0;

	iparm[0] = 1;                              /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
					  /* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;//omp_get_max_threads();
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 16; /* Default logical fortran unit number for output */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Not in use */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 0; /* Don't print statistical information in file */
	error = 0; /* Initialize error flag */
	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
		pt[i] = 0;

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;                           //only reordering and symbolic factorization
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&matrix_size, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&matrix_size, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;    // only factorization
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */

	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&matrix_size, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, b, x, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	m_solute.resize(matrix_size);
	for (i = 0; i < matrix_size; i++)
	{
		m_solute[i] = x[i];
	}

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1; /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&matrix_size, &ddum, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
}

void FunctionBSpline::SolveLinearEquation(int mtype, int GlobalSize, std::vector<int>& rows, std::vector<int>& cols,
	std::vector<double>& elems, std::vector<double>& load_vect, std::vector<double>& m_solute)
{
	/*rows = { 0,2,3,4 };
	cols = { 0,1,1,2 };
	elems = { 1,2,1,1 };
	load_vect = { 5.0,4.0,3.0 };
	int matrix_size = 3;*/
	int matrix_size = GlobalSize;
	int* ia = new int[matrix_size + 1];   //next row start index
	std::vector<int> jav;
	std::vector<double> val;

	//ia[0] = 1;
	for (int k = 0; k <= matrix_size; k++)
		ia[k] = rows[k];

	int n = elems.size();      //non zero elements

	int* ja = new int[n];
	double* a = new double[n];
	double* b = new double[matrix_size];
	double* x = new double[matrix_size];   /* RHS and solution vectors. */

	for (int i = 0; i < matrix_size; ++i)
		b[i] = load_vect[i];

	for (int i = 0; i < n; ++i)
	{
		/* ja[i] = jav[i];
		 a[i] = val[i];*/
		ja[i] = cols[i];
		a[i] = elems[i];
	}

	int nrhs = 1;                       /* Number of right hand sides. */

	void* pt[64];                           /* Internal solver memory pointer pt, */
	int iparm[64];                          /* Pardiso control parameters. */
	int maxfct, mnum, phase, error, msglvl;
	//	int i;
	double ddum; /* double dummy *//* Auxiliary variables. */
	int idum; /* Integer dummy. */
	int i;

	for (i = 0; i < 64; i++)
		iparm[i] = 0;

	iparm[0] = 1;                              /* No solver default */
	iparm[1] = 2;
	//iparm[1] = 2; /* Fill-in reordering from METIS */
					  /* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1; //omp_get_max_threads();
	// 2023.5.24  change iparm[3] from 0 to 62
	iparm[3] = 62; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 16; /* Default logical fortran unit number for output */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Not in use */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */
	iparm[34] = 1;
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 0; /* Don't print statistical information in file */
	error = 0; /* Initialize error flag */
	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
		pt[i] = 0;

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;                           //only reordering and symbolic factorization
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&matrix_size, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&matrix_size, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;    // only factorization
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */

	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&matrix_size, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, b, x, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	m_solute.resize(matrix_size);
	for (i = 0; i < matrix_size; i++)
	{
		m_solute[i] = x[i];
	}
	//LogsolveErr(m_solute, rows, cols, elems, load_vect);
	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1; /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&matrix_size, &ddum, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	std::vector<int>().swap(rows);
	std::vector<int>().swap(cols);
	std::vector<double>().swap(elems);
	delete[] ia;
	delete[] ja;
	delete[] a;
	delete[] b;
	delete[] x;
}

void FunctionBSpline::Error(std::vector<double> m_coef, std::vector<double>& m_errors)
{

}

double FunctionBSpline::ExactSol(double x, double y, double z)
{
	return x * (1 - x) * y * (1 - y) * (1 + z) * (1 - z);
	//return x * x + y * y + z * z;
}

double FunctionBSpline::RHS(double x, double y, double z)
{
	//return 2 * x*y*(x - 1)*(y - 1) + 2 * x*(x - 1)*(z - 1)*(z + 1) + 2 * y*(y - 1)*(z - 1)*(z + 1);
	return 1;//f = 1加一个负号
	//return -8;

}

double FunctionBSpline::RHS1(double u, double v, double w)
{
	double x;
	x = sqrt((u - 0.3) * (u - 0.3) + (v - 0.3) * (v - 0.3) + (w - 0.3) * (w - 0.3));
	return (1 - x) * 1;
}

void FunctionBSpline::BdyFunFace(std::vector< std::vector<double> >& m_coef)
{
	std::vector< std::vector<double> > m_coeftest = m_coef;
	for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		//第一个面
		int uvw = U, uvw_id = 0;
		facePDE(i, uvw, uvw_id, m_coef, m_coeftest);
		//第二个面
		uvw = U, uvw_id = 1;
		facePDE(i, uvw, uvw_id, m_coef, m_coeftest);
		//第三个面
		uvw = V, uvw_id = 0;
		facePDE(i, uvw, uvw_id, m_coef, m_coeftest);
		//第四个面
		uvw = V, uvw_id = 1;
		facePDE(i, uvw, uvw_id, m_coef, m_coeftest);
		//第五个面
		uvw = W, uvw_id = 0;
		facePDE(i, uvw, uvw_id, m_coef, m_coeftest);
		//第六个面
		uvw = W, uvw_id = 1;
		facePDE(i, uvw, uvw_id, m_coef, m_coeftest);
	}
}

void FunctionBSpline::Global_MatrixAndRhsFace(int i, int uvw, int uvw_id, Matrix& m_bigmatrix, std::vector<double>& m_bigrhs)
{
	int m_GlobalSize = m_solids.Global_BSplineSolid.vertices_.size();
	//int m_GlobalSize = m_solids.Local_BSplineSolid[i].vertices_.size();
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.resize(m_GlobalSize, 0);

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	std::vector<double> m_1Intervel, m_2Intervel, m_3Intervel;
	std::vector<std::vector<double>> m_jacobiInverseMatrix;
	std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
	double ders_1[2], ders_2[2];
	std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
	double m_jacobidet;
	int m_BaseIndex, m_1Location, m_2Location;
	std::vector<double> m_1, m_2;
	int m_1Ptnum, m_2Ptnum, m_3Ptnum, m_MatrixSize;
	std::vector<std::vector<double>> m_BaseDValue1, m_BaseDValue2;//存储样条基关于u和v在高斯点的导数，
	std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵

	m_1Ptnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
	m_2Ptnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
	m_3Ptnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
	m_MatrixSize = m_1Ptnum * m_1Ptnum * m_3Ptnum;
	int uvw_size;
	if (uvw_id == 0)
		uvw_size = 0;
	else if (uvw_id == 1)
		uvw_size = m_1Ptnum - 1;

	m_BaseDValue1.resize(m_MatrixSize); m_BaseDValue2.resize(m_MatrixSize);
	m_1.resize(m_GaussNum); m_2.resize(m_GaussNum);
	m_solids.IntegralArea(i, m_1Intervel, m_2Intervel, m_3Intervel);
	/***************************构造刚度矩阵按照单元循环****************************************/
	for (int u_interval = 0; u_interval < m_1Intervel.size() - 1; u_interval++)
	{
		for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
			m_1[gaussnum] = (m_1Intervel[u_interval] + m_1Intervel[u_interval + 1] + (-m_1Intervel[u_interval]
				+ m_1Intervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
		//确定u方向这些高斯点所在区间的位置
		m_1Location = m_solids.X_LOCATION(i, 0, m_1[0]);

		for (int v_intervel = 0; v_intervel < m_2Intervel.size() - 1; v_intervel++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_2[gaussnum] = (m_2Intervel[v_intervel] + m_2Intervel[v_intervel + 1] + (-m_2Intervel[v_intervel]
					+ m_2Intervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定v方向这些高斯点所在区间的位置
			m_2Location = m_solids.X_LOCATION(i, 1, m_2[0]);

			m_det.clear(); m_InverseMatrix.clear();
			/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
			m_solids.JacobiMatrixFace(i, uvw, uvw_id, m_1, m_2, m_det, m_InverseMatrix);
			std::vector<int> m_baseindex;
			for (int ii = m_1Location - m_solids.Local_Order[i][0]; ii <= m_1Location; ii++)
				for (int j = m_2Location - m_solids.Local_Order[i][1]; j <= m_2Location; j++)
				{
					double m_tempvalue = 0;
					if (uvw == U)
						m_BaseIndex = m_solids.Control_Point_id[i][ii][j][uvw_size];
					else if (uvw == V)
						m_BaseIndex = m_solids.Control_Point_id[i][ii][uvw_size][j];
					else if (uvw == W)
						m_BaseIndex = m_solids.Control_Point_id[i][uvw_size][ii][j];

					m_baseindex.push_back(m_BaseIndex);
					m_BaseDValue1[m_BaseIndex].clear(); m_BaseDValue2[m_BaseIndex].clear();
					for (int m = 0; m < m_GaussNum; m++)
					{
						m_solids.DersOneBasisFun(i, 1, j, m_2[m], 1, ders_2);
						for (int l = 0; l < m_GaussNum; l++)
						{
							m_solids.DersOneBasisFun(i, 0, j, m_1[l], 1, ders_1);
							V3f gasv;
							if (uvw == U)
								gasv = m_solids.Geometry(i, uvw_id, m_1[l], m_2[m]);
							else if (uvw == V)
								gasv = m_solids.Geometry(i, m_1[l], uvw_id, m_2[m]);
							else if (uvw == W)
								gasv = m_solids.Geometry(i, m_1[l], m_2[m], uvw_id);

							double m_temp = ders_1[1] * ders_2[0];
							m_BaseDValue1[m_BaseIndex].push_back(m_temp);
							m_temp = ders_1[0] * ders_2[1];
							m_BaseDValue2[m_BaseIndex].push_back(m_temp);
							int m_index = m * m_GaussNum + l;
							m_tempvalue += ders_1[0] * ders_2[0] * RHS(gasv.x, gasv.y, gasv.z) *
								fabs(m_det[m_index]) * m_Weight[l] * m_Weight[m];
						}
					}
					/*m_bigrhs[m_solids.CP_Local_to_Global[i][m_BaseIndex]] += m_tempvalue * (m_1Intervel[u_interval + 1] -
						m_1Intervel[u_interval])*(m_2Intervel[v_intervel + 1] - m_2Intervel[v_intervel]) / 4;*/
					m_bigrhs[m_BaseIndex] += m_tempvalue * (m_1Intervel[u_interval + 1] -
						m_1Intervel[u_interval]) * (m_2Intervel[v_intervel + 1] - m_2Intervel[v_intervel]) / 4;
				}
			int m_ptNum = m_BaseDValue1[m_BaseIndex].size();
			for (int ii = 0; ii < m_baseindex.size(); ii++)
			{
				for (int j = 0; j < m_baseindex.size(); j++)
				{
					for (int k = 0; k < m_ptNum; k++)
					{
						int m_uptIndex = k % m_GaussNum;
						int m_vptIndex = k / m_GaussNum;
						int m_index;
						double m_Rx, m_Ry, m_Qx, m_Qy;               //存储基的梯度在区域上变换后的量。
						m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValue1[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValue2[m_baseindex[ii]][k];
						m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValue1[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValue2[m_baseindex[ii]][k];
						m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValue1[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValue2[m_baseindex[j]][k];
						m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValue1[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValue2[m_baseindex[j]][k];
						if (m_solids.CP_Local_to_Global[i][m_baseindex[ii]] == 0)
							int asdasd = 0;
						double test = (m_Rx * m_Qx + m_Ry * m_Qy) * fabs(m_det[k])
							* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * (m_1Intervel[u_interval + 1] -
								m_1Intervel[u_interval]) * (m_2Intervel[v_intervel + 1] - m_2Intervel[v_intervel]) / 4;
						m_bigmatrix[m_solids.CP_Local_to_Global[i][m_baseindex[ii]]][m_solids.CP_Local_to_Global[i][m_baseindex[j]]] += (m_Rx * m_Qx + m_Ry * m_Qy) * fabs(m_det[k])
							* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * (m_1Intervel[u_interval + 1] -
								m_1Intervel[u_interval]) * (m_2Intervel[v_intervel + 1] - m_2Intervel[v_intervel]) / 4;
						/*m_bigmatrix[m_baseindex[ii]][m_baseindex[j]] += (m_Rx*m_Qx + m_Ry * m_Qy)*fabs(m_det[k])
							*m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * (m_1Intervel[u_interval + 1] -
								m_1Intervel[u_interval])*(m_2Intervel[v_intervel + 1] - m_2Intervel[v_intervel]) / 4;*/
					}
				}
			}
		}
	}
}

void FunctionBSpline::facePDE(int solid_id, int uvw, int uvw_id, std::vector< std::vector<double> >& m_coef,
	std::vector< std::vector<double> > m_coeftest)
{
	/*-----------------构造刚度矩阵和右端项-----------------*/
	Matrix m_bigmatrix;
	std::vector<double> m_bigrhs;
	Global_MatrixAndRhsFace(solid_id, uvw, uvw_id, m_bigmatrix, m_bigrhs);
	/*-----------------处理刚度矩阵和右端项-----------------*/
	int m_GlobalSize = m_solids.Global_BSplineSolid.vertices_.size();
	std::vector<int> bdyedge;
	int ptnum = m_solids.Local_Node[solid_id][0].size() - m_solids.Local_Order[solid_id][0] - 1;
	int uvw_size;
	if (uvw_id == 0)
		uvw_size = 0;
	else if (uvw_id == 1)
		uvw_size = ptnum - 1;
	if (uvw == U)
	{
		for (int i = 0; i < ptnum; i++)
		{
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][0][i][uvw_size]);//第一条边
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][i][0][uvw_size]);//第二条边
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][i][1][uvw_size]);//第三条边
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][1][i][uvw_size]);//第四条边
		}

	}
	else if (uvw == V)
	{
		for (int i = 0; i < ptnum; i++)
		{
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][0][uvw_size][i]);//第一条边
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][i][uvw_size][0]);//第二条边
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][i][uvw_size][1]);//第三条边
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][1][uvw_size][i]);//第四条边
		}
	}
	else if (uvw == W)
	{
		for (int i = 0; i < ptnum; i++)
		{
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][uvw_size][0][i]);//第一条边
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][uvw_size][i][0]);//第二条边
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][uvw_size][i][1]);//第三条边
			bdyedge.push_back(m_solids.Control_Point_id[solid_id][uvw_size][1][i]);//第四条边
		}
	}
	for (int k = 0; k < bdyedge.size(); k++)
	{
		for (int j = 0; j < m_GlobalSize; j++)
		{
			m_bigrhs[j] = m_bigrhs[j] - m_coeftest[solid_id][bdyedge[k]] * m_bigmatrix[j][m_solids.CP_Local_to_Global[solid_id][bdyedge[k]]];
			//m_bigrhs[j] = m_bigrhs[j] - m_coeftest[solid_id][k] * m_bigmatrix[j][k];
		}
		for (int j = 0; j < m_GlobalSize; j++)
		{
			m_bigmatrix[j][m_solids.CP_Local_to_Global[solid_id][bdyedge[k]]] = 0;
			m_bigmatrix[m_solids.CP_Local_to_Global[solid_id][bdyedge[k]]][j] = 0;
			/*m_bigmatrix[j][k] = 0;
			m_bigmatrix[k][j] = 0;*/
		}
		m_bigmatrix[m_solids.CP_Local_to_Global[solid_id][bdyedge[k]]][m_solids.CP_Local_to_Global[solid_id][bdyedge[k]]] = 1;
		m_bigrhs[m_solids.CP_Local_to_Global[solid_id][bdyedge[k]]] = m_coeftest[solid_id][bdyedge[k]];
		/*m_bigmatrix[k][k] = 1;
		m_bigrhs[k] = m_coeftest[solid_id][k];*/
	}
	/*-----------------求解线型方程组-----------------*/
	std::vector<double> m_face;
	SolveLinearEquation(-2, m_bigmatrix, m_bigrhs, m_face);
	for (int v_id = 0; v_id < m_solids.Local_BSplineSolid[solid_id].vertices_.size(); v_id++)
	{
		m_coef[solid_id][v_id] = m_face[m_solids.CP_Local_to_Global[solid_id][v_id]];
		//m_coef[solid_id][v_id] = m_face[v_id];
	}

}

double FunctionBSpline::uTKu(std::vector<double> u, Matrix K0)
{
	int size = u.size();
	std::vector<double> cal_uTK(size, 0);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cal_uTK[i] += u[j] * K0[i][j];
		}
	}
	double cal_uTKu = 0;
	for (int i = 0; i < size; i++)
	{
		cal_uTKu += cal_uTK[i] * u[i];
	}
	return cal_uTKu;
}

double FunctionBSpline::PTu(std::vector<double> P, std::vector<double> u)
{
	double cal_PTu = 0;
	int size = P.size();
	for (int i = 0; i < size; i++)
	{
		cal_PTu += P[i] * u[i];
	}
	return cal_PTu;
}

void FunctionBSpline::Global_MatrixAndRhs1_new(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix
	, std::map<int, int>& base_check)
{
	int m_GlobalSize = m_solids.Global_size_without_bdy;

	local_matrix.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_matrix[i].resize(size, std::vector<double>(size, 0));
	}
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.resize(m_GlobalSize, 0);


	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	int base_num = 0;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵


		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								bool base_bdy_check = false;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								if (base_check.find(m_solids.Local_to_Global[i][m_BaseIndex]) == base_check.end()
									&& !m_solids.local_bdy[i][m_BaseIndex])
								{
									base_check.insert(std::pair<int, int>(m_solids.Local_to_Global[i][m_BaseIndex], base_num));
									base_num++;
									base_bdy_check = true;
								}
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								if (base_bdy_check == true)
									m_bigrhs[base_check[m_solids.Local_to_Global[i][m_BaseIndex]]] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double cal_temp = (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii]][m_baseindex[j]] += cal_temp;
								if (base_check.find(m_solids.Local_to_Global[i][m_baseindex[ii]]) != base_check.end()
									&& base_check.find(m_solids.Local_to_Global[i][m_baseindex[j]]) != base_check.end())
									m_bigmatrix[base_check[m_solids.Local_to_Global[i][m_baseindex[ii]]]][base_check[m_solids.Local_to_Global[i][m_baseindex[j]]]] +=
									cal_temp * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
						}
					}

				}
			}
		}
	}
}

void FunctionBSpline::Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix)
{
	int m_GlobalSize = m_solids.Global_size;

	local_matrix.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_matrix[i].resize(size, std::vector<double>(size, 0));
	}
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.resize(m_GlobalSize, 0);

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								m_bigrhs[m_solids.Local_to_Global[i][m_BaseIndex]] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double cal_temp = (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii]][m_baseindex[j]] += cal_temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]]][m_solids.Local_to_Global[i][m_baseindex[j]]] +=
									cal_temp * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
						}
					}
				}
			}
		}
	}
}

void FunctionBSpline::Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix)
{
	int m_GlobalSize = m_solids.Global_size;
	m_bigmatrix.clear();
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	for (int i = 0; i < local_matrix.size(); i++)
	{
		std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << local_matrix.size() << std::endl;
		for (int j = 0; j < local_matrix[i].size(); j++)
		{
			for (int w = 0; w < local_matrix[i][j].size(); w++)
			{
				m_bigmatrix[m_solids.Local_to_Global[i][j]][m_solids.Local_to_Global[i][w]] += local_matrix[i][j][w]
					* (umin + (1 - umin) * pow(material_distribution[i], p));
			}
		}
	}
}

void FunctionBSpline::Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix, std::vector<Matrix>& local_matrix1)
{
	int m_GlobalSize = m_solids.Global_size;
	m_bigmatrix.clear();
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	local_matrix1.clear();
	local_matrix1.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_matrix1[i].resize(size, std::vector<double>(size, 0));
	}
	for (int i = 0; i < local_matrix1.size(); i++)
	{
		std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << local_matrix1.size() << std::endl;
		for (int j = 0; j < local_matrix1[i].size(); j++)
		{
			for (int w = 0; w < local_matrix1[i][j].size(); w++)
			{
				local_matrix1[i][j][w] = local_matrix[i][j][w] * (umin + (1 - umin) * pow(material_distribution[i], p));
				m_bigmatrix[m_solids.Local_to_Global[i][j]][m_solids.Local_to_Global[i][w]] += local_matrix[i][j][w]
					* (umin + (1 - umin) * pow(material_distribution[i], p));
			}
		}
	}
}

void FunctionBSpline::Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix
	, std::vector<double>& m_bigrhs1)
{
	int m_GlobalSize = m_solids.Global_size;
	local_matrix.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_matrix[i].resize(size, std::vector<double>(size, 0));
	}
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.resize(m_GlobalSize, 0);
	m_bigrhs1.resize(m_GlobalSize, 0);


	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								m_bigrhs[m_solids.Local_to_Global[i][m_BaseIndex]] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								m_bigrhs1[m_solids.Local_to_Global[i][m_BaseIndex]] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double cal_temp = (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii]][m_baseindex[j]] += cal_temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]]][m_solids.Local_to_Global[i][m_baseindex[j]]] +=
									cal_temp * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
						}
					}
				}
			}
		}
	}
}

void FunctionBSpline::Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs
	, std::vector<double>& m_bigrhs1)
{
	int m_GlobalSize = m_solids.Global_size;
	m_bigmatrix.clear();
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.clear();
	m_bigrhs.resize(m_GlobalSize, 0);
	for (int i = 0; i < local_matrix.size(); i++)
	{
		std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << local_matrix.size() << std::endl;
		for (int j = 0; j < local_matrix[i].size(); j++)
		{
			for (int w = 0; w < local_matrix[i][j].size(); w++)
			{
				m_bigmatrix[m_solids.Local_to_Global[i][j]][m_solids.Local_to_Global[i][w]] += local_matrix[i][j][w]
					* (umin + (1 - umin) * pow(material_distribution[i], p));
			}
		}
	}
	for (int i = 0; i < m_GlobalSize; i++)
	{
		std::cout << "Quickly schedule of bigrhs " << i + 1 << ":" << m_GlobalSize << std::endl;
		m_bigrhs[i] = m_bigrhs1[i];
	}
}

void FunctionBSpline::Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix
	, std::vector<std::vector<double>>& local_bigrhs)
{
	int m_GlobalSize = m_solids.Global_size;
	local_matrix.resize(m_solids.Local_V3f.size());
	local_bigrhs.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_matrix[i].resize(size, std::vector<double>(size, 0));
		local_bigrhs[i].resize(size, 0);
	}
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.resize(m_GlobalSize, 0);



	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								local_bigrhs[i][m_BaseIndex] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								m_bigrhs[m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8) * (umin + (1 - umin) * pow(material_distribution[i], p));
								// m_bigrhs[m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval])*(m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
								// *(m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8);
							}
					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double cal_temp = (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii]][m_baseindex[j]] += cal_temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]]][m_solids.Local_to_Global[i][m_baseindex[j]]] +=
									cal_temp * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
						}
					}
				}
			}
		}
	}
}

void FunctionBSpline::Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs
	, std::vector<std::vector<double>>& local_bigrhs)
{
	int m_GlobalSize = m_solids.Global_size;
	m_bigmatrix.clear();
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.clear();
	m_bigrhs.resize(m_GlobalSize, 0);
	int local_matrix_size = local_matrix.size();
	for (int i = 0; i < local_matrix_size; i++)
	{
		std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << local_matrix_size << std::endl;
		int local_matrixi_size = local_matrix[i].size();
		for (int j = 0; j < local_matrixi_size; j++)
		{
			int local_matrixij_size = local_matrix[i][j].size();
			for (int w = 0; w < local_matrixij_size; w++)
			{
				m_bigmatrix[m_solids.Local_to_Global[i][j]][m_solids.Local_to_Global[i][w]] += local_matrix[i][j][w]
					* (umin + (1 - umin) * pow(material_distribution[i], p));
			}
		}
	}
	for (int i = 0; i < local_bigrhs.size(); i++)
	{
		for (int j = 0; j < local_bigrhs[i].size(); j++)
			m_bigrhs[m_solids.Local_to_Global[i][j]] += local_bigrhs[i][j] * (umin + (1 - umin) * pow(material_distribution[i], p));
		// m_bigrhs[m_solids.Local_to_Global[i][j]] += local_bigrhs[i][j];
	}
}

void FunctionBSpline::Global_MatrixAndRhs(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs
	, std::vector<std::vector<double>>& local_bigrhs)
{
	int m_GlobalSize = m_solids.Global_BSplineSolid.vertices_.size();
	m_bigmatrix.clear();
	m_bigmatrix.resize(m_GlobalSize, std::vector<double>(m_GlobalSize, 0));
	m_bigrhs.clear();
	m_bigrhs.resize(m_GlobalSize, 0);
	int local_matrix_size = local_matrix.size();
	for (int i = 0; i < local_matrix_size; i++)
	{
		std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << local_matrix_size << std::endl;
		int local_matrixi_size = local_matrix[i].size();
		for (int j = 0; j < local_matrixi_size; j++)
		{
			int local_matrixij_size = local_matrix[i][j].size();
			for (int w = 0; w < local_matrixij_size; w++)
			{
				m_bigmatrix[m_solids.CP_Local_to_Global[i][j]][m_solids.CP_Local_to_Global[i][w]] += local_matrix[i][j][w]
					* (umin + (1 - umin) * pow(material_distribution[i], p));
			}
		}
	}
	for (int i = 0; i < local_bigrhs.size(); i++)
	{
		for (int j = 0; j < local_bigrhs[i].size(); j++)
			m_bigrhs[m_solids.CP_Local_to_Global[i][j]] += local_bigrhs[i][j] * (umin + (1 - umin) * pow(material_distribution[i], p));
	}
}

void FunctionBSpline::Global_MatrixAndRhs_Sparse(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<double>& m_bigrhs, std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals)
{
	int m_GlobalSize = m_solids.Global_size;

	/*---------------存储每一行中的非零元下标------------*/
	std::vector<std::vector<int>> m_colids;
	m_colids.resize(m_GlobalSize);

	m_bigrhs.resize(m_GlobalSize);
	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			m_bigmatrix[m_solids.Local_to_Global[i][j]].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2), 0);
		}
	}

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵
		std::vector<Matrix> m_LocalMatrixes;

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		m_LocalMatrixes.resize((m_uIntervel.size() - 1) * (m_vIntervel.size() - 1) * (m_wIntervel.size() - 1));

		Matrix m_LocalMatrix;
		m_LocalMatrix.resize(m_MatrixSize);
		for (int k = 0; k < m_MatrixSize; k++)
		{
			m_LocalMatrix[k].resize(m_MatrixSize, 0);
		}

		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								m_bigrhs[m_solids.Local_to_Global[i][m_BaseIndex]] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double cal_temp = (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								m_LocalMatrix[m_baseindex[ii]][m_baseindex[j]] += cal_temp;
							}
						}
					}
				}
			}
		}
		//矩阵中的非零元素全部存起来
		for (int ii = 0; ii < m_MatrixSize; ii++)
		{
			for (int j = 0; j < m_MatrixSize; j++)
			{
				if (fabs(m_LocalMatrix[ii][j]) > 1.0e-7)
				{
					//如果m_faces[i].m_cptsGlobelID[j]不在m_faces[i].m_cptsGlobelID[ii]行的非零元下标中，
					//则把下标放入m_colids[m_faces[i].m_cptsGlobelID[ii]]，元素值放入矩阵中
					if (!InVecID(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]])) //全矩阵存储
					{
						m_colids[m_solids.Local_to_Global[i][ii]].push_back(m_solids.Local_to_Global[i][j]);
						int m_size = m_colids[m_solids.Local_to_Global[i][ii]].size();
						m_bigmatrix[m_solids.Local_to_Global[i][ii]][m_size - 1] += m_LocalMatrix[ii][j];
					}
					else
					{
						int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]]);
						m_bigmatrix[m_solids.Local_to_Global[i][ii]][tmp_id] += m_LocalMatrix[ii][j];
					}
				}
			}
		}
	}
	//处理边界条件
	//处理刚度矩阵和右端项
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		for (int j = 0; j < m_GlobalSize; j++)
		{
			for (int k = 0; k < m_colids[j].size(); k++)
				if ((m_colids[j][k]) == m_globalid)
				{
					m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid] * m_bigmatrix[j][k];
					m_bigmatrix[j][k] = 0;
					break;
				}

		}
		for (int k = 0; k < m_colids[m_globalid].size(); k++)
		{
			if (m_colids[m_globalid][k] == m_globalid)
			{
				m_bigmatrix[m_globalid][k] = 1;
			}
			else
				m_bigmatrix[m_globalid][k] = 0;
		}
		m_bigrhs[m_globalid] = m_bdycoef[m_solidid][m_baseid];
	}
	std::cout << "边界条件处理完毕" << std::endl;
	r_rowids.push_back(1);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first + 1);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num + 1);
	}
}

void FunctionBSpline::Global_MatrixAndRhs_SparseOld(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<double>& m_bigrhs, std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals)
{
	int m_GlobalSize = m_solids.Global_BSplineSolid.vertices_.size();

	/*---------------存储每一行中的非零元下标------------*/
	std::vector<std::vector<int>> m_colids;
	m_colids.resize(m_GlobalSize);

	m_bigrhs.resize(m_GlobalSize);
	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize);
	for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_BSplineSolid[i].vertices_.size(); j++)
		{
			m_bigmatrix[m_solids.CP_Local_to_Global[i][j]].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2), 0);
		}
	}

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	std::vector<int> rowids, colids;
	std::vector<double> elems;

	for (int i = 0; i < m_solids.Local_BSplineSolid.size(); i++)
	{
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵
		std::vector<Matrix> m_LocalMatrixes;

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		m_LocalMatrixes.resize((m_uIntervel.size() - 1) * (m_vIntervel.size() - 1) * (m_wIntervel.size() - 1));

		Matrix m_LocalMatrix;
		m_LocalMatrix.resize(m_MatrixSize);
		for (int k = 0; k < m_MatrixSize; k++)
		{
			m_LocalMatrix[k].resize(m_MatrixSize, 0);
		}

		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								m_BaseIndex = m_solids.Control_Point_id[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								m_bigrhs[m_solids.CP_Local_to_Global[i][m_BaseIndex]] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double cal_temp = (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								m_LocalMatrix[m_baseindex[ii]][m_baseindex[j]] += cal_temp;
							}
						}
					}
				}
			}
		}
		//矩阵中的非零元素全部存起来
		for (int ii = 0; ii < m_MatrixSize; ii++)
		{
			for (int j = 0; j < m_MatrixSize; j++)
			{
				if (fabs(m_LocalMatrix[ii][j]) > 1.0e-7)
				{
					//如果m_faces[i].m_cptsGlobelID[j]不在m_faces[i].m_cptsGlobelID[ii]行的非零元下标中，
					//则把下标放入m_colids[m_faces[i].m_cptsGlobelID[ii]]，元素值放入矩阵中
					if (!InVecID(m_solids.CP_Local_to_Global[i][j], m_colids[m_solids.CP_Local_to_Global[i][ii]])) //全矩阵存储
					{
						m_colids[m_solids.CP_Local_to_Global[i][ii]].push_back(m_solids.CP_Local_to_Global[i][j]);
						int m_size = m_colids[m_solids.CP_Local_to_Global[i][ii]].size();
						m_bigmatrix[m_solids.CP_Local_to_Global[i][ii]][m_size - 1] += m_LocalMatrix[ii][j];
					}
					else
					{
						int tmp_id = InVecIndex(m_solids.CP_Local_to_Global[i][j], m_colids[m_solids.CP_Local_to_Global[i][ii]]);
						m_bigmatrix[m_solids.CP_Local_to_Global[i][ii]][tmp_id] += m_LocalMatrix[ii][j];
					}
				}
			}
		}
	}
	//处理边界条件
	//处理刚度矩阵和右端项
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		for (int j = 0; j < m_GlobalSize; j++)
		{
			for (int k = 0; k < m_colids[j].size(); k++)
				if ((m_colids[j][k]) == m_globalid)
				{
					m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid] * m_bigmatrix[j][k];
					m_bigmatrix[j][k] = 0;
					break;
				}

		}
		for (int k = 0; k < m_colids[m_globalid].size(); k++)
		{
			if (m_colids[m_globalid][k] == m_globalid)
			{
				m_bigmatrix[m_globalid][k] = 1;
			}
			else
				m_bigmatrix[m_globalid][k] = 0;
		}
		m_bigrhs[m_globalid] = m_bdycoef[m_solidid][m_baseid];
	}
	std::cout << "边界条件处理完毕" << std::endl;
	r_rowids.push_back(1);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first + 1);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num + 1);
	}
}

void FunctionBSpline::Global_MatrixAndRhs_Sparse(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
	std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
	std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick)
{
	int m_GlobalSize = m_solids.Global_size;

	/*---------------存储每一行中的非零元下标------------*/
	m_colids.resize(m_GlobalSize);

	m_bigrhs.resize(m_GlobalSize);
	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			m_bigmatrix[m_solids.Local_to_Global[i][j]].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2), 0);
		}
	}
	local_matrix.resize(m_solids.Local_V3f.size());
	local_bigrhs.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_matrix[i].resize(size, std::vector<double>(size, 0));
		local_bigrhs[i].resize(size, 0);
	}

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	std::vector<int> rowids, colids;
	std::vector<double> elems;

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵
		std::vector<Matrix> m_LocalMatrixes;

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		m_LocalMatrixes.resize((m_uIntervel.size() - 1) * (m_vIntervel.size() - 1) * (m_wIntervel.size() - 1));

		// Matrix m_LocalMatrix;
		// m_LocalMatrix.resize(m_MatrixSize);
		// for (int k = 0; k < m_MatrixSize; k++)
		// {
			// m_LocalMatrix[k].resize(m_MatrixSize, 0);
		// }

		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								local_bigrhs[i][m_BaseIndex] += m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								m_bigrhs[m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8) * (umin + (1 - umin) * pow(material_distribution[i], p));
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double cal_temp = (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k])
									* m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii]][m_baseindex[j]] += cal_temp;
							}
						}
					}
				}
			}
		}
		// local_matrix.push_back(m_LocalMatrix);
		//矩阵中的非零元素全部存起来
		for (int ii = 0; ii < m_MatrixSize; ii++)
		{
			for (int j = 0; j < m_MatrixSize; j++)
			{
				if (fabs(local_matrix[i][ii][j]) > 1.0e-7)
				{
					//如果m_faces[i].m_cptsGlobelID[j]不在m_faces[i].m_cptsGlobelID[ii]行的非零元下标中，
					//则把下标放入m_colids[m_faces[i].m_cptsGlobelID[ii]]，元素值放入矩阵中
					if (!InVecID(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]])) //全矩阵存储
					{
						m_colids[m_solids.Local_to_Global[i][ii]].push_back(m_solids.Local_to_Global[i][j]);
						int m_size = m_colids[m_solids.Local_to_Global[i][ii]].size();
						m_bigmatrix[m_solids.Local_to_Global[i][ii]][m_size - 1] += local_matrix[i][ii][j] * (umin + (1 - umin) * pow(material_distribution[i], p));
					}
					else
					{
						int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]]);
						m_bigmatrix[m_solids.Local_to_Global[i][ii]][tmp_id] += local_matrix[i][ii][j] * (umin + (1 - umin) * pow(material_distribution[i], p));
					}
				}
			}
		}
	}
	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		std::vector<std::pair<int, int>> temp_bdyquick;
		std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		for (int j = 0; j < m_GlobalSize; j++)
		{
			for (int k = 0; k < m_colids[j].size(); k++)
				if ((m_colids[j][k]) == m_globalid)
				{
					temp_bdyquick.push_back(std::pair<int, int>(j, k));
					m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid] * m_bigmatrix[j][k];
					m_bigmatrix[j][k] = 0;
					break;
				}

		}
		m_bdyquick.push_back(temp_bdyquick);
		for (int k = 0; k < m_colids[m_globalid].size(); k++)
		{
			if (m_colids[m_globalid][k] == m_globalid)
			{
				m_bigmatrix[m_globalid][k] = 1;
			}
			else
				m_bigmatrix[m_globalid][k] = 0;
		}
		m_bigrhs[m_globalid] = m_bdycoef[m_solidid][m_baseid];
	}
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(1);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first + 1);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num + 1);
	}
}

void FunctionBSpline::Global_MatrixAndRhs_Sparse_quick(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
	std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
	std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick)
{
	int m_GlobalSize = m_solids.Global_size;
	m_bigrhs.clear();
	m_bigrhs.resize(m_GlobalSize, 0);

	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			m_bigmatrix[m_solids.Local_to_Global[i][j]].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2), 0);
		}
	}

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		int m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		int m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		int m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		int m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
		//矩阵中的非零元素全部存起来
		for (int ii = 0; ii < m_MatrixSize; ii++)
		{
			for (int j = 0; j < m_MatrixSize; j++)
			{
				if (fabs(local_matrix[i][ii][j]) > 1.0e-7)
				{
					//如果m_faces[i].m_cptsGlobelID[j]不在m_faces[i].m_cptsGlobelID[ii]行的非零元下标中，
					//则把下标放入m_colids[m_faces[i].m_cptsGlobelID[ii]]，元素值放入矩阵中
					if (!InVecID(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]])) //全矩阵存储
					{
						m_colids[m_solids.Local_to_Global[i][ii]].push_back(m_solids.Local_to_Global[i][j]);
						int m_size = m_colids[m_solids.Local_to_Global[i][ii]].size();
						m_bigmatrix[m_solids.Local_to_Global[i][ii]][m_size - 1] += local_matrix[i][ii][j] * (umin + (1 - umin) * pow(material_distribution[i], p));
					}
					else
					{
						int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]]);
						m_bigmatrix[m_solids.Local_to_Global[i][ii]][tmp_id] += local_matrix[i][ii][j] * (umin + (1 - umin) * pow(material_distribution[i], p));
					}
				}
			}
		}
	}
	for (int i = 0; i < local_bigrhs.size(); i++)
	{
		for (int j = 0; j < local_bigrhs[i].size(); j++)
			m_bigrhs[m_solids.Local_to_Global[i][j]] += local_bigrhs[i][j] * (umin + (1 - umin) * pow(material_distribution[i], p));
	}
	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		std::vector<std::pair<int, int>> temp_bdyquick = m_bdyquick[i];
		for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
		{
			int j = (*temp_bdyquick_it).first;
			int k = (*temp_bdyquick_it).second;
			m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid] * m_bigmatrix[j][k];
			m_bigmatrix[j][k] = 0;
		}
		/*for (int j = 0; j < m_GlobalSize; j++)
		{
			for (int k = 0; k < m_colids[j].size(); k++)
				if ((m_colids[j][k]) == m_globalid)
				{
					m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid] * m_bigmatrix[j][k];
					m_bigmatrix[j][k] = 0;
					break;
				}

		}*/
		for (int k = 0; k < m_colids[m_globalid].size(); k++)
		{
			if (m_colids[m_globalid][k] == m_globalid)
			{
				m_bigmatrix[m_globalid][k] = 1;
			}
			else
				m_bigmatrix[m_globalid][k] = 0;
		}
		m_bigrhs[m_globalid] = m_bdycoef[m_solidid][m_baseid];
	}
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(1);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first + 1);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num + 1);
	}
}

void FunctionBSpline::Global_MatrixAndRhs_liner(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix,
	std::vector<std::vector<double>>& local_bigrhs)
{
	int m_GlobalSize = m_solids.Global_size;
	local_matrix.resize(m_solids.Local_V3f.size());
	local_bigrhs.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_matrix[i].resize(size * 3, std::vector<double>(size * 3, 0));
		local_bigrhs[i].resize(size * 3, 0);
	}
	m_bigmatrix.resize(m_GlobalSize * 3, std::vector<double>(m_GlobalSize * 3, 0));
	m_bigrhs.resize(m_GlobalSize * 3, 0);

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue1 = 0, m_tempvalue2 = 0, m_tempvalue3 = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue1 += u_ders[0] * v_ders[0] * w_ders[0] * f1(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
											m_tempvalue2 += u_ders[0] * v_ders[0] * w_ders[0] * f2(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
											m_tempvalue3 += u_ders[0] * v_ders[0] * w_ders[0] * f3(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								double temp = (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_bigrhs[i][3 * m_BaseIndex] += m_tempvalue1 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue1 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
								local_bigrhs[i][3 * m_BaseIndex + 1] += m_tempvalue2 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 1] += (m_tempvalue2 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
								local_bigrhs[i][3 * m_BaseIndex + 2] += m_tempvalue3 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 2] += (m_tempvalue3 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double temp = fabs(m_det[k]) * m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1]
									- m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii] * 3][m_baseindex[j] * 3] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
								local_matrix[i][m_baseindex[ii] * 3][m_baseindex[j] * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
								local_matrix[i][m_baseindex[ii] * 3][m_baseindex[j] * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]] * 3][m_solids.Local_to_Global[i][m_baseindex[j]] * 3] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]] * 3][m_solids.Local_to_Global[i][m_baseindex[j]] * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]] * 3][m_solids.Local_to_Global[i][m_baseindex[j]] * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]] * 3 + 1][m_solids.Local_to_Global[i][m_baseindex[j]] * 3] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]] * 3 + 1][m_solids.Local_to_Global[i][m_baseindex[j]] * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]] * 3 + 1][m_solids.Local_to_Global[i][m_baseindex[j]] * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]] * 3 + 2][m_solids.Local_to_Global[i][m_baseindex[j]] * 3] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]] * 3 + 2][m_solids.Local_to_Global[i][m_baseindex[j]] * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
								m_bigmatrix[m_solids.Local_to_Global[i][m_baseindex[ii]] * 3 + 2][m_solids.Local_to_Global[i][m_baseindex[j]] * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;
							}
						}
					}
				}
			}
		}
	}
}

void FunctionBSpline::Global_MatrixAndRhs_liner_quick(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix,
	std::vector<std::vector<double>>& local_bigrhs)
{
	int m_GlobalSize = m_solids.Global_size;
	Matrix().swap(m_bigmatrix);
	m_bigmatrix.resize(m_GlobalSize * 3, std::vector<double>(m_GlobalSize * 3, 0));
	std::vector<double>().swap(m_bigrhs);
	m_bigrhs.resize(m_GlobalSize * 3, 0);
	int local_matrix_size = local_matrix.size();
	for (int i = 0; i < local_matrix_size; i++)
	{
		std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << local_matrix_size << std::endl;
		for (int j = 0; j < local_matrix[i].size(); j++)
		{
			for (int w = 0; w < local_matrix[i][j].size(); w++)
			{
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++)
						m_bigmatrix[3 * m_solids.Local_to_Global[i][j / 3] + k][3 * m_solids.Local_to_Global[i][w / 3] + l] += local_matrix[i][j / 3 * 3 + k][w / 3 * 3 + l]
						* (umin + (1 - umin) * pow(material_distribution[i], p));
			}
		}
	}
	for (int i = 0; i < local_bigrhs.size(); i++)
		for (int j = 0; j < local_bigrhs[i].size(); j++)
			for (int k = 0; k < 3; k++)
				m_bigrhs[3 * m_solids.Local_to_Global[i][j / 3] + k] += local_bigrhs[i][j / 3 * 3 + k] * (umin + (1 - umin) * pow(material_distribution[i], p));
}

void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids, std::vector<int>& r_colids,
	std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick)
{
	int m_GlobalSize = m_solids.Global_size;        // total control size   < 64 * patches
	
	m_colids.resize(m_GlobalSize * 3);
	m_bigrhs.resize(m_GlobalSize * 3);

	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize * 3);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			for (int k = 0; k < 3; k++)
				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2) * martix_size, 0);
		}
	}
	local_bigrhs.resize(m_solids.Local_V3f.size());
	local_matrix.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_bigrhs[i].resize(size * 3, 0);
		local_matrix[i].resize(size * 3, std::vector<double>(size * 3, 0));
	}
	// local_matrix.reserve(m_solids.Local_V3f.size());

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;

//#pragma omp parallel for num_threads(8)
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		if (i % 100 == 0)    std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;


		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);

		// Matrix m_LocalMatrix;
		// m_LocalMatrix.resize(m_MatrixSize * 3);
		// for (int k = 0; k < m_MatrixSize * 3; k++)
		// {
			// m_LocalMatrix[k].resize(m_MatrixSize * 3, 0);
		// }

		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;

					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue1 = 0, m_tempvalue2 = 0, m_tempvalue3 = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue1 += u_ders[0] * v_ders[0] * w_ders[0] * f1(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
											m_tempvalue2 += u_ders[0] * v_ders[0] * w_ders[0] * f2(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
											m_tempvalue3 += u_ders[0] * v_ders[0] * w_ders[0] * f3(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								double temp = (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_bigrhs[i][3 * m_BaseIndex] += m_tempvalue1 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue1 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
								local_bigrhs[i][3 * m_BaseIndex + 1] += m_tempvalue2 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 1] += (m_tempvalue2 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
								local_bigrhs[i][3 * m_BaseIndex + 2] += m_tempvalue3 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 2] += (m_tempvalue3 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;

								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double temp = fabs(m_det[k]) * m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1]
									- m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii] * 3][m_baseindex[j] * 3] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
								local_matrix[i][m_baseindex[ii] * 3][m_baseindex[j] * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
								local_matrix[i][m_baseindex[ii] * 3][m_baseindex[j] * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;
							}
						}
					}
				}
			}
		}
		// local_matrix.push_back(m_LocalMatrix);
		//矩阵中的非零元素全部存起来
		for (int ii = 0; ii < m_MatrixSize; ii++)
		{
			for (int j = 0; j < m_MatrixSize; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						if (fabs(local_matrix[i][ii * 3 + k][j * 3 + l]) > 1.0e-7)
						{
							if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
							{
								m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
								int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//#pragma omp atomic
								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrix[i][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
							else
							{
								int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//#pragma omp atomic
								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrix[i][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
						}
					}
				}
			}
		}
	}
	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		std::vector<std::pair<int, int>> temp_bdyquick;
		if (i % 100 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		for (int j = 0; j < m_GlobalSize * 3; j++)
		{
			for (int k = 0; k < m_colids[j].size(); k++)
			{
				if (bdybase_count < subnum)
				{
					if ((m_colids[j][k]) == m_globalid * 3 + forceorder)
					{
						temp_bdyquick.push_back(std::pair<int, int>(j, k));
						m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
						m_bigmatrix[j][k] = 0;
						break;
					}
				}
				else
				{
					for (int l = 0; l < 3; l++)
					{
						if ((m_colids[j][k]) == m_globalid * 3 + l)
						{
							temp_bdyquick.push_back(std::pair<int, int>(j, k));
							m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
							m_bigmatrix[j][k] = 0;
							break;
						}
					}
				}
			}
			// for (int k = 0; k < m_colids[j].size(); k++)
			// {
				// if (bdybase_count < 4)
				// {
					// if ((m_colids[j][k]) == m_globalid * 3 + 1)
					// {
						// temp_bdyquick.push_back(std::pair<int, int>(j, k));
						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + 1] * m_bigmatrix[j][k];
						// m_bigmatrix[j][k] = 0;
						// break;
					// }
				// }
				// else if (bdybase_count < 8)
				// {
					// if ((m_colids[j][k]) == m_globalid * 3)
					// {
						// temp_bdyquick.push_back(std::pair<int, int>(j, k));
						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3] * m_bigmatrix[j][k];
						// m_bigmatrix[j][k] = 0;
						// break;
					// }
				// }
				// else
				// {
					// for (int l = 0; l < 3; l++)
					// {
						// if ((m_colids[j][k]) == m_globalid * 3 + l)
						// {
							// temp_bdyquick.push_back(std::pair<int, int>(j, k));
							// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
							// m_bigmatrix[j][k] = 0;
							// break;
						// }
					// }
				// }
			// }
		}
		m_bdyquick.push_back(temp_bdyquick);
		if (bdybase_count < subnum)
		{
			for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
			{
				if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
				{
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
				}
				else
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
			}
			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
		}
		else
		{
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}
		}
		// if (bdybase_count < 4)
		// {
			// for (int k = 0; k < m_colids[m_globalid * 3 + 1].size(); k++)
			// {
				// if (m_colids[m_globalid * 3 + 1][k] == m_globalid * 3 + 1)
				// {
					// m_bigmatrix[m_globalid * 3 + 1][k] = 1;
				// }
				// else
					// m_bigmatrix[m_globalid * 3 + 1][k] = 0;
			// }
			// m_bigrhs[m_globalid * 3 + 1] = m_bdycoef[m_solidid][m_baseid * 3 + 1];
		// }
		// else if (bdybase_count < 8)
		// {
			// for (int k = 0; k < m_colids[m_globalid * 3].size(); k++)
			// {
				// if (m_colids[m_globalid * 3][k] == m_globalid * 3)
				// {
					// m_bigmatrix[m_globalid * 3 ][k] = 1;
				// }
				// else
					// m_bigmatrix[m_globalid * 3][k] = 0;
			// }
			// m_bigrhs[m_globalid * 3] = m_bdycoef[m_solidid][m_baseid * 3];
		// }
		// else
		// {
			// for (int m = 0; m < 3; m++)
			// {
				// for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				// {
					// if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						// || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					// {
						// m_bigmatrix[m_globalid * 3 + m][k] = 1;
					// }
					// else
						// m_bigmatrix[m_globalid * 3 + m][k] = 0;
				// }
				// m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			// }
		// }
		bdybase_count++;
	}
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(0);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize * 3; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first );
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num);
	}
	Matrix().swap(m_bigmatrix);
	std::vector<double>().swap(m_Weight);
	std::vector<double>().swap(m_IntePt);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	os << "[";
	for (size_t i = 0; i < v.size(); ++i) {
		os << v[i];
		if (i != v.size() - 1)
			os << ", ";
	}
	os << "]";
	return os;
}

// Linear_Elastic_CC_IGA_Solver

/// <summary>
/// 参照王宁改的版本
/// 效率待更新 先保证对 暂时未保证
/// </summary>
/// <param name="m_bdycoef"></param>
/// <param name="bdybaseid">边界条件</param>
/// <param name="local_matrix">局部刚度矩阵</param>
/// <param name="m_bigrhs">右端项</param>
/// <param name="local_bigrhs">局部右端项</param>
/// <param name="r_rowids">未使用</param>
/// <param name="r_colids">未使用</param>
/// <param name="r_vals">未使用</param>
/// <param name="m_colids"></param>
/// <param name="m_bdyquick">未使用</param>
/// <param name="u"></param>
void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_V2(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids, std::vector<int>& r_colids,
	std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick, std::vector<double>& u)
{


	const int global_basis_cnt = m_solids.Global_size;          // 全局基函数数量
	const int global_dof_cnt = m_solids.Global_size * 3;        // 全局自由度数量
	const int patch_cnt = m_solids.Local_V3f.size();            // Bezier patch数量
	const int local_basis_cnt = 64;
	const int local_dof_cnt = local_basis_cnt * 3;              // 局部自由度数量
	std::cout << "global_dof_cnt = " << global_dof_cnt << "\n";


	// 初始化右端项
	m_bigrhs.resize(global_dof_cnt, 0);
	m_colids.resize(global_dof_cnt);

	// 初始化刚度矩阵
	Matrix m_bigmatrix;
	m_bigmatrix.resize(global_dof_cnt);
	local_matrix.resize(patch_cnt);
	for (int i = 0; i < patch_cnt; i++)
	{
		local_matrix[i].resize(local_dof_cnt, std::vector<double>(local_dof_cnt, 0));
	}

	// 初始化高斯点值与权重
	std::vector<double> gauss_x4 = { -0.8611363116, -0.3399810436, 0.3399810436,  0.8611363116 };
	std::vector<double> gauss_w4 = { 0.3478548451 ,0.6521451549, 0.6521451549, 0.3478548451 };
	int gauss_cnt = 4;
	int patch_gauss_cnt = gauss_cnt * gauss_cnt * gauss_cnt;

	std::vector<double> para_gauss{ 0.0694318442,0.3300094782,0.6699905218,0.9305681558 };

	std::vector<std::vector<std::vector<double>>> pre_calc_ders;
	pre_calc_ders.resize(4);
	for (int i = 0; i < 4; ++i) {
		pre_calc_ders[i].resize(gauss_cnt);
		if (i == 0) {
			auto func1 = [](double x)->double {
				return -3 * x * x + 6 * x - 3;
			};
			auto func0 = [](double x)->double {
				return (1 - x) * (1 - x) * (1 - x);
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
		if (i == 1) {
			auto func1 = [](double x)->double {
				return 9 * x * x - 12 * x + 3;
			};
			auto func0 = [](double x)->double {
				return 3 * x * x * x - 6 * x * x + 3 * x;
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
		if (i == 2) {
			auto func1 = [](double x)->double {
				return -9 * x * x + 6 * x;
			};
			auto func0 = [](double x)->double {
				return (1 - x) * 3 * x * x;
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
		if (i == 3) {
			auto func1 = [](double x)->double {
				return 3 * x * x;
			};
			auto func0 = [](double x)->double {
				return (x) * (x) * (x);
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
	}

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;

	//#pragma omp parallel for num_threads(8)
	for (int patch_id = 0; patch_id < patch_cnt; patch_id++)
	{
		if (patch_id % 100 == 0)    std::cout << "Schedule of stiffness matrix " << patch_id + 1 << ":" << patch_cnt << std::endl;

		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;                                           //每个基对应一组在各个高斯点上的雅可比逆矩阵

		std::vector<double> m_det;                                                     //记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;


		int u_basis_cnt = 4;
		int v_basis_cnt = 4;
		int w_basis_cnt = 4;



		std::vector<double> u_Intervel{ 0, 1 };
		std::vector<double> v_Intervel{ 0, 1 };
		std::vector<double> w_Intervel{ 0, 1 };

		int elem_u_cnt = u_Intervel.size() - 1;
		int elem_v_cnt = v_Intervel.size() - 1;
		int elem_w_cnt = w_Intervel.size() - 1;



		/***************************构造刚度矩阵按照单元循环****************************************/
		m_det.clear(); m_InverseMatrix.clear();
		m_solids.JacobiMatrix_new(patch_id, para_gauss, para_gauss, para_gauss, m_det, m_InverseMatrix);

		std::vector<int> patch_global_base_index;

		for (int u_id = 0; u_id < elem_u_cnt; u_id++)
		{
			for (int v_id = 0; v_id < elem_v_cnt; v_id++)
			{
				for (int w_id = 0; w_id < elem_w_cnt; w_id++)
				{

					//std::cout << "Calc Derivative \n";
					// 第i个基第j个高斯点的导数
					std::vector<std::vector<double>> u_basis_derivative(local_basis_cnt, std::vector<double>(patch_gauss_cnt, 0));
					std::vector<std::vector<double>> v_basis_derivative(local_basis_cnt, std::vector<double>(patch_gauss_cnt, 0));
					std::vector<std::vector<double>> w_basis_derivative(local_basis_cnt, std::vector<double>(patch_gauss_cnt, 0));



					int cur_basis_id = 0;
					for (int u_basis_id = 0; u_basis_id < u_basis_cnt; u_basis_id++) {
						for (int v_basis_id = 0; v_basis_id < v_basis_cnt; v_basis_id++) {
							for (int w_basis_id = 0; w_basis_id < w_basis_cnt; w_basis_id++) {

								patch_global_base_index.push_back(m_solids.Control_Point[patch_id][w_basis_id][v_basis_id][u_basis_id]);

								int cur_gauss_id = 0;
								for (int w_gauss_id = 0; w_gauss_id < gauss_cnt; w_gauss_id++)
								{
									auto w_ders = pre_calc_ders[w_basis_id][w_gauss_id];

									for (int v_gauss_id = 0; v_gauss_id < gauss_cnt; v_gauss_id++)
									{
										auto v_ders = pre_calc_ders[v_basis_id][v_gauss_id];
										for (int u_gauss_id = 0; u_gauss_id < gauss_cnt; u_gauss_id++)
										{
											auto u_ders = pre_calc_ders[u_basis_id][u_gauss_id];

											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											u_basis_derivative[cur_basis_id][cur_gauss_id] = m_temp;
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											v_basis_derivative[cur_basis_id][cur_gauss_id] = m_temp;
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											w_basis_derivative[cur_basis_id][cur_gauss_id] = m_temp;
											cur_gauss_id++;
										}
									}
								}
								cur_basis_id++;
							}
						}
					}


					//std::cout << "Local Matrix \n";
					for (int i = 0; i < patch_global_base_index.size(); i++)
					{
						int ii = patch_global_base_index[i];
						for (int j = 0; j < patch_global_base_index.size(); j++)
						{
							int jj = patch_global_base_index[j];
							for (int k = 0; k < patch_gauss_cnt; k++)       // 64
							{
								int m_uptIndex = (k % (gauss_cnt * gauss_cnt)) % gauss_cnt;
								int m_vptIndex = (k % (gauss_cnt * gauss_cnt)) / gauss_cnt;
								int m_wptIndex = k / (gauss_cnt * gauss_cnt);
								int m_index;

								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;
								m_Rx = m_InverseMatrix[k][0][0] * u_basis_derivative[ii][k] + m_InverseMatrix[k][0][1] * v_basis_derivative[ii][k] + m_InverseMatrix[k][0][2] * w_basis_derivative[ii][k];
								m_Ry = m_InverseMatrix[k][1][0] * u_basis_derivative[ii][k] + m_InverseMatrix[k][1][1] * v_basis_derivative[ii][k] + m_InverseMatrix[k][1][2] * w_basis_derivative[ii][k];
								m_Rz = m_InverseMatrix[k][2][0] * u_basis_derivative[ii][k] + m_InverseMatrix[k][2][1] * v_basis_derivative[ii][k] + m_InverseMatrix[k][2][2] * w_basis_derivative[ii][k];
								m_Qx = m_InverseMatrix[k][0][0] * u_basis_derivative[jj][k] + m_InverseMatrix[k][0][1] * v_basis_derivative[jj][k] + m_InverseMatrix[k][0][2] * w_basis_derivative[jj][k];
								m_Qy = m_InverseMatrix[k][1][0] * u_basis_derivative[jj][k] + m_InverseMatrix[k][1][1] * v_basis_derivative[jj][k] + m_InverseMatrix[k][1][2] * w_basis_derivative[jj][k];
								m_Qz = m_InverseMatrix[k][2][0] * u_basis_derivative[jj][k] + m_InverseMatrix[k][2][1] * v_basis_derivative[jj][k] + m_InverseMatrix[k][2][2] * w_basis_derivative[jj][k];
								double temp = fabs(m_det[k])
									* gauss_w4[m_uptIndex] * gauss_w4[m_vptIndex] * gauss_w4[m_wptIndex]
									* 1
									* 1
									* 1 / 8;

								local_matrix[patch_id][ii * 3 + 0][jj * 3 + 0] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
								local_matrix[patch_id][ii * 3 + 0][jj * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
								local_matrix[patch_id][ii * 3 + 0][jj * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
								local_matrix[patch_id][ii * 3 + 1][jj * 3 + 0] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
								local_matrix[patch_id][ii * 3 + 1][jj * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
								local_matrix[patch_id][ii * 3 + 1][jj * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
								local_matrix[patch_id][ii * 3 + 2][jj * 3 + 0] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
								local_matrix[patch_id][ii * 3 + 2][jj * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
								local_matrix[patch_id][ii * 3 + 2][jj * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;
							}
						}
					}
				}
			}
		}
		//std::cout << "Global Matrix \n";
		// local_matrix.push_back(m_LocalMatrix);
		//矩阵中的非零元素全部存起来
		for (int i = 0; i < patch_global_base_index.size(); i++)
		{
			for (int j = 0; j < patch_global_base_index.size(); j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						//std::cout << i << " " << j << " " << k << " " << l << "\n";
						int local_id1 = i * 3 + k;
						int local_id2 = j * 3 + l;
						int global_id1 = m_solids.Local_to_Global[patch_id][patch_global_base_index[i]] * 3 + k;
						int global_id2 = m_solids.Local_to_Global[patch_id][patch_global_base_index[j]] * 3 + l;

						if (fabs(local_matrix[patch_id][local_id1][local_id2]) > 1.0e-12)
						{

							if (!InVecID(global_id2, m_colids[global_id1]))
							{
								m_colids[global_id1].push_back(global_id2);
								m_bigmatrix[global_id1].push_back(0);
								m_bigmatrix[global_id1].back() += local_matrix[patch_id][local_id1][local_id2];
							}
							else
							{
								int tmp_id = InVecIndex(global_id2, m_colids[global_id1]);
								assert(tmp_id < m_bigmatrix[global_id1].size() && tmp_id >= 0);
								m_bigmatrix[global_id1][tmp_id] += local_matrix[patch_id][local_id1][local_id2];
							}
						}
					}
				}
			}
		}
		//std::cout << "Global Matrix Done once\n";
	}

	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;

	//std::cout << "bdybaseid.size() = " << bdybaseid.size() << std::endl;
	std::unordered_set<int> fix_base;

	for (int i = 0; i < bdybaseid.size(); ++i) {

		int m_globalid = bdybaseid[i].global_id;
		int m_bdytype = bdybaseid[i].bdy_type;
		if (m_bdytype != 1) continue;
		for (int l = 0; l < 3; ++l) {
			fix_base.insert(m_globalid * 3 + l);
		}
	}

	for (int j = 0; j < global_dof_cnt; j++)
	{
		//std::cout << j << std::endl;
		for (int k = 0; k < m_colids[j].size(); k++)
		{
			if (fix_base.count(m_colids[j][k]))
			{
				m_bigmatrix[j][k] = 0;
			}
		}
	}

	for (int i = 0; i < bdybaseid.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		int m_bdytype = bdybaseid[i].bdy_type;
		//std::cout << "bdybaseid[i].bdy_type = " << bdybaseid[i].bdy_type << std::endl;
		// 
		//for (int j = 0; j < m_GlobalSize * 3; j++)
		//{
		//	for (int k = 0; k < m_colids[j].size(); k++)
		//	{
		//		//if (bdybase_count < subnum)
		//		//{
		//		//	if ((m_colids[j][k]) == m_globalid * 3 + forceorder)
		//		//	{
		//		//		temp_bdyquick.push_back(std::pair<int, int>(j, k));
		//		//		m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
		//		//		m_bigmatrix[j][k] = 0;
		//		//		break;
		//		//	}
		//		//}
		//		//else
		//		//{
		//		//	for (int l = 0; l < 3; l++)
		//		//	{
		//		//		if ((m_colids[j][k]) == m_globalid * 3 + l)
		//		//		{
		//		//			temp_bdyquick.push_back(std::pair<int, int>(j, k));
		//		//			m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
		//		//			m_bigmatrix[j][k] = 0;
		//		//			break;
		//		//		}
		//		//	}
		//		//}
		//		if (m_bdytype == 1) {
		//			// 支撑
		//			for (int l = 0; l < 3; l++)
		//			{
		//				if ((m_colids[j][k]) == m_globalid * 3 + l)
		//				{
		//					m_bigmatrix[j][k] = 0;
		//					break;
		//				}
		//			}

		//		}
		//	}
		//}

		if (m_bdytype == 2)
		{
			/*for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
			{
				if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
				{
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
				}
				else
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
			}*/
			// 施加力
			for (int dir = 0; dir < 3; ++dir) {
				m_bigrhs[m_globalid * 3 + dir] = m_bdycoef[m_solidid][m_baseid * 3 + dir];
			}
		}
		else if (m_bdytype == 1)
		{
			/*for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}*/
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + m)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = 0;
			}
		}
		bdybase_count++;
	}

	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	std::vector<T> coeffs;
	for (int i = 0; i < global_dof_cnt; ++i) {
		for (int k = 0; k < m_colids[i].size(); k++)
		{
			if (fabs(m_bigmatrix[i][k]) > 1.0e-12)
				coeffs.push_back(T(i, m_colids[i][k], m_bigmatrix[i][k]));
		}
	}
	{
		Timer timer;
		std::cout << "prepare solving...\n";
		Eigen::setNbThreads(8);  // 线程池大小
		Eigen::initParallel();   // 初始化多线程


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

		u = std::vector<double>(eigen_u.data(), eigen_u.data() + eigen_u.size());

		std::cout << "end solving...\n";
	}


	/*std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(0);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize * 3; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-9)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num);
	}*/
	Matrix().swap(m_bigmatrix);
}

void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_V1(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids, std::vector<int>& r_colids,
	std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick, std::vector<double>& u)
{
	const int global_basis_cnt = m_solids.Global_size;          // 全局基函数数量
	const int global_dof_cnt = m_solids.Global_size * 3;        // 全局自由度数量
	const int patch_cnt = m_solids.Local_V3f.size();            // Bezier patch数量
	const int local_basis_cnt = 64;
	const int local_dof_cnt = local_basis_cnt * 3;              // 局部自由度数量
	std::cout << "global_dof_cnt = " << global_dof_cnt << "\n";


	// 初始化右端项
	m_bigrhs.resize(global_dof_cnt, 0);
	m_colids.resize(global_dof_cnt);

	// 初始化刚度矩阵
	Matrix m_bigmatrix;
	m_bigmatrix.resize(global_dof_cnt);
	local_matrix.resize(patch_cnt);
	for (int i = 0; i < patch_cnt; i++)
	{
		local_matrix[i].resize(local_dof_cnt, std::vector<double>(local_dof_cnt, 0));
	}

	// 初始化高斯点值与权重
	std::vector<double> gauss_x4 = { -0.8611363116, -0.3399810436, 0.3399810436,  0.8611363116 };
	std::vector<double> gauss_w4 = { 0.3478548451 ,0.6521451549, 0.6521451549, 0.3478548451 };
	int gauss_cnt = 4;
	int patch_gauss_cnt = gauss_cnt * gauss_cnt * gauss_cnt;

	std::vector<double> para_gauss{ 0.0694318442,0.3300094782,0.6699905218,0.9305681558 };

	std::vector<std::vector<std::vector<double>>> pre_calc_ders;

	pre_calc_ders.resize(4);


	for (int i = 0; i < 4; ++i) {
		pre_calc_ders[i].resize(gauss_cnt);
		if (i == 0) {
			auto func1 = [](double x)->double{
				return -3 * x * x + 6 * x - 3;
			};
			auto func0 = [](double x)->double {
				return (1-x)*(1-x)*(1-x);
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
		if (i == 1) {
			auto func1 = [](double x)->double {
				return 9 * x * x - 12 * x + 3;
			};
			auto func0 = [](double x)->double {
				return 3 * x * x * x - 6 * x * x + 3 * x ;
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
		if (i == 2) {
			auto func1 = [](double x)->double {
				return -9 * x * x + 6 * x;
			};
			auto func0 = [](double x)->double {
				return (1 - x) * 3 * x * x;
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
		if (i == 3) {
			auto func1 = [](double x)->double {
				return 3*x*x;
			};
			auto func0 = [](double x)->double {
				return (x) * (x) * (x);
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
	}

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;

	//#pragma omp parallel for num_threads(8)
	for (int patch_id = 0; patch_id < patch_cnt; patch_id++)
	{
		if (patch_id % 100 == 0)    std::cout << "Schedule of stiffness matrix " << patch_id + 1 << ":" << patch_cnt << std::endl;
		
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;                                           //每个基对应一组在各个高斯点上的雅可比逆矩阵
		
		std::vector<double> m_det;                                                     //记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;


		int u_basis_cnt = 4;
		int v_basis_cnt = 4;
		int w_basis_cnt = 4;
		
		

		std::vector<double> u_Intervel{ 0, 1 };
		std::vector<double> v_Intervel{ 0, 1 };
		std::vector<double> w_Intervel{ 0, 1 };
		
		int elem_u_cnt = u_Intervel.size() - 1;
		int elem_v_cnt = v_Intervel.size() - 1;
		int elem_w_cnt = w_Intervel.size() - 1;



		/***************************构造刚度矩阵按照单元循环****************************************/
		m_det.clear(); m_InverseMatrix.clear();
		m_solids.JacobiMatrix_new(patch_id, para_gauss, para_gauss, para_gauss, m_det, m_InverseMatrix);

		std::vector<int> patch_global_base_index;

		for (int u_id = 0; u_id < elem_u_cnt; u_id++)
		{
			for (int v_id = 0; v_id < elem_v_cnt; v_id++)
			{
				for (int w_id = 0; w_id < elem_w_cnt; w_id++)
				{

					//std::cout << "Calc Derivative \n";
					// 第i个基第j个高斯点的导数
					std::vector<std::vector<double>> u_basis_derivative(local_basis_cnt, std::vector<double>(patch_gauss_cnt, 0));
					std::vector<std::vector<double>> v_basis_derivative(local_basis_cnt, std::vector<double>(patch_gauss_cnt, 0));
					std::vector<std::vector<double>> w_basis_derivative(local_basis_cnt, std::vector<double>(patch_gauss_cnt, 0));
					


					int cur_basis_id = 0;
					for (int u_basis_id = 0; u_basis_id < u_basis_cnt; u_basis_id++) {
						for (int v_basis_id = 0; v_basis_id < v_basis_cnt; v_basis_id++) {
							for (int w_basis_id = 0; w_basis_id < w_basis_cnt; w_basis_id++) {

								patch_global_base_index.push_back(m_solids.Control_Point[patch_id][w_basis_id][v_basis_id][u_basis_id]);

								int cur_gauss_id = 0;
								for (int w_gauss_id = 0; w_gauss_id < gauss_cnt; w_gauss_id++)
								{
									auto w_ders = pre_calc_ders[w_basis_id][w_gauss_id];
									
									for (int v_gauss_id = 0; v_gauss_id < gauss_cnt; v_gauss_id++)
									{
										auto v_ders = pre_calc_ders[v_basis_id][v_gauss_id];
										for (int u_gauss_id = 0; u_gauss_id < gauss_cnt; u_gauss_id++)
										{
											auto u_ders = pre_calc_ders[u_basis_id][u_gauss_id];

											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											u_basis_derivative[cur_basis_id][cur_gauss_id] = m_temp;
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											v_basis_derivative[cur_basis_id][cur_gauss_id] = m_temp;
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											w_basis_derivative[cur_basis_id][cur_gauss_id] = m_temp;
											cur_gauss_id++;
										}
									}
								}
								cur_basis_id++;
							}
						}
					}
					

					//std::cout << "Local Matrix \n";
					for (int i = 0; i < patch_global_base_index.size(); i++)
					{
						int ii = patch_global_base_index[i];
						for (int j = 0; j < patch_global_base_index.size(); j++)
						{
							int jj = patch_global_base_index[j];
							for (int k = 0; k < patch_gauss_cnt; k++)       // 64
							{
								int m_uptIndex = (k % (gauss_cnt * gauss_cnt)) % gauss_cnt;
								int m_vptIndex = (k % (gauss_cnt * gauss_cnt)) / gauss_cnt;
								int m_wptIndex = k / (gauss_cnt * gauss_cnt);
								int m_index;

								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;
								m_Rx = m_InverseMatrix[k][0][0] * u_basis_derivative[ii][k] + m_InverseMatrix[k][0][1] * v_basis_derivative[ii][k] + m_InverseMatrix[k][0][2] * w_basis_derivative[ii][k];
								m_Ry = m_InverseMatrix[k][1][0] * u_basis_derivative[ii][k] + m_InverseMatrix[k][1][1] * v_basis_derivative[ii][k] + m_InverseMatrix[k][1][2] * w_basis_derivative[ii][k];
								m_Rz = m_InverseMatrix[k][2][0] * u_basis_derivative[ii][k] + m_InverseMatrix[k][2][1] * v_basis_derivative[ii][k] + m_InverseMatrix[k][2][2] * w_basis_derivative[ii][k];
								m_Qx = m_InverseMatrix[k][0][0] * u_basis_derivative[jj][k] + m_InverseMatrix[k][0][1] * v_basis_derivative[jj][k] + m_InverseMatrix[k][0][2] * w_basis_derivative[jj][k];
								m_Qy = m_InverseMatrix[k][1][0] * u_basis_derivative[jj][k] + m_InverseMatrix[k][1][1] * v_basis_derivative[jj][k] + m_InverseMatrix[k][1][2] * w_basis_derivative[jj][k];
								m_Qz = m_InverseMatrix[k][2][0] * u_basis_derivative[jj][k] + m_InverseMatrix[k][2][1] * v_basis_derivative[jj][k] + m_InverseMatrix[k][2][2] * w_basis_derivative[jj][k];
								double temp = fabs(m_det[k]) 
									* gauss_w4[m_uptIndex] * gauss_w4[m_vptIndex] * gauss_w4[m_wptIndex] 
									* 1 
									* 1
									* 1 / 8;

								local_matrix[patch_id][ii * 3 + 0][jj * 3 + 0] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
								local_matrix[patch_id][ii * 3 + 0][jj * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
								local_matrix[patch_id][ii * 3 + 0][jj * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
								local_matrix[patch_id][ii * 3 + 1][jj * 3 + 0] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
								local_matrix[patch_id][ii * 3 + 1][jj * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
								local_matrix[patch_id][ii * 3 + 1][jj * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
								local_matrix[patch_id][ii * 3 + 2][jj * 3 + 0] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
								local_matrix[patch_id][ii * 3 + 2][jj * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
								local_matrix[patch_id][ii * 3 + 2][jj * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;
							}
						}
					}
				}
			}
		}
		//std::cout << "Global Matrix \n";
		// local_matrix.push_back(m_LocalMatrix);
		//矩阵中的非零元素全部存起来
		for (int i = 0; i < patch_global_base_index.size(); i++)
		{
			for (int j = 0; j < patch_global_base_index.size(); j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						//std::cout << i << " " << j << " " << k << " " << l << "\n";
						int local_id1 = i * 3 + k;
						int local_id2 = j * 3 + l;
						int global_id1 = m_solids.Local_to_Global[patch_id][patch_global_base_index[i]] * 3 + k;
						int global_id2 = m_solids.Local_to_Global[patch_id][patch_global_base_index[j]] * 3 + l;

						if (fabs(local_matrix[patch_id][local_id1][local_id2]) > 1.0e-12)
						{
							
							if (!InVecID(global_id2, m_colids[global_id1 ])) 
							{
								m_colids[global_id1].push_back(global_id2 );
								m_bigmatrix[global_id1].push_back(0);
								m_bigmatrix[global_id1].back() += local_matrix[patch_id][local_id1][local_id2];
							}
							else
							{
								int tmp_id = InVecIndex(global_id2 , m_colids[global_id1 ]);
								assert(tmp_id < m_bigmatrix[global_id1].size() && tmp_id >= 0);
								m_bigmatrix[global_id1][tmp_id] += local_matrix[patch_id][local_id1][local_id2];
							}
						}
					}
				}
			}
		}
		//std::cout << "Global Matrix Done once\n";
	}

	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;

	//std::cout << "bdybaseid.size() = " << bdybaseid.size() << std::endl;
	std::unordered_set<int> fix_base;

	for (int i = 0; i < bdybaseid.size(); ++i) {
		
		int m_globalid = bdybaseid[i].global_id;
		int m_bdytype = bdybaseid[i].bdy_type;
		if (m_bdytype != 1) continue;
		for (int l = 0; l < 3; ++l) {
			fix_base.insert(m_globalid * 3 + l);
		}
	}

	for (int j = 0; j <global_dof_cnt; j++)
	{
		//std::cout << j << std::endl;
		for (int k = 0; k < m_colids[j].size(); k++)
		{
			if (fix_base.count(m_colids[j][k]))
			{
				m_bigmatrix[j][k] = 0;
			}
		}
	}

	for (int i = 0; i < bdybaseid.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		int m_bdytype = bdybaseid[i].bdy_type;
		//std::cout << "bdybaseid[i].bdy_type = " << bdybaseid[i].bdy_type << std::endl;
		// 
		//for (int j = 0; j < m_GlobalSize * 3; j++)
		//{
		//	for (int k = 0; k < m_colids[j].size(); k++)
		//	{
		//		//if (bdybase_count < subnum)
		//		//{
		//		//	if ((m_colids[j][k]) == m_globalid * 3 + forceorder)
		//		//	{
		//		//		temp_bdyquick.push_back(std::pair<int, int>(j, k));
		//		//		m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
		//		//		m_bigmatrix[j][k] = 0;
		//		//		break;
		//		//	}
		//		//}
		//		//else
		//		//{
		//		//	for (int l = 0; l < 3; l++)
		//		//	{
		//		//		if ((m_colids[j][k]) == m_globalid * 3 + l)
		//		//		{
		//		//			temp_bdyquick.push_back(std::pair<int, int>(j, k));
		//		//			m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
		//		//			m_bigmatrix[j][k] = 0;
		//		//			break;
		//		//		}
		//		//	}
		//		//}
		//		if (m_bdytype == 1) {
		//			// 支撑
		//			for (int l = 0; l < 3; l++)
		//			{
		//				if ((m_colids[j][k]) == m_globalid * 3 + l)
		//				{
		//					m_bigmatrix[j][k] = 0;
		//					break;
		//				}
		//			}

		//		}
		//	}
		//}

		if (m_bdytype == 2)
		{
			/*for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
			{
				if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
				{
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
				}
				else
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
			}*/
			// 施加力
			for (int dir = 0; dir < 3; ++dir) {
				m_bigrhs[m_globalid * 3 + dir] = m_bdycoef[m_solidid][m_baseid * 3 + dir];
			}
		}
		else if(m_bdytype == 1)
		{
			/*for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}*/
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + m)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = 0;
			}
		}
		bdybase_count++;
	}

	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	std::vector<T> coeffs;
	for (int i = 0; i < global_dof_cnt; ++i) {
		for (int k = 0; k < m_colids[i].size(); k++)
		{
			if (fabs(m_bigmatrix[i][k]) > 1.0e-12) 
				coeffs.push_back(T(i, m_colids[i][k], m_bigmatrix[i][k]));
		}
	}
	{
		Timer timer;
		std::cout << "prepare solving...\n";
		Eigen::setNbThreads(8);  // 线程池大小
		Eigen::initParallel();   // 初始化多线程


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

		u = std::vector<double>(eigen_u.data(), eigen_u.data() + eigen_u.size());

		std::cout << "end solving...\n";
	}
	

	/*std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(0);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize * 3; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-9)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num);
	}*/
	Matrix().swap(m_bigmatrix);
}


void FunctionBSpline::Single_Resolution_TopOpt_Init_Step(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs,std::vector<std::vector<int>>& m_colids, std::vector<double>& uu)
{
	const int global_basis_cnt = m_solids.Global_size;          // 全局基函数数量
	const int global_dof_cnt = m_solids.Global_size * 3;        // 全局自由度数量
	const int patch_cnt = m_solids.Local_V3f.size();            // Bezier patch数量
	const int local_basis_cnt = 64;
	const int local_dof_cnt = local_basis_cnt * 3;              // 局部自由度数量
	std::cout << "global_dof_cnt = " << global_dof_cnt << "\n";


	// 初始化右端项
	m_bigrhs.resize(global_dof_cnt, 0);
	m_colids.resize(global_dof_cnt);

	// 初始化刚度矩阵
	Matrix m_bigmatrix;
	m_bigmatrix.resize(global_dof_cnt);
	local_matrix.resize(patch_cnt);
	for (int i = 0; i < patch_cnt; i++)
	{
		local_matrix[i].resize(local_dof_cnt, std::vector<double>(local_dof_cnt, 0));
	}

	// 初始化高斯点值与权重
	std::vector<double> gauss_x4 = { -0.8611363116, -0.3399810436, 0.3399810436,  0.8611363116 };
	std::vector<double> gauss_w4 = { 0.3478548451 ,0.6521451549, 0.6521451549, 0.3478548451 };
	int gauss_cnt = 4;
	int patch_gauss_cnt = gauss_cnt * gauss_cnt * gauss_cnt;

	std::vector<double> para_gauss{ 0.0694318442,0.3300094782,0.6699905218,0.9305681558 };

	std::vector<std::vector<std::vector<double>>> pre_calc_ders;
	pre_calc_ders.resize(4);
	for (int i = 0; i < 4; ++i) {
		pre_calc_ders[i].resize(gauss_cnt);
		if (i == 0) {
			auto func1 = [](double x)->double {
				return -3 * x * x + 6 * x - 3;
			};
			auto func0 = [](double x)->double {
				return (1 - x) * (1 - x) * (1 - x);
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
		if (i == 1) {
			auto func1 = [](double x)->double {
				return 9 * x * x - 12 * x + 3;
			};
			auto func0 = [](double x)->double {
				return 3 * x * x * x - 6 * x * x + 3 * x;
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
		if (i == 2) {
			auto func1 = [](double x)->double {
				return -9 * x * x + 6 * x;
			};
			auto func0 = [](double x)->double {
				return (1 - x) * 3 * x * x;
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
		if (i == 3) {
			auto func1 = [](double x)->double {
				return 3 * x * x;
			};
			auto func0 = [](double x)->double {
				return (x) * (x) * (x);
			};
			for (int gauss_id = 0; gauss_id < gauss_cnt; ++gauss_id) {
				pre_calc_ders[i][gauss_id].resize(2);
				pre_calc_ders[i][gauss_id][0] = func0(para_gauss[gauss_id]);
				pre_calc_ders[i][gauss_id][1] = func1(para_gauss[gauss_id]);
			}
		}
	}

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;

	//#pragma omp parallel for num_threads(8)
	for (int patch_id = 0; patch_id < patch_cnt; patch_id++)
	{
		if (patch_id % 100 == 0)    std::cout << "Schedule of stiffness matrix " << patch_id + 1 << ":" << patch_cnt << std::endl;

		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;                                           //每个基对应一组在各个高斯点上的雅可比逆矩阵

		std::vector<double> m_det;                                                     //记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;


		int u_basis_cnt = 4;
		int v_basis_cnt = 4;
		int w_basis_cnt = 4;



		std::vector<double> u_Intervel{ 0, 1 };
		std::vector<double> v_Intervel{ 0, 1 };
		std::vector<double> w_Intervel{ 0, 1 };

		int elem_u_cnt = u_Intervel.size() - 1;
		int elem_v_cnt = v_Intervel.size() - 1;
		int elem_w_cnt = w_Intervel.size() - 1;



		/***************************构造刚度矩阵按照单元循环****************************************/
		m_det.clear(); m_InverseMatrix.clear();
		m_solids.JacobiMatrix_new(patch_id, para_gauss, para_gauss, para_gauss, m_det, m_InverseMatrix);

		std::vector<int> patch_global_base_index;

		for (int u_id = 0; u_id < elem_u_cnt; u_id++)
		{
			for (int v_id = 0; v_id < elem_v_cnt; v_id++)
			{
				for (int w_id = 0; w_id < elem_w_cnt; w_id++)
				{

					//std::cout << "Calc Derivative \n";
					// 第i个基第j个高斯点的导数
					std::vector<std::vector<double>> u_basis_derivative(local_basis_cnt, std::vector<double>(patch_gauss_cnt, 0));
					std::vector<std::vector<double>> v_basis_derivative(local_basis_cnt, std::vector<double>(patch_gauss_cnt, 0));
					std::vector<std::vector<double>> w_basis_derivative(local_basis_cnt, std::vector<double>(patch_gauss_cnt, 0));



					int cur_basis_id = 0;
					for (int u_basis_id = 0; u_basis_id < u_basis_cnt; u_basis_id++) {
						for (int v_basis_id = 0; v_basis_id < v_basis_cnt; v_basis_id++) {
							for (int w_basis_id = 0; w_basis_id < w_basis_cnt; w_basis_id++) {

								patch_global_base_index.push_back(m_solids.Control_Point[patch_id][w_basis_id][v_basis_id][u_basis_id]);

								int cur_gauss_id = 0;
								for (int w_gauss_id = 0; w_gauss_id < gauss_cnt; w_gauss_id++)
								{
									auto w_ders = pre_calc_ders[w_basis_id][w_gauss_id];

									for (int v_gauss_id = 0; v_gauss_id < gauss_cnt; v_gauss_id++)
									{
										auto v_ders = pre_calc_ders[v_basis_id][v_gauss_id];
										for (int u_gauss_id = 0; u_gauss_id < gauss_cnt; u_gauss_id++)
										{
											auto u_ders = pre_calc_ders[u_basis_id][u_gauss_id];

											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											u_basis_derivative[cur_basis_id][cur_gauss_id] = m_temp;
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											v_basis_derivative[cur_basis_id][cur_gauss_id] = m_temp;
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											w_basis_derivative[cur_basis_id][cur_gauss_id] = m_temp;
											cur_gauss_id++;
										}
									}
								}
								cur_basis_id++;
							}
						}
					}


					//std::cout << "Local Matrix \n";
					for (int i = 0; i < patch_global_base_index.size(); i++)
					{
						int ii = patch_global_base_index[i];
						for (int j = 0; j < patch_global_base_index.size(); j++)
						{
							int jj = patch_global_base_index[j];
							for (int k = 0; k < patch_gauss_cnt; k++)       // 64
							{
								int m_uptIndex = (k % (gauss_cnt * gauss_cnt)) % gauss_cnt;
								int m_vptIndex = (k % (gauss_cnt * gauss_cnt)) / gauss_cnt;
								int m_wptIndex = k / (gauss_cnt * gauss_cnt);
								int m_index;

								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;
								m_Rx = m_InverseMatrix[k][0][0] * u_basis_derivative[ii][k] + m_InverseMatrix[k][0][1] * v_basis_derivative[ii][k] + m_InverseMatrix[k][0][2] * w_basis_derivative[ii][k];
								m_Ry = m_InverseMatrix[k][1][0] * u_basis_derivative[ii][k] + m_InverseMatrix[k][1][1] * v_basis_derivative[ii][k] + m_InverseMatrix[k][1][2] * w_basis_derivative[ii][k];
								m_Rz = m_InverseMatrix[k][2][0] * u_basis_derivative[ii][k] + m_InverseMatrix[k][2][1] * v_basis_derivative[ii][k] + m_InverseMatrix[k][2][2] * w_basis_derivative[ii][k];
								m_Qx = m_InverseMatrix[k][0][0] * u_basis_derivative[jj][k] + m_InverseMatrix[k][0][1] * v_basis_derivative[jj][k] + m_InverseMatrix[k][0][2] * w_basis_derivative[jj][k];
								m_Qy = m_InverseMatrix[k][1][0] * u_basis_derivative[jj][k] + m_InverseMatrix[k][1][1] * v_basis_derivative[jj][k] + m_InverseMatrix[k][1][2] * w_basis_derivative[jj][k];
								m_Qz = m_InverseMatrix[k][2][0] * u_basis_derivative[jj][k] + m_InverseMatrix[k][2][1] * v_basis_derivative[jj][k] + m_InverseMatrix[k][2][2] * w_basis_derivative[jj][k];
								double temp = fabs(m_det[k])
									* gauss_w4[m_uptIndex] * gauss_w4[m_vptIndex] * gauss_w4[m_wptIndex]
									* 1
									* 1
									* 1 / 8;

								local_matrix[patch_id][ii * 3 + 0][jj * 3 + 0] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
								local_matrix[patch_id][ii * 3 + 0][jj * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
								local_matrix[patch_id][ii * 3 + 0][jj * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
								local_matrix[patch_id][ii * 3 + 1][jj * 3 + 0] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
								local_matrix[patch_id][ii * 3 + 1][jj * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
								local_matrix[patch_id][ii * 3 + 1][jj * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
								local_matrix[patch_id][ii * 3 + 2][jj * 3 + 0] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
								local_matrix[patch_id][ii * 3 + 2][jj * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
								local_matrix[patch_id][ii * 3 + 2][jj * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;
							}
						}
					}
				}
			}
		}
		//std::cout << "Global Matrix \n";
		// local_matrix.push_back(m_LocalMatrix);
		//矩阵中的非零元素全部存起来
		for (int i = 0; i < patch_global_base_index.size(); i++)
		{
			for (int j = 0; j < patch_global_base_index.size(); j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						//std::cout << i << " " << j << " " << k << " " << l << "\n";
						int local_id1 = i * 3 + k;
						int local_id2 = j * 3 + l;
						int global_id1 = m_solids.Local_to_Global[patch_id][patch_global_base_index[i]] * 3 + k;
						int global_id2 = m_solids.Local_to_Global[patch_id][patch_global_base_index[j]] * 3 + l;

						if (fabs(local_matrix[patch_id][local_id1][local_id2]) > 1.0e-12)
						{

							if (!InVecID(global_id2, m_colids[global_id1]))
							{
								m_colids[global_id1].push_back(global_id2);
								m_bigmatrix[global_id1].push_back(0);
								m_bigmatrix[global_id1].back() += local_matrix[patch_id][local_id1][local_id2];
							}
							else
							{
								int tmp_id = InVecIndex(global_id2, m_colids[global_id1]);
								assert(tmp_id < m_bigmatrix[global_id1].size() && tmp_id >= 0);
								m_bigmatrix[global_id1][tmp_id] += local_matrix[patch_id][local_id1][local_id2];
							}
						}
					}
				}
			}
		}
		//std::cout << "Global Matrix Done once\n";
	}

	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;

	//std::cout << "bdybaseid.size() = " << bdybaseid.size() << std::endl;
	std::unordered_set<int> fix_base;

	for (int i = 0; i < bdybaseid.size(); ++i) {

		int m_globalid = bdybaseid[i].global_id;
		int m_bdytype = bdybaseid[i].bdy_type;
		if (m_bdytype != 1) continue;
		for (int l = 0; l < 3; ++l) {
			fix_base.insert(m_globalid * 3 + l);
		}
	}

	for (int j = 0; j < global_dof_cnt; j++)
	{
		//std::cout << j << std::endl;
		for (int k = 0; k < m_colids[j].size(); k++)
		{
			if (fix_base.count(m_colids[j][k]))
			{
				m_bigmatrix[j][k] = 0;
			}
		}
	}

	for (int i = 0; i < bdybaseid.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		int m_bdytype = bdybaseid[i].bdy_type;
		//std::cout << "bdybaseid[i].bdy_type = " << bdybaseid[i].bdy_type << std::endl;
		// 
		//for (int j = 0; j < m_GlobalSize * 3; j++)
		//{
		//	for (int k = 0; k < m_colids[j].size(); k++)
		//	{
		//		//if (bdybase_count < subnum)
		//		//{
		//		//	if ((m_colids[j][k]) == m_globalid * 3 + forceorder)
		//		//	{
		//		//		temp_bdyquick.push_back(std::pair<int, int>(j, k));
		//		//		m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
		//		//		m_bigmatrix[j][k] = 0;
		//		//		break;
		//		//	}
		//		//}
		//		//else
		//		//{
		//		//	for (int l = 0; l < 3; l++)
		//		//	{
		//		//		if ((m_colids[j][k]) == m_globalid * 3 + l)
		//		//		{
		//		//			temp_bdyquick.push_back(std::pair<int, int>(j, k));
		//		//			m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
		//		//			m_bigmatrix[j][k] = 0;
		//		//			break;
		//		//		}
		//		//	}
		//		//}
		//		if (m_bdytype == 1) {
		//			// 支撑
		//			for (int l = 0; l < 3; l++)
		//			{
		//				if ((m_colids[j][k]) == m_globalid * 3 + l)
		//				{
		//					m_bigmatrix[j][k] = 0;
		//					break;
		//				}
		//			}

		//		}
		//	}
		//}

		if (m_bdytype == 2)
		{
			/*for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
			{
				if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
				{
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
				}
				else
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
			}*/
			// 施加力
			for (int dir = 0; dir < 3; ++dir) {
				m_bigrhs[m_globalid * 3 + dir] = m_bdycoef[m_solidid][m_baseid * 3 + dir];
			}
		}
		else if (m_bdytype == 1)
		{
			/*for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}*/
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + m)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = 0;
			}
		}
		bdybase_count++;
	}

	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	std::vector<T> coeffs;
	for (int i = 0; i < global_dof_cnt; ++i) {
		for (int k = 0; k < m_colids[i].size(); k++)
		{
			if (fabs(m_bigmatrix[i][k]) > 1.0e-12)
				coeffs.push_back(T(i, m_colids[i][k], m_bigmatrix[i][k]));
		}
	}
	{
		Timer timer;
		std::cout << "prepare solving...\n";
		Eigen::setNbThreads(8);  // 线程池大小
		Eigen::initParallel();   // 初始化多线程


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

		uu = std::vector<double>(eigen_u.data(), eigen_u.data() + eigen_u.size());

		std::cout << "end solving...\n";
	}

	Matrix().swap(m_bigmatrix);
}


void FunctionBSpline::Single_Resolution_TopOpt_Loop_Step(std::vector<Matrix>& local_matrixs, std::vector<double>& m_bigrhs,
	std::vector<std::vector<int>>& m_colids, std::vector<BdyBase>& bdybaseid, std::vector<double>&u)
{
	int m_GlobalSize = m_solids.Global_size;        // total control size   < 64 * patches


	Matrix m_bigmatrix;
	std::cout << "m_solids.Global_size() = " << m_solids.Global_size << std::endl;
	std::cout << "m_colids.size() = " << m_colids.size() << std::endl;
	m_bigmatrix.resize(m_GlobalSize * 3);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			for (int k = 0; k < 3; k++)
				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize(m_colids[m_solids.Local_to_Global[i][j] * 3 + k].size(), 0);
		}
	}
	std::cout << "m_bigmatrix.size() = " << m_bigmatrix.size() << std::endl;

	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{

		if (i % 100 == 0)    std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		
		int m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		int m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		int m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		int m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		

		double rho = (umin + (1 - umin) * pow(material_distribution[i], p));
		//std::cout << "rho = " << rho << std::endl;

		//矩阵中的非零元素全部存起来
		for (int ii = 0; ii < m_MatrixSize; ii++)
		{
			for (int j = 0; j < m_MatrixSize; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						if (fabs(local_matrixs[i][ii * 3 + k][j * 3 + l]) > 1.0e-12)
						{
							int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
							if (tmp_id >= m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size()) continue;
							
							m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrixs[i][ii * 3 + k][j * 3 + l] * rho;

						}
					}
				}
			}
		}	
	}



	int bdybase_count = 0;

	//std::cout << "bdybaseid.size() = " << bdybaseid.size() << std::endl;
	std::unordered_set<int> fix_base;

	for (int i = 0; i < bdybaseid.size(); ++i) {

		int m_globalid = bdybaseid[i].global_id;
		int m_bdytype = bdybaseid[i].bdy_type;
		if (m_bdytype != 1) continue;
		for (int l = 0; l < 3; ++l) {
			fix_base.insert(m_globalid * 3 + l);
		}
	}

	for (int j = 0; j < m_GlobalSize * 3; j++)
	{
		for (int k = 0; k < m_colids[j].size(); k++)
		{
			if (fix_base.count(m_colids[j][k]))
			{
				m_bigmatrix[j][k] = 0;
			}
		}
	}

	for (int i = 0; i < bdybaseid.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		int m_bdytype = bdybaseid[i].bdy_type;
		

		if (m_bdytype == 2)
		{
		
			// 施加力

		}
		else if (m_bdytype == 1)
		{
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + m)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = 0;
			}
		}
		bdybase_count++;
	}

	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	std::vector<T> coeffs;
	for (int i = 0; i < m_GlobalSize * 3; ++i) {
		for (int k = 0; k < m_colids[i].size(); k++)
		{
			if (fabs(m_bigmatrix[i][k]) > 1.0e-12)
				coeffs.push_back(T(i, m_colids[i][k], m_bigmatrix[i][k]));
		}
	}
	{
		Timer timer;
		std::cout << "prepare solving...\n";
		Eigen::setNbThreads(8);  // 线程池大小
		Eigen::initParallel();   // 初始化多线程


		std::cout << "Current threads: " << Eigen::nbThreads() << std::endl;
		SpMat A(m_GlobalSize * 3, m_GlobalSize * 3);
		A.setFromTriplets(coeffs.begin(), coeffs.end());
		Eigen::VectorXd b = Eigen::VectorXd::Map(m_bigrhs.data(), m_bigrhs.size());

		Eigen::PardisoLU<SpMat> solver;
		std::cout << "start solving...\n";
		solver.compute(A);
		//solver.compute(A.selfadjointView<Eigen::Lower>());
		/*solver.analyzePattern(A);
		solver.factorize(A)*/;

		Eigen::VectorXd eigen_u = solver.solve(b);

		u = std::vector<double>(eigen_u.data(), eigen_u.data() + eigen_u.size());

		std::cout << "end solving...\n";
	}
	
	Matrix().swap(m_bigmatrix);
}


void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_quick(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
	std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
	std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick, std::vector<double>& u)
{
	int m_GlobalSize = m_solids.Global_size;
	m_bigrhs.clear();
	m_bigrhs.resize(m_GlobalSize * 3, 0);

	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize * 3);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			for (int k = 0; k < 3; k++)
				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2) * martix_size, 0);
		}
	}

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		int m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		int m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		int m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		int m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
		//矩阵中的非零元素全部存起来
		for (int ii = 0; ii < m_MatrixSize; ii++)
		{
			for (int j = 0; j < m_MatrixSize; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						if (fabs(local_matrix[i][ii * 3 + k][j * 3 + l]) > 1.0e-7)
						{
							if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
							{
								m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
								int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrix[i][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
							else
							{
								int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrix[i][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < local_bigrhs.size(); i++)
	{
		for (int k = 0; k < 3; k++)
			for (int j = 0; j < local_bigrhs[i].size() / 3; j++)
				m_bigrhs[m_solids.Local_to_Global[i][j] * 3 + k] += local_bigrhs[i][j * 3 + k] * (umin + (1 - umin) * pow(material_distribution[i], p));
	}
	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		if (i % 1000 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		std::vector<std::pair<int, int>> temp_bdyquick = m_bdyquick[i];
		if (bdybase_count < subnum)
		{
			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			{
				int j = (*temp_bdyquick_it).first;
				int k = (*temp_bdyquick_it).second;
				m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
				m_bigmatrix[j][k] = 0;
			}
		}
		else
		{
			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			{
				int j = (*temp_bdyquick_it).first;
				int k = (*temp_bdyquick_it).second;
				for (int l = 0; l < 3; l++)
				{
					if ((m_colids[j][k]) == m_globalid * 3 + l)
					{
						m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
						break;
					}
				}
				m_bigmatrix[j][k] = 0;
			}
		}
		// if (bdybase_count < 4)
		// {
			// for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			// {
				// int j = (*temp_bdyquick_it).first;
				// int k = (*temp_bdyquick_it).second;
				// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + 1] * m_bigmatrix[j][k];
				// m_bigmatrix[j][k] = 0;
			// }
		// }
		// else if (bdybase_count < 8)
		// {
			// for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			// {
				// int j = (*temp_bdyquick_it).first;
				// int k = (*temp_bdyquick_it).second;
				// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3] * m_bigmatrix[j][k];
				// m_bigmatrix[j][k] = 0;
			// }
		// }
		// else
		// {
			// for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			// {
				// int j = (*temp_bdyquick_it).first;
				// int k = (*temp_bdyquick_it).second;
				// for (int l = 0; l < 3; l++)
				// {
					// if ((m_colids[j][k]) == m_globalid * 3 + l)
					// {
						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
						// break;
					// }
				// }
				// m_bigmatrix[j][k] = 0;
			// }
		// }
		if (bdybase_count < subnum)
		{
			for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
			{
				if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
				{
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
				}
				else
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
			}
			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
		}
		else
		{
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}
		}
		// if (bdybase_count < 4)
		// {
			// for (int k = 0; k < m_colids[m_globalid * 3 + 1].size(); k++)
			// {
				// if (m_colids[m_globalid * 3 + 1][k] == m_globalid * 3 + 1)
				// {
					// m_bigmatrix[m_globalid * 3 + 1][k] = 1;
				// }
				// else
					// m_bigmatrix[m_globalid * 3 + 1][k] = 0;
			// }
			// m_bigrhs[m_globalid * 3 + 1] = m_bdycoef[m_solidid][m_baseid * 3 + 1];
		// }
		// else if (bdybase_count < 8)
		// {
			// for (int k = 0; k < m_colids[m_globalid * 3].size(); k++)
			// {
				// if (m_colids[m_globalid * 3][k] == m_globalid * 3)
				// {
					// m_bigmatrix[m_globalid * 3][k] = 1;
				// }
				// else
					// m_bigmatrix[m_globalid * 3][k] = 0;
			// }
			// m_bigrhs[m_globalid * 3] = m_bdycoef[m_solidid][m_baseid * 3];
		// }
		// else
		// {
			// for (int m = 0; m < 3; m++)
			// {
				// for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				// {
					// if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						// || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					// {
						// m_bigmatrix[m_globalid * 3 + m][k] = 1;
					// }
					// else
						// m_bigmatrix[m_globalid * 3 + m][k] = 0;
				// }
				// m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			// }
		// }
		bdybase_count++;
	}
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(1);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize * 3; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first + 1);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num + 1);
	}



	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	std::vector<T> coeffs;
	for (int i = 0; i < m_GlobalSize * 3; ++i) {
		for (int k = 0; k < m_colids[i].size(); k++)
		{
			if (fabs(m_bigmatrix[i][k]) > 1.0e-12)
				coeffs.push_back(T(i, m_colids[i][k], m_bigmatrix[i][k]));
		}
	}
	{
		Timer timer;
		std::cout << "prepare solving...\n";
		Eigen::setNbThreads(8);  // 线程池大小
		Eigen::initParallel();   // 初始化多线程


		std::cout << "Current threads: " << Eigen::nbThreads() << std::endl;
		SpMat A(m_GlobalSize * 3, m_GlobalSize * 3);
		A.setFromTriplets(coeffs.begin(), coeffs.end());
		Eigen::VectorXd b = Eigen::VectorXd::Map(m_bigrhs.data(), m_bigrhs.size());

		Eigen::PardisoLU<SpMat> solver;
		std::cout << "start solving...\n";
		solver.compute(A);
		//solver.compute(A.selfadjointView<Eigen::Lower>());
		/*solver.analyzePattern(A);
		solver.factorize(A)*/;

		Eigen::VectorXd eigen_u = solver.solve(b);

		u = std::vector<double>(eigen_u.data(), eigen_u.data() + eigen_u.size());

		std::cout << "end solving...\n";
	}

	Matrix().swap(m_bigmatrix);
}

void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids, std::vector<int>& r_colids,
	std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick,
	std::string name)
{
	int m_GlobalSize = m_solids.Global_size;
	/*---------------存储每一行中的非零元下标------------*/
	m_colids.resize(m_GlobalSize * 3);

	m_bigrhs.resize(m_GlobalSize * 3);
	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize * 3);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			for (int k = 0; k < 3; k++)
				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2) * martix_size, 0);
		}
	}
	local_bigrhs.resize(m_solids.Local_V3f.size());
	local_matrix.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_bigrhs[i].resize(size * 3, 0);
		local_matrix[i].resize(size * 3, std::vector<double>(size * 3, 0));
	}
	// local_matrix.reserve(m_solids.Local_V3f.size());

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);

		// Matrix m_LocalMatrix;
		// m_LocalMatrix.resize(m_MatrixSize * 3);
		// for (int k = 0; k < m_MatrixSize * 3; k++)
		// {
			// m_LocalMatrix[k].resize(m_MatrixSize * 3, 0);
		// }

		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue1 = 0, m_tempvalue2 = 0, m_tempvalue3 = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue1 += u_ders[0] * v_ders[0] * w_ders[0] * f1(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
											m_tempvalue2 += u_ders[0] * v_ders[0] * w_ders[0] * f2(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
											m_tempvalue3 += u_ders[0] * v_ders[0] * w_ders[0] * f3(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								double temp = (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_bigrhs[i][3 * m_BaseIndex] += m_tempvalue1 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue1 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
								local_bigrhs[i][3 * m_BaseIndex + 1] += m_tempvalue2 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 1] += (m_tempvalue2 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
								local_bigrhs[i][3 * m_BaseIndex + 2] += m_tempvalue3 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 2] += (m_tempvalue3 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double temp = fabs(m_det[k]) * m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1]
									- m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_matrix[i][m_baseindex[ii] * 3][m_baseindex[j] * 3] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
								local_matrix[i][m_baseindex[ii] * 3][m_baseindex[j] * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
								local_matrix[i][m_baseindex[ii] * 3][m_baseindex[j] * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
								local_matrix[i][m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;
							}
						}
					}
				}
			}
		}
		// local_matrix.push_back(m_LocalMatrix);
		//矩阵中的非零元素全部存起来
		for (int ii = 0; ii < m_MatrixSize; ii++)
		{
			for (int j = 0; j < m_MatrixSize; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						if (fabs(local_matrix[i][ii * 3 + k][j * 3 + l]) > 1.0e-7)
						{
							if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
							{
								m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
								int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrix[i][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
							else
							{
								int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrix[i][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
						}
					}
				}
			}
		}
	}
	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		std::vector<std::pair<int, int>> temp_bdyquick;
		if (i % 100 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		for (int j = 0; j < m_GlobalSize * 3; j++)
		{
			for (int k = 0; k < m_colids[j].size(); k++)
			{
				if (bdybase_count < subnum)
				{
					if ((m_colids[j][k]) == m_globalid * 3 + forceorder)
					{
						temp_bdyquick.push_back(std::pair<int, int>(j, k));
						m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
						m_bigmatrix[j][k] = 0;
						break;
					}
				}
				else
				{
					for (int l = 0; l < 3; l++)
					{
						if ((m_colids[j][k]) == m_globalid * 3 + l)
						{
							temp_bdyquick.push_back(std::pair<int, int>(j, k));
							m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
							m_bigmatrix[j][k] = 0;
							break;
						}
					}
				}
			}
			// for (int k = 0; k < m_colids[j].size(); k++)
			// {
				// if (bdybase_count < 4)
				// {
					// if ((m_colids[j][k]) == m_globalid * 3 + 1)
					// {
						// temp_bdyquick.push_back(std::pair<int, int>(j, k));
						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + 1] * m_bigmatrix[j][k];
						// m_bigmatrix[j][k] = 0;
						// break;
					// }
				// }
				// else if (bdybase_count < 8)
				// {
					// if ((m_colids[j][k]) == m_globalid * 3)
					// {
						// temp_bdyquick.push_back(std::pair<int, int>(j, k));
						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3] * m_bigmatrix[j][k];
						// m_bigmatrix[j][k] = 0;
						// break;
					// }
				// }
				// else
				// {
					// for (int l = 0; l < 3; l++)
					// {
						// if ((m_colids[j][k]) == m_globalid * 3 + l)
						// {
							// temp_bdyquick.push_back(std::pair<int, int>(j, k));
							// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
							// m_bigmatrix[j][k] = 0;
							// break;
						// }
					// }
				// }
			// }
		}
		m_bdyquick.push_back(temp_bdyquick);
		if (bdybase_count < subnum)
		{
			for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
			{
				if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
				{
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
				}
				else
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
			}
			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
		}
		else
		{
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}
		}
		// if (bdybase_count < 4)
		// {
			// for (int k = 0; k < m_colids[m_globalid * 3 + 1].size(); k++)
			// {
				// if (m_colids[m_globalid * 3 + 1][k] == m_globalid * 3 + 1)
				// {
					// m_bigmatrix[m_globalid * 3 + 1][k] = 1;
				// }
				// else
					// m_bigmatrix[m_globalid * 3 + 1][k] = 0;
			// }
			// m_bigrhs[m_globalid * 3 + 1] = m_bdycoef[m_solidid][m_baseid * 3 + 1];
		// }
		// else if (bdybase_count < 8)
		// {
			// for (int k = 0; k < m_colids[m_globalid * 3].size(); k++)
			// {
				// if (m_colids[m_globalid * 3][k] == m_globalid * 3)
				// {
					// m_bigmatrix[m_globalid * 3 ][k] = 1;
				// }
				// else
					// m_bigmatrix[m_globalid * 3][k] = 0;
			// }
			// m_bigrhs[m_globalid * 3] = m_bdycoef[m_solidid][m_baseid * 3];
		// }
		// else
		// {
			// for (int m = 0; m < 3; m++)
			// {
				// for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				// {
					// if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						// || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					// {
						// m_bigmatrix[m_globalid * 3 + m][k] = 1;
					// }
					// else
						// m_bigmatrix[m_globalid * 3 + m][k] = 0;
				// }
				// m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			// }
		// }
		bdybase_count++;
	}
	FileHandler filenewxx = FileHandler();
	filenewxx.writebdy(name, m_bdyquick);
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(1);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize * 3; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first + 1);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num + 1);
	}
	Matrix().swap(m_bigmatrix);
	std::vector<double>().swap(m_Weight);
	std::vector<double>().swap(m_IntePt);
}

void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner1(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids, std::vector<int>& r_colids,
	std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick)
{
	int m_GlobalSize = m_solids.Global_size;
	/*---------------存储每一行中的非零元下标------------*/
	m_colids.resize(m_GlobalSize * 3);

	m_bigrhs.resize(m_GlobalSize * 3);
	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize * 3);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			for (int k = 0; k < 3; k++)
				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2) * martix_size, 0);
		}
	}
	local_bigrhs.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		int size = m_solids.Local_V3f[i].size();
		local_bigrhs[i].resize(size * 3, 0);
	}
	local_matrix.reserve(m_solids.Local_V3f.size());

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	std::vector<int> rowids, colids;
	std::vector<double> elems;

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		std::vector<std::vector<double>> m_jacobiInverseMatrix;
		std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
		double u_ders[2], v_ders[2], w_ders[2];
		std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
		double m_jacobidet;
		int m_BaseIndex, m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
		std::vector<double> m_baseIntegralforboundary;//每个基在区域的积分，对于矩阵奇异的情形要修改系数矩阵
		std::vector<Matrix> m_LocalMatrixes;

		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
		m_LocalMatrixes.resize((m_uIntervel.size() - 1) * (m_vIntervel.size() - 1) * (m_wIntervel.size() - 1));

		Matrix m_LocalMatrix;
		m_LocalMatrix.resize(m_MatrixSize * 3);
		for (int k = 0; k < m_MatrixSize * 3; k++)
		{
			m_LocalMatrix[k].resize(m_MatrixSize * 3, 0);
		}

		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue1 = 0, m_tempvalue2 = 0, m_tempvalue3 = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue1 += u_ders[0] * v_ders[0] * w_ders[0] * f1(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
											m_tempvalue2 += u_ders[0] * v_ders[0] * w_ders[0] * f2(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
											m_tempvalue3 += u_ders[0] * v_ders[0] * w_ders[0] * f3(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								double temp = (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								local_bigrhs[i][3 * m_BaseIndex] += m_tempvalue1 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue1 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
								local_bigrhs[i][3 * m_BaseIndex + 1] += m_tempvalue2 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 1] += (m_tempvalue2 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
								local_bigrhs[i][3 * m_BaseIndex + 2] += m_tempvalue3 * temp;
								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 2] += (m_tempvalue3 * temp) * (umin + (1 - umin) * pow(material_distribution[i], p));
							}

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double temp = fabs(m_det[k]) * m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1]
									- m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;
							}
						}
					}
				}
			}
		}
		local_matrix.push_back(m_LocalMatrix);
		//矩阵中的非零元素全部存起来
		for (int ii = 0; ii < m_MatrixSize; ii++)
		{
			for (int j = 0; j < m_MatrixSize; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						if (fabs(m_LocalMatrix[ii * 3 + k][j * 3 + l]) > 1.0e-7)
						{
							if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
							{
								m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
								int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += m_LocalMatrix[ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
							else
							{
								int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += m_LocalMatrix[ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution[i], p));
							}
						}
					}
				}
			}
		}
	}
	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		if (i % 1000 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		std::vector<std::pair<int, int>> temp_bdyquick = m_bdyquick[i];
		if (bdybase_count < subnum)
		{
			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			{
				int j = (*temp_bdyquick_it).first;
				int k = (*temp_bdyquick_it).second;
				m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
				m_bigmatrix[j][k] = 0;
			}
		}
		else
		{
			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			{
				int j = (*temp_bdyquick_it).first;
				int k = (*temp_bdyquick_it).second;
				for (int l = 0; l < 3; l++)
				{
					if ((m_colids[j][k]) == m_globalid * 3 + l)
					{
						m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
						break;
					}
				}
				m_bigmatrix[j][k] = 0;
			}
		}
		/*if (bdybase_count == 0)
		{
			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			{
				int j = (*temp_bdyquick_it).first;
				int k = (*temp_bdyquick_it).second;
				m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + 1] * m_bigmatrix[j][k];
				m_bigmatrix[j][k] = 0;
			}
		}
		else if (bdybase_count == 1)
		{
			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			{
				int j = (*temp_bdyquick_it).first;
				int k = (*temp_bdyquick_it).second;
				m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3] * m_bigmatrix[j][k];
				m_bigmatrix[j][k] = 0;
			}
		}
		else
		{
			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			{
				int j = (*temp_bdyquick_it).first;
				int k = (*temp_bdyquick_it).second;
				for (int l = 0; l < 3; l++)
				{
					if ((m_colids[j][k]) == m_globalid * 3 + l)
					{
						m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
						break;
					}
				}
				m_bigmatrix[j][k] = 0;
			}
		}*/
		if (bdybase_count < subnum)
		{
			for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
			{
				if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
				{
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
				}
				else
					m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
			}
			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
		}
		else
		{
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}
		}
		/*if (bdybase_count == 0)
		{
			for (int k = 0; k < m_colids[m_globalid * 3 + 1].size(); k++)
			{
				if (m_colids[m_globalid * 3 + 1][k] == m_globalid * 3 + 1)
				{
					m_bigmatrix[m_globalid * 3 + 1][k] = 1;
				}
				else
					m_bigmatrix[m_globalid * 3 + 1][k] = 0;
			}
			m_bigrhs[m_globalid * 3 + 1] = m_bdycoef[m_solidid][m_baseid * 3 + 1];
		}
		else if (bdybase_count == 1)
		{
			for (int k = 0; k < m_colids[m_globalid * 3].size(); k++)
			{
				if (m_colids[m_globalid * 3][k] == m_globalid * 3)
				{
					m_bigmatrix[m_globalid * 3][k] = 1;
				}
				else
					m_bigmatrix[m_globalid * 3][k] = 0;
			}
			m_bigrhs[m_globalid * 3] = m_bdycoef[m_solidid][m_baseid * 3];
		}
		else
		{
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}
		}*/
		bdybase_count++;
	}
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(1);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize * 3; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first + 1);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num + 1);
	}
}

//void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_mutires(std::vector<std::vector<double>> &m_bdycoef, std::vector<BdyBase> &bdybaseid,
//	std::vector<std::vector<Matrix>> &local_matrixs, std::vector<double> &m_bigrhs, std::vector<std::vector<std::vector<double>>> &local_bigrhs,
//	std::vector<int> &r_rowids, std::vector<int> &r_colids, std::vector<double> &r_vals,
//	std::vector<std::vector<int>> &m_colids, std::vector<std::vector<std::pair<int, int>>> &m_bdyquick,
//	std::vector<std::vector<V3f>> &centers,std::vector<std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>>> &matrixs_quick,
//	std::vector<std::vector<std::vector<int>>> &matrixs_quick_check)
//{
//	/*---------------初始化全局刚度矩阵和右端项-----------*/
//	int m_GlobalSize = m_solids.Global_size;
//	m_colids.resize(m_GlobalSize * 3);
//
//	m_bigrhs.resize(m_GlobalSize * 3);
//	Matrix m_bigmatrix;
//	m_bigmatrix.resize(m_GlobalSize * 3);
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
//		{
//			for (int k = 0; k < 3; k++)
//				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize((2 * 3 + 2)*(2 * 3 + 2)*(2 * 3 + 2) * martix_size, 0);
//		}
//	}
//
//	std::vector<double> m_Weight, m_IntePt;
//	int m_GaussNum = 4;
//
//	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4
//
//	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		if(i%100 == 0)
//			std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
//		std::vector<double> m_uIntervel{ 0,0.5,1 }, m_vIntervel{ 0,0.5,1 }, m_wIntervel{ 0,0.5,1 };
//		// std::vector<double> m_uIntervel{ 0,0.25,0.5,0.75,1 }, m_vIntervel{ 0,0.25,0.5,0.75,1 }, m_wIntervel{ 0,0.25,0.5,0.75,1 };
//		int m_uLocation, m_vLocation, m_wLocation;
//		std::vector<double> m_u, m_v, m_w;
//		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
//		//m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
//
//		std::vector<Matrix> local_matrix;
//		std::vector<std::vector<double>> local_bigrh;
//		std::vector<std::vector<std::pair<int, int>>> m_bdyquick;
//		std::vector<V3f> center;
//
//		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
//		int m_BaseIndex;
//		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
//		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
//		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
//		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
//		std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> temtemmatris;
//		std::vector<std::vector<int>> temtemmatrischeck;
//		/***************************构造刚度矩阵按照单元循环****************************************/
//		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
//		{
//			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
//				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
//					+ m_uIntervel[u_interval + 1])*m_IntePt[gaussnum])*0.5;
//			//确定u方向这些高斯点所在区间的位置
//			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);
//
//			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
//			{
//				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
//					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
//						+ m_vIntervel[v_intervel + 1])*m_IntePt[gaussnum])*0.5;
//				//确定v方向这些高斯点所在区间的位置
//				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);
//
//				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
//				{
//					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
//						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
//							+ m_wIntervel[w_intervel + 1])*m_IntePt[gaussnum])*0.5;
//					//确定w方向这些高斯点所在区间的位置
//					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);
//
//					/*-------初始化这部分矩阵的信息------*/
//					V3f m_center(0,0,0);
//					int center_num = 0;
//					std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
//					double u_ders[2], v_ders[2], w_ders[2];
//					std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
//					std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
//					
//					m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
//					Matrix m_LocalMatrix;
//					m_LocalMatrix.resize(m_MatrixSize * 3);
//					for (int k = 0; k < m_MatrixSize * 3; k++)
//					{
//						m_LocalMatrix[k].resize(m_MatrixSize * 3, 0);
//					}
//					std::vector<double> m_bigrh(m_MatrixSize * 3, 0);
//					/*--------------初始化结束---------------*/
//
//					m_det.clear(); m_InverseMatrix.clear();
//					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
//					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
//					std::vector<int> m_baseindex;
//					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
//						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
//							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
//							{
//								double m_tempvalue1 = 0, m_tempvalue2 = 0, m_tempvalue3 = 0;
//								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
//								m_baseindex.push_back(m_BaseIndex);
//								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
//								for (int m = 0; m < m_GaussNum; m++)
//								{
//									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
//									for (int l = 0; l < m_GaussNum; l++)
//									{
//										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
//										for (int k = 0; k < m_GaussNum; k++)
//										{
//											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);
//
//											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
//											center_num++;
//											m_center = m_center + gasv;
//											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
//											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
//											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
//											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
//											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
//											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
//											int m_index = m * m_GaussNum*m_GaussNum + l * m_GaussNum + k;
//											m_tempvalue1 += u_ders[0] * v_ders[0] * w_ders[0] * f1(gasv.x, gasv.y, gasv.z)*
//												fabs(m_det[m_index])*m_Weight[k] * m_Weight[l] * m_Weight[m];
//											m_tempvalue2 += u_ders[0] * v_ders[0] * w_ders[0] * f2(gasv.x, gasv.y, gasv.z)*
//												fabs(m_det[m_index])*m_Weight[k] * m_Weight[l] * m_Weight[m];
//											m_tempvalue3 += u_ders[0] * v_ders[0] * w_ders[0] * f3(gasv.x, gasv.y, gasv.z)*
//												fabs(m_det[m_index])*m_Weight[k] * m_Weight[l] * m_Weight[m];
//										}
//									}
//								}
//								double temp = (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval])*(m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
//									*(m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
//								m_bigrh[3 * m_BaseIndex] += m_tempvalue1 * temp;
//								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue1 * temp) * (umin + (1 - umin)*pow(1, p));
//								m_bigrh[3 * m_BaseIndex + 1] += m_tempvalue2 * temp;
//								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 1] += (m_tempvalue2 * temp) * (umin + (1 - umin)*pow(1, p));
//								m_bigrh[3 * m_BaseIndex + 2] += m_tempvalue3 * temp;
//								m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 2] += (m_tempvalue3 * temp) * (umin + (1 - umin)*pow(1, p));
//							}
//					local_bigrh.push_back(m_bigrh);
//					m_center = m_center / center_num;
//					center.push_back(m_center);
//
//					int m_ptNum = m_BaseDValueu[0].size();
//					for (int ii = 0; ii < m_baseindex.size(); ii++)
//					{
//						for (int j = 0; j < m_baseindex.size(); j++)
//						{
//							for (int k = 0; k < m_ptNum; k++)
//							{
//								int m_uptIndex = (k % (m_GaussNum*m_GaussNum)) % m_GaussNum;
//								int m_vptIndex = (k % (m_GaussNum*m_GaussNum)) / m_GaussNum;
//								int m_wptIndex = k / (m_GaussNum*m_GaussNum);
//								int m_index;
//								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
//								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
//								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
//								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
//								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
//								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
//								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
//								double temp = fabs(m_det[k])*m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1]
//									- m_uIntervel[u_interval])*(m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
//									*(m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
//								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3] += ((lbda + 2 * mu)*m_Rx*m_Qx + mu * (m_Ry*m_Qy + m_Rz * m_Qz))*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 1] += (lbda*m_Rx*m_Qy + mu * m_Ry*m_Qx)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 2] += (lbda*m_Rx*m_Qz + mu * m_Rz*m_Qx)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3] += (lbda*m_Ry*m_Qx + mu * m_Rx*m_Qy)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 1] += ((lbda + 2 * mu)*m_Ry*m_Qy + mu * (m_Rx*m_Qx + m_Rz * m_Qz))*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 2] += (lbda*m_Ry*m_Qz + mu * m_Rz*m_Qy)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3] += (lbda*m_Rz*m_Qx + mu * m_Rx*m_Qz)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 1] += (lbda*m_Rz*m_Qy + mu * m_Ry*m_Qz)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 2] += ((lbda + 2 * mu)*m_Rz*m_Qz + mu * (m_Rx*m_Qx + m_Ry * m_Qy))*temp;
//							}
//						}
//					}
//					local_matrix.push_back(m_LocalMatrix);
//					
//					//矩阵中的非零元素全部存起来
//					std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> temp_Matrixquick;
//					std::vector<int> temp_Matrixquickcheck;
//					for (int ii = 0; ii < m_MatrixSize; ii++)
//					{
//						for (int j = 0; j < m_MatrixSize; j++)
//						{
//							for (int k = 0; k < 3; k++)
//							{
//								for (int l = 0; l < 3; l++)
//								{
//									if (fabs(local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l]) > 1.0e-7)
//									{
//										temp_Matrixquick.push_back({{ii,k},{j,l}});
//										if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
//										{
//											temp_Matrixquickcheck.push_back(0);
//											m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
//											int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//											m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(1, p));
//										}
//										else
//										{
//											temp_Matrixquickcheck.push_back(1);
//											int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//											m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(1, p));
//										}
//									}
//								}
//							}
//						}
//					}
//					temtemmatris.push_back(temp_Matrixquick);
//					temtemmatrischeck.push_back(temp_Matrixquickcheck);
//				}
//				
//			}
//		}
//		local_matrixs.push_back(local_matrix);
//		local_bigrhs.push_back(local_bigrh);
//		centers.push_back(center);
//		matrixs_quick.push_back(temtemmatris);
//		matrixs_quick_check.push_back(temtemmatrischeck);
//	}
//	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
//	//处理边界条件
//	//处理刚度矩阵和右端项
//	std::cout << "Start processing boundary conditions" << std::endl;
//	int bdybase_count = 0;
//	for (int i = 0; i < bdybaseid.size(); i++)
//	{
//		std::vector<std::pair<int, int>> temp_bdyquick;
//		if(i%100 == 0)
//			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
//		int m_baseid = bdybaseid[i].local_id;
//		int m_solidid = bdybaseid[i].solid_id;
//		int m_globalid = bdybaseid[i].global_id;
//		for (int j = 0; j < m_GlobalSize * 3; j++)
//		{
//			for (int k = 0; k < m_colids[j].size(); k++)
//			{
//				if (bdybase_count < subnum)
//				{
//					if ((m_colids[j][k]) == m_globalid * 3 + forceorder)
//					{
//						temp_bdyquick.push_back(std::pair<int, int>(j, k));
//						m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
//						m_bigmatrix[j][k] = 0;
//						break;
//					}
//				}
//				else
//				{
//					for (int l = 0; l < 3; l++)
//					{
//						if ((m_colids[j][k]) == m_globalid * 3 + l)
//						{
//							temp_bdyquick.push_back(std::pair<int, int>(j, k));
//							m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//							m_bigmatrix[j][k] = 0;
//							break;
//						}
//					}
//				}
//			}
//			// for (int k = 0; k < m_colids[j].size(); k++)
//			// {
//				// if (bdybase_count < 16)
//				// {
//					// if ((m_colids[j][k]) == m_globalid * 3 + 1)
//					// {
//						// temp_bdyquick.push_back(std::pair<int, int>(j, k));
//						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + 1] * m_bigmatrix[j][k];
//						// m_bigmatrix[j][k] = 0;
//						// break;
//					// }
//				// }
//				// else if (bdybase_count < 32)
//				// {
//					// if ((m_colids[j][k]) == m_globalid * 3)
//					// {
//						// temp_bdyquick.push_back(std::pair<int, int>(j, k));
//						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3] * m_bigmatrix[j][k];
//						// m_bigmatrix[j][k] = 0;
//						// break;
//					// }
//					// for (int l = 0; l < 2; l++)
//					// {
//						// if ((m_colids[j][k]) == m_globalid * 3 + l)
//						// {
//							// temp_bdyquick.push_back(std::pair<int, int>(j, k));
//							// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//							// m_bigmatrix[j][k] = 0;
//							// break;
//						// }
//					// }
//				// }
//				// else
//				// {
//					// for (int l = 0; l < 3; l++)
//					// {
//						// if ((m_colids[j][k]) == m_globalid * 3 + l)
//						// {
//							// temp_bdyquick.push_back(std::pair<int, int>(j, k));
//							// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//							// m_bigmatrix[j][k] = 0;
//							// break;
//						// }
//					// }
//				// }
//			// }
//		}
//		m_bdyquick.push_back(temp_bdyquick);
//		if (bdybase_count < subnum)
//		{
//			for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
//			{
//				if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
//				{
//					m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
//				}
//				else
//					m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
//			}
//			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
//		}
//		else
//		{
//			for (int m = 0; m < 3; m++)
//			{
//				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				{
//					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
//					{
//						m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					}
//					else
//						m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				}
//				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			}
//		}
//		// if (bdybase_count < 16)
//		// {
//			// for (int k = 0; k < m_colids[m_globalid * 3 + 1].size(); k++)
//			// {
//				// if (m_colids[m_globalid * 3 + 1][k] == m_globalid * 3 + 1)
//				// {
//					// m_bigmatrix[m_globalid * 3 + 1][k] = 1;
//				// }
//				// else
//					// m_bigmatrix[m_globalid * 3 + 1][k] = 0;
//			// }
//			// m_bigrhs[m_globalid * 3 + 1] = m_bdycoef[m_solidid][m_baseid * 3 + 1];
//		// }
//		// else if (bdybase_count < 32)
//		// {
//			// for (int k = 0; k < m_colids[m_globalid * 3].size(); k++)
//			// {
//				// if (m_colids[m_globalid * 3][k] == m_globalid * 3)
//				// {
//					// m_bigmatrix[m_globalid * 3 ][k] = 1;
//				// }
//				// else
//					// m_bigmatrix[m_globalid * 3][k] = 0;
//			// }
//			// m_bigrhs[m_globalid * 3] = m_bdycoef[m_solidid][m_baseid * 3];
//			// for (int m = 0; m < 2; m++)
//			// {
//				// for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				// {
//					// if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						// || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
//					// {
//						// m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					// }
//					// else
//						// m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				// }
//				// m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			// }
//		// }
//		// else
//		// {
//			// for (int m = 0; m < 3; m++)
//			// {
//				// for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				// {
//					// if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						// || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
//					// {
//						// m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					// }
//					// else
//						// m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				// }
//				// m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			// }
//		// }
//		bdybase_count++;
//	}
//	std::cout << "Processing boundary conditions completely" << std::endl;
//	r_rowids.push_back(1);
//	int ele_num = 0;
//	for (int i = 0; i < m_GlobalSize * 3; i++)
//	{
//		std::map<int, double> ele_colid;
//		for (int j = 0; j < m_bigmatrix[i].size(); j++)
//		{
//			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
//			{
//				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
//			}
//		}
//		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
//		{
//			if ((*ele_colid_it).first >= i)
//			{
//				r_colids.push_back((*ele_colid_it).first + 1);
//				r_vals.push_back((*ele_colid_it).second);
//				ele_num++;
//			}
//		}
//		r_rowids.push_back(ele_num + 1);
//	}
//}
//
//void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_mutires_quick(std::vector<std::vector<double>> &m_bdycoef, std::vector<BdyBase> &bdybaseid,
//	std::vector<std::vector<Matrix>> &local_matrixs, std::vector<double> &m_bigrhs, std::vector<std::vector<std::vector<double>>> &local_bigrhs,
//	std::vector<int> &r_rowids, std::vector<int> &r_colids, std::vector<double> &r_vals,
//	std::vector<std::vector<int>> &m_colids, std::vector<std::vector<std::pair<int, int>>> &m_bdyquick,std::vector<std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>>>  &matrixs_quick,
//	std::vector<std::vector<std::vector<int>>> &matrixs_quick_check)
//{
//	int m_GlobalSize = m_solids.Global_size;
//	m_bigrhs.clear();
//	m_bigrhs.resize(m_GlobalSize * 3, 0);
//	m_colids.clear();
//	m_colids.resize(m_GlobalSize * 3);
//	
//	Matrix m_bigmatrix;
//	m_bigmatrix.resize(m_GlobalSize * 3);
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
//		{
//			for (int k = 0; k < 3; k++)
//				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize((2 * 3 + 2)*(2 * 3 + 2)*(2 * 3 + 2) * martix_size, 0);
//		}
//	}
//
//	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		if(i%100 == 0)
//			std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
//		// int m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
//		// int m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
//		// int m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
//		// int m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
//
//		for (int lm = 0; lm < local_matrixs[i].size(); ++lm)
//		{
//			//矩阵中的非零元素全部存起来
//			for(int x = 0;x<matrixs_quick[i][lm].size();++x)
//			{
//				int ii = matrixs_quick[i][lm][x].first.first;
//				int k = matrixs_quick[i][lm][x].first.second;
//				int j = matrixs_quick[i][lm][x].second.first;
//				int l = matrixs_quick[i][lm][x].second.second;
//				int check = matrixs_quick_check[i][lm][x];
//				if (check == 0) //全矩阵存储
//				{
//					m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
//					int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//					m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(material_distribution_mutires[i][lm], p));
//				}
//				else
//				{
//					int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//					m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(material_distribution_mutires[i][lm], p));
//				}
//			}
//			//矩阵中的非零元素全部存起来
//			// for (int ii = 0; ii < m_MatrixSize; ii++)
//			// {
//				// for (int j = 0; j < m_MatrixSize; j++)
//				// {
//					// for (int k = 0; k < 3; k++)
//					// {
//						// for (int l = 0; l < 3; l++)
//						// {
//							// if (fabs(local_matrixs[i][lm][ii * 3 + k][j * 3 + l]) > 1.0e-7)
//							// {
//								// if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
//								// {
//									// m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
//									// int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//									// m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(material_distribution_mutires[i][lm], p));
//								// }
//								// else
//								// {
//									// int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//									// m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(material_distribution_mutires[i][lm], p));
//								// }
//							// }
//						// }
//					// }
//				// }
//			// }
//		}
//		
//	}
//	for (int i = 0; i < local_bigrhs.size(); i++)
//	{
//		for (int lb = 0; lb < local_bigrhs[i].size(); ++lb)
//		{
//			for (int k = 0; k < 3; k++)
//				for (int j = 0; j < local_bigrhs[i].size() / 3; j++)
//					m_bigrhs[m_solids.Local_to_Global[i][j] * 3 + k] += local_bigrhs[i][lb][j * 3 + k] /** (umin + (1 - umin)*pow(material_distribution_mutires[i][lb], p))*/;
//		}	
//	}
//	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
//	//处理边界条件
//	//处理刚度矩阵和右端项
//	std::cout << "Start processing boundary conditions" << std::endl;
//	int bdybase_count = 0;
//	for (int i = 0; i < bdybaseid.size(); i++)
//	{
//		if(i%1000 == 0)
//			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
//		int m_baseid = bdybaseid[i].local_id;
//		int m_solidid = bdybaseid[i].solid_id;
//		int m_globalid = bdybaseid[i].global_id;
//		std::vector<std::pair<int, int>> temp_bdyquick = m_bdyquick[i];
//		if (bdybase_count < subnum)
//		{
//			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			{
//				int j = (*temp_bdyquick_it).first;
//				int k = (*temp_bdyquick_it).second;
//				m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
//				m_bigmatrix[j][k] = 0;
//			}
//		}
//		else
//		{
//			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			{
//				int j = (*temp_bdyquick_it).first;
//				int k = (*temp_bdyquick_it).second;
//				for (int l = 0; l < 3; l++)
//				{
//					if ((m_colids[j][k]) == m_globalid * 3 + l)
//					{
//						m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//						break;
//					}
//				}
//				m_bigmatrix[j][k] = 0;
//			}
//		}
//		// if (bdybase_count < 16)
//		// {
//			// for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			// {
//				// int j = (*temp_bdyquick_it).first;
//				// int k = (*temp_bdyquick_it).second;
//				// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + 1] * m_bigmatrix[j][k];
//				// m_bigmatrix[j][k] = 0;
//			// }
//		// }
//		// else if (bdybase_count < 32)
//		// {
//			// for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			// {
//				// int j = (*temp_bdyquick_it).first;
//				// int k = (*temp_bdyquick_it).second;
//				// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3] * m_bigmatrix[j][k];
//				// m_bigmatrix[j][k] = 0;
//			// }
//			// for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			// {
//				// int j = (*temp_bdyquick_it).first;
//				// int k = (*temp_bdyquick_it).second;
//				// for (int l = 0; l < 2; l++)
//				// {
//					// if ((m_colids[j][k]) == m_globalid * 3 + l)
//					// {
//						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//						// break;
//					// }
//				// }
//				// m_bigmatrix[j][k] = 0;
//			// }
//		// }
//		// else
//		// {
//			// for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			// {
//				// int j = (*temp_bdyquick_it).first;
//				// int k = (*temp_bdyquick_it).second;
//				// for (int l = 0; l < 3; l++)
//				// {
//					// if ((m_colids[j][k]) == m_globalid * 3 + l)
//					// {
//						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//						// break;
//					// }
//				// }
//				// m_bigmatrix[j][k] = 0;
//			// }
//		// }
//		if (bdybase_count < subnum)
//		{
//			for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
//			{
//				if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
//				{
//					m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
//				}
//				else
//					m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
//			}
//			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
//		}
//		else
//		{
//			for (int m = 0; m < 3; m++)
//			{
//				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				{
//					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
//					{
//						m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					}
//					else
//						m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				}
//				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			}
//		}
//		// if (bdybase_count < 16)
//		// {
//			// for (int k = 0; k < m_colids[m_globalid * 3 + 1].size(); k++)
//			// {
//				// if (m_colids[m_globalid * 3 + 1][k] == m_globalid * 3 + 1)
//				// {
//					// m_bigmatrix[m_globalid * 3 + 1][k] = 1;
//				// }
//				// else
//					// m_bigmatrix[m_globalid * 3 + 1][k] = 0;
//			// }
//			// m_bigrhs[m_globalid * 3 + 1] = m_bdycoef[m_solidid][m_baseid * 3 + 1];
//		// }
//		// else if (bdybase_count < 32)
//		// {
//			// for (int k = 0; k < m_colids[m_globalid * 3].size(); k++)
//			// {
//				// if (m_colids[m_globalid * 3][k] == m_globalid * 3)
//				// {
//					// m_bigmatrix[m_globalid * 3][k] = 1;
//				// }
//				// else
//					// m_bigmatrix[m_globalid * 3][k] = 0;
//			// }
//			// m_bigrhs[m_globalid * 3] = m_bdycoef[m_solidid][m_baseid * 3];
//			// for (int m = 0; m < 2; m++)
//			// {
//				// for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				// {
//					// if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						// || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
//					// {
//						// m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					// }
//					// else
//						// m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				// }
//				// m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			// }
//		// }
//		// else
//		// {
//			// for (int m = 0; m < 3; m++)
//			// {
//				// for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				// {
//					// if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						// || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
//					// {
//						// m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					// }
//					// else
//						// m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				// }
//				// m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			// }
//		// }
//		bdybase_count++;
//	}
//	std::cout << "Processing boundary conditions completely" << std::endl;
//	r_rowids.push_back(1);
//	int ele_num = 0;
//	for (int i = 0; i < m_GlobalSize * 3; i++)
//	{
//		std::map<int, double> ele_colid;
//		for (int j = 0; j < m_bigmatrix[i].size(); j++)
//		{
//			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
//			{
//				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
//			}
//		}
//		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
//		{
//			if ((*ele_colid_it).first >= i)
//			{
//				r_colids.push_back((*ele_colid_it).first + 1);
//				r_vals.push_back((*ele_colid_it).second);
//				ele_num++;
//			}
//		}
//		r_rowids.push_back(ele_num + 1);
//	}
//}


//void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_mutires(std::vector<std::vector<double>> &m_bdycoef, std::vector<BdyBase> &bdybaseid,
//	std::vector<Matrix> &local_matrix, std::vector<double> &m_bigrhs, std::vector<std::vector<std::vector<double>>> &local_bigrhs,
//	std::vector<int> &r_rowids, std::vector<int> &r_colids, std::vector<double> &r_vals,
//	std::vector<std::vector<int>> &m_colids, std::vector<std::vector<std::pair<int, int>>> &m_bdyquick,
//	std::vector<std::vector<V3f>> &centers)
//{	
//	std::cout << "11_18 1:10 new file check" << std::endl;
//	/*---------------初始化全局刚度矩阵和右端项-----------*/
//	int m_GlobalSize = m_solids.Global_size;
//	m_colids.resize(m_GlobalSize * 3);
//
//	Matrix m_bigmatrix;
//
//	m_bigmatrix.resize(m_GlobalSize * 3);
//	 
//
//	std::vector<double> m_Weight, m_IntePt;
//	int m_GaussNum = 4;
//
//	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4
//
//
//	std::vector<double> m_BaseDValueu_global, m_BaseDValuev_global, m_BaseDValuew_global;
//	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		if (i % 100 == 0)
//			std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
//		std::vector<double> m_uIntervel{ 0,1 }, m_vIntervel{ 0,1 }, m_wIntervel{ 0,1 };
//		// std::vector<double> m_uIntervel{ 0,0.25,0.5,0.75,1 }, m_vIntervel{ 0,0.25,0.5,0.75,1 }, m_wIntervel{ 0,0.25,0.5,0.75,1 };
//		//std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
//		int m_uLocation, m_vLocation, m_wLocation;
//		std::vector<double> m_u, m_v, m_w;
//		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
//		//m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
//
//		std::vector<std::vector<double>> local_bigrh;
//		// std::vector<std::vector<std::pair<int, int>>> m_bdyquick;
//		std::vector<V3f> center;
//
//		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
//		int m_BaseIndex;
//		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
//		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
//		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
//		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
//		/***************************构造刚度矩阵按照单元循环****************************************/
//		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
//		{
//			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
//				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
//					+ m_uIntervel[u_interval + 1])*m_IntePt[gaussnum])*0.5;
//			//确定u方向这些高斯点所在区间的位置
//			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);
//
//			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
//			{
//				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
//					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
//						+ m_vIntervel[v_intervel + 1])*m_IntePt[gaussnum])*0.5;
//				//确定v方向这些高斯点所在区间的位置
//				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);
//
//				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
//				{
//					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
//						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
//							+ m_wIntervel[w_intervel + 1])*m_IntePt[gaussnum])*0.5;
//					//确定w方向这些高斯点所在区间的位置
//					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);
//
//					/*-------初始化这部分矩阵的信息------*/
//					V3f m_center(0, 0, 0);
//					int center_num = 0;
//					std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
//					double u_ders[2], v_ders[2], w_ders[2];
//					std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
//					std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
//
//					m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
//					Matrix m_LocalMatrix;
//					m_LocalMatrix.resize(m_MatrixSize * 3);
//					for (int k = 0; k < m_MatrixSize * 3; k++)
//					{
//						m_LocalMatrix[k].resize(m_MatrixSize * 3, 0);
//					}
//					std::vector<double> m_bigrh(m_MatrixSize * 3, 0);
//					/*--------------初始化结束---------------*/
//
//					m_det.clear(); m_InverseMatrix.clear();
//					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
//					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
//					std::vector<int> m_baseindex;
//					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
//						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
//							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
//							{
//								double m_tempvalue1 = 0, m_tempvalue2 = 0, m_tempvalue3 = 0;
//								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
//								m_baseindex.push_back(m_BaseIndex);
//								if (m_BaseDValueu_global.size() == 0) {
//									m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
//									for (int m = 0; m < m_GaussNum; m++)
//									{
//										m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
//										for (int l = 0; l < m_GaussNum; l++)
//										{
//											m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
//											for (int k = 0; k < m_GaussNum; k++)
//											{
//												m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);
//
//												V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
//												center_num++;
//												m_center = m_center + gasv;
//												double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
//												m_BaseDValueu[m_BaseIndex].push_back(m_temp);
//												m_temp = u_ders[0] * v_ders[1] * w_ders[0];
//												m_BaseDValuev[m_BaseIndex].push_back(m_temp);
//												m_temp = u_ders[0] * v_ders[0] * w_ders[1];
//												m_BaseDValuew[m_BaseIndex].push_back(m_temp);
//												/*int m_index = m * m_GaussNum*m_GaussNum + l * m_GaussNum + k;*/
//												/*m_tempvalue1 += u_ders[0] * v_ders[0] * w_ders[0] * fx(i,m_BaseIndex)*
//													fabs(m_det[m_index])*m_Weight[k] * m_Weight[l] * m_Weight[m];
//												m_tempvalue2 += u_ders[0] * v_ders[0] * w_ders[0] * fy(i,m_BaseIndex)*
//													fabs(m_det[m_index])*m_Weight[k] * m_Weight[l] * m_Weight[m];
//												m_tempvalue3 += u_ders[0] * v_ders[0] * w_ders[0] * fz(i,m_BaseIndex)*
//													fabs(m_det[m_index])*m_Weight[k] * m_Weight[l] * m_Weight[m];*/
//											}
//										}
//									}
//									m_BaseDValueu_global.resize(m_BaseDValueu[m_BaseIndex].size());
//									m_BaseDValuev_global.resize(m_BaseDValuev[m_BaseIndex].size());
//									m_BaseDValuew_global.resize(m_BaseDValuew[m_BaseIndex].size());
//									m_BaseDValueu_global = m_BaseDValueu[m_BaseIndex];
//									m_BaseDValuev_global = m_BaseDValuev[m_BaseIndex];
//									m_BaseDValuew_global = m_BaseDValuew[m_BaseIndex];
//								}
//								else {
//									m_BaseDValueu[m_BaseIndex].resize(m_BaseDValueu_global.size());
//									m_BaseDValuev[m_BaseIndex].resize(m_BaseDValuev_global.size());
//									m_BaseDValuew[m_BaseIndex].resize(m_BaseDValuew_global.size());
//									m_BaseDValueu[m_BaseIndex] = m_BaseDValueu_global;
//									m_BaseDValuev[m_BaseIndex] = m_BaseDValuev_global;
//									m_BaseDValuew[m_BaseIndex] = m_BaseDValuew_global;
//									//std::cout << "m_BaseDvalueu.size(): " << m_BaseDValueu[m_BaseIndex].size() << std::endl;
//									for (int m = 0; m < m_GaussNum; m++)
//									{
//										for (int l = 0; l < m_GaussNum; l++)
//										{
//											for (int k = 0; k < m_GaussNum; k++)
//											{
//
//												V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
//												center_num++;
//												m_center = m_center + gasv;
//												
//											}
//										}
//									}
//								}
//								
//								/*double temp = (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval])*(m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
//									*(m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;*/
//								
//								//右端项直接通过处理刚度矩阵的时候进行设置 这边就不设置了 全部注释
//								// m_bigrh[3 * m_BaseIndex] = fx(i,m_BaseIndex);
//								// m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex]] = fx(i,m_BaseIndex);
//								// m_bigrh[3 * m_BaseIndex + 1] = fy(i,m_BaseIndex);
//								// m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 1] = fy(i,m_BaseIndex);
//								// m_bigrh[3 * m_BaseIndex + 2] = fz(i,m_BaseIndex);
//								// m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 2] = fz(i,m_BaseIndex);
//								
//								//m_bigrh[3 * m_BaseIndex] += m_tempvalue1 * temp;
//								//m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue1 * temp) * (umin + (1 - umin)*pow(1, p));
//								//m_bigrh[3 * m_BaseIndex + 1] += m_tempvalue2 * temp;
//								//m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 1] += (m_tempvalue2 * temp) * (umin + (1 - umin)*pow(1, p));
//								//m_bigrh[3 * m_BaseIndex + 2] += m_tempvalue3 * temp;
//								//m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 2] += (m_tempvalue3 * temp) * (umin + (1 - umin)*pow(1, p));
//							}
//					local_bigrh.push_back(m_bigrh);
//					m_center = m_center / center_num;
//					center.push_back(m_center);
//
//					int m_ptNum = m_BaseDValueu[0].size();
//					for (int ii = 0; ii < m_baseindex.size(); ii++)
//					{
//						for (int j = 0; j < m_baseindex.size(); j++)
//						{
//							for (int k = 0; k < m_ptNum; k++)
//							{
//								int m_uptIndex = (k % (m_GaussNum*m_GaussNum)) % m_GaussNum;
//								int m_vptIndex = (k % (m_GaussNum*m_GaussNum)) / m_GaussNum;
//								int m_wptIndex = k / (m_GaussNum*m_GaussNum);
//								int m_index;
//								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
//								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
//								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
//								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
//								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
//								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
//								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
//								double temp = fabs(m_det[k])*m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1]
//									- m_uIntervel[u_interval])*(m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
//									*(m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
//								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3] += ((lbda + 2 * mu)*m_Rx*m_Qx + mu * (m_Ry*m_Qy + m_Rz * m_Qz))*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 1] += (lbda*m_Rx*m_Qy + mu * m_Ry*m_Qx)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 2] += (lbda*m_Rx*m_Qz + mu * m_Rz*m_Qx)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3] += (lbda*m_Ry*m_Qx + mu * m_Rx*m_Qy)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 1] += ((lbda + 2 * mu)*m_Ry*m_Qy + mu * (m_Rx*m_Qx + m_Rz * m_Qz))*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 2] += (lbda*m_Ry*m_Qz + mu * m_Rz*m_Qy)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3] += (lbda*m_Rz*m_Qx + mu * m_Rx*m_Qz)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 1] += (lbda*m_Rz*m_Qy + mu * m_Ry*m_Qz)*temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 2] += ((lbda + 2 * mu)*m_Rz*m_Qz + mu * (m_Rx*m_Qx + m_Ry * m_Qy))*temp;
//							}
//						}
//					}
//					local_matrix.push_back(m_LocalMatrix);
//					
//					//矩阵中的非零元素全部存起来
//					for (int ii = 0; ii < m_MatrixSize; ii++)
//					{
//						for (int j = 0; j < m_MatrixSize; j++)
//						{
//							for (int k = 0; k < 3; k++)
//							{
//								for (int l = 0; l < 3; l++)
//								{
//									if (m_solids.Local_to_Global[i][j] * 3 + l >= m_solids.Local_to_Global[i][ii] * 3 + k && fabs(local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l]) > 1.0e-7)
//									{
//										//if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
//										//{
//										//	m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
//										//	int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//										//	if(m_size - 1 >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//										//	{
//										//		m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(m_size,0);
//										//	}
//										//	m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(1, p));
//										//}
//										//else
//										//{
//										//	int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//										//	if(tmp_id >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//										//	{
//										//		m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(tmp_id + 1,0);
//										//	}
//										//	m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(1, p));
//										//}
//										if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
//										{
//											m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
//											m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(0);
//											int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//											/*if (m_size - 1 >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//											{
//												m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(m_size, 0);
//											}*/
//											m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(1, p));
//										}
//										else
//										{
//											int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//											/*if (tmp_id >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//											{
//												m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(tmp_id + 1, 0);
//											}*/
//											m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(1, p));
//										}
//									}
//								}
//							}
//						}
//					}
//				}
//				
//			}
//		}
//		local_bigrhs.push_back(local_bigrh);
//		centers.push_back(center);
//	}
//	// int checknum = 0;
//	// for(int i = 0;i<m_bigrhs.size();++i)
//		// if(m_bigrhs[i]<0)
//			// checknum++;
//	// std::cout<<"checknum: "<<checknum<<std::endl;
//	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
//	//处理边界条件
//	//处理刚度矩阵和右端项
//	std::cout << "Start processing boundary conditions" << std::endl;
//	m_bigrhs.clear();
//	m_bigrhs.resize(m_GlobalSize * 3,0);
//	int bdybase_count = 0;
//	for (int i = 0; i < bdybaseid.size(); i++)
//	{
//		std::vector<std::pair<int, int>> temp_bdyquick;
//		if (i % 100 == 0)
//			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
//		int m_baseid = bdybaseid[i].local_id;
//		int m_solidid = bdybaseid[i].solid_id;
//		int m_globalid = bdybaseid[i].global_id;
//		for (int j = 0; j < m_GlobalSize * 3; j++)
//		{
//			for (int k = 0; k < m_colids[j].size(); k++)
//			{
//				if (bdybase_count < subnum)
//				{
//					//if (m_colids[j][k] == m_globalid * 3 + forceorder)
//					//{
//					//	//temp_bdyquick.push_back(std::pair<int, int>(j, k));
//					//	// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
//					//	// m_bigmatrix[j][k] = 0;
//					//	//break;
//					//}
//				}
//				else
//				{
//					for (int l = 0; l < 3; l++)
//					{
//						if (m_colids[j][k] == m_globalid * 3 + l)
//						{
//							temp_bdyquick.push_back(std::pair<int, int>(j, k));
//							// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//							m_bigmatrix[j][k] = 0;
//							break;
//						}
//					}
//				}
//			}
//		}
//		m_bdyquick.push_back(temp_bdyquick);
//		if (bdybase_count < subnum)
//		{
//			// for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
//			// {
//				// if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
//				// {
//					// m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
//				// }
//				// else
//					// m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
//			// }
//			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
//		}
//		else
//		{
//			for (int m = 0; m < 3; m++)
//			{
//				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				{
//					/*if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)*/
//					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3+m)
//					{
//						m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					}
//					else
//						m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				}
//				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			}
//		}
//		bdybase_count++;
//	}
//	int checknum = 0;
//	//for(int i = 0;i< m_bigrhs.size();++i)
//	//	if(m_bigrhs[i]<0)
//	//		checknum++;
//	//std::cout<<"checknum: "<<checknum<<std::endl;
//	std::cout << "Processing boundary conditions completely" << std::endl;
//	r_rowids.push_back(1);
//	int ele_num = 0;
//	for (int i = 0; i < m_GlobalSize * 3; i++)
//	{
//		// if (i % 100 == 0)
//			// std::cout << "Quickly K " << i + 1 << ":" << m_GlobalSize * 3 << std::endl;
//		std::map<int, double> ele_colid;
//		for (int j = 0; j < m_bigmatrix[i].size(); j++)
//		{
//			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
//			{
//				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
//			}
//		}
//		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
//		{
//			if ((*ele_colid_it).first >= i)
//			{
//				r_colids.push_back((*ele_colid_it).first + 1);
//				r_vals.push_back((*ele_colid_it).second);
//				ele_num++;
//			}
//		}
//		r_rowids.push_back(ele_num + 1);
//	}
//	std::cout<<"end1"<<std::endl;
//}
//
//void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_mutires_quick(std::vector<std::vector<double>> &m_bdycoef, std::vector<BdyBase> &bdybaseid,
//	std::vector<Matrix> &local_matrix, std::vector<double> &m_bigrhs, std::vector<std::vector<std::vector<double>>> &local_bigrhs,
//	std::vector<int> &r_rowids, std::vector<int> &r_colids, std::vector<double> &r_vals,
//	std::vector<std::vector<int>> &m_colids, std::vector<std::vector<std::pair<int, int>>> &m_bdyquick)
//{
//	int m_GlobalSize = m_solids.Global_size;
//	m_bigrhs.clear();
//	m_bigrhs.resize(m_GlobalSize * 3, 0);
//
//	Matrix m_bigmatrix;
//	m_bigmatrix.resize(m_GlobalSize * 3);
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
//		{
//			for (int k = 0; k < 3; k++)
//				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize(m_colids[m_solids.Local_to_Global[i][j] * 3 + k].size(), 0);
//		}
//	}
//
//	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		if (i % 100 == 0)
//			std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
//		int m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
//		int m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
//		int m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
//		int m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
//
//
//			//矩阵中的非零元素全部存起来
//		for (int ii = 0; ii < m_MatrixSize; ii++)
//		{
//			for (int j = 0; j < m_MatrixSize; j++)
//			{
//				for (int k = 0; k < 3; k++)
//				{
//					for (int l = 0; l < 3; l++)
//					{
//						if (m_solids.Local_to_Global[i][j] * 3 + l >= m_solids.Local_to_Global[i][ii] * 3 + k && fabs(local_matrix[i][ii * 3 + k][j * 3 + l]) > 1.0e-7)
//						{
//							//if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
//							//{
//							//	m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
//							//	int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//							//	if(m_size - 1 >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//							//	{
//							//		m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(m_size,0);
//							//	}
//							//	m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(material_distribution_mutires[i][lm], p));
//							//}
//							//else
//							//{
//							//	int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//							//	if(tmp_id >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//							//	{
//							//		m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(tmp_id + 1,0);
//							//	}
//							//	m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * (umin + (1 - umin)*pow(material_distribution_mutires[i][lm], p));
//							//}
//							if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
//							{
//								//m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
//								//m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(0);
//								//int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//								//*if (m_size - 1 >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//								//{
//									//m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(m_size, 0);
//								//}*/
//								//m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution_mutires[i][lm], p));
//								;
//							}
//							else
//							{
//								
//
//								int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//								/*if (tmp_id >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//								{
//									m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(tmp_id + 1, 0);
//								}*/
//								
//								/*m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] = m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] 
//								> fabs(local_matrix[i][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution_mutires1[i], p)))?
//									m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id]: 
//									local_matrix[i][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution_mutires1[i], p));*/
//								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrix[i][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution_mutires1[i], p));
//
//							}
//						}
//					}
//				}
//			}
//		} 
//
//	}
//	// for (int i = 0; i < local_bigrhs.size(); i++)
//	// {
//		// for (int lb = 0; lb < local_bigrhs[i].size(); ++lb)
//		// {
//			// for (int k = 0; k < 3; k++)
//				// for (int j = 0; j < local_bigrhs[i].size() / 3; j++)
//					// m_bigrhs[m_solids.Local_to_Global[i][j] * 3 + k] += local_bigrhs[i][lb][j * 3 + k] /** (umin + (1 - umin)*pow(material_distribution_mutires[i][lb], p))*/;
//		// }
//	// }
//	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
//	//处理边界条件
//	//处理刚度矩阵和右端项
//	std::cout << "Start processing boundary conditions" << std::endl;
//	int bdybase_count = 0;
//	for (int i = 0; i < bdybaseid.size(); i++)
//	{
//		if (i % 1000 == 0)
//			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
//		int m_baseid = bdybaseid[i].local_id;    //该控制点对应的局部id
//		int m_solidid = bdybaseid[i].solid_id;   //该控制点对应的体序号
//		int m_globalid = bdybaseid[i].global_id; //该控制点对应的全局id
//		std::vector<std::pair<int, int>> temp_bdyquick = m_bdyquick[i];
//		if (bdybase_count < subnum)   //对于受力部分
//		{
//			
//			// for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			// {
//				// int j = (*temp_bdyquick_it).first;
//				// int k = (*temp_bdyquick_it).second;
//				// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
//				// m_bigmatrix[j][k] = 0;
//				// for (int l = 0; l < 3; l++)
//				// {
//					// if ((m_colids[j][k]) == m_globalid * 3 + l)
//					// {
//						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//						// break;
//					// }
//				// }
//				// m_bigmatrix[j][k] = 0;
//			// }
//		}
//		else    //对于固定部分
//		{
//			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			{
//				int j = (*temp_bdyquick_it).first;
//				int k = (*temp_bdyquick_it).second;
//				//for (int l = 0; l < 3; l++)
//				//{
//				//	if ((m_colids[j][k]) == m_globalid * 3 + l)
//				//	{
//				//		// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//				//		break;
//				//	}
//				//}
//				m_bigmatrix[j][k] = 0;
//			}
//		}
//		if (bdybase_count < subnum)
//		{
//			// for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
//			// {
//				// if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
//				// {
//					// m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
//				// }
//				// else
//					// m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
//			// }
//			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];  //受力的右端项
//		}
//		else
//		{
//			for (int m = 0; m < 3; m++)
//			{
//				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				{
//					/*if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)*/
//					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3+m)
//					{
//						m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					}
//					else
//						m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				}
//				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			}
//		}
//		bdybase_count++;
//	}
//	std::cout << "Processing boundary conditions completely" << std::endl;
//	r_rowids.push_back(1);
//	int ele_num = 0;
//	for (int i = 0; i < m_GlobalSize * 3; i++)
//	{
//		std::map<int, double> ele_colid;  //可自动排序
//		for (int j = 0; j < m_bigmatrix[i].size(); j++)
//		{
//			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
//			{
//				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
//			}
//		}
//		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
//		{
//			if ((*ele_colid_it).first >= i)
//			{
//				r_colids.push_back((*ele_colid_it).first + 1);
//				r_vals.push_back((*ele_colid_it).second);
//				ele_num++;
//			}
//		}
//		r_rowids.push_back(ele_num + 1);
//	}
//}


 //下面是正确的细分之后的求解器刚度矩阵装配
void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_mutires(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
	std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
	std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick,
	std::vector<std::vector<V3f>>& centers)
{
	/*---------------初始化全局刚度矩阵和右端项-----------*/
	int m_GlobalSize = m_solids.Global_size;       // 总控制顶点个数
	m_colids.resize(m_GlobalSize * 3);             // 每个控制顶点有三维
	colsID.resize(m_GlobalSize * 3);
	colsize.resize(m_GlobalSize * 3, 0);
	Matrix m_bigmatrix;

	m_bigmatrix.resize(m_GlobalSize * 3);

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4
	//for (int ii = 0; ii < 4; ++ii) std::cout << "Guass weight: " << m_Weight[ii] << std::endl;

	centers.resize(m_solids.Local_V3f.size());
	local_matrixs.resize(m_solids.Local_V3f.size());
	PreCalDerBasisFunction.clear();
	PreCalDerBasisFunction.resize(8);
	for (int a = 0; a < 2; ++a) {
		std::vector<double> m_u(m_GaussNum);
		for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++) {
			m_u[gaussnum] = (0.5 * a + 0.5 * (a + 1) + (-0.5 * a + 0.5 * (a + 1)) * m_IntePt[gaussnum]) * 0.5;
		}
		for (int b = 0; b < 2; ++b) {

			std::vector<double> m_v(m_GaussNum);
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++) {
				m_v[gaussnum] = (0.5 * b + 0.5 * (b + 1) + (-0.5 * b + 0.5 * (b + 1)) * m_IntePt[gaussnum]) * 0.5;
			}
			for (int c = 0; c < 2; ++c) {
				std::vector<double> m_w(m_GaussNum);
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++) {
					m_w[gaussnum] = (0.5 * c + 0.5 * (c + 1) + (-0.5 * c + 0.5 * (c + 1)) * m_IntePt[gaussnum]) * 0.5;
				}

				size_t ElementOrder = a * 1 + b * 2 + c * 4;
				std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
				std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式

				m_solids.JacobiMatrix_new(0, m_u, m_v, m_w, m_det, m_InverseMatrix, PreCalDerBasisFunction[ElementOrder], false);
			}
		}
	}

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++) {
		for (int ii = 0; ii < 64; ii++)
		{
			size_t globalidx1 = m_solids.Local_to_Global[i][ii] * 3;
			for (int k = 0; k < 3; k++)
			{
				size_t index1 = globalidx1 + k;
				for (int j = 0; j < 64; j++)
				{
					size_t globalidx2 = m_solids.Local_to_Global[i][j] * 3;
					for (int l = 0; l < 3; l++)
					{
						size_t index2 = globalidx2 + l;
						int temp_id = InVecIndex(index2, m_colids[index1]);

						if (temp_id == m_colids[index1].size())
						{
							m_colids[index1].push_back(index2);
							m_bigmatrix[index1].push_back(0);
						}

					}
				}
			}
		}
	}
	std::cout << "Start constructing stiffness matrix and the right end term old" << std::endl;
#pragma omp parallel for num_threads(8)
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{

		if (i % 50 == 0)
			std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		//std::vector<double> m_uIntervel{ 0,1 }, m_vIntervel{ 0,1 }, m_wIntervel{ 0,1 };
		std::vector<double> m_uIntervel{ 0,0.5,1 }, m_vIntervel{ 0,0.5,1 }, m_wIntervel{ 0,0.5,1 };
		// std::vector<double> m_uIntervel{ 0,0.25,0.5,0.75,1 }, m_vIntervel{ 0,0.25,0.5,0.75,1 }, m_wIntervel{ 0,0.25,0.5,0.75,1 };
		//std::vector<double> m_uIntervel, m_vIntervel, m_wIntervel;
		int m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		//m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);

		std::vector<Matrix> local_matrix;
		std::vector<std::vector<double>> local_bigrh;
		// std::vector<std::vector<std::pair<int, int>>> m_bdyquick;
		std::vector<V3f> center;

		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		int m_BaseIndex;
		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
		/***************************构造刚度矩阵按照单元循环****************************************/

		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;


					bool isPreCalFlag = false;
					size_t ElementOrder = u_interval * 1 + v_intervel * 2 + w_intervel * 4;

					if (PreCalDerBasisFunction[ElementOrder].size() != 0) isPreCalFlag = true;


					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					/*-------初始化这部分矩阵的信息------*/
					V3f m_center(0, 0, 0);
					int center_num = 0;
					std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
					double u_ders[2], v_ders[2], w_ders[2];
					std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
					std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，

					m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
					Matrix m_LocalMatrix;
					m_LocalMatrix.resize(m_MatrixSize * 3);
					for (int k = 0; k < m_MatrixSize * 3; k++)
					{
						m_LocalMatrix[k].resize(m_MatrixSize * 3, 0);
					}
					std::vector<double> m_bigrh(m_MatrixSize * 3, 0);
					/*--------------初始化结束---------------*/

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix, PreCalDerBasisFunction[ElementOrder], isPreCalFlag);
					//m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue1 = 0, m_tempvalue2 = 0, m_tempvalue3 = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									//m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									w_ders[0] = PreCalDerBasisFunction[ElementOrder][2 * 32 + (ww * 1 + m * 4) * 2];
									w_ders[1] = PreCalDerBasisFunction[ElementOrder][2 * 32 + (ww * 1 + m * 4) * 2 + 1];
									for (int l = 0; l < m_GaussNum; l++)
									{
										//m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										v_ders[0] = PreCalDerBasisFunction[ElementOrder][1 * 32 + (j * 1 + l * 4) * 2];
										v_ders[1] = PreCalDerBasisFunction[ElementOrder][1 * 32 + (j * 1 + l * 4) * 2 + 1];
										for (int k = 0; k < m_GaussNum; k++)
										{
											//m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);
											u_ders[0] = PreCalDerBasisFunction[ElementOrder][0 * 32 + (ii * 1 + k * 4) * 2];
											u_ders[1] = PreCalDerBasisFunction[ElementOrder][0 * 32 + (ii * 1 + k * 4) * 2 + 1];

											//V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m], PreCalDerBasisFunction[ElementOrder]);
											V3f gasv = m_solids.Geometry_new(i, k, l, m, PreCalDerBasisFunction[ElementOrder]);
											//V3f gasv = m_solids.Geometry_new(i, k, l, m);
											center_num++;
											m_center = m_center + gasv;
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
										}
									}
								}
							}
					local_bigrh.push_back(m_bigrh);
					m_center = m_center / center_num;
					center.push_back(m_center);

					int m_ptNum = m_BaseDValueu[0].size();
					int preuvw[64][3];
					double pretemp[64];
					for (int k = 0; k < 64; ++k) {
						preuvw[k][0] = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
						preuvw[k][1] = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
						preuvw[k][2] = k / (m_GaussNum * m_GaussNum);
						pretemp[k] = fabs(m_det[k]) * m_Weight[preuvw[k][0]] * m_Weight[preuvw[k][1]] * m_Weight[preuvw[k][2]] * (m_uIntervel[u_interval + 1]
							- m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
							* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
					}
					double XXX[64][64];
					double YYY[64][64];
					double ZZZ[64][64];

					for (int k = 0; k < m_ptNum; k++)
					{
						int m_uptIndex = preuvw[k][0];
						int m_vptIndex = preuvw[k][1];
						int m_wptIndex = preuvw[k][2];
						for (int ii = 0; ii < m_baseindex.size(); ii++)
						{
							XXX[k][ii] = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
							YYY[k][ii] = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
							ZZZ[k][ii] = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
						}
					}

					for (int k = 0; k < m_ptNum; k++)
					{
						int m_uptIndex = preuvw[k][0];
						int m_vptIndex = preuvw[k][1];
						int m_wptIndex = preuvw[k][2];
						for (int ii = 0; ii < m_baseindex.size(); ii++)
						{
							double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
							m_Rx = XXX[k][ii];
							m_Ry = YYY[k][ii];
							m_Rz = ZZZ[k][ii];

							for (int j = 0; j < m_baseindex.size(); j++)
							{
								//int m_index;
								m_Qx = XXX[k][j];
								m_Qy = YYY[k][j];
								m_Qz = ZZZ[k][j];
								double temp = pretemp[k];
								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;

							}
						}
					}
					//std::cout << "m_LocalMatrix[0][0] " << std::setprecision(15) << m_LocalMatrix[0][0] << std::endl;
					local_matrix.push_back(m_LocalMatrix);

					//矩阵中的非零元素全部存起来
					for (int ii = 0; ii < m_MatrixSize; ii++)
					{
						size_t globalidx1 = m_solids.Local_to_Global[i][ii] * 3;
						for (int k = 0; k < 3; k++)
						{
							size_t index1 = globalidx1 + k;
							size_t localidx1 = ii * 3 + k;
							for (int j = 0; j < m_MatrixSize; j++)
							{

								size_t globalidx2 = m_solids.Local_to_Global[i][j] * 3;
								//if (globalidx1 > globalidx2) continue;
								for (int l = 0; l < 3; l++)
								{
									size_t index2 = globalidx2 + l;
									size_t localidx2 = j * 3 + l;
									double localval = m_LocalMatrix[localidx1][localidx2];
									if (fabs(localval) > 1e-9)
									{
										int temp_id = InVecIndex(index2, m_colids[index1]);
#pragma omp atomic
										m_bigmatrix[index1][temp_id] += localval;

										//std::cout << "2" << std::endl;
									}

								}
							}
						}
					}
				}

			}
		}
		local_matrixs[i] = local_matrix;
		centers[i] = center;
	}


	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	m_bigrhs.clear();
	m_bigrhs.resize(m_GlobalSize * 3, 0);

	int bdybase_count = 0;

	std::unordered_map<int, int> fixed;
	std::unordered_map<int, int> force;

	m_bdyquick.resize(bdybaseid.size());

	for (int i = 0; i < bdybaseid.size(); i++)
	{
		//std::vector<std::pair<int, int>> temp_bdyquick;
		if (i % 100 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;

		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;

		//m_bdyquick.push_back(temp_bdyquick);
		if (bdybase_count < subnum)
		{
			force[m_globalid * 3 + forceorder] = i;
			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
		}
		else
			// fix
		{
			for (int m = 0; m < 3; m++) 
			{
				fixed[m_globalid * 3 + m] = i;

				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}
		}
		bdybase_count++;
	}
	//#pragma omp parallel for num_threads(8)
	for (int j = 0; j < m_GlobalSize * 3; j++)
	{
		for (int k = 0; k < m_colids[j].size(); k++)
		{
			int cur = m_colids[j][k];
			if (fixed.count(cur)) {
				m_bigmatrix[j][k] = 0;
				m_bdyquick[fixed[cur]].push_back(std::pair<int, int>(j, k));
			}
			else if (force.count(cur)) {
				m_bdyquick[fixed[cur]].push_back(std::pair<int, int>(j, k));
			}
			/*if (bdybase_count < subnum)
			{
				if (m_colids[j][k] == m_globalid * 3 + forceorder)
				{
					temp_bdyquick.push_back(std::pair<int, int>(j, k));

					break;
				}
			}
			else
			{
				for (int l = 0; l < 3; l++)
				{
					if (m_colids[j][k] == m_globalid * 3 + l)
					{
						temp_bdyquick.push_back(std::pair<int, int>(j, k));
						m_bigmatrix[j][k] = 0;
						break;
					}
				}
			}*/
		}
	}

	int checknum = 0;
	for (int i = 0; i < m_bigrhs.size(); ++i)
		if (m_bigrhs[i] < 0)
			checknum++;
	std::cout << "checknum: " << checknum << std::endl;
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(0);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize * 3; i++)
	{
		// if (i % 100 == 0)
			// std::cout << "Quickly K " << i + 1 << ":" << m_GlobalSize * 3 << std::endl;
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num);
	}
	std::cout << "end1" << std::endl;
}



void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_mutires_quick(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
	std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
	std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick)
{
	int m_GlobalSize = m_solids.Global_size;
	m_bigrhs.clear();
	m_bigrhs.resize(m_GlobalSize * 3, 0);

	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize * 3);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			for (int k = 0; k < 3; k++)
				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize(m_colids[m_solids.Local_to_Global[i][j] * 3 + k].size(), 0);
		}
	}

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
#pragma omp parallel for num_threads(8)
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		int m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		int m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		int m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		int m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		for (int lm = 0; lm < local_matrixs[i].size(); ++lm)
			//for(int lm = 0; lm < 1; ++lm)
		{
			double rho = (umin + (1 - umin) * pow(material_distribution_mutires[i][lm], p));
			//矩阵中的非零元素全部存起来
			for (int ii = 0; ii < m_MatrixSize; ii++)
			{
				for (int j = 0; j < m_MatrixSize; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						for (int l = 0; l < 3; l++)
						{
							if (m_solids.Local_to_Global[i][ii] * 3 + k > m_solids.Local_to_Global[i][j] * 3 + l) {}
							else if (fabs(local_matrixs[i][lm][ii * 3 + k][j * 3 + l]) > 1.0e-7)
							{
								int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
								if (tmp_id >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size()) {
									std::cout << "m_bigmatrix 越界:  FunctionBSpline.cpp:6310行" << std::endl;
								}
#pragma omp atomic
								m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * rho;

							}
						}
					}
				}
			}
		}

	}
	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		if (i % 1000 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;


		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;

		std::vector<std::pair<int, int>> temp_bdyquick = m_bdyquick[i];

		if (bdybase_count < subnum)
		{
		}
		else
		{
			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
			{
				int j = (*temp_bdyquick_it).first;
				int k = (*temp_bdyquick_it).second;
				for (int l = 0; l < 3; l++)
				{
					if ((m_colids[j][k]) == m_globalid * 3 + l)
					{
						break;
					}
				}
				m_bigmatrix[j][k] = 0;
			}
		}
		if (bdybase_count < subnum)
		{
			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
		}
		else
		{
			for (int m = 0; m < 3; m++)
			{
				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
				{
					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
					{
						m_bigmatrix[m_globalid * 3 + m][k] = 1;
					}
					else
						m_bigmatrix[m_globalid * 3 + m][k] = 0;
				}
				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
			}
		}
		bdybase_count++;
	}
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(0);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize * 3; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			//if (fabs(m_bigmatrix[i][j]) > 1.0e-9)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{

			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num);
	}

}

//void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_mutires(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
//	std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
//	std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
//	std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick,
//	std::vector<std::vector<V3f>>& centers)
//{
//	/*---------------初始化全局刚度矩阵和右端项-----------*/
//	int m_GlobalSize = m_solids.Global_size;
//	m_colids.resize(m_GlobalSize * 3);
//
//	// m_bigrhs.resize(m_GlobalSize * 3);
//	Matrix m_bigmatrix;
//	m_bigmatrix.resize(m_GlobalSize * 3, std::vector<double>((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2) * martix_size, 0));
//	// for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	// {
//		// for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
//		// {
//			// for (int k = 0; k < 3; k++)
//				// m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize((2 * 3 + 2)*(2 * 3 + 2)*(2 * 3 + 2) * martix_size, 0);
//		// }
//	// }
//
//	std::vector<double> m_Weight, m_IntePt;
//	int m_GaussNum = 4;
//
//	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4
//
//	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		if (i % 100 == 0)
//			std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
//		std::vector<double> m_uIntervel{ 0,0.5,1 }, m_vIntervel{ 0,0.5,1 }, m_wIntervel{ 0,0.5,1 };
//		// std::vector<double> m_uIntervel{ 0,0.25,0.5,0.75,1 }, m_vIntervel{ 0,0.25,0.5,0.75,1 }, m_wIntervel{ 0,0.25,0.5,0.75,1 };
//		int m_uLocation, m_vLocation, m_wLocation;
//		std::vector<double> m_u, m_v, m_w;
//		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
//		//m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);
//
//		std::vector<Matrix> local_matrix;
//		std::vector<std::vector<double>> local_bigrh;
//		// std::vector<std::vector<std::pair<int, int>>> m_bdyquick;
//		std::vector<V3f> center;
//
//		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
//		int m_BaseIndex;
//		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
//		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
//		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
//		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
//		/***************************构造刚度矩阵按照单元循环****************************************/
//		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
//		{
//			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
//				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
//					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
//			//确定u方向这些高斯点所在区间的位置
//			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);
//
//			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
//			{
//				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
//					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
//						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
//				//确定v方向这些高斯点所在区间的位置
//				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);
//
//				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
//				{
//					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
//						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
//							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
//					//确定w方向这些高斯点所在区间的位置
//					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);
//
//					/*-------初始化这部分矩阵的信息------*/
//					V3f m_center(0, 0, 0);
//					int center_num = 0;
//					std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
//					double u_ders[2], v_ders[2], w_ders[2];
//					std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
//					std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，
//
//					m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
//					Matrix m_LocalMatrix;
//					m_LocalMatrix.resize(m_MatrixSize * 3);
//					for (int k = 0; k < m_MatrixSize * 3; k++)
//					{
//						m_LocalMatrix[k].resize(m_MatrixSize * 3, 0);
//					}
//					std::vector<double> m_bigrh(m_MatrixSize * 3, 0);
//					/*--------------初始化结束---------------*/
//
//					m_det.clear(); m_InverseMatrix.clear();
//					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
//					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
//					std::vector<int> m_baseindex;
//					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
//						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
//							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
//							{
//								double m_tempvalue1 = 0, m_tempvalue2 = 0, m_tempvalue3 = 0;
//								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
//								m_baseindex.push_back(m_BaseIndex);
//								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
//								for (int m = 0; m < m_GaussNum; m++)
//								{
//									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
//									for (int l = 0; l < m_GaussNum; l++)
//									{
//										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
//										for (int k = 0; k < m_GaussNum; k++)
//										{
//											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);
//
//											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
//											center_num++;
//											m_center = m_center + gasv;
//											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
//											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
//											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
//											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
//											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
//											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
//											/*int m_index = m * m_GaussNum*m_GaussNum + l * m_GaussNum + k;*/
//											/*m_tempvalue1 += u_ders[0] * v_ders[0] * w_ders[0] * fx(i,m_BaseIndex)*
//												fabs(m_det[m_index])*m_Weight[k] * m_Weight[l] * m_Weight[m];
//											m_tempvalue2 += u_ders[0] * v_ders[0] * w_ders[0] * fy(i,m_BaseIndex)*
//												fabs(m_det[m_index])*m_Weight[k] * m_Weight[l] * m_Weight[m];
//											m_tempvalue3 += u_ders[0] * v_ders[0] * w_ders[0] * fz(i,m_BaseIndex)*
//												fabs(m_det[m_index])*m_Weight[k] * m_Weight[l] * m_Weight[m];*/
//										}
//									}
//								}
//								/*double temp = (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval])*(m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
//									*(m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;*/
//
//									//右端项直接通过处理刚度矩阵的时候进行设置 这边就不设置了 全部注释
//									// m_bigrh[3 * m_BaseIndex] = fx(i,m_BaseIndex);
//									// m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex]] = fx(i,m_BaseIndex);
//									// m_bigrh[3 * m_BaseIndex + 1] = fy(i,m_BaseIndex);
//									// m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 1] = fy(i,m_BaseIndex);
//									// m_bigrh[3 * m_BaseIndex + 2] = fz(i,m_BaseIndex);
//									// m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 2] = fz(i,m_BaseIndex);
//
//									//m_bigrh[3 * m_BaseIndex] += m_tempvalue1 * temp;
//									//m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue1 * temp) * (umin + (1 - umin)*pow(1, p));
//									//m_bigrh[3 * m_BaseIndex + 1] += m_tempvalue2 * temp;
//									//m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 1] += (m_tempvalue2 * temp) * (umin + (1 - umin)*pow(1, p));
//									//m_bigrh[3 * m_BaseIndex + 2] += m_tempvalue3 * temp;
//									//m_bigrhs[3 * m_solids.Local_to_Global[i][m_BaseIndex] + 2] += (m_tempvalue3 * temp) * (umin + (1 - umin)*pow(1, p));
//							}
//					local_bigrh.push_back(m_bigrh);
//					m_center = m_center / center_num;
//					center.push_back(m_center);
//
//					int m_ptNum = m_BaseDValueu[0].size();
//					for (int ii = 0; ii < m_baseindex.size(); ii++)
//					{
//						for (int j = 0; j < m_baseindex.size(); j++)
//						{
//							for (int k = 0; k < m_ptNum; k++)
//							{
//								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
//								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
//								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
//								int m_index;
//								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
//								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
//								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
//								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
//								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
//								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
//								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
//								double temp = fabs(m_det[k]) * m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1]
//									- m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
//									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
//								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3] += ((lbda + 2 * mu) * m_Rx * m_Qx + mu * (m_Ry * m_Qy + m_Rz * m_Qz)) * temp;
//								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 1] += (lbda * m_Rx * m_Qy + mu * m_Ry * m_Qx) * temp;
//								m_LocalMatrix[m_baseindex[ii] * 3][m_baseindex[j] * 3 + 2] += (lbda * m_Rx * m_Qz + mu * m_Rz * m_Qx) * temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3] += (lbda * m_Ry * m_Qx + mu * m_Rx * m_Qy) * temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 1] += ((lbda + 2 * mu) * m_Ry * m_Qy + mu * (m_Rx * m_Qx + m_Rz * m_Qz)) * temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 1][m_baseindex[j] * 3 + 2] += (lbda * m_Ry * m_Qz + mu * m_Rz * m_Qy) * temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3] += (lbda * m_Rz * m_Qx + mu * m_Rx * m_Qz) * temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 1] += (lbda * m_Rz * m_Qy + mu * m_Ry * m_Qz) * temp;
//								m_LocalMatrix[m_baseindex[ii] * 3 + 2][m_baseindex[j] * 3 + 2] += ((lbda + 2 * mu) * m_Rz * m_Qz + mu * (m_Rx * m_Qx + m_Ry * m_Qy)) * temp;
//							}
//						}
//					}
//					local_matrix.push_back(m_LocalMatrix);
//
//					//矩阵中的非零元素全部存起来
//					for (int ii = 0; ii < m_MatrixSize; ii++)
//					{
//						for (int j = 0; j < m_MatrixSize; j++)
//						{
//							for (int k = 0; k < 3; k++)
//							{
//								for (int l = 0; l < 3; l++)
//								{
//									if (fabs(local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l]) > 1.0e-7)
//									{
//										if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
//										{
//											m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
//											int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//											if (m_size - 1 >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//											{
//												m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(m_size, 0);
//											}
//											m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(1, p));
//										}
//										else
//										{
//											int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//											if (tmp_id >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//											{
//												m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(tmp_id + 1, 0);
//											}
//											m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrix[local_matrix.size() - 1][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(1, p));
//										}
//									}
//								}
//							}
//						}
//					}
//				}
//
//			}
//		}
//		local_matrixs.push_back(local_matrix);
//		local_bigrhs.push_back(local_bigrh);
//		centers.push_back(center);
//	}
//	// int checknum = 0;
//	// for(int i = 0;i<m_bigrhs.size();++i)
//		// if(m_bigrhs[i]<0)
//			// checknum++;
//	// std::cout<<"checknum: "<<checknum<<std::endl;
//	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
//	//处理边界条件
//	//处理刚度矩阵和右端项
//	std::cout << "Start processing boundary conditions" << std::endl;
//	m_bigrhs.clear();
//	m_bigrhs.resize(m_GlobalSize * 3, 0);
//	int bdybase_count = 0;
//	for (int i = 0; i < bdybaseid.size(); i++)
//	{
//		std::vector<std::pair<int, int>> temp_bdyquick;
//		if (i % 100 == 0)
//			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
//		int m_baseid = bdybaseid[i].local_id;
//		int m_solidid = bdybaseid[i].solid_id;
//		int m_globalid = bdybaseid[i].global_id;
//		for (int j = 0; j < m_GlobalSize * 3; j++)
//		{
//			for (int k = 0; k < m_colids[j].size(); k++)
//			{
//				if (bdybase_count < subnum)
//				{
//					if (m_colids[j][k] == m_globalid * 3 + forceorder)
//					{
//						temp_bdyquick.push_back(std::pair<int, int>(j, k));
//						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
//						// m_bigmatrix[j][k] = 0;
//						break;
//					}
//				}
//				else
//				{
//					for (int l = 0; l < 3; l++)
//					{
//						if (m_colids[j][k] == m_globalid * 3 + l)
//						{
//							temp_bdyquick.push_back(std::pair<int, int>(j, k));
//							// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//							m_bigmatrix[j][k] = 0;
//							break;
//						}
//					}
//				}
//			}
//		}
//		m_bdyquick.push_back(temp_bdyquick);
//		if (bdybase_count < subnum)
//		{
//			// for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
//			// {
//				// if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
//				// {
//					// m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
//				// }
//				// else
//					// m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
//			// }
//			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
//		}
//		else
//		{
//			for (int m = 0; m < 3; m++)
//			{
//				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				{
//					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
//					{
//						m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					}
//					else
//						m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				}
//				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			}
//		}
//		bdybase_count++;
//	}
//	int checknum = 0;
//	for (int i = 0; i < m_bigrhs.size(); ++i)
//		if (m_bigrhs[i] < 0)
//			checknum++;
//	std::cout << "checknum: " << checknum << std::endl;
//	std::cout << "Processing boundary conditions completely" << std::endl;
//	r_rowids.push_back(1);
//	int ele_num = 0;
//	for (int i = 0; i < m_GlobalSize * 3; i++)
//	{
//		// if (i % 100 == 0)
//			// std::cout << "Quickly K " << i + 1 << ":" << m_GlobalSize * 3 << std::endl;
//		std::map<int, double> ele_colid;
//		for (int j = 0; j < m_bigmatrix[i].size(); j++)
//		{
//			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
//			{
//				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
//			}
//		}
//		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
//		{
//			if ((*ele_colid_it).first >= i)
//			{
//				r_colids.push_back((*ele_colid_it).first + 1);
//				r_vals.push_back((*ele_colid_it).second);
//				ele_num++;
//			}
//		}
//		r_rowids.push_back(ele_num + 1);
//	}
//	std::cout << "end1" << std::endl;
//}
//
//void FunctionBSpline::Global_MatrixAndRhs_Sparse_liner_mutires_quick(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
//	std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
//	std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
//	std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick)
//{
//	int m_GlobalSize = m_solids.Global_size;
//	m_bigrhs.clear();
//	m_bigrhs.resize(m_GlobalSize * 3, 0);
//
//	Matrix m_bigmatrix;
//	m_bigmatrix.resize(m_GlobalSize * 3);
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
//		{
//			for (int k = 0; k < 3; k++)
//				m_bigmatrix[m_solids.Local_to_Global[i][j] * 3 + k].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2) * martix_size, 0);
//		}
//	}
//
//	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
//	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//	{
//		if (i % 100 == 0)
//			std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
//		int m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
//		int m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
//		int m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
//		int m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
//
//		for (int lm = 0; lm < local_matrixs[i].size(); ++lm)
//		{
//			//矩阵中的非零元素全部存起来
//			for (int ii = 0; ii < m_MatrixSize; ii++)
//			{
//				for (int j = 0; j < m_MatrixSize; j++)
//				{
//					for (int k = 0; k < 3; k++)
//					{
//						for (int l = 0; l < 3; l++)
//						{
//							if (fabs(local_matrixs[i][lm][ii * 3 + k][j * 3 + l]) > 1.0e-7)
//							{
//								if (!InVecID(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k])) //全矩阵存储
//								{
//									m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].push_back(m_solids.Local_to_Global[i][j] * 3 + l);
//									int m_size = m_colids[m_solids.Local_to_Global[i][ii] * 3 + k].size();
//									if (m_size - 1 >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//									{
//										m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(m_size, 0);
//									}
//									m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][m_size - 1] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution_mutires[i][lm], p));
//								}
//								else
//								{
//									int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j] * 3 + l, m_colids[m_solids.Local_to_Global[i][ii] * 3 + k]);
//									if (tmp_id >= m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].size())
//									{
//										m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k].resize(tmp_id + 1, 0);
//									}
//									m_bigmatrix[m_solids.Local_to_Global[i][ii] * 3 + k][tmp_id] += local_matrixs[i][lm][ii * 3 + k][j * 3 + l] * (umin + (1 - umin) * pow(material_distribution_mutires[i][lm], p));
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//
//	}
//	// for (int i = 0; i < local_bigrhs.size(); i++)
//	// {
//		// for (int lb = 0; lb < local_bigrhs[i].size(); ++lb)
//		// {
//			// for (int k = 0; k < 3; k++)
//				// for (int j = 0; j < local_bigrhs[i].size() / 3; j++)
//					// m_bigrhs[m_solids.Local_to_Global[i][j] * 3 + k] += local_bigrhs[i][lb][j * 3 + k] /** (umin + (1 - umin)*pow(material_distribution_mutires[i][lb], p))*/;
//		// }
//	// }
//	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
//	//处理边界条件
//	//处理刚度矩阵和右端项
//	std::cout << "Start processing boundary conditions" << std::endl;
//	int bdybase_count = 0;
//	for (int i = 0; i < bdybaseid.size(); i++)
//	{
//		if (i % 1000 == 0)
//			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
//		int m_baseid = bdybaseid[i].local_id;
//		int m_solidid = bdybaseid[i].solid_id;
//		int m_globalid = bdybaseid[i].global_id;
//		std::vector<std::pair<int, int>> temp_bdyquick = m_bdyquick[i];
//		if (bdybase_count < subnum)
//		{
//			// for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			// {
//				// int j = (*temp_bdyquick_it).first;
//				// int k = (*temp_bdyquick_it).second;
//				// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + forceorder] * m_bigmatrix[j][k];
//				// m_bigmatrix[j][k] = 0;
//				// for (int l = 0; l < 3; l++)
//				// {
//					// if ((m_colids[j][k]) == m_globalid * 3 + l)
//					// {
//						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//						// break;
//					// }
//				// }
//				// m_bigmatrix[j][k] = 0;
//			// }
//		}
//		else
//		{
//			for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
//			{
//				int j = (*temp_bdyquick_it).first;
//				int k = (*temp_bdyquick_it).second;
//				for (int l = 0; l < 3; l++)
//				{
//					if ((m_colids[j][k]) == m_globalid * 3 + l)
//					{
//						// m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid * 3 + l] * m_bigmatrix[j][k];
//						break;
//					}
//				}
//				m_bigmatrix[j][k] = 0;
//			}
//		}
//		if (bdybase_count < subnum)
//		{
//			// for (int k = 0; k < m_colids[m_globalid * 3 + forceorder].size(); k++)
//			// {
//				// if (m_colids[m_globalid * 3 + forceorder][k] == m_globalid * 3 + forceorder)
//				// {
//					// m_bigmatrix[m_globalid * 3 + forceorder][k] = 1;
//				// }
//				// else
//					// m_bigmatrix[m_globalid * 3 + forceorder][k] = 0;
//			// }
//			m_bigrhs[m_globalid * 3 + forceorder] = m_bdycoef[m_solidid][m_baseid * 3 + forceorder];
//		}
//		else
//		{
//			for (int m = 0; m < 3; m++)
//			{
//				for (int k = 0; k < m_colids[m_globalid * 3 + m].size(); k++)
//				{
//					if (m_colids[m_globalid * 3 + m][k] == m_globalid * 3 || m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 1
//						|| m_colids[m_globalid * 3 + m][k] == m_globalid * 3 + 2)
//					{
//						m_bigmatrix[m_globalid * 3 + m][k] = 1;
//					}
//					else
//						m_bigmatrix[m_globalid * 3 + m][k] = 0;
//				}
//				m_bigrhs[m_globalid * 3 + m] = m_bdycoef[m_solidid][m_baseid * 3 + m];
//			}
//		}
//		bdybase_count++;
//	}
//	std::cout << "Processing boundary conditions completely" << std::endl;
//	r_rowids.push_back(1);
//	int ele_num = 0;
//	for (int i = 0; i < m_GlobalSize * 3; i++)
//	{
//		std::map<int, double> ele_colid;
//		for (int j = 0; j < m_bigmatrix[i].size(); j++)
//		{
//			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
//			{
//				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
//			}
//		}
//		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
//		{
//			if ((*ele_colid_it).first >= i)
//			{
//				r_colids.push_back((*ele_colid_it).first + 1);
//				r_vals.push_back((*ele_colid_it).second);
//				ele_num++;
//			}
//		}
//		r_rowids.push_back(ele_num + 1);
//	}
//}

void FunctionBSpline::Global_MatrixAndRhs_Sparse_hot_mutires(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
	std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
	std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick,
	std::vector<std::vector<V3f>>& centers)
{
	/*---------------初始化全局刚度矩阵和右端项-----------*/
	int m_GlobalSize = m_solids.Global_size;
	m_colids.resize(m_GlobalSize);

	m_bigrhs.resize(m_GlobalSize);
	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			m_bigmatrix[m_solids.Local_to_Global[i][j]].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2) * martix_size, 0);
		}
	}

	std::vector<double> m_Weight, m_IntePt;
	int m_GaussNum = 4;

	GaussIntegalPoint(m_GaussNum, m_Weight, m_IntePt);//m_Weight和m_IntePt的size是4

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		std::vector<double> m_uIntervel{ 0,0.5,1 }, m_vIntervel{ 0,0.5,1 }, m_wIntervel{ 0,0.5,1 };
		// std::vector<double> m_uIntervel{ 0,0.25,0.5,0.75,1 }, m_vIntervel{ 0,0.25,0.5,0.75,1 }, m_wIntervel{ 0,0.25,0.5,0.75,1 };
		int m_uLocation, m_vLocation, m_wLocation;
		std::vector<double> m_u, m_v, m_w;
		m_u.resize(m_GaussNum); m_v.resize(m_GaussNum); m_w.resize(m_GaussNum);
		//m_solids.IntegralArea(i, m_uIntervel, m_vIntervel, m_wIntervel);

		std::vector<Matrix> local_matrix;
		std::vector<std::vector<double>> local_bigrh;
		std::vector<std::vector<std::pair<int, int>>> m_bdyquick;
		std::vector<V3f> center;

		int m_uPtnum, m_vPtnum, m_wPtnum, m_MatrixSize;
		int m_BaseIndex;
		m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;
		/***************************构造刚度矩阵按照单元循环****************************************/
		for (int u_interval = 0; u_interval < m_uIntervel.size() - 1; u_interval++)
		{
			for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
				m_u[gaussnum] = (m_uIntervel[u_interval] + m_uIntervel[u_interval + 1] + (-m_uIntervel[u_interval]
					+ m_uIntervel[u_interval + 1]) * m_IntePt[gaussnum]) * 0.5;
			//确定u方向这些高斯点所在区间的位置
			m_uLocation = m_solids.X_LOCATION(i, 0, m_u[0]);

			for (int v_intervel = 0; v_intervel < m_vIntervel.size() - 1; v_intervel++)
			{
				for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
					m_v[gaussnum] = (m_vIntervel[v_intervel] + m_vIntervel[v_intervel + 1] + (-m_vIntervel[v_intervel]
						+ m_vIntervel[v_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
				//确定v方向这些高斯点所在区间的位置
				m_vLocation = m_solids.X_LOCATION(i, 1, m_v[0]);

				for (int w_intervel = 0; w_intervel < m_wIntervel.size() - 1; w_intervel++)
				{
					for (int gaussnum = 0; gaussnum < m_GaussNum; gaussnum++)
						m_w[gaussnum] = (m_wIntervel[w_intervel] + m_wIntervel[w_intervel + 1] + (-m_wIntervel[w_intervel]
							+ m_wIntervel[w_intervel + 1]) * m_IntePt[gaussnum]) * 0.5;
					//确定w方向这些高斯点所在区间的位置
					m_wLocation = m_solids.X_LOCATION(i, 2, m_w[0]);

					/*-------初始化这部分矩阵的信息------*/
					V3f m_center(0, 0, 0);
					int center_num = 0;
					std::vector<Matrix> m_InverseMatrix;//每个基对应一组在各个高斯点上的雅可比逆矩阵
					double u_ders[2], v_ders[2], w_ders[2];
					std::vector<double> m_det;//记录每个单元上高斯点上的jacobi行列式
					std::vector<std::vector<double>> m_BaseDValueu, m_BaseDValuev, m_BaseDValuew;//存储样条基关于u和v在高斯点的导数，

					m_BaseDValueu.resize(m_MatrixSize); m_BaseDValuev.resize(m_MatrixSize); m_BaseDValuew.resize(m_MatrixSize);
					Matrix m_LocalMatrix;
					m_LocalMatrix.resize(m_MatrixSize);
					for (int k = 0; k < m_MatrixSize; k++)
					{
						m_LocalMatrix[k].resize(m_MatrixSize, 0);
					}
					std::vector<double> m_bigrh(m_MatrixSize, 0);
					/*--------------初始化结束---------------*/

					m_det.clear(); m_InverseMatrix.clear();
					/****m_det存储当前单元中高斯点处的Jacobi矩阵的行列式，m_InverseMatrix存储高斯点处的Jacobi矩阵的转置的逆****/
					m_solids.JacobiMatrix_new(i, m_u, m_v, m_w, m_det, m_InverseMatrix);
					std::vector<int> m_baseindex;
					for (int ii = m_uLocation - m_solids.Local_Order[i][0]; ii <= m_uLocation; ii++)
						for (int j = m_vLocation - m_solids.Local_Order[i][1]; j <= m_vLocation; j++)
							for (int ww = m_wLocation - m_solids.Local_Order[i][2]; ww <= m_wLocation; ww++)
							{
								double m_tempvalue = 0;
								m_BaseIndex = m_solids.Control_Point[i][ww][j][ii];
								m_baseindex.push_back(m_BaseIndex);
								m_BaseDValueu[m_BaseIndex].clear(); m_BaseDValuev[m_BaseIndex].clear(); m_BaseDValuew[m_BaseIndex].clear();
								for (int m = 0; m < m_GaussNum; m++)
								{
									m_solids.DersOneBasisFun(i, 2, ww, m_w[m], 1, w_ders);
									for (int l = 0; l < m_GaussNum; l++)
									{
										m_solids.DersOneBasisFun(i, 1, j, m_v[l], 1, v_ders);
										for (int k = 0; k < m_GaussNum; k++)
										{
											m_solids.DersOneBasisFun(i, 0, ii, m_u[k], 1, u_ders);

											V3f gasv = m_solids.Geometry_new(i, m_u[k], m_v[l], m_w[m]);
											center_num++;
											m_center = m_center + gasv;
											double m_temp = u_ders[1] * v_ders[0] * w_ders[0];
											m_BaseDValueu[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[1] * w_ders[0];
											m_BaseDValuev[m_BaseIndex].push_back(m_temp);
											m_temp = u_ders[0] * v_ders[0] * w_ders[1];
											m_BaseDValuew[m_BaseIndex].push_back(m_temp);
											int m_index = m * m_GaussNum * m_GaussNum + l * m_GaussNum + k;
											m_tempvalue += u_ders[0] * v_ders[0] * w_ders[0] * RHS(gasv.x, gasv.y, gasv.z) *
												fabs(m_det[m_index]) * m_Weight[k] * m_Weight[l] * m_Weight[m];
										}
									}
								}
								double temp = (m_uIntervel[u_interval + 1] - m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								m_bigrh[m_BaseIndex] += m_tempvalue * temp;
								m_bigrhs[m_solids.Local_to_Global[i][m_BaseIndex]] += (m_tempvalue * temp) * (umin + (1 - umin) * pow(1, p));
							}
					local_bigrh.push_back(m_bigrh);
					m_center = m_center / center_num;
					center.push_back(m_center);

					int m_ptNum = m_BaseDValueu[0].size();
					for (int ii = 0; ii < m_baseindex.size(); ii++)
					{
						for (int j = 0; j < m_baseindex.size(); j++)
						{
							for (int k = 0; k < m_ptNum; k++)
							{
								int m_uptIndex = (k % (m_GaussNum * m_GaussNum)) % m_GaussNum;
								int m_vptIndex = (k % (m_GaussNum * m_GaussNum)) / m_GaussNum;
								int m_wptIndex = k / (m_GaussNum * m_GaussNum);
								int m_index;
								double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;               //存储基的梯度在区域上变换后的量。
								m_Rx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Ry = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Rz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[ii]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[ii]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[ii]][k];
								m_Qx = m_InverseMatrix[k][0][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][0][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][0][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qy = m_InverseMatrix[k][1][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][1][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][1][2] * m_BaseDValuew[m_baseindex[j]][k];
								m_Qz = m_InverseMatrix[k][2][0] * m_BaseDValueu[m_baseindex[j]][k] + m_InverseMatrix[k][2][1] * m_BaseDValuev[m_baseindex[j]][k] + m_InverseMatrix[k][2][2] * m_BaseDValuew[m_baseindex[j]][k];
								double temp = fabs(m_det[k]) * m_Weight[m_uptIndex] * m_Weight[m_vptIndex] * m_Weight[m_wptIndex] * (m_uIntervel[u_interval + 1]
									- m_uIntervel[u_interval]) * (m_vIntervel[v_intervel + 1] - m_vIntervel[v_intervel])
									* (m_wIntervel[w_intervel + 1] - m_wIntervel[w_intervel]) / 8;
								m_LocalMatrix[m_baseindex[ii]][m_baseindex[j]] += (m_Rx * m_Qx + m_Ry * m_Qy + m_Rz * m_Qz) * fabs(m_det[k]) * temp;

							}
						}
					}
					local_matrix.push_back(m_LocalMatrix);

					//矩阵中的非零元素全部存起来
					for (int ii = 0; ii < m_MatrixSize; ii++)
					{
						for (int j = 0; j < m_MatrixSize; j++)
						{
							if (fabs(local_matrix[local_matrix.size() - 1][ii][j]) > 1.0e-7)
							{
								if (!InVecID(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]])) //全矩阵存储
								{
									m_colids[m_solids.Local_to_Global[i][ii]].push_back(m_solids.Local_to_Global[i][j]);
									int m_size = m_colids[m_solids.Local_to_Global[i][ii]].size();
									m_bigmatrix[m_solids.Local_to_Global[i][ii]][m_size - 1] += local_matrix[local_matrix.size() - 1][ii][j] * (umin + (1 - umin) * pow(1, p));
								}
								else
								{
									int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]]);
									m_bigmatrix[m_solids.Local_to_Global[i][ii]][tmp_id] += local_matrix[local_matrix.size() - 1][ii][j] * (umin + (1 - umin) * pow(1, p));
								}
							}
						}
					}
				}

			}
		}
		local_matrixs.push_back(local_matrix);
		local_bigrhs.push_back(local_bigrh);
		centers.push_back(center);
	}
	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		std::vector<std::pair<int, int>> temp_bdyquick;
		if (i % 1000 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		for (int j = 0; j < m_GlobalSize; j++)
		{
			for (int k = 0; k < m_colids[j].size(); k++)
			{
				if ((m_colids[j][k]) == m_globalid)
				{
					temp_bdyquick.push_back(std::pair<int, int>(j, k));
					m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid] * m_bigmatrix[j][k];
					m_bigmatrix[j][k] = 0;
					break;
				}
			}
		}
		m_bdyquick.push_back(temp_bdyquick);
		for (int k = 0; k < m_colids[m_globalid].size(); k++)
		{
			if (m_colids[m_globalid][k] == m_globalid)
			{
				m_bigmatrix[m_globalid][k] = 1;
			}
			else
				m_bigmatrix[m_globalid][k] = 0;
		}
		m_bigrhs[m_globalid] = m_bdycoef[m_solidid][m_baseid];
		bdybase_count++;
	}
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(1);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first + 1);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num + 1);
	}
}

void FunctionBSpline::Global_MatrixAndRhs_Sparse_hot_mutires_quick(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
	std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
	std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
	std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick)
{
	int m_GlobalSize = m_solids.Global_size;
	m_bigrhs.clear();
	m_bigrhs.resize(m_GlobalSize, 0);

	Matrix m_bigmatrix;
	m_bigmatrix.resize(m_GlobalSize);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
		{
			m_bigmatrix[m_solids.Local_to_Global[i][j]].resize((2 * 3 + 2) * (2 * 3 + 2) * (2 * 3 + 2) * martix_size, 0);
		}
	}

	std::cout << "Start constructing stiffness matrix and the right end term" << std::endl;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		if (i % 100 == 0)
			std::cout << "Quickly schedule of stiffness matrix " << i + 1 << ":" << m_solids.Local_V3f.size() << std::endl;
		int m_uPtnum = m_solids.Local_Node[i][0].size() - m_solids.Local_Order[i][0] - 1;
		int m_vPtnum = m_solids.Local_Node[i][1].size() - m_solids.Local_Order[i][1] - 1;
		int m_wPtnum = m_solids.Local_Node[i][2].size() - m_solids.Local_Order[i][2] - 1;
		int m_MatrixSize = m_uPtnum * m_vPtnum * m_wPtnum;

		for (int lm = 0; lm < local_matrixs[i].size(); ++lm)
		{
			//矩阵中的非零元素全部存起来
			for (int ii = 0; ii < m_MatrixSize; ii++)
			{
				for (int j = 0; j < m_MatrixSize; j++)
				{

					if (fabs(local_matrixs[i][lm][ii][j]) > 1.0e-7)
					{
						if (!InVecID(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]])) //全矩阵存储
						{
							m_colids[m_solids.Local_to_Global[i][ii]].push_back(m_solids.Local_to_Global[i][j]);
							int m_size = m_colids[m_solids.Local_to_Global[i][ii]].size();
							m_bigmatrix[m_solids.Local_to_Global[i][ii]][m_size - 1] += local_matrixs[i][lm][ii][j] * (umin + (1 - umin) * pow(material_distribution_mutires[i][lm], p));
						}
						else
						{
							int tmp_id = InVecIndex(m_solids.Local_to_Global[i][j], m_colids[m_solids.Local_to_Global[i][ii]]);
							m_bigmatrix[m_solids.Local_to_Global[i][ii]][tmp_id] += local_matrixs[i][lm][ii][j] * (umin + (1 - umin) * pow(material_distribution_mutires[i][lm], p));
						}
					}

				}
			}
		}

	}
	for (int i = 0; i < local_bigrhs.size(); i++)
	{
		for (int lb = 0; lb < local_bigrhs[i].size(); ++lb)
		{
			for (int j = 0; j < local_bigrhs[i].size(); j++)
				m_bigrhs[m_solids.Local_to_Global[i][j]] += local_bigrhs[i][lb][j] /** (umin + (1 - umin)*pow(material_distribution_mutires[i][lb], p))*/;
		}
	}
	std::cout << "Constructing stiffness matrix and the right end term completely" << std::endl;
	//处理边界条件
	//处理刚度矩阵和右端项
	std::cout << "Start processing boundary conditions" << std::endl;
	int bdybase_count = 0;
	for (int i = 0; i < bdybaseid.size(); i++)
	{
		std::vector<std::pair<int, int>> temp_bdyquick = m_bdyquick[i];
		int m_baseid = bdybaseid[i].local_id;
		int m_solidid = bdybaseid[i].solid_id;
		int m_globalid = bdybaseid[i].global_id;
		if (i % 1000 == 0)
			std::cout << "Schedule of processing boundary conditions " << i + 1 << ":" << bdybaseid.size() << std::endl;
		for (auto temp_bdyquick_it = temp_bdyquick.begin(); temp_bdyquick_it != temp_bdyquick.end(); temp_bdyquick_it++)
		{
			int j = (*temp_bdyquick_it).first;
			int k = (*temp_bdyquick_it).second;
			m_bigrhs[j] -= m_bdycoef[m_solidid][m_baseid] * m_bigmatrix[j][k];
			m_bigmatrix[j][k] = 0;
		}
		m_bdyquick.push_back(temp_bdyquick);
		for (int k = 0; k < m_colids[m_globalid].size(); k++)
		{
			if (m_colids[m_globalid][k] == m_globalid)
			{
				m_bigmatrix[m_globalid][k] = 1;
			}
			else
				m_bigmatrix[m_globalid][k] = 0;
		}
		m_bigrhs[m_globalid] = m_bdycoef[m_solidid][m_baseid];
		bdybase_count++;
	}
	std::cout << "Processing boundary conditions completely" << std::endl;
	r_rowids.push_back(1);
	int ele_num = 0;
	for (int i = 0; i < m_GlobalSize; i++)
	{
		std::map<int, double> ele_colid;
		for (int j = 0; j < m_bigmatrix[i].size(); j++)
		{
			if (fabs(m_bigmatrix[i][j]) > 1.0e-7)
			{
				ele_colid[m_colids[i][j]] = m_bigmatrix[i][j];
			}
		}
		for (auto ele_colid_it = ele_colid.begin(); ele_colid_it != ele_colid.end(); ele_colid_it++)
		{
			if ((*ele_colid_it).first >= i)
			{
				r_colids.push_back((*ele_colid_it).first + 1);
				r_vals.push_back((*ele_colid_it).second);
				ele_num++;
			}
		}
		r_rowids.push_back(ele_num + 1);
	}
}
