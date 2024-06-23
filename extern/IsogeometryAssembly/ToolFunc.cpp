#include"ToolFunc.h"
void ToolFunc::transformGradients(Eigen::MatrixXd jac, int q, const Eigen::MatrixXd& allGrads, Eigen::MatrixXd& trfGradsK, int dim)
{
	if (allGrads.size() % dim != 0) cout << "Invalid size of gradient matrix";
	const int numGrads = allGrads.rows() / dim;
	const Eigen::Map<const Eigen::MatrixXd> grads_k(allGrads.col(q).data(), dim, numGrads);
	trfGradsK.noalias() = jac.inverse() * grads_k;
}

//基函数梯度+原值
void ToolFunc::transformG(std::vector<CVector> jac, int q, const Eigen::MatrixXd& allGrads, Eigen::MatrixXd& trfGradsK, Eigen::MatrixXd& G_hat, int dim)
{
	if (allGrads.size() % dim != 0) cout << "Invalid size of gradient matrix";

	const int numGrads = allGrads.rows() / dim;
	const Eigen::Map<const Eigen::MatrixXd> grads_k(allGrads.col(q).data(), dim, numGrads);
	Eigen::Matrix2d jacobian;
	G_hat = grads_k;
	jacobian << jac[0].getX(), jac[0].getY(), jac[1].getX(), jac[1].getY();
	trfGradsK.noalias() = jacobian.inverse() * grads_k;
}

//高斯点获取
void ToolFunc::gaussIntegalPoint(int N, double lowpoint, double uppoint, std::vector<double>& point, std::vector<double>& weight)
{
	if (N == 2) {
		point = { -0.5773502691896250, 0.5773502691896250 };
		weight = { 1,1 };
	}
	else if (N == 3) {
		point = { -0.7745966692414830, 0.0000000000, 0.7745966692414830 };
		weight = { 0.5555555555555550, 0.88888888880, 0.5555555555555550 };
	}
	else if (N == 4) {
		point = { -0.8611363115940520, -0.3399810435848560, 0.3399810435848560, 0.8611363115940520 };
		weight = { 0.3478548451374530, 0.6521451548625460, 0.6521451548625460, 0.3478548451374530 };
	}
	else
		std::cout << "高斯点数量超出范围" << "\n";
	// 积分：0 - 1 f(t)dt
	// 令 x = t * (b - a) / 2 + (b + a) / 2
	// 积分：a - b f(x)dx = a - b f(t * (b - a) / 2 + (b + a) / 2)dt * (b - a) / 2
	double hrop = fabs(uppoint - lowpoint) < 1e-8 ? double(0.5) : (uppoint - lowpoint) / 2;
	for (int i = 0; i < N; ++i) {
		point[i] = (uppoint + lowpoint) / 2 - point[i] * (lowpoint - uppoint) / 2;
		weight[i] = weight[i] * hrop;
	}
}

//材料矩阵
void ToolFunc::setC(Eigen::MatrixXd& C, const Eigen::MatrixXd& R, double lambda, double mu)
{
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

//B
void ToolFunc::setB(Eigen::MatrixXd& B, const Eigen::MatrixXd& F, const Eigen::VectorXd& bGrad)
{
	short a_dim = F.cols();
	short dimTensor = a_dim * (a_dim + 1) / 2;
	B.resize(dimTensor, a_dim);
	//B.setZero(dimTensor, a_dim);

	//B(0, 0) = bGrad(0);
	//B(1, 1) = bGrad(1);
	//B(2, 2) = bGrad(2);
	//B(3, 0) = bGrad(1); B(3, 1) = bGrad(0);
	//B(4, 1) = bGrad(2); B(4, 2) = bGrad(1);
	//B(5, 0) = bGrad(2); B(5, 2) = bGrad(0);

	for (short j = 0; j < a_dim; ++j)
	{
		for (short i = 0; i < a_dim; ++i)
			B(i, j) = F(j, i) * bGrad(i);

		if (a_dim == 2)
			B(2, j) = F(j, 0) * bGrad(1) + F(j, 1) * bGrad(0);

		if (a_dim == 3)
			for (short i = 0; i < a_dim; ++i)
			{
				short k = (i + 1) % a_dim;
				B(i + a_dim, j) = F(j, i) * bGrad(k) + F(j, k) * bGrad(i);
			}
	}
}

//删除矩阵第i,j后计算determine
void ToolFunc::firstMinor(Eigen::MatrixXd& input, int i, int j, Eigen::MatrixXd& result) {
	int r = input.rows(), c = input.cols();
	if (i < r && j < c) {
		result.block(0, 0, i, j) = input.block(0, 0, i, j);
		result.block(0, j, i, c - j - 1) = input.block(0, j + 1, i, c - j - 1);
		result.block(i, 0, r - i - 1, j) = input.block(i + 1, 0, r - i - 1, j);
		result.block(i, j, r - i - 1, c - j - 1) = input.block(i + 1, j + 1, r - i - 1, c - j - 1);
	}
}

//右端项的雅克比计算
void ToolFunc::outerNormal(const Eigen::MatrixXd& jac, boxSide s, Eigen::VectorXd& result) {
	int m_orientation = jac.determinant() >= 0 ? 1 : -1;
	const double sgn = s.sideOrientation() * m_orientation;
	const int dir = s.direction();
	result.resize(3);
	Eigen::MatrixXd JK = jac.transpose();
	double alt_sgn = sgn;
	Eigen::MatrixXd minor;
	minor.resize(JK.rows() - 1, JK.cols() - 1);
	for (short i = 0; i != 3; ++i) // for all components of the normal
	{
		firstMinor(JK, i, dir, minor);
		result[i] = alt_sgn * minor.determinant();
		alt_sgn = -alt_sgn;
	}
}

void ToolFunc::outCsv(const Eigen::MatrixXd& info, string filename)
{
	std::ofstream myfile;
	myfile.open("D:/Resource/FileDataInput/threeD/"+filename+".csv");
	for (int i = 0; i < info.rows(); ++i) {
		for (int j = 0; j < info.cols(); ++j)
			myfile << info(i, j) << ",";
		myfile << "\n";
	}
	myfile.close();
}

void ToolFunc::readFileXml(const string filename, std::vector<CBSplineSolid>& geometry)
{
	//const char* filepath = (char*)filename.data();
	//TiXmlDocument doc(filepath);
	//bool loadOkay = doc.LoadFile();

	//if (!loadOkay) {
	//	printf("Could not load test file %s. Error='%s'. Exiting.\n", filepath, doc.ErrorDesc());
	//	exit(1);
	//}

	//// get dom root of 'phonebookdata.xml', here root should be 'phonebook'.  
	//TiXmlElement* root = doc.RootElement();
	//// trace every items below root.  
	//for (TiXmlNode* item = root->FirstChild("Geometry");item;item = item->NextSibling("Geometry")) {
	//	
	//	

	//	int dim = 3;
	//	TiXmlElement* basis = item->FirstChildElement("Basis");
	//	//dim = stoi(basis->Attribute("parDim"));
	//	vector<int> degree;
	//	vector<vector<int>> mint(dim);
	//	vector<vector<double>> uvw_knot(dim);
	//	int iter = 0;

	//	for (TiXmlNode* basis_item = basis->FirstChild("Basis"); basis_item; basis_item = basis_item->NextSibling("Basis")) {
	//		
	//		TiXmlElement* knot = basis_item->FirstChild("KnotVector")->ToElement();
	//		degree.emplace_back(stoi(knot->Attribute("degree")));
	//		string value = knot->GetText();
	//		std::regex pattern(" ");
	//		vector<string> result(sregex_token_iterator(value.begin(),value.end(),pattern,-1),std::sregex_token_iterator());
	//		for (int be = degree.back(); be < result.size()-degree.back(); ++be) {
	//			uvw_knot[iter].emplace_back(stod(result[be]));
	//		}
	//		for (int mut = 0; mut < result.size() - 2 * degree.back(); ++mut) mint[iter].emplace_back(1);
	//		mint[iter][0] = degree.back() + 1;
	//		mint[iter].back() = degree.back() + 1;
	//		iter++;
	//	}
	//	TiXmlNode* coefs = item->FirstChild("coefs");
	//	string value = coefs->ToElement()->GetText();
	//	vector<vector<vector<CPoint>>> contropoints;
	//	std::regex pattern(" ");
	//	vector<string> result(sregex_token_iterator(value.begin(), value.end(), pattern, -1), std::sregex_token_iterator());
	//	int uCoef = uvw_knot[0].size() + degree[0] - 1;
	//	int vCoef = uvw_knot[1].size() + degree[1] - 1;
	//	int wCoef = uvw_knot[2].size() + degree[2] - 1;
	//	if (result.size() != uCoef * vCoef * wCoef*3) {
	//		cout << "Input file Coefs Error" << endl;
	//		exit(1);
	//	}
	//	contropoints.resize(uCoef, vector<vector<CPoint>>(vCoef, vector<CPoint>(wCoef)));
	//	iter = 0;
	//	for (int k = 0; k < wCoef; ++k) {
	//		for (int j = 0; j < vCoef; ++j) {
	//			for (int i = 0; i < uCoef; ++i) {
	//				contropoints[i][j][k].setX(stod(result[iter * 3]));
	//				contropoints[i][j][k].setY(stod(result[iter * 3+1]));
	//				contropoints[i][j][k].setZ(stod(result[iter * 3+2]));
	//				iter++;
	//			}
	//		}
	//	}

	//	CBSplineSolid solid(degree[0], degree[1], degree[2], contropoints, uvw_knot[0], uvw_knot[1], uvw_knot[2], mint[0], mint[1], mint[2]);
	//	geometry.push_back(solid);
	//}
}

void ToolFunc::voigtStressToK(const Eigen::MatrixXd& S, Eigen::VectorXd& Svec)
{
	size_t dim = S.cols();
	size_t dimTensor = dim * (dim + 1) / 2;
	Svec.resize(dimTensor);
	if (dim == 2)
	{
		Svec(0) = S(0, 0);
		Svec(1) = S(1, 1);
		Svec(2) = S(0, 1);
	}
	if (dim == 3)
	{
		Svec(0) = S(0, 0);
		Svec(1) = S(1, 1);
		Svec(2) = S(2, 2);
		Svec(3) = S(0, 1);
		Svec(4) = S(1, 2);
		Svec(5) = S(0, 2);
	}
}