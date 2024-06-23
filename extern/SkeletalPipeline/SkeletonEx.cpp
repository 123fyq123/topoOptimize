#include "SkeletonEx.h"


double get180RotateAngle(double x1, double y1, double x2, double y2)
{
	const double epsilon = 1.0e-6;

	double dist, dot, degree, angle;

	// normalize
	dist = sqrt(x1 * x1 + y1 * y1);
	x1 /= dist;
	y1 /= dist;
	dist = sqrt(x2 * x2 + y2 * y2);
	x2 /= dist;
	y2 /= dist;
	// dot product
	dot = x1 * x2 + y1 * y2;
	if (std::fabs(dot - 1.0) <= epsilon)
		angle = 0.0;
	else if (std::fabs(dot + 1.0) <= epsilon)
		angle = M_PI;
	else {
		double cross;

		angle = std::acos(dot);
		//cross product
		cross = x1 * y2 - x2 * y1;
		// vector p2 is clockwise from vector p1
		// with respect to the origin (0.0)
		if (cross < 0) {
			angle = 2 * M_PI - angle;
		}
	}
	//std::cout << angle;
	

	return angle;
}


bool compare_pd(const std::pair<Eigen::Vector2i, double> &a, const std::pair<Eigen::Vector2i, double> &b)
{
	
	if (a.second == b.second)
	{
		if (a.first(0) == b.first(0))
		{
			return std::abs(a.first(1) - 360.5) < std::abs(b.first(1) - 360.5);
		}
		return a.first(0) < b.first(0);
	}
	return a.second < b.second;
	
	/*
	if (a.second == b.second)
	{
		if (a.first(1) == b.first(1))
		{
			return std::abs(a.first(0) - 80.5) >= std::abs(b.first(0) - 80.5);
			//return a.first(0)>b.first(0);
		}
		else
			return a.first(1) < b.first(1);
	}else
		return a.second < b.second;
	*/
}
void SkeletonEx::ReadTopMesh(std::string mesh_path, bool need_c)
{
	OpenMesh::IO::Options opt;
	opt += OpenMesh::IO::Options::FaceColor;
	mesh.request_face_colors();
	if (!OpenMesh::IO::read_mesh(mesh, mesh_path, opt))
	{
		readed = false;
		std::cerr << "Error loading mesh from file: " << mesh_path << std::endl;
		return;
	}
	if (!opt.check(OpenMesh::IO::Options::FaceColor))
	{
		std::cerr << "Error reading mesh from file, it doesn't has FaceColor properties!";
	}
	readed = true;
	//Assembly the matrix within the range of width*height
	uint32_t i = 0;
	for (auto f_i : mesh.faces())
	{
		indexToXYZ(i, grid_xy[i]);
		grid_index(grid_xy[i][0], grid_xy[i][1]) = i;
		++i;
	}

	std::cout << "mesh readed\nwidth:" << width << " height:" << height << "\n";
}

void SkeletonEx::InitBoundary(bool is_8_nei)
{

	OpenMesh::Vec3uc color_bdy(0, 0, 255);
	OpenMesh::Vec3uc color_inner(255, 0, 0);

	//auto propsBdy = OpenMesh::makeTemporaryProperty<OpenMesh::FaceHandle,int>(mesh);
	auto propsBdy = OpenMesh::FProp<int>(mesh, bdy_flag.c_str());

	//TODO: verify the 8 nei_c and 4 nei_c
	if (is_8_nei)
	{ //find all inner grid by 8 neiboring grids
		auto propsVBdy = OpenMesh::VProp<int>(mesh);
		for (auto i : mesh.vertices())
		{
			propsVBdy[i] = 0;
		}
	

		for (const auto &face : mesh.faces())
		{
			if (mesh.color(face) == color_bdy)
			{
				propsBdy[face] = -1;
				for (const auto &v : mesh.fv_range(face))
				{
					propsVBdy[v] = 1;
				}
			}
			else
			{
				propsBdy[face] = 1;
			}
		}

		std::vector<OpenMesh::FaceHandle> updateElem;
		for (const auto &face : mesh.faces())
		{
			if (propsBdy[face] == -1)
				continue;
			for (const auto &v4 : mesh.fv_range(face))
			{
				if (propsVBdy[v4] == 1)
				{
					updateElem.push_back(face);
					break;
				}
			}
		}
		for (auto fi : updateElem)
		{
			propsBdy[fi] = 0;
		}
	}
	else
	{ //find all inner grid by 4 neiboring grids
		for (const auto &face : mesh.faces())
		{
			if (mesh.color(face) == color_bdy)
			{
				propsBdy[face] = -1;
			}
			else
			{
				propsBdy[face] = 1;
			}
		}

		std::vector<OpenMesh::FaceHandle> updateElem;
		for (const auto &face : mesh.faces())
		{
			if (propsBdy[face] == -1)
				continue;
			for (const auto &f_nei : mesh.ff_range(face))
			{
				if (propsBdy[f_nei] == -1)
				{
					updateElem.push_back(face);
					break;
				}
			}
		}
		for (auto fi : updateElem)
		{
			propsBdy[fi] = 0;
		}
	}

	//Assembly the matrix within the range of width*height
	uint32_t i = 0;
	for (uint32_t row = 0; row < height; ++row)
	{
		for (uint32_t col = 0; col < width; ++col)
		{
			//set scalar value
			grid_value(row, col) = propsBdy[OpenMesh::FaceHandle(grid_index(row, col))];
		}
	}
}

void SkeletonEx::ComDistanceField()
{
	DistanceMap dc(width, height, grid_value);
	dc.Generate_Bdy_inner(true, true);
	//std::cout << dc.grid2.dis.row(1);
	//dc.grid1.dis = -dc.grid1.dis;
	//dc.grid1.cp = dc.grid2.cp;
	dc.grid2.dis = -dc.grid1.dis;
	distance_map = dc.grid2;
	//std::cout << "calculate distance field\n";
	std::cout<<distance_map.dis.block(0,0,10,10);
}

void SkeletonEx::SampleSphere2D(Eigen::MatrixXd &samples)
{
	double alpha = 2.00 * M_PI / (double)NUMBER_SAMPLES;
	for (uint32_t i = 0; i < NUMBER_SAMPLES; ++i)
	{
		samples(0, i) = std::cos(alpha * i) * epsilon;
		samples(1, i) = std::sin(alpha * i) * epsilon;
	}
}

void SkeletonEx::ComGradField()
{
	//calculate the outter flux field
	Eigen::MatrixXd samples;
	samples.resize(2, NUMBER_SAMPLES);

	SampleSphere2D(samples);

	flux.resize(height, width);
	flux.setZero();

	double px, py, dis_ij;
	uint32_t CI, CJ;
	Eigen::Vector2i closet_p;
	Eigen::Vector2d dis_p;

	for (uint32_t row = 0; row < height; ++row)
		for (uint32_t col = 0; col < width; ++col)
		{
			if ((row <= (uint32_t)epsilon) || (col <= (uint32_t)epsilon) || (row >= height - (uint32_t)epsilon) || (col >= width - (uint32_t)epsilon))
				continue;
			double flux_value = 0.00;
			for (uint32_t i = 0; i < NUMBER_SAMPLES; ++i)
			{

				//sample point in the sphere with radius
				px = row + samples(0, i) + 0.5;
				py = col + samples(1, i) + 0.5;

				//index of cell which sample point fall into
				CI = std::floor(px);
				CJ = std::floor(py);

				//closest point on the boundary to that cell
				//if (distance_map.cp[CI][CJ].cloestp.size() == 0)
				//	continue;
				closet_p = *distance_map.cp[CI][CJ].cloestp.begin();
				//TODO: check error
				closet_p = (Eigen::Vector2i(px, py)-closet_p);
				//closet_p = (Eigen::Vector2i(CI, CJ)-closet_p);
				dis_ij = std::sqrt(closet_p[0] * closet_p[0] + closet_p[1] * closet_p[1]);
				if (dis_ij != 0.00)
				{
					dis_p = 1.00/dis_ij*Eigen::Vector2d(closet_p(0), closet_p(1)) ;
				}
				else
				{
					dis_p = Eigen::Vector2d(0.00, 0.00);
				}
				flux_value = flux_value + dis_p.dot(samples.col(i));
			}

			flux(row, col) = flux_value / (double)NUMBER_SAMPLES;
		}
}

void SkeletonEx::BindCloestPoints()
{
}

void SkeletonEx::ExThickSkeleton(bool is_8_nei)
{
	//double threshold = -1.00 / M_PI * (std::asin((double)epsilon / (double)(fluxThreshold + epsilon)) + ((double)epsilon / (double)(fluxThreshold + epsilon)));
	double threshold = fluxThreshold;//-0.41;
	std::cout << "Accuracy threshold under discrete elements:" << threshold << "\n";
	bool have_next = true;
	std::vector<Eigen::Vector2i> p_resort;

	newm = grid_value;

	std::vector<std::pair<Eigen::Vector2i, double>> p_d_pairs(height * width);
	Eigen::MatrixXi sort_mat(height, width);
	uint32_t i_rc = 0;
	for (uint32_t row = 0; row < height; ++row)
		for (uint32_t col = 0; col < width; ++col)
		{
			p_d_pairs[i_rc] = std::make_pair(Eigen::Vector2i(row, col), distance_map.dis(row, col));
			++i_rc;
		}
	std::sort(p_d_pairs.begin(), p_d_pairs.end(), compare_pd);
	//按照同伦薄化步骤进行逐步删除
	p_resort.resize(height * width);
	for (uint32_t i = 0; i < height * width; ++i)
	{
		p_resort[i] = p_d_pairs[i].first;
		sort_mat((p_d_pairs[i].first)(0), (p_d_pairs[i].first)(1)) = i;
	}
	std::vector<bool> deleted(p_resort.size(), false);

	for (uint32_t i = 0; i < height * width; ++i)
	{
		if (newm(p_resort[i][0], p_resort[i][1]) < 0)
		{
			deleted[i] = true;
		}
	}

	while (have_next)
	{
		have_next = false;
		for (uint32_t i = 0; i < height * width; ++i)
		{

			if ((newm(p_resort[i][0], p_resort[i][1]) == 0) && (!deleted[i]))
			{
				if (flux(p_resort[i][0], p_resort[i][1]) > threshold)
				{
					//判断删除是否造成亏格变化
					//K = F-E+V 求△K即可
					int32_t F_n = -1, E_n = 0, V_n = 0;

					//△E
					OpenMesh::FaceHandle nowf(grid_index(p_resort[i][0], p_resort[i][1]));
					//std::cout << "bdy:" << mesh.is_boundary(nowf);
					for (auto nei_f : mesh.ff_range(nowf))
					{
						++E_n;
						if (newm(grid_xy[nei_f.idx()](0), grid_xy[nei_f.idx()](1)) >= 0)
						{
							--E_n;
						}
					}
					E_n = -E_n;
					//△V
					uint32_t up0 = 0;
					for (auto v : mesh.fv_range(nowf))
					{
						up0 = 0;
						for (auto nei_f : mesh.vf_range(v))
						{
							if (newm(grid_xy[nei_f.idx()](0), grid_xy[nei_f.idx()](1)) >= 0)
							{
								++up0;
							}
						}
						if (up0 == 1)
							--V_n;
					}

					//未发生亏格变化 进行删除操作和边界修改
					if (F_n - E_n + V_n == 0)
					{
						//std::cout << "Deleted index: " << i<<" \n";
						deleted[i] = true;
						newm(p_resort[i][0], p_resort[i][1]) = -1;

						//按照四邻域或八邻域更新边界
						if (is_8_nei)
						{
							for (auto v4 : mesh.fv_range(nowf))
							{
								for (auto fv : mesh.vf_range(v4))
								{
									if ((fv != nowf) && (newm(grid_xy[fv.idx()](0), grid_xy[fv.idx()](1)) >= 0))
									{
										newm(grid_xy[fv.idx()](0), grid_xy[fv.idx()](1)) = 0;
									}
								}
							}
						}
						else
						{
							for (auto nei_f4 : mesh.ff_range(nowf))
							{
								if (newm(grid_xy[nei_f4.idx()](0), grid_xy[nei_f4.idx()](1)) >= 0)
								{
									newm(grid_xy[nei_f4.idx()](0), grid_xy[nei_f4.idx()](1)) = 0;
								}
							}
						}

						have_next = true;
						break;
					}
				}
			}
		}
		if (!have_next)
		{
			break;
		}
	};

	del_elem.resize(height, width);
	for (uint32_t i = 0; i < height * width; ++i)
	{
		if (deleted[i])
		{
			del_elem(p_resort[i](0), p_resort[i](1)) = 0;
		}
		else
		{
			del_elem(p_resort[i](0), p_resort[i](1)) = 1;
		}
	}

	//DisplayField(FieldType::AOF_BDY);
	//DisplayField(FieldType::Ske_AOF);

	//get thin skeleton
	have_next = true;

	//find end point
	std::vector<uint32_t> end_points(p_resort.size(), 0);
	std::map<uint32_t, std::vector<uint32_t>> p_type3_nei;
	std::vector<OpenMesh::FaceHandle> neic8(8);
	uint32_t count_de = 0, end_all = 0;
	bool count_ie = false;
	for (uint32_t i = 0; i < height * width; ++i)
	{
		//uint32_t ac = p_resort[i](0),ab=p_resort[i](1);
		//if(ab==3) std::cout<<"\ncoutc:"<<ac<<","<<ab<<","<<newm(p_resort[i](0), p_resort[i](1));
		if (newm(p_resort[i](0), p_resort[i](1))>=0)
		{
			count_de = 0;
			OpenMesh::FaceHandle nowf(grid_index(p_resort[i][0], p_resort[i][1]));
			Eigen::Vector2i index_nowf = p_resort[i];
			//bool nei8_contain_end = false;
			for (int32_t rowa = -1; rowa < 2; ++rowa)
			{
				for (int32_t cola = -1; cola < 2; ++cola)
				{
					if ((rowa == 0) && (cola == 0))
						continue;
					if ((p_resort[i](0) + rowa < 0) || (p_resort[i](1) + cola < 0) || (p_resort[i](0) + rowa >= height) || (p_resort[i](1) + cola >= width))
						continue;
					//if (end_points[sort_mat(p_resort[i][0] + rowa, p_resort[i][1] + cola)]>0) nei8_contain_end = true;
					if (newm(p_resort[i](0) + rowa, p_resort[i](1) + cola) >= 0)
					{
						++count_de;
						neic8[count_de - 1] = OpenMesh::FaceHandle(grid_index(p_resort[i](0) + rowa, p_resort[i](1) + cola));
					}
				}
			}
			//if (nei8_contain_end) continue;
			
			if (count_de == 1)
			{
				//std::cout<<"\ncount1:"<<index_nowf;
				end_points[i] = 1;
				++end_all;
				//break;
			}
			else if (count_de == 2)
			{
				//std::cout<<"\ncount2:"<<index_nowf;
				//count_ie = false;
				for (auto ff : mesh.ff_range(neic8[0]))
				{
					if (ff == neic8[1])
					{
						//count_ie = true;
						end_points[i] = 2;
						++end_all;
						break;
					}
				}
			}
			else if (count_de == 3)
			{
				//std::cout<<"\ncount3:"<<index_nowf;
				count_ie = false;

				std::vector<uint32_t> temp3;
				for (uint32_t kk = 0; kk < 3; ++kk)
				{
					for (auto ff : mesh.ff_range(neic8[kk]))
					{
						if (ff == neic8[(kk + 1) % 3])
						{
							//count_ie = true;
							end_points[i] = 3;
							count_ie = true;
							++end_all;

							auto axy = grid_xy[neic8[(kk + 2) % 3].idx()];
							temp3.push_back(sort_mat(axy(0), axy(1)));

							break;
						}
						else if (ff == neic8[(kk + 2) % 3])
						{
							//count_ie = true;
							end_points[i] = 3;
							count_ie = true;
							++end_all;

							auto axy = grid_xy[neic8[(kk + 1) % 3].idx()];
							temp3.push_back(sort_mat(axy(0), axy(1)));

							break;
						}
					}

					//if (count_ie)break;
				}
				p_type3_nei[i] = temp3;
			}
		}
	}
	std::cout << "Initial end points: " << end_all << "\n";

	for (uint32_t i = 0; i < height * width; ++i)
	{
		if (!deleted[i])
		{
			del_elem(p_resort[i](0), p_resort[i](1)) = end_points[i];
		}
		else
		{
			del_elem(p_resort[i](0), p_resort[i](1)) = -1;
		}
		//else {
		//	del_elem(p_resort[i](0), p_resort[i](1)) = 0;
		//}
	}
	//DisplayField(FieldType::Ske_End);
	/*
	for (uint32_t i = 0; i < height*width; ++i) {
		if (end_points[i] == 3) {
			uint32_t add_i = 0;
			for (int32_t rowa = -1; rowa < 2; ++rowa) {
				for (int32_t cola = -1; cola < 2; ++cola) {
					if ((rowa == 0) && (cola == 0))continue;
					if ((p_resort[i][0] + rowa < 0) || (p_resort[i][1] + cola < 0) || (p_resort[i][0] + rowa >= height) || (p_resort[i][1] + cola >= width)) continue;
					if (end_points[sort_mat(p_resort[i][0] + rowa, p_resort[i][1] + cola)] == 3) {
						
						end_points[sort_mat(p_resort[i][0] + rowa, p_resort[i][1] + cola)] = 0;
						end_all--;
						++add_i;
					}
				}
			}
			if (add_i == 0) end_points[i] = 0;
		}
	}*/
	for (int32_t i = (height * width) - 1; i >= 0; --i)
	{
		//std::cout << i << " ";
		if (end_points[i] == 3)
		{
			if (p_type3_nei.find(i) != p_type3_nei.end())
			{
				bool finded3 = false;
				for (auto a : p_type3_nei[i])
				{
					if (end_points[a] == 3)
					{
						//check topology conectivity

						//K = F-E+V 求△K即可
						int32_t F_n = -1, E_n = 0, V_n = 0;

						//△E
						OpenMesh::FaceHandle nowf(grid_index(p_resort[a][0], p_resort[a][1]));

						for (auto nei_f : mesh.ff_range(nowf))
						{
							++E_n;
							if (newm(grid_xy[nei_f.idx()](0), grid_xy[nei_f.idx()](1)) >= 0)
							{
								--E_n;
							}
						}
						E_n = -E_n;
						//△V
						uint32_t up0 = 0;
						for (auto v : mesh.fv_range(nowf))
						{
							up0 = 0;
							for (auto nei_f : mesh.vf_range(v))
							{
								if (newm(grid_xy[nei_f.idx()](0), grid_xy[nei_f.idx()](1)) >= 0)
								{
									++up0;
								}
							}
							if (up0 == 1)
								--V_n;
						}

						//未发生亏格变化 进行删除操作和边界修改
						if (F_n - E_n + V_n != 0)
						{
							continue;
						}

						finded3 = true;
						end_points[a] = 0;
						p_type3_nei.erase(a);
						end_all--;
					}
				}
				if (!finded3)
				{
					end_points[i] = 0;
					end_all--;
				}
			}
		}
	}
	std::cout << "Final end points: " << end_all << "\n";

	for (uint32_t i = 0; i < height * width; ++i)
	{
		if (!deleted[i])
		{
			del_elem(p_resort[i](0), p_resort[i](1)) = end_points[i];
		}
		else
		{
			del_elem(p_resort[i](0), p_resort[i](1)) = -1;
		}
		//else {
		//	del_elem(p_resort[i](0), p_resort[i](1)) = 0;
		//}
	}
	//DisplayField(FieldType::Ske_End);

	while (have_next)
	{
		have_next = false;
		for (uint32_t i = 0; i < height * width; ++i)
		{

			if ((newm(p_resort[i][0], p_resort[i][1]) == 0) && (!deleted[i]))
			{
				if (end_points[i] == 0)
				{
					//判断删除是否造成亏格变化
					//K = F-E+V 求△K即可
					int32_t F_n = -1, E_n = 0, V_n = 0;

					//△E
					OpenMesh::FaceHandle nowf(grid_index(p_resort[i][0], p_resort[i][1]));
					//std::cout << "bdy:" << mesh.is_boundary(nowf);
					for (auto nei_f : mesh.ff_range(nowf))
					{
						++E_n;
						if (newm(grid_xy[nei_f.idx()](0), grid_xy[nei_f.idx()](1)) >= 0)
						{
							--E_n;
						}
					}
					E_n = -E_n;
					//△V
					uint32_t up0 = 0;
					for (auto v : mesh.fv_range(nowf))
					{
						up0 = 0;
						for (auto nei_f : mesh.vf_range(v))
						{
							if (newm(grid_xy[nei_f.idx()](0), grid_xy[nei_f.idx()](1)) >= 0)
							{
								++up0;
							}
						}
						if (up0 == 1)
							--V_n;
					}

					//未发生亏格变化 进行删除操作和边界修改
					if (F_n - E_n + V_n == 0)
					{
						//std::cout << "Deleted index: " << i<<" \n";
						deleted[i] = true;
						newm(p_resort[i][0], p_resort[i][1]) = -1;

						//按照四邻域或八邻域更新边界
						if (is_8_nei)
						{
							for (auto v4 : mesh.fv_range(nowf))
							{
								for (auto fv : mesh.vf_range(v4))
								{
									if ((fv != nowf) && (newm(grid_xy[fv.idx()](0), grid_xy[fv.idx()](1)) >= 0))
									{
										newm(grid_xy[fv.idx()](0), grid_xy[fv.idx()](1)) = 0;
									}
								}
							}
						}
						else
						{
							for (auto nei_f4 : mesh.ff_range(nowf))
							{
								if (newm(grid_xy[nei_f4.idx()](0), grid_xy[nei_f4.idx()](1)) >= 0)
								{
									newm(grid_xy[nei_f4.idx()](0), grid_xy[nei_f4.idx()](1)) = 0;
								}
							}
						}

						have_next = true;
						break;
					}
				}
			}
		}
		if (!have_next)
		{
			break;
		}
	};

	for (uint32_t i = 0; i < height * width; ++i)
	{
		if (deleted[i])
		{
			del_elem(p_resort[i](0), p_resort[i](1)) = 0;
		}
		else
		{
			del_elem(p_resort[i](0), p_resort[i](1)) = 1;
		}
	}
	//DisplayField(FieldType::Ske_AOF_D);
}

void SkeletonEx::ExThinSkeleton(double pis)
{
}

void SkeletonEx::ShiftSkeleton()
{
}

void SkeletonEx::BuildBranches()
{
	//collect all skeleton points
	for (uint32_t row = 0; row < height; ++row)
	{
		for (uint32_t col = 0; col < width; ++col)
		{
			if (del_elem(row, col) == 0)
				continue;
			uint32_t count = 0;
			std::vector<Eigen::Vector2i> neic;
			for (int32_t rowa = -1; rowa < 2; ++rowa)
			{
				for (int32_t cola = -1; cola < 2; ++cola)
				{
					if ((rowa == 0) && (cola == 0))
						continue;
					if ((row + rowa < 0) || (col + cola < 0) || (row + rowa >= height) || (col + cola >= width))
						continue;
					if (del_elem(row + rowa, col + cola) == 1)
					{
						++count;
						neic.push_back(Eigen::Vector2i(rowa, cola));
					}
				}
			}
			SkeletonPoint ske;

			if (count == 1)
			{ //end point
				ske.type = SkPType::End;
				//check the ske_point is boundary point or not
			}
			else if (count == 2)
			{ //regular point
				ske.type = SkPType::Regular;
				
			}
			else
			{ //junction point
				bool valid = true;
				for (uint32_t ii = 0; ii < neic.size(); ++ii)
				{
					for (uint32_t jj = 0; jj < neic.size(); ++jj)
					{
						if (ii == jj)
							continue;
						if (neic[ii](0) == neic[jj](0))
						{
							if (neic[ii](1) + neic[jj](1) != 0)
							{
								valid = false;
								goto HE_VALID;
							}
						}
						else if (neic[ii](1) == neic[jj](1))
						{
							if (neic[ii](0) + neic[jj](0) != 0)
							{
								valid = false;
								goto HE_VALID;
							}
						}
					}
				}
			HE_VALID:
				if (valid)
				{
					ske.type = SkPType::Junction;
				}
				else
				{
					ske.type = SkPType::Regular;
				}
			}

			ske.index = ske_points.size();
			ske.coor = Eigen::Vector2i(row, col);
			auto c = mesh.calc_face_centroid(OpenMesh::FaceHandle(grid_index(row, col)));
			ske.coord = Eigen::Vector2d(c[0], c[1]);
			ske_points.push_back(ske);
			ske_p_index[Eigen::Vector2i(row, col)] = ske.index;
		}
	}
	std::cout << "Number of skeleton points:" << ske_points.size() << "\n";
	//DisplayField(FieldType::Ske_BRANCH_P);

	

	std::vector<bool> is_checked(ske_points.size(), false);
	uint32_t idb = 0;
	for (uint32_t i = 0; i < ske_points.size(); ++i)
	{
		if (is_checked[i])
			continue;
		if (ske_points[i].type == SkPType::Regular)
			continue;

		//Iterated finding next point
		bool have_next = true;
		uint32_t row, col;
		uint32_t lasti = i, nowi = i;
		std::map<Eigen::Vector2i, uint32_t>::iterator auto_i;

		std::vector<uint32_t> branch_all;
		//find all branch link this point
		row = ske_points[nowi].coor[0];
		col = ske_points[nowi].coor[1];
		for (int32_t rowa = -1; rowa < 2; ++rowa)
		{
			for (int32_t cola = -1; cola < 2; ++cola)
			{
				if ((rowa == 0) && (cola == 0))
					continue;
				if ((row + rowa < 0) || (col + cola < 0) || (row + rowa >= height) || (col + cola >= width))
					continue;
				if ((auto_i = ske_p_index.find(Eigen::Vector2i(row + rowa, col + cola))) == ske_p_index.end())
					continue;
				if (is_checked[(*auto_i).second])
					continue;
				branch_all.push_back((*auto_i).second);
			}
		}

		is_checked[i] = true;
		if (branch_all.size() == 0)
			continue;

		for (uint32_t bi = 0; bi < branch_all.size(); ++bi)
		{
			lasti = i;
			nowi = branch_all[bi];
			SkeletonBranch branch;
			branch.id = idb++;
			branch.vlist.push_back(i);
			branch.vlist.push_back(nowi);

			do
			{
				row = ske_points[nowi].coor[0];
				col = ske_points[nowi].coor[1];
				have_next = false;
				is_checked[nowi] = true;

				//check include end or junction
				for (int32_t rowa = -1; rowa < 2; ++rowa)
				{
					for (int32_t cola = -1; cola < 2; ++cola)
					{
						if ((rowa == 0) && (cola == 0))
							continue;
						if ((row + rowa < 0) || (col + cola < 0) || (row + rowa >= height) || (col + cola >= width))
							continue;
						if ((auto_i = ske_p_index.find(Eigen::Vector2i(row + rowa, col + cola))) == ske_p_index.end())
							continue;
						if (is_checked[(*auto_i).second])
							continue;

						if (ske_points[(*auto_i).second].type == SkPType::End)
						{
							is_checked[(*auto_i).second] = true;
							branch.vlist.push_back((*auto_i).second);
							have_next = false;
							goto OUTTER_H;
						}
						else if (ske_points[(*auto_i).second].type == SkPType::Junction)
						{
							have_next = false;
							branch.vlist.push_back((*auto_i).second);
							goto OUTTER_H;
						}
					}
				}

				for (int32_t rowa = -1; rowa < 2; ++rowa)
				{
					for (int32_t cola = -1; cola < 2; ++cola)
					{
						if ((rowa == 0) && (cola == 0))
							continue;
						if ((row + rowa < 0) || (col + cola < 0) || (row + rowa >= height) || (col + cola >= width))
							continue;
						if ((auto_i = ske_p_index.find(Eigen::Vector2i(row + rowa, col + cola))) == ske_p_index.end())
							continue;
						if (is_checked[(*auto_i).second])
							continue;

						lasti = nowi;
						nowi = (*auto_i).second;

						branch.vlist.push_back(nowi);

						if (ske_points[nowi].type == SkPType::Regular)
						{
							is_checked[nowi] = true;
							have_next = true;
							goto OUTTER_H;
						}
					}
				}
			OUTTER_H:
			{
			}

			} while (have_next);
			branches.push_back(branch);
		}

		//all the branch link this junction point
		is_checked[i] = true;
	}

	//DisplayField(FieldType::Ske_BRANCH);
	std::cout << "Number of branches:" << branches.size() << "\n";

	//DisplayField(FieldType::Ske_SEGLINE);



	uint32_t from1,from2,th1,to1,to2;
	uint32_t kl,kr;
	Eigen::Vector2d h2,h1,h12,h3,h32,ave;

	//smooth skeleton line
	std::vector<Eigen::Vector2d> new_skep;
	new_skep.resize(ske_points.size()) ;
	uint32_t kkk = 0;
	while (kkk++<2)
	{
		/* code */
	
	
	for(uint32_t i = 0 ; i<branches.size() ; ++i){	
		//std::cout<<"``````branch:"<<i<<"\n";
		for(kl = 0; kl<branches[i].vlist.size() ; ++kl){
			if((kl == 0)||(kl == branches[i].vlist.size()-1)){
				new_skep[branches[i].vlist[kl]] = ske_points[branches[i].vlist[kl]].coord;
				continue;
			}else if((kl == 1)||(kl == branches[i].vlist.size()-2)){
				from1 = branches[i].vlist[kl-1] ;
				th1 = branches[i].vlist[kl] ;
				to1 = branches[i].vlist[kl+1];
				h1 = ske_points[from1].coord ;
				h2 = ske_points[th1].coord ;
				h3 = ske_points[to1].coord ;
				auto brans_temp = (1.00/2.00*(h3+h1));
				new_skep[branches[i].vlist[kl]] = brans_temp;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,brans_temp(0),brans_temp(1))) ;
				continue;
				
			}else{
				from2 = branches[i].vlist[kl-2] ;
				from1 = branches[i].vlist[kl-1] ;
				th1 = branches[i].vlist[kl] ;
				to1 = branches[i].vlist[kl+1];
				to2 = branches[i].vlist[kl+2];
				h12 = ske_points[from2].coord ;
				h1 = ske_points[from1].coord ;
				h2 = ske_points[th1].coord ;
				h3 = ske_points[to1].coord ;
				h32 = ske_points[to2].coord ;
				auto brans_temp = (1.00/4.00*(h3+h1+h12+h32));
				new_skep[branches[i].vlist[kl]] = brans_temp;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,brans_temp(0),brans_temp(1))) ;
				continue;
			}
		}
		
	}
	for(uint32_t i = 0; i<ske_points.size() ; ++i){
		ske_points[i].coord = new_skep[i];
	}

	}
//DisplayField(FieldType::Ske_SEGLINE);
	
	//simple normal
	for(uint32_t i = 0 ; i<branches.size() ; ++i){	
		//std::cout<<"``````branch:"<<i<<"\n";
		for(kl = 0; kl<branches[i].vlist.size() ; ++kl){
			if(kl == 0){
				th1 = branches[i].vlist[kl] ;
				to1 = branches[i].vlist[kl+1];
				h2 = ske_points[to1].coord ; 
				h1 = ske_points[th1].coord ;
				auto brans_temp = (h2-h1).normalized();
				ske_points[th1].normal_brans.push_back(std::make_pair(branches[i].id,brans_temp)) ;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,brans_temp(0),brans_temp(1))) ;
				continue;
			}else if((kl==1)||(kl == branches[i].vlist.size()-2)){

				th1 = branches[i].vlist[kl] ;
				to1 = branches[i].vlist[kl+1];
				h1 = ske_points[th1].coord ;
				h2 = ske_points[to1].coord ;
				auto brans_temp = (h2-h1).normalized();
				ske_points[th1].normal_brans.push_back(std::make_pair(branches[i].id,brans_temp)) ;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,brans_temp(0),brans_temp(1))) ;
				continue;
				
			}else if(kl == branches[i].vlist.size()-1){
				th1 = branches[i].vlist[kl] ;
				from1 = branches[i].vlist[kl-1];
				h2 = ske_points[th1].coord ; 
				h1 = ske_points[from1].coord ;
				auto brans_temp = (h2-h1).normalized();
				ske_points[th1].normal_brans.push_back(std::make_pair(branches[i].id,brans_temp)) ;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,brans_temp(0),brans_temp(1))) ;
				continue;
			}else{
				
				th1 = branches[i].vlist[kl];
				to1 = branches[i].vlist[kl+1];
				
				h1 = ske_points[th1].coord ;
				h2 = ske_points[to1].coord ;

				auto brans_temp = (h2-h1).normalized();
				ske_points[th1].normal_brans.push_back(std::make_pair(branches[i].id,brans_temp)) ;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,ave1(0),ave1(1))) ;
				//std::cout<<th1<<" normal:("<<ave1(0)<<","<<ave1(1)<<")\n";
				continue;
			}
		}
		
	}

	/*
	for(uint32_t i = 0 ; i<branches.size() ; ++i){	
		//std::cout<<"``````branch:"<<i<<"\n";
		for(kl = 0; kl<branches[i].vlist.size() ; ++kl){
			if(kl == 0){
				th1 = branches[i].vlist[kl] ;
				to1 = branches[i].vlist[kl+1];
				h2 = ske_points[to1].coord ; 
				h1 = ske_points[th1].coord ;
				auto brans_temp = (h2-h1).normalized();
				ske_points[th1].normal_brans.push_back(std::make_pair(branches[i].id,brans_temp)) ;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,brans_temp(0),brans_temp(1))) ;
				continue;
			}else if((kl==1)||(kl == branches[i].vlist.size()-2)){
				
				from1 = branches[i].vlist[kl-1] ;
				th1 = branches[i].vlist[kl] ;
				to1 = branches[i].vlist[kl+1];
				h2 = ske_points[th1].coord ; 
				h1 = ske_points[from1].coord ;
				h3 = ske_points[to1].coord ;
				auto brans_temp = (h3-h1).normalized();
				ske_points[th1].normal_brans.push_back(std::make_pair(branches[i].id,brans_temp)) ;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,brans_temp(0),brans_temp(1))) ;
				continue;
				
			}else if(kl == branches[i].vlist.size()-1){
				th1 = branches[i].vlist[kl] ;
				from1 = branches[i].vlist[kl-1];
				h2 = ske_points[th1].coord ; 
				h1 = ske_points[from1].coord ;
				auto brans_temp = (h2-h1).normalized();
				ske_points[th1].normal_brans.push_back(std::make_pair(branches[i].id,brans_temp)) ;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,brans_temp(0),brans_temp(1))) ;
				continue;
			}else{
				from2 = branches[i].vlist[kl-2] ;
				from1 = branches[i].vlist[kl-1] ;
				th1 = branches[i].vlist[kl];
				to1 = branches[i].vlist[kl+1];
				to2 = branches[i].vlist[kl+2] ;
				h12 = ske_points[from2].coord ; 
				h1 = ske_points[from1].coord ; 
				h2 = ske_points[th1].coord ;
				h3 = ske_points[to1].coord ;
				h32 = ske_points[to2].coord ; 


				auto index_th1 = ske_points[th1].coor;
				auto cloestp1 = *(distance_map.cp[index_th1(0)][index_th1(1)].cloestp.begin());
				auto cloestv1_ = mesh.calc_face_centroid(OpenMesh::FaceHandle(grid_index(cloestp1(0), cloestp1(1))));
				auto cloestv1 = Eigen::Vector2d(cloestv1_[0],cloestv1_[1]) ;
				auto p11 = (h32-h12).normalized();
				auto p22 = (h3-h1).normalized();
				Eigen::Vector2d ave1 = (1.00/3.00*p11)+(2.00/3.00*p22) ;
				ave1 = ave1.normalized();
				//120 TODO: end direction
				
				//std::cout<<ave<<"\n";
				ske_points[th1].normal_brans.push_back(std::make_pair(branches[i].id,ave1)) ;
				//ske_points[th1].rotation_ans.push_back(getRotateAngle(1.00,0.00,ave1(0),ave1(1))) ;
				//std::cout<<th1<<" normal:("<<ave1(0)<<","<<ave1(1)<<")\n";
				continue;
			}
		}
		
	}*/

	//calculate all rotation angle
	uint32_t k,k1 ;
	for(uint32_t i = 0; i<branches.size() ; ++i){
		for(kl = 0 ; kl<branches[i].vlist.size() ; ++kl){
			if(kl==0){
				k1 = branches[i].vlist[kl] ;
				uint32_t j = 0 ;
				for(j = 0 ; j<ske_points[k1].normal_brans.size() ; ++j){
					if(ske_points[k1].normal_brans[j].first == i) break ;
				}
				ske_points[k1].rotation_ans.push_back(get180RotateAngle(1.00,0.00,ske_points[k1].normal_brans[j].second(0),ske_points[k1].normal_brans[j].second(1))) ;

			}else{
				k = branches[i].vlist[kl-1] ;
				k1 = branches[i].vlist[kl] ;
				uint32_t j = 0,j1 = 0;
				for(j = 0 ; j<ske_points[k].normal_brans.size() ; ++j){
					if(ske_points[k].normal_brans[j].first == i) break ;
				}
				for(j1 = 0 ; j1<ske_points[k1].normal_brans.size() ; ++j1){
					if(ske_points[k1].normal_brans[j1].first == i) break ;
				}
				auto rot = ske_points[k].rotation_ans[j] ;
				auto rot1 = get180RotateAngle(ske_points[k].normal_brans[j].second(0),ske_points[k].normal_brans[j].second(1),ske_points[k1].normal_brans[j1].second(0),ske_points[k1].normal_brans[j1].second(1));
				if(rot1>M_PI){
					rot1 = rot1-2*M_PI ;
				}
				ske_points[k1].rotation_ans.push_back(rot+rot1) ;
			}
			
		}
	}
	
	//calculate all alphas
	//TODO: 需要修改α
	for(uint32_t i = 0 ; i<branches.size() ; ++i){	
		for(kl = 0; kl<branches[i].vlist.size() ; ++kl){
			auto i_ = branches[i].vlist[kl];
			//TODO::目前为尺寸相关，后续更改
			//Eigen::Vector2i index_th1 = ske_points[i_].coor;
			Eigen::Vector2i index_th1 ;
			index_th1(0) = (int)(std::abs(ske_points[i_].coord(1)));
			index_th1(1) = width -1 - (int)(std::abs(ske_points[i_].coord(0)));
			//重新定位位置
			auto cloestp1 = *(distance_map.cp[index_th1(0)][index_th1(1)].cloestp.begin());
			auto cloestv1_ = mesh.calc_face_centroid(OpenMesh::FaceHandle(grid_index(cloestp1(0), cloestp1(1))));
			auto cloestv1 = Eigen::Vector2d(cloestv1_[0],cloestv1_[1]) ;
			ske_points[i_].cloest_p = cloestv1 ;
			ske_points[i_].dis = (cloestv1 - ske_points[i_].coord).norm();
			auto link_ = (cloestv1-ske_points[i_].coord).normalized();

			for(uint32_t j = 0 ; j<ske_points[i_].normal_brans.size() ; ++j){
				if(ske_points[i_].normal_brans[j].first != i)continue;
				
				ske_points[i_].alphas.push_back(std::acos(link_.dot(ske_points[i_].normal_brans[j].second))) ;
				
			}
			
		}
	}

	//adjust alphas and rotation
	for(uint32_t i = 0 ; i<branches.size() ; ++i){	
		for(kl = 0; kl<branches[i].vlist.size() ; ++kl){
			uint32_t j = 0;
			auto k1 = branches[i].vlist[kl];
			
			for(j = 0 ; j<ske_points[k1].normal_brans.size() ; ++j){
				if(ske_points[k1].normal_brans[j].first != i)continue;
				break ;
			}
			if(kl==0){
				auto k2 = branches[i].vlist[kl+1] ;
				if(std::abs(ske_points[k1].alphas[j]-ske_points[k2].alphas[0])>0.4){
					ske_points[k1].alphas[j] = ske_points[k2].alphas[0];
				}
				if(std::abs(ske_points[k1].rotation_ans[j]-ske_points[k2].rotation_ans[0])>0.4){
					ske_points[k1].rotation_ans[j] = ske_points[k2].rotation_ans[0];
				}
			}else if(kl==branches[i].vlist.size()-1){
				auto k0 = branches[i].vlist[kl-1] ;
				if(std::abs(ske_points[k1].alphas[j]-ske_points[k0].alphas[0])>0.4){
					ske_points[k1].alphas[j] = ske_points[k0].alphas[0];
				}

				if(std::abs(ske_points[k1].rotation_ans[j]-ske_points[k0].rotation_ans[0])>0.4){
					ske_points[k1].rotation_ans[j] = ske_points[k0].rotation_ans[0];
				}
			}

			
		}
	}
	
	//DisplayField(FieldType::Ske_NOR_VECTOR);

}

void SkeletonEx::Display(uint32_t type_sk)
{
}
/*
void SkeletonEx::DisplayField(FieldType type_field)
{
#ifdef GUI_ON
	if (type_field == FieldType::Flag)
	{
		DataVisualization data;
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());
		auto propsBdy = OpenMesh::FProp<int>(mesh, bdy_flag.c_str());

		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				color[i * grid_index.cols() + j] = propsBdy[OpenMesh::FaceHandle(grid_index(i, j))];
			}
		}
		data.ShowGrayPic2D("BoundaryFlag", grid_index.rows(), grid_index.cols(), color);
	}
	else if (type_field == FieldType::Distance)
	{
		DataVisualization data;
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());

		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				color[i * grid_index.cols() + j] = distance_map.dis(i, j);
			}
		}
		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("DistanceField", grid_index.rows(), grid_index.cols(), color);
	}
	else if (type_field == FieldType::AOF)
	{
		DataVisualization data;
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());

		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				color[i * grid_index.cols() + j] = flux(i, j);
			}
		}
		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("AOFField", grid_index.rows(), grid_index.cols(), color);
	}
	else if (type_field == FieldType::AOF_BDY)
	{
		DataVisualization data;
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());

		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				color[i * grid_index.cols() + j] = newm(i, j);
			}
		}
		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("AOFBoundary", grid_index.rows(), grid_index.cols(), color);
	}
	else if (type_field == FieldType::Ske_End)
	{
		DataVisualization data;
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());

		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				color[i * grid_index.cols() + j] = del_elem(i, j);
			}
		}
		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("AOFEndPoints", grid_index.rows(), grid_index.cols(), color);
	}
	else if (type_field == FieldType::Ske_AOF)
	{
		DataVisualization data;
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());

		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				color[i * grid_index.cols() + j] = del_elem(i, j);
			}
		}
		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("AOFThickSkeleton", grid_index.rows(), grid_index.cols(), color);
	}
	else if (type_field == FieldType::Ske_AOF_D)
	{
		DataVisualization data;
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());

		auto propsBdy = OpenMesh::FProp<int>(mesh, bdy_flag.c_str());

		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				color[i * grid_index.cols() + j] = propsBdy[OpenMesh::FaceHandle(grid_index(i, j))];
			}
		}
		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				if(del_elem(i, j)!=0)
					color[i * grid_index.cols() + j] = 5;
			}
		}
		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("AOFThinSkeleton", grid_index.rows(), grid_index.cols(), color);
	}
	else if (type_field == FieldType::Ske_BRANCH_P)
	{
		DataVisualization data;
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());

		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				color[i * grid_index.cols() + j] = 0; // del_elem(i, j);
			}
		}

		for (uint32_t i = 0; i < ske_points.size(); ++i)
		{
			if (ske_points[i].type == SkPType::Regular)
			{
				color[ske_points[i].coor(0) * grid_index.cols() + ske_points[i].coor(1)] = 1;
			}
			else if (ske_points[i].type == SkPType::Junction)
			{
				color[ske_points[i].coor(0) * grid_index.cols() + ske_points[i].coor(1)] = 2;
			}
			else if (ske_points[i].type == SkPType::End)
			{
				color[ske_points[i].coor(0) * grid_index.cols() + ske_points[i].coor(1)] = 3;
			}
		}
		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("BranchPoints", grid_index.rows(), grid_index.cols(), color);
	}
	else if (type_field == FieldType::Ske_BRANCH)
	{
		DataVisualization data;
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());

		for (uint32_t i = 0; i < ske_points.size(); ++i)
		{
			if (ske_points[i].type == SkPType::Junction)
			{
				color[ske_points[i].coor(0) * grid_index.cols() + ske_points[i].coor(1)] = 1;
			}
			else if (ske_points[i].type == SkPType::End)
			{
				color[ske_points[i].coor(0) * grid_index.cols() + ske_points[i].coor(1)] = 2;
			}
		}

		for (uint32_t i = 0; i < branches.size(); ++i)
		{
			for (uint32_t j = 0; j < branches[i].vlist.size(); ++j)
			{
				auto t = ske_points[branches[i].vlist[j]];
				color[t.coor(0) * grid_index.cols() + t.coor(1)] = 3 * (i + 1);
			}
		}

		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("Branches", grid_index.rows(), grid_index.cols(), color);
	}
	else if (type_field == FieldType::Ske_SEGLINE)
	{
		DataVisualization data;
		std::vector<std::vector<float>> x, y;
		for (uint32_t i = 0; i < branches.size(); ++i)
		{
			std::vector<float> x_(branches[i].vlist.size()), y_(branches[i].vlist.size());
			for (uint32_t j = 0; j < branches[i].vlist.size(); ++j)
			{
				x_[j] = ske_points[branches[i].vlist[j]].coord[0];
				y_[j] = ske_points[branches[i].vlist[j]].coord[1];
			}
			x.push_back(x_);
			y.push_back(y_);
		}
		data.ShowSegLines("SegLines", x, y);
	}
	else if (type_field == FieldType::Ske_NOR_VECTOR)
	{
		DataVisualization data;
		
		std::vector<std::vector<float>> all;
		all.resize(4);
		
		for (uint32_t i = 0; i < ske_points.size(); ++i)
		{
			//std::cout<<"\ncccccccccc:"<<i<<"\n";
			
			if(ske_points[i].type == SkPType::Regular){
				all[0].push_back(ske_points[i].coord(0)) ;
				all[1].push_back(ske_points[i].coord(1)) ;
				//std::cout<<"re";
				all[2].push_back(ske_points[i].normal_brans[0].second(0)) ;
				all[3].push_back(ske_points[i].normal_brans[0].second(1)) ;
				//std::cout<<i<<" "<<ske_points[i].coord(0)<<" "<<ske_points[i].coord(1)<<" "<<ske_points[i].normal_brans[0].second[0]<<" "<<ske_points[i].normal_brans[0].second[1]<<"\n";
			}else if(ske_points[i].type == SkPType::Junction){
				for(uint32_t j = 0 ; j<ske_points[i].normal_brans.size() ; ++j){
					all[0].push_back(ske_points[i].coord(0)) ;
					all[1].push_back(ske_points[i].coord(1)) ;
					//std::cout<<"re";
					all[2].push_back(ske_points[i].normal_brans[j].second(0)) ;
					all[3].push_back(ske_points[i].normal_brans[j].second(1)) ;
				}
			}
		}
		data.ShowVectorField("Normal", all);
		
	
		std::vector<float> color;
		color.resize(grid_index.rows() * grid_index.cols());

		auto propsBdy = OpenMesh::FProp<int>(mesh, bdy_flag.c_str());

		for (uint32_t i = 0; i < grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < grid_index.cols(); ++j)
			{
				color[i * grid_index.cols() + j] = 0.1*propsBdy[OpenMesh::FaceHandle(grid_index(i, j))];
			}
		}
		for (uint32_t i = 0; i < ske_points.size(); ++i)
		{
			
			if(ske_points[i].type==SkPType::Regular)
				color[ske_points[i].coor(0) * grid_index.cols() + ske_points[i].coor(1)] = ske_points[i].alphas[0];
			
		}
		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("Alpha", grid_index.rows(), grid_index.cols(), color);
	}
#endif 
}*/
