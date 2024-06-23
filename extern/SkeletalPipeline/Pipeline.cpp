#include "Pipeline.h"

polyMesh& Pipeline::BuildSkeletonToFile(std::string meshPath, uint32_t width, uint32_t height, uint32_t samplesCir, double fluxThreshold, double epsilon){
	ske = new SkeletonEx(samplesCir,fluxThreshold,epsilon,width,height) ;
	ske->ReadTopMesh(meshPath,true);
	ske->InitBoundary(false);
	ske->ComDistanceField();
	ske->ComGradField();
	ske->ExThickSkeleton(false);
	ske->BuildBranches();
	auto propsBdy = OpenMesh::FProp<int>(ske->mesh, SKE_STR.c_str());

	for (uint32_t i = 0; i < ske->ske_points.size(); ++i)
	{
		Eigen::Vector2i index_th1l;
		index_th1l(0) = width - 1 - ske->ske_points[i].coor(1);
		index_th1l(1) = ske->ske_points[i].coor(0);
		//std::cout<<index_th1l(0) * ske->grid_index.cols() + index_th1l(1)<<" ";
		for (uint32_t i = 0; i < ske->ske_points.size(); ++i)
		{
			Eigen::Vector2i index_th1l;
			index_th1l(0) = width - 1 - ske->ske_points[i].coor(1);
			index_th1l(1) = ske->ske_points[i].coor(0);
			OpenMesh::FaceHandle h(ske->grid_index(index_th1l(1), index_th1l(0)));
			//std::cout<<index_th1l(0) * ske->grid_index.cols() + index_th1l(1)<<" ";
			if (ske->ske_points[i].type == SkPType::Regular)
			{
				propsBdy[h] = 1;
			} else if (ske->ske_points[i].type == SkPType::Junction)
			{
				propsBdy[h] = 2;
			} else if (ske->ske_points[i].type == SkPType::End)
			{
				propsBdy[h] = 3;
			}
		}
		//std::cout << "FaceHandle: " << index_th1l(0) * ske->grid_index.cols() + index_th1l(1) << std::endl;
	}
	
	return ske->mesh;
}
polyMesh& Pipeline::SparseOptimization(uint32_t iterationTimes, double dMaxError){
	pc = new PointClusters(*ske,dMaxError) ;
	pc->ExtractOrderedVList() ;
	//ske.DisplayField(FieldType::Ske_AOF);
	pc->init_sparse_optimizer();
	//p.set_stop_it(std::pow(10,-4),std::pow(10,-4),std::pow(10,-7)) ;
	for(uint32_t i =0  ;i <iterationTimes ; ++i){
		pc->set_initial_val();
	}
	return pc->se.mesh;
}
void Pipeline::GenerateSymAxis(std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>>& axis) {
	ca = new ComponentSymmetryAxi(pc->se,pc->pos_split,pc->branchClusteres) ;
	ca->CalOBBSymAxis();
	ca->UpdateHalfLens();
	for (uint32_t i = 0; i < ca->sym_axi_arr.size(); ++i)
	{
		for(uint32_t j = 0 ; j<ca->sym_axi_arr[i].size() ; ++j){
		
			Eigen::Vector2d p1 = ca->sym_axi_arr[i][j].center+ca->sym_axi_arr[i][j].half_len*ca->sym_axi_arr[i][j].ori;
			Eigen::Vector2d p2 = ca->sym_axi_arr[i][j].center-ca->sym_axi_arr[i][j].half_len*ca->sym_axi_arr[i][j].ori;
			axis.push_back(std::make_pair(p1,p2)) ;
		}
		
	}
	
}
polyMesh& Pipeline::ShapeApproximation(uint32_t maxIterations) {
	sat = new ShapeApproximate(pc->se,pc->pos_split,ca->sym_axi_arr,ca->sym_axi_ends,"./para.txt","./mmcmesh",1) ;
	sat->InitMMCList();
	sat->init_sparse_optimizer();
	//sat.mmc_list.display_all() ;
	sat->set_initial_val(maxIterations) ;
	return ske->mesh;
}
void Pipeline::ExtractMMCBoundary(std::vector<std::vector<Eigen::Vector2d>>& bdy,double sampleLen){
	sat->mmc_list.DisplayBoundary(sampleLen,bdy);
}