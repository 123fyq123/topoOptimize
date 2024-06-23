#include "PointClustersCal.h"
//TODO:更新
void PointClusters::ExtractOrderedVList()
{
	//se
	double alpha;
	
	Eigen::Vector2d ori;
	for(uint32_t bi = 0 ; bi<se.branches.size() ; ++bi){
		std::vector<Eigen::Vector2d> temp_vi;
		for(uint32_t vi = 0 ; vi<se.branches[bi].vlist.size() ; ++vi){
			auto& nowp = se.ske_points[se.branches[bi].vlist[vi]];
			//if(nowp.type != SkPType::Regular)continue;
			//calculate 
			for(uint32_t ni = 0; ni<nowp.normal_brans.size() ; ++ni){
				auto branindex = nowp.normal_brans[ni].first;
				if(branindex != bi)continue;

				ori = nowp.normal_brans[ni].second;
				
				//std::vector<std::pair<uint32_t,std::vector<Eigen::Vector2d>>> branchClusteres;
				//temp_vi.push_back((nowp.dis*std::cos(nowp.alphas[ni]))*(nowp.normal_brans[ni].second)+nowp.coord);
				temp_vi.push_back(nowp.dis*std::cos(nowp.alphas[ni])*(nowp.normal_brans[ni].second)+nowp.coord) ;
				//temp_vi.push_back(((nowp.cloest_p-nowp.coord).dot(ori))*(nowp.normal_brans[ni].second)+nowp.coord);
			}


		}
		branchClusteres.push_back(std::make_pair(bi,temp_vi)) ;  
	}

	//加入双边对称点重排
    std::vector<Eigen::Vector2d> mainOris;
    for(uint32_t i = 0 ; i<se.branches.size(); ++i){
        Eigen::Vector2d ori (0,0) ;
        for(uint32_t j = 1 ; j<se.branches[i].vlist.size()-1; ++j){
            ori+=se.ske_points[se.branches[i].vlist[j]].normal_brans[0].second ;
        }
        ori = (1.00/(se.branches[i].vlist.size()-2.00))*ori;
        mainOris.push_back(ori) ;
    }

    for(uint32_t i = 0 ; i<se.branches.size(); ++i){
        auto a = Comp_v2d(mainOris[i]);
        std::sort(branchClusteres[i].second.begin(), branchClusteres[i].second.end(), a);
        
    }
    
	//DisplayField( FieldType::Com_SYM_P_VEC);
	//滤波器
	uint32_t times = 0;
	while(times++<3)
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		std::vector<Eigen::Vector2d> smooth = branchClusteres[i].second;
		for(uint32_t j = 1; j<branchClusteres[i].second.size()-1 ; ++j){
			smooth[j] = (smooth[j-1]+smooth[j+1]+smooth[j])/3;
		}
		branchClusteres[i].second = smooth;
	}


	//DisplayField( FieldType::Com_SYM_P_VEC);

}
/*
void PointClusters::DisplayField(FieldType type_field){
	#ifdef GUI_OPT_ON 
	if(type_field == FieldType::Sym_POINTS){
		DataVisualization data;
		std::vector<std::vector<float>> xy(2);
		
		for (uint32_t i = 0; i < branchClusteres.size(); ++i)
		{
			
			for (uint32_t j = 0; j < branchClusteres[i].second.size(); ++j)
			{
				xy[0].push_back(branchClusteres[i].second[j](0)) ;
				xy[1].push_back(branchClusteres[i].second[j](1));
			}
		}
		data.Show2DPlot("Plot2D", xy);
	}else if(type_field == FieldType::Sym_OPT){
		DataVisualization data;
		
		std::vector<std::vector<float>> all;
		all.resize(4);
		
		for (uint32_t i = 0; i < se.ske_points.size(); ++i)
		{
			//std::cout<<"\ncccccccccc:"<<i<<"\n";
			
			if(se.ske_points[i].type == SkPType::Regular){

				auto dc = Eigen::Vector2d(std::cos(copy_ves[i].rotation_ans[0]) ,std::sin(copy_ves[i].rotation_ans[0])) ;

				all[0].push_back(se.ske_points[i].coord(0)) ;
				all[1].push_back(se.ske_points[i].coord(1)) ;
				//std::cout<<"re";
				all[2].push_back(dc(0)) ;
				all[3].push_back(dc(1)) ;
				//std::cout<<i<<" "<<ske_points[i].coord(0)<<" "<<ske_points[i].coord(1)<<" "<<ske_points[i].normal_brans[0].second[0]<<" "<<ske_points[i].normal_brans[0].second[1]<<"\n";
			}else if(se.ske_points[i].type == SkPType::Junction){
				for(uint32_t j = 0 ; j<se.ske_points[i].normal_brans.size() ; ++j){
					all[0].push_back(se.ske_points[i].coord(0)) ;
					all[1].push_back(se.ske_points[i].coord(1)) ;
					//std::cout<<"re";
					all[2].push_back(copy_ves[i].normal_brans[j].second(0)) ;
					all[3].push_back(copy_ves[i].normal_brans[j].second(1)) ;
				}
			}
		}
		data.ShowVectorField("OPTNormal", all);
	}else if(type_field == FieldType::Sym_BROLINES){
		DataVisualization data;
		
		for(uint32_t i = 0 ; i<result_show.size(); ++i){
			if(result_show[i]>2)result_show[i] = 0.00;
		}
		data.ShowBrokenLine("OPTNormal",result_show) ;
	}else if(type_field==FieldType::sym_SEGLINE_CO){
		DataVisualization data;
		std::vector<float> color;
		color.resize(se.grid_index.rows() * se.grid_index.cols());

		auto propsBdy = OpenMesh::FProp<int>(se.mesh, se.bdy_flag.c_str());

		for (uint32_t i = 0; i < se.grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < se.grid_index.cols(); ++j)
			{
				color[i * se.grid_index.cols() + j] = propsBdy[OpenMesh::FaceHandle(se.grid_index(i, j))];
			}
		}
		//3 5
		
		for(uint32_t i = 0 ; i<se.branches.size(); ++i){
			for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
				auto k = se.branches[i].vlist[j];
				if(se.ske_points[k].type != SkPType::Regular){
					color[se.ske_points[k].coor(0) * se.grid_index.cols() + se.ske_points[k].coor(1)] = 10;
					continue;
				} 
				//se.ske_points[k].normal_brans[0].second(0) = array_init[2*ci+0];
				//se.ske_points[k].normal_brans[0].second(1) = array_init[2*ci+1];
				
				if(j==se.branches[i].vlist.size()-2){
					
					color[se.ske_points[k].coor(0) * se.grid_index.cols() + se.ske_points[k].coor(1)] = 3;
					continue;
				}
				auto k_c = se.branches[i].vlist[j+1];
				auto re_a = (copy_ves[k].normal_brans[0].second).normalized();
				auto re_c = (copy_ves[k_c].normal_brans[0].second).normalized();
				auto rs = (re_a(0)-re_c(0))*(re_a(0)-re_c(0))+(re_a(1)-re_c(1))*(re_a(1)-re_c(1)) ;

				auto errc = std::abs(copy_ves[k].rotation_ans[0] - copy_ves[k_c].rotation_ans[0]);
				if(errc> del_pp){
					color[se.ske_points[k].coor(0) * se.grid_index.cols() + se.ske_points[k].coor(1)] = 10;
				}else{
					color[se.ske_points[k].coor(0) * se.grid_index.cols() + se.ske_points[k].coor(1)] = 3;
				}
				
			}
		}
		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("InitialSegment", se.grid_index.rows(), se.grid_index.cols(), color);
	}else if(type_field == FieldType::Com_SYM_PS){
		DataVisualization data;
		std::vector<float> color;
		color.resize(se.grid_index.rows() * se.grid_index.cols());

		auto propsBdy = OpenMesh::FProp<int>(se.mesh, se.bdy_flag.c_str());

		for (uint32_t i = 0; i < se.grid_index.rows(); ++i)
		{
			for (uint32_t j = 0; j < se.grid_index.cols(); ++j)
			{
				color[i * se.grid_index.cols() + j] = propsBdy[OpenMesh::FaceHandle(se.grid_index(i, j))];
			}
		}
		//3 5
		
		for(uint32_t i = 0 ; i<se.branches.size(); ++i){
			for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
				auto k = se.branches[i].vlist[j];
				if(se.ske_points[k].type != SkPType::Regular){
					color[se.ske_points[k].coor(0) * se.grid_index.cols() + se.ske_points[k].coor(1)] = 10;
					continue;
				} 
				//se.ske_points[k].normal_brans[0].second(0) = array_init[2*ci+0];
				//se.ske_points[k].normal_brans[0].second(1) = array_init[2*ci+1];
				
				
				color[se.ske_points[k].coor(0) * se.grid_index.cols() + se.ske_points[k].coor(1)] = 3;
				
				
			}
		}

		for(uint32_t i = 0 ; i<se.branches.size(); ++i){
			for(uint32_t j = 0 ; j<pos_split[i].size() ; ++j){
				auto k = se.branches[i].vlist[pos_split[i][j]];
				
				//se.ske_points[k].normal_brans[0].second(0) = array_init[2*ci+0];
				//se.ske_points[k].normal_brans[0].second(1) = array_init[2*ci+1];
				
				
				color[se.ske_points[k].coor(0) * se.grid_index.cols() + se.ske_points[k].coor(1)] = 10;
				
				
			}
		}

		//std::cout << distance_map.dis;
		data.ShowGrayPic2D("LastSegment", se.grid_index.rows(), se.grid_index.cols(), color);
	}else if(type_field == FieldType::Com_SYM_P_VEC){
		DataVisualization data;
		std::vector<std::vector<float>> x, y;
		for (uint32_t i = 0; i < branchClusteres.size(); ++i)
		{
			std::vector<float> x_(branchClusteres[i].second.size()), y_(branchClusteres[i].second.size());
            for(uint32_t j = 0 ; j<branchClusteres[i].second.size() ; ++j){
                x_[j] = branchClusteres[i].second[j](0) ;
                y_[j] = branchClusteres[i].second[j](1) ;   
            }
			x.push_back(x_);
			y.push_back(y_);
		}
		data.ShowSegLines("SegLines", x, y);
		
	}
	#endif
}*/
void PointClusters::CalSparseCluster()
{
}

void PointClusters::GetLeastSquareLine()
{
}

void PointClusters::CalSymmetryAxis()
{
}

inline double PointClusters::GetSymmetryError(uint32_t from, uint32_t to)
{
	return 0.0;
}

void PointClusters::init_sparse_optimizer(){
	uint32_t count_regular = 0;
	uint32_t count_wi = 0;
	/*
	for(uint32_t i = 0 ; i<se.ske_points.size() ; ++i){
		if(se.ske_points[i].type == SkPType::Regular){
			++count_regular ;
		}
	}*/
	delta_error.resize(se.branches.size()) ;
	for(uint32_t i = 0 ; i<se.branches.size() ; ++i){
		delta_error[i] = 0;
		if(se.branches[i].vlist.size() != 3){
			
			for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
				if(se.ske_points[se.branches[i].vlist[j]].type == SkPType::Regular){
					++count_regular;
				}
				if((se.ske_points[se.branches[i].vlist[j]].type == SkPType::Regular)&&
				(j!=se.branches[i].vlist.size()-2)){
					++count_wi;
				}
				if(j>=1&&j<=se.branches[i].vlist.size()-3){
					Eigen::Vector2d nb1 = se.ske_points[se.branches[i].vlist[j+1]].normal_brans[0].second;//.normalized();
					Eigen::Vector2d nb2 = se.ske_points[se.branches[i].vlist[j]].normal_brans[0].second;//.normalized();
					double err = nb1.dot(nb2)/(nb1.norm()*nb2.norm());
					if(err>=1) err = 1;
					delta_error[i]+=std::abs(std::acos(err));
					//std::cout<<"d:"<<std::abs(std::acos(err))<<" ";
				}
			}
		}
		std::cout<<"i:"<<i<<" del:"<<delta_error[i]<<"\n";
	}
	
	
	pos_split.resize(se.branches.size()) ;
	val_num = count_regular;
	result_show.resize(val_num-1) ;
	weight_w.resize(count_wi,1.00) ;
	copy_ves = se.ske_points ;

	std::ofstream oof("./weight.txt") ;

	uint32_t wi = 0 ;
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		if(se.branches[i].vlist.size()==3)continue;
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
			auto k = se.branches[i].vlist[j];
			if(se.ske_points[k].type != SkPType::Regular) continue;

			if((se.branches[i].vlist.size()!=3)&&(j!=se.branches[i].vlist.size()-2)){
				auto k1 = se.branches[i].vlist[j+1];
				double temp_sum1 = 0;//std::abs(se.ske_points[k].rotation_ans[0]-se.ske_points[k1].rotation_ans[0]  )/des_c;
				//if(temp_sum1 < 1e-14) temp_sum1 = 1e-14;
				double error_fe = delta_error[i] ;//1/std::abs(se.ske_points[se.branches[i].vlist[1]].rotation_ans[0]-se.ske_points[se.branches[i].vlist[se.branches[i].vlist.size()-2]].rotation_ans[0]);
				//if(error_fe>1e14)error_fe = 1e14;
				weight_w[wi] = error_fe*std::exp(-temp_sum1*temp_sum1*temp_sum1*temp_sum1) ;
				//oof<<"\nbrani:"<<i<<" vlistj:"<<j<<" wi:"<<wi<<" temp_sum2 = "<<error_fe<<" temp_sum1 = "<<temp_sum1<<" weight:"<<weight_w[wi];
				wi++;
			}

		}
	}
	oof.close();
	
	//std::cout<<weight_w.size();
	//weight_w.assign(weight_w.size() , 1.00) ;
	//NLOPT_LN_PRAXIS//
	//local_opter = nlopt_create(NLOPT_LD_LBFGS,val_num*DIM_PRO) ;
	//nlopt_set_ftol_abs(local_opter,0.000001) ;//0.0000001
	global_opter = nlopt_create(NLOPT_LD_LBFGS,val_num) ;
	//nlopt_set_local_optimizer(global_opter,local_opter) ; 
	
}

double WarpperFunction( unsigned n, const double *x, double *grad, void *pdata ){
    
    return reinterpret_cast<PointClusters*>(pdata)->appro_term_fun(n,x,grad) ;
}
double WrapperConstraints(unsigned n, const double *x, double *grad, void *pdata){
    return reinterpret_cast<PointClusters*>(pdata)->boudary_constraints(n,x,grad) ;
}
void WrapperMConstraints(unsigned m, double* result, unsigned n, const double* x, double* grad, void* pdata){
	return reinterpret_cast<PointClusters*>(pdata)->multi_boudary_constraints(m,result,n,x,grad) ;
}

void PointClusters::multi_boudary_constraints(unsigned m, double *result, unsigned n, const double* x, double* grad){
	/*
	the index of grad:
	i*n+j 
	i: the index of constraint
	j: the index of value
	n: dimension of x
	*/
	uint32_t ci = 0;
	//calculate grad
	if(grad){
		//dis
		ci = 0;
		for (size_t i = 0; i < val_num; i++)
		for(size_t j = 0; j < 2*val_num ; j++)
		{
			grad[i*n+j] = 0.00;
		}
		
		for(uint32_t i = 0 ; i<se.branches.size(); ++i){
			if(se.branches[i].vlist.size()==3) continue;
			for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
				auto k = se.branches[i].vlist[j];
				if(se.ske_points[k].type != SkPType::Regular) continue;
				
				grad[2*val_num*ci+2*ci] = 2*x[2*ci];
				grad[2*val_num*ci+2*ci+1] = 2*x[2*ci+1];
				
				++ci;
				
			}
		}


	}
	
	//calculate value
	ci = 0;
	double new_value = 0.00;
	double all = 0.00,all_c = 0.00;
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		if(se.branches[i].vlist.size()==3) continue;
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
			auto k = se.branches[i].vlist[j];
			if(se.ske_points[k].type != SkPType::Regular) continue;
			
			new_value=(x[2*ci]-se.ske_points[k].normal_brans[0].second(0))*(x[2*ci]-se.ske_points[k].normal_brans[0].second(0))+
			(x[2*ci+1]-se.ske_points[k].normal_brans[0].second(1))*(x[2*ci+1]-se.ske_points[k].normal_brans[0].second(1));
			
			result[ci] = new_value-0.2;
			all+=(result[ci]>0.00)?1.00:0.00;
			all_c+=new_value-0.3;
			++ci;
			
		}
	}
	std::cout<<" cons:"<<all<<" val:"<<all_c<<"\n";
	return ;
}

//approximation error term + sparse term
double PointClusters::appro_term_fun(unsigned n, const double *x, double *grad){
	//std::cout<<"Iteration times: "<<iteration_times<<"\n";
	/*
	if(iteration_times>0){
		for(uint32_t i = 0 ; i<se.ske_points.size() ; ++i){
			if(se.ske_points[i].type != SkPType::Regular) continue;
			se.ske_points[i].normal_brans[0].second(0) = array_init[2*i+0];
			se.ske_points[i].normal_brans[0].second(1) = array_init[2*i+1];
		}
	}*/
	double alpha = 0.5;//50.00;
	uint32_t ci = 0 ;
	uint32_t wi = 0 ;

	//std::ofstream oof("./weightCal.txt") ;

	if(grad){
		//re-set grad value
		
		double grad_temp1 = 0.00,grad_temp2;
		for(uint32_t i = 0 ; i<se.branches.size(); ++i){
			if(se.branches[i].vlist.size()==3) continue;
			for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
				auto k = se.branches[i].vlist[j];
				if(se.ske_points[k].type != SkPType::Regular) continue;
				
				//auto k_1 = se.branches[i].vlist[j+1];
				grad_temp1 = 0.00;
				grad_temp2 = 0.00;
				//TODO: 暂时移除
				//grad_temp2+=2*(x[ci]-se.ske_points[k].rotation_ans[0]);

				//start point + 1
				if(j==1){
					double tempvc;
					if((x[ci+1]-x[ci])>0.00){
						tempvc = -alpha*weight_w[wi];
						grad_temp1+=-alpha*weight_w[wi];
					}else if((x[ci+1]-x[ci])<0.00){
						tempvc = -alpha*weight_w[wi];
						grad_temp1+=alpha*weight_w[wi];
					}

					//oof<<"\nbrani:"<<i<<" vlistj:"<<j<<" ci:"<<ci<<" wi:"<<wi<<" g1:"<<tempvc<<" x[ci]:"<<x[ci]<<" x[ci+1]:"<<x[ci+1];
					wi+=1;
					goto OUTTER_C;
				//end point - 1
				}else if(j==se.branches[i].vlist.size()-2){
					double tempvc;
					if((x[ci]-x[ci-1])>0.00){
						tempvc = alpha*weight_w[wi-1];
						grad_temp1+=alpha*weight_w[wi-1];
					}else if((x[ci]-x[ci-1])<0.00){
						tempvc = -alpha*weight_w[wi-1];
						grad_temp1+=-alpha*weight_w[wi-1];
					}
					//oof<<"\nbrani:"<<i<<" vlistj:"<<j<<" ci:"<<ci<<" wi:"<<wi<<" g1:"<<tempvc<<" x[ci]:"<<x[ci]<<" x[ci-1]:"<<x[ci-1];
					goto OUTTER_C;
				
				}else{
					double tempvc = 0.00;
					if((x[ci]-x[ci-1])>0.00){
						tempvc += alpha*weight_w[wi-1];
						grad_temp1+=alpha*weight_w[wi-1];
					}else if((x[ci]-x[ci-1])<0.00){
						tempvc += -alpha*weight_w[wi-1];
						grad_temp1+=-alpha*weight_w[wi-1];
					}

					if((x[ci+1]-x[ci])>0.00){
						tempvc+=-alpha*weight_w[wi];
						grad_temp1+=-alpha*weight_w[wi];
					}else if((x[ci+1]-x[ci])<0.00){
						tempvc+=alpha*weight_w[wi];
						grad_temp1+=alpha*weight_w[wi];
					}
					//oof<<"\nbrani:"<<i<<" vlistj:"<<j<<" ci:"<<ci<<" wi:"<<wi<<" g1:"<<tempvc<<" x[ci]:"<<x[ci]<<" x[ci-1]:"<<x[ci-1]<<" x[ci+1]:"<<x[ci+1];
					//std::cout<<"bran:"<<i<<" 3ci:"<<ci+3<<" wi:"<<wi+1<<" g1:"<<grad_temp1<<" g2:"<<grad_temp2<<"\n";
					wi+=1;
				}
				OUTTER_C:
				//std::cout<<"wi:"<<weight_w[wi-1]<<"g1:"<<grad_temp1<<" g2:"<<grad_temp2<<"\n";
				grad[ci] = grad_temp1+grad_temp2;
				
				
				ci+=1;
				
					
			}
		//std::cout<<"-------------bran:"<<i<<" ci:"<<ci<<" wi:"<<wi<<"\n";
		}
	}
	//oof.close();

	//re-calculate
	ci = 0;
	wi = 0;
	double new_value = 0.00;
	double value1 = 0.00,value2 = 0.00;
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		if(se.branches[i].vlist.size()==3) continue;
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
			auto k = se.branches[i].vlist[j];
			if(se.ske_points[k].type != SkPType::Regular) continue;
			
			//auto k_1 = se.branches[i].vlist[j+1];
			//TODO: 暂时移除
			//new_value+=(x[ci]-se.ske_points[k].rotation_ans[0])*(x[ci]-se.ske_points[k].rotation_ans[0]) ;

			//value1+=(x[ci]-se.ske_points[k].rotation_ans[0])*(x[ci]-se.ske_points[k].rotation_ans[0]) ;

			if(j!=se.branches[i].vlist.size()-2){
				new_value+=alpha*(weight_w[wi]*std::abs(x[ci]-x[ci+1]));
				value2+=alpha*(weight_w[wi]*std::abs(x[ci]-x[ci+1]));
				//std::cout<<"ci:"<<ci+3<<" wi:"<<wi+1<<" v:"<<new_value<<"\n";
				wi+=1;
			}

			ci+=1;
			
		}
	}

	if(iteration_times%100 == 0)
	std::cout<<"Iteration times: "<<iteration_times<<" value:"<<new_value<<" value1:"<<value1<<" value2:"<<value2<<"\n";
	++iteration_times;
	return new_value;
}
//constraints 
double PointClusters::boudary_constraints(unsigned n, const double *x , double *grad){
	//std::cout<<"Iteration times: "<<iteration_times<<"\n";
	uint32_t ci = 0;
	//if(!grad.empty()){
		for(uint32_t i = 0 ; i<se.branches.size(); ++i){
			for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
				auto k = se.branches[i].vlist[j];
				if(se.ske_points[k].type != SkPType::Regular) continue;
				
				//auto k_1 = se.branches[i].vlist[j+1];

				grad[ci] = 4.00*(x[ci]*x[ci]+x[ci+1]*x[ci+1]-1.00)*x[ci];
				grad[ci+1] = 4.00*(x[ci]*x[ci]+x[ci+1]*x[ci+1]-1.00)*x[ci+1];
				ci+=2;
		
			}
		}
	//}

	ci = 0;
	double new_value = 0.00;
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
			auto k = se.branches[i].vlist[j];
			if(se.ske_points[k].type != SkPType::Regular) continue;
			
			//auto k_1 = se.branches[i].vlist[j+1];

			new_value+=std::pow(x[ci]*x[ci]+x[ci+1]*x[ci+1]-1.00,2);

			ci+=2;
		}
	}
	std::cout<<"Iteration times: "<<iteration_times<<" cons:"<<new_value-0.000000001<<"\n";
	return new_value-0.000000001;
}
// lower bounds and upper bounds of (x_ori,y_ori)
void PointClusters::set_lower_upper_bounds(){
	double* array_upper = new double[val_num] ;
	double* array_lower = new double[val_num] ;
	uint32_t ci = 0;
	
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		if(se.branches[i].vlist.size()==3) continue;
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
			auto k = se.branches[i].vlist[j];
			if(se.ske_points[k].type != SkPType::Regular) continue;
			
			if((j==1)||(j==se.branches[i].vlist.size()-2)/*||(j==2)||(j==se.branches[i].vlist.size()-3)*/){
				array_lower[ci] = se.ske_points[k].rotation_ans[0]-0.1;
				array_upper[ci] = se.ske_points[k].rotation_ans[0]+0.1;
			}else{
				//array_lower[ci] = se.ske_points[k].rotation_ans[0]-0.6;
				//array_upper[ci] = se.ske_points[k].rotation_ans[0]+0.6;
				double dsita = 2*std::asin( 0.5*delta_d/se.ske_points[k].dis );
				array_lower[ci] = se.ske_points[k].rotation_ans[0]-dsita;
				array_upper[ci] = se.ske_points[k].rotation_ans[0]+dsita;
			}
			
			ci+=1;
			
		}
	}

	nlopt_set_lower_bounds(global_opter,array_lower) ;
    nlopt_set_upper_bounds(global_opter,array_upper) ;
}
//set stop threshold
void PointClusters::set_stop_it(double xtol, double ftol, double force){
	
	
	nlopt_set_ftol_abs(global_opter,force) ;

}
void PointClusters::set_initial_val() {
	/*
	uint32_t count_regular = 0;
	for(uint32_t i = 0 ; i<se.ske_points.size() ; ++i){
		if(se.ske_points[i].type == SkPType::Regular){
			++count_regular ;
		}
	}*/
	double* array_init = new double[val_num] ;
	
	uint32_t ci = 0;
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		if(se.branches[i].vlist.size()==3)continue;
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
			auto k = se.branches[i].vlist[j];
			if(se.ske_points[k].type != SkPType::Regular) continue;
			array_init[ci] = se.ske_points[k].rotation_ans[0];
			ci++;
		}
	}
	ci = 0;
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		if(se.branches[i].vlist.size()==3)continue;
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
			auto k = se.branches[i].vlist[j];
			if(se.ske_points[k].type != SkPType::Regular) continue;
			
			if(ci<val_num-1){
				auto rs = (array_init[ci+1]-array_init[ci])*(array_init[ci+1]-array_init[ci]);
				result_show[ci] = rs;
			}
			ci++;
		}
	}
	


	//DisplayField(FieldType::Sym_BROLINES);
	//DisplayField(FieldType::sym_SEGLINE_CO);
	double f_max = INF;
	
	nlopt_set_min_objective(global_opter,WarpperFunction,this) ;
	nlopt_set_maxeval(global_opter, 4000);
	set_lower_upper_bounds() ;
	//nlopt_add_inequality_mconstraint(global_opter,val_num,WrapperMConstraints,this,NULL) ;

	//nlopt_add_inequality_constraint(global_opter,WrapperConstraints,this,0.000001);
	nlopt_result result = nlopt_optimize(global_opter, array_init, &f_max) ;
	std::cout<<"result code:"<<result;
	ci = 0;
	uint32_t wi = 0 ;
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		if(se.branches[i].vlist.size()==3)continue;
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
			auto k = se.branches[i].vlist[j];
			if(se.ske_points[k].type != SkPType::Regular) continue;
			
			copy_ves[k].rotation_ans[0] = array_init[ci];

			if((se.branches[i].vlist.size()!=3)&&(j!=se.branches[i].vlist.size()-2)){
				double temp_sum1 = std::abs(array_init[ci+1]-array_init[ci])/des_c;
				double error_fe = delta_error[i];//1/std::abs(se.ske_points[se.branches[i].vlist[1]].rotation_ans[0]-se.ske_points[se.branches[i].vlist[se.branches[i].vlist.size()-2]].rotation_ans[0]);
				//if(temp_sum1 < 1e-14) temp_sum1 = 1e-14;
				//if(error_fe>1e14)error_fe = 1e14;
				weight_w[wi] = error_fe*std::exp(-temp_sum1*temp_sum1*temp_sum1*temp_sum1) ;
				
				wi++;
			}
			if(ci<val_num-1){
				result_show[ci] = (array_init[ci+1]-array_init[ci])*(array_init[ci+1]-array_init[ci]);
				
			}
			ci++;

		}
	}
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		std::vector<uint32_t>().swap(pos_split[i]);
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size(); ++j){
			auto k = se.branches[i].vlist[j];
			if(copy_ves[k].type!=SkPType::Regular){
				pos_split[i].push_back(j) ;
				continue;
			}
			if((j == 1) || (j==se.branches[i].vlist.size()-2)) continue ;
			auto k_1 = se.branches[i].vlist[j+1];
			
			auto rs = std::abs(copy_ves[k_1].rotation_ans[0]-copy_ves[k].rotation_ans[0]);
			if(rs>del_pp){
				pos_split[i].push_back(j) ;
			}

		}
	}

	//update pos_split
	
	//DisplayField(FieldType::sym_SEGLINE_CO);
    //for(uint32_t kc = 0 ; kc<2 ; ++kc)
	//{
		/* code */
	
	
    for(uint32_t i = 0; i<se.branches.size(); ++i){
        if(se.branches[i].vlist.size()<=3)continue;
        std::vector<uint32_t> leftall ;
        leftall.push_back(0);
        for(uint32_t j = 1; j<pos_split[i].size()-1; ++j){
            if((pos_split[i][j]==1)||(pos_split[i][j]==se.branches[i].vlist.size()-2)||
            (pos_split[i][j]==2)||(pos_split[i][j]==se.branches[i].vlist.size()-3)/*||
			(pos_split[i][j]==3)||(pos_split[i][j]==se.branches[i].vlist.size()-4)*/){
                continue;
            }

			
            std::vector<int32_t> j_collect;
            j_collect.push_back(pos_split[i][j]) ;
            
            for(uint32_t k = j+1; k<=pos_split[i].size()-1; ++k){
				if(pos_split[i][k]>=se.branches[i].vlist.size()-3){
					j = pos_split[i].size()-1;
					//std::cout<<"aj:"<<j;
					break ;
				}
                if(pos_split[i][k]-j_collect[j_collect.size()-1]<=4){
                    j_collect.push_back(pos_split[i][k]) ;
					
					//std::cout<<"bj:"<<j;
                }else{
                    j = k-1;
					//std::cout<<"cj:"<<j;
                    break;
                }
            }
            uint32_t selected = pos_split[i][j];
            double max_er = 0.00;
            for(uint32_t c = 0 ; c<j_collect.size() ; ++c){
				/*
                auto vnow_0 = se.branches[i].vlist[j_collect[c]-1];
                auto vnow = se.branches[i].vlist[j_collect[c]];
                auto vnow_1 = se.branches[i].vlist[j_collect[c]+1];

                double er_temp = std::abs(se.ske_points[vnow].rotation_ans[0]-se.ske_points[vnow_0].rotation_ans[0])+
                std::abs(se.ske_points[vnow].rotation_ans[0]-se.ske_points[vnow_1].rotation_ans[0]);
                if(er_temp>max_er){
                    selected = j_collect[c] ;
                    max_er = er_temp ;
                }*/
				if(std::abs((j_collect[c]-se.branches[i].vlist.size())/2.00)>max_er){
                    selected = j_collect[c] ;
                    max_er = std::abs((j_collect[c]-se.branches[i].vlist.size())/2.00) ;
                }
            }
            leftall.push_back(selected) ;

        }

        leftall.push_back(se.branches[i].vlist.size()-1);
        pos_split[i] = leftall;
    }
	//}
	//DisplayField(FieldType::Com_SYM_PS);
	//DisplayField(FieldType::Sym_OPT);
	//DisplayField(FieldType::Sym_BROLINES);
	
	auto propsBdy = OpenMesh::FProp<int>(se.mesh, se.bdy_flag.c_str());
	auto propsSP = OpenMesh::FProp<int>(se.mesh, std::string("ske_sp_flag").c_str());
	
	for (uint32_t i = 0; i < se.grid_index.rows(); ++i)
	{
		for (uint32_t j = 0; j < se.grid_index.cols(); ++j)
		{
			propsSP[OpenMesh::FaceHandle(se.grid_index(i, j))] = propsBdy[OpenMesh::FaceHandle(se.grid_index(i, j))];
		}
	}
	//3 5
	
	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		for(uint32_t j = 0 ; j<se.branches[i].vlist.size() ; ++j){
			auto k = se.branches[i].vlist[j];
			Eigen::Vector2i index_th1l ;
			index_th1l(0) = se.width -1- se.ske_points[k].coor(1);
			index_th1l(1) = se.ske_points[k].coor(0);
			OpenMesh::FaceHandle h(se.grid_index(index_th1l(1),index_th1l(0))) ;

			if(se.ske_points[k].type != SkPType::Regular){
				propsSP[h] = 10;
				continue;
			} 
			//se.ske_points[k].normal_brans[0].second(0) = array_init[2*ci+0];
			//se.ske_points[k].normal_brans[0].second(1) = array_init[2*ci+1];
			
			
			propsSP[h] = 3;
			
			
		}
	}

	for(uint32_t i = 0 ; i<se.branches.size(); ++i){
		for(uint32_t j = 0 ; j<pos_split[i].size() ; ++j){
			auto k = se.branches[i].vlist[pos_split[i][j]];
			Eigen::Vector2i index_th1l ;
			index_th1l(0) = se.width -1- se.ske_points[k].coor(1);
			index_th1l(1) = se.ske_points[k].coor(0);
			OpenMesh::FaceHandle h(se.grid_index(index_th1l(1),index_th1l(0))) ;
			//se.ske_points[k].normal_brans[0].second(0) = array_init[2*ci+0];
			//se.ske_points[k].normal_brans[0].second(1) = array_init[2*ci+1];
			
			
			propsSP[h] = 10;
			
			
		}
	}
	return ;
}