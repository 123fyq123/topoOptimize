#include "ShapeApproximate.h"



double WarpperFunctionMMC( unsigned n, const double *x, double *grad, void *pdata ){
    
    return reinterpret_cast<ShapeApproximate*>(pdata)->appro_term_fun(n,x,grad) ;
}
double WrapperConstraintsMMC(unsigned n, const double *x, double *grad, void *pdata){
    return reinterpret_cast<ShapeApproximate*>(pdata)->boudary_constraints(n,x,grad) ;
}
void WrapperMConstraintsMMC(unsigned m, double* result, unsigned n, const double* x, double* grad, void* pdata){
	return reinterpret_cast<ShapeApproximate*>(pdata)->multi_boudary_constraints(m,result,n,x,grad) ;
}

void ShapeApproximate::InitMMCList(){
    //init MMC list
    
    uint32_t count_mmc = 0;
    for(uint32_t i = 0 ; i<axis.size() ; ++i){
        for(uint32_t j = 0 ; j<axis[i].size() ; ++j){
            if(axis[i][j].half_len>0.00) count_mmc++;
        }
        //count_mmc+=axis[i].size() ;
    }
    std::cout<<"\ncout_of_mmc:"<<count_mmc<<"\n";


    mmc_list.num_of_mmc = count_mmc ;
    mmc_list.elem_list.resize(count_mmc) ;

   
    uint32_t count_i = 0;
    for(uint32_t i = 0 ; i<axis.size() ; ++i){
        for(uint32_t j = 0 ; j<axis[i].size() ; ++j){
            if(axis[i][j].half_len>0.00){ 


	        //(uint32_t i, uint32_t i_bran, uint32_t i_part, std::vector<double> ti, Eigen::Vector2d vo, double hl, Eigen::Vector2d an);
            std::vector<double> ti ;
            if(axis[i][j].inner){
                ti.resize(4) ;
            }else{
                ti.resize(3) ;
            }
            auto brani = axis[i][j].index_bran;
            auto indexi = axis[i][j].index_part;
            //auto starti = pos_split[i][indexi];
            //auto endi = pos_split[i][indexi+1];
            //auto dis_se = endi - starti+1;
            //参数位置
            double para_ti = 0.00;
            Eigen::Vector2d para_xy;
            for(uint32_t k = 0; k<ti.size() ; ++k){
                //TODO: 改进延伸出对称轴上末端宽度初始值
                para_ti = (double)k/(double)(ti.size()-1);
                para_xy = (1-para_ti)*axis[i][j].start_coord+para_ti*axis[i][j].end_coord ;
                
                //左端最近网格单元
                Eigen::Vector2i index_th1l ;
                index_th1l(0) = (int)( std::abs(para_xy(1)) );
                index_th1l(1) = se.width -1 - (int)( std::abs( para_xy(0) ));
                //std::cout<<"\n"<<para_xy<<"\n"<<index_th1l<<"\n";
                if(index_th1l(0)<0.00){
                    index_th1l(0) = 0 ;
                }else if(index_th1l(0)>(double)(se.height-1)){
                    index_th1l(0) = se.height-1;
                }
                if(index_th1l(1)<0.00){
                    index_th1l(1) = 0;
                }else if(index_th1l(1)>(double)(se.width-1)){
                    index_th1l(1) = se.width-1;
                }
                
                //std::cout<<"\n"<<para_xy<<"\n"<<index_th1l<<"\n";
                //重新定位位置
                auto cloestp1l = *(se.distance_map.cp[index_th1l(0)][index_th1l(1)].cloestp.begin());
                auto cloestv1l_ = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(cloestp1l(0), cloestp1l(1))));
                auto cloestv1l = Eigen::Vector2d(cloestv1l_[0],cloestv1l_[1]) ;

                //std::cout<<"\n i:"<<i<<" j:"<<j<<" ti[k]:"<<(cloestv1l-para_xy).norm();
                ti[k] = (cloestv1l-para_xy).norm();//std::abs(now_skp.dis*std::sin(now_skp.alphas[brani_]));

                
            }
            
            /*
            std::cout<<"\nbran:"<<brani<<" indexj:"<<indexi<<" ti:";
            for(uint32_t cc = 0 ; cc<4 ; ++cc){
                std::cout<<ti[cc]<<" ";
            }*/
            /*
            std::vector<uint32_t> left_mmc,right_mmc ;

            if(axis[i].size()!=1){

                if(j==0){
                    
                    right_mmc.push_back(count_i+1) ;
                }else if(j==axis[i].size()-1){
                    left_mmc.push_back(count_i-1) ;
                }else{
                    left_mmc.push_back(count_i-1) ;
                    right_mmc.push_back(count_i+1) ;
                }

            }*/
            std::cout<<"\nindex:"<<axis[i][j].index;
            mmc_list.initialize_one_component(axis[i][j].index,axis[i][j].index_bran,axis[i][j].index_part,ti,axis[i][j].center,axis[i][j].half_len,axis[i][j].ori) ;
            mmc_list.elem_list[axis[i][j].index].nei_mmcs_l = axis[i][j].l_neis ;
            mmc_list.elem_list[axis[i][j].index].nei_mmcs_r = axis[i][j].r_neis ;
            mmc_list.elem_list[axis[i][j].index].is_inner = (int)axis[i][j].inner ;
            mmc_list.elem_list[axis[i][j].index].cons_nei = axis[i][j].cons_neis;
            ++count_i ;
            }
        }

    }
    //约束连接
    for(uint32_t i = 0; i<mmc_list.elem_list.size(); ++i){
        for(uint32_t j = 0; j<mmc_list.elem_list[i].nei_mmcs_l.size(); ++j){
            auto nei = mmc_list.elem_list[i].nei_mmcs_l[j];
            if(mmc_list.elem_list[nei.first].index_branch != mmc_list.elem_list[i].index_branch)continue;
            if(nei.second==0){
                double temp = (mmc_list.elem_list[nei.first].ti[0]+mmc_list.elem_list[i].ti[0])/2;
                mmc_list.elem_list[nei.first].ti[0] = temp;
                mmc_list.elem_list[i].ti[0] = temp;
            }else{
                double temp = (mmc_list.elem_list[nei.first].ti.back()+mmc_list.elem_list[i].ti[0])/2;
                mmc_list.elem_list[nei.first].ti[mmc_list.elem_list[nei.first].ti.size()-1] = temp;
                mmc_list.elem_list[i].ti[0] = temp;
            }
        }
        for(uint32_t j = 0; j<mmc_list.elem_list[i].nei_mmcs_r.size(); ++j){
            auto nei = mmc_list.elem_list[i].nei_mmcs_r[j];
            if(mmc_list.elem_list[nei.first].index_branch != mmc_list.elem_list[i].index_branch)continue;
            if(nei.second==0){
                double temp = (mmc_list.elem_list[nei.first].ti[0]+mmc_list.elem_list[i].ti.back())/2;
                mmc_list.elem_list[nei.first].ti[0] = temp;
                mmc_list.elem_list[i].ti[mmc_list.elem_list[i].ti.size()-1] = temp;
            }else{
                double temp = (mmc_list.elem_list[nei.first].ti.back()+mmc_list.elem_list[i].ti.back())/2;
                mmc_list.elem_list[nei.first].ti[mmc_list.elem_list[nei.first].ti.size()-1] = temp;
                mmc_list.elem_list[i].ti[mmc_list.elem_list[i].ti.size()-1] = temp;
            }
        }
    }
    /*
    std::vector<std::pair<double,double>> l_r_update ;
    l_r_update.resize(mmc_list.elem_list.size()) ;
    //update half-lenght
    count_i = 0;
    for(uint32_t i = 0 ; i<axis.size() ; ++i){
        for(uint32_t j = 0 ; j<axis[i].size() ; ++j){
            if(axis[i][j].half_len>0.00){ 
                double left_add = 0.00,right_add = 0.00;
                if(j == 0){
                    auto p1 = ;
                    auto p2 = ;
                    auto p3 = ;
                }else if(j == axis[i].size()-1){
                    
                }
                axis[i][j].end_p2
            }
        }
    }*/
    

    //建立组件左右邻接组件关系


    /*
    //build left and right mmcs
    for(uint32_t i = 0 ; i<mmc_list.elem_list.size() ; ++i){
        
        //for junction components
        if(mmc_list.elem_list[i].index_part == 999){
            for(uint32_t j = 0 ; j < mmc_list.elem_list.size() ; ++j){
                if(mmc_list.elem_list[j].index_part != 999){
                    if(mmc_list.elem_list[i].nei_mmcs_l[0] == mmc_list.elem_list[j].index_branch){
                        mmc_list.elem_list[i].nei_mmcs_l[0] = j;
                    }
                    if(mmc_list.elem_list[i].nei_mmcs_r[0] == mmc_list.elem_list[j].index_branch){
                        mmc_list.elem_list[i].nei_mmcs_r[0] = j;
                    }

                    //build left and right mmcs between junction component and inner component
                    auto index_b = mmc_list.elem_list[j].index_branch;
                    auto index_r = mmc_list.elem_list[j].index_part;
                    
                    auto p2_lr = axis[index_b][index_r].end_p2;
                    auto p2_dis_1 = (se.ske_points[p2_lr[0]].coord - Eigen::Vector2d(mmc_list.elem_list[i].x0,mmc_list.elem_list[i].y0)).norm();
                    auto p2_dis_2 = (se.ske_points[p2_lr[1]].coord - Eigen::Vector2d(mmc_list.elem_list[i].x0,mmc_list.elem_list[i].y0)).norm();

                    if(p2_dis_1<p2_dis_2){
                        mmc_list.elem_list[j].nei_mmcs_l.push_back(i) ;
                    }else{
                        mmc_list.elem_list[j].nei_mmcs_r.push_back(i) ;
                    }

                }
            }
        //for inner components    
        }else{
            for(uint32_t j = 0 ; j < mmc_list.elem_list.size() ; ++j){
                if(mmc_list.elem_list[j].index_part != 999){
                    if(mmc_list.elem_list[j].index_branch == mmc_list.elem_list[i].index_branch){
                        if(mmc_list.elem_list[j].index_part == mmc_list.elem_list[i].index_part-1){
                            mmc_list.elem_list[i].nei_mmcs_l.push_back(j) ;
                        }else if(mmc_list.elem_list[j].index_part == mmc_list.elem_list[i].index_part+1){
                            mmc_list.elem_list[i].nei_mmcs_r.push_back(j) ;
                            
                        }
                    }
                }
            }
        }
    
        
    }*/

    //TODO:添加半长更新
    //update half length
    //std::vector<std::vector<double>> hl_update_l;
    //std::vector<std::vector<double>> hl_update_r;
    //hl_update_l.resize(mmc_list.elem_list.size()) ;
    //hl_update_r.resize(mmc_list.elem_list.size()) ;
    /*
    Eigen::MatrixXi is_checked;
    is_checked.resize(mmc_list.elem_list.size(),mmc_list.elem_list.size()) ;
    is_checked.setZero() ;
    for(uint32_t i = 0; i<hl_update.size() ; ++i){
        double left_add = 0.00,right_add = 0.00;
        for(uint32_t j = 0 ; j<mmc_list.elem_list[i].nei_mmcs_l.size() ; ++j){
            auto now_mmc_i = mmc_list.elem_list[i].nei_mmcs_l[j] ;
            if(is_checked(i,j)==1) continue ;

            is_checked(i,j) = 1;
            is_checked(j,i) = 1;
        }

        for(uint32_t j = 0 ; j<mmc_list.elem_list[i].nei_mmcs_r.size() ; ++j){
            auto now_mmc_i = mmc_list.elem_list[i].nei_mmcs_l[j] ;
            if(is_checked(i,j)==1) continue ;

            is_checked(i,j) = 1;
            is_checked(j,i) = 1;
        }
        
    }*/


}

void ShapeApproximate::searchNLayers(uint32_t N) {
    std::vector<int16_t> f_is_inner(se.mesh.n_faces(),0) ;
    auto propsBdy = OpenMesh::FProp<int>(se.mesh, se.bdy_flag.c_str());
    

    //find N layers of elements on boundary or inner regions
    for(uint32_t rowi = 0 ; rowi<se.distance_map.dis.rows() ; ++rowi){
        for(uint32_t coli = 0; coli<se.distance_map.dis.cols(); ++coli){
            //OpenMesh::FaceHandle nowf(grid_index(p_resort[i][0], p_resort[i][1]));
            //mesh.calc_face_centroid(OpenMesh::FaceHandle(grid_index(row, col)));
            if((se.distance_map.dis(rowi,coli)>0.00)&&(se.distance_map.dis(rowi,coli)<(double)(N+1))){
                auto p = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(rowi,coli)));    
                inner_p.push_back(Eigen::Vector2d(p[0],p[1]));
                f_is_inner[se.grid_index(rowi,coli)] = 1;
            }else if((se.distance_map.dis(rowi,coli)<0.00)&&(se.distance_map.dis(rowi,coli)>-(double)(N+1))){
                auto p = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(rowi,coli)));    
                outter_p.push_back(Eigen::Vector2d(p[0],p[1]));
                f_is_inner[se.grid_index(rowi,coli)] = -1;

            }
        }
    }
    

    
    for(uint32_t i = 0; i<se.ske_points.size() ; ++i){
        if(se.ske_points[i].type!=SkPType::Junction) continue;
        //find N neiboring cells
        
        auto center_now = se.ske_points[i].coor ;
        Eigen::Vector2i move_xy;
        int N_ = 5;
        //std::cout<<"f:"<<se.height<<" "<<se.width;
        for(int32_t x_i = -N_ ; x_i<=N_; ++x_i){
            for(int32_t y_i = -N_; y_i<=N_; ++y_i){
                move_xy(0) = center_now(0) + y_i;
                move_xy(1) = center_now(1) + x_i;
                
                if((move_xy(0)<0)||(move_xy(0)>(int)(se.height-1))||(move_xy(1)<0)||(move_xy(1)>(int)(se.width-1))) continue;
                auto f_index = OpenMesh::FaceHandle(se.grid_index(move_xy(0),move_xy(1))) ;
               
                if(propsBdy[f_index]<0) continue;
                
                if(f_is_inner[f_index.idx()]==0){
                    //std::cout<<"f:"<<f_index.idx() <<" "<<move_xy(0)<<" "<<move_xy(1);
                    f_is_inner[f_index.idx()] = 1;
                    auto p = se.mesh.calc_face_centroid(f_index);
                    inner_p.push_back(Eigen::Vector2d(p[0],p[1])) ;
                }

            }
            
        }
    }
    
    //init grid points
    //auto propsBdy = OpenMesh::FProp<int>(se.mesh,se.bdy_flag.c_str()) ;
    for(const auto& face : se.mesh.faces()){
        if(propsBdy[face] >0){
            /*
            auto cc = se.mesh.calc_face_centroid(face);
            inner_p.push_back(Eigen::Vector2d(cc[0],cc[1])) ;
            */
            //std::cout<<"inner-i:"<<face.idx()<<" d:"<<Eigen::Vector2d(cc[0],cc[1])<<"\n" ;
        }else if(propsBdy[face] ==0){
            auto cc = se.mesh.calc_face_centroid(face);
            //inner_p.push_back(Eigen::Vector2d(cc[0],cc[1])) ;
            for(auto v : se.mesh.fv_range(face)){
                for(auto f : se.mesh.vf_range(v)){
                    if(propsBdy[f]<0){
                        bdy_p.push_back(Eigen::Vector2d(se.mesh.point(v)[0],se.mesh.point(v)[1]));
                        break;
                    }
                }
                
            }
        }else{/*
            auto cc = se.mesh.calc_face_centroid(face);
            outter_p.push_back(Eigen::Vector2d(cc[0],cc[1])) ;*/
            //std::cout<<"outter-i:"<<face.idx()<<" d:"<<Eigen::Vector2d(cc[0],cc[1])<<"\n" ;
        }
    }

    //DataVisualization data;
    std::vector<float> color;
    color.resize(se.grid_index.rows() * se.grid_index.cols());

    //auto propsBdy = OpenMesh::FProp<int>(se.mesh, se.bdy_flag.c_str());

    for (uint32_t i = 0; i < se.grid_index.rows(); ++i)
    {
        for (uint32_t j = 0; j < se.grid_index.cols(); ++j)
        {
            color[i * se.grid_index.cols() + j] = propsBdy[OpenMesh::FaceHandle(se.grid_index(i, j))];
        }
    }
    for(uint32_t i = 0; i<f_is_inner.size() ; ++i){
        auto c = se.grid_xy[i];
        if(f_is_inner[i]==1){
            color[c(0) * se.grid_index.cols() + c(1)] = 4;
        }else if(f_is_inner[i]==-1){
            color[c(0) * se.grid_index.cols() + c(1)] = 6;
        }
    }


    //std::cout << distance_map.dis;
    //data.ShowGrayPic2D("flag", se.grid_index.rows(), se.grid_index.cols(), color);

}
//init optimizer
void ShapeApproximate::init_sparse_optimizer(){
    count_x = 0;
    for(uint32_t i = 0; i<mmc_list.elem_list.size() ; ++i){
        count_x+=mmc_list.elem_list[i].ti.size() ;
    }

    
    searchNLayers(2) ;

    global_opter = nlopt_create(NLOPT_LN_SBPLX,count_x) ;
}
//approximation error term + sparse term
double ShapeApproximate::appro_term_fun(unsigned n, const double *x, double *grad){
    //reset mmc_list
    uint32_t ci = 0;
    if(iteration_times>0){
        for(uint32_t i = 0; i<mmc_list.num_of_mmc; ++i){
            for(uint32_t j = 0; j<mmc_list.elem_list[i].ti.size(); ++j){
                mmc_list.elem_list[i].ti[j] = x[ci] ;
                ci++;
            }
        }
    }
    if(grad){

    }

    double result1 = 0.00;
    
    for(uint32_t i = 0 ; i<inner_p.size() ; ++i){
        double global_val = Heaviside( mmc_list.get_EM_xy(inner_p[i]) );
        //std::cout<<"1-i:"<<i<<"global_val:"<<global_val<<"\n";
        result1-=global_val;
    }

    
    double result2 = 0.00 ;
    /*
    for(uint32_t i = 0 ; i<bdy_p.size() ; ++i){
        //std::cout<<"\nsize:"<<bdy_p.size()<<"\n"<<bdy_p[i];
        double global_val = mmc_list.get_EM_xy(bdy_p[i]);
        //std::cout<<"1-i:"<<i<<"global_val:"<<global_val<<"\n";
        result2+=global_val*global_val;
    }*/
    

    double result3 = 0.00;
    
    for(uint32_t i = 0 ; i<outter_p.size() ; ++i){
        double global_val = Heaviside( mmc_list.get_EM_xy(outter_p[i]) );
       
        //std::cout<<"2-i:"<<i<<"global_val:"<<global_val<<"\n";
        result3+=global_val;
    }
    double ratioc = 1;//(double)inner_p.size()/(double)outter_p.size();
    if(iteration_times%100==0){
        std::cout<<"\niteration_times: "<<iteration_times<<" ratio:"<<ratioc<<" result1:"<<result1<<" result2:"<<result2;
        std::cout<<" result3:"<<ratioc*result3<<" result:"<<result1+result2+ratioc*result3 <<"\n";
    }
    
    if(iteration_times == max_iterations){//10000
        //std::string string_path = "./mmc_"+std::to_string((int)iteration_times)+".txt";
        //mmc_list.write_vtk(string_path,temp_m);
        //mmc_list.write_to_file(string_path);
        mmc_list.write_to_file(paraMMC);
        //mmc_list.display_all() ;
        
    }
    iteration_times++;
    return result1+result2+ratioc*result3;
}
//constraints 
double ShapeApproximate::boudary_constraints(unsigned n, const double *x , double *grad){
    return 0.0;
}

//multi-constraints
void ShapeApproximate::multi_boudary_constraints(unsigned m, double *result, unsigned n, const double* x, double* grad){

}

// lower bounds and upper bounds of (x_ori,y_ori)
void ShapeApproximate::set_lower_upper_bounds(std::pair<double,double> x_ori, std::pair<double,double> y_ori){

}
//set stop threshold
void ShapeApproximate::set_stop_it(double xtol, double ftol, double force){
    nlopt_set_ftol_abs(global_opter,force) ;
}

void ShapeApproximate::set_initial_val(uint32_t max_iterations){
    this->max_iterations = max_iterations;
    double* array_init = new double[count_x] ;
    uint32_t ci = 0;
    for(uint32_t i = 0 ; i<mmc_list.num_of_mmc; ++i){
        for(uint32_t j = 0 ; j<mmc_list.elem_list[i].ti.size() ; ++j){
            array_init[ci] = mmc_list.elem_list[i].ti[j];
            ++ci;
        }
    }

    double* array_lower = new double[count_x] ;
    double* array_upper = new double[count_x] ;
    ci = 0;
    double max_len = 0.00;
    for(uint32_t i = 0 ; i<mmc_list.num_of_mmc; ++i){
        max_len = 1.2*(*std::max_element(mmc_list.elem_list[i].ti.begin(),mmc_list.elem_list[i].ti.end())) ;
        for(uint32_t j = 0 ; j<mmc_list.elem_list[i].ti.size() ; ++j){
            array_lower[ci] = 0;
            array_upper[ci] = max_len;
            if(j==0||j==mmc_list.elem_list[i].ti.size()-1){
                array_lower[ci] = mmc_list.elem_list[i].ti[j]-0.01;
                array_upper[ci] = mmc_list.elem_list[i].ti[j]+0.01;
            }
            ++ci;
        }
    }
    
    for(uint32_t i = 0 ; i<count_x; ++i){
        std::cout<<array_init[i]<<" ";
    }
    double f_max = INF;
    nlopt_set_maxeval(global_opter,max_iterations+1) ;
    nlopt_set_min_objective(global_opter,WarpperFunctionMMC,this) ;
    nlopt_set_lower_bounds(global_opter,array_lower) ;
    nlopt_set_upper_bounds(global_opter,array_upper) ;
    std::cout<<"\n inner_p size:"<<inner_p.size()<<" outter_p size:"<<outter_p.size()<<"\n";
    
	nlopt_result result = nlopt_optimize(global_opter, array_init, &f_max) ;
	std::cout<<"result code:"<<result;

    auto propsSP = OpenMesh::VProp<double>(se.mesh, std::string("mmc_p_val").c_str());
	for(auto v : se.mesh.vertices()){
        auto vp = se.mesh.point(v);
        Eigen::Vector2d vc(vp[0],vp[1]) ;
        propsSP[v] = mmc_list.get_EM_xy(vc);
    }
	
    // free
 
    nlopt_destroy(global_opter);
    return;
}