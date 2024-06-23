#include "ComponentSymmetryAxi.h"

double get2RotateAngle(bool is_clockwise, double x1, double y1, double x2, double y2)
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
    if(is_clockwise)
	    return 2*M_PI-angle;
    else
        return angle;
}
//t1 t2同号代表为同向节点相接
bool Compute2LinesIntersection(Eigen::Vector2d& ise, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& rot1, const Eigen::Vector2d& rot2, double& t1, double& t2){
    
    if(rot1 == rot2) return false;
    if(rot1 == -rot2) return false;
    
    Eigen::Vector2d d1 = rot1.normalized();
    Eigen::Vector2d d2 = rot2.normalized();
    
    double offset1 = (d2(0)*(p1(1)-p2(1))-d2(1)*(p1(0)-p2(0)))/(d2(1)*d1(0)-d1(1)*d2(0));
    t1 = offset1;
    double offset2 = (d1(0)*(p2(1)-p1(1))-d1(1)*(p2(0)-p1(0)))/( d2(0)*d1(1)-d1(0)*d2(1) );
    //std::cout<<"\np1"<<p1(0)<<" "<<p1(1)<<" p2"<<p2(0)<<" "<<p2(1)<<" rot1"<<rot1(0)<<" "<<rot1(1)<<" rot2"<<rot2(0)<<" "<<rot2(1);
    t2 = offset2;
    ise = p1+offset1*d1;
    //std::cout<<" t1 "<<t1<<" t2 "<<t2;
    //std::cout<<" mid"<<ise(0)<<" "<<ise(1)<<"\n";

    
    return true ;
}
//为每个分支线计算主成分分析 i分支 i_now起始对称轴序号
void ComponentSymmetryAxi::OBBBranch(uint32_t i, uint32_t i_now){
        std::vector<Eigen::Vector2d> mid_p_vec;
        //从骨架分支线开始的计算
        if(se.ske_points[se.branches[i].vlist[0]].type == SkPType::Junction&&se.ske_points[se.branches[i].vlist.back()].type == SkPType::Junction
        &&se.branches[i].vlist.size()<15){
            mid_p_vec.resize(se.branches[i].vlist.size()) ;
            for(uint32_t j = 0; j<se.branches[i].vlist.size(); ++j){
                mid_p_vec[j] = se.ske_points[se.branches[i].vlist[j]].coord;
            }
        }else{//从双边对称线开始的计算
            mid_p_vec = branchClusteres[i].second ;
        }
        
        
        
        for(uint32_t j = 0 ; j<pos_split[i].size()-1; ++j){
            auto dis_12 = pos_split[i][j+1] - pos_split[i][j] + 1;
            Eigen::MatrixXd temp_mat;
            temp_mat.resize(2,dis_12) ;
            temp_mat.setZero() ;
            uint32_t st = 0 ; 
            Eigen::Vector2d v_ave(0.00,0.00);
            for(uint32_t k = pos_split[i][j]; k<=pos_split[i][j+1] ; ++k){
                auto cloest_p1 = se.ske_points[se.branches[i].vlist[k]].cloest_p ;
                Eigen::Vector2d norm_ori;
                uint32_t bi = 0;
                for(bi = 0 ; bi<se.ske_points[se.branches[i].vlist[k]].normal_brans.size(); ++bi){
                    if(se.ske_points[se.branches[i].vlist[k]].normal_brans[bi].first != i) continue;
                    norm_ori = se.ske_points[se.branches[i].vlist[k]].normal_brans[bi].second ;
                    break;
                }
                
                
                v_ave+=se.ske_points[se.branches[i].vlist[k]].normal_brans[bi].second ;
                
                temp_mat.col(st) = mid_p_vec[k] ;
                st+=1;
                 
            }
            v_ave = (1.00/dis_12)*v_ave ;
            v_ave.normalize() ;
            //std::cout<<"\ni_now:"<<i_now++<<"\n";
            //std::cout<<"\nid_i:"<<i<<" id_from:"<<pos_split[i][j]<<" id_to:"<<pos_split[i][j+1]<<"\n";
            //std::cout<<"v_ave:"<<v_ave(0)<<" "<<v_ave(1)<<"\n" ;
            //std::cout<<"\n"<<temp_mat<<"\n";
            
            auto cov = obb.Compute_covariance_matrix(temp_mat) ;
            obb.Compute_eigen_value(cov,temp_mat) ;
            
            SymAxi s;
            s.index_bran = i ;
            s.index_part = j ;
            s.center = obb.center ;
            s.end_p2(0) = pos_split[i][j];
            s.end_p2(1) = pos_split[i][j+1];
            s.end_p2_brani(0) = i;
            s.end_p2_brani(1) = i;
            s.end_p2_skpi(0) = se.branches[i].vlist[s.end_p2(0)];
            s.end_p2_skpi(1) = se.branches[i].vlist[s.end_p2(1)];
            s.index = i_now++;
            
            if(std::abs(obb.eigen_value(0))>std::abs(obb.eigen_value(1))){
                s.ori = obb.eigen_vectors.col(0); 
                if(v_ave.dot(obb.eigen_vectors.col(0))<-0.5) s.ori = -obb.eigen_vectors.col(0); 
                s.half_len = obb.half_axi_len[0];
            }else{
                s.ori = obb.eigen_vectors.col(1); 
                if(v_ave.dot(obb.eigen_vectors.col(1))<-0.5) s.ori = -obb.eigen_vectors.col(1); 
                s.half_len = obb.half_axi_len[1];
            }

            s.start_coord = s.center - s.half_len*s.ori;
            s.end_coord = s.center + s.half_len*s.ori;
            sym_axi_arr[i].push_back(s) ;
        }
    //}
}
bool ComponentSymmetryAxi::CheckBranchSymAxis(uint32_t i, bool only_e2e){
    //判断OBB出来的结果是否合理，不合理即合并掉稀疏优化的部分结果
    std::vector<uint32_t> resort_sp;
    uint32_t oSize = pos_split[i].size();
    resort_sp.push_back(pos_split[i][0]);
    double lens_sum = 0.00;
    for(uint32_t j = 0; j<sym_axi_arr[i].size(); ++j){
        lens_sum+=sym_axi_arr[i][j].half_len;
    }
    for(uint32_t j = 1; j<sym_axi_arr[i].size(); ++j){
        //相邻单位
        Eigen::Vector2d ise(0.00,0.00);
        double t1,t2;
        if(!Compute2LinesIntersection(ise,sym_axi_arr[i][j-1].center,sym_axi_arr[i][j].center,sym_axi_arr[i][j-1].ori,sym_axi_arr[i][j].ori,t1,t2)){
            continue;
        }
        //尾结点相连接
        if(t1*t2>0.00) continue;
        if(only_e2e){
            resort_sp.push_back(pos_split[i][j]) ;
        }else{
            if(std::abs(sym_axi_arr[i][j].ori.dot(sym_axi_arr[i][j-1].ori))<0.98){
                resort_sp.push_back(pos_split[i][j]) ;
            }else if(sym_axi_arr[i][j].half_len>lens_sum*0.2){
                resort_sp.push_back(pos_split[i][j]) ;
            }
        }
    }
    resort_sp.push_back(pos_split[i].back());
    pos_split[i] = resort_sp;
    if(oSize!=pos_split[i].size()) return true;
    else return false;
}
void ComponentSymmetryAxi::CalOBBSymAxis(){
    sym_axi_arr.resize(se.branches.size()) ;

    //当前序号
    uint32_t i_now = 0 ;
    
    for(uint32_t i = 0; i<se.branches.size(); ++i){
        //初次对双边对称线分段生成OBB
        OBBBranch(i,i_now) ;
        //检查是否存在近似平行对称轴，合并掉，重新OBB
        if(CheckBranchSymAxis(i,true)){
            //清理无效结果
            std::vector<SymAxi>().swap(sym_axi_arr[i]) ;
            //重生成
            OBBBranch(i,i_now) ;
        }
        while(CheckBranchSymAxis(i,true)){
            //清理无效结果
            std::vector<SymAxi>().swap(sym_axi_arr[i]) ;
            //重生成
            OBBBranch(i,i_now) ;
        }
        i_now += pos_split[i].size()-1; 
    }


    //调整连接点(先确定左右终点)
    //DisplayField(FieldType::Com_SYM_LINES) ;
    for(uint32_t i=0; i<se.branches.size() ; ++i){
        for(uint32_t j = 0 ; j<pos_split[i].size()-2; ++j){
            //相邻单位
            Eigen::Vector2d ise(0.00,0.00);
            double t1,t2;
            if(!Compute2LinesIntersection(ise,sym_axi_arr[i][j].center,sym_axi_arr[i][j+1].center,sym_axi_arr[i][j].ori,sym_axi_arr[i][j+1].ori,t1,t2)){
                std::cout<<"error: 平行对称轴\n";
                continue;
            }
            
            //尾结点相连接
            if(t1*t2>0.00){
                if((sym_axi_arr[i][j+1].half_len>sym_axi_arr[i][j].half_len)&&(std::abs(t1)>2*sym_axi_arr[i][j+1].half_len)&&(std::abs(t1)>5*sym_axi_arr[i][j].half_len)) continue;
                if((sym_axi_arr[i][j].half_len>sym_axi_arr[i][j+1].half_len)&&(std::abs(t2)>2*sym_axi_arr[i][j].half_len)&&(std::abs(t2)>5*sym_axi_arr[i][j+1].half_len)) continue;
                Eigen::Vector2d p11;
                Eigen::Vector2d p12;
                if(t1>0){
                    p11 = (sym_axi_arr[i][j].center-sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori);
                    p12 = (sym_axi_arr[i][j].center+t1*sym_axi_arr[i][j].ori);
                }else{
                    p11 = (sym_axi_arr[i][j].center+sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori);
                    p12 = (sym_axi_arr[i][j].center+t1*sym_axi_arr[i][j].ori);
                }
                

                //Eigen::Vector2d p12 = ise;
                sym_axi_arr[i][j].center = (p11+p12)/2.00;
                sym_axi_arr[i][j].half_len = (p11-p12).norm()/2.00;
                sym_axi_arr[i][j].start_coord = sym_axi_arr[i][j].center - sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori;
                sym_axi_arr[i][j].end_coord = sym_axi_arr[i][j].center + sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori;
                
                

                Eigen::Vector2d p21,p22;
                
                if(t2>0){
                    p21 = (sym_axi_arr[i][j+1].center-sym_axi_arr[i][j+1].half_len*sym_axi_arr[i][j+1].ori);
                    p22 = (sym_axi_arr[i][j+1].center+t2*sym_axi_arr[i][j+1].ori);
                }else{
                    p21 = (sym_axi_arr[i][j+1].center+sym_axi_arr[i][j+1].half_len*sym_axi_arr[i][j+1].ori);
                    p22 = (sym_axi_arr[i][j+1].center+t2*sym_axi_arr[i][j+1].ori);
                }
                
                //Eigen::Vector2d p22 = ise;
                sym_axi_arr[i][j+1].center = (p21+p22)/2.00;
                sym_axi_arr[i][j+1].half_len = (p21-p22).norm()/2.00;
                sym_axi_arr[i][j].r_update = true;
                sym_axi_arr[i][j+1].start_coord = sym_axi_arr[i][j+1].center - sym_axi_arr[i][j+1].half_len*sym_axi_arr[i][j+1].ori;
                sym_axi_arr[i][j+1].end_coord = sym_axi_arr[i][j+1].center + sym_axi_arr[i][j+1].half_len*sym_axi_arr[i][j+1].ori;
                
                
                if(t1>0){
                    sym_axi_arr[i][j].r_neis.push_back(std::make_pair( sym_axi_arr[i][j+1].index,1)) ;
                    sym_axi_arr[i][j+1].r_neis.push_back(std::make_pair( sym_axi_arr[i][j].index,1)) ;
                }else{
                    sym_axi_arr[i][j].l_neis.push_back(std::make_pair( sym_axi_arr[i][j+1].index,0)) ;
                    sym_axi_arr[i][j+1].l_neis.push_back(std::make_pair( sym_axi_arr[i][j].index,0)) ;
                }
                
                sym_axi_arr[i][j].inner = true;
                sym_axi_arr[i][j+1].inner = true;
                
            //首尾结点相连接
            }else{
                
                Eigen::Vector2d p11 = sym_axi_arr[i][j].center-sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori;
                Eigen::Vector2d p12 = sym_axi_arr[i][j].center+t1*sym_axi_arr[i][j].ori;;
                sym_axi_arr[i][j].center = (p11+p12)/2.00;
                sym_axi_arr[i][j].half_len = (p12-p11).norm()/2.00;

                sym_axi_arr[i][j].start_coord = sym_axi_arr[i][j].center - sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori;
                sym_axi_arr[i][j].end_coord = sym_axi_arr[i][j].center + sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori;

                Eigen::Vector2d p21 = sym_axi_arr[i][j+1].center+t2*sym_axi_arr[i][j+1].ori;;
                Eigen::Vector2d p22 = sym_axi_arr[i][j+1].center+sym_axi_arr[i][j+1].half_len*sym_axi_arr[i][j+1].ori;
                sym_axi_arr[i][j+1].center = (p21+p22)/2.00;
                sym_axi_arr[i][j+1].half_len = (p22-p21).norm()/2.00;

                sym_axi_arr[i][j+1].start_coord = sym_axi_arr[i][j+1].center - sym_axi_arr[i][j+1].half_len*sym_axi_arr[i][j+1].ori;
                sym_axi_arr[i][j+1].end_coord = sym_axi_arr[i][j+1].center + sym_axi_arr[i][j+1].half_len*sym_axi_arr[i][j+1].ori;
                

                sym_axi_arr[i][j].r_update = true;
                sym_axi_arr[i][j+1].l_update = true;
                
                sym_axi_arr[i][j].r_neis.push_back(std::make_pair( sym_axi_arr[i][j+1].index,0)) ;
                sym_axi_arr[i][j+1].l_neis.push_back(std::make_pair( sym_axi_arr[i][j].index,1)) ;
                sym_axi_arr[i][j].inner = true;
                sym_axi_arr[i][j+1].inner = true;
            }

        }
    }
    //骨架分支线首尾终结点以及对齐方向
    //TODO: 在分支密集的情况下可能存在很大隐患,
    //sym_axi_end_lr: 左右终点上与对称轴哪端对齐，左0右1
    //sym_axi_ends: 分支线首尾对称轴坐标
    //sym_axi_end_i.resize(se.branches.size()) ;
    sym_axi_end_lr.resize(se.branches.size()) ;
    sym_axi_ends.resize(se.branches.size()) ;
    for(uint32_t i = 0 ; i<se.branches.size() ; ++i){
        //std::cout<<" " ;
        auto p_si = se.branches[i].vlist[0] ;
        auto p_si1 = se.branches[i].vlist[1];
        auto junction_p_s = se.ske_points[p_si].coord ;
        auto junction_p_s1 = se.ske_points[p_si1].coord ;
        auto p1_start = sym_axi_arr[i][0].center+sym_axi_arr[i][0].half_len*sym_axi_arr[i][0].ori;
        auto p2_start = sym_axi_arr[i][0].center-sym_axi_arr[i][0].half_len*sym_axi_arr[i][0].ori;
        //std::cout<<"\n*******from:"<<p_si <<" "<<junction_p_s(0)<<" "<<junction_p_s(1);
        auto p_ei = *(se.branches[i].vlist.rbegin()) ;
        auto p_ei1 = *(se.branches[i].vlist.rbegin()+1);
        auto junction_p_e = se.ske_points[p_ei].coord ;
        auto junction_p_e1 = se.ske_points[p_ei1].coord ;
        auto p1_end = sym_axi_arr[i][sym_axi_arr[i].size()-1].center+sym_axi_arr[i][sym_axi_arr[i].size()-1].half_len*sym_axi_arr[i][sym_axi_arr[i].size()-1].ori;
        auto p2_end = sym_axi_arr[i][sym_axi_arr[i].size()-1].center-sym_axi_arr[i][sym_axi_arr[i].size()-1].half_len*sym_axi_arr[i][sym_axi_arr[i].size()-1].ori;
        //std::cout<<" to:"<<p_ei <<" "<<junction_p_e(0)<<" "<<junction_p_e(1);
        //std::cout<<"\n p1s:"<<p1_start(0)<<" "<<p1_start(1);
        //std::cout<<"\n p2s:"<<p2_start(0)<<" "<<p2_start(1);
        //std::cout<<"\n p1e:"<<p1_end(0)<<" "<<p1_end(1);
        //std::cout<<"\n p2e:"<<p2_end(0)<<" "<<p2_end(1);

        Eigen::Vector2d p1s;
        uint32_t l_link_flag ;
        if((junction_p_s-junction_p_s1).dot(p1_start-p2_start)>=0.00){
            p1s = p1_start ;
            l_link_flag = 1;
        }else{
            p1s = p2_start;
            l_link_flag = 0;
        }

       
        //std::cout<<"\n start:"<<p1s(0)<<" "<<p1s(1);

        Eigen::Vector2d p1e;
        uint32_t r_link_flag ;
        if((junction_p_e-junction_p_e1).dot(p1_end-p2_end)>=0.00){
            p1e = p1_end ;
            r_link_flag = 1;
        }else{
            p1e = p2_end;
            r_link_flag = 0;
        }

        sym_axi_end_lr[i] = std::make_pair(l_link_flag,r_link_flag) ;
        //std::cout<<"\n end:"<<p1e(0)<<" "<<p1e(1);
        sym_axi_ends[i] = std::make_pair(p1s,p1e) ;
    }

    //设置反向映射关系    
    for(uint32_t i = 0; i<sym_axi_arr.size() ; ++i){
        for(uint32_t j = 0; j<sym_axi_arr[i].size(); ++j){
            index_mapping.push_back(std::make_pair(i,j)) ;
        }
    }
    //DisplayField(FieldType::Com_SYM_LINES) ;
    //TODO: 有误
    //添加交叉点附近连接关系
    for(uint32_t i = 0; i<se.ske_points.size(); ++i){
        if(se.ske_points[i].type != SkPType::Junction)continue;
        std::cout<<"\njunction i:"<<i ;
        auto& ske_jun = se.ske_points[i];
        //bran_i-symAxi_i bool:left right
        std::unordered_map<uint32_t,std::pair<uint32_t,bool>> index_sym;
        //逆时针旋转角 rot-brani
        std::map<double,uint32_t> rot_bran;
        //收集交叉点附近组件
        for(uint32_t j = 0; j<se.ske_points[i].normal_brans.size(); ++j){
            auto bran_index = se.ske_points[i].normal_brans[j].first ;
            Eigen::Vector2d ori_temp ;
            //填充序号
            if(i==se.branches[bran_index].vlist[0]){//首连接
                
                index_sym[bran_index] = std::make_pair(sym_axi_arr[bran_index][0].index,
                sym_axi_end_lr[bran_index].first) ;
                ori_temp = (sym_axi_end_lr[bran_index].first? -sym_axi_arr[bran_index][0].ori: sym_axi_arr[bran_index][0].ori);
                //填充旋转角
                rot_bran[get2RotateAngle(false,1,0,ori_temp(0),ori_temp(1))] =  bran_index;
            }else{//尾连接
                
                index_sym[bran_index] = std::make_pair(sym_axi_arr[bran_index].back().index,
                sym_axi_end_lr[bran_index].second) ;
                ori_temp = (sym_axi_end_lr[bran_index].second? -sym_axi_arr[bran_index].back().ori: sym_axi_arr[bran_index].back().ori);
                //填充旋转角
                rot_bran[get2RotateAngle(false,1,0,ori_temp(0),ori_temp(1))] =  bran_index;
            }
            
            
        }
        std::vector<uint32_t> sym_resort;
        std::vector<std::pair<uint32_t,bool>> sym_all;
        //收集
        for(auto pairi : rot_bran){
            auto index_b = pairi.second;
            sym_resort.push_back(index_b) ;
            sym_all.push_back(index_sym[index_b]) ;
        }
        syms_jun_all[i] = sym_all;
        //加入symaxi的邻接关系
        for(uint32_t k = 0; k<sym_resort.size(); ++k){
            auto index_b = sym_resort[k] ;
            auto index_s = index_mapping[index_sym[index_b].first];

            //左右邻居
            auto index_l_b = sym_resort[(k+sym_resort.size()-1)%sym_resort.size()];
            auto index_l_s = index_mapping[index_sym[index_l_b].first];
            auto bool_l_s = index_sym[index_l_b].second;
            auto index_r_b = sym_resort[(k+1)%sym_resort.size()];
            auto index_r_s = index_mapping[index_sym[index_r_b].first];
            auto bool_r_s = index_sym[index_r_b].second;

            if(index_sym[index_b].second){//右端
                sym_axi_arr[index_s.first][index_s.second].r_neis.push_back(std::make_pair(index_sym[index_l_b].first,bool_l_s));
                sym_axi_arr[index_s.first][index_s.second].r_neis.push_back(std::make_pair(index_sym[index_r_b].first,bool_r_s));
                
            }else{//左端
                sym_axi_arr[index_s.first][index_s.second].l_neis.push_back(std::make_pair(index_sym[index_l_b].first,bool_l_s));
                sym_axi_arr[index_s.first][index_s.second].l_neis.push_back(std::make_pair(index_sym[index_r_b].first,bool_r_s));

            }       
        }
        
    }
    


}
void CalPointMappingToLine(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& m, Eigen::Vector2d& mapping){
    mapping(0) = ( m(0)*(p2(0)-p1(0))*(p2(0)-p1(0))+m(1)*(p2(1)-p1(1))*(p2(0)-p1(0))+(p1(0)*p2(1)-p2(0)*p1(1))*(p2(1)-p1(1)) )/((p2(0)-p1(0))*(p2(0)-p1(0))+(p2(1)-p1(1))*(p2(1)-p1(1)));
    mapping(1) = ( m(0)*(p2(0)-p1(0))*(p2(1)-p1(1))+m(1)*(p2(1)-p1(1))*(p2(1)-p1(1))+(p2(0)*p1(1)-p1(0)*p2(1))*(p2(0)-p1(0)) )/((p2(0)-p1(0))*(p2(0)-p1(0))+(p2(1)-p1(1))*(p2(1)-p1(1)));
}
//计算在固定方向上点st2i到边界最小距离
//st2i(y,x)左上角起点 st2d(x,y)原点起点(st2d空间) ori单位向量(st2d空间) bdyp最近边界点(st2d空间) 
double ComponentSymmetryAxi::CalOrientationDistance(int min_dis, const Eigen::Vector2i& st2i, const Eigen::Vector2d& st2d, const Eigen::Vector2d& ori, Eigen::Vector2d& bdyp) {
    /*Eigen::Vector2i index_th1(
                    (int)std::abs(sym_axi_arr[i][j].end_coord(1)),
                    se.width -1 - (int)std::abs(sym_axi_arr[i][j].end_coord(0))
                );*/
    //+-min_dis范围内取点
    //int min_dis = 0;
    Eigen::Vector2i st2i_update;
    double min_dis_pp = MAX_LEN;
    double min_dis_pl = MAX_LEN;
    for(int add_h = -min_dis; add_h<=min_dis; ++add_h){
        for(int add_w = -min_dis; add_w<=min_dis; ++add_w){
            st2i_update(0) = st2i(0)+add_h ;
            st2i_update(1) = st2i(1)+add_w ;
            if((st2i_update(0)<0)||(st2i_update(0)>=this->se.height)) continue;
            if((st2i_update(1)<0)||(st2i_update(1)>=this->se.width)) continue;
            if(se.grid_value(st2i_update(0),st2i_update(1))!=0) continue;
            
            Eigen::Vector3d st3i_update(-(se.width-1-st2i_update(1)+0.5),-(st2i_update(0)+0.5),0) ;
            //Eigen::Vector3d st2i_(-(se.width-1-st2i_update(1)+0.5),-(st2i_update(0)+0.5),0) ;st2d
            Eigen::Vector3d st3d(st2d(0),st2d(1),0) ;
            Eigen::Vector3d ori_(ori(0),ori(1),0) ;
            Eigen::Vector3d link = (st3i_update-st3d) ;

            double sinphi = std::abs((link.cross(ori_))(2))/link.norm() ;
            double distance_normal = link.norm()*sinphi ;
            double distance_link = link.norm();
            double distance_cosl = link.dot(ori_) ;

            if((distance_cosl>0)&&(distance_normal<min_dis_pl)){
                bdyp(0) = st3i_update(0);
                bdyp(1) = st3i_update(1);

                min_dis_pp = distance_cosl ;
                min_dis_pl = distance_normal ;
            }
            
            
        }
    }
    
    return min_dis_pp;
}

Eigen::Vector2d ComponentSymmetryAxi::FindOrientation(uint32_t in_, const std::vector<uint32_t>& outs_, const Eigen::Vector2d& in_vec, const std::vector<Eigen::Vector2d>& outs_vec, std::vector<uint32_t>& selected){
    //Eigen::Vector2d mid;
    //co-c
    double an_in_out1 = get2RotateAngle(false,in_vec(0),in_vec(1),outs_vec[0](0),outs_vec[0](1)) ;
    double an_in_out2 = get2RotateAngle(true,in_vec(0),in_vec(1),outs_vec[1](0),outs_vec[1](1)) ;

    //c-co
    if((an_in_out1+an_in_out2)>2*M_PI){
        an_in_out1 = 2*M_PI-an_in_out1 ;
        an_in_out2 = 2*M_PI-an_in_out2 ;
    }

    //判断是否需要延伸长度
    //延伸长度
    if(an_in_out1>M_PI){
        selected.push_back(0) ;
        return (-(outs_vec[0]+in_vec)/2).normalized();
    }else if(an_in_out2>M_PI){
        selected.push_back(1) ;
        return (-(outs_vec[1]+in_vec)/2).normalized();
    }
    return Eigen::Vector2d(0.00,0.00) ;
}
std::pair<double,double> ComponentSymmetryAxi::CalUpdateLenOfInOut(uint32_t outb, uint32_t outa, uint32_t inb, uint32_t ina, bool is_left_out, bool is_left_in){
    if(is_left_out){
        //重新定位位置
        Eigen::Vector2i index_th1(
            (int)std::abs(sym_axi_arr[outb][outa].start_coord(1)),
            se.width -1 - (int)std::abs(sym_axi_arr[outb][outa].start_coord(0))
        );
        auto cloestp1 = *(se.distance_map.cp[index_th1(0)][index_th1(1)].cloestp.begin());
        auto cloestv1_ = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(cloestp1(0), cloestp1(1))));
        auto cloestv1 = Eigen::Vector2d(cloestv1_[0],cloestv1_[1]) ;
        //边界最近点
        auto distance_clo_end = (cloestv1-sym_axi_arr[outb][outa].start_coord).norm();
        //计算对称轴夹角
        Eigen::Vector2d from_vec = (sym_axi_arr[outb][outa].center-sym_axi_arr[outb][outa].start_coord).normalized();
        Eigen::Vector2d to_vec = (sym_axi_arr[inb][ina].center-sym_axi_arr[outb][outa].start_coord).normalized();
        Eigen::Vector2d mid_vec = -((from_vec+to_vec)/2).normalized();

        //夹角重计算
        auto alpha_c = std::acos(from_vec.dot(to_vec)) ;
        auto new_alp1 = alpha_c*0.5;//(cos_l_alp/(cos_l_alp+cos_2_alp));
        auto new_alp2 = alpha_c*0.5;//(cos_2_alp/(cos_l_alp+cos_2_alp));

        Eigen::Vector2d cloestv2(0.00,0.00);
        double distance_clo_end2 = CalOrientationDistance(distance_clo_end,index_th1, sym_axi_arr[outb][outa].start_coord, mid_vec, cloestv2);
        if(distance_clo_end2 == MAX_LEN){
            distance_clo_end2 = distance_clo_end ;
        }

        //半长延伸
        auto r_1_add = distance_clo_end2*std::abs(std::cos(new_alp1)) ;
        auto l_2_add = distance_clo_end2*std::abs(std::cos(new_alp2)) ;
        return std::make_pair(r_1_add, l_2_add) ;
    }else{
        //重新定位位置
        Eigen::Vector2i index_th1(
            (int)std::abs(sym_axi_arr[outb][outa].end_coord(1)),
            se.width -1 - (int)std::abs(sym_axi_arr[outb][outa].end_coord(0))
        );
        auto cloestp1 = *(se.distance_map.cp[index_th1(0)][index_th1(1)].cloestp.begin());
        auto cloestv1_ = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(cloestp1(0), cloestp1(1))));
        auto cloestv1 = Eigen::Vector2d(cloestv1_[0],cloestv1_[1]) ;
        //边界最近点
        auto distance_clo_end = (cloestv1-sym_axi_arr[outb][outa].end_coord).norm();
        //计算对称轴夹角
        Eigen::Vector2d from_vec = (sym_axi_arr[outb][outa].center-sym_axi_arr[outb][outa].end_coord).normalized();
        Eigen::Vector2d to_vec = (sym_axi_arr[inb][ina].center-sym_axi_arr[outb][outa].end_coord).normalized();
        Eigen::Vector2d mid_vec = -((from_vec+to_vec)/2).normalized();

        
        //夹角重计算
        auto alpha_c = std::acos(from_vec.dot(to_vec)) ;
        auto new_alp1 = alpha_c*0.5;//(cos_l_alp/(cos_l_alp+cos_2_alp));
        auto new_alp2 = alpha_c*0.5;//(cos_2_alp/(cos_l_alp+cos_2_alp));

        Eigen::Vector2d cloestv2(0.00,0.00);
        double distance_clo_end2 = CalOrientationDistance(distance_clo_end,index_th1, sym_axi_arr[outb][outa].end_coord, mid_vec, cloestv2);
        if(distance_clo_end2 == MAX_LEN){
            distance_clo_end2 = distance_clo_end ;
        }

        //半长延伸
        auto add1_l = distance_clo_end2*std::abs(std::cos(new_alp1)) ;
        auto add2_lr = distance_clo_end2*std::abs(std::cos(new_alp2)) ;
        return std::make_pair(add1_l,add2_lr) ;
    }
}
//针对外部组件间的半长修正
std::pair<double,double> ComponentSymmetryAxi::CalUpdateLenOf2InOut(uint32_t outb, uint32_t outa, uint32_t inb, uint32_t ina, bool is_left_out, bool is_left_in){
    if(is_left_out){
        //重新定位位置
        Eigen::Vector2i index_th1(
            (int)std::abs(sym_axi_arr[outb][outa].start_coord(1)),
            se.width -1 - (int)std::abs(sym_axi_arr[outb][outa].start_coord(0))
        );
        auto cloestp1 = *(se.distance_map.cp[index_th1(0)][index_th1(1)].cloestp.begin());
        auto cloestv1_ = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(cloestp1(0), cloestp1(1))));
        auto cloestv1 = Eigen::Vector2d(cloestv1_[0],cloestv1_[1]) ;
        //边界最近点
        auto distance_clo_end = (cloestv1-sym_axi_arr[outb][outa].start_coord).norm();
        auto link_to_bdy = (cloestv1-sym_axi_arr[outb][outa].start_coord).normalized();
        //计算对称轴夹角
        Eigen::Vector2d from_vec = (sym_axi_arr[outb][outa].center-sym_axi_arr[outb][outa].start_coord).normalized();
        Eigen::Vector2d to_vec = (sym_axi_arr[inb][ina].center-sym_axi_arr[outb][outa].start_coord).normalized();
        auto cos_1 = std::abs((from_vec).dot(link_to_bdy)) ;
        auto cos_l_alp =  std::acos(cos_1);
        auto cos_2 = std::abs((to_vec).dot(link_to_bdy)) ;
        auto cos_2_alp = std::acos(cos_2);
        auto alpha_c = std::acos(from_vec.dot(to_vec)) ;
        //夹角重计算
        auto new_alp1 = alpha_c*0.5;//(cos_l_alp/(cos_l_alp+cos_2_alp));
        auto new_alp2 = alpha_c*0.5;//(cos_2_alp/(cos_l_alp+cos_2_alp));
        //半长延伸
        auto add1_l = distance_clo_end*std::abs(std::cos(new_alp1)) ;
        auto add2_lr = distance_clo_end*std::abs(std::cos(new_alp2)) ;
        return std::make_pair(add1_l, add2_lr) ;
    }else{
        //重新定位位置
        Eigen::Vector2i index_th1(
            (int)std::abs(sym_axi_arr[outb][outa].end_coord(1)),
            se.width -1 - (int)std::abs(sym_axi_arr[outb][outa].end_coord(0))
        );
        auto cloestp1 = *(se.distance_map.cp[index_th1(0)][index_th1(1)].cloestp.begin());
        auto cloestv1_ = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(cloestp1(0), cloestp1(1))));
        auto cloestv1 = Eigen::Vector2d(cloestv1_[0],cloestv1_[1]) ;
        //边界最近点
        auto distance_clo_end = (cloestv1-sym_axi_arr[outb][outa].end_coord).norm();
        auto link_to_bdy = (cloestv1-sym_axi_arr[outb][outa].end_coord).normalized();
        //计算对称轴夹角
        Eigen::Vector2d from_vec = (sym_axi_arr[outb][outa].center-sym_axi_arr[outb][outa].end_coord).normalized();
        Eigen::Vector2d to_vec = (sym_axi_arr[inb][ina].center-sym_axi_arr[outb][outa].end_coord).normalized();
        auto cos_1 = std::abs((from_vec).dot(link_to_bdy)) ;
        auto cos_l_alp =  std::acos(cos_1);
        auto cos_2 = std::abs((to_vec).dot(link_to_bdy)) ;
        auto cos_2_alp = std::acos(cos_2);
        auto alpha_c = std::acos(from_vec.dot(to_vec)) ;
        //夹角重计算
        auto new_alp1 = alpha_c*0.5;//(cos_l_alp/(cos_l_alp+cos_2_alp));
        auto new_alp2 = alpha_c*0.5;//(cos_2_alp/(cos_l_alp+cos_2_alp));
        //半长延伸
        auto add1_l = distance_clo_end*std::abs(std::cos(new_alp1)) ;
        auto add2_lr = distance_clo_end*std::abs(std::cos(new_alp2)) ;
        return std::make_pair(add1_l,add2_lr) ;
    }
}
//针对内部组件间的半长修正
std::pair<double,double> ComponentSymmetryAxi::CalUpdateLenOf2In(uint32_t b1, uint32_t a1, uint32_t b2, uint32_t a2){

    Eigen::Vector2i index_th1 ;
    index_th1(0) = (int)( std::abs( sym_axi_arr[b1][a1].end_coord(1) ));
    index_th1(1) = se.width -1 - (int)( std::abs( sym_axi_arr[b1][a1].end_coord(0) ));
    //重新定位位置
    auto cloestp1 = *(se.distance_map.cp[index_th1(0)][index_th1(1)].cloestp.begin());
    auto cloestv1_ = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(cloestp1(0), cloestp1(1))));
    auto cloestv1 = Eigen::Vector2d(cloestv1_[0],cloestv1_[1]) ;

    double distance_clo_end = (cloestv1-sym_axi_arr[b1][a1].end_coord).norm();

    
    //计算对称轴夹角
    Eigen::Vector2d from_vec = (sym_axi_arr[b1][a1].center-sym_axi_arr[b1][a1].end_coord).normalized();
    Eigen::Vector2d to_vec = (sym_axi_arr[b2][a2].center-sym_axi_arr[b2][a2].start_coord).normalized();
    Eigen::Vector2d mid_vec = -((from_vec+to_vec)/2.00).normalized();
    double alpha_c;
    alpha_c = std::acos(from_vec.dot(to_vec)) ;
    

    double new_alp1,new_alp2;
    new_alp1 = alpha_c/2.00;
    new_alp2 = alpha_c/2.00;

    Eigen::Vector2d cloestv2(0.00,0.00);
    double distance_clo_end2 = CalOrientationDistance(distance_clo_end,index_th1, sym_axi_arr[b1][a1].end_coord, mid_vec, cloestv2);
    if(distance_clo_end2 == MAX_LEN){
        distance_clo_end2 = distance_clo_end ;
    } 

    /*
    if(cos_1_alp<M_PI){
        alpha_c = M_PI - (cos_1_alp)/2;
    }else{
        alpha_c = (cos_1_alp)/2;
    }*/
    
    //alpha_c = (cos_1_alp)/2.00;
    
    //alpha_c = std::cos(alpha_c) ;

    auto r_1_add = distance_clo_end2*std::abs(std::cos(new_alp1)) ;
    auto l_2_add = distance_clo_end2*std::abs(std::cos(new_alp2)) ;
    //std::cout<<"\n cos:"<<cos_1<<" alpha:"<<cos_1_alp/M_PI*180<<" half_cos:"<<alpha_c<<" r_1_add:"<<r_1_add<<" l_2_add:"<<l_2_add;
    //std::cout<<"\n alp1:"<<new_alp1/M_PI*180<<" alp2:"<<new_alp2/M_PI*180<<" alphac:"<<alpha_c/M_PI*180<<" r_1_add:"<<r_1_add<<" l_2_add:"<<l_2_add;
            
    return std::make_pair(r_1_add,l_2_add) ;
}
void ComponentSymmetryAxi::UpdateHalfLens(){
    
    
    //更新带标示的内部对称轴的半长
    for(uint32_t i=0; i<se.branches.size() ; ++i){
        std::cout<<"\nupdate branch:"<<i;
        for(uint32_t j = 0 ; j<pos_split[i].size()-2; ++j){
            //检查相邻两对称轴是否有标记
            if(!(sym_axi_arr[i][j].r_update&&sym_axi_arr[i][j+1].l_update)){
                //调整近似平行组件的半长
                if((!sym_axi_arr[i][j].r_update)&&(!sym_axi_arr[i][j+1].l_update))
                {
                    sym_axi_arr[i][j].center = sym_axi_arr[i][j].center+4*sym_axi_arr[i][j].ori;
                    sym_axi_arr[i][j+1].center = sym_axi_arr[i][j+1].center-4*sym_axi_arr[i][j+1].ori;
                    
                    sym_axi_arr[i][j].half_len = sym_axi_arr[i][j].half_len+(4);
                    sym_axi_arr[i][j+1].half_len = sym_axi_arr[i][j+1].half_len+(4);

                    sym_axi_arr[i][j].end_coord = sym_axi_arr[i][j].center+sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori;
                    sym_axi_arr[i][j+1].start_coord = sym_axi_arr[i][j+1].center-sym_axi_arr[i][j+1].half_len*sym_axi_arr[i][j+1].ori;
                }
                continue;
            }
            
            auto pair_add = CalUpdateLenOf2In(i,j,i,j+1) ;
            auto r_1_add = pair_add.first;
            auto l_2_add = pair_add.second;

            sym_axi_arr[i][j].center = sym_axi_arr[i][j].center+(r_1_add/2)*sym_axi_arr[i][j].ori;
            sym_axi_arr[i][j+1].center = sym_axi_arr[i][j+1].center-(l_2_add/2)*sym_axi_arr[i][j+1].ori;
            
            sym_axi_arr[i][j].half_len = sym_axi_arr[i][j].half_len+(r_1_add/2);
            sym_axi_arr[i][j+1].half_len = sym_axi_arr[i][j+1].half_len+(l_2_add/2);

            sym_axi_arr[i][j].end_coord = sym_axi_arr[i][j].center+sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori;
            sym_axi_arr[i][j+1].start_coord = sym_axi_arr[i][j+1].center-sym_axi_arr[i][j+1].half_len*sym_axi_arr[i][j+1].ori;
            
            //添加存在约束的邻近组件
            sym_axi_arr[i][j].cons_neis.insert(sym_axi_arr[i][j+1].index) ;
            sym_axi_arr[i][j+1].cons_neis.insert(sym_axi_arr[i][j].index) ;

        }
    }
    
    //修正交叉点附近组件对称轴半长（包括pi以下夹角与pi以上夹角的处理）
    for(auto& syms_pi : syms_jun_all){
        auto& index_p = syms_pi.first;
        auto& index_lr = syms_pi.second;
        //需要直接延长求交的对称轴 加上约束标示
        std::pair<uint32_t,uint32_t> dull2;
        std::pair<bool,bool> dull2lr;
        std::unordered_map<uint32_t,Eigen::Vector2d> all_oris;
        bool find_dull2 = false;
        //收集方向
        for(auto& iter_i : index_lr){
            auto index_s = iter_i.first;
            auto lr_ = iter_i.second;
            auto re_s = index_mapping[index_s];
            all_oris[index_s] = (lr_?-sym_axi_arr[re_s.first][re_s.second].ori:sym_axi_arr[re_s.first][re_s.second].ori);
        }
        //寻找钝角
        for(auto& iter_i : index_lr){
            auto index_s = iter_i.first;
            auto lr_ = iter_i.second;
            auto remap_s = index_mapping[index_s];
            //邻接2
            auto& neis = (lr_?sym_axi_arr[remap_s.first][remap_s.second].r_neis:sym_axi_arr[remap_s.first][remap_s.second].l_neis) ;
            std::vector<Eigen::Vector2d> nei2ori;
            std::vector<uint32_t> nei2in;
            for(auto neii:neis){
                nei2ori.push_back(all_oris[neii.first]) ;
                nei2in.push_back(neii.first);
            }
            std::vector<uint32_t> selected;
            FindOrientation(index_s,nei2in,all_oris[index_s],nei2ori,selected);
            if(selected.size()==0){
                find_dull2 = false;
                continue;
            }else{
                dull2.first = index_s;
                dull2.second = neis[selected[0]].first;
                dull2lr.first = lr_;

                for(auto& c : index_lr){
                    if(c.first == dull2.second){
                        dull2lr.second = c.second;
                        break;
                    }
                }
                
                find_dull2 = true;
                break;
            }
        }
        
        //需要交叉点投影的对称轴
        auto jun_co = se.ske_points[index_p].coord ;
        std::cout<<"\nall oris:"<<all_oris.size()<<"\n" ;
        double add_len = 0;
        for(auto iiter : all_oris){
            auto ps1 = index_mapping[iiter.first];
            auto len1 = (sym_axi_arr[ps1.first][ps1.second].start_coord-jun_co).norm();
            auto len2 = (sym_axi_arr[ps1.first][ps1.second].end_coord-jun_co).norm();
            add_len += (len2>len1?len1:len2) ;
        }
        add_len/=all_oris.size() ;
        add_len*=0.2;
        for(auto iiter : all_oris){
            auto ps1 = index_mapping[iiter.first];
            //投影点
            Eigen::Vector2d mapping_p;
            CalPointMappingToLine(sym_axi_arr[ps1.first][ps1.second].start_coord,sym_axi_arr[ps1.first][ps1.second].end_coord,jun_co,mapping_p) ;
            
            std::cout<<"before--- index:"<<iiter.first<<" ps1:"<<ps1.first<<","<<ps1.second<<" center:"<<sym_axi_arr[ps1.first][ps1.second].center(0)<<","<<sym_axi_arr[ps1.first][ps1.second].center(1);
            std::cout<<" hl:"<<sym_axi_arr[ps1.first][ps1.second].half_len;
            std::cout<<" s:"<<sym_axi_arr[ps1.first][ps1.second].start_coord(0)<<","<<sym_axi_arr[ps1.first][ps1.second].start_coord(1);
            std::cout<<" e:"<<sym_axi_arr[ps1.first][ps1.second].end_coord(0)<<","<<sym_axi_arr[ps1.first][ps1.second].end_coord(1)<<"\n";

            Eigen::Vector2d start_p = sym_axi_arr[ps1.first][ps1.second].center+sym_axi_arr[ps1.first][ps1.second].half_len*iiter.second;
            std::cout<<"start_p:"<<start_p(0)<<","<<start_p(1)<<" mapping:"<<mapping_p(0)<<","<<mapping_p(1)
            <<" newori:"<< (mapping_p-start_p).normalized()(0)<<","<<(mapping_p-start_p).normalized()(1)<<"\n";

            //更新分支
            sym_axi_arr[ps1.first][ps1.second].center = (start_p+mapping_p)/2;
            sym_axi_arr[ps1.first][ps1.second].half_len = (start_p-mapping_p).norm()/2;
            sym_axi_arr[ps1.first][ps1.second].start_coord = sym_axi_arr[ps1.first][ps1.second].center-sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;
            sym_axi_arr[ps1.first][ps1.second].end_coord = sym_axi_arr[ps1.first][ps1.second].center+sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;
            std::cout<<"after---"<<" center:"<<sym_axi_arr[ps1.first][ps1.second].center(0)<<","<<sym_axi_arr[ps1.first][ps1.second].center(1);
            std::cout<<" hl:"<<sym_axi_arr[ps1.first][ps1.second].half_len;
            std::cout<<" s:"<<sym_axi_arr[ps1.first][ps1.second].start_coord(0)<<","<<sym_axi_arr[ps1.first][ps1.second].start_coord(1);
            std::cout<<" e:"<<sym_axi_arr[ps1.first][ps1.second].end_coord(0)<<","<<sym_axi_arr[ps1.first][ps1.second].end_coord(1)<<"\n";
        
            //对无边界约束的组件对称轴为了鲁棒性延长部分长度
            if(find_dull2&&(iiter.first==dull2.first||iiter.first==dull2.second))continue;
            
            sym_axi_arr[ps1.first][ps1.second].center = sym_axi_arr[ps1.first][ps1.second].center-add_len*iiter.second;
            sym_axi_arr[ps1.first][ps1.second].half_len = sym_axi_arr[ps1.first][ps1.second].half_len+add_len;
            sym_axi_arr[ps1.first][ps1.second].start_coord = sym_axi_arr[ps1.first][ps1.second].center-sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;
            sym_axi_arr[ps1.first][ps1.second].end_coord = sym_axi_arr[ps1.first][ps1.second].center+sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;
        }
        

        //处理交叉点附近可能存在边界连接的对称轴
        if(find_dull2){
            auto pa1 = dull2.first;
            auto pa2 = dull2.second;
            auto ps1 = index_mapping[pa1];
            auto ps2 = index_mapping[pa2];
            Eigen::Vector2d ise(0,0);
            double t1,t2;
            //TODO:检查t1 t2符号
            /*
            Compute2LinesIntersection(ise,sym_axi_arr[ps1.first][ps1.second].center,sym_axi_arr[ps2.first][ps2.second].center,all_oris[pa1],all_oris[pa2],t1,t2);
            if(t1*t2>=0){//可直接连接后计算
                Eigen::Vector2d p11 = sym_axi_arr[ps1.first][ps1.second].center+sym_axi_arr[ps1.first][ps1.second].half_len*all_oris[pa1];
                Eigen::Vector2d p12 = sym_axi_arr[ps1.first][ps1.second].center+t1*all_oris[pa1];
                sym_axi_arr[ps1.first][ps1.second].center = (p11+p12)/2.00;
                sym_axi_arr[ps1.first][ps1.second].half_len = (p12-p11).norm()/2.00;
                sym_axi_arr[ps1.first][ps1.second].start_coord = sym_axi_arr[ps1.first][ps1.second].center-sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;
                sym_axi_arr[ps1.first][ps1.second].end_coord = sym_axi_arr[ps1.first][ps1.second].center+sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;


                Eigen::Vector2d p21 = sym_axi_arr[ps2.first][ps2.second].center+t2*all_oris[pa2];
                Eigen::Vector2d p22 = sym_axi_arr[ps2.first][ps2.second].center+sym_axi_arr[ps2.first][ps2.second].half_len*all_oris[pa2];
                sym_axi_arr[ps2.first][ps2.second].center = (p21+p22)/2.00;
                sym_axi_arr[ps2.first][ps2.second].half_len = (p22-p21).norm()/2.00;
                sym_axi_arr[ps2.first][ps2.second].start_coord = sym_axi_arr[ps2.first][ps2.second].center-sym_axi_arr[ps2.first][ps2.second].half_len*sym_axi_arr[ps2.first][ps2.second].ori;
                sym_axi_arr[ps2.first][ps2.second].end_coord = sym_axi_arr[ps2.first][ps2.second].center+sym_axi_arr[ps2.first][ps2.second].half_len*sym_axi_arr[ps2.first][ps2.second].ori;

                sym_axi_arr[ps1.first][ps1.second].cons_neis.insert(pa2) ;
                sym_axi_arr[ps2.first][ps2.second].cons_neis.insert(pa1) ;


                //计算延长长度
                Eigen::Vector2i index_th1 ;
                index_th1(0) = (int)( std::abs( ise(1) ));
                index_th1(1) = se.width -1 - (int)( std::abs( ise(0) ));
                //重新定位位置
                auto cloestp1 = *(se.distance_map.cp[index_th1(0)][index_th1(1)].cloestp.begin());
                auto cloestv1_ = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(cloestp1(0), cloestp1(1))));
                auto cloestv1 = Eigen::Vector2d(cloestv1_[0],cloestv1_[1]) ;
                double distance_clo_end = (cloestv1-ise).norm();
                Eigen::Vector2d from_vec = all_oris[pa1];
                Eigen::Vector2d to_vec = all_oris[pa2];
                Eigen::Vector2d mid_vec = -((from_vec+to_vec)/2.00).normalized();
                Eigen::Vector2d cloestv2(0.00,0.00);
                double distance_clo_end2 = CalOrientationDistance(distance_clo_end,index_th1, ise, mid_vec, cloestv2);
                if(distance_clo_end2 == MAX_LEN){
                    distance_clo_end2 = distance_clo_end ;
                } 
                auto alpha_c = std::acos(from_vec.dot(to_vec)) ;
                auto new_alp1 = alpha_c*0.5;//(cos_l_alp/(cos_l_alp+cos_2_alp));
                auto new_alp2 = alpha_c*0.5;//(cos_2_alp/(cos_l_alp+cos_2_alp));
                auto r_1_add = distance_clo_end2*std::abs(std::cos(new_alp1)) ;
                auto l_2_add = distance_clo_end2*std::abs(std::cos(new_alp2)) ;

                //更新分支
                sym_axi_arr[ps1.first][ps1.second].center = sym_axi_arr[ps1.first][ps1.second].center - r_1_add/2.00*all_oris[pa1];
                sym_axi_arr[ps1.first][ps1.second].half_len = sym_axi_arr[ps1.first][ps1.second].half_len + r_1_add/2.00;
                sym_axi_arr[ps1.first][ps1.second].start_coord = sym_axi_arr[ps1.first][ps1.second].center-sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;
                sym_axi_arr[ps1.first][ps1.second].end_coord = sym_axi_arr[ps1.first][ps1.second].center+sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;

                sym_axi_arr[ps2.first][ps2.second].center = sym_axi_arr[ps2.first][ps2.second].center - l_2_add/2.00*all_oris[pa2];
                sym_axi_arr[ps2.first][ps2.second].half_len = sym_axi_arr[ps2.first][ps2.second].half_len + l_2_add/2.00;
                sym_axi_arr[ps2.first][ps2.second].start_coord = sym_axi_arr[ps2.first][ps2.second].center-sym_axi_arr[ps2.first][ps2.second].half_len*sym_axi_arr[ps2.first][ps2.second].ori;
                sym_axi_arr[ps2.first][ps2.second].end_coord = sym_axi_arr[ps2.first][ps2.second].center+sym_axi_arr[ps2.first][ps2.second].half_len*sym_axi_arr[ps2.first][ps2.second].ori;


                //删除已经处理过的两个分支
                auto iter = all_oris.find(pa1) ;
                all_oris.erase(iter) ;
                auto iter2 = all_oris.find(pa2) ;
                all_oris.erase(iter2) ;
            }else{//不可直接连接时，先做投影，继而计算
               */
                sym_axi_arr[ps1.first][ps1.second].cons_neis.insert(pa2) ;
                sym_axi_arr[ps2.first][ps2.second].cons_neis.insert(pa1) ;

                //计算延长长度
                Eigen::Vector2i index_th1 ;
                index_th1(0) = (int)( std::abs( jun_co(1) ));
                index_th1(1) = se.width -1 - (int)( std::abs( jun_co(0) ));
                //重新定位位置
                auto cloestp1 = *(se.distance_map.cp[index_th1(0)][index_th1(1)].cloestp.begin());
                auto cloestv1_ = se.mesh.calc_face_centroid(OpenMesh::FaceHandle(se.grid_index(cloestp1(0), cloestp1(1))));
                auto cloestv1 = Eigen::Vector2d(cloestv1_[0],cloestv1_[1]) ;
                double distance_clo_end = (cloestv1-jun_co).norm();
                Eigen::Vector2d from_vec = all_oris[pa1];
                Eigen::Vector2d to_vec = all_oris[pa2];
                Eigen::Vector2d mid_vec = -((from_vec+to_vec)/2.00).normalized();
                Eigen::Vector2d cloestv2(0.00,0.00);
                double distance_clo_end2 = CalOrientationDistance(distance_clo_end,index_th1, jun_co, mid_vec, cloestv2);
                if(distance_clo_end2 == MAX_LEN){
                    distance_clo_end2 = distance_clo_end ;
                } 
                auto alpha_c = std::acos(from_vec.dot(to_vec)) ;
                auto new_alp1 = alpha_c*0.5;//(cos_l_alp/(cos_l_alp+cos_2_alp));
                auto new_alp2 = alpha_c*0.5;//(cos_2_alp/(cos_l_alp+cos_2_alp));
                auto r_1_add = distance_clo_end2*std::abs(std::cos(new_alp1)) ;
                auto l_2_add = distance_clo_end2*std::abs(std::cos(new_alp2)) ;

                //更新分支
                sym_axi_arr[ps1.first][ps1.second].center = sym_axi_arr[ps1.first][ps1.second].center - r_1_add/2.00*all_oris[pa1];
                sym_axi_arr[ps1.first][ps1.second].half_len = sym_axi_arr[ps1.first][ps1.second].half_len + r_1_add/2.00;
                sym_axi_arr[ps1.first][ps1.second].start_coord = sym_axi_arr[ps1.first][ps1.second].center-sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;
                sym_axi_arr[ps1.first][ps1.second].end_coord = sym_axi_arr[ps1.first][ps1.second].center+sym_axi_arr[ps1.first][ps1.second].half_len*sym_axi_arr[ps1.first][ps1.second].ori;

                sym_axi_arr[ps2.first][ps2.second].center = sym_axi_arr[ps2.first][ps2.second].center - l_2_add/2.00*all_oris[pa2];
                sym_axi_arr[ps2.first][ps2.second].half_len = sym_axi_arr[ps2.first][ps2.second].half_len + l_2_add/2.00;
                sym_axi_arr[ps2.first][ps2.second].start_coord = sym_axi_arr[ps2.first][ps2.second].center-sym_axi_arr[ps2.first][ps2.second].half_len*sym_axi_arr[ps2.first][ps2.second].ori;
                sym_axi_arr[ps2.first][ps2.second].end_coord = sym_axi_arr[ps2.first][ps2.second].center+sym_axi_arr[ps2.first][ps2.second].half_len*sym_axi_arr[ps2.first][ps2.second].ori;

                /*
                //删除已经处理过的两个分支
                auto iter = all_oris.find(pa1) ;
                all_oris.erase(iter) ;
                auto iter2 = all_oris.find(pa2) ;
                all_oris.erase(iter2) ;*/
            //}
        }
        
    }
    
    //DisplayField(FieldType::Com_SYM_LINES) ;
}
/*
void ComponentSymmetryAxi::DisplayField(FieldType type_field){
    #ifdef GUI_CSA_ON
    if(type_field==FieldType::Com_SYM_LINES){
        DataVisualization data;
		std::vector<std::vector<float>> x, y;
		for (uint32_t i = 0; i < sym_axi_arr.size(); ++i)
		{
            for(uint32_t j = 0 ; j<sym_axi_arr[i].size() ; ++j){
                std::vector<float> x_(2), y_(2);
			
                Eigen::Vector2d p1 = sym_axi_arr[i][j].center+sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori;
                Eigen::Vector2d p2 = sym_axi_arr[i][j].center-sym_axi_arr[i][j].half_len*sym_axi_arr[i][j].ori;
                
                x_[0] = p1(0) ;
                x_[1] = p2(0) ;
                y_[0] = p1(1) ;
                y_[1] = p2(1) ;
                
                x.push_back(x_);
                y.push_back(y_);
            }
			
		}
		data.ShowSegLines("SegLines", x, y);
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