#ifndef __COMPONENTSYMMETRYAXI_H__
#define __COMPONENTSYMMETRYAXI_H__
#include"common.h"
#include"DistanceMap.h"
#include"SkeletonEx.h"
#include"pointcloud.h"


struct SymAxi
{
    public:
    uint32_t index;
    uint32_t index_bran ;
    uint32_t index_part ;
    Eigen::Vector2d center;
    Eigen::Vector2d ori;
    double half_len;
    //首尾节点对应骨架分支线和双边对称线上节点序号（每个分支0起始）
    Eigen::Vector2i end_p2;
    //首尾节点所在分支
    Eigen::Vector2i end_p2_brani;
    //首尾节点对应骨架分支线上骨架点序号
    Eigen::Vector2i end_p2_skpi;
    //首尾节点坐标
    Eigen::Vector2d start_coord;
    Eigen::Vector2d end_coord;
    //是否需要设置延长半径
    bool l_update = false;
    bool r_update = false;
    //是否内部对称轴
    bool inner = true;

    //左右邻接组件(pair-<组件号,邻接对称轴的首点或尾点（非本对称轴）>)
    std::vector<std::pair<uint32_t,uint32_t>> l_neis;
    std::vector<std::pair<uint32_t,uint32_t>> r_neis;
    //左右邻接组件的约束关系
    std::set<uint32_t> cons_neis;
};


class ComponentSymmetryAxi
{
private:
    /* data */
    SkeletonEx& se;
    std::vector<std::vector<uint32_t>>& pos_split;

    std::vector<std::pair<uint32_t,std::vector<Eigen::Vector2d>>>& branchClusteres;
    std::vector<Eigen::Vector2d> pClusteres;
    OBBBoundingBox_Extend obb;
    
public:
    //完备的对称轴提取
    std::vector<std::vector<SymAxi>> sym_axi_arr;
    //序号对应 index-bran_i+index_i
    std::vector<std::pair<uint32_t,uint32_t>> index_mapping;
    //******************************针对通过分解分支线得到的对称轴**************************
    //分支线对应分解对称轴的首尾点
    std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>> sym_axi_ends;
    //连接点在首尾对称轴的起始或终止点
    std::vector<std::pair<uint32_t,uint32_t>> sym_axi_end_lr;
    //std::vector<std::pair<uint32_t,uint32_t>> sym_axi_end_i;
    
    //交叉点附近对称轴信息 skepi-<<symi,l or r>,>
    std::unordered_map<uint32_t,std::vector<std::pair<uint32_t,bool>>> syms_jun_all;

    ComponentSymmetryAxi(SkeletonEx& se, std::vector<std::vector<uint32_t>>& pos_split, std::vector<std::pair<uint32_t,std::vector<Eigen::Vector2d>>>& branchClusteres):se(se),pos_split(pos_split),branchClusteres(branchClusteres){
        
    };
    ~ComponentSymmetryAxi(){

    };
    void OBBBranch(uint32_t i, uint32_t i_now);
    bool CheckBranchSymAxis(uint32_t i, bool only_e2e);
    void CalOBBSymAxis() ;
    double CalOrientationDistance(int min_dis, const Eigen::Vector2i& st2i, const Eigen::Vector2d& st2d, const Eigen::Vector2d& ori,Eigen::Vector2d& bdyp) ;
    Eigen::Vector2d FindOrientation(uint32_t in_, const std::vector<uint32_t>& outs_, const Eigen::Vector2d& in_vec, const std::vector<Eigen::Vector2d>& outs_vec, std::vector<uint32_t>& selected);
    std::pair<double,double> CalUpdateLenOf2In(uint32_t b1, uint32_t a1, uint32_t b2, uint32_t a2) ;
    std::pair<double,double> CalUpdateLenOfInOut(uint32_t outb, uint32_t outa, uint32_t inb, uint32_t ina, bool is_left_out, bool is_left_in);
    std::pair<double,double> CalUpdateLenOf2InOut(uint32_t outb, uint32_t outa, uint32_t inb, uint32_t ina, bool is_left_out, bool is_left_in);
    void UpdateHalfLens() ;
    //void DisplayField(FieldType type_field);
};



#endif