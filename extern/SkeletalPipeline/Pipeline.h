#ifndef __LIB_PIPELINE_H__
#define __LIB_PIPELINE_H__
//#ifdef Handle
#undef Handle


#include "common.h"
//#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS
//#define _SILENCE_CXX17_NEGATORS_DEPRECATION_WARNING

#include "SkeletonEx.h"
#include "PointClustersCal.h"
#include "ComponentSymmetryAxi.h"
#include "mmc.h"
#include "ShapeApproximate.h"
//namespace plt = matplotlibcpp;

class Pipeline
{
private:
	/* data */
	SkeletonEx* ske;
	PointClusters* pc;
	ComponentSymmetryAxi* ca;
	ShapeApproximate* sat ;
	//面属性 1-3各类型骨架点
	const std::string SKE_STR = "ske_flag";
	//面属性 0-2一般网格单元 3规则骨架点 10分段点
	const std::string SKE_SP_STR = "ske_sp_flag";
	//点属性 double隐函数值
	const std::string MMC_VAL_STR = "mmc_p_val" ;
public:
	/*BuildSkeletonToFile接口建议默认参数，不要瞎改····*/
	Pipeline(){};

	// 网格文件地址(全路径带后缀.obj) 极小圆采样数 通量场阀值 边界阀值 背景网格宽长
	polyMesh& BuildSkeletonToFile(std::string meshPath, uint32_t width, uint32_t height, uint32_t samplesCir=60, double fluxThreshold=-0.41, double epsilon=1);
	//提取分段点
	polyMesh& SparseOptimization(uint32_t iterationTimes=1, double dMaxError=1.0 );
	//生成组件对称轴
	void GenerateSymAxis(std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>>& axis) ;
	//形状逼近
	//设置为0可以显示初值，演示的话不优化直接用0省时间，平均要10000左右数值上有很好的结果，4000尚可
	polyMesh& ShapeApproximation(uint32_t maxIterations = 4000) ;
	
	//提取组件精确边界
	void ExtractMMCBoundary(std::vector<std::vector<Eigen::Vector2d>>& bdy,double sampleLen=0.1);
	~Pipeline(){

	};
};
//#endif // Handle
#endif

