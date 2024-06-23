#pragma once
#include"Hyperelastic.h"
#include"massAssembler.h"
namespace Poster_Needed {
    class timeIntegrator : public stiffAssembler {
        massAssembler &mass;           //质量矩阵计算
        stiffAssembler &stiff;
        double tStep;                //时间步长
        Eigen::MatrixXd dispVector, velVector, accVector;     //位移、速度、加速度
        Eigen::MatrixXd pre_dispVector, pre_velVector, pre_accVector;     //上一步的位移、速度、加速度
        Eigen::SparseMatrix<double> matrix;
        Eigen::VectorXd rhs;
        double Beta, Gamma;
        int Scheme, Verbosity;
        //bool more_fixed = false;
        //Eigen::MatrixXd newRhs;
        //Eigen::SparseMatrix<double> newMatrix;
    public:
        timeIntegrator(stiffAssembler &stiff_, massAssembler &_mass);

        //设置初始位移与速度
        void setInitialDisplacement(const Eigen::MatrixXd &initialDisplacement);

        void setInitialVelocity(const Eigen::MatrixXd &initialVelocity);

        /// 获得初始状态加速度
        void initialize();

        /// 设置时间步长
        void makeTimeStep(double timeStep);

        /// 装配线性系统);
        void assemble(const Eigen::VectorXd &solution, bool isAssemblerMatrix = true);

        const Eigen::MatrixXd &displacementVector() const; //位移向量

        const Eigen::MatrixXd &velocityVector() const; //加速度向量

        void defaultOptions();     //设置默认参数


    protected:
        Eigen::MatrixXd implicitLinear();

        Eigen::MatrixXd implicitNonlinear();

        double alpha1() { return 1. / Beta / pow(tStep, 2); }

        double alpha2() { return 1. / Beta / tStep; }

        double alpha3() { return (1 - 2 * Beta) / 2 / Beta; }

        double alpha4() { return Gamma / Beta / tStep; }

        double alpha5() { return 1 - Gamma / Beta; }

        double alpha6() { return (1 - Gamma / Beta / 2) * tStep; }
    };
}