#include "TimeIntegrator.h"
namespace Poster_Needed {
    timeIntegrator::timeIntegrator(stiffAssembler &stiff_, massAssembler &_mass)
            : stiff(stiff_), mass(_mass) {
        defaultOptions();
    }

    void timeIntegrator::setInitialDisplacement(const Eigen::MatrixXd &initialDisplacement) {
        dispVector = initialDisplacement;
        pre_dispVector = dispVector;
    }

    void timeIntegrator::setInitialVelocity(const Eigen::MatrixXd &initialVelocity) {
        velVector = initialVelocity;
        pre_velVector = velVector;
    }

    void timeIntegrator::initialize() {
        mass.assemble();
        if (dispVector.rows() != numDofs())
            dispVector.setZero(stiff.numDofs(), 1);
        if (velVector.rows() != numDofs())
            velVector.setZero(stiff.numDofs(), 1);
        stiff.assemble(dispVector, false);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(mass.matrix());
        accVector = solver.solve(stiff.getRhs());
        pre_accVector = accVector;
    }

    void timeIntegrator::makeTimeStep(double timeStep) {
        tStep = timeStep;
        Eigen::MatrixXd newDispVector;
        if (Scheme == time_integration::implicit_linear)
            newDispVector = implicitLinear();
        if (Scheme == time_integration::implicit_nonlinear)
            newDispVector = implicitNonlinear();

        Eigen::MatrixXd tempVelVector = velVector;
        velVector = alpha4() * (newDispVector - dispVector) + alpha5() * tempVelVector + alpha6() * accVector;
        accVector = alpha1() * (newDispVector - dispVector) - alpha2() * tempVelVector - alpha3() * accVector;
        dispVector = newDispVector;

    }


    void timeIntegrator::assemble(const Eigen::VectorXd &solution, bool isAssemblerMatrix) {
        stiff.assemble(solution, isAssemblerMatrix);
        if (isAssemblerMatrix) {
            m_system.matrix() = alpha1() * mass.matrix() + stiff.getMatrix();
            m_system.matrix().makeCompressed();
        }
        m_system.rhs() = stiff.getRhs() + mass.matrix() * (alpha1() * (dispVector - solution) + alpha2() * velVector +
                                                           alpha3() * accVector);
        //if (more_fixed) {
        //	m_system.matrix() += newMatrix;
        //	m_system.rhs() += newRhs;
        //}
    }

    const Eigen::MatrixXd &timeIntegrator::displacementVector() const {
        return dispVector;
    }

    const Eigen::MatrixXd &timeIntegrator::velocityVector() const {
        return velVector;
    }

    void timeIntegrator::defaultOptions() {
        Beta = 0.25;
        Gamma = 0.5;
        Scheme = time_integration::implicit_nonlinear;
        Verbosity = newton_verbosity::none;
    }

    Eigen::MatrixXd timeIntegrator::implicitLinear() {
        m_system.matrix() = alpha1() * mass.matrix() + stiff.getMatrix();
        matrix.makeCompressed();
        m_system.rhs() =
                mass.matrix() * (alpha1() * dispVector + alpha2() * velVector + alpha3() * accVector) + stiff.getRhs();
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(matrix);
        return solver.solve(rhs);
    }

    Eigen::MatrixXd timeIntegrator::implicitNonlinear() {
        Hyperelastic solver(this, dispVector);
        solver.setInfoNone();
        solver.solve();
        return solver.solution();
    }
}