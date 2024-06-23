#include "GP_BoundaryCPoint.h"

void helper_provide_neumann_dirichlet_points(std::vector<boundary_point>& neumann_points, std::vector<boundary_point>& dirichlet_points)
{

    CVector traction(0, 1500, -1500);
    CVector traction1(10000, 0, 0);
    CVector traction2(3000, -3000, 0);

    std::vector<boundary_point> neumann_tmp, dirichlet_tmp;

    int mode = 1;
    switch (mode) {
        case 1:
            //neumann
            for (int i = 0; i <  100; ++i)
            {
                neumann_tmp.push_back(boundary_point(condition_type::neumann, traction1, i));
            }


            //dirichlet
            for (int i = 100; i < 200; ++i)
            {
                dirichlet_tmp.push_back(boundary_point(condition_type::dirichlet, CVector(0,0,0), i));
            }

            break;
        case 2:
            break;
        case 3:
            break;
        default:
            break;
    };

    neumann_points.swap(neumann_tmp);
    dirichlet_points.swap(dirichlet_tmp);
}