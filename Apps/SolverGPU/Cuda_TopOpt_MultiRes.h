#pragma once

#include <string>
#include <execution>
#include <numeric>
#include "Cuda_Linear_Elastic_Solver.h"
#include "Hyperelasticity/CBSplineSolid.h"
#include<atomic>

struct FineElement {
    double energy;     // 应变能密度
    double volume;     // 体积
    double sensitivity; // 灵敏度

    FineElement() : energy(0.0), volume(1.0), sensitivity(0.0) {}
};

struct MPoint{
    double x, y, z;
};

template <typename T>
struct PointCloud
{
    using coord_t = T;
    std::vector<MPoint> pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline T kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if (dim == 0)
            return pts[idx].x;
        else if (dim == 1)
            return pts[idx].y;
        else
            return pts[idx].z;
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
        return false;
    }
};

struct MultiResMapping {
    std::vector<std::vector<int>> coarse_to_fine;
    std::vector<int> fine_to_coarse;
    std::vector<double> fine_volumes;
    int coarse_size, fine_size, subdiv_factor;

    bool init(int coarse_sz, int subdiv_level) {
        if(coarse_sz <= 0 || subdiv_level <= 0 || subdiv_level > 5) {
            std::cerr << "Invalid parameters: coarse_sz=" << coarse_sz
                      << ", subdiv_level=" << subdiv_level << std::endl;
            return false;
        }

        coarse_size = coarse_sz;
        subdiv_factor = std::pow(8, subdiv_level);
        fine_size = coarse_size * subdiv_factor;

        try {
            coarse_to_fine.clear();
            coarse_to_fine.resize(coarse_size);
            fine_to_coarse.clear();
            fine_to_coarse.resize(fine_size);
            fine_volumes.clear();
            fine_volumes.resize(fine_size, 1.0);

            for(int c = 0; c < coarse_size; ++c) {
                coarse_to_fine[c].reserve(subdiv_factor);
                for(int i = 0; i < subdiv_factor; ++i) {
                    int fine_id = c * subdiv_factor + i;
                    if(fine_id >= fine_size) {
                        std::cerr << "fine_id overflow: " << fine_id << " >= " << fine_size << std::endl;
                        return false;
                    }
                    coarse_to_fine[c].push_back(fine_id);
                    fine_to_coarse[fine_id] = c;
                }
            }
        } catch(const std::exception& e) {
            std::cerr << "Memory allocation failed: " << e.what() << std::endl;
            return false;
        }

        std::cout << "Mapping initialized: " << coarse_size << " coarse -> "
                  << fine_size << " fine (factor=" << subdiv_factor << ")" << std::endl;
        return true;
    }

    void update_fine_volumes_from_elements(const std::vector<FineElement>& fine_elements) {
        if(fine_elements.size() != fine_size) return;
        for(int i = 0; i < fine_size; ++i) {
            fine_volumes[i] = fine_elements[i].volume;
        }
    }

    bool isValid() const {
        return coarse_size > 0 && fine_size > 0 &&
               coarse_to_fine.size() == coarse_size &&
               fine_to_coarse.size() == fine_size &&
               fine_volumes.size() == fine_size;
    }
};



void GPU_OCC(std::vector<double>& Density, const double TargetRate, const std::vector<double>& dc,
             const std::vector<int>& isPicked, const std::vector<double>& volumes);

void MultiRes_TopOpt_Pipeline(
        const std::string& input_model_file,
        const std::string& bdy_file,
        const std::vector<int>& protect_cell_ids,
        const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
        double target_rate,
        double max_iter,
        double Yangs,
        double Poissons,
        int subdivision_level,
        std::vector<double>& density,
        std::atomic_bool& is_need_update_model);

void ComputeFineElementEnergies(
        //const std::vector<Eigen::MatrixXd>& coarse_local_matrix,
        const Eigen::MatrixXd& coarse_UU,
        const std::vector<CBSplineSolid>& coarse_patch,
        const double Youngs, const double Poissons,
        MultiResMapping& mapping,
        std::vector<FineElement>& fine_elements);

void ComputeSensitivitiesParallel(
        std::vector<FineElement>& fine_elements,
        const std::vector<double>& density);