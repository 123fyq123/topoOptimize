#include "Poster/source/nanoflann.h"

namespace gpmesh {
    using namespace nanoflann;
    struct Point {
        double x, y, z;
    };

    template<typename T>
    struct PointCloud {

        using coord_t = T;  //!< The type of each coordinate

        std::vector<Point> pts;

        // Must return the number of data points
        inline size_t kdtree_get_point_count() const { return pts.size(); }

        // Returns the dim'th component of the idx'th point in the class:
        // Since this is inlined and the "dim" argument is typically an immediate
        // value, the
        //  "if/else's" are actually solved at compile time.
        inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
            if (dim == 0)
                return pts[idx].x;
            else if (dim == 1)
                return pts[idx].y;
            else
                return pts[idx].z;
        }

        // Optional bounding-box computation: return false to default to a standard
        // bbox computation loop.
        //   Return true if the BBOX was already computed by the class and returned
        //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
        //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
        template<class BBOX>
        inline bool kdtree_get_bbox(BBOX & /* bb */) const {
            return false;
        }
    };

    typedef KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<double, PointCloud<double>>,
            PointCloud<double>,
            3 /* dim */
    > KDTree;

    inline void GetNearestIds(std::vector<Point> &p, std::vector<Point> &q, std::vector<int> &ans) {
        ans.resize(p.size());
        PointCloud<double> A, B;
        A.pts = p;
        B.pts = q;
        KDTree index(3, B, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
        index.buildIndex();
        int a_index = 0;
        for (const auto &a: A.pts) {
            size_t ret_index;
            double out_dist_sqr;
            nanoflann::KNNResultSet<double> resultSet(1);
            resultSet.init(&ret_index, &out_dist_sqr);
            index.findNeighbors(resultSet, &a.x);
            ans[a_index++] = ret_index;
        }
    }

}