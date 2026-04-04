#pragma once

#include <vector>
#include <Eigen/Core>

namespace path_planning {

/// @brief 2-dimensional path
struct Path2D {
    std::vector<Eigen::Vector2d> waypoints;

    /// @brief `true` if a solution was found, `false` otherwise
    bool valid = false;
};

/// @brief N-dimensional path
struct Path {
    std::vector<Eigen::VectorXd> waypoints;

    /// @brief `true` if a solution was found, `false` otherwise
    bool valid = false;

    /// @brief Project the waypoints to 2D workspace
    /// @return 2D waypoints using the first two dimensions of the full state (x, y)
    std::vector<Eigen::Vector2d> getWaypoints2D() const {
        std::vector<Eigen::Vector2d> waypoints_2d;
        for (const Eigen::VectorXd& wp : waypoints) {
            waypoints_2d.push_back(Eigen::Vector2d(wp(0), wp(1)));
        }
        return waypoints_2d;
    }

    /// @brief Project the waypoints to 3D workspace
    /// @return 3D waypoints using the first three dimensions of the full state (x, y, z)
    std::vector<Eigen::Vector3d> getWaypoints3D() const {
        std::vector<Eigen::Vector3d> waypoints_3d;
        for (const Eigen::VectorXd& wp : waypoints) {
            waypoints_3d.push_back(Eigen::Vector3d(wp(0), wp(1), wp(2)));
        }
        return waypoints_3d;
    }
};
}
