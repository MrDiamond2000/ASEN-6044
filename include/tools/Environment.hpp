#pragma once

#include <Eigen/Core>

#include "tools/Obstacle.hpp" 
#include "tools/MyPath.hpp" 

namespace path_planning { 

/// @brief 2D workspace with rectangular bounds and obstacles
struct Environment2D {
    double x_min = 0.0; 
    double x_max = 10.0; 
    double y_min = 0.0; 
    double y_max = 10.0; 
    std::vector<Obstacle2D> obstacles;
};

/// @brief Environment with initial state and goal state
struct Problem2D : Environment2D {
    /// @brief Mobile robot: location of reference point on robot
    /// Manipulator: end effector location
    Eigen::Vector2d q_init;
    /// @brief Mobile robot: location of reference point on robot
    /// Manipulator: end effector location
    Eigen::Vector2d q_goal;
};

}