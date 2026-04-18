#pragma once

#include <Eigen/Core>

#include "tools/Obstacle.hpp" 
#include "tools/MyPath.hpp" 
#include "tools/Serializer.hpp"

namespace path_planning { 

/// @brief 2D workspace with rectangular bounds and obstacles
struct Environment2D {
    double x_min = 0.0; 
    double x_max = 10.0; 
    double y_min = 0.0; 
    double y_max = 10.0; 
    std::vector<Obstacle2D> obstacles;

    /// @brief Read and setup problemenviroumet from yaml file
    /// @param yaml_file Path to yaml file containing problem description
    void deserialize(const std::string& yaml_file){
        Deserializer deserializer(yaml_file);
        if (!deserializer){
            return;
        }
        const YAML::Node& node = deserializer.get();
        x_min = node["x_min"].as<double>();
        x_max = node["x_max"].as<double>();
        y_min = node["y_min"].as<double>();
        y_max = node["y_max"].as<double>();
        obstacles.clear();
        for (const auto& obs_node : node["obstacles"]){
            std::vector<Eigen::Vector2d> vertices_ccw;
            for (const auto& vertex_node : obs_node["vertices_ccw"]){
                double x;
                double y;
                if (vertex_node.IsSequence()){
                    x = vertex_node[0].as<double>();
                    y = vertex_node[1].as<double>();
                } else {
                    x = vertex_node["x"].as<double>();
                    y = vertex_node["y"].as<double>();
                }
                vertices_ccw.push_back(Eigen::Vector2d(x, y));
            }
            obstacles.push_back(vertices_ccw);
        }
    };
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