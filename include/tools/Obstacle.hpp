#pragma once

#include <vector>
#include <Eigen/Core>

// #include "tools/Serializer.h"

namespace path_planning {

class Polygon {
    public:
        /// @brief Used for deserialization
        Polygon() = default;
        
        /// @brief Construct from a set of vertices ASSUMED to be in counter clock wise order
        /// @param vertices_ccw Ordered counter clock wise order
        Polygon(std::vector<Eigen::Vector2d>& vertices_ccw){
            m_vertices_ccw = vertices_ccw;
        }
        
        /// @brief Access the vertices in memory of the obstacle (counter clock wise coordinates)
        /// @return Reference to stored CCW vertices
        std::vector<Eigen::Vector2d>& verticesCCW() {
            return m_vertices_ccw;
        } 

        /// @brief Access the vertices in memory of the obstacle (counter clock wise coordinates)
        /// @return Reference to stored CCW vertices
        const std::vector<Eigen::Vector2d>& verticesCCW() const{
            return m_vertices_ccw;
        }; 

    private:
        std::vector<Eigen::Vector2d> m_vertices_ccw;
};

using Obstacle2D = Polygon;

}
