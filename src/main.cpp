#include <iostream>
#include "core/MySamplingBasedPlanners.hpp"

int main(){
    // build environment
    path_planning::Environment2D env;
    std::vector<Eigen::Vector2d> vertices_ccw = {Eigen::Vector2d(2.0, 3.0), Eigen::Vector2d(8.0, 3.0), Eigen::Vector2d(8.0, 4.0), Eigen::Vector2d(2.0, 4.0)};
    path_planning::Obstacle2D obs1(vertices_ccw);
    vertices_ccw = {Eigen::Vector2d(2.0, 6.0), Eigen::Vector2d(8.0, 6.0), Eigen::Vector2d(8.0, 7.0), Eigen::Vector2d(2.0, 7.0)};
    path_planning::Obstacle2D obs2(vertices_ccw);

    env.obstacles.push_back(obs1);
    env.obstacles.push_back(obs2);

    // test collision checker
    Point2DCollisionChecker collision_checker(env);
    LOG(collision_checker.isCollide(Eigen::Vector2d(3, 3.5)));
}