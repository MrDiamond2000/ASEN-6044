#pragma once

// This includes all of the necessary header files in the toolbox
// #include "HelpfulClass.h"
#include <time.h>
#include <cmath>
#include "tools/Graph.hpp"
#include "tools/MyPath.hpp"
#include "tools/HelpfulClass.hpp"
#include <Eigen/Core>
#include <map>

class MyGenericRRT {
    public:
        MyGenericRRT(double bias_, int iteration_, double step_size_) : bias(bias_), iteration(iteration_), step_size(step_size_) {}

        std::shared_ptr<path_planning::Graph<double>> getGraphPtr() { return graphPtr; }

        std::map<path_planning::Node , Eigen::VectorXd> getNodes() { return nodes; }

        path_planning::Path planND(Eigen::VectorXd init_, Eigen::VectorXd goal_, BaseCollisionChecker<Eigen::VectorXd>& collision_checker_);

        // path_planning::Path planNDDecen(Eigen::VectorXd init_, Eigen::VectorXd goal_, MultiAgentDisk2DCollisionCheckerDecen& collision_checker_);

    private:
        std::shared_ptr<path_planning::Graph<double>> graphPtr = std::make_shared<path_planning::Graph<double>>();
        std::map<path_planning::Node , Eigen::VectorXd> nodes;
        Eigen::VectorXd generatePoint(const std::vector<std::pair<double, double>>& bounds);
        path_planning::Node  closestPoint(const Eigen::VectorXd& point);
        Eigen::VectorXd extendRRT(const Eigen::VectorXd& point, BaseCollisionChecker<Eigen::VectorXd>& collision_checker_);
        // Eigen::VectorXd extendRRTDecen(const Eigen::VectorXd& point, MultiAgentDisk2DCollisionCheckerDecen& collision_checker_);
        bool checkDistance(Eigen::VectorXd direction, double requirement);
        double magnitude(Eigen::VectorXd vec);
        double bias;
        int iteration;
        double step_size;
};

class MyRRT : public MyGenericRRT {
    public:
        MyRRT(double bias_, int iteration_, double step_size_) : MyGenericRRT(bias_, iteration_, step_size_) {}
        void set_environment(path_planning::Environment2D environment_) {environment = environment_;}

        path_planning::Path2D plan(Eigen::Vector2d init_, Eigen::Vector2d goal_); 
        std::map<path_planning::Node, Eigen::Vector2d> getNodes2D();

    private:
        path_planning::Environment2D environment;
};
