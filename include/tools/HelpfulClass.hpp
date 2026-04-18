#pragma once
#include <cstdint>
#include "tools/MyPath.hpp"
#include "tools/ConfigurationSpace.hpp"
#include "core/UsefulMacros.hpp"
#include "tools/Environment.hpp"

// #define GAP 0.05
// #define GAPDECEN 0.1

template <typename P>
class BaseCollisionChecker {
    public:
        virtual bool isCollide(const P& point_) {
            (void)point_;
            return false;
        }
        virtual bool isCollide2P(const P& point1_, const P& point2_) {
            (void)point1_;
            (void)point2_;
            return false;
        }

        std::vector<std::pair<double,double>> getBounds() {return bounds;}
        void setBounds(std::vector<std::pair<double,double>> bounds_) {bounds = bounds_;}

    private:
        std::vector<std::pair<double,double>> bounds;
};

class Point2DCollisionChecker : public BaseCollisionChecker<Eigen::VectorXd>{
    public:
        Point2DCollisionChecker(const path_planning::Environment2D& enviroument_);

        bool isCollide(const Eigen::VectorXd& point_) override;
        bool isCollide2P(const Eigen::VectorXd& point1_, const Eigen::VectorXd& point2_) override;

        bool thisObstacle(int ob_idx_, const Eigen::VectorXd& point_);

    private:
        const path_planning::Environment2D& env;
};

class Point2DCollisionCheckerGrid : public BaseCollisionChecker<Eigen::VectorXd>{
    public:
        Point2DCollisionCheckerGrid(const path_planning::GridCSpace2D_T<int8_t>& map_);

        bool isCollide(const Eigen::VectorXd& point_) override;
        bool isCollide2P(const Eigen::VectorXd& point1_, const Eigen::VectorXd& point2_) override;

    private:
        const path_planning::GridCSpace2D_T<int8_t>& map;
};

class MultiAgentPoint2DCollisionCheckerGrid{
    public:
        MultiAgentPoint2DCollisionCheckerGrid(const std::vector<path_planning::GridCSpace2D_T<int8_t>>& maps_)
        : maps(maps_){}

        bool isCollide(int agent_idx_, const Eigen::VectorXd& point_);
        bool isCollide2P(int agent_idx_, const Eigen::VectorXd& point1_, const Eigen::VectorXd& point2_);

    private:
        const std::vector<path_planning::GridCSpace2D_T<int8_t>>& maps;
};

template <typename T>
class QuickSort{
    public:
        void run(std::vector<int>& arr, const std::unordered_map<int, T>& map, int low, int high){
            if (low < high) {
                // pi is the partition return index of pivot
                int pi = partition(arr, map, low, high);

                // recursion calls for smaller elements
                // and greater or equals elements
                run(arr, map, low, pi - 1);
                run(arr, map, pi + 1, high);
            }
        }

    private:
        int partition(std::vector<int>& arr, const std::unordered_map<int, T>& map, int low, int high){
            // choose the pivot
            double pivot = map.at(arr[high]);
        
            // undex of smaller element and indicates 
            // the right position of pivot found so far
            int i = low - 1;

            // Traverse arr[low..high] and move all smaller
            // elements on left side. Elements from low to 
            // i are smaller after every iteration
            for (int j = low; j <= high - 1; j++) {
                if (map.at(arr[j]) < pivot) {
                    i++;
                    std::swap(arr[i], arr[j]);
                }
            }
            
            // move pivot after smaller elements and
            // return its position
            std::swap(arr[i + 1], arr[high]);  
            return i + 1;
        }
};
