#include <iostream>
#include "core/MySamplingBasedPlanners.hpp"
#include "PF/ParticleFilter.hpp"
#include <fstream>
#include <yaml-cpp/yaml.h>

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
    LOG("Is (3, 3.5) in collision: " << collision_checker.isCollide(Eigen::Vector2d(3, 3.5)));

    // load in parameters from the parameter yaml file
    YAML::Node config = YAML::LoadFile("src/PF/ParticleFilterParams.yaml");

    double maxX = config["environment"]["maxX"].as<double>();
    double maxY = config["environment"]["maxY"].as<double>();

    double rrtBias = config["pathPlanner"]["bias"].as<double>();
    double rrtIteration = config["pathPlanner"]["iteration"].as<double>();
    double rrtStepSize = config["pathPlanner"]["stepSize"].as<double>();

    int Ns = config["particleFilter"]["Ns"].as<int>();
    double stepSize = config["particleFilter"]["stepSize"].as<double>();
    double resampleThreshold = config["particleFilter"]["resampleThreshold"].as<double>();
    int gridSizeX = config["particleFilter"]["gridSizeX"].as<int>();
    int gridSizeY = config["particleFilter"]["gridSizeY"].as<int>();
    int maxPropogationAttempts = config["particleFilter"]["maxPropogationAttempts"].as<int>();
    double maxObservedScore = config["particleFilter"]["maxObservedScore"].as<double>();
    double observedDecayRate = config["particleFilter"]["observedDecayRate"].as<double>();

    double vehicleInitialX = config["vehicle"]["vehicleInitialX"].as<double>();
    double vehicleInitialY = config["vehicle"]["vehicleInitialY"].as<double>();
    double vehicleInitialHeading = config["vehicle"]["vehicleInitialHeading"].as<double>();

    double fovFractionOfPi = config["sensor"]["fovFractionOfPi"].as<double>();
    double range = config["sensor"]["range"].as<double>();

    int maxIterations = config["misc"]["maxIterations"].as<int>();
    std::string outputFile = config["misc"]["outputFile"].as<std::string>();

    // test path planner
    MyRRT rrt(rrtBias, rrtIteration, rrtStepSize);
    rrt.set_environment(env);
    //path_planning::Path2D path_output = rrt.plan(Eigen::Vector2d(1,1), Eigen::Vector2d(5,5));

    //LOG("Path valid: " << path_output.valid);
    //if (path_output.valid){
    //    for (size_t i = 0; i < path_output.waypoints.size(); i++){
    //        LOG(path_output.waypoints[i].transpose());
    //    }
    //}

    // test particle filter
    BaseCollisionChecker<Eigen::VectorXd>& lineOfSightChecker = collision_checker;

    // create particle filter object and pass in the number of particles and environment size
    ParticleFilter pf(Ns, maxX, maxY, gridSizeX, gridSizeY);

    // generate initial particles, target, and estimate
    pf.setMaxPropogationAttempts(maxPropogationAttempts);
    pf.setMaxObservedScore(maxObservedScore);
    pf.setObservedDecayRate(observedDecayRate);
    Eigen::Vector2d target = pf.initializeParticles(collision_checker);
    Eigen::Vector2d estimate = pf.estimate();
    Eigen::Vector2d direction;
    Eigen::Vector2d headingVector;

    // open output csv file
    std::ofstream file(outputFile, std::ios::out);

    // initialize a search vehicle
    Eigen::Vector2d vehicle(vehicleInitialX, vehicleInitialY);

    // define sensor parameters
    double fov = M_PI/fovFractionOfPi;
    double fovCosine = cos(fov/2.0);

    // run particle filter iterations and save data to a csv
    bool found = false;
    double currentHeading = vehicleInitialHeading;
    for (int i = 0; i < maxIterations; i++) {

        // calculate the vehicle heading
        direction = estimate - vehicle;

        if (direction.norm() > 1e-6) {
            direction.normalize();
            currentHeading = atan2(direction(1), direction(0));
        }

        headingVector = {cos(currentHeading), sin(currentHeading)};

        // output the vehicle, estimate, and target positions to the csv file
        file << vehicle(0) << "," << vehicle(1) << "," << headingVector(0) << ","  << headingVector(1) << "," << estimate(0) << "," << estimate(1) << "," << target(0) << "," << target(1);
        
        // output the particle positions to the csv file
        for (auto& p : pf.particles) {
            file << "," << p.position(0) << "," << p.position(1);
        }
        file << "\n";

        // run the particle filter
        pf.step(collision_checker, lineOfSightChecker, vehicle, headingVector, fovCosine, range, stepSize, resampleThreshold);
        estimate = pf.estimate();

        // check if the target has been detected
        if (pf.detectTarget(target, vehicle, headingVector, lineOfSightChecker, fovCosine, range)) {
            std::cout << "Target FOUND at iteration " << i << std::endl;
            std::cout << "True target:      " << target.transpose() << std::endl;
            std::cout << "Estimate:         " << pf.estimate().transpose() << std::endl;
            std::cout << "Vehicle Position: " << vehicle.transpose() << std::endl;
            found = true;
            break;
        }

        // plan a path towards the filter estimate
        path_planning::Path2D path = rrt.plan(vehicle, estimate);

        if (path.valid && path.waypoints.size() > 1) {
            // move one step along the path, but normalize to the step size (velocity) of the vehicle to ensure consistency
            Eigen::Vector2d next = path.waypoints[1];
            Eigen::Vector2d delta = next - vehicle;

            if (delta.norm() > stepSize) {
                delta.normalize();
                delta *= stepSize;
            }

            vehicle += delta;
        }

        // propogate the target position for the next iteration
        pf.propagate(target, collision_checker, stepSize);
    }

    // output the final particle filter estimate to the csv
    file << vehicle(0) << "," << vehicle(1) << "," << headingVector(0) << ","  << headingVector(1) << "," << estimate(0) << "," << estimate(1) << "," << target(0) << "," << target(1);

    // output the final particle positions to the csv
    for (auto& p : pf.particles) {
            file << "," << p.position(0) << "," << p.position(1);
        }
        file << "\n";

    // close the output csv file
    file.close();

    // if the filter did not locate the target, print it's final estimate
    if (!found) {
        std::cout << "Target NOT found.\n";
        std::cout << "Final estimate: " << estimate.transpose() << std::endl;
    }

    return 0;
}