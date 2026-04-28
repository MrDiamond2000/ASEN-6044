#include <iostream>
#include "core/MySamplingBasedPlanners.hpp"
#include "PF/ParticleFilter.hpp"
#include "GM/GaussianMixtureFilter.hpp"
#include <fstream>
#include <yaml-cpp/yaml.h>
#include "tools/Logging.hpp"
#include <chrono>

//     // test collision checker
//     Point2DCollisionChecker collision_checker(env);
//     // LOG("Is (3, 3.5) in collision: " << collision_checker.isCollide(Eigen::Vector2d(3, 3.5)));

//     // test path planner
//     MyRRT rrt(rrtBias, rrtIteration, rrtStepSize);
//     rrt.set_environment(env);
//     //path_planning::Path2D path_output = rrt.plan(Eigen::Vector2d(1,1), Eigen::Vector2d(5,5));

//     //LOG("Path valid: " << path_output.valid);
//     //if (path_output.valid){
//     //    for (size_t i = 0; i < path_output.waypoints.size(); i++){
//     //        LOG(path_output.waypoints[i].transpose());
//     //    }
//     //}

int main(){

    // Create environment and collision checkers
    path_planning::Environment2D env;
    env.deserialize("input/problem.yaml");

    Point2DCollisionChecker collision_checker(env);
    BaseCollisionChecker<Eigen::VectorXd>& lineOfSightChecker = collision_checker;

    double maxX = env.x_max;
    double maxY = env.y_max;

    // Load config files
    YAML::Node pfConfig = YAML::LoadFile("src/PF/ParticleFilterParams.yaml");
    YAML::Node gmConfig = YAML::LoadFile("src/GM/GaussianMixtureParams.yaml");

    // Define shared parameters
    double vehicleInitialX = pfConfig["vehicle"]["vehicleInitialX"].as<double>();
    double vehicleInitialY = pfConfig["vehicle"]["vehicleInitialY"].as<double>();
    double vehicleInitialHeading = pfConfig["vehicle"]["vehicleInitialHeading"].as<double>();

    double fovFractionOfPi = pfConfig["sensor"]["fovFractionOfPi"].as<double>();
    double range = pfConfig["sensor"]["range"].as<double>();
    int maxIterations = pfConfig["misc"]["maxIterations"].as<int>();

    double fov = M_PI/fovFractionOfPi;
    double fovCosine = cos(fov/2.0);

    // Generate target trajectory in advance
    LOG("Generating shared target trajectory...");

    ParticleFilter pfTemp(pfConfig["Filter"]["Ns"].as<int>(), maxX, maxY, pfConfig["Filter"]["gridSizeX"].as<int>(), pfConfig["Filter"]["gridSizeY"].as<int>());

    // Create initial target location
    Eigen::Vector2d target = pfTemp.initializeParticles(collision_checker);
    std::vector<Eigen::Vector2d> targetTrajectory;
    targetTrajectory.push_back(target);

    // Propagate the target trajectory through the max number of iterations
    double stepSizeTarget = pfConfig["Filter"]["stepSize"].as<double>();
    for(int i = 0; i < maxIterations; i++) {
        pfTemp.propagate(target, collision_checker, stepSizeTarget);
        targetTrajectory.push_back(target);
    }

    // Run the particle filter
    LOG("RUNNING PARTICLE FILTER");

    // Create particle filter object and pass in the number of particles and environment size
    ParticleFilter pf(pfConfig["Filter"]["Ns"].as<int>(), maxX, maxY, pfConfig["Filter"]["gridSizeX"].as<int>(), pfConfig["Filter"]["gridSizeY"].as<int>());

    // Generate initial particles
    pf.setMaxPropogationAttempts(pfConfig["Filter"]["maxPropogationAttempts"].as<int>());
    pf.setMaxObservedScore(pfConfig["Filter"]["maxObservedScore"].as<double>());
    pf.setObservedDecayRate(pfConfig["Filter"]["observedDecayRate"].as<double>());
    pf.initializeParticles(collision_checker); // ignore target

    Eigen::Vector2d pfVehicle(vehicleInitialX, vehicleInitialY);
    double pfHeading = vehicleInitialHeading;

    // Initialize path planner
    MyRRT pfRRT(pfConfig["pathPlanner"]["bias"].as<double>(), pfConfig["pathPlanner"]["iteration"].as<double>(), pfConfig["pathPlanner"]["stepSize"].as<double>());
    pfRRT.set_environment(env);
    path_planning::Path2D pfPath;
    size_t pfPathCounter = 0;

    // Open output csv file
    std::ofstream pfFile(pfConfig["misc"]["outputFile"].as<std::string>());

    int pfFinalIteration = -1;
    double pfTime = 0.0;

    // Run particle filter iterations
    for(int i = 0; i < maxIterations; i++) {

        // Find current target, estimate, and density estimate
        Eigen::Vector2d currTarget = targetTrajectory[i];
        Eigen::Vector2d estimate = pf.estimate();
        auto density = pf.estimate_density(collision_checker);

        // Calculate the vehicle heading
        Eigen::Vector2d direction;
        if(pfPath.valid && pfPathCounter+2 < pfPath.waypoints.size()) {
            direction = pfPath.waypoints[pfPathCounter+2] - pfVehicle;
        } else {
            direction = estimate - pfVehicle;
        }

        if(direction.norm() > 1e-6) {
            direction.normalize();
            pfHeading = atan2(direction(1), direction(0));
        }

        Eigen::Vector2d headingVec(cos(pfHeading), sin(pfHeading));

        // Output the vehicle, estimate, and target positions to the csv file
        pfFile << pfVehicle(0) << "," << pfVehicle(1) << "," << headingVec(0) << "," << headingVec(1) << "," << estimate(0) << "," << estimate(1) << "," << currTarget(0) << "," << currTarget(1) << "," << density.second(0) << "," << density.second(1);

        // Output the particle positions to the csv file
        for(auto& p : pf.particles) {
            pfFile << "," << p.position(0) << "," << p.position(1);
        }
        pfFile << "\n";

        // Run the particle filter and record runtime
        auto time1 = std::chrono::high_resolution_clock::now();
        pf.step(collision_checker, lineOfSightChecker, pfVehicle, headingVec, fovCosine, range, pfConfig["Filter"]["stepSize"].as<double>(), pfConfig["Filter"]["resampleThreshold"].as<double>());
        auto time2 = std::chrono::high_resolution_clock::now();

        pfTime += std::chrono::duration<double>(time2 - time1).count();

        // Check if the target has been detected
        if(pfFinalIteration == -1 && pf.detectTarget(currTarget, pfVehicle, headingVec, lineOfSightChecker, fovCosine, range)) {
            LOG("[PF] FOUND at iteration " << i);
            pfFinalIteration = i;
            break;
        }

        // Plan a path towards the filter estimate
        if (!pfPath.valid || pfPathCounter > 20 || pfPathCounter+2 >= pfPath.waypoints.size()) {

            Eigen::Vector2d goal = estimate;

            // Get estimate using density grid method if above threshold
            if(density.first > pfConfig["Filter"]["densityThreshold"].as<double>()) {
                LOG("[PF] Using density estimate | Density = " << density.first);
                goal = density.second;
            }
            else {
                LOG("[PF] Using weighting and exploration estimate | Density = " << density.first);
            }

            // Plan the path
            pfPath = pfRRT.plan(pfVehicle, goal);
            if(!pfPath.valid) {
                LOG("[PF] No path found");
                break;
            }

            pfPathCounter = 0;
        } 
        else {
            pfPathCounter++;
        }

        // Move one step along the path
        pfVehicle = pfPath.waypoints[pfPathCounter+1];
    }

    pfFile.close();

    // Run the gaussian mixture filter
    LOG("RUNNING GAUSSIAN MIXTURE FILTER");

    // Create gaussian mixture filter object and pass in the number of particles and environment size
    GaussianMixtureFilter gm(gmConfig["Filter"]["M"].as<int>(), maxX, maxY, gmConfig["Filter"]["gridSizeX"].as<int>(), gmConfig["Filter"]["gridSizeY"].as<int>());

    // Generate initial components
    gm.setMaxPropogationAttempts(gmConfig["Filter"]["maxPropogationAttempts"].as<int>());
    gm.setMaxObservedScore(gmConfig["Filter"]["maxObservedScore"].as<double>());
    gm.setObservedDecayRate(gmConfig["Filter"]["observedDecayRate"].as<double>());
    gm.initializeObjects(collision_checker);

    Eigen::Vector2d gmVehicle(vehicleInitialX, vehicleInitialY);
    double gmHeading = vehicleInitialHeading;

    // Initialize path planner
    MyRRT gmRRT(gmConfig["pathPlanner"]["bias"].as<double>(), gmConfig["pathPlanner"]["iteration"].as<double>(), gmConfig["pathPlanner"]["stepSize"].as<double>());
    gmRRT.set_environment(env);
    path_planning::Path2D gmPath;
    size_t gmPathCounter = 0;

    // Open output csv file
    std::ofstream gmFile(gmConfig["misc"]["outputFile"].as<std::string>());

    int gmFinalIteration = -1;
    double gmTime = 0.0;

    // Run gaussian mixture filter iterations
    for(int i = 0; i < maxIterations; i++){

        // Find current target, estimate, and density estimate
        Eigen::Vector2d currTarget = targetTrajectory[i];
        Eigen::Vector2d estimate = gm.estimate(collision_checker);
        auto density = gm.estimate_density(collision_checker);

        // Calculate the vehicle heading
        Eigen::Vector2d direction;
        if(gmPath.valid && gmPathCounter+2 < gmPath.waypoints.size()) {
            direction = gmPath.waypoints[gmPathCounter+2] - gmVehicle;
        } else {
            direction = estimate - gmVehicle;
        }

        if(direction.norm() > 1e-6){
            direction.normalize();
            gmHeading = atan2(direction(1), direction(0));
        }

        Eigen::Vector2d headingVec(cos(gmHeading), sin(gmHeading));

        // Output the vehicle, estimate, and target positions to the csv file
        gmFile << gmVehicle(0) << "," << gmVehicle(1) << "," << headingVec(0) << "," << headingVec(1) << "," << estimate(0) << "," << estimate(1) << "," << currTarget(0) << "," << currTarget(1) << "," << density.second(0) << "," << density.second(1);

        // Output the component positions, covariances, and weights to the csv file
        for(auto& g : gm.components) {
            gmFile << "," << g.position(0) << "," << g.position(1) << "," << g.covariance(0,0) << "," << g.covariance(0,1) << "," << g.covariance(1,0) << "," << g.covariance(1,1) << "," << g.weight;
        }
        gmFile << "\n";

        // Run the gaussian mixture filter and record runtime
        Eigen::Matrix2d randWalkCov = 0.1*gmConfig["Filter"]["stepSize"].as<double>()*gmConfig["Filter"]["stepSize"].as<double>()*Eigen::Matrix2d::Identity();

        auto time1 = std::chrono::high_resolution_clock::now();
        gm.step(collision_checker, lineOfSightChecker, gmVehicle, headingVec, fov, fovCosine, range, randWalkCov, gmConfig["Filter"]["stepSize"].as<double>());
        auto time2 = std::chrono::high_resolution_clock::now();

        gmTime += std::chrono::duration<double>(time2 - time1).count();

        // Check if the target has been detected
        if(gmFinalIteration == -1 && gm.detectTarget(currTarget, gmVehicle, headingVec, lineOfSightChecker, fovCosine, range)) {
            LOG("[GM] FOUND at iteration " << i);
            gmFinalIteration = i;
            break;
        }

        // Plan a path towards the filter estimate
        if (!gmPath.valid || gmPathCounter > 20 || gmPathCounter+2 >= gmPath.waypoints.size()) {

            Eigen::Vector2d goal = estimate;

            // Get estimate using density grid method if above threshold
            if(density.first > gmConfig["Filter"]["densityThreshold"].as<double>()) {
                LOG("[GM] Using density estimate | Density = " << density.first);
                goal = density.second;
            }
            else {
                LOG("[GM] Using weighting and exploration estimate | Density = " << density.first);
            }

            // Plan the path
            gmPath = gmRRT.plan(gmVehicle, goal);
            if(!gmPath.valid){
                LOG("[GM] No path found");
                break;
            }

            gmPathCounter = 0;
        } 
        else {
            gmPathCounter++;
        }

        // Move one step along the path
        gmVehicle = gmPath.waypoints[gmPathCounter+1];
    }

    gmFile.close();

    // Print results
    std::cout << "\nFINAL RESULTS\n";

    std::cout << "Particle Filter:\n";
    std::cout << "Iteration: " << pfFinalIteration << "\n";
    std::cout << "Time: " << pfTime << " sec\n";

    std::cout << "\nGaussian Mixture:\n";
    std::cout << "Iteration: " << gmFinalIteration << "\n";
    std::cout << "Time: " << gmTime << " sec\n";

    std::cout << "\nSpeed ratio (GM / PF): " << (gmTime/pfTime) << "\n";

    return 0;
}