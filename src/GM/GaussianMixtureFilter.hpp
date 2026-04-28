#include <time.h>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <map>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include "tools/HelpfulClass.hpp"

class GaussianMixtureFilter {
    public:
        GaussianMixtureFilter(int M_, int maxX_, int maxY_, int gridSizeX_, int gridSizeY_)
        :   M(M_), 
            gen(std::random_device{}()), 
            angleDist(0, 2*M_PI), 
            gridSizeX(gridSizeX_), 
            gridSizeY(gridSizeY_), 
            maxPropogationAttempts(20), 
            maxObservedScore(10), 
            observedDecayRate(0.995), 
            maxX(maxX_), maxY(maxY_), 
            distX(0.0, maxX), 
            distY(0.0, maxY) {
                observedGrid.resize(gridSizeX, std::vector<double>(gridSizeY, 0.0));
                densityGrid.resize(gridSizeX, std::vector<double>(gridSizeY, 0.0));
                surroundingDensityGrid.resize(gridSizeX, std::vector<double>(gridSizeY, 0.0));
        }

        // Define GM structure
        struct GaussianObject {
            Eigen::Vector2d position;
            Eigen::Matrix2d covariance;
            double weight;

            GaussianObject() = default;
            GaussianObject(const GaussianObject&) = default;
            GaussianObject& operator=(const GaussianObject&) = default;
        };

        // Define the list of components
        std::vector<GaussianObject> components;
        std::vector<GaussianObject> updatedComponents;

        // Randomely generate initial component positions and add corresponding components, then return a randomely generated target
        Eigen::Vector2d initializeObjects(Point2DCollisionChecker& checker) {
            for (int i = 0; i < M; i++) {
                Eigen::Vector2d pos;
                do { pos = {distX(gen), distY(gen)}; } while (checker.isCollide(pos));
                addComponent(pos, 1.0/M);
            }

            Eigen::Vector2d target;
            do { target = {distX(gen), distY(gen)}; } while (checker.isCollide(target));
            return target;
        }

        // Add a component to the filter
        void addComponent(Eigen::Vector2d& position, double weight) {
            Eigen::Matrix2d covariance = Eigen::Matrix2d::Identity()*0.1; // Not totally sure what a good initial covariance is
            components.push_back({position, covariance, weight});
        }

        // Add a component to the filter
        // void addComponent(Eigen::Vector2d& position, double weight, Eigen::Matrix2d covariance) {
        //     Eigen::Matrix2d covariance = Eigen::Matrix2d::Identity()*0.1; // Not totally sure what a good initial covariance is
        //     components.push_back({position, covariance, weight});
        // }

        // Setters
        void setMaxPropogationAttempts(double val) {
            maxPropogationAttempts = val;
        }

        void setMaxObservedScore(double val) {
            maxObservedScore = val;
        }

        void setObservedDecayRate(double val) {
            observedDecayRate = val;
        }

        
        // Updates the "observed" score of each point in the grid based on the current observer position
        // Decays the observed score of each grid cell slightly, and then increases the score of the currently observed grid cells
        void updateObservedGrid(const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, BaseCollisionChecker<Eigen::VectorXd>& checker, double fovCosine, double range)
        {
            // Decay the observed score
            for (int i = 0; i < gridSizeX; i++) {
                for (int j = 0; j < gridSizeY; j++) {
                    observedGrid[i][j] *= observedDecayRate;
                }
            }

            // Increase observed score for currently observable grid cells
            for (int i = 0; i < gridSizeX; i++) {
                for (int j = 0; j < gridSizeY; j++) {
                    Eigen::Vector2d cellPosition((i + 0.5)*maxX/static_cast<double>(gridSizeX), (j + 0.5)*maxY/static_cast<double>(gridSizeY));

                    if (observed(cellPosition, observerPosition, observerHeading, fovCosine, range) && !checker.isCollide2P(Eigen::VectorXd(observerPosition), Eigen::VectorXd(cellPosition))) {
                        // Do not increase observed score above 5 to avoid blowing up the score
                        observedGrid[i][j] = std::min(observedGrid[i][j] + 1.0, maxObservedScore);
                    }
                }
            }
        }
        
        // Get the observed score of a given position (from nearest grid point)
        double getObservedScore(const Eigen::Vector2d& pos) {
            int ix = std::min(std::max(int((pos(0)/maxX)*gridSizeX), 0), gridSizeX - 1);
            int iy = std::min(std::max(int((pos(1)/maxY)*gridSizeY), 0), gridSizeY - 1);
            return observedGrid[ix][iy];
        }

        // Check whether the component (or target) is observed in a defined cone infront of the observer
        bool observed(const Eigen::Vector2d& position, const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, double fovCosine, double range) {

            double dx = position(0) - observerPosition(0);
            double dy = position(1) - observerPosition(1);

            if (dx > range || dx < -range || dy > range || dy < -range) return false;

            double distanceSquared = dx*dx + dy*dy;

            if (distanceSquared > range*range) return false;

            double dot = dx*observerHeading(0) + dy*observerHeading(1);
            return dot > sqrt(distanceSquared)*fovCosine;
        }

        // Set the likelihood of any particles in the conical observable area low, otherwise one
        double likelihood(const GaussianObject& g, const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, BaseCollisionChecker<Eigen::VectorXd>& checker, double fovCosine, double range) {
            if (observed(g.position, observerPosition, observerHeading, fovCosine, range) && !checker.isCollide2P(Eigen::VectorXd(observerPosition), Eigen::VectorXd(g.position))) {

                double mahalanobisDistance = (g.position - observerPosition).transpose()*g.covariance.inverse()*(g.position - observerPosition);

                // Apply negative information
                double probabilityDetected = std::exp(-0.5*mahalanobisDistance);
                return std::max(1.0 - 0.7*probabilityDetected, 0.05); // minimum likelihood at 0.05 if observed, maximum 0.3 seems to work
            }

            return 1.0;
        }

        // Check whether the target is observed in the conical observable area
        bool detectTarget(const Eigen::Vector2d& truthTargetPosition, const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, BaseCollisionChecker<Eigen::VectorXd>& checker, double fovCosine, double range) {
            if (observed(truthTargetPosition, observerPosition, observerHeading, fovCosine, range) && !checker.isCollide2P(Eigen::VectorXd(observerPosition), Eigen::VectorXd(truthTargetPosition))) {
                return true;
            }
            else {
                return false;
            }
        }

        // Normalize weights, and if all weights are zero then reinitialize weights as uniform
        void normalize() {
            double sum = 0.0;
            
            for (auto& p : components) sum += p.weight;

            if (sum == 0.0) {
                double uniformWeight = 1.0/components.size();
                for (auto& p : components) p.weight = uniformWeight;
                return;
            }

            for (auto& p : components) p.weight /= sum;
        }

        // Execute a random walk of distance 0.1 in a random direction (provided there are no obstacles)
        // Serves to update the state of the true target
        void propagate(Eigen::Vector2d& position, Point2DCollisionChecker& checker, double stepSize) {
            for (int i = 0; i < maxPropogationAttempts; i++) {
                double movementAngle = angleDist(gen);
                Eigen::Vector2d step(stepSize*cos(movementAngle), stepSize*sin(movementAngle));
                Eigen::Vector2d newPos = position + step;

                if ((!checker.isCollide(newPos)) && (newPos(0) > 0 && newPos(0) < maxX) && (newPos(1) > 0 && newPos(1) < maxY)) {
                    position = newPos;
                    return;
                }
            }
        }

        // Finds the position of the combined highest weight and most unobserved component
        // Eigen::Vector2d estimate(Point2DCollisionChecker& checker) {
        //     double maxScore = 0.0;
        //     Eigen::Vector2d bestComponentPosition;

        //     for (auto& p : components) {
        //         double observedScore = getObservedScore(p.position);
        //         double score = p.weight*exp(-observedScore) + 0.5*exp(-observedScore);

        //         if (score > maxScore) {
        //             maxScore = score;
        //             bestComponentPosition = p.position;
        //         }
        //     }

        //     // Compared to the particle filter, produced a very jumpy estimate, so apply smoothing
        //     const double alpha = 0.15; // smoothing factor

        //     if (!hasSmoothedEstimate) {
        //         smoothedEstimate = bestComponentPosition;
        //         hasSmoothedEstimate = true;
        //     }
        //     else {
        //         Eigen::Vector2d candidate = (1.0 - alpha)*smoothedEstimate + alpha*bestComponentPosition;

        //         // If the candidate is valid, use it
        //         if (!checker.isCollide(candidate)) {
        //             smoothedEstimate = candidate;
        //         }
        //         else {
        //             // If not valid, move toward bestComponent instead
        //             smoothedEstimate = projectToFreeSpace(bestComponentPosition, checker);
        //         }
        //     }

        //     return smoothedEstimate;
        // }

        // // Finds the highest density grid cell (in terms of number of components), as long as it is not highly observed
        // std::pair<double, Eigen::Vector2d> estimate_density(Point2DCollisionChecker& checker) {
        //     const int radius = 3;

        //     // Reset grids
        //     for (int i = 0; i < gridSizeX; i++) {
        //         for (int j = 0; j < gridSizeY; j++) {
        //             densityGrid[i][j] = 0.0;
        //             surroundingDensityGrid[i][j] = 0.0;
        //         }
        //     }

        //     // Add component weights to density grid
        //     for (auto& p : components) {
        //         int ix = std::min(std::max(int((p.position(0)/maxX)*gridSizeX), 0), gridSizeX - 1);
        //         int iy = std::min(std::max(int((p.position(1)/maxY)*gridSizeY), 0), gridSizeY - 1);
        //         densityGrid[ix][iy] += p.weight;
        //     }

        //     // Calculate the number of cells in the circular area
        //     int cellCount = 0;
        //     for (int di = -radius; di <= radius; di++) {
        //         for (int dj = -radius; dj <= radius; dj++) {
        //             if (di*di + dj*dj <= radius*radius) {
        //                 cellCount++;
        //             }
        //         }
        //     }

        //     double cellArea = (maxX/gridSizeX)*(maxY/gridSizeY);
        //     double area = cellCount*cellArea;

        //     // Circular smoothing of density
        //     for (int i = 0; i < gridSizeX; i++) {
        //         for (int j = 0; j < gridSizeY; j++) {

        //             double weightedSum = 0.0;

        //             for (int di = -radius; di <= radius; di++) {
        //                 for (int dj = -radius; dj <= radius; dj++) {

        //                     // Circular mask
        //                     if (di*di + dj*dj > radius*radius) continue;
        //                     int ni = i + di;
        //                     int nj = j + dj;

        //                     if (ni >= 0 && ni < gridSizeX && nj >= 0 && nj < gridSizeY) {
        //                         double obsPenalty = std::exp(-observedGrid[ni][nj]);
        //                         weightedSum += densityGrid[ni][nj] * obsPenalty;
        //                     }
        //                 }
        //             }

        //             surroundingDensityGrid[i][j] = weightedSum;
        //         }
        //     }

        //     // Find cell with highest density
        //     double maxDensity = 0.0;
        //     Eigen::Vector2d bestComponentPosition;

        //     for (int i = 0; i < gridSizeX; i++) {
        //         for (int j = 0; j < gridSizeY; j++) {
        //             if (surroundingDensityGrid[i][j] <= maxDensity) continue;

        //             Eigen::Vector2d pos((i + 0.5)*maxX/static_cast<double>(gridSizeX), (j + 0.5)*maxY/static_cast<double>(gridSizeY));
        //             if (checker.isCollide(pos)) continue;

        //             maxDensity = surroundingDensityGrid[i][j];
        //             bestComponentPosition = pos;
        //         }
        //     }

        //     // Convert density to weight per unit area
        //     double densityPerArea = maxDensity/area;
        //     return {densityPerArea, bestComponentPosition};
        // }

        // Generate estimate based on component sampling density
        // Essentially constructs the posterior using Monte Carlo sampling and chooses the point of the highest probability density as the estimate
        Eigen::Vector2d estimate(Point2DCollisionChecker& checker) {

            // Build distribution based on component weights
            std::vector<double> weights;
            weights.reserve(components.size());
            for (const auto& g : components) weights.push_back(g.weight);

            std::discrete_distribution<int> weightDistribution(weights.begin(), weights.end());

            // Sample from the gaussian mixture
            std::vector<Eigen::Vector2d> samples;
            samples.reserve(10*M);

            std::normal_distribution<double> normal(0.0, 1.0);

            for (int k = 0; k < M; k++) {
                int componentIndex = weightDistribution(gen);
                const auto& g = components[componentIndex];

                Eigen::Vector2d standardNormal;
                standardNormal << normal(gen), normal(gen);

                Eigen::Matrix2d choleskyCovFactor = g.covariance.llt().matrixL();
                Eigen::Vector2d sample = g.position + choleskyCovFactor*standardNormal;

                if (!checker.isCollide(sample)) {
                    samples.push_back(sample);
                }
            }

            if (samples.empty()) {
                return smoothedEstimate;
            }

            // Build density grid

            // Reset grids
            for (int i = 0; i < gridSizeX; i++) {
                for (int j = 0; j < gridSizeY; j++) {
                    densityGrid[i][j] = 0.0;
                    surroundingDensityGrid[i][j] = 0.0;
                }
            }

            // Add uniform sample weights to build a histogram approximation of the pdf
            for (auto& s : samples) {
                int ix = std::min(std::max(int((s(0)/maxX)*gridSizeX), 0), gridSizeX - 1);
                int iy = std::min(std::max(int((s(1)/maxY)*gridSizeY), 0), gridSizeY - 1);
                densityGrid[ix][iy] += 1.0;
            }

            const int radius = 3;

            // Count circular cells
            // int cellCount = 0;
            // for (int di = -radius; di <= radius; di++) {
            //     for (int dj = -radius; dj <= radius; dj++) {
            //         if (di*di + dj*dj <= radius*radius) {
            //             cellCount++;
            //         }
            //     }
            // }

            //double cellArea = (maxX/gridSizeX)*(maxY/gridSizeY);
            //double area = cellCount * cellArea;

            // Circular density smoothing based on grid cells within the defined radius
            for (int i = 0; i < gridSizeX; i++) {
                for (int j = 0; j < gridSizeY; j++) {

                    double weightedSum = 0.0;

                    for (int di = -radius; di <= radius; di++) {
                        for (int dj = -radius; dj <= radius; dj++) {

                            if (di*di + dj*dj > radius*radius) continue;

                            int ni = i + di;
                            int nj = j + dj;

                            if (ni >= 0 && ni < gridSizeX && nj >= 0 && nj < gridSizeY) {
                                weightedSum += densityGrid[ni][nj];
                            }
                        }
                    }

                    surroundingDensityGrid[i][j] = weightedSum;
                }
            }

            // Scan through all cells to find the highest density cell
            double maxDensity = 0.0;
            Eigen::Vector2d highestDensityPosition = smoothedEstimate;

            for (int i = 0; i < gridSizeX; i++) {
                for (int j = 0; j < gridSizeY; j++) {

                    if (surroundingDensityGrid[i][j] <= maxDensity) continue;

                    Eigen::Vector2d pos((i + 0.5)*maxX/static_cast<double>(gridSizeX), (j + 0.5)*maxY/static_cast<double>(gridSizeY));

                    if (checker.isCollide(pos)) continue;

                    maxDensity = surroundingDensityGrid[i][j];
                    highestDensityPosition = pos;
                }
            }

            // Smooth the density estimate so it doesn't move locations as quickly
            const double alpha = 0.05;

            if (!hasSmoothedEstimate) {
                smoothedEstimate = highestDensityPosition;
                hasSmoothedEstimate = true;
            } 
            else {
                Eigen::Vector2d candidate = (1.0 - alpha)*smoothedEstimate + alpha*highestDensityPosition;

                if (!checker.isCollide(candidate)) {
                    smoothedEstimate = candidate;
                } 
                else {
                    smoothedEstimate = projectToFreeSpace(highestDensityPosition, checker);
                }
            }

            return smoothedEstimate;
        }

        void choppedGaussian(GaussianObject& g, std::pair<Eigen::Vector2d,Eigen::Vector2d> seg) {
            // Get the line segment in the form a^T x = b
            Eigen::Vector2d d = seg.second - seg.first;
            Eigen::Vector2d a(-d(1), d(0)); // normal vector
            double b = a.dot(seg.second);

            Eigen::Vector2d mu = g.position;

            if (a.dot(mu) > b) return;

            Eigen::Matrix2d Sigma = g.covariance;

            double m = a.dot(mu);
            double s2 = a.transpose()*Sigma*a;
            double s = std::sqrt(s2);

            double alpha = (b - m)/s;

            double cdf_alpha = 0.5*(1.0 + std::erf(alpha/std::sqrt(2.0)));
            double pdf_alpha = std::exp(-0.5*alpha*alpha)/std::sqrt(2.0*M_PI);

            // Left-side truncation: a^T x <= b
            double lambda = pdf_alpha/cdf_alpha;

            g.position = mu - lambda*(Sigma*a)/s;
            g.covariance = Sigma - (alpha*lambda + lambda*lambda)*(Sigma*a)*(a.transpose()*Sigma)/s2;
        }

        GaussianObject choppedGaussianAndWeight(GaussianObject& g, std::pair<Eigen::Vector2d,Eigen::Vector2d> seg) {
            // Get the line segment in the form a^T x = b
            Eigen::Vector2d d = seg.second - seg.first;
            Eigen::Vector2d a(-d(1), d(0)); // normal vector
            double b = a.dot(seg.second);

            Eigen::Vector2d mu = g.position;

            // if (a.dot(mu) > b) return;

            Eigen::Matrix2d Sigma = g.covariance;

            double m = a.dot(mu);
            double s2 = a.transpose()*Sigma*a;
            double s = std::sqrt(s2);

            double alpha = (b - m)/s;

            double cdf_alpha = 0.5*(1.0 + std::erf(alpha/std::sqrt(2.0)));
            double pdf_alpha = std::exp(-0.5*alpha*alpha)/std::sqrt(2.0*M_PI);

            // Left-side truncation: a^T x <= b
            double lambda = pdf_alpha/cdf_alpha;

            GaussianObject newGaussian;

            newGaussian.position = mu - lambda*(Sigma*a)/s;

            // g.position =
            //     mu - lambda * (Sigma * a) / s;

            newGaussian.covariance = Sigma - (alpha*lambda + lambda*lambda)*(Sigma*a)*(a.transpose()*Sigma)/s2;
            // g.covariance =
            //     Sigma - (alpha * lambda + lambda * lambda)
            //     * (Sigma * a) * (a.transpose() * Sigma) / s2;

            newGaussian.weight = g.weight*cdf_alpha;
            // g.weight *= cdf_alpha; // Adjust weight based on the probability mass that remains after chopping

            return newGaussian;
        }

        // Linear GSF prediction step
        void predictionStep(const Eigen::Matrix2d& randWalkCov, Point2DCollisionChecker& checker) {
            updatedComponents = components;
            
            for (size_t i = 0; i < components.size(); i++) {

                // update covariance based on random walk
                updatedComponents[i].covariance += randWalkCov; // Add some process noise

                // std::normal_distribution<double> normalDist(0.0, stepSize);

                
                std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>> collision_segments = checker.isCollideEllipse(components[i].position, components[i].covariance, 0.8, true);
                
                // If the predicted Gaussian intersects with an obstacle, seperate into two Gaussians and keep the one that is not colliding
                for (const auto& seg : collision_segments) {
                    choppedGaussian(updatedComponents[i], seg);
                    // // Get the line segment in the form a^T x = b
                    // Eigen::Vector2d d = seg.second - seg.first;
                    // Eigen::Vector2d a(-d(1), d(0)); // normal vector
                    // double b = a.dot(seg.first);

                    // // Use the current updated Gaussian, not always the original one
                    // Eigen::Vector2d mu = updatedComponents[i].position;

                    // if (a.dot(mu) > b) continue; // If the mean is on the left side of the line, then we are good and can skip to the next segment

                    // Eigen::Matrix2d Sigma = updatedComponents[i].covariance;

                    // // // Ensure updated components do not move outside the environment
                    // // updatedComponents[i].position(0) = std::clamp(updatedComponents[i].position(0), 0.0, maxX);
                    // // updatedComponents[i].position(1) = std::clamp(updatedComponents[i].position(1), 0.0, maxY);

                    // double m = a.dot(mu);
                    // double s2 = a.transpose() * Sigma * a;
                    // double s = std::sqrt(s2);

                    // double alpha = (b - m) / s;

                    // double cdf_alpha = 0.5 * (1.0 + std::erf(alpha / std::sqrt(2.0)));
                    // double pdf_alpha = std::exp(-0.5 * alpha * alpha) / std::sqrt(2.0 * M_PI);

                    // // Left-side truncation: a^T x <= b
                    // double lambda = pdf_alpha / cdf_alpha;

                    // updatedComponents[i].position =
                    //     mu - lambda * (Sigma * a) / s;

                    // updatedComponents[i].covariance =
                    //     Sigma - (alpha * lambda + lambda * lambda)
                    //     * (Sigma * a) * (a.transpose() * Sigma) / s2;
                }

                // // Propogate random walk dynamics
                // Eigen::Vector2d newPos;
                // bool valid = false;

                // for (int k = 0; k < 20; k++) {
                //     Eigen::Vector2d noise;
                //     noise << normalDist(gen), normalDist(gen);

                //     newPos = components[i].position + noise;

                //     if (!checker.isCollide(newPos)) {
                //         valid = true;
                //         break;
                //     }
                // }

                // // Check that new component positions are valid
                // if (valid) {
                //     updatedComponents[i].position = newPos;
                // } else {
                //     // Keep old position as a fallback
                //     updatedComponents[i].position = components[i].position;
                // }

            }



            components = updatedComponents;

        }
        std::vector<path_planning::Obstacle2D> getFOVObstacle(const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, Point2DCollisionChecker& checker, double fov, double range) {
            // create point within FOV and check if that collid from center out. Then use the out most point to create obstacle
            
            double startAngle = std::atan2(observerHeading(1), observerHeading(0)) - fov/2.0;
            size_t numAngleSegments = 10; // Number of segments to discretize the FOV into
            double angleInterval = fov/numAngleSegments; // 10 segments in the FOV, can be adjusted for more accuracy at the cost of computation time
            
            size_t numRangeSegments = 20; // Number of segments to discretize the range into
            double rangeInterval = range/numRangeSegments;

            std::vector<size_t> collision_points; // Store the indices of the angle and range segments where collisions occur
            for (size_t i = 0; i < numAngleSegments+1; i++) {
                double angle = startAngle + i*angleInterval;

                for (size_t j = 1; j < numRangeSegments+1; j++) {
                    if (j == numRangeSegments) {
                        collision_points.push_back(j);
                        break;
                    }

                    double r = j*rangeInterval;
                    Eigen::Vector2d point = observerPosition + r*Eigen::Vector2d(std::cos(angle), std::sin(angle));

                    if (checker.isCollide(point)) {
                        collision_points.push_back(j);
                        break; // Stop checking further along this angle segment after the first collision
                    }
                }
            }

            std::vector<path_planning::Obstacle2D> fovObstacles;
            for (size_t i = 0; i < numAngleSegments; i++) {
                double angle1 = startAngle + i*angleInterval;
                double angle2 = startAngle + (i+1)*angleInterval;

                double r = std::min(collision_points[i]*rangeInterval, collision_points[i+1]*rangeInterval);
                Eigen::Vector2d point1 = observerPosition + r*Eigen::Vector2d(std::cos(angle1), std::sin(angle1));
                Eigen::Vector2d point2 = observerPosition + r*Eigen::Vector2d(std::cos(angle2), std::sin(angle2));
                std::vector<Eigen::Vector2d> fovVertices = {observerPosition, point1, point2};
                fovObstacles.push_back(fovVertices);
            }

            return fovObstacles;
        }

        // Linear GSF measurement step
        void measurementStep(const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, Point2DCollisionChecker& checker, double fov, double range) {
            // Eigen::Vector2d leftFOVPoint = observerPosition + range*Eigen::Vector2d(std::cos(std::atan2(observerHeading(1), observerHeading(0)) + fov/2.0), std::sin(std::atan2(observerHeading(1), observerHeading(0)) + fov/2.0));
            // Eigen::Vector2d rightFOVPoint = observerPosition + range*Eigen::Vector2d(std::cos(std::atan2(observerHeading(1), observerHeading(0)) - fov/2.0), std::sin(std::atan2(observerHeading(1), observerHeading(0)) - fov/2.0));
            // std::vector<Eigen::Vector2d> fovVertices = {observerPosition, rightFOVPoint, leftFOVPoint};

            path_planning::Environment2D env = checker.getEnvironment();
            env.obstacles.clear(); // Clear existing obstacles to only consider FOV polygon for collision checking in this step
            
            env.obstacles = getFOVObstacle(observerPosition, observerHeading, checker, fov, range); // Get the FOV polygon as an obstacle

            Point2DCollisionChecker temp_checker(env); // Create a temporary collision checker with the FOV polygon as an obstacle to check whether components are in the FOV, without modifying the original collision checker which is needed for other checks in the filter. Could be optimized by just adding and removing the FOV polygon from the original checker each step, but this is simpler to implement and debug for now. Just make sure not to use the original checker for any checks that need to consider the FOV as an obstacle during the measurement step, and use
            for (auto& g : components) {
                std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>> collision_segments = temp_checker.isCollideEllipse(g.position, g.covariance, 0.6, false);
                std::vector<GaussianObject> newGaussians;
                for (const auto& seg : collision_segments) {
                    newGaussians.push_back(choppedGaussianAndWeight(g, seg));
                }

                // assign g with the new gaussians with highest weight
                if (!newGaussians.empty()) {
                    auto maxIt = std::max_element(newGaussians.begin(), newGaussians.end(), [](const GaussianObject& a, const GaussianObject& b) {
                        return a.weight < b.weight;
                    });
                    g = *maxIt;
                }

                // g.weight *= likelihood(g, observerPosition, observerHeading, checker, fovCosine, range);
            }
        }

        // Find nearest free position for a position which has collided with an obstacle
        Eigen::Vector2d projectToFreeSpace(const Eigen::Vector2d& pos, Point2DCollisionChecker& checker) {
            if (!checker.isCollide(pos)) return pos;

            const double step = 0.05; // search resolution
            const double maxRadius = maxX/5; // search limit

            // Search for open positions
            for (double r = step; r <= maxRadius; r += step) {
                for (double theta = 0; theta < 2*M_PI; theta += M_PI/16.0) {
                    Eigen::Vector2d candidate = pos + r*Eigen::Vector2d(std::cos(theta), std::sin(theta));

                    if (!checker.isCollide(candidate) && candidate(0) >= 0 && candidate(0) <= maxX && candidate(1) >= 0 && candidate(1) <= maxY) {
                        return candidate;
                    }
                }
            }

            // Random valid point as a backup
            Eigen::Vector2d fallback;
            do {
                fallback = {distX(gen), distY(gen)};
            } while (checker.isCollide(fallback));

            return fallback;
        }

        // I tried this method which was in the slides but something is wrong with it and I am not sure what, it breaks the spread of components throughout the environment
        // Merge mixand pairs that are the most statistically similar using a pseudo Mahalabobis distance
        // void mergeComponents(double mahalanobisThreshold) {
        //     std::vector<bool> merged(components.size(), false);
        //     std::vector<GaussianObject> newComponents;

        //     for (size_t i = 0; i < components.size(); i++) {
        //         if (merged[i]) continue;

        //         GaussianObject g = components[i];

        //         for (size_t j = i+1; j < components.size(); j++) {
        //             if (merged[j]) continue;

        //             const auto& gComparison = components[j];

        //             // Position difference between mixands
        //             Eigen::Vector2d posDifference = g.position - gComparison.position;

        //             // Combined covariance
        //             Eigen::Matrix2d covCombined = g.covariance + gComparison.covariance;

        //             // Clamp eigenvalues
        //             Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(covCombined);
        //             Eigen::Vector2d evals = solver.eigenvalues().cwiseMax(1e-2).cwiseMin(1.0);
        //             covCombined = solver.eigenvectors() * evals.asDiagonal() * solver.eigenvectors().transpose();

        //             // Regularization
        //             covCombined += 1e-6*Eigen::Matrix2d::Identity();

        //             double mahalanobisDist = posDifference.transpose()*covCombined.inverse()*posDifference;

        //             // If mixands are similar, merge them into one component
        //             if (mahalanobisDist < mahalanobisThreshold) {
        //                 double weight1 = g.weight;
        //                 double weight2 = gComparison.weight;

        //                 // Merge mean and covariance
        //                 Eigen::Vector2d newPosition = (weight1*g.position + weight2*gComparison.position)/weightSum;
        //                 Eigen::Matrix2d newCovariance = (weight1*(g.covariance + (g.position - newPosition)*(g.position - newPosition).transpose()) + weight2*(gComparison.covariance + (gComparison.position - newPosition)*(gComparison.position - newPosition).transpose()))/weightSum;

        //                 g.position = newPosition;
        //                 g.covariance = newCovariance;
        //                 g.weight = weight1 + weight2;

        //                 merged[j] = true;
        //             }
        //         }

        //         newComponents.push_back(g);
        //     }

        //     components = newComponents;
        // }

        // Merge gaussian components if their means are close together
        void mergeComponents(double distThresh) {
            std::vector<bool> merged(components.size(), false);
            std::vector<GaussianObject> newComponents;

            for (size_t i = 0; i < components.size(); i++) {
                if (merged[i]) continue;

                GaussianObject g = components[i];

                for (size_t j = i+1; j < components.size(); j++) {
                    if (merged[j]) continue;

                    // If mixands are similar, merge them into one component
                    if ((components[i].position - components[j].position).norm() < distThresh) {
                        double weight1 = g.weight;
                        double weight2 = components[j].weight;

                        // Merge mean and covariance
                        g.position = (weight1*g.position + weight2*components[j].position)/(weight1 + weight2);
                        g.covariance = (weight1*g.covariance + weight2*components[j].covariance)/(weight1 + weight2);
                        g.weight = weight1 + weight2;

                        merged[j] = true;
                    }
                }

                newComponents.push_back(g);
            }

            components = newComponents;
        }

        // Generate a random position jitter so split components will not be identical
        Eigen::Vector2d jitterNoise() {
            static std::normal_distribution<double> dist(0.0, 0.1);
            return Eigen::Vector2d(dist(gen), dist(gen));
        }

        // Split components when the number of components drops below the number of expected components
        void splitComponents(size_t targetNumComponents, Point2DCollisionChecker& checker) {
            while (components.size() < targetNumComponents) {

                // Find highest weight component
                auto maxComponent = std::max_element(components.begin(), components.end(), [](const GaussianObject& comp1, const GaussianObject& comp2) { return comp1.weight < comp2.weight; });
                GaussianObject g = *maxComponent;

                // Find the largest eigenvector and spread
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(g.covariance);
                Eigen::Vector2d direction = solver.eigenvectors().col(1); // largest eigenvector
                double spread = std::sqrt(solver.eigenvalues()(1));

                // Create new component objects
                GaussianObject g1 = g;
                GaussianObject g2 = g;

                // Generate new positions with jitter, weights, and covariances
                Eigen::Vector2d noise = jitterNoise();
                g1.position += 0.5*spread*direction + noise;
                g2.position -= 0.5*spread*direction - noise;
                g1.position = projectToFreeSpace(g1.position, checker);
                g2.position = projectToFreeSpace(g2.position, checker);

                g1.weight = g.weight*0.5;
                g2.weight = g.weight*0.5;

                g1.covariance *= 0.8;
                g2.covariance *= 0.8;

                // Replace the original component
                *maxComponent = g1;
                components.push_back(g2);
            }
        }

        // Performs a Linear Gaussian Sum Filter step
        void step(Point2DCollisionChecker& checker, BaseCollisionChecker<Eigen::VectorXd>& lineOfSightChecker, const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, double fov, double fovCosine, double range, const Eigen::Matrix2d& randWalkCov) {

            // Prediction step 
            predictionStep(randWalkCov, checker);

            // Measurement step
            measurementStep(observerPosition, observerHeading, checker, fov, range);

            // Update spatial memory
            updateObservedGrid(observerPosition, observerHeading, lineOfSightChecker, fovCosine, range);

            // Normalize
            normalize();

            // // Re-seed low weight components into free, unobserved, grid cells
            // std::vector<double> cellWeights;
            // std::vector<std::pair<int,int>> cells;
            // cellWeights.reserve(gridSizeX*gridSizeY);
            // cells.reserve(gridSizeX*gridSizeY);

            // for (int i = 0; i < gridSizeX; i++) {
            //     for (int j = 0; j < gridSizeY; j++) {

            //         double weight = std::exp(-observedGrid[i][j]);
            //         Eigen::Vector2d pos((i + 0.5)*maxX/gridSizeX, (j + 0.5)*maxY/gridSizeY);

            //         // Avoid obstacles
            //         if (checker.isCollide(pos)) continue;

            //         cellWeights.push_back(weight);
            //         cells.emplace_back(i,j);
            //     }
            // }

            // std::discrete_distribution<int> cellDist(cellWeights.begin(), cellWeights.end());

            // // Weight threshold to prune
            // const double pruneThresh = 1e-3;

            // for (auto& g : components) {
            //     if (g.weight < pruneThresh) {

            //         int cellIndex = cellDist(gen);
            //         auto [i, j] = cells[cellIndex];

            //         do {
            //             // Add jitter inside the cell
            //             double dx = (double)gen()/gen.max();
            //             double dy = (double)gen()/gen.max();

            //             g.position = Eigen::Vector2d((i + dx)*maxX/gridSizeX, (j + dy)*maxY/gridSizeY);
            //         } while (checker.isCollide(g.position));

            //         // Reset covariance
            //         g.covariance = Eigen::Matrix2d::Identity()*1.0;

            //         // Give uniform small weight
            //         g.weight = 1.0/M;
            //     }
            // }

            // // Spatial hashing grid for component repulsion
            // const double cellSize = 0.3;
            // std::unordered_map<int, std::vector<int>> grid;

            // // Random hash function
            // auto hash = [&](const Eigen::Vector2d& p) {
            //     int x = int(p(0)/cellSize);
            //     int y = int(p(1)/cellSize);
            //     return x*73856093^y*19349663;
            // };

            // // Build the grid
            // for (size_t i = 0; i < components.size(); i++) {
            //     grid[hash(components[i].position)].push_back(i);
            // }

            // // Create repulsion between neighbouring components if they are too close
            // for (auto& [key, indices] : grid) {
            //     for (int i : indices) {
            //         for (int j : indices) {
            //             if (i >= j) continue;

            //             Eigen::Vector2d diff = components[i].position - components[j].position;
            //             double dist = diff.norm();

            //             if (dist < 0.3) {
            //                 Eigen::Vector2d push = 0.05*diff.normalized();
            //                 components[i].position += push;
            //                 components[j].position -= push;
            //             }
            //         }
            //     }
            // }

            // // Maintain the number of components
            // size_t target = M;

            // // Merge components
            // mergeComponents(0.5);

            // // If too many components, prune the lowest weights
            // if (components.size() > target) {
            //     std::sort(components.begin(), components.end(), [](const GaussianObject& comp1, const GaussianObject& comp2) { return comp1.weight > comp2.weight; });
            //     components.resize(target);
            // }

            // // If too few components, split them up
            // while (components.size() < target) {
            //     splitComponents(target, checker);
            // }

            // // Normalize again
            // normalize();

            // // Final check that the components are in open areas
            // for (auto& g : components) {
            //     g.position = projectToFreeSpace(g.position, checker);
            // }
        }

    private:
        int M;
        std::mt19937 gen;
        std::uniform_real_distribution<double> angleDist;
        int gridSizeX;
        int gridSizeY;
        int maxPropogationAttempts;
        double maxObservedScore;
        double observedDecayRate;
        std::vector<std::vector<double>> observedGrid;
        double maxX;
        double maxY;
        std::uniform_real_distribution<double> distX;
        std::uniform_real_distribution<double> distY;
        std::vector<std::vector<double>> densityGrid;
        std::vector<std::vector<double>> surroundingDensityGrid;
        Eigen::Vector2d smoothedEstimate;
        bool hasSmoothedEstimate = false;
};
