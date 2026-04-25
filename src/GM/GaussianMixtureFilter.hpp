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
        GaussianMixtureFilter(int Ns_, int maxX_, int maxY_, int gridSizeX_, int gridSizeY_)
        :   Ns(Ns_), 
            gen(std::random_device{}()), 
            angleDist(0, 2*M_PI), 
            gridSizeX(gridSizeX_), 
            gridSizeY(gridSizeY_), 
            maxX(maxX_), maxY(maxY_), 
            distX(0.0, maxX), 
            distY(0.0, maxY) {
                densityGrid.resize(gridSizeX, std::vector<double>(gridSizeY, 0.0));
                surroundingDensityGrid.resize(gridSizeX, std::vector<double>(gridSizeY, 0.0));
        }

        // Define GM structure
        struct GaussianObject {
            Eigen::Vector2d position;
            Eigen::Matrix2d covariance;
            double weight;

            // define = as copying everything
            GaussianObject& operator=(const GaussianObject& other) {
                position = other.position;
                covariance = other.covariance;
                weight = other.weight;
                return *this;
            }
        };

        // Define the list of particles
        std::vector<GaussianObject> gmObjects;
        std::vector<GaussianObject> updatedgmObjects;

        // Randomely generate initial particle positions and add corresponding particles, then return a randomely generated target
        Eigen::Vector2d initializeParticles(Point2DCollisionChecker& checker) {
            for (int i = 0; i < Ns; i++) {
                Eigen::Vector2d pos;
                do { pos = {distX(gen), distY(gen)}; } while (checker.isCollide(pos));
                addParticle(pos, 1.0/Ns);
            }

            Eigen::Vector2d target;
            do { target = {distX(gen), distY(gen)}; } while (checker.isCollide(target));
            return target;
        }

        // Add a particle to the filter
        void addParticle(Eigen::Vector2d position, double weight) {
            particles.push_back({position, weight});
        }

        // Check whether the particle (or target) is observed in a defined cone infront of the observer
        bool observed(Eigen::Vector2d position, Eigen::Vector2d observerPosition, Eigen::Vector2d observerHeading, double fovCosine, double range) {

            double dx = position(0) - observerPosition(0);
            double dy = position(1) - observerPosition(1);

            if (dx > range || dx < -range || dy > range || dy < -range) return false;

            double distanceSquared = dx*dx + dy*dy;

            if (distanceSquared > range*range) return false;

            double dot = dx*observerHeading(0) + dy*observerHeading(1);
            return dot > sqrt(distanceSquared)*fovCosine;
        }

        // Set the weights of any particles in the conical observable area to zero, otherwise one
        void updateWeights(const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, BaseCollisionChecker<Eigen::VectorXd>& checker, double fovCosine, double range) {
            for (auto& p : particles) {
                if (observed(p.position, observerPosition, observerHeading, fovCosine, range) && !checker.isCollide2P(Eigen::VectorXd(observerPosition), Eigen::VectorXd(p.position))) {
                    p.weight = 0;
                }
            }
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
            
            for (auto& p : particles) sum += p.weight;

            if (sum == 0.0) {
                double uniformWeight = 1.0/particles.size();
                for (auto& p : particles) p.weight = uniformWeight;
                return;
            }

            for (auto& p : particles) p.weight /= sum;
        }

        // Finds the position of the combined highest weight and most unobserved particle
        Eigen::Vector2d estimate() {
            double maxScore = 0.0;
            Eigen::Vector2d bestParticlePosition;

            for (auto& p : particles) {
                double observedScore = getObservedScore(p.position);
                double score = p.weight*exp(-observedScore) + 0.5*exp(-observedScore);

                if (score > maxScore) {
                    maxScore = score;
                    bestParticlePosition = p.position;
                }
            }

            return bestParticlePosition;
        }

        std::pair<double, Eigen::Vector2d> estimate_density(Point2DCollisionChecker& checker) {
            // Reset density grid
            for (int i = 0; i < gridSizeX; i++) {
                for (int j = 0; j < gridSizeY; j++) {
                    densityGrid[i][j] = 0.0;
                    surroundingDensityGrid[i][j] = 0.0;
                }
            }

            // Add particle weights to density grid
            for (auto& p : particles) {
                int ix = std::min(std::max(int((p.position(0)/maxX)*gridSizeX), 0), gridSizeX - 1);
                int iy = std::min(std::max(int((p.position(1)/maxY)*gridSizeY), 0), gridSizeY - 1);
                densityGrid[ix][iy] += p.weight;
            }

            // Add surrounding cells to density grid to create a smoother estimate
            for (int i = 0; i < gridSizeX; i++) {
                for (int j = 0; j < gridSizeY; j++) {
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -10; dj <= 10; dj++) {
                            int ni = i + di;
                            int nj = j + dj;
                            if (ni >= 0 && ni < gridSizeX && nj >= 0 && nj < gridSizeY) {
                                surroundingDensityGrid[i][j] += densityGrid[ni][nj];
                            }
                        }
                    }
                }
            }

            // Find cell with highest density
            double maxDensity = 0.0;
            Eigen::Vector2d bestParticlePosition;
            Eigen::Vector2d tempBestParticlePosition;


            for (int i = 0; i < gridSizeX; i++) {
                for (int j = 0; j < gridSizeY; j++) {
                    if (surroundingDensityGrid[i][j] > maxDensity) {
                        tempBestParticlePosition = Eigen::Vector2d((i + 0.5)*maxX/static_cast<double>(gridSizeX), (j + 0.5)*maxY/static_cast<double>(gridSizeY));
                        if (checker.isCollide(tempBestParticlePosition)) continue; // Skip if cell is in an obstacle   

                        maxDensity = surroundingDensityGrid[i][j];
                        bestParticlePosition = Eigen::Vector2d((i + 0.5)*maxX/static_cast<double>(gridSizeX), (j + 0.5)*maxY/static_cast<double>(gridSizeY));
                    }
                }
            }

            return {maxDensity*gridSizeX*gridSizeY/(10*10), bestParticlePosition};
        }

        void predict(Eigen::Matrix2d randWalkCov, Point2DCollisionChecker& checker) {
            for (size_t i = 0; i < gmObjects.size(); i++) {
                updatedgmObjects[i] = gmObjects[i];

                // update covariance base on random walk
                updatedgmObjects[i].covariance += randWalkCov; // Add some process noise

                
                std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>> collition_sagments = checker.isCollideEllipse(gmObjects[i].position, gmObjects[i].covariance, 0.95);
                
                // If the predicted Gaussian intersects with an obstacle, seperate into two Gaussians and keep the one that is not colliding
                for (const auto& seg : collition_sagments) {
                    // Get the line segment in the form a^T x = b
                    Eigen::Vector2d d = seg.second - seg.first;
                    Eigen::Vector2d a(-d(1), d(0));   // normal vector
                    double b = a.dot(seg.first);

                    // Use the current updated Gaussian, not always the original one
                    Eigen::Vector2d mu = updatedgmObjects[i].position;
                    Eigen::Matrix2d Sigma = updatedgmObjects[i].covariance;

                    double m = a.dot(mu);
                    double s2 = a.transpose() * Sigma * a;
                    double s = std::sqrt(s2);

                    double alpha = (b - m) / s;

                    double cdf_alpha = 0.5 * (1.0 + std::erf(alpha / std::sqrt(2.0)));
                    double pdf_alpha = std::exp(-0.5 * alpha * alpha) / std::sqrt(2.0 * M_PI);

                    // Left-side truncation: a^T x <= b
                    double lambda = pdf_alpha / cdf_alpha;

                    updatedgmObjects[i].position =
                        mu - lambda * (Sigma * a) / s;

                    updatedgmObjects[i].covariance =
                        Sigma - (alpha * lambda + lambda * lambda)
                        * (Sigma * a) * (a.transpose() * Sigma) / s2;
                }

            }

        }

        // Performs a Sequential Importance Sampling (SIS) Particle Filter step
        void step(Point2DCollisionChecker& checker, BaseCollisionChecker<Eigen::VectorXd>& lineOfSightChecker, const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, double fovCosine, double range, Eigen::Matrix2d randWalkCov, double resampleThreshold) {

            // Prediction step: 
            predict(randWalkCov, checker);

            // Propogate particles
            propagateAll(checker, stepSize);

            // Update weights
            updateWeights(observerPosition, observerHeading, lineOfSightChecker, fovCosine, range);

            // Update spatial memory
            updateObservedGrid(observerPosition, observerHeading, lineOfSightChecker, fovCosine, range);

            // Normalize
            normalize();

            // Compute Ness
            double Ness = computeNess();

            // Resample and add jitter
            if (Ness < resampleThreshold*Ns) {
                resample();
            }
        }

    private:
        int Ns;
        std::mt19937 gen;
        std::uniform_real_distribution<double> angleDist;
        int gridSizeX;
        int gridSizeY;
        double maxX;
        double maxY;
        std::uniform_real_distribution<double> distX;
        std::uniform_real_distribution<double> distY;
        std::vector<std::vector<double>> densityGrid;
        std::vector<std::vector<double>> surroundingDensityGrid;
};