#include <time.h>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <map>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>

class ParticleFilter {
    public:
        ParticleFilter(int Ns_, int maxX_, int maxY_, int gridSizeX_, int gridSizeY_)
        :   Ns(Ns_), 
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

        // Define particle structure
        struct Particle {
            Eigen::Vector2d position;
            double weight;
        };

        // Define the list of particles
        std::vector<Particle> particles;

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

        // Check whether the particle (or target) is observed in a defined cone infront of the observer
        bool observed(Eigen::Vector2d position, Eigen::Vector2d observerPosition, Eigen::Vector2d observerHeading, double fovCosine, double range) {

            double dx = position(0) - observerPosition(0);
            double dy = position(1) - observerPosition(1);
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

        // Check whether the target is observed in the conical observable area
        bool detectTarget(const Eigen::Vector2d& truthTargetPosition, const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, BaseCollisionChecker<Eigen::VectorXd>& checker, double fovCosine, double range) {
            if (observed(truthTargetPosition, observerPosition, observerHeading, fovCosine, range) && !checker.isCollide2P(Eigen::VectorXd(observerPosition), Eigen::VectorXd(truthTargetPosition))) {
                return true;
            }
            else {
                return false;
            }
        }

        // Execute a random walk of distance 0.1 in a random direction (provided there are no obstacles)
        // Serves to update the state of each particle and for the true target
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

        // Propogate position of all particles
        void propagateAll(Point2DCollisionChecker& checker, double stepSize) {
            for (auto& p : particles) {
                propagate(p.position, checker, stepSize);
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

        // Compute effective sample size
        double computeNess() {
            double sum = 0.0;
            for (auto& p : particles) {
                sum += p.weight*p.weight;
            }

            return 1.0/sum;
        }

        // Resample particles
        void resample() {
            std::vector<Particle> new_particles;
            new_particles.resize(Ns);
            double step = 1.0/Ns;

            // Pick a random offset
            std::uniform_real_distribution<double> dist(0.0, step);
            double start = dist(gen);

            // Initialize index and starting sum of particle weights
            double sum = particles[0].weight;
            int index = 0;

            // Loop through equally spaced targets
            for (int i = 0; i < Ns; i++) {
                double target = start + i*step;

                // Move through the sum of normalized particle weights until the target is passed
                while (target > sum && index < Ns - 1) {
                    index++;
                    sum += particles[index].weight;
                }

                // Copy the current particle, particles of higher weights will be selected more often
                new_particles[i] = particles[index];
                new_particles[i].weight = step;
            }

            particles = new_particles;
        }

        // Calculate the MMSE estimate based on the particles
        //Eigen::Vector2d estimate() {
        //    Eigen::Vector2d mmse(0,0);

        //    for (auto& p : particles) {
        //        mmse += p.weight*p.position;
        //    }

        //    return mmse;
        //}

        // Calculate the target position estimate based on the highest weight particle and least searched area
        //Eigen::Vector2d estimate(double lambda = 1.0) {
        //    Eigen::Vector2d bestWeight = bestParticle();
        //    Eigen::Vector2d bestLocation = bestUnexplored();

        //    return (1.0 - lambda)*bestWeight + lambda*bestLocation;
        //}

        // Finds the position of the highest weight particle
        //Eigen::Vector2d bestParticle() {
        //    double maxW = 0.0;
        //    Eigen::Vector2d bestParticlePosition;

        //    for (auto& p : particles) {
        //        if (p.weight > maxW) {
        //            maxW = p.weight;
        //            bestParticlePosition = p.position;
        //        }
        //    }

        //    return bestParticlePosition;
        //}

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

        std::pair<double, Eigen::Vector2d> estimate_density() {
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

            // for (int i = 0; i < gridSizeX; i++) {
            //     for (int j = 0; j < gridSizeY; j++) {
            //         if (densityGrid[i][j] > maxDensity) {
            //             maxDensity = densityGrid[i][j];
            //             bestParticlePosition = Eigen::Vector2d((i + 0.5)*maxX/static_cast<double>(gridSizeX), (j + 0.5)*maxY/static_cast<double>(gridSizeY));
            //         }
            //     }
            // }

            for (int i = 0; i < gridSizeX; i++) {
                for (int j = 0; j < gridSizeY; j++) {
                    if (surroundingDensityGrid[i][j] > maxDensity) {
                        maxDensity = surroundingDensityGrid[i][j];
                        bestParticlePosition = Eigen::Vector2d((i + 0.5)*maxX/static_cast<double>(gridSizeX), (j + 0.5)*maxY/static_cast<double>(gridSizeY));
                    }
                }
            }

            return {maxDensity*gridSizeX*gridSizeY/(10*10), bestParticlePosition};
        }

        // Performs a Sequential Importance Sampling (SIS) Particle Filter step
        void step(Point2DCollisionChecker& checker, BaseCollisionChecker<Eigen::VectorXd>& lineOfSightChecker, const Eigen::Vector2d& observerPosition, const Eigen::Vector2d& observerHeading, double fovCosine, double range, double stepSize, double resampleThreshold) {

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
};