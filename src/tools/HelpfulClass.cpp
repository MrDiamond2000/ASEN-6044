#include "tools/HelpfulClass.hpp"

Point2DCollisionChecker::Point2DCollisionChecker(const path_planning::Environment2D& enviroument_)
: env(enviroument_){
    std::vector<std::pair<double,double>> bounds;
    bounds.push_back({enviroument_.x_min, enviroument_.x_max});
    bounds.push_back({enviroument_.y_min, enviroument_.y_max});
    setBounds(bounds);
}

bool Point2DCollisionChecker::isCollide(const Eigen::VectorXd& point) {
    for (size_t ob_idx = 0; ob_idx < env.obstacles.size(); ++ob_idx) {
        if (thisObstacle(ob_idx, point)){
            return true; // Point is inside this polygons
        }
    }
    return false; // Point is outside all polygons
    // Implementation
}

bool Point2DCollisionChecker::isCollide2P(const Eigen::VectorXd& point1_, const Eigen::VectorXd& point2_){
    Eigen::VectorXd one_2_two = point2_ - point1_;
    double distance = one_2_two.norm();
    int sec_num = 1;
    while(distance/sec_num > GAP){
        for(int i = 0; i < sec_num; i++){
            Eigen::VectorXd check_location = point1_ + one_2_two*(1+2*i)/(sec_num*2);
            if (isCollide(check_location))
            {
                return true;
            }
        }
        sec_num = sec_num * 2;
    }
    return false;
}

std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>> Point2DCollisionChecker::isCollideEllipse(const Eigen::Vector2d& center, const Eigen::Matrix2d& covariance, double confidence, bool with_bound) const {
    std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>> collision_segments;

    double s = -2 * log(1 - confidence);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(covariance * s);
    Eigen::Vector2d eigenvalues = es.eigenvalues();
    Eigen::Matrix2d eigenvectors = es.eigenvectors();

    // find angle the eigenvector are from the x-axis
    double angle = atan2(eigenvectors(1, 0), eigenvectors(0, 0));

    for (size_t ob_idx = 0; ob_idx < env.obstacles.size(); ++ob_idx) {
        const path_planning::Obstacle2D& obstacle = env.obstacles[ob_idx];
        int n = obstacle.verticesCCW().size();
        for (int vertix_idx = 0; vertix_idx < n; ++vertix_idx) {
            const Eigen::VectorXd& v1 = obstacle.verticesCCW()[vertix_idx];
            const Eigen::VectorXd& v2 = obstacle.verticesCCW()[(vertix_idx + 1) % n];

            // transform the line segment to the ellipse's coordinate frame
            Eigen::Matrix2d rotation;
            rotation << cos(-angle), -sin(-angle),
                        sin(-angle), cos(-angle);

            // Transform segment endpoints into ellipse frame
            Eigen::Vector2d p1 = rotation * (v1 - center);
            Eigen::Vector2d p2 = rotation * (v2 - center);

            Eigen::Vector2d d = p2 - p1;

            double a2 = eigenvalues(0);
            double b2 = eigenvalues(1);

            // Degenerate segment: v1 == v2
            if (d.squaredNorm() < 1e-12) {
                double value = p1(0) * p1(0) / a2 + p1(1) * p1(1) / b2;
                if (value <= 1.0) {
                    collision_segments.push_back({v1, v2});
                }
            }

            // Quadratic coefficients for:
            // ((p1_x + t d_x)^2 / a2) + ((p1_y + t d_y)^2 / b2) = 1
            double A = d(0) * d(0) / a2 + d(1) * d(1) / b2;
            double B = 2.0 * (p1(0) * d(0) / a2 + p1(1) * d(1) / b2);
            double C = p1(0) * p1(0) / a2 + p1(1) * p1(1) / b2 - 1.0;

            double discriminant = B * B - 4.0 * A * C;

            if (discriminant < 0.0) {
                // No intersection with the infinite line
                continue;
            }

            double sqrt_disc = std::sqrt(discriminant);

            double t1 = (-B - sqrt_disc) / (2.0 * A);
            double t2 = (-B + sqrt_disc) / (2.0 * A);

            // Collision if either intersection point lies on the segment
            if ((t1 >= 0.0 && t1 <= 1.0) || (t2 >= 0.0 && t2 <= 1.0)) {
                collision_segments.push_back({v1, v2});
            }


            // Eigen::Vector2d transformed_v1 = rotation * (v1 - center);
            // Eigen::Vector2d transformed_v2 = rotation * (v2 - center);

            // // create a line for y = mx + c form
            // double m = (transformed_v2(1) - transformed_v1(1)) / (transformed_v2(0) - transformed_v1(0));
            // double c = transformed_v1(1) - m * transformed_v1(0);

            // if (c*c < eigenvalues(1) + m*m*eigenvalues(0)) {
            //     collision_segments.push_back({v1, v2});
            // }
        }
    }

    if (!with_bound) {
        return collision_segments;
    }

    std::vector<Eigen::Vector2d> temp_bounds = {Eigen::Vector2d(env.x_min, env.y_min), Eigen::Vector2d(env.x_min, env.y_max), Eigen::Vector2d(env.x_max, env.y_max), Eigen::Vector2d(env.x_max, env.y_min), Eigen::Vector2d(env.x_min, env.y_min)};

    for (size_t idx = 0; idx < 4; ++idx) {
        const Eigen::VectorXd& v1 = temp_bounds[idx];
        const Eigen::VectorXd& v2 = temp_bounds[idx+1];

        // transform the line segment to the ellipse's coordinate frame
        Eigen::Matrix2d rotation;
        rotation << cos(-angle), -sin(-angle),
                    sin(-angle), cos(-angle);

        // Transform segment endpoints into ellipse frame
        Eigen::Vector2d p1 = rotation * (v1 - center);
        Eigen::Vector2d p2 = rotation * (v2 - center);

        Eigen::Vector2d d = p2 - p1;

        double a2 = eigenvalues(0);
        double b2 = eigenvalues(1);

        // Degenerate segment: v1 == v2
        if (d.squaredNorm() < 1e-12) {
            double value = p1(0) * p1(0) / a2 + p1(1) * p1(1) / b2;
            if (value <= 1.0) {
                collision_segments.push_back({v1, v2});
            }
        }

        // Quadratic coefficients for:
        // ((p1_x + t d_x)^2 / a2) + ((p1_y + t d_y)^2 / b2) = 1
        double A = d(0) * d(0) / a2 + d(1) * d(1) / b2;
        double B = 2.0 * (p1(0) * d(0) / a2 + p1(1) * d(1) / b2);
        double C = p1(0) * p1(0) / a2 + p1(1) * p1(1) / b2 - 1.0;

        double discriminant = B * B - 4.0 * A * C;

        if (discriminant < 0.0) {
            // No intersection with the infinite line
            continue;
        }

        double sqrt_disc = std::sqrt(discriminant);

        double t1 = (-B - sqrt_disc) / (2.0 * A);
        double t2 = (-B + sqrt_disc) / (2.0 * A);

        // Collision if either intersection point lies on the segment
        if ((t1 >= 0.0 && t1 <= 1.0) || (t2 >= 0.0 && t2 <= 1.0)) {
            collision_segments.push_back({v1, v2});
        }
    }

    return collision_segments;
}

bool Point2DCollisionChecker::thisObstacle(int ob_idx_, const Eigen::VectorXd& point_){
    const path_planning::Obstacle2D& obstacle = env.obstacles[ob_idx_];
    double cross_product;

    int n = obstacle.verticesCCW().size();
    for (int vertix_idx = 0; vertix_idx < n; ++vertix_idx) {
        const Eigen::VectorXd& v1 = obstacle.verticesCCW()[vertix_idx];
        const Eigen::VectorXd& v2 = obstacle.verticesCCW()[(vertix_idx + 1) % n];
        Eigen::VectorXd edge = v1 - v2;
        Eigen::VectorXd point_to_v1 = point_ - v1;
        cross_product = edge(0) * point_to_v1(1) - edge(1) * point_to_v1(0);
        if (cross_product > 0){
            return false; // Point is outside of this obstacle
        }
    }
    return true; // Point is inside this obstacle
}

Point2DCollisionCheckerGrid::Point2DCollisionCheckerGrid(const path_planning::GridCSpace2D_T<int8_t>& map_)
: map(map_){
    std::vector<std::pair<double,double>> bounds;
    bounds.push_back({map.x0Bounds().first, map.x0Bounds().second});
    bounds.push_back({map.x1Bounds().first, map.x1Bounds().second});
    setBounds(bounds);
}

bool Point2DCollisionCheckerGrid::isCollide(const Eigen::VectorXd& point_){
    return(map.inCollision(point_(0), point_(1)));
}

bool Point2DCollisionCheckerGrid::isCollide2P(const Eigen::VectorXd& point1_, const Eigen::VectorXd& point2_){
    Eigen::VectorXd one_2_two = point2_ - point1_;
    double distance = one_2_two.norm();
    int sec_num = 1;
    while(distance/sec_num > GAP){
        for(int i = 0; i < sec_num; i++){
            Eigen::VectorXd check_location = point1_ + one_2_two*(1+2*i)/(sec_num*2);
            if (isCollide(check_location))
            {
                return true;
            }
        }
        sec_num = sec_num * 2;
    }
    return false;
}

bool MultiAgentPoint2DCollisionCheckerGrid::isCollide(int agent_idx_, const Eigen::VectorXd& point_){
    return(maps[agent_idx_].inCollision(point_(0), point_(1)));
}

bool MultiAgentPoint2DCollisionCheckerGrid::isCollide2P(int agent_idx_, const Eigen::VectorXd& point1_, const Eigen::VectorXd& point2_){
    Eigen::VectorXd one_2_two = point2_ - point1_;
    double distance = one_2_two.norm();
    int sec_num = 1;
    while(distance/sec_num > GAP){
        for(int i = 0; i < sec_num; i++){
            Eigen::VectorXd check_location = point1_ + one_2_two*(1+2*i)/(sec_num*2);
            if (isCollide(agent_idx_, check_location))
            {
                return true;
            }
        }
        sec_num = sec_num * 2;
    }
    return false;
}
