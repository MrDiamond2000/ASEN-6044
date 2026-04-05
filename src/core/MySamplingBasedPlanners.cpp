# include "core/MySamplingBasedPlanners.hpp"

// MyGenericRRT
path_planning::Path MyGenericRRT::planND(Eigen::VectorXd init_, Eigen::VectorXd goal_, BaseCollisionChecker<Eigen::VectorXd>& collision_checker_){
    graphPtr->clear();
    nodes.clear();
    nodes[0] = init_;
    bool success = false;

    // int temp = 0;
    for (int num_it = 0; num_it < iteration; num_it++){
        int step_num = 0;
        if(double(rand())/RAND_MAX < bias){
            Eigen::VectorXd new_point;
            do{
                new_point = extendRRT(goal_, collision_checker_);
                if(checkDistance(new_point, pow(10,-6))){
                    break;
                }
                step_num++;
            }while(!checkDistance(goal_ - new_point, pow(10,-6)));
            
            if(checkDistance(goal_ - new_point, pow(10,-5))){
                success = true;
                break;
            }

            continue;
        }

        Eigen::VectorXd rand_point = generatePoint(collision_checker_.getBounds());
        Eigen::VectorXd new_point = extendRRT(rand_point, collision_checker_);
        // if (new_point.isZero() || step_num == 0){
        //     num_it--;
        // }
        // if (num_it != temp){
        //     LOG(num_it);
        // }
        // temp = num_it;
    }
    // LOG(temp);

    path_planning::Path path;
    path.valid = false;
    if(success){
        path.valid = true;
        path.waypoints.push_back(goal_);
        path_planning::Node current_node = (graphPtr->parents(nodes.size()-1))[0];
        while(current_node != 0){
            path.waypoints.insert(path.waypoints.begin(), nodes[current_node]);
            current_node = (graphPtr->parents(current_node))[0];
        }
        path.waypoints.insert(path.waypoints.begin(), init_);
    }
    return path;
}

Eigen::VectorXd MyGenericRRT::generatePoint(const std::vector<std::pair<double, double>>& bounds){
    Eigen::VectorXd point;
    point.resize(bounds.size());
    for(size_t i = 0; i < bounds.size(); i++){
        double lower = bounds[i].first;
        double upper = bounds[i].second;
        point(i) = lower + (upper - lower)*double(rand())/RAND_MAX;
    }
    return point;
}

path_planning::Node MyGenericRRT::closestPoint(const Eigen::VectorXd& point){
    double closest_distance = -1;
    path_planning::Node closest_node;
    for (const auto& pair : getNodes()){
        if (closest_distance > magnitude(point-pair.second) || closest_distance == -1){
            closest_distance = magnitude(point-pair.second);
            closest_node = pair.first;
        }
        // if (closest_distance > (point-pair.second).norm() || closest_distance == -1){
        //     closest_distance = (point-pair.second).norm();
        //     closest_node = pair.first;
        // }
    }
    return closest_node;
}

Eigen::VectorXd MyGenericRRT::extendRRT(const Eigen::VectorXd& point, BaseCollisionChecker<Eigen::VectorXd>& collision_checker_){
    path_planning::Node closest_node = closestPoint(point);
    Eigen::VectorXd one_2_two = point - nodes[closest_node];
    Eigen::VectorXd step = one_2_two;
    
    if(checkDistance(one_2_two, pow(10,-6))){
        return step;
    }

    for(int i = 0 ; i < one_2_two.size()/2; i++){
        Eigen::Vector2d temp = {one_2_two(2*i), one_2_two(2*i+1)};
        if(temp.norm() < step_size){
            continue;
        }
        temp.normalize();
        step(2*i) = temp(0)*step_size;
        step(2*i+1) = temp(1)*step_size;
    }

    if(!collision_checker_.isCollide2P(nodes[closest_node], nodes[closest_node]+step)){
        nodes[nodes.size()] = nodes[closest_node]+step;
        graphPtr->connect(closest_node, nodes.size()-1, step.norm());
        return nodes[closest_node]+step;
    }
    step.setZero();
    return step;
}

bool MyGenericRRT::checkDistance(Eigen::VectorXd one_2_two, double requirement){
    std::vector<Eigen::Vector2d> direction;

    for (int i = 0; i < one_2_two.size()/2; i++){
        direction.push_back(Eigen::Vector2d(one_2_two(2*i), one_2_two(2*i+1)));
    }

    for(size_t i = 0; i < direction.size(); i++){
        if(direction[i].norm()>requirement){
            return false;
        }
    }
    return true;
}

double MyGenericRRT::magnitude(Eigen::VectorXd vec){
    double magnitude = 0.0;
    for(int i = 0; i < vec.size()/2; i++){
        magnitude += pow(vec(2*i)*vec(2*i) + vec(2*i+1)*vec(2*i+1),0.5);
    }
    // DEBUG(magnitude);
    return magnitude;
}

// MyRRT
path_planning::Path2D MyRRT::plan(Eigen::Vector2d init_, Eigen::Vector2d goal_) {
    path_planning::Path2D path;
    Point2DCollisionChecker collision_checker(environment);
    path_planning::Path ori_path = planND(init_, goal_, collision_checker);
    path.waypoints = ori_path.getWaypoints2D();
    // path.waypoints.push_back(problem_.q_init);
    // path.waypoints.push_back(problem_.q_goal);
    path.valid = ori_path.valid;
    return path;
}

// std::map<amp::Node, Eigen::Vector2d> MyRRT::getNodes2D(){
//     std::map<amp::Node, Eigen::Vector2d> nodes_2d;
//     for(const auto& pair : getNodes()){
//         nodes_2d.insert({pair.first, Eigen::Vector2d(pair.second(0), pair.second(1))});
//     }
//     return nodes_2d;
// }
