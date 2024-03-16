#pragma once

#include <iostream>
#include <chrono>
#include <Eigen/Dense>

#include "clipperplus/clique_corenumber.h"
#include "clipperplus/clique_optimization.h"
#include "clipperplus/utils.h"


namespace clipperplus {

int clipperplus_clique(const Eigen::MatrixXd& adj,
                       int& clique_size,
                       std::vector<int>& clique,
                       int& certificate);

} // namespace clipperplus
