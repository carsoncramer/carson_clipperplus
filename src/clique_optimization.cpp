/**
compute core number of graph verticeis, and quickly find a miximal clique.
Code based on Parallel Maximum Clique Algorithm by Ryan A. Rossi (http://ryanrossi.com/pmc)

Author: kaveh fathian (kavehfathian@gmail.com)
 */

#include "clipperplus/clique_optimization.h"

namespace clipperplus {

int clique_optimization(const Eigen::MatrixXd& affinity_matrix, 
                        const Eigen::VectorXd& u0,
                        int& clique_size,
                        std::vector<int>& clique) {

    // run clipper on pruned graph
    clipper::invariants::EuclideanDistance::Params iparams; 
    clipper::invariants::EuclideanDistancePtr invariant =
            std::make_shared<clipper::invariants::EuclideanDistance>(iparams); // instantiate the invariant function
    clipper::Params params; // instantiate parameters     
    clipper::CLIPPER clipper(invariant, params); // clipper object
    // clipper.setMatrixData(M, M);
    clipper.setMatrixData_old(affinity_matrix, affinity_matrix); // identical to matlab

    // find the densest clique of the previously constructed consistency graph
    // clipper.solve(u0);
    // clipper.solveBinary(u0);
    clipper.solveBinary_old(u0);
    clipper::Solution sol = clipper.getSolution();
    
    clique_size = sol.nodes.size(); // clique size
    #ifdef DEBUG
        std::cout << "clipper's clique size: " << clique_size << std::endl;
    #endif
    
    clique = sol.nodes; // indices of clique
    #ifdef DEBUG
        std::cout << "clipper's solution:";
        for (int i=0; i<clique_size; i++) {
            std::cout << clique[i] << " ";
        }
        std::cout << std::endl;
    #endif

    return 1;
}

} // namespace clipperplus
