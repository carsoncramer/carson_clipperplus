/*
c++ tests for clipper+ library
*/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <Eigen/Dense>

// for computating core numbers and PMC-heuristic clique 
#include "clipperplus/clipperplus_clique.h"
#include "pmc/pmc_heuristic.h"

#include "clipperplus/utils.h"


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void test_optimization() {
    std::cout << "\ntesting clipper with affinity matrix:" << std::endl;

    // graph affinity matrix (adjacency matrix+ identity)
    Eigen::MatrixXd M(5,5);
    M << 1, 0, 1, 1, 1,
         0, 1, 0, 1, 1,
         1, 0, 1, 1, 1,
         1, 1, 1, 1, 1,
         1, 1, 1, 1, 1;       
    std::cout << M << std::endl;

    // instantiate the invariant function
    clipper::invariants::EuclideanDistance::Params iparams;
    clipper::invariants::EuclideanDistancePtr invariant =
            std::make_shared<clipper::invariants::EuclideanDistance>(iparams);
    clipper::Params params;
    // clipper object 
    clipper::CLIPPER clipper(invariant, params);
    clipper.setMatrixData(M, M);

    // find the densest clique of the previously constructed consistency graph
    clipper.solve();

    clipper::Solution sol = clipper.getSolution();
    std::cout << "clipper's solution:";
    for (int i=0; i<sol.nodes.size(); i++) {
        std::cout << sol.nodes[i] << " ";
    }
    std::cout << std::endl;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
int test_pmcheu() {
    std::cout << "\ntesting pmc heuristic with adjacency matrix:" << std::endl;
    Eigen::MatrixXd adj(5,5); // graph affinity matrix (adjacency matrix+ identity)
    adj << 0, 0, 1, 1, 1,
           0, 0, 0, 1, 1,
           1, 0, 0, 1, 1,
           1, 1, 1, 0, 1,
           1, 1, 1, 1, 0; 
                 
    std::cout << adj << std::endl;

    const int nnodes = adj.rows(); // number of graph nodes
    // run pmu heuristic    
    int clique_size = 0;
    std::vector<int> clique;
    pmc::pmc_heuristic(adj, clique_size, clique);

    std::cout << "pmc heuristic clique size: " << clique_size << std::endl;
    std::cout << "pmc heuristic clique: ";
    for (int i : clique) {std::cout << i << " ";}
    std::cout << std::endl;    
    return 1;
}



/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
int test_clipperplus(const Eigen::MatrixXd& adj) {
    std::cout << "\ntesting clipper-plus." << std::endl;

    const int nnodes = adj.rows(); // number of graph nodes
    // Eigen::MatrixXd M = adj + Eigen::MatrixXd::Identity(nnodes,nnodes); // affinity matrix

    // run clipperplus_clique    
    int clique_size = 0;
    std::vector<int> clique;
    Certificate certificate = NONE;
    clipperplus::clipperplus_clique(adj, clique_size, clique, certificate);

    std::cout << "clipperplus clique size: " << clique_size << std::endl;
    std::cout << "clipperplus clique: ";
    for (int i : clique) {std::cout << i << " ";}
    std::cout << std::endl;
    std::cout << "clipperplus certification type: " << cert_to_string(certificate) << std::endl;
    
    return 1;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
    std::cout << "adjacency matrix:" << std::endl;

    // Eigen::MatrixXd adj(10,10); // graph affinity matrix 
    // adj << 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
    //        0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
    //        1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
    //        1, 1, 1, 0, 1, 1, 1, 0, 1, 1,
    //        1, 1, 1, 1, 0, 1, 0, 1, 1, 0,
    //        1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
    //        1, 1, 1, 1, 0, 1, 0, 1, 1, 1,
    //        1, 1, 1, 0, 1, 1, 1, 0, 1, 1,
    //        1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
    //        1, 1, 1, 1, 0, 1, 1, 1, 1, 0;       

    // Eigen::MatrixXd adj(10,10);
    // adj << 0, 0, 1, 1, 1, 1, 1, 0, 1, 0,
    //        0, 0, 1, 1, 1, 0, 1, 1, 1, 1,
    //        1, 1, 0, 1, 0, 1, 1, 1, 0, 1,
    //        1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
    //        1, 1, 0, 1, 0, 0, 1, 1, 1, 1,
    //        1, 0, 1, 1, 0, 0, 1, 1, 1, 1,
    //        1, 1, 1, 1, 1, 1, 0, 1, 1, 0,
    //        0, 1, 1, 1, 1, 1, 1, 0, 1, 1,
    //        1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
    //        0, 1, 1, 1, 1, 1, 0, 1, 1, 0;

    Eigen::MatrixXd adj(20,20); // graph affinity matrix
    adj << 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
           0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1,
           1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0,
           0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1,
           0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1,
           1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
           0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0,
           0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
           1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1,
           0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1,
           1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1,
           1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
           1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0,
           1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
           0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
           1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0;

    std::cout << adj << std::endl;

    // test_optimization(); // test old clipper optimization clique
    // int clique_size_heu = test_pmcheu(); // test pmc heuristic algorithm
    int clique_size = test_clipperplus(adj); // test clipper_plus
    
};
