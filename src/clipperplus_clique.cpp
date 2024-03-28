/*
computes a maximal clique in graph, and certifies if it's maximum clique

Author: kaveh fathian (kavehfathian@gmail.com)
 */

#include "clipperplus/clipperplus_clique.h"
#include "clipperplus/utils.h"

namespace clipperplus {

int clipperplus_clique(const Eigen::MatrixXd& adj,
                       int& clique_size,
                       std::vector<int>& clique,
                       Certificate& certificate) {
    
    certificate = NONE; // initialize to 0

    const int nnodes = adj.rows(); // number of graph nodes
    const int nedges = adj.sum()/2; // number of graph edges    
 
    // affinity matrix (adjacency matrix+ identity)
    Eigen::MatrixXd affinity_matrix = adj + Eigen::MatrixXd::Identity(nnodes,nnodes); 

    #ifdef DEBUG_TIMING
        const auto t1 = std::chrono::high_resolution_clock::now(); // timer
    #endif

    // calculate a maximal clique via core numbers
    std::vector<int> clique_core; // initialize index vector of clique
    clique_core.reserve(nnodes); // allocate memory
    std::vector<int> core_numbers(nnodes,0); // initialize vector of core numbers
    int core_bound = 0; // initialize max clique upper bound based on max kcore
    int chromatic_bound = 0;; // initialize chromatic number upper bound
    std::vector<int> node_colors(nnodes, 0); // initialize graph node coloring
    
    int clique_size_core = clipperplus::clique_corenumber(adj, clique_core, 
                                            core_numbers, core_bound, 
                                            node_colors, chromatic_bound);
    
    #ifdef DEBUG_TIMING
        const auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
        double elapsed = static_cast<double>(duration.count()) / 1e6;
        std::cout << "time for heuristic clique & core numbers: " << elapsed << std::endl;
    #endif
    
    if (clique_size_core == core_bound) {
        clique = clique_core;
        clique_size = clique_size_core;
        certificate = PRUNING; // certified based on pruning heuristic clique
        return clique_size_core;
    }     

    // prune graph: index of the nodes to prune or keep
    std::vector<int> idx_prune, idx_keep;
    idx_prune.reserve(nnodes); // reserve memory
    idx_keep.reserve(nnodes); // reserve memory
    for (int i=0; i<nnodes; i++) {
        if (core_numbers[i] < clique_size_core) { 
            idx_prune.push_back(i); // node with core number lower than the clique size
        } else {
            idx_keep.push_back(i);
        }
    }

    // pruned graph affinity matirx
    Eigen::MatrixXd affinity_matrix_pruned = affinity_matrix(idx_keep, idx_keep);

    // create initial vector for optimization
    Eigen::VectorXd u0(idx_keep.size()); ///< initial vector used for local solver
    u0.setOnes(); // set all elements to 1

    // remove heuristic clique indices
    for(int i = 0; i < clique_size_core; i++){
        // find index of the element in idx_keep that is equal to clique_core_array[i]
        int idx = clipperplus::find_index(idx_keep, clique_core[i]);
        u0[idx] = 0; // remove heuristic clique 
    }

    // Normalize the vector
    u0.normalize();

    #ifdef DEBUG_TIMING
        const auto t3 = std::chrono::high_resolution_clock::now(); // timer
        duration = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
        elapsed = static_cast<double>(duration.count()) / 1e6;
        std::cout << "time for graph pruning: " << elapsed << std::endl;
    #endif

    // calculate a maximal clique on pruned graph via optimization
    std::vector<int> clique_optim_pruned; // indices of clique
    int clique_size_optim = 0; // initialize size of clique
    
    // find clique on pruned graph via optimization
    clipperplus::clique_optimization(affinity_matrix_pruned, u0, clique_size_optim, clique_optim_pruned);

    #ifdef DEBUG_TIMING
        const auto t4 = std::chrono::high_resolution_clock::now(); // timer
        duration = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);
        elapsed = static_cast<double>(duration.count()) / 1e6;
        std::cout << "time for clipper optimization: " << elapsed << std::endl;
    #endif

    // map pruned idicies back to original idicies
    std::vector<int> clique_optim;
    clique_optim.reserve(nnodes); // reserve memory
    for (int i : clique_optim_pruned) {
        if (i >= 0 && i < idx_keep.size()) {  // Check if the index is valid
            clique_optim.push_back(idx_keep[i]);
        } else {
            std::cerr << "index " << i << " is out of bounds." << std::endl;
        }
    }

    // check if clique is certifiable as max clique
    clique_size = 0;
    if (clique_size_optim >= clique_size_core) {
        clique_size = clique_size_optim;
        clique = clique_optim;
        #ifdef DEBUG
            std::cout << "optimization gave better or equal clique\n" << std::endl;
        #endif
    } else {
        clique_size = clique_size_core;
        clique = clique_core;
        #ifdef DEBUG
            std::cout << "optimization returned a worst clique!\n" << std::endl;
        #endif
    }
    
    if (clique_size == core_bound) {
        certificate = KCORE; // certified based on max kcore 
        #ifdef DEBUG
            std::cout << "max clique found; certified by max kcore." << std::endl;
        #endif
    } else if (clique_size == chromatic_bound) {
        certificate = CHROMATIC; // certified based on chromatic number
        #ifdef DEBUG
            std::cout << "max clique found; certified by chromatic number." << std::endl;
        #endif
    } else {
        #ifdef DEBUG
            std::cout << "could not certify max clique" << std::endl; 
        #endif
    }
    
    #ifdef DEBUG_TIMING
        const auto t5 = std::chrono::high_resolution_clock::now(); // timer
        duration = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4);
        elapsed = static_cast<double>(duration.count()) / 1e6;
        std::cout << "time for certification: " << elapsed << std::endl;
    #endif

    // return results
    return clique_size; // return best clique
}

} // namespace clipperplus
