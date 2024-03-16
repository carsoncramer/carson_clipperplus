/*
return a maximal clique in graph by PMC heuristic algorithm
Author: kaveh fathian (kavehfathian@gmail.com)
*/

#include "pmc/pmc_heuristic.h"

namespace pmc {

int pmc_heuristic(const Eigen::MatrixXd& adj,
                  int& clique_size,
                  std::vector<int>& clique) {
    
    const int nnodes = adj.rows(); // number of graph nodes
    const int nedges = adj.sum()/2; // number of graph edges

    #ifdef DEBUG
        std::cout << "\nrunning PMC heuristic algorithm..." << std::endl;
        std::cout << "number of graph nodes: " << nnodes << std::endl;
        std::cout << "number of graph edges: " << nedges << std::endl;

        // std::cout << "\nadjacency matrix:" << std::endl;
        // std::cout << adj << std::endl;
    #endif

    // create a PMC graph from the adjacency matrix
    std::vector<int> edges;
    std::vector<long long> vertices;
    vertices.push_back(0);

    size_t total_num_edges = 0;
    for (size_t i=0; i<nnodes; i++) {
        for (size_t j=0; j<nnodes; j++) {
        if (adj(i,j)==1) {
            edges.push_back(j);
            total_num_edges++;
        }
        }
        vertices.push_back(total_num_edges);
    }

    // construct PMC graph
    pmc::pmc_graph G(vertices, edges);    
    pmc::input in;
   
    // compute k-cores 
    double seconds = get_time();
    G.compute_cores(); // compute core numbers
    in.ub = G.get_max_core() + 1; // max clique upper bound as max k-core + 1
    int kcore_bound = in.ub; // output

    // std::vector<int> kcores = G.kcore; // k-core values for each node
    // core number ouput
    // for(int i = 0; i < nnodes; i++) {
    //     core_numbers[i] = kcores[i] - 1;
    // }

    #ifdef DEBUG
        // std::cout << "G.get_max_core(): " << G.get_max_core() << std::endl;
        std::cout << "max clique upper bound: " << kcore_bound << std::endl;
        std::cout << "k-cores compute time: " << get_time()-seconds << std::endl;

        // std::cout << "core numbers: ";
        // for (int i=0; i < nnodes; i++){std::cout << kcores[i]-1 << " ";}
        // std::cout << std::endl;
    #endif

    // find a heuristic maximal clique using core numbers
    pmc::pmc_heu maxclique(G,in);
    std::vector<int> C; // init maximal clique
    in.lb = maxclique.search(G, C); // find maximal clique
    if (C.size() == 0) {C.push_back(0);} //if graph is disconnected, return vertex 0 as clique
    clique_size = C.size(); // size of maximal clique

    #ifdef DEBUG
        std::cout << "Heuristic found clique of size " << clique_size;
        std::cout << " in " << get_time() - seconds << " seconds" << std::endl;
        
        std::cout << "Heuristic clique: ";
        for(int i = 0; i < clique_size; i++){std::cout << C[i] << " ";}
        std::cout << std::endl;
    
        if (in.lb == in.ub) {
            std::cout << "Heuristic found max clique." << std::endl; 
        }
    #endif
   
    // the heuristic clique ouput
    for(int i = 0; i < clique_size && i < nnodes; i++) {
        clique.push_back(C[i]);
    }
    
    return clique_size; 
}

} // namespace pmc
