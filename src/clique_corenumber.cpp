/*
compute core number of graph nodes, and use them to find a miximal clique.
Code inspired by PMC Algorithm by Ryan A. Rossi (http://ryanrossi.com/pmc)
Author: kaveh fathian (fathian@ariarobotics.com)
*/

#include "clipperplus/clique_corenumber.h"

namespace clipperplus {

///////////////////////////////////////////////////////////////////////////////////

// ALGORITHM 1: DEGENERACY-ORDERED GREEDY MAXIMAL CLIQUE HEURISTIC
int find_heuristic_clique(const Eigen::MatrixXd& adj, 
                          const std::vector<int>& core_numbers,
                          std::vector<int>& clique) {
    
    const int nnodes = adj.cols(); // number of nodes

    // if graph has no edges, return trivial clique
    if (adj.sum() == 0) { 
        clique.push_back(0);
        return 1;
    }
    // Line 3, initialize max value to 0
    int max_val = 0; 
    std::vector<std::pair<int,int>> nodes; // .first = id, .second = core number
    nodes.reserve(nnodes);
    for(int i = 0; i < nnodes; i++){
        nodes.emplace_back(i, core_numbers[i]);
    }

    //Line 5, sort vertices in decending order by core number
    std::stable_sort(nodes.begin(),nodes.end(), [&](std::pair<int,int> a, std::pair<int,int> b) {
        return a.second > b.second; });

    //Line 6
    for(std::pair<int,int> node : nodes){
        //Line 7
        if(node.second > max_val){
            // Line 9, neighbors of node i with a core number > max_val
            std::vector<std::pair<int,int>> nodes_subset; // .first = id, .second = core number
            for(std::pair<int,int> n : nodes){
                if(n.second >= max_val && adj(node.first, n.first) != 0){
                    nodes_subset.push_back(n);
                }
            } if(nodes_subset.empty()){continue;}

            //Line 11
            std::vector<std::pair<int,int>> clq;
            clq.reserve(nnodes);
            clq.push_back(node);
            //Line 12
            for(std::pair<int,int> n : nodes){
                bool allElementsTrue = true;
                //Line 13
                for(std::pair<int,int> clq_node : clq){
                    if(adj(n.first, clq_node.first) == 0){
                        allElementsTrue = false;
                        break;
                    }
                }
                //Line 14
                if (allElementsTrue) {clq.push_back(n);}
            }

            //Line 15-16
            if (clq.size() > max_val) { // store better clique
                clique.clear();
                clique.reserve(clq.size());
                std::transform(clq.begin(), clq.end(), std::back_inserter(clique), [](const std::pair<int,int>& clq_e){ return clq_e.first; });
                max_val = clq.size();
            }    
        }
    } 
    return 1;
} //find_heuristic_clique

///////////////////////////////////////////////////////////////////////////////////

int clique_corenumber(const Eigen::MatrixXd& adj,
                      std::vector<int>& clique,
                      std::vector<int>& core_numbers,
                      int& core_bound,
                      std::vector<int>& node_colors,
                      int& chromatic_bound) {
    
    const int nnodes = adj.rows(); // number of graph nodes

    // create a PMC graph from the adjacency matrix
    std::vector<int> edges;
    std::vector<long long> vertices;
    vertices.push_back(0);

    size_t total_num_edges = 0; 
    for (size_t i=0; i<nnodes; i++) { 
        for (size_t j=0; j<nnodes; j++) { //NOTE: Was going to try and take advantage of the parallelism in adjacency matrices in unweighted graphs, but there is no way to create the graph for the pmc_graph constructor in a more efficient manner, because their constructor functions basically run this code regardless. 
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

    ////////////////////////////////////////////////////    
    // compute k-cores 
    G.compute_cores(); // compute k-cores
    in.ub = G.get_max_core() + 1; // max clique upper bound as max k-core + 1
    core_bound = in.ub; // output
    std::vector<int> kcores = G.kcore; // k-core values for each node

    // core number ouput
    for(int i = 0; i < nnodes; i++) {
        core_numbers[i] = kcores[i] - 1;
    }
    

    ////////////////////////////////////////////////////    
    // find a maximal clique using core numbers
    find_heuristic_clique(adj, core_numbers, clique);
     
    const int clique_size = clique.size();
   
    ////////////////////////////////////////////////////
    // graph chromatic number (estimate)
    std::vector<int> nodeOrder(nnodes); // vector of node indices sorted by kcores in descending order
    std::iota(nodeOrder.begin(), nodeOrder.end(), 0); // set elements with sequential integers from 0 to n-1
    std::sort(nodeOrder.begin(), nodeOrder.end(), [&](int a, int b) {
        return kcores[a] > kcores[b];
    });

    // std::vector<int> node_colors(nnodes, 0);
    for (int i = 0; i < nnodes; i++) {
        int node = nodeOrder[i];
        std::set<int> neighborColors;

        // find the colors used by neighboring nodes
        for (int j = 0; j < nnodes; j++) {
            if (adj(node,j)==1) {
                neighborColors.insert(node_colors[j]);
            }
        }

        // find the minimum available color
        int color = 1;
        while (neighborColors.count(color) > 0) {
            color++;
        }

        node_colors[node] = color;
    }

    chromatic_bound = *std::max_element(node_colors.begin(), node_colors.end());

    ////////////////////////////////////////////////////
    // return output
    return clique_size;


} //clique_corenumber

} // namespace clipperplus
