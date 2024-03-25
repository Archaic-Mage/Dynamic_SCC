#pragma once
#include <iostream>
#include <vector>
#include <string>
#include "graph.hpp"
#include <boost/mpi.hpp>

//printing some information with format strings
void dInfo(const boost::mpi::communicator& world, std::string message) {
    std::cout << "Rank " << world.rank() << ": " << message << std::endl;
}

// print Edge
void dEdge(const Edge &edge) {
    std::cout << "Edge: " << edge.from << "->" << edge.to << std::endl;
}

// Print SCC
void dScc(const std::vector<std::vector<long int>> &scc) {
    std::cout << "Strongly connected components are: " << std::endl;
    for (const auto &component : scc) {
        for (const auto &node : component) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
}