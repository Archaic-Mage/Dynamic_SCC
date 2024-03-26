#pragma once
#include <iostream>
#include "graph.hpp"
#include <boost/mpi.hpp>

//printing some information with format strings
void dInfo(const boost::mpi::communicator& world, std::string message) {
    std::cout << "Rank " << world.rank() << ": " << message << std::endl;
}

// print Message
void dPrint(const std::string &message) {
    std::cout << message << std::endl;
}

// print Edge
void dEdge(const Edge &edge) {
    std::cout << "Edge: " << edge.from << "->" << edge.to << std::endl;
}

// Print SCC
void dScc(const std::unordered_map<long int, long int> &sccs) {
    std::unordered_map<long int, std::vector<long int>> scc;
    for (auto &i : sccs) {
        scc[i.second].push_back(i.first);
    }
    std::cout << "Strongly connected components are: " << std::endl;
    for (const auto &component : scc) {
        for (const auto &node : component.second) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
}

void dSccTree(const SccTree &sccTree) {
    std::cout << "SCC Tree: " << std::endl;
    std::function<void(const TreeNode&, int)> dfs = [&](const TreeNode& node, int depth) {
        for (int i = 0; i < depth; i++) {
            std::cout << "  ";
        }
        std::cout << node.label << std::endl;
        for (auto& child : node.contains) {
            dfs(sccTree.nodes.at(child), depth + 1);
        }
    };
    dfs(sccTree.nodes.at(sccTree.root), 0);
}