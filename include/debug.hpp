#pragma once
#include <iostream>
#include "graph.hpp"
#include <boost/mpi.hpp>

// Debug flag
#ifndef NDEBUG
#define DEBUG true
#endif

#ifdef NDEBUG
#define DEBUG false
#endif



using namespace SCC;
//printing some information with format strings
void dInfo(const boost::mpi::communicator& world, std::string message) {
    if(!DEBUG) return;
    std::cout << "Rank " << world.rank() << ": " << message << std::endl;
}

// print Message
void dPrint(const std::string &message) {
    if(!DEBUG) return;
    std::cout << message << std::endl;
}

// print Edge
void dEdge(const Edge &edge) {
    if(!DEBUG) return;
    std::cout << "Edge: " << edge.from << "->" << edge.to << std::endl;
}

// Print SCC
void dScc(const std::unordered_map<long int, long int> &sccs) {
    if(!DEBUG) return;
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

void dTreeNode(const TreeNode &node) {
    if(!DEBUG) return;
    std::cout << "TreeNode: " << std::endl;
    std::cout << "Label: " << node.label << std::endl;
    std::cout << "Parent: " << node.parent << std::endl;
    std::cout << "Dept: " << node.dept << std::endl;
    std::cout << "Contains: ";
    for (auto &i : node.contains) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    std::cout << "Corresponds to: ";
    int max_threads = omp_get_max_threads();
    for(int i = 0; i<max_threads; i++)
        for (auto &i : node.corresponds_to[i]) {
            std::cout << "(" << i.first.from << "->" << i.first.to << ", " << i.second.from << "->" << i.second.to << ") ";
        }
    std::cout << std::endl;
}

void dSccTree(const long int &root, const std::unordered_map<long int, TreeNode> &nodes) {
    if(!DEBUG) return;
    std::cout << "SCC Tree: " << std::endl;
    std::function<void(const TreeNode&, int)> dfs = [&](const TreeNode& node, int depth) {
        for (int i = 0; i < depth; i++) {
            std::cout << "  ";
        }
        std::cout << node.label << std::endl;
        for (auto& child : node.contains) {
            dfs(nodes.at(child), depth + 1);
        }
    };
    dfs(nodes.at(root), 0);
}