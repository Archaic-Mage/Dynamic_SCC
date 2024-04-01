#pragma once
#include <map>
#include <set>
#include <stack>
#include <vector>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

class Edge {
public:
    long int from;
    long int to;
    Edge() {}
    Edge(int from , int to) : from(from), to(to) {}
    friend boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & from;
        ar & to;
    }
    bool operator<(const Edge& other) const {
        if (from == other.from) {
            return to < other.to;
        }
        return from < other.from;
    }
    bool operator==(const Edge& other) const {
        return from == other.from && to == other.to;
    }
};

class TreeNode {
public:
    long int label;                             // label for each SCC-Tree Node
    long int parent;
    std::vector<std::pair<Edge, Edge>> corresponds_to;
    std::unordered_set<long int> contains;
    long int dept;
    friend boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & label;
        ar & parent;
        ar & corresponds_to;
        ar & contains;
        ar & dept;
    }
    void condenseFill(std::vector<Edge>& edges, std::unordered_map<long int, long int>& sccs);
    void checkUnreachable(std::unordered_set<long int>& unreachable);
    void getNewLabels(std::unordered_set<long int>& unreachable, std::unordered_map<long int, long int>& new_labels);
    void updateLabels(std::unordered_map<long int, long int>& new_labels, long int child_label);
    void exposeToParent(TreeNode& parent, std::unordered_set<long int> unreachable);
    void removeEdge(const Edge& edge);
};

class SccTree {
public:
    long int root;
    std::unordered_map<long int, TreeNode> nodes;
    SccTree() { root = -1; }
    friend boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & root;
        ar & nodes;
    }
};