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
};

class TreeNode {
public:
    long int label;                             // label for each SCC-Tree Node
    long int parent;
    std::map<Edge, Edge> corresponds_to;
    std::set<long int> contains;
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
    void checkUnreachable(std::vector<long int>& unreachable) {
        long int splitOn = -1;
        for(auto& edge: corresponds_to) {
            if(edge.second.to < 0) {
                splitOn = -edge.first.to;
                edge.second.to = -edge.second.to;
            }
        }

    }
};

class SccTree {
public:
    long int root;
    std::map<long int, TreeNode> nodes;
    SccTree() { root = -1; }
    friend boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & root;
        ar & nodes;
    }
};