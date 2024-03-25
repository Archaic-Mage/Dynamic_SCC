#pragma once
#include <map>
#include <set>

class Edge {
public:
    long int from;
    long int to;
};

class TreeNode {
public:
    long int label;                             // label for each SCC-Tree Node
    TreeNode* parent;
    std::map<Edge, Edge> corresponds_to;
    std::set<TreeNode*> contains;
    long int dept;
};