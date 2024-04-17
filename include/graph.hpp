#pragma once
#include <map>
#include <set>
#include <stack>
#include <vector>
#include <queue>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>

class Edge {
public:
    long int        from;
    long int        to;
    friend          boost::serialization::access;

    Edge() {}
    Edge(int from , int to) : from(from), to(to) {}


    template < class Archive >
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
    long int                            label;
    long int                            parent;
    std::vector<std::pair<Edge, Edge>>  corresponds_to;
    std::unordered_set<long int>        contains;
    long int                            dept;
    friend                              boost::serialization::access;

    template < class Archive >
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
    void updateLabels(std::unordered_map<long int, long int>& new_labels);
    void exposeToParent(TreeNode& parent, std::unordered_set<long int> unreachable);
    void removeEdge(const Edge& edge);

    //copy constructor
    TreeNode(const TreeNode& other) {
        label = other.label;
        parent = other.parent;
        corresponds_to = other.corresponds_to;
        contains = other.contains;
        dept = other.dept;
    }
    //assignment operator
    TreeNode& operator=(const TreeNode& other) {
        label = other.label;
        parent = other.parent;
        corresponds_to = other.corresponds_to;
        contains = other.contains;
        dept = other.dept;
        return *this;
    }
    //default constructor
    TreeNode() {
        label = -1;
        parent = -1;
        dept = 0;
    }
};

class Cache {
public:
    long int        dept;
    long int        label;

    bool operator<(const Cache& other) const {
        return dept < other.dept;
    }
    bool operator==(const Cache& other) const {
        return dept == other.dept && label == other.label;
    }
};

class DecCache : public Cache {
public:
    bool operator() (const Cache& a, const Cache& b) const {
        return a.dept > b.dept;
    }
};

class IncCache : public Cache {
public:
    bool operator() (const Cache& a, const Cache& b) const {
        return a.dept < b.dept;
    }
};

class MaintainSCC {

    long int                                SCC_LABEL = 0;
    const long int                          MOD = 1e10;
    int                                     world_size;
    std::vector<long int>                   roots;
    std::unordered_map<long int, TreeNode>  scc_tree_nodes;
    std::set<Cache, DecCache>               delete_cache;
    std::set<Cache, IncCache>               insert_cache;
    std::unordered_map<long int, int>       which_rank;
    std::unordered_map<long int, long int>  sccs;

    enum MessageType
    {
        SCC_TREE,
        EDGE_DELETE,
        EDGE_INSERT,
        EDGE_QUERY,
        SCC_TRANSFER,
        CLR_DEC_CACHE,
        CLR_INC_CACHE,
        EXIT
    };

    enum STATUS
    {
        DONE_NO_NEW,
        DONE_NEW,
        TRANSFER,
        ERROR
    };

    int getLabel();
    // void findScc(const std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs);
    void traverseNode(int node, std::unordered_map<long int, long int> &new_sccs);
    void splitGraphOnNode(std::vector<Edge> &edges, long int node);
    void divideEdgesByScc(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs, std::unordered_map<long int, std::vector<Edge>> &sccEdges);
    void makeSccTreeInternals(std::vector<Edge> &edge, TreeNode &currentNode);
    void makeSccTree(std::vector<Edge> &edges, std::vector<long int> &nodes);
    void changeSccLabels(std::unordered_map<long int, long int> &sccs, std::unordered_map<long int, std::vector<long int>> &sccNodes);
    void constructMasterNode(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs);
    void updateMasterNode(int type, std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs);
    int checkAndRemoveUnreachable(TreeNode &curr_node);
    int checkAndDetachUnreachable(TreeNode &root_node);
    int getLCANode(int node1, int node2);
    void deleteEdge(Edge edge);
    bool clearDeleteCache(TreeNode &node);
    void deleteEdgeFromMaster(Edge edge);
    void insertEdge(Edge edge);
    void clearInsertCache(TreeNode &node);
    void processMessage();

public:
    bool query(long int v1, long int v2);
    void deleteEdges(std::vector<Edge> &decrement);
    void insertEdges(std::vector<Edge> &increament);
    void endAll();

    MaintainSCC(long int n, std::vector<Edge> &edges);
    ~MaintainSCC();
};