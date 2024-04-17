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

namespace SCC {

/**
 * @brief Edge class to represent an edge in the graph
 */
class Edge {
public:
    /**
     * @brief from is the source node of the edge
     *
     */
    long int        from;
    /**
     * @brief to is the destination node of the edge
     *
     */
    long int        to;

    /**
     * @brief Construct a new Edge object (DEFAULT CONSTRUCTOR)
     * 
     */
    Edge() {}
    /**
     * @brief Construct a new Edge object
     * 
     * @param from specifies the source node
     * @param to specifies the destination node
     */
    Edge(int from , int to) : from(from), to(to) {}

    /**
     * @brief Serialize the Edge object
     * 
     * Making the boost::serialization friend to access the private members
     * @tparam Archive the type of the archive
     * 
     */
    friend          boost::serialization::access;
    template < class Archive >
    /**
     * @brief This function is used to serialize the Edge object
     * 
     * @param ar The archive object
     * @param version The version of the serialization
     */
    void serialize(Archive & ar, const unsigned int version) {
        ar & from;
        ar & to;
    }

    /**
     * @brief Comparator == for the Edge objects 
     * 
     * This comparator is essential for the set operations to work
     * the unordered_set is used to store the edges in the graph
     * 
     * @param other Compare the current object with this object
     * @return true if both the from and to nodes are equal
     * @return false otherwise
     */
    bool operator==(const Edge& other) const {
        return from == other.from && to == other.to;
    }
};

/**
 * @brief This class is used to represent a single node in the SCC tree. 
 * This contains information about the original graph nodes
 * and structure which is used to dynamically maintain the SCCs
 * 
 */
class TreeNode {
public:
    /**
     * @brief this would represent the label of the node in the SCC tree
     * the label is unique for each node in the SCC tree
     * 
     */
    long int                            label;
    /**
     * @brief this represents the parent of the current node in the SCC tree
     *  and is used to traverse the tree upwords.
     * 
     */
    long int                            parent;
    /**
     * @brief stores the actual edge to labelled edge pair,
     * nodes forming a scc are condesed to one labelled node, and so are the edges,
     * the corresponds_to is used to store the actual edges to labelled edges
     * relationship.
     * 
     */
    std::vector<std::pair<Edge, Edge>>  corresponds_to;
    /**
     * @brief contains the nodes which are part of the current node in the SCC tree.
     * They store the label of the child nodes and is used to travese the tree downwards
     * 
     */
    std::unordered_set<long int>        contains;
    /**
     * @brief stores the depth of the current node in the SCC tree and is
     * used to find the LCA of two nodes (Lowest Common Ancestor), 
     * important for update operations.
     * 
     */
    long int                            dept;
    /**
     * @brief The following function is used to serialize the TreeNode object,
     * this is essential for the boost::serialization to work. The TreeNodes are 
     * transferred during update queries.
     * 
     */
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

}