#include <boost/mpi.hpp>
#include <iostream>
#include "graph.hpp"
#include "debug.hpp"

//broadcasted by the root process
int SCC_LABEL = 0;

enum RET_INT {
  OK, ERR
};

/*** MPI ***/

// Function to get the free label from the root process
int getFreeLabel(const int &who) {
  boost::mpi::communicator world;
  int rank = world.rank();
  int ret = RET_INT::ERR;
  if (rank == 0) {
    world.send(who, 0, SCC_LABEL);
    SCC_LABEL++;
    ret = RET_INT::OK;
  } else if (rank == who) {
    world.recv(0, 0, ret);
  }
  return ret;
}


//Send SCC tree from one processor
void sendSccTree(const int &to, const int &tag, SccTree &sccTree, boost::mpi::communicator &world) {
  world.send(to, tag, sccTree);
}

//Receive SCC tree from another processor
void recieveSccTree(const int &from, const int &tag, SccTree &sccTree, boost::mpi::communicator &world) {
  world.recv(from, tag, sccTree);
}

// Function returns the free unused label in the process
int getLabel() {
  int ret = SCC_LABEL;
  SCC_LABEL++;
  return ret;
}

/*** GRAPH RELATED FUNCS ***/

void TreeNode::condenseFill(std::vector<Edge>& edges, std::unordered_map<long int, long int>& sccs) {
  for(auto& edge: edges) {
    if(sccs[edge.from] != sccs[edge.to]) {
      Edge e;
      e.from = sccs[edge.from];
      e.to = sccs[edge.to];
      corresponds_to[edge] = e;
    }
  }
}

// Function to find the strongly connected components O(n+m)
void findScc(const std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs);

// Function to get the graph from the input (O(m))
void getGraph(std::vector<Edge> &edges, long int &n) {
  int m;
  std::cin >> n >> m;
  for (int i = 0; i < m; i++) {
    Edge edge;
    std::cin >> edge.from >> edge.to;
    edges.push_back(edge);
  }
}

// split a node's edges into two sets (node_in and node_out) creating new graph (O(m))
void splitGraphOnNode(std::vector<Edge> &edges, long int node) {
  for(auto& edge: edges) {
    if(edge.to == node) {
      edge.to = -node;
    }
  }
}

// Divide edges by groups based on SCCs (O(mlog(n)))
void divideEdgesByScc(std::vector<Edge> &edges, std::unordered_map<long int, long int> sccs, std::unordered_map<long int, std::vector<Edge>> &sccEdges) {
  for(auto& edge: edges) {
    if(sccs[edge.from] == sccs[edge.to]) {
      sccEdges[sccs[edge.from]].push_back(edge);
    }
  }
}


/*** Functions to make SCC Tree ***/
void makeSccTreeInternals(std::vector<Edge> &edge, SccTree &sccTree, TreeNode& currentNode) {
  dPrint("Current Node: " + std::to_string(currentNode.label));
  splitGraphOnNode(edge, edge[0].from);
  // actual node -> temp label
  std::unordered_map<long int, long int> sccs;
  findScc(edge, sccs);
  // temp label -> edges in scc
  std::unordered_map<long int, std::vector<Edge>> sccEdges;
  divideEdgesByScc(edge, sccs, sccEdges);
  // temp label -> label
  std::unordered_map<long int, long int> label_mapping;
  for(auto &scc: sccEdges) {
    label_mapping[scc.first] = getLabel();
  }
  //changing sccs mappings to labels
  for(auto& scc: sccs) {
    int temp_label = scc.second;
    if(label_mapping.find(temp_label) != label_mapping.end()) {
      scc.second = label_mapping[temp_label];
    } else {
      scc.second = scc.first;
      label_mapping[temp_label] = scc.first;
    }
  }
  currentNode.condenseFill(edge, sccs);
  dPrint("Condensed Fill");
  for(auto& scc: currentNode.corresponds_to) {
    dEdge(scc.second);
  }
  dPrint("Creating Child Nodes");
  for(auto& scc: label_mapping) {
    //we don't create a node for the split node
    if (scc.second == -edge[0].from)
      continue;
    TreeNode& child = sccTree.nodes[scc.second];
    child.label = scc.second;
    child.parent = currentNode.label;
    child.dept = currentNode.dept + 1;
    currentNode.contains.insert(child.label);
    if(sccEdges.find(scc.first) != sccEdges.end()) {
      makeSccTreeInternals(sccEdges[scc.first], sccTree, child);
    }
  }
}

// Construct the SCC tree
void makeSccTree(std::vector<Edge> &edges, SccTree &sccTree) {
  //Create a TreeNode
  int label = getLabel();
  sccTree.root = label;
  TreeNode& root = sccTree.nodes[label];
  root.label = label;
  root.parent = label;
  root.dept = 0;
  dPrint("SCC Tree Root: " + std::to_string(root.label));
  makeSccTreeInternals(edges, sccTree, root);
}

/*** Functions to update the graph ***/

TreeNode& getLCANode(TreeNode& node1, TreeNode& node2, const SccTree &sccTree) {
  while(node1.dept > node2.dept) {
    node1 = sccTree.nodes.at(node1.parent);
  }
  while(node2.dept > node1.dept) {
    node2 = sccTree.nodes.at(node2.parent);
  }
  while(node1.label != node2.label) {
    node1 = sccTree.nodes.at(node1.parent);
    node2 = sccTree.nodes.at(node2.parent);
  }
  return node1;
}

void deleteEdge(Edge edge, SccTree &sccTree) {
  TreeNode& from = sccTree.nodes.at(edge.from);
  TreeNode& to = sccTree.nodes.at(edge.to);
  TreeNode& lca = getLCANode(from, to, sccTree);
  //delete the edge from the lca
  lca.corresponds_to.erase(edge);
  //do unreachability and reachability checks
}

//query if two nodes are in the same SCC


//TODO: Build a dispatch system for the update/queries
//TODO: Implement the update and query functions
//TODO: Implement SCC (finding) in parallel.

int main(int argc, char *argv[])
{
  // setting up fast input output
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(NULL);
  std::cout.tie(NULL);

  // setting up mpi (boost)
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;

  int rank = world.rank();
  if (rank == 0) {
    std::vector<Edge> edges;
    long int n;
    getGraph(edges, n);
    //set the SCC label to n
    SCC_LABEL = n+1;
    std::unordered_map<long int, long int> sccs;
    findScc(edges, sccs);
    dScc(sccs);
    SccTree sccTree;
    makeSccTree(edges, sccTree);
    sendSccTree(1, 0, sccTree, world);
    // dSccTree(sccTree);
  } 
  else if (rank == 1) {
    SccTree sccTree;
    recieveSccTree(0, 0, sccTree, world);
    dSccTree(sccTree);
  }

  return 0;
}


// change to the adj list type representation for easy transfer
// give set of unique number to generate from them
// bits/std should be made to use individually