#include <boost/mpi.hpp>
#include <iostream>
#include <chrono>
#include "graph.hpp"
#include "debug.hpp"

// broadcasted by the root process
int SCC_LABEL = 0;

/*** MPI ***/

// Function returns the free unused label in the process
int getLabel()
{
  int ret = SCC_LABEL;
  SCC_LABEL++;
  return ret;
}

/*** GRAPH RELATED FUNCS ***/

// Function to find the strongly connected components O(n+m)
void findScc(const std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs);

void TreeNode::removeEdge(const Edge &edge)
{
  for (auto it = corresponds_to.begin(); it != corresponds_to.end(); it++)
  {
    if (it->first == edge)
    {
      corresponds_to.erase(it);
      break;
    }
  }
}

void TreeNode::condenseFill(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs)
{
  for (auto &edge : edges)
  {
    if (sccs[edge.from] != sccs[edge.to])
    {
      Edge e;
      e.from = sccs[edge.from];
      e.to = sccs[edge.to];
      //removing the split node
      if(edge.to < 0) 
        edge.to = -edge.to;
      corresponds_to.emplace_back(std::make_pair(edge, e));
    }
  }
}

void TreeNode::checkUnreachable(std::unordered_set<long int> &unreachable)
{
  //get edge list 
  std::vector<Edge> edges;
  long int split_on = -1;
  for (auto edge: corresponds_to) {
    if (edge.second.to < 0) {
      split_on = -edge.second.to;
      edge.second.to = split_on;
    }
    edges.emplace_back(edge.second);
  }
  for(auto &edge: edges) {
    dEdge(edge);
  }
  // get the nodes in the current scc
  std::unordered_map<long int, long int> sccs;
  for(auto &node: contains) {
    sccs[node] = node;
  }
  findScc(edges, sccs);
  dScc(sccs);
  int split_on_label = sccs[split_on];
  for(const auto &node: contains) {
    if (sccs[node] != split_on_label) {
      unreachable.insert(node);
    }
  }
}

void TreeNode::getNewLabels(std::unordered_set<long int> &unreachable, std::unordered_map<long int, long int> &new_labels)
{
  for(auto &correspond_to : corresponds_to) {
    const Edge &actual_edge = correspond_to.first;
    const Edge &label_edge = correspond_to.second;
    if (unreachable.find(label_edge.from) != unreachable.end())
      new_labels[actual_edge.from] = label_edge.from;
    else 
      new_labels[actual_edge.from] = label;
    if (unreachable.find(label_edge.to) != unreachable.end())
      new_labels[actual_edge.to] = label_edge.to;
    else
      new_labels[actual_edge.to] = label;
  }
}

void TreeNode::updateLabels(std::unordered_map<long int, long int> &new_labels, long int child_label) {
  for(auto &correspond_to : corresponds_to) {
    Edge &label_edge = correspond_to.second;
    Edge &actual_edge = correspond_to.first;
    if (label_edge.from == child_label) {
      label_edge.from = new_labels[actual_edge.from];
    }
    if(label_edge.to == child_label) {
      label_edge.to = new_labels[actual_edge.to];
    }
  }
}

void TreeNode::exposeToParent(TreeNode &parent, std::unordered_set<long int> unreachable) {
  std::vector<std::pair<Edge, Edge>> temp;
  while(!corresponds_to.empty()) {
    std::pair<Edge, Edge> &correspond_to = corresponds_to.back();
    corresponds_to.pop_back();
    if (unreachable.find(correspond_to.second.from) != unreachable.end() || unreachable.find(correspond_to.second.to) != unreachable.end()) {
      if(unreachable.find(correspond_to.second.from) == unreachable.end()) {
        correspond_to.second.from = label;
      }
      if(unreachable.find(correspond_to.second.to) == unreachable.end()) {
        correspond_to.second.to = label;
      }
      parent.corresponds_to.push_back(correspond_to);
    } else {
      temp.push_back(correspond_to);
    }
  }
  corresponds_to = temp;
}

// Function to get the graph from the input (O(m))
void getGraph(std::vector<Edge> &edges, long int &n)
{
  std::cout << "started taking input" << std::endl;
  long int m;
  std::cin >> n >> m;
  for (long unsigned int i = 0; i < m; i++)
  {
    long int from, to;
    std::cin >> from >> to;
    from++, to++;
    edges.emplace_back(from, to);
  }
}

// split a node's edges into two sets (node_in and node_out) creating new graph (O(m))
void splitGraphOnNode(std::vector<Edge> &edges, long int node)
{
  for (auto &edge : edges)
  {
    if (edge.to == node)
    {
      edge.to = -node;
    }
  }
}

// Divide edges by groups based on SCCs (O(mlog(n)))
void divideEdgesByScc(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs, std::unordered_map<long int, std::vector<Edge>> &sccEdges)
{
  for (auto &edge : edges)
  {
    if (sccs[edge.from] == sccs[edge.to])
    {
      sccEdges[sccs[edge.from]].push_back(edge);
    }
  }
}

/*** Functions to make SCC Tree ***/
void makeSccTreeInternals(std::vector<Edge> &edge, SccTree &sccTree, TreeNode &currentNode)
{
  // dPrint("Current Node: " + std::to_string(currentNode.label));
  splitGraphOnNode(edge, edge[0].from);
  // actual node -> temp label
  std::unordered_map<long int, long int> sccs;
  findScc(edge, sccs);
  // temp label -> edges in scc
  std::unordered_map<long int, std::vector<Edge>> sccEdges;
  divideEdgesByScc(edge, sccs, sccEdges);
  // temp label -> label
  std::unordered_map<long int, long int> label_mapping;
  for (auto &scc : sccEdges)
  {
    label_mapping[scc.first] = getLabel();
  }
  // changing sccs mappings to labels
  for (auto &scc : sccs)
  {
    int temp_label = scc.second;
    if (label_mapping.find(temp_label) != label_mapping.end())
    {
      scc.second = label_mapping[temp_label];
    }
    else
    {
      scc.second = scc.first;
      label_mapping[temp_label] = scc.first;
    }
  }
  currentNode.condenseFill(edge, sccs);
  // dPrint("Condensed Fill");
  // for(auto& scc: currentNode.corresponds_to) {
  //   dEdge(scc.second);
  // }
  // dPrint("Creating Child Nodes");
  for (auto &scc : label_mapping)
  {
    // we don't create a node for the split node
    if (scc.second == -edge[0].from)
      continue;
    TreeNode &child = sccTree.nodes[scc.second];
    child.label = scc.second;
    child.parent = currentNode.label;
    child.dept = currentNode.dept + 1;
    currentNode.contains.insert(child.label);
    if (sccEdges.find(scc.first) != sccEdges.end())
    {
      makeSccTreeInternals(sccEdges[scc.first], sccTree, child);
    }
  }
}

// Construct the SCC tree
void makeSccTree(std::vector<Edge> &edges, SccTree &sccTree)
{
  // Create a TreeNode
  int label = getLabel();
  sccTree.root = label;
  TreeNode &root = sccTree.nodes[label];
  root.label = label;
  root.parent = label;
  root.dept = 0;
  // dPrint("SCC Tree Root: " + std::to_string(root.label));
  makeSccTreeInternals(edges, sccTree, root);
}

void checkAndRemoveUnreachable(TreeNode &curr_node, SccTree &sccTree) {
  // do unreachability and reachability checks
  std::unordered_set<long int> unreachable;
  curr_node.checkUnreachable(unreachable);
  for (auto &node : unreachable)
  {
    dPrint("Unreachable: " + std::to_string(node));
  }
  //remove all unreachable nodes from the current node
  for(auto &unreachable_node: unreachable) {
    curr_node.contains.erase(unreachable_node);
  }
  //add the unreachable nodes to the parent
  TreeNode &parent = sccTree.nodes.at(curr_node.parent);
  for(auto &unreachable_node: unreachable) {
    parent.contains.insert(unreachable_node);
  }
  //new_labeling
  std::unordered_map<long int, long int> new_labels;
  //change label of node if only one node present O(1)
  int curr_node_label = curr_node.label;
  if (curr_node.contains.size() == 1) {
    long int new_label = *(curr_node.contains.begin());
    curr_node.label = new_label;
    parent.contains.erase(curr_node.label);
    parent.contains.insert(new_label);
    sccTree.nodes[new_label] = curr_node;
  }
  //fill the new_labels with new mappings that would go in the parent
  curr_node.getNewLabels(unreachable, new_labels);
  for(auto &new_label: new_labels) {
    dPrint("New Label: " + std::to_string(new_label.first) + " -> " + std::to_string(new_label.second));
  }
  //update the parent with the new labels
  parent.updateLabels(new_labels, curr_node_label);


  //exposing internal structure to parent node
  curr_node.exposeToParent(parent, unreachable);
  dTreeNode(curr_node);
  dTreeNode(parent);
}

/*** Functions to update the graph ***/

TreeNode &getLCANode(TreeNode &node1, TreeNode &node2, const SccTree &sccTree)
{
  while (node1.dept > node2.dept)
  {
    node1 = sccTree.nodes.at(node1.parent);
  }
  while (node2.dept > node1.dept)
  {
    node2 = sccTree.nodes.at(node2.parent);
  }
  while (node1.label != node2.label)
  {
    node1 = sccTree.nodes.at(node1.parent);
    node2 = sccTree.nodes.at(node2.parent);
  }
  return node1;
}

void deleteEdge(Edge edge, SccTree &sccTree)
{
  TreeNode &from = sccTree.nodes.at(edge.from);
  TreeNode &to = sccTree.nodes.at(edge.to);
  TreeNode &lca = getLCANode(from, to, sccTree);
  // delete the edge from the lca
  lca.removeEdge(edge);
  checkAndRemoveUnreachable(lca, sccTree);
}

// query if two nodes are in the same SCC

// TODO: Build a dispatch system for the update/queries
// TODO: Implement the update and query functions
// TODO: Implement SCC (finding) in parallel.

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
  // int rank = 0;
  if (rank == 0)
  {
    std::vector<Edge> edges;
    long int n;
    getGraph(edges, n);
    std::cout << "Graph taken" << std::endl;
    // set the SCC label to n
    SCC_LABEL = n + 1;
    std::unordered_map<long int, long int> sccs;
    // start timer
    auto start = std::chrono::high_resolution_clock::now();
    findScc(edges, sccs);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = stop - start;
    double time_taken = (double)duration.count() / 1e9;
    printf("Time taken for SCC: %.9lf s\n", time_taken);
    SccTree sccTree;
    start = std::chrono::high_resolution_clock::now();
    makeSccTree(edges, sccTree);
    stop = std::chrono::high_resolution_clock::now();
    duration = stop - start;
    time_taken = (double)duration.count() / 1e9;
    printf("Time taken for SCC Tree: %.9lf s\n", time_taken);

    dSccTree(sccTree);

    deleteEdge(Edge(4, 5), sccTree);

    // sendSccTree(1, 0, sccTree, world);
  }
  else if (rank == 1)
  {
    // SccTree sccTree;
    // recieveSccTree(0, 0, sccTree, world);
    // dSccTree(sccTree);
  }

  return 0;
}

// give set of unique number to generate from them