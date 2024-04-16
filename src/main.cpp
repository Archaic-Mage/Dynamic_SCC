#include "graph.hpp"
#include "debug.hpp"
#include <queue>
#include <chrono>
#include <fstream>
#include <climits>
#include <iostream>
#include <boost/mpi.hpp>

long int SCC_LABEL = 0;
const long int MOD = 1e10;
int world_size;
std::vector<long int> roots;
std::unordered_map<long int, TreeNode> scc_tree_nodes;
std::set<Cache, DecCache> deleteCache;
std::set<Cache, IncCache> insertCache;

/*** MPI ***/

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

// Function returns the free unused label in the process
int getLabel()
{
  int ret = SCC_LABEL;
  SCC_LABEL+=world_size;
  return ret;
}

/*** GRAPH RELATED FUNCS ***/

// Function to find the strongly connected components O(n+m)
void findScc(const std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs);

void TreeNode::removeEdge(const Edge &edge)
{
  for (auto it = this->corresponds_to.begin(); it != this->corresponds_to.end(); it++)
  {
    if (it->first == edge)
    {
      this->corresponds_to.erase(it);
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
  // get the nodes in the current scc
  std::unordered_map<long int, long int> sccs;
  for(auto &node: contains) {
    sccs[node] = node;
  }
  findScc(edges, sccs);
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

//Function to get the list of updates
void getUpdates(int &updates, std::vector<Edge> &decrement, std::vector<Edge> &increament) {
  std::cin >> updates;
  for(int i = 0; i<updates; i++) {
    int type; std::cin >> type;
    long int from, to;
    std::cin >> from >> to;
    from++, to++;
    if(type == 0) decrement.emplace_back(from, to);
    else increament.emplace_back(from, to);
  }
}

void getQueries(int &q, std::vector<std::pair<long int,long int>> &queries) {
  std::cin >> q;
  for(int i = 0; i<q; i++) {
    int v1, v2; std::cin >> v1 >> v2;
    v1++, v2++;
    queries.emplace_back(v1, v2);
  }
}

// traverse the Node and all the containing nodes
void traverseNode(int node, std::vector<std::pair<long int, long int>> &new_sccs) {
  std::queue<long int> q;
  q.push(node);
  while(!q.empty()) {
    long int curr = q.front();
    q.pop();
    new_sccs.emplace_back(std::make_pair(curr, node));
    for(auto &child: scc_tree_nodes[curr].contains) {
      if (child != curr) {
        q.push(child);
      }
    }
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
void makeSccTreeInternals(std::vector<Edge> &edge, TreeNode &currentNode)
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
    TreeNode &child = scc_tree_nodes[scc.second];
    child.label = scc.second;
    child.parent = currentNode.label;
    child.dept = currentNode.dept + 1;
    currentNode.contains.insert(child.label);
    if (sccEdges.find(scc.first) != sccEdges.end())
    {
      makeSccTreeInternals(sccEdges[scc.first], child);
    }
  }
}

// Construct the SCC tree
void makeSccTree(std::vector<Edge> &edges, std::vector<long int> &nodes)
{
  //case when only one node is present in the scc
  if (nodes.size() == 1) {
    TreeNode &root = scc_tree_nodes[nodes[0]];
    root.label = nodes[0];
    root.parent = nodes[0];
    root.dept = 0;
    return;
  }
  // Create a TreeNode
  int label = getLabel();
  roots.push_back(label);
  TreeNode &root = scc_tree_nodes[label];
  root.label = label;
  root.parent = label;
  root.dept = 0;
  // dPrint("SCC Tree Root: " + std::to_string(root.label));
  makeSccTreeInternals(edges, root);
}

void changeSccLabels(std::unordered_map<long int, long int> &sccs, std::unordered_map<long int, std::vector<long int>> &sccNodes) {
  std::unordered_map<long int, long int> new_labels;
  for(auto &scc: sccs) {
    if(new_labels.find(scc.second) == new_labels.end())
      if(sccNodes[scc.second].size() == 1)
        new_labels[scc.second] = sccNodes[scc.second][0];
      else
        new_labels[scc.second] = getLabel();
  }
  for(auto &scc: sccs) {
    scc.second = new_labels[scc.second];
  }
}

// Construct Master Node
void constructMasterNode(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs)
{
  // Create a TreeNode
  int label = 0;
  roots.push_back(label);
  TreeNode &root = scc_tree_nodes[label];
  root.label = label;
  root.parent = label;
  root.dept = 0;
  for(const auto &edge: edges) {
    if(sccs[edge.from] != sccs[edge.to]) {
      Edge e;
      e.from = sccs[edge.from];
      e.to = sccs[edge.to];
      root.corresponds_to.emplace_back(std::make_pair(edge, e));
    }
  }
  for(auto &scc: sccs) {
    root.contains.insert(scc.second);
  }
}

//update master node 0-new scc creation
void updateMasterNode(int type, std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs) {
  if(type == 0) {
    for (auto &correspond_to: scc_tree_nodes[0].corresponds_to) {
      Edge &actual_edge = correspond_to.first;
      Edge &label_edge = correspond_to.second;
      long int from = actual_edge.from;
      if (sccs.find(from) != sccs.end()) {
        label_edge.from = sccs[from];
      }
      long int to = actual_edge.to;
      if (sccs.find(to) != sccs.end()) {
        label_edge.to = sccs[to];
      }
    }
    for(const auto &edge: edges) {
      if(sccs[edge.from] != sccs[edge.to]) {
        Edge e;
        e.from = sccs[edge.from];
        e.to = sccs[edge.to];
        scc_tree_nodes[0].corresponds_to.emplace_back(std::make_pair(edge, e));
      }
    }
    for(auto &scc: sccs) {
      scc_tree_nodes[0].contains.insert(scc.second);
    }
  }
}

int checkAndRemoveUnreachable(TreeNode &curr_node) {
  // do unreachability and reachability checks
  std::unordered_set<long int> unreachable;
  curr_node.checkUnreachable(unreachable);
  //if no unreachable nodes then return
  if (unreachable.empty()) {
    return 0;
  }
  for (auto &node : unreachable)
  {
    dPrint("Unreachable: " + std::to_string(node));
  }
  //remove all unreachable nodes from the current node
  for(auto &unreachable_node: unreachable) {
    curr_node.contains.erase(unreachable_node);
  }
  //add the unreachable nodes to the parent
  TreeNode &parent = scc_tree_nodes.at(curr_node.parent);
  //change parent label of unreachable nodes
  for(auto &unreachable_node: unreachable) {
    scc_tree_nodes[unreachable_node].parent = parent.label;
    parent.contains.insert(unreachable_node);
  }
  //new_labeling
  std::unordered_map<long int, long int> new_labels;
  //change label of node if only one node present O(1)
  int curr_node_label = curr_node.label;
  if (curr_node.contains.size() == 1) {
    long int new_label = *(curr_node.contains.begin());
    parent.contains.erase(curr_node.label);
    curr_node.label = new_label;
    parent.contains.erase(curr_node.label);
    parent.contains.insert(new_label);
    scc_tree_nodes[new_label] = curr_node;
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
  return 1;
}
//TODO: check and detach unreachable nodes
int checkAndDetachUnreachable(TreeNode &root_node, boost::mpi::communicator &world) {
  std::unordered_set<long int> unreachable;
  root_node.checkUnreachable(unreachable);
  if (unreachable.empty()) {
    return 0;
  }
  for (auto &node : unreachable)
  {
    dPrint("Unreachable: " + std::to_string(node));
  }
  for(auto &unreachable_node: unreachable) {
    root_node.contains.erase(unreachable_node);
  }
  TreeNode temp;
  root_node.exposeToParent(temp, unreachable);
  std::vector<Edge> edges;
  for(auto &correspond_to: temp.corresponds_to) {
    edges.emplace_back(correspond_to.first);
  }
  std::vector<std::pair<long int,long int>> new_sccs;
  for(auto &unreachable_node: unreachable) {
    roots.push_back(unreachable_node);
    traverseNode(unreachable_node, new_sccs);
  }
  //transfer to master node for creating new scc
  world.send(0, 0, STATUS::DONE_NEW);
  world.send(0, 1, edges);
  world.send(0, 2, new_sccs);

  return 0;
}

/*** Functions to update the graph ***/

int getLCANode(int node1, int node2)
{
  while (node1 != node2)
  {
    if (scc_tree_nodes[node1].dept > scc_tree_nodes[node2].dept)
    {
      node1 = scc_tree_nodes[node1].parent;
    }
    else
    {
      node2 = scc_tree_nodes[node2].parent;
    }
  }
  return node1;
}

//DELETE EDGE

void deleteEdge(Edge edge) {
  TreeNode &lca = scc_tree_nodes.at(getLCANode(edge.from, edge.to));
  // delete the edge from the lca
  lca.removeEdge(edge);
  //add the lca to the delete cache
  deleteCache.insert({lca.dept, lca.label});
}

bool clearDeleteCache(TreeNode &node, boost::mpi::communicator &world)
{
  if(node.dept != 0) return checkAndRemoveUnreachable(node);
  return checkAndDetachUnreachable(node, world);
}

void deleteEdgeFromMaster(Edge edge) {
  TreeNode &root = scc_tree_nodes[0];
  root.removeEdge(edge);
}

//INSERT EDGE
void insertEdge(Edge edge) {
  TreeNode &lca = scc_tree_nodes.at(getLCANode(edge.from, edge.to));
  int label_from = edge.from;
  while(lca.label != scc_tree_nodes[label_from].parent) {
    label_from = scc_tree_nodes[label_from].parent;
  }
  int label_to = edge.to;
  while(lca.label != scc_tree_nodes[label_to].parent) {
    label_to = scc_tree_nodes[label_to].parent;
  }
  Edge e(label_from, label_to);
  // add the edge to the lca
  lca.corresponds_to.emplace_back(std::make_pair(edge, e));
  //add the lca to the insert cache
  insertCache.insert({lca.dept, lca.label});
}

void clearInsertCache(TreeNode &node) {
  std::vector<Edge> edges;
  for(auto &correspond_to: node.corresponds_to) {
    edges.emplace_back(correspond_to.second);
  }
  std::unordered_map<long int, long int> sccs;
  findScc(edges, sccs);

  std::unordered_map<long int, std::vector<long int>> sccNodes;
  for(auto &scc: sccs) {
    sccNodes[scc.second].push_back(scc.first);
  }

  bool to_end = true;
  for(auto &scc_nodes: sccNodes) {
    if(scc_nodes.second.size() > 1) {
      to_end = false;
    }
  }
  if(to_end) return;

  changeSccLabels(sccs, sccNodes);

  sccNodes.clear();
  for(auto &scc: sccs) {
    sccNodes[scc.second].push_back(scc.first);
  }

  std::vector<long int> new_tree_nodes;
  for(auto &scc: sccs) {
    if(scc.second < 0) continue;
    if(scc_tree_nodes.find(scc.second) == scc_tree_nodes.end()) {
      TreeNode &child = scc_tree_nodes[scc.second];
      child.label = scc.second;
      child.parent = node.label;
      child.dept = node.dept + 1;
      node.contains.insert(sccNodes[scc.second].begin(), sccNodes[scc.second].end());
      new_tree_nodes.push_back(child.label);
    }
  }

  std::vector<std::pair<Edge, Edge>> temp;
  while(!node.corresponds_to.empty()) {
    std::pair<Edge, Edge> &correspond_to = node.corresponds_to.back();
    node.corresponds_to.pop_back();
    if(sccs[correspond_to.second.from] != sccs[correspond_to.second.to]) {
      temp.emplace_back(correspond_to);
    } else {
      scc_tree_nodes[sccs[correspond_to.second.from]].corresponds_to.emplace_back(correspond_to);
    }
  }
  node.corresponds_to = temp;

  //changing connections in node
  for(auto &correspond_to: node.corresponds_to) {
    Edge &label_edge = correspond_to.second;
    Edge &actual_edge = correspond_to.first;
    label_edge.from = sccs[label_edge.from];
    label_edge.to = sccs[label_edge.to];
  }

  //update the contains set
  node.contains.clear();
  for(auto &scc: sccs) {
    if(scc.second < 0) continue;
    node.contains.insert(scc.second);
  }
  dTreeNode(node);
  for(auto &new_tree_node: new_tree_nodes) {
    TreeNode &child = scc_tree_nodes[new_tree_node];
    dTreeNode(child);
    // insertCache.insert({child.dept, child.label});
  }
}

// query if two nodes are in the same SCC
bool query(long int v1, long int v2, std::unordered_map<long int, long int> sccs) {
  return sccs[v1] == sccs[v2];
}

// TODO: Implement the update functions

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
  world_size = world.size();

  if (rank == 0)
  {
    std::unordered_map<long int, int> which_rank;
    std::vector<Edge> edges;
    long int n;

    getGraph(edges, n);

    std::vector<Edge> decrement;
    std::vector<Edge> increament;
    int updates;

    getUpdates(updates, decrement, increament);

    int q;
    std::vector<std::pair<long int, long int>> queries;

    getQueries(q, queries);

    auto time_con = std::chrono::high_resolution_clock::now();

    SCC_LABEL = n + 1;
    std::unordered_map<long int, long int> sccs;

    for(long int i = 1; i<=n; i++) {
      sccs[i] = i;
    }

    findScc(edges, sccs);

    std::unordered_map<long int, std::vector<Edge>> sccEdges;
    divideEdgesByScc(edges, sccs, sccEdges);

    std::unordered_map<long int, std::vector<long int>> sccNodes;
    for(auto &scc: sccs) {
      sccNodes[scc.second].push_back(scc.first);
      which_rank[scc.first] = scc.second % (world.size()-1) + 1;
    }

    for(auto &scc_edges: sccEdges) {
      std::vector<long int> &nodes = sccNodes[scc_edges.first];
      int to_rank = which_rank[nodes[0]];
      world.send(to_rank, 0, MessageType::SCC_TREE);
      world.send(to_rank, 1, SCC_LABEL+to_rank);
      world.send(to_rank, 2, nodes);
      world.send(to_rank, 3, scc_edges.second);
    }

    changeSccLabels(sccs, sccNodes);
    constructMasterNode(edges, sccs);

    auto time_con_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_con = time_con_end - time_con;
    std::cout << "Time taken for construction: " << elapsed_con.count() << "s" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    //updates
    //deleting edges
    for(auto &edge: decrement) {
      if(sccs[edge.from] == sccs[edge.to]) {
        int to_rank = which_rank[edge.from];
        world.send(to_rank, 0, MessageType::EDGE_DELETE);
        world.send(to_rank, 1, edge);
      } else {
        deleteEdgeFromMaster(edge);
      }
    }
    for(int i = 0; i<world.size(); i++) {
      world.send(i, 0, MessageType::CLR_DEC_CACHE);
    }

    //waiting for signal of completion
    boost::mpi::request req[world.size()];
    STATUS status[world.size()];
    for(int i = 1; i<world.size(); i++) {
      req[i] = world.irecv(i, 0, status[i]);
    }
    boost::mpi::wait_all(req+1, req+world.size());

    for(int i = 1; i<world.size(); i++) {
      if(status[i] == STATUS::DONE_NEW) {
        std::vector<Edge> edges;
        std::vector<std::pair<long int,long int>> node_list;
        // std::vector<std::pair<long int, int>> transfer_list;
        world.recv(i, 1, edges);
        world.recv(i, 2, node_list);
        for(auto edge: edges) {
          dEdge(edge);
        }
        for(auto &node: node_list) {
          sccs[node.first] = node.second;
          int to_rank = node.second % (world.size()-1) + 1;
          which_rank[node.first] = to_rank;
          // transfer_list.emplace_back(std::make_pair(node.first, to_rank));
        }
        updateMasterNode(0, edges, sccs);
        // world.send(i, 0, STATUS::TRANSFER);
        // world.send(i, 1, transfer_list);
        dTreeNode(scc_tree_nodes[0]);
      }
    }

    //inserting edges
    for(auto &edge: increament) {
      if(sccs[edge.from] == sccs[edge.to]) {
        int to_rank = which_rank[edge.from];
        world.send(to_rank, 0, MessageType::EDGE_INSERT);
        world.send(to_rank, 1, edge);
      } else {
        //insert edge in master node
      }
    }

    for(int i = 1; i<world.size(); i++) {
      world.send(i, 0, MessageType::CLR_INC_CACHE);
    }

    for(int i = 1; i<world.size(); i++) {
      req[i] = world.irecv(i, 0, status[i]);
    }
    boost::mpi::wait_all(req+1, req+world.size());

    //answering query
    std::ofstream out("output.txt");
    for(int i = 0; i<q; i++) {
      long int v1 = queries[i].first;
      long int v2 = queries[i].second;
      if(query(v1, v2, sccs)) {
        out << "YES" << std::endl;
      } else {
        out << "NO" << std::endl;
      }
    }
    out.close();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken for updates/queries: " << elapsed.count() << "s" << std::endl;

    for(int i = 1; i<world.size(); i++) {
      world.send(i, 0, MessageType::EXIT);
    }
  }
  else
  {
    MessageType type;
    while(true) {
      world.recv(0, 0, type);
      if(type == MessageType::EXIT) 
        break;
      else if (type == MessageType::SCC_TREE) 
      {
        std::vector<long int> nodes;
        std::vector<Edge> edgeList;
        world.recv(0, 1, SCC_LABEL);
        world.recv(0, 2, nodes);
        world.recv(0, 3, edgeList);
        makeSccTree(edgeList, nodes);
        for(auto &root: roots) {
          dSccTree(root, scc_tree_nodes);
        }
      }
      else if (type == MessageType::EDGE_DELETE) 
      {
        Edge edge;
        world.recv(0, 1, edge);
        deleteEdge(edge);
        std::cout << "Deleted Edge: " << edge.from << " " << edge.to << std::endl;
      } else if (type == MessageType::CLR_DEC_CACHE) 
      {
        while(!deleteCache.empty()) {
          dInfo(world, "Deleting Cache");
          Cache c = *deleteCache.begin();
          deleteCache.erase(deleteCache.begin());
          TreeNode &node = scc_tree_nodes.at(c.label);
          dTreeNode(node);
          if (clearDeleteCache(node, world)) {
            deleteCache.insert({node.dept-1, node.parent});
          }
        }
        world.send(0, 0, STATUS::DONE_NO_NEW);
      } else if (type == MessageType::EDGE_INSERT) 
      {
        Edge edge;
        world.recv(0, 1, edge);
        std::cout << "Inserted Edge: " << edge.from << " " << edge.to << std::endl;
        insertEdge(edge);
      } else if(type == MessageType::CLR_INC_CACHE) {
        while(!insertCache.empty()) {
          dInfo(world, "Inserting Cache");
          Cache c = *insertCache.begin();
          insertCache.erase(insertCache.begin());
          TreeNode &node = scc_tree_nodes.at(c.label);
          clearInsertCache(node);
        }
        world.send(0, 0, STATUS::DONE_NO_NEW);
      }
    }
  }

  return 0;
}

// give set of unique number to generate from them