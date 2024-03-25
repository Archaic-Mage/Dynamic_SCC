#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <iostream>
#include <vector>
#include <map>
#include "graph.hpp"
#include "debug.hpp"

// Function to find the strongly connected components
void findScc(const std::vector<Edge> &edges, std::vector<std::vector<long int>> &sccs);

// Function to get the graph from the input
void getGraph(std::vector<Edge> &edges, long int &n) {
  int m;
  std::cin >> n >> m;
  for (int i = 0; i < m; i++) {
    Edge edge;
    std::cin >> edge.from >> edge.to;
    edges.push_back(edge);
  }
}

// split a node's edges into two sets (node_in and node_out) creating new graph
void splitGraphOnNode(std::vector<Edge> &edges, long int node) {
  for(auto& edge: edges) {
    if(edge.to == node) {
      edge.to = -node;
    }
  }
}

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
    std::vector<std::vector<long int>> sccs;
    std::vector<std::vector<long int>> sccs2;
    findScc(edges, sccs);
    dScc(sccs);
    splitGraphOnNode(edges, 2);
    for(auto& edge: edges) {
      dEdge(edge);
    }
    findScc(edges, sccs2);
    dScc(sccs2);
  }

  return 0;
}