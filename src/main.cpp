#include "graph.hpp"
#include <chrono>
#include <fstream>
#include <climits>
#include <iostream>
#include <boost/mpi.hpp>

// Function to get the graph from the input (O(m))
void getGraph(std::vector<SCC::Edge> &edges, long int &n)
{
  long int m;
  std::cin >> n >> m;
  for (long unsigned int i = 0; i < m; i++)
  {
    long int from, to;
    std::cin >> from >> to;
    from++, to++;
    if(from == to)
      continue;
    edges.emplace_back(from, to);
  }
}

namespace mt = boost::mpi::threading;

int main(int argc, char *argv[])
{
  // setting up fast input output
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(NULL);
  std::cout.tie(NULL);

  // setting up mpi (boost)
  boost::mpi::environment env(argc, argv, mt::funneled);
  boost::mpi::communicator world;

  // checking if the mpi supports the required threading level
  if(env.thread_level() < mt::funneled) {
    if(world.rank() == 0) {
      std::cerr << "This program requires MPI_THREAD_FUNNELED support from MPI. Exiting..." << std::endl;
    }
    return 1;
  }

  long int n;
  std::vector<SCC::Edge> edges;
  std::vector<SCC::Edge> decrement;
  std::vector<SCC::Edge> increament;
  std::vector<std::pair<long int, long int>> queries;

  if (world.rank() == 0)
  {
    for(int i = 0; i<7; i++) {
      std::string tmp;
      getline(std::cin, tmp);
    }
    getGraph(edges, n);

    //static algorithm
    auto s = std::chrono::high_resolution_clock::now();
    std::unordered_map<long int, long int> sccs;
    for (long int i = 0; i < n; i++)
      sccs[i] = i;
    SCC::findScc(edges, sccs);

    // emulating the dynamic updates
    float ratio = std::stof(argv[1]);
    int m = edges.size();
    int updates = m * ratio;
    srand(time(0));
    for(int i = 0; i<updates; i++) {
      int k = rand() % m;
      decrement.push_back(edges[k]);
      std::swap(edges[k], edges[m-1]);
      m--;
    }

    // emulating the re-run of the algorithm
    SCC::findScc(edges, sccs);
    auto e = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = e - s;
    std::cout << "STATIC: " << elapsed.count() << "s" << std::endl;
  }

  auto s = std::chrono::high_resolution_clock::now();
  // this structrue is maintained by the master but is shared with all the workers
  SCC::MaintainSCC scc(n, edges);

  if(world.rank() == 0) {
    auto e = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = e - s;
    std::cout << "INIT: " << elapsed.count() << "s" << std::endl;

    s = std::chrono::high_resolution_clock::now();
    // sending the updates to all the workers
    scc.deleteEdges(decrement);
    e = std::chrono::high_resolution_clock::now();
    elapsed = e - s;
    std::cout << "UPDATES: " << elapsed.count() << "s" << std::endl;
  }


  return 0;
}

// give set of unique number to generate from them