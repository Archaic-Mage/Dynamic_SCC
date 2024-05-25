#include "graph.hpp"
#include <chrono>
#include <fstream>
#include <climits>
#include <iostream>
#include <boost/mpi.hpp>

// Function to get the graph from the input (O(m))
void getGraph(std::ifstream &file ,std::vector<SCC::Edge> &edges, long int &n)
{
  long int m;
  file >> n >> m;
  for (long unsigned int i = 0; i < m; i++)
  {
    long int from, to;
    file >> from >> to;
    from++, to++;
    if(from == to)
      continue;
    edges.emplace_back(from, to);
  }
}

void getUpdates(std::ifstream &file, std::vector<SCC::Edge> &decrement, std::vector<SCC::Edge> &increament)
{
  long int m;
  file >> m;
  for (long unsigned int i = 0; i < m; i++)
  {
    int type;
    file >> type;
    long int from, to;
    file >> from >> to;
    from++, to++;
    if(type)
      increament.emplace_back(from, to);
    else
      decrement.emplace_back(from, to);
  }
}

void getQueries(std::ifstream &file, std::vector<std::pair<long int, long int>> &queries)
{
  long int m;
  file >> m;
  for (long unsigned int i = 0; i < m; i++)
  {
    long int from, to;
    file >> from >> to;
    from++, to++;
    queries.emplace_back(from, to);
  }
}

namespace mt = boost::mpi::threading;

int main(int argc, char *argv[])
{
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

  std::string filename = argv[2];
  std::ifstream file(filename);

  for(int i = 0; i<7; i++) {
    std::string tmp;
    getline(file, tmp);
  }
  getGraph(file,edges, n);

  // getUpdates(file, decrement, increament);

  // getQueries(file, queries);

  if (world.rank() == 0)
  {
    // the dynamic updates
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
  }

  boost::mpi::broadcast(world, decrement, 0);

  world.barrier();

  auto s = std::chrono::high_resolution_clock::now();
  // this structrue is maintained by the master but is shared with all the workers
  SCC::MaintainSCC scc(n, edges);

  world.barrier();
  if(world.rank() == 0) {
    auto e = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = e - s;
    std::cout << "INIT: " << elapsed.count() << "s" << std::endl;
  }

  world.barrier();
  int num_of_sccs = scc.getNumberOfSCCs();
  // if(world.rank() == 0) {
  //   std::cout << "Number of SCCs: " << num_of_sccs << std::endl;
  // }


  world.barrier();
  s = std::chrono::high_resolution_clock::now();
  // sending the updates to all the workers
  scc.deleteEdges(decrement);
  world.barrier();
  // scc.insertEdges(increament);
  // world.barrier();
  if(world.rank() == 0) {
    auto e = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = e - s;
    std::cout << "UPD: " << elapsed.count() << "s" << std::endl;
  }

  world.barrier();

  // std::vector<bool> ans;
  // for(auto &i : queries) {
  //   ans.push_back(scc.query(i.first, i.second));
  // }
  // world.barrier();
  // if(world.rank() == 0) {
  //   std::ofstream out("output.tmp");
  //   for(const auto &val: ans) {
  //     if(val) {
  //       out << "YES" << std::endl;
  //     } else {
  //       out << "NO" << std::endl;
  //     }
  //   }
  // }
  // world.barrier();
  return 0;
}

// give set of unique number to generate from them