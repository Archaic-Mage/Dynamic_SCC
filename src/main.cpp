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
    edges.emplace_back(from, to);
  }
}

// Function to get the list of updates
void getUpdates(int &updates, std::vector<SCC::Edge> &decrement, std::vector<SCC::Edge> &increament)
{
  std::cin >> updates;
  for (int i = 0; i < updates; i++)
  {
    int type;
    std::cin >> type;
    long int from, to;
    std::cin >> from >> to;
    from++, to++;
    if (type == 0)
      decrement.emplace_back(from, to);
    else
      increament.emplace_back(from, to);
  }
}

void getQueries(int &q, std::vector<std::pair<long int, long int>> &queries)
{
  std::cin >> q;
  for (int i = 0; i < q; i++)
  {
    int v1, v2;
    std::cin >> v1 >> v2;
    v1++, v2++;
    queries.emplace_back(v1, v2);
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
    for(int i = 0; i<4; i++) {
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

    // int updates;
    // getUpdates(updates, decrement, increament);
    // int q;
    // getQueries(q, queries);
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

    SCC::findScc(edges, sccs);
    auto e = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = e - s;
    std::cout << "STATIC: " << elapsed.count() << "s" << std::endl;
  }

  auto start = std::chrono::high_resolution_clock::now();
  SCC::MaintainSCC scc(n, edges);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  if (world.rank() == 0)
    std::cout << "INIT: " << elapsed.count() << "s\n";

  start = std::chrono::high_resolution_clock::now();  
  scc.deleteEdges(decrement);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  if (world.rank() == 0)
    std::cout << "UPD: " << elapsed.count() << "s\n";

  if (world.rank() == 0) {
    // std::cout << "Number of SCCs: " << scc.getNumberOfSCCs() << std::endl;
    std::ofstream out("output.tmp");
    for (auto query : queries)
    {
      if(scc.query(query.first, query.second))
        out << "YES\n";
      else
        out << "NO\n";
    }
  }

  return 0;
}

// give set of unique number to generate from them