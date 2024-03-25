#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include "graph.hpp"

void findScc(const std::vector<Edge>& edges, std::vector<std::vector<long int>> &sccs) {
    std::unordered_map<long int, std::vector<long int>> adj;
    std::unordered_map<long int, std::vector<long int>> adj_rev;
    std::unordered_map<long int, bool> visited;
    std::unordered_set<long int> labels;
    std::vector<long int> order;
    std::vector<long int> component;

    std::function<void(void)> mark_unvisited = [&]() {
        for (const auto& label : labels) {
            visited[label] = false;
        }
    };
    for (Edge edge : edges) {
        labels.insert(edge.from);
        labels.insert(edge.to);
        adj[edge.from].push_back(edge.to);
        adj_rev[edge.to].push_back(edge.from);
    }
    std::function<void(long int)> dfs1 = [&](long int u) {
        visited[u] = 1;
        for (long int v : adj[u]) {
            if (!visited[v]) {
                dfs1(v);
            }
        }
        order.push_back(u);
    };
    std::function<void(long int)> dfs2 = [&](long int u) {
        visited[u] = 1;
        component.push_back(u);
        for (long int v : adj_rev[u]) {
            if (!visited[v]) {
                dfs2(v);
            }
        }
    };
    for (auto& i: labels) {
        if (!visited[i]) {
            dfs1(i);
        }
    }
    mark_unvisited();
    for (long int i = order.size()-1; i >= 0; i--) {
        long int u = order[i];
        if (!visited[u]) {
            dfs2(u);
            sccs.push_back(component);
            component.clear();
        }
    }
}