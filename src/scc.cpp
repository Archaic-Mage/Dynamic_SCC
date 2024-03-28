#include "graph.hpp"

void findScc(const std::vector<Edge>& edges, std::unordered_map<long int, long int> &sccs) {
    std::unordered_map<long int, std::vector<long int>> adj;
    std::unordered_map<long int, std::vector<long int>> adj_rev;
    std::unordered_map<long int, bool> visited;
    std::unordered_set<long int> labels;
    std::vector<long int> order;
    std::stack<long int> st;
    int mx = 0;

    for (Edge edge : edges) {
        labels.insert(edge.from);
        labels.insert(edge.to);
        adj[edge.from].push_back(edge.to);
        adj_rev[edge.to].push_back(edge.from);
    }
    std::function<void(void)> mark_unvisited = [&]() {
        for (const auto& label : labels) {
            if (mx < label) mx = label;
            visited[label] = false;
        }
    };
    std::function<void(long int)> dfs1 = [&](long int u) {
        st.push(u);
        visited[u] = 1;
        while(!st.empty()) {
            long int v = st.top();
            st.pop();
            order.push_back(v);
            for (long int w : adj[v]) {
                if (!visited[w]) {
                    st.push(w);
                    visited[w] = 1;
                }
            }
        }
    };
    std::function<void(long int)> dfs2 = [&](long int u) {
        st.push(u);
        visited[u] = 1;
        sccs[u] = mx;
        while(!st.empty()) {
            long int v = st.top();
            st.pop();
            for (long int w : adj_rev[v]) {
                if (!visited[w]) {
                    st.push(w);
                    visited[w] = 1;
                    sccs[w] = mx;
                }
            }
        }
    };
    mark_unvisited();
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
            mx++;
        }
    }
}