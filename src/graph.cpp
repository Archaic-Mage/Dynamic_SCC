#include "graph.hpp"
#include "debug.hpp"

namespace SCC {

namespace COMPUTE_SCC {

int getVerticesDone(std::vector<int> & trim, int max_threads) {
    int vertices_done = 0;
    #pragma omp parallel for num_threads(max_threads) reduction(+:vertices_done)
    for(int i = 0; i<trim.size(); i++) {
        vertices_done += trim[i];
    }
    return vertices_done;
}

void tarjan_scc(std::vector< std::vector < int >> & adj_list, std::vector<int> & sccs, std::vector<int> & trim, int v, std::vector<int> & low, std::vector<int> & disc, std::stack<int> & st, std::vector<bool> & stack_member) {
    static int time = 0;
    disc[v] = low[v] = ++time;
    st.push(v);
    stack_member[v] = true;

    for(int c: adj_list[v]) {
        if(trim[c] == 1) continue;
        if(disc[c] == -1) {
            tarjan_scc(adj_list, sccs, trim, c, low, disc, st, stack_member);
            low[v] = std::min(low[v], low[c]);
        } else if(stack_member[c] == true) {
            low[v] = std::min(low[v], disc[c]);
        }
    }

    int w = 0;
    if(low[v] == disc[v]) {
        while(st.top() != v) {
            w = st.top();
            sccs[w] = v;
            trim[w] = 1;
            stack_member[w] = false;
            st.pop();
        }
        w = st.top();
        sccs[w] = v;
        trim[w] = 1;
        stack_member[w] = false;
        st.pop();
    }
}

void init_tarjan(std::vector< std::vector < int >> & adj_list, std::vector<int> & sccs, std::vector<int> & trim, int n) {
    std::vector<int> low(n, -1);
    std::vector<int> disc(n, -1);
    std::stack<int> st;
    std::vector<bool> stack_member(n, false);

    for(int i = 0; i<n; i++) {
        if(trim[i] == 0 && disc[i] == -1) {
            tarjan_scc(adj_list, sccs, trim, i, low, disc, st, stack_member);
        }
    }
}

void fwbw(std::vector< std::vector < int >> & adj_list, std::vector< std::vector < int >> & adj_list_rev, std::vector<int> & sccs, std::vector<int> & trim, int pivot, int n, int max_threads) {
    std::vector<int> visited(n, 0);
    std::vector<int> queue;
    std::vector<std::vector<int>> thread_queue(max_threads);
    queue.push_back(pivot);
    // Forward traversal
    while(!queue.empty()) {
        #pragma omp parallel for num_threads(max_threads)
        for(int v: queue) {
            int tid = omp_get_thread_num();
            for(int c: adj_list[v]) {
                if(trim[c] == 1) continue;
                if(!visited[c]) {
                    visited[c] = 1;
                    thread_queue[tid].push_back(c);
                }
            }
        }
        queue.clear();
        for(int i = 0; i<max_threads; i++) {
            for(int v: thread_queue[i]) {
                queue.push_back(v);
            }
            thread_queue[i].clear();
        }
    }
  

    // Backward traversal
    queue.push_back(pivot);
    sccs[pivot] = pivot;
    trim[pivot] = 1;
    while(!queue.empty()) {
        #pragma omp parallel for num_threads(max_threads)
        for(int v: queue) {
            int tid = omp_get_thread_num();
            for(int c: adj_list_rev[v]) {
                if(trim[c] == 1) continue;
                if(visited[c] == 1) {
                    visited[c] = 2;
                    sccs[c] = pivot;
                    trim[c] = 1;
                    thread_queue[tid].push_back(c);
                }
            }
        }
        queue.clear();
        for(int i = 0; i<max_threads; i++) {
            for(int v: thread_queue[i]) {
                queue.push_back(v);
            }
            thread_queue[i].clear();
        }
    }

}
    
void colorGraph(std::vector< std::vector < int >> & adj_list, std::vector< std::vector < int >> & adj_list_rev, std::vector<int> & sccs, std::vector<int> & trim, int n, int max_threads) {
    // Coloring
    std::vector<int> color(n, -1);
    std::vector<int> queue;
    std::vector<std::vector<int>> thread_queue(max_threads);
    std::vector<bool> in_next_queue(n,0);


    for(int i = 0; i<n; i++) {
        if(trim[i] == 0) {
            color[i] = i;
            queue.push_back(i);
            in_next_queue[i] = false;
        }
    }

    while(!queue.empty()) {
        #pragma omp parallel for num_threads(max_threads)
        for(int v: queue) {
            bool any_color_change = false;
            int tid = omp_get_thread_num();
            for(int c: adj_list[v]) {
                if(trim[c] == 1) continue;
                if(color[v] > color[c]) {
                    color[c] = color[v];
                    any_color_change = true;
                    if(!in_next_queue[c]) {
                        thread_queue[tid].push_back(c);
                        in_next_queue[c] = true;
                    }
                }
            }
            if(any_color_change) {
                if(!in_next_queue[v]) {
                    thread_queue[tid].push_back(v);
                    in_next_queue[v] = true;
                }
            }
        }
        queue.clear();
        for(int i = 0; i<max_threads; i++) {
            for(int v: thread_queue[i]) {
                queue.push_back(v);
                in_next_queue[v] = false;
            }
            thread_queue[i].clear();
        }
    }


    for(int i = 0; i<n; i++) {
        if(trim[i] == 0) {
            if(!in_next_queue[color[i]]) {
                queue.push_back(color[i]);
                sccs[color[i]] = color[i];
                trim[color[i]] = 1;
                in_next_queue[color[i]] = true;
            }
        }
    }

    // #pragma omp parallel for num_threads(max_threads)
    while(!queue.empty()) {
        #pragma omp parallel for num_threads(max_threads)
        for(int v: queue) {
            int tid = omp_get_thread_num();
            for(int c: adj_list_rev[v]) {
                if(!trim[c] && color[c] == color[v]) {
                    sccs[c] = sccs[v];
                    trim[c] = 1;
                    if(!in_next_queue[c]) {
                        thread_queue[tid].push_back(c);
                        in_next_queue[c] = true;
                    }
                }
            }
        }
        queue.clear();
        for(int i = 0; i<max_threads; i++) {
            for(int v: thread_queue[i]) {
                queue.push_back(v);
                in_next_queue[v] = false;
            }
            thread_queue[i].clear();
        }
    }
}

// requires continuous vertex labels
void find_scc(std::vector< std::vector < int >> & adj_list, std::vector< std::vector < int >> &adj_list_rev, std::vector<int> & sccs, int n, int m, int max_threads) {

    // degree of each vertex
    std::vector<int> in_degree(n, 0);
    std::vector<int> out_degree(n, 0);
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i < n; i++) {
        out_degree[i] = adj_list[i].size();
        in_degree[i] = adj_list_rev[i].size();
    }

        // trim the graph - remove vertices with degree(in/out) = 0 one-step
    std::vector<int> trim(n, 0);

    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<n; i++) {
        if(in_degree[i] == 0 || out_degree[i] == 0) {
            trim[i] = 1;
            sccs[i] = i;
        }
    }

    // int trimmed_vertices = getVerticesDone(trim, max_threads);
    // std::cout << "Number of vertices trimmed: " << trimmed_vertices << std::endl;

    //forward backward traversal
    int pivot = 0;
    long int mx_deg_mul = 0;

    for(int i = 0; i<n; i++) {
        if((long int) in_degree[i]*out_degree[i] > mx_deg_mul) {
            mx_deg_mul = in_degree[i]*out_degree[i];
            pivot = i;
        }
    }

    fwbw(adj_list, adj_list_rev, sccs, trim, pivot, n, max_threads);

    // int fw_bw_vertices = getVerticesDone(trim, max_threads);
    // std::cout << "Number of vertices done in FW-BW: " << fw_bw_vertices << std::endl;

    //coloring
    int vertex_done = getVerticesDone(trim, max_threads);
    while(vertex_done < n && (n-vertex_done) > 100000) {
        // vertex_done = 0;
        colorGraph(adj_list, adj_list_rev, sccs, trim, n, max_threads);
        
        vertex_done = getVerticesDone(trim, max_threads);
    }

    // int colored_vertices = getVerticesDone(trim, max_threads);
    // std::cout << "Number of vertices done in coloring: " << colored_vertices << std::endl;

    // tarjans
    init_tarjan(adj_list, sccs, trim, n);

    // int scc_vertices = getVerticesDone(trim, max_threads);
    // std::cout << "Number of vertices done in Tarjan's: " << scc_vertices << std::endl;

}

}


int MaintainSCC::getLabel()
{
    int ret = SCC_LABEL;
    SCC_LABEL += world_size;
    return ret;
}

void dfs(long int u, std::unordered_map<long int, std::vector<long int>> &adj, std::vector<long int> &order, std::unordered_map<long int, int> &visited)
{
    visited[u] = true;
    for (long int v : adj[u])
    {
        if (!visited[v])
        {
            dfs(v, adj, order, visited);
        }
    }
    order.push_back(u);
}


void findScc(const std::vector<Edge>& edges, std::unordered_map<long int, long int> &sccs) {
    int n = sccs.size();
    int m = edges.size();
    int max_threads = omp_get_max_threads();
    std::vector<Edge> edges_copy(m);
    std::vector< std::vector < int >> adj_list_new(n);
    std::vector< std::vector < int >> adj_list_rev_new(n);
    std::vector<int> sccs_new(n, -1);
    std::unordered_map<int, int> vertex_hash;
    std::vector<int> vertex_unhash(n,0);
    
    vertex_hash.reserve(n);
    //filling the vertex_map
    for(auto it = sccs.begin(); it != sccs.end(); it++) {
        int at = vertex_hash.size();
        vertex_hash[it->first] = at;
        vertex_unhash[at] = it->first;
    }

    // change the vertex labels in the edges
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<m; i++) {
        edges_copy[i].from = vertex_hash[edges[i].from];
        edges_copy[i].to = vertex_hash[edges[i].to];
    }

    //locks required for adj_list_new
    omp_lock_t** lock = new omp_lock_t*[n];
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<n; i++) {
        lock[i] = new omp_lock_t;
        omp_init_lock(lock[i]);
    }

    //filling the adj_list_new and adj_list_rev_new
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<m; i++) {
        int u = edges_copy[i].from;
        int v = edges_copy[i].to;
        omp_set_lock(lock[u]);
        adj_list_new[u].push_back(v);
        omp_unset_lock(lock[u]);
    }
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<m; i++) {
        int u = edges_copy[i].from;
        int v = edges_copy[i].to;
        omp_set_lock(lock[v]);
        adj_list_rev_new[v].push_back(u);
        omp_unset_lock(lock[v]);
    }

    COMPUTE_SCC::find_scc(adj_list_new, adj_list_rev_new, sccs_new, n, m, max_threads);
    
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<n; i++) {
        sccs[vertex_unhash[i]] = sccs_new[i];
    }
}

void TreeNode::removeEdge(const Edge &edge)
{
    //erase using key
    int max_threads = omp_get_max_threads();
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<max_threads; i++) {
        corresponds_to[i].erase(edge);
    }
}

void TreeNode::condenseFill(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs)
{
    int max_threads = omp_get_max_threads();
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<edges.size(); i++) {
        int tid = omp_get_thread_num();
        if(sccs[edges[i].from] != sccs[edges[i].to]) {
            Edge e;
            e.from = sccs[edges[i].from];
            e.to = sccs[edges[i].to];
            if(edges[i].to < 0) {
                e.to = -e.to;
            }
            corresponds_to[tid][edges[i]] = e;
        }
    }
}

void TreeNode::checkUnreachable(std::unordered_set<long int> &unreachable)
{
    int max_threads = omp_get_max_threads();
    // get edge list
    std::vector<Edge> edges;
    std::vector<int> prev_edges_size(max_threads, 0);
    for(int i = 1; i<max_threads; i++) {
        prev_edges_size[i] = prev_edges_size[i-1] + corresponds_to[i-1].size();
    }
    int edges_size = prev_edges_size[max_threads-1] + corresponds_to[max_threads-1].size();
    edges.resize(edges_size);

    long int split_on = -1;
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<max_threads; i++) {
        int j = prev_edges_size[i];
        for(auto edge : corresponds_to[i]) {
            if(edge.second.to < 0) {
                split_on = -edge.second.to;
                edge.second.to = split_on;
            }
            edges[j++] = edge.second;
        }
    }

    // get the nodes in the current scc
    std::unordered_map<long int, long int> sccs;
    for (auto &node : contains)
    {
        sccs[node] = node;
    }
    findScc(edges, sccs);
    int split_on_label = sccs[split_on];
    for (const auto &node : contains)
    {
        if (sccs[node] != split_on_label)
        {
            unreachable.insert(node);
        }
    }
}

void TreeNode::updateLabels(std::unordered_map<long int, long int> &new_labels)
{
    int max_threads = omp_get_max_threads();
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<max_threads; i++) 
        for (auto &correspond_to : corresponds_to[i])
        {
            Edge &label_edge = correspond_to.second;
            const Edge &actual_edge = correspond_to.first;

            if ( new_labels.find ( actual_edge.from ) != new_labels.end() )
                label_edge.from = new_labels[actual_edge.from];
            if ( new_labels.find ( actual_edge.to ) != new_labels.end() )
                label_edge.to = new_labels[actual_edge.to];
        }
}

void TreeNode::exposeToParent(TreeNode &parent, std::unordered_set<long int> unreachable)
{
    int max_threads = omp_get_max_threads();
    #pragma omp parallel for num_threads(max_threads)
    for(int i= 0; i<max_threads; i++) {
        std::unordered_map<Edge, Edge, EdgeHash> temp;
        for(auto &correspond_to: corresponds_to[i])
        {
            if (unreachable.find(correspond_to.second.from) != unreachable.end() || unreachable.find(correspond_to.second.to) != unreachable.end())
            {
                if (unreachable.find(correspond_to.second.from) == unreachable.end())
                {
                    correspond_to.second.from = label;
                }
                if (unreachable.find(correspond_to.second.to) == unreachable.end())
                {
                    correspond_to.second.to = label;
                }
                parent.corresponds_to[i][correspond_to.first] = correspond_to.second;
            }
            else
            {
                temp[correspond_to.first] = correspond_to.second;
            }
        }
        corresponds_to[i].swap(temp);
    }
}

// traverse the Node and all the containing nodes
void MaintainSCC::traverseNode(int node, std::unordered_map<long int, long int> &sccs)
{
    std::queue<long int> q;
    q.push(node);
    while (!q.empty())
    {
        long int curr = q.front();
        q.pop();
        sccs[curr] = node;
        for (auto &child : scc_tree_nodes[curr].contains)
        {
            if (child != curr)
            {
                q.push(child);
            }
        }
    }
}

// split a node's edges into two sets (node_in and node_out) creating new graph (O(m))
void MaintainSCC::splitGraphOnNode(std::vector<Edge> &edges, long int node)
{
    int max_threads = omp_get_max_threads();
    #pragma omp parallel for num_threads(max_threads)
    for (int i = 0; i<edges.size(); i++) 
    {
        if (edges[i].to == node)
        {
            edges[i].to = -node;
        }
    }
}

// Divide edges by groups based on SCCs (O(mlog(n)))
void MaintainSCC::divideEdgesByScc(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs, std::unordered_map<long int, std::vector<Edge>> &sccEdges, std::vector<Edge> &inter_scc_edges, bool split = false)
{
    int max_threads = omp_get_max_threads();
    std::unordered_map<int, std::vector<std::vector<Edge>>> thread_scc_edges;
    std::vector<std::vector<Edge>> thread_inter_scc_edges(max_threads);

    for(auto it = sccs.begin(); it != sccs.end(); it++) {
        thread_scc_edges[it->second].resize(max_threads);
    }

    #pragma omp parallel for num_threads(max_threads)
    for (int i = 0; i<edges.size(); i++)
    {
        int tid = omp_get_thread_num();
        if (sccs[edges[i].from] == sccs[edges[i].to])
        {
            thread_scc_edges[sccs[edges[i].from]][tid].push_back(edges[i]);
        } else if(split) {
            thread_inter_scc_edges[tid].push_back(edges[i]);
        }
    }

    std::vector<int> prev_edges_size(max_threads, 0);
    for(int i = 1; i<max_threads; i++) {
        prev_edges_size[i] = prev_edges_size[i-1] + thread_inter_scc_edges[i-1].size();
    }
    int edges_size = prev_edges_size[max_threads-1] + thread_inter_scc_edges[max_threads-1].size();
    inter_scc_edges.resize(edges_size);
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<max_threads; i++) {
        int j = prev_edges_size[i];
        for(auto edge: thread_inter_scc_edges[i]) {
            inter_scc_edges[j++] = edge;
        }
    }

    for(auto &scc_edges: thread_scc_edges) {
        for(int i = 1; i<max_threads; i++) {
            prev_edges_size[i] = prev_edges_size[i-1] + scc_edges.second[i-1].size();
        }
        edges_size = prev_edges_size[max_threads-1] + scc_edges.second[max_threads-1].size();
        sccEdges[scc_edges.first].resize(edges_size);
        #pragma omp parallel for num_threads(max_threads)
        for(int i = 0; i<max_threads; i++) {
            int j = prev_edges_size[i];
            for(auto edge: scc_edges.second[i]) {
                sccEdges[scc_edges.first][j++] = edge;
            }
        }
    }
}

int MaintainSCC::getSplitNode(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs)
{
    // int n = sccs.size();
    // int m = edges.size();
    // int max_threads = omp_get_max_threads();
    // std::unordered_map<int , int> label_hash;
    // std::vector<int> label_unhash(sccs.size(), 0);
    // for(const auto &[node, _]: sccs) {
    //     label_hash[node] = label_hash.size();
    //     label_unhash[label_hash.size()-1] = node;
    // }
    // std::vector<int> in_degree(sccs.size(), 0);
    // std::vector<int> out_degree(sccs.size(), 0);
    // #pragma omp parallel for num_threads(max_threads)
    // for(int i = 0; i<m; i++) {
    //     out_degree[label_hash[sccs[edges[i].from]]]++;
    //     in_degree[label_hash[sccs[edges[i].to]]]++;
    // }

    // int deg_mul = 0, pivot = 0;
    // #pragma omp parallel for num_threads(max_threads) reduction(max:deg_mul)
    // for (int i = 0; i < sccs.size(); i++)
    // {
    //     deg_mul = std::max(deg_mul, in_degree[i] * out_degree[i]);
    // }
    // #pragma omp parallel for num_threads(max_threads)
    // for(int i = 0; i<n; i++) {
    //     if(in_degree[i]*out_degree[i] == deg_mul) {
    //         pivot = label_unhash[i];
    //     }
    // }
    int pivot = edges[0].from;
    return pivot;
}

/*** Functions to make SCC Tree ***/
void MaintainSCC::makeSccTreeInternals(std::vector<Edge> &edge, std::vector<long int> &nodes, TreeNode &currentNode)
{
    // actual node -> temp label
    std::unordered_map<long int, long int> sccs;
    for(auto &node: nodes) {
        sccs[node] = node;
    }

    int n = sccs.size();
    int m = edge.size();
    int max_threads = omp_get_max_threads();

    scc_tree_labels[currentNode.label] = currentNode.label;


    //get split on vertex (we choose vertex with highest degree (total))
    long int split_on = getSplitNode(edge, sccs);

    // O(m)
    splitGraphOnNode(edge, split_on);

    sccs[-split_on] = -split_on;

    // Optimizing timing by reducing the depts (capping at MAX_DEPT)
    int next_dept = currentNode.dept + 1;
    if( next_dept > MAX_DEPT ) {
        for(auto &node: sccs) {
            if(node.first == -split_on) continue;
            scc_tree_labels[node.first] = currentNode.label;
        }
        currentNode.condenseFill(edge, sccs);
        return;
    }

    findScc(edge, sccs);


    // O(m)
    // temp label -> edges in scc
    std::unordered_map<long int, std::vector<Edge>> sccEdges;
    std::vector<Edge> temp;
    divideEdgesByScc(edge, sccs, sccEdges, temp);

    // O(n)
    std::unordered_map<long int, std::vector<long int>> sccNodes;
    for (auto &scc : sccs)
    {
        sccNodes[scc.second].push_back(scc.first);
    }

    // O(n)
    // changing sccs mappings to labels
    for (auto &node : sccNodes)
    {
        bool single_node = (node.second.size() == 1);
        int label = single_node ? node.second[0] : getLabel();
        for(auto &nd: node.second)
            sccs[nd] = label;
        if(label == -split_on) continue;
        TreeNode &child = scc_tree_nodes[label];
        child.label = label;
        child.parent = currentNode.label;
        child.dept = currentNode.dept + 1;
        currentNode.contains.insert(child.label);
        if (!single_node)
            makeSccTreeInternals(sccEdges[node.first], node.second, child);
        else
            scc_tree_labels[node.second[0]] = label;
    }

    // O(m)
    currentNode.condenseFill(edge, sccs);
}

// Construct the SCC tree
void MaintainSCC::makeSccTree(std::vector<Edge> &edges, std::vector<long int> &nodes)
{
    // case when only one node is present in the scc
    if (nodes.size() == 1)
    {
        int label = nodes[0];
        roots.push_back(label);
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
    makeSccTreeInternals(edges, nodes, root);
}

void MaintainSCC::changeSccLabels(std::unordered_map<long int, long int> &sccs, std::unordered_map<long int, std::vector<long int>> &sccNodes)
{
    std::unordered_map<long int, long int> new_labels;
    for (auto &scc : sccs)
    {
        if (new_labels.find(scc.second) == new_labels.end())
            if (sccNodes[scc.second].size() == 1)
                new_labels[scc.second] = sccNodes[scc.second][0];
            else
                new_labels[scc.second] = getLabel();
    }
    for (auto &scc : sccs)
    {
        scc.second = new_labels[scc.second];
    }
}

// Construct Master Node
void MaintainSCC::constructMasterNode(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs)
{
    int max_threads = omp_get_max_threads();
    // Create a TreeNode
    int label = 0;
    roots.push_back(label);
    TreeNode &root = scc_tree_nodes[label];
    root.label = label;
    root.parent = label;
    root.dept = 0;
    #pragma omp parallel for num_threads(max_threads)
    for (int i = 0; i<edges.size(); i++)
    {
        int tid = omp_get_thread_num();
        if (sccs[edges[i].from] != sccs[edges[i].to])
        {
            Edge e;
            e.from = sccs[edges[i].from];
            e.to = sccs[edges[i].to];
            root.corresponds_to[tid][edges[i]] = e;
        }
    }
    for (auto &scc : sccs)
    {
        root.contains.insert(scc.second);
    }
}

// update master node 0-new scc creation
void MaintainSCC::updateMasterNode(int type, std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs)
{
    int max_threads = omp_get_max_threads();
    if (type == 0)
    {
        #pragma omp parallel for num_threads(max_threads)
        for(int i = 0; i<max_threads; i++)
            for (auto &correspond_to : scc_tree_nodes[0].corresponds_to[i])
            {
                const Edge &actual_edge = correspond_to.first;
                Edge &label_edge = correspond_to.second;
                long int from = actual_edge.from;
                if (sccs.find(from) != sccs.end())
                {
                    label_edge.from = sccs[from];
                }
                long int to = actual_edge.to;
                if (sccs.find(to) != sccs.end())
                {
                    label_edge.to = sccs[to];
                }
            }
        #pragma omp parallel for num_threads(max_threads)
        for (int i = 0; i<edges.size(); i++)
        {
            int tid = omp_get_thread_num();
            if (sccs[edges[i].from] != sccs[edges[i].to])
            {
                Edge e;
                e.from = sccs[edges[i].from];
                e.to = sccs[edges[i].to];
                scc_tree_nodes[0].corresponds_to[tid][edges[i]] = e;
            }
        }
        for (auto &scc : sccs)
        {
            scc_tree_nodes[0].contains.insert(scc.second);
        }
    }
}

int MaintainSCC::checkAndRemoveUnreachable(TreeNode &curr_node)
{
    int max_threads = omp_get_max_threads();
    // do unreachability and reachability checks
    std::unordered_set<long int> unreachable;
    curr_node.checkUnreachable(unreachable);

    // if no unreachable nodes then return
    if (unreachable.empty())
        return 0;

    for (auto &node : unreachable)
    {
        dPrint("Unreachable: " + std::to_string(node));
    }

    // remove all unreachable nodes from the current node
    for (auto &unreachable_node : unreachable)
    {
        curr_node.contains.erase(unreachable_node);
    }
    // add the unreachable nodes to the parent
    TreeNode &parent = scc_tree_nodes.at(curr_node.parent);
    // change parent label of unreachable nodes
    for (auto &unreachable_node : unreachable)
    {
        scc_tree_nodes[unreachable_node].parent = parent.label;
        parent.contains.insert(unreachable_node);
    }
    std::unordered_map<long int, long int> new_labels;
    if (curr_node.contains.size() == 1)
    {
        long int new_label = *(curr_node.contains.begin());
        parent.contains.erase(curr_node.label);
        #pragma omp parallel for num_threads(max_threads)
        for(int i = 0; i<max_threads; i++)
            for(auto &edge_pair: parent.corresponds_to[i]) {
                Edge &label_edge = edge_pair.second;
                if(label_edge.from == curr_node.label) {
                    label_edge.from = new_label;
                }
                if(label_edge.to == curr_node.label) {
                    label_edge.to = new_label;
                }
            }
        curr_node.label = new_label;
        parent.contains.insert(new_label);
    }
    for(const auto &nd: unreachable) {
        traverseNode(nd, new_labels);
    }
    for (auto &new_label : new_labels)
    {
        dPrint("New Label: " + std::to_string(new_label.first) + " -> " + std::to_string(new_label.second));
    }
    dTreeNode(parent);
    // update the parent with the new labels
    parent.updateLabels(new_labels);

    // exposing internal structure to parent node
    curr_node.exposeToParent(parent, unreachable);
    return 1;
}

int MaintainSCC::checkAndDetachUnreachable(TreeNode &root_node)
{
    boost::mpi::communicator world;
    int max_threads = omp_get_max_threads();
    std::unordered_set<long int> unreachable;
    root_node.checkUnreachable(unreachable);
    if (unreachable.empty())
    {
        return 0;
    }
    for (auto &node : unreachable)
    {
        dPrint("Unreachable: " + std::to_string(node));
    }
    for (auto &unreachable_node : unreachable)
    {
        root_node.contains.erase(unreachable_node);
    }
    TreeNode temp;
    root_node.exposeToParent(temp, unreachable);
    std::vector<Edge> edges;
    std::vector<int> prev_edges_size(max_threads, 0);
    int edges_size = 0;
    for(int i = 1; i<max_threads; i++) {
        prev_edges_size[i] = prev_edges_size[i-1] + temp.corresponds_to[i-1].size();
    }
    edges_size = prev_edges_size[max_threads-1] + temp.corresponds_to[max_threads-1].size();
    edges.resize(edges_size);
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<max_threads; i++) {
        int j = prev_edges_size[i];
        for(auto &correspond_to: temp.corresponds_to[i]) {
            edges[j++] = correspond_to.first;
        }
    }
    std::unordered_map<long int, long int> new_sccs;
    for (auto &unreachable_node : unreachable)
    {
        roots.push_back(unreachable_node);
        traverseNode(unreachable_node, new_sccs);
    }
    // transfer to master node for creating new scc
    world.send(0, 0, STATUS::DONE_NEW);
    world.send(0, 1, edges);
    world.send(0, 2, new_sccs);

    return 0;
}

/*** Functions to update the graph ***/

int MaintainSCC::getLCANode(int node1, int node2)
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

// DELETE EDGE
void MaintainSCC::deleteEdge(Edge edge)
{   
    long int from = scc_tree_labels[edge.from];
    long int to = scc_tree_labels[edge.to];
    TreeNode &lca = scc_tree_nodes.at(getLCANode(from, to));
    // delete the edge from the lca
    lca.removeEdge(edge);
    // add the lca to the delete cache
    delete_cache.insert({lca.dept, lca.label});
}

bool MaintainSCC::clearDeleteCache(TreeNode &node)
{
    if (node.dept != 0)
        return checkAndRemoveUnreachable(node);
    return checkAndDetachUnreachable(node);
}

void MaintainSCC::deleteEdgeFromMaster(Edge edge)
{
    TreeNode &root = scc_tree_nodes[0];
    root.removeEdge(edge);
}

void MaintainSCC::insertEdgesInMaster(std::vector< Edge > &edges) 
{
    int max_threads = omp_get_max_threads();
    TreeNode &root = scc_tree_nodes[0];

    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<edges.size(); i++) {
        int tid = omp_get_thread_num();
        Edge label_edge(sccs[edges[i].from], sccs[edges[i].to]);
        root.corresponds_to[tid][edges[i]] = label_edge;
    }
    edges.clear();
    std::vector<int> prev_edges_size(max_threads, 0);
    int edges_size = 0;
    for(int i = 1; i<max_threads; i++) {
        prev_edges_size[i] = prev_edges_size[i-1] + root.corresponds_to[i-1].size();
    }
    edges_size = prev_edges_size[max_threads-1] + root.corresponds_to[max_threads-1].size();
    edges.resize(edges_size);
    // filling the edges for computing sccs
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<max_threads; i++) {
        int j = prev_edges_size[i];
        for(auto &correspond_to: root.corresponds_to[i]) {
            edges[j++] = correspond_to.second;
        }
    }
    std::unordered_map<long int, long int> sccs;
    for(const auto &node: root.contains)
        sccs[node] = node;
    findScc(edges, sccs);

    dScc(sccs);
    dTreeNode(root);
    std::unordered_map<long int, std::vector<long int>> sccNodes;
    for(auto &scc: sccs)
        sccNodes[scc.second].push_back(scc.first);
    std::unordered_map<long int, long int> actual_node;
    //dividing the edges based on the new sccs formed (restrucring the root master node)
    // std::unordered_map<long int, std::vector<std::pair<const Edge, Edge>>> new_scc_root;
    // std::unordered_map<Edge, Edge, EdgeHash> temp;
    // for(auto &correspond_to: root.corresponds_to) {
    //     actual_node[correspond_to.second.from] = correspond_to.first.from;
    //     actual_node[correspond_to.second.to] = correspond_to.first.to;
    //     if(sccs[correspond_to.second.from] == sccs[correspond_to.second.to]) {
    //         new_scc_root[sccs[correspond_to.second.from]].emplace_back(correspond_to);
    //     } else {
    //         temp.emplace_hint(temp.end(), correspond_to.first, correspond_to.second);
    //     }
    // }
    // root.corresponds_to.swap(temp);
}

// INSERT EDGE
void MaintainSCC::insertEdge(Edge edge)
{
    TreeNode &lca = scc_tree_nodes.at(getLCANode(edge.from, edge.to));
    int label_from = edge.from;
    while (lca.label != scc_tree_nodes[label_from].parent)
    {
        label_from = scc_tree_nodes[label_from].parent;
    }
    int label_to = edge.to;
    while (lca.label != scc_tree_nodes[label_to].parent)
    {
        label_to = scc_tree_nodes[label_to].parent;
    }
    Edge e(label_from, label_to);

    int min_thread = 0;
    int min_thread_size = INT_MAX;
    for (int i = 0; i < omp_get_max_threads(); i++)
    {
        if (lca.corresponds_to[i].size() < min_thread_size)
        {
            min_thread = i;
            min_thread_size = lca.corresponds_to[i].size();
        }
    }

    // add the edge to the lca
    lca.corresponds_to[min_thread][edge] = e;
    // add the lca to the insert cache
    insert_cache.insert({lca.dept, lca.label});
}

void MaintainSCC::clearInsertCache(TreeNode &node)
{
    int max_threads = omp_get_max_threads();
    std::vector<Edge> edges;
    std::vector<int> prev_edges_size(omp_get_max_threads(), 0);
    int edges_size = 0;
    for(int i = 1; i<omp_get_max_threads(); i++) {
        prev_edges_size[i] = prev_edges_size[i-1] + node.corresponds_to[i-1].size();
    }
    edges_size = prev_edges_size[omp_get_max_threads()-1] + node.corresponds_to[omp_get_max_threads()-1].size();
    edges.resize(edges_size);
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<max_threads; i++) {
        int j = prev_edges_size[i];
        for(auto &correspond_to: node.corresponds_to[i]) {
            edges[j++] = correspond_to.second;
        }
    }

    std::unordered_map<long int, long int> sccs;
    for (auto &node : node.contains)
        sccs[node] = node;
    findScc(edges, sccs);

    std::unordered_map<long int, std::vector<long int>> sccNodes;
    for (auto &scc : sccs)
        sccNodes[scc.second].push_back(scc.first);

    bool to_end = true;
    for (auto &scc_nodes : sccNodes)
        if (scc_nodes.second.size() > 1)
            to_end = false;
    if (to_end)
        return;

    changeSccLabels(sccs, sccNodes);

    sccNodes.clear();
    for (auto &scc : sccs)
        sccNodes[scc.second].push_back(scc.first);

    std::vector<long int> new_tree_nodes;
    std::unordered_map<long int, long int> split_on;
    for (auto &scc : sccs)
    {
        if (scc.second < 0)
            continue;
        if (scc_tree_nodes.find(scc.second) == scc_tree_nodes.end())
        {
            TreeNode &child = scc_tree_nodes[scc.second];
            child.label = scc.second;
            child.parent = node.label;
            child.dept = node.dept + 1;
            child.contains.insert(sccNodes[scc.second].begin(), sccNodes[scc.second].end());
            split_on[scc.second] = sccNodes[scc.second][0];
            new_tree_nodes.push_back(child.label);
        }
    }

    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<max_threads; i++) {
        std::unordered_map<Edge, Edge, EdgeHash> temp;
        for(auto &correspond_to: node.corresponds_to[i])
        {
            int scc_label_from = sccs[correspond_to.second.from], scc_label_to = sccs[correspond_to.second.to];
            if (scc_label_from != scc_label_to)
                temp[correspond_to.first] = correspond_to.second;
            else {
                if(correspond_to.second.to == split_on[scc_label_to])
                    correspond_to.second.to = -correspond_to.second.to;
                scc_tree_nodes[sccs[correspond_to.second.from]].corresponds_to[i][correspond_to.first] = correspond_to.second;
            }
        }
        node.corresponds_to[i].swap(temp);
    }

    // changing connections in node
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i< max_threads; i++)
        for (auto &correspond_to : node.corresponds_to[i])
        {
            Edge &label_edge = correspond_to.second;
            const Edge &actual_edge = correspond_to.first;
            label_edge.from = sccs[label_edge.from];
            label_edge.to = sccs[label_edge.to];
        }

    // update the contains set
    node.contains.clear();
    for (auto &scc : sccs)
    {
        if (scc.second < 0)
            continue;
        node.contains.insert(scc.second);
    }
    for (auto &new_tree_node : new_tree_nodes)
    {
        TreeNode &child = scc_tree_nodes[new_tree_node];
        insert_cache.insert({child.dept, child.label});
    }
}

// query if two nodes are in the same SCC
bool MaintainSCC::query(long int v1, long int v2)
{
    return sccs[v1] == sccs[v2];
}

void MaintainSCC::processMessage()
{
    boost::mpi::communicator world;
    MessageType type;
    while (true)
    {
        world.recv(0, 0, type);
        if (type == MessageType::EXIT)
            break;
        else if(type == MessageType::LABELS) {
            world.recv(0, 1, SCC_LABEL);
        }
        else if (type == MessageType::SCC_TREE)
        {
            std::vector<long int> nodes;
            std::vector<Edge> edgeList;
            world.recv(0, 1, nodes);
            world.recv(0, 2, edgeList);
            makeSccTree(edgeList, nodes);
            world.send(0, 0, STATUS::DONE_NO_NEW);
        }
        else if (type == MessageType::EDGE_DELETE)
        {
            Edge edge;
            world.recv(0, 1, edge);
            // std::cout << "Deleting Edge: " << edge.from << " " << edge.to << std::endl;
            deleteEdge(edge);
            // std::cout << "Deleted Edge: " << edge.from << " " << edge.to << std::endl;
        }
        else if (type == MessageType::CLR_DEC_CACHE)
        {
            // std::cout << "Clearing Delete Cache : " << world.rank() << std::endl;
            while (!delete_cache.empty())
            {
                Cache c = *delete_cache.begin();
                delete_cache.erase(delete_cache.begin());
                TreeNode &node = scc_tree_nodes.at(c.label);
                dTreeNode(node);
                if (clearDeleteCache(node))
                {
                    delete_cache.insert({node.dept - 1, node.parent});
                }
            }
            world.send(0, 0, STATUS::DONE_NO_NEW);
            // std::cout << "Cleared Delete Cache : " << world.rank() << std::endl;
        }
        else if (type == MessageType::EDGE_INSERT)
        {
            Edge edge;
            world.recv(0, 1, edge);
            // std::cout << "Inserted Edge: " << edge.from << " " << edge.to << std::endl;
            insertEdge(edge);
        }
        else if (type == MessageType::CLR_INC_CACHE)
        {
            while (!insert_cache.empty())
            {
                dInfo(world, "Inserting Cache");
                Cache c = *insert_cache.begin();
                insert_cache.erase(insert_cache.begin());
                TreeNode &node = scc_tree_nodes.at(c.label);
                clearInsertCache(node);
                dTreeNode(node);
            }
            world.send(0, 0, STATUS::DONE_NO_NEW);
        }
    }
}

MaintainSCC::MaintainSCC(long int n, std::vector<Edge> &edges)
{
    boost::mpi::communicator world;

    int world_rank = world.rank();
    int max_threads = omp_get_max_threads();
    int m = edges.size();

    world_size = world.size();
    SCC_LABEL = n + 1 + world_rank;
    
    for(int i = 1; i<=n; i++) {
        sccs[i] = i;
    }

    findScc(edges, sccs);

    std::unordered_map<long int, std::vector<long int>> sccNodes;
    for(const auto &scc: sccs) 
        sccNodes[scc.second].push_back(scc.first);

    changeSccLabels(sccs, sccNodes);

    std::unordered_map<long int, std::vector<Edge>> sccEdges;
    std::vector<Edge> inter_scc_edges;
    divideEdgesByScc(edges, sccs, sccEdges, inter_scc_edges, true);

    //update the sccNodes
    sccNodes.clear();
    for(const auto &scc: sccs) 
        sccNodes[scc.second].push_back(scc.first);


    std::cout << "Number of SCCs: " << sccNodes.size() << std::endl;


    int temp = 0;
    //check and distribute the sccs
    for(auto it = sccEdges.begin(); it != sccEdges.end(); it++) {
        int scc_size = it->second.size();
        if(scc_size > 0.1 * m) {
            int to_rank = temp % (world_size - 1) + 1;
            which_rank[it->first] = to_rank;
            if(to_rank == world_rank){
                makeSccTree(it->second, sccNodes[it->first]);
            }
            temp++;
        } else {
            int to_rank = 0;
            which_rank[it->first] = to_rank;
            if(to_rank == world_rank)
                makeSccTree(it->second, sccNodes[it->first]);
        }
    }

    world.barrier();

    // if (world_rank == 0)
    // {
    //     SCC_LABEL = n + 1;
    //     sccs.reserve(n+1);
    //     for (long int i = 1; i <= n; i++)
    //         sccs[i] = i;
    //     findScc(edges, sccs);

    //     std::cout << "Done FIND SCC\n" << std::endl;

    //     std::unordered_map<long int, std::vector<Edge>> sccEdges;
    //     std::vector<Edge> inter_scc_edges;
    //     divideEdgesByScc(edges, sccs, sccEdges, inter_scc_edges, true);

    //     std::cout << "Done DIVIDE EDGES\n" << std::endl;

    //     std::unordered_map<long int, std::vector<long int>> sccNodes;
    //     for (auto &scc : sccs)
    //     {
    //         sccNodes[scc.second].push_back(scc.first);
    //         which_rank[scc.first] = scc.second % (world.size() - 1) + 1;
    //     }

    //     std::cout << "Done DIVIDE NODES\n" << std::endl;

    //     std::cout << "Number of SCCs: " << sccNodes.size() << std::endl;
    //     std::cout << "Set of Edge: " << sccEdges.size() << std::endl;

    //     for(int i = 1; i<world_size; i++) {
    //         world.send(i, 0, MessageType::LABELS);
    //         world.send(i, 1, SCC_LABEL+i);
    //     }

    //     for(auto &scc_edge: sccEdges) {
    //         int to_rank = which_rank[scc_edge.second[0].from];
    //         world.send(to_rank, 0, MessageType::SCC_TREE);
    //         world.send(to_rank, 1, sccNodes[scc_edge.first]);
    //         world.send(to_rank, 2, scc_edge.second);
    //     }

    //     // for(auto &scc_node: sccNodes) {
    //     //     int to_rank = which_rank[scc_node.second[0]];
    //     //     world.send(to_rank, 0, MessageType::SCC_TREE);
    //     //     world.send(to_rank, 1, scc_node.second);
    //     //     world.send(to_rank, 2, sccEdges[scc_node.first]);
    //     // }

    //     changeSccLabels(sccs, sccNodes);
    //     constructMasterNode(edges, sccs);

    //     std::cout << "Done CONSTRUCT MASTER NODE\n" << std::endl;

    //     // wait for the completion of the process
    //     boost::mpi::request req[world.size()];
    //     STATUS status[world.size()];
    //     for (int i = 1; i < world.size(); i++)
    //     {
    //         req[i] = world.irecv(i, 0, status[i]);
    //     }
    //     boost::mpi::wait_all(req + 1, req + world.size());
    //     std::cout << "Done WAITING\n" << std::endl;
    // }
    // else
    //     processMessage();
}

void MaintainSCC::deleteEdges(std::vector<Edge> &decrement)
{
    boost::mpi::communicator world;
    int world_rank = world.rank();
    if (world_rank == 0)
    {
        for (auto &edge : decrement)
        {
            if (sccs[edge.from] == sccs[edge.to])
            {
                int to_rank = which_rank[edge.from];
                world.send(to_rank, 0, MessageType::EDGE_DELETE);
                world.send(to_rank, 1, edge);
            }
            else
            {
                deleteEdgeFromMaster(edge);
            }
        }
        std::cout << "Deleted Edges" << std::endl;
        for (int i = 1; i < world.size(); i++)
        {
            world.send(i, 0, MessageType::CLR_DEC_CACHE);
        }

        // waiting for signal of completion
        boost::mpi::request req[world.size()];
        STATUS status[world.size()];
        for (int i = 1; i < world.size(); i++)
        {
            req[i] = world.irecv(i, 0, status[i]);
        }
        boost::mpi::wait_all(req + 1, req + world.size());
        
        std::cout << "Done Delete-Cache Request" << std::endl;

        for (int i = 1; i < world.size(); i++)
        {
            if (status[i] == STATUS::DONE_NEW)
            {
                std::vector<Edge> edges;
                std::unordered_map<long int, long int> node_list;
                // std::vector<std::pair<long int, int>> transfer_list;
                world.recv(i, 1, edges);
                world.recv(i, 2, node_list);
                
                for (auto &node : node_list)
                {
                    sccs[node.first] = node.second;
                    int to_rank = node.second % (world.size() - 1) + 1;
                    which_rank[node.first] = to_rank;
                    // transfer_list.emplace_back(std::make_pair(node.first, to_rank));
                }
                updateMasterNode(0, edges, sccs);
                world.recv(i, 0, status[i]);
                std::cout << "Waited for New SCC completion" << std::endl;
            }
        }

        std::cout << "Done DELETE EDGES\n" << std::endl;
    }
}

void MaintainSCC::insertEdges(std::vector<Edge> &increament)
{
    boost::mpi::communicator world;
    int world_rank = world.rank();
    if (world_rank == 0)
    {
        std::vector<Edge> insert_in_master;
        for (auto &edge : increament)
        {
            if (sccs[edge.from] == sccs[edge.to])
            {
                int to_rank = which_rank[edge.from];
                world.send(to_rank, 0, MessageType::EDGE_INSERT);
                world.send(to_rank, 1, edge);
            }
            else
            {
                insert_in_master.push_back(edge);
            }
        }
        if (!insert_in_master.empty()) 
            insertEdgesInMaster(insert_in_master);
        for (int i = 0; i < world.size(); i++)
        {
            world.send(i, 0, MessageType::CLR_INC_CACHE);
        }

        // waiting for signal of completion
        boost::mpi::request req[world.size()];
        STATUS status[world.size()];
        for (int i = 1; i < world.size(); i++)
        {
            req[i] = world.irecv(i, 0, status[i]);
        }
        boost::mpi::wait_all(req + 1, req + world.size());
    }
}

void MaintainSCC::endAll()
{
    boost::mpi::communicator world;
    if (world.rank() == 0)
    {
        for (int i = 1; i < world.size(); i++)
        {
            world.send(i, 0, MessageType::EXIT);
        }
    }
}

MaintainSCC::~MaintainSCC()
{
    endAll();
}

}