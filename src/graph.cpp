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
                edges[i].to = -edges[i].to;
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
    std::unordered_map<long int, long int> nwsccs;
    for (auto &node : this->contains)
        nwsccs[node] = node;
    findScc(edges, nwsccs);
    int split_on_label = nwsccs[split_on];
    for (const auto &node : this->contains)
    {
        if (nwsccs[node] != split_on_label)
            unreachable.insert(node);
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
void MaintainSCC::traverseNode(int node, int to_set, std::unordered_map<long int, long int> &update)
{
    int max_threads = omp_get_max_threads();
    std::vector<int> queue;
    std::unordered_map<int,bool> in_next_queue;
    std::vector<std::vector<int>> thread_queue(max_threads);
    queue.push_back(node);
    in_next_queue[node] = true;
    while(!queue.empty()) {
        #pragma omp parallel for num_threads(max_threads)
        for(int v: queue) {
            int tid = omp_get_thread_num();
            update[v] = to_set;
            for(auto &child : scc_tree_nodes[v].contains) {
                if(!in_next_queue[child]) {
                    thread_queue[tid].push_back(child);
                    in_next_queue[child] = true;
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
            currentNode.contains.insert(node.first);
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
void MaintainSCC::makeSccTree(std::vector<Edge> &edges, std::vector<long int> &nodes, int label)
{
    // case when only one node is present in the scc
    if (nodes.size() == 1)
    {
        roots.push_back(label);
        TreeNode &root = scc_tree_nodes[nodes[0]];
        root.label = nodes[0];
        root.parent = nodes[0];
        root.dept = 0;
        return;
    }
    // Create a TreeNode
    roots.push_back(label);
    TreeNode &root = scc_tree_nodes[label];
    root.label = label;
    root.parent = label;
    root.dept = 0;
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
    root.dept = -1;
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
void MaintainSCC::updateMasterNode()
{
    int max_threads = omp_get_max_threads();
    //re-run the scc algorithm on the master node

    //1 - get the edges from the master node
    std::vector<Edge> edges;
    std::vector<int> prev_edges_size(max_threads, 0);
    int edges_size = 0;
    for(int i = 1; i<max_threads; i++) {
        prev_edges_size[i] = prev_edges_size[i-1] + scc_tree_nodes[0].corresponds_to[i-1].size();
    }
    edges_size = prev_edges_size[max_threads-1] + scc_tree_nodes[0].corresponds_to[max_threads-1].size();
    edges.resize(edges_size);
    #pragma omp parallel for num_threads(max_threads)
    for(int i = 0; i<max_threads; i++) {
        int j = prev_edges_size[i];
        for(auto &correspond_to: scc_tree_nodes[0].corresponds_to[i]) {
            edges[j++] = correspond_to.second;
        }
    }

    //2 - get the nodes from the master node
    std::unordered_map<long int, long int> new_sccs;
    for(auto &node: scc_tree_nodes[0].contains) {
        new_sccs[node] = node;
    }

    //3 - run the scc algorithm
    findScc(edges, new_sccs);

    //4 - divide the nodes based on the new sccs formed
    std::unordered_map<long int, std::vector<long int>> sccNodes;
    for(auto &scc: new_sccs) {
        sccNodes[scc.second].push_back(scc.first);
    }

    //5 - change the sccs mappings to labels
    changeSccLabels(new_sccs, sccNodes);

    std::unordered_map<long int, std::vector<long int>> init_sccNodes;
    for(auto &scc: sccs) {
        init_sccNodes[scc.second].push_back(scc.first);
    }

    std::unordered_map<long int, long int> updated_labels;
    for(auto &scc: new_sccs) {
        if(init_sccNodes.find(scc.second) == init_sccNodes.end()) {
            traverseNode(scc.first, scc.second, updated_labels);
            for(auto &nd: init_sccNodes[scc.first]) {
                updated_labels[nd] = scc.second;
            }
        }
    }

    for(auto &scc: updated_labels) {
        if(sccs.find(scc.first) != sccs.end()) {
            sccs[scc.first] = scc.second;
        }
    }
}

int MaintainSCC::checkAndRemoveUnreachable(TreeNode &curr_node)
{
    int max_threads = omp_get_max_threads();
    // do unreachability and reachability checks O(n+m)
    std::unordered_set<long int> unreachable;
    curr_node.checkUnreachable(unreachable);

    // if no unreachable nodes then return
    if (unreachable.empty())
        return 0;

    // remove all unreachable nodes from the current node O(n)
    for (auto &unreachable_node : unreachable)
        curr_node.contains.erase(unreachable_node);

    // add the unreachable nodes to the parent
    int parent_label = curr_node.parent;
    if(curr_node.dept == 0) parent_label = 0;
    TreeNode &parent = scc_tree_nodes.at(parent_label);
    // change parent label of unreachable nodes O(n)
    for (auto &unreachable_node : unreachable)
    {
        scc_tree_nodes[unreachable_node].parent = parent.label;
        scc_tree_labels[unreachable_node] = parent.label;
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
        scc_tree_labels[new_label] = parent.label;
        parent.contains.insert(new_label);
    }

    for(const auto &nd: unreachable) {
        traverseNode(nd, nd, new_labels);
    }

    // update the parent with the new labels
    parent.updateLabels(new_labels);

    // exposing internal structure to parent node
    curr_node.exposeToParent(parent, unreachable);

    return 1;
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


void MaintainSCC::deleteEdgeFromMaster(Edge edge)
{
    TreeNode &root = scc_tree_nodes[0];
    root.removeEdge(edge);
}

void MaintainSCC::insertEdgeInMaster(Edge edge)
{
    TreeNode &root = scc_tree_nodes[0];
    Edge label_edge(sccs[edge.from], sccs[edge.to]);
    int min_thread = 0;
    int min_thread_size = INT_MAX;
    for (int i = 0; i < omp_get_max_threads(); i++)
    {
        if (root.corresponds_to[i].size() < min_thread_size)
        {
            min_thread = i;
            min_thread_size = root.corresponds_to[i].size();
        }
    }
    root.corresponds_to[min_thread][edge] = label_edge;
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


//clear the delete cache
bool MaintainSCC::clearDeleteCache()
{
    bool new_scc = false;
    int max_threads = omp_get_max_threads();

    while (!delete_cache.empty())
    {
        Cache c = *delete_cache.begin();
        delete_cache.erase(delete_cache.begin());
        TreeNode &node = scc_tree_nodes.at(c.label);
        dTreeNode(node);
        if (checkAndRemoveUnreachable(node)) {
            if(node.dept != 0)
                delete_cache.insert({node.dept - 1, node.parent});
            else
                new_scc = true;
        }
    }

    return new_scc;
}

// query if two nodes are in the same SCC
bool MaintainSCC::query(long int v1, long int v2)
{
    boost::mpi::communicator world;
    bool ret = (sccs[v1] == sccs[v2]);
    bool ans;
    boost::mpi::reduce(world, ret, ans, std::logical_and<bool>(), 0);
    boost::mpi::broadcast(world, ans, 0);
    return ans;
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

    int temp = 0;
    //check and distribute the sccs
    for(auto it = sccEdges.begin(); it != sccEdges.end(); it++) {
        int scc_size = it->second.size();
        if(scc_size > 0.1 * m) {
            int to_rank = temp % (world_size - 1) + 1;
            which_rank[it->first] = to_rank;
            if(to_rank == world_rank){
                makeSccTree(it->second, sccNodes[it->first], it->first);
            }
            temp++;
        } else {
            int to_rank = 0;
            which_rank[it->first] = to_rank;
            if(to_rank == world_rank)
                makeSccTree(it->second, sccNodes[it->first], it->first);
        }
    }

    //maitained by every node
    constructMasterNode(inter_scc_edges, sccs);
    world.barrier();
    
    for(auto root: roots) {
        if(root == 0) continue;
        dSccTree(root, scc_tree_nodes);
    }

    world.barrier();

    dTreeNode(scc_tree_nodes[0]);

    world.barrier();
}

void MaintainSCC::deleteEdges(std::vector<Edge> &decrement)
{
    boost::mpi::communicator world;
    int world_rank = world.rank();
    int max_threads = omp_get_max_threads();

    //assuming every core has access to decrement edges
    for(int i = 0; i<decrement.size(); i++) {
        if(sccs[decrement[i].from] == sccs[decrement[i].to]) {
            int scc_label = sccs[decrement[i].from];
            int to_rank = which_rank[scc_label];
            if(to_rank == world_rank) {
                deleteEdge(decrement[i]);
            }
        } else {
            deleteEdgeFromMaster(decrement[i]);
        }
    }

    bool new_scc = clearDeleteCache();
    
    if(new_scc) {
        updateMasterNode();
    }

    world.barrier();
}

void MaintainSCC::insertEdges(std::vector<Edge> &increament)
{
    boost::mpi::communicator world;
    int world_rank = world.rank();

    for(int i = 0; i<increament.size(); i++) {
        if(sccs[increament[i].from] == sccs[increament[i].to]) {
            int scc_label = sccs[increament[i].from];
            int to_rank = which_rank[scc_label];
            if(to_rank == world_rank) {
                insertEdge(increament[i]);
            }
        } else {
            insertEdgeInMaster(increament[i]);
        }
    }

    updateMasterNode();
    world.barrier();
}

int MaintainSCC::getNumberOfSCCs()
{
    boost::mpi::communicator world;
    int n = scc_tree_nodes[0].contains.size();
    int ret = -1;
    // find the max number of sccs
    boost::mpi::reduce(world, n, ret, boost::mpi::maximum<int>(), 0);
    boost::mpi::broadcast(world, ret, 0);
    return ret;
}

}