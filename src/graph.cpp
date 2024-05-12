#include "graph.hpp"
#include "debug.hpp"

namespace SCC {

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
    std::unordered_map<long int, std::vector<long int>> adj;
    std::unordered_map<long int, std::vector<long int>> adj_rev;
    std::unordered_map<long int, int> visited;
    std::vector<long int> order;
    std::stack<long int> st;
    int mx = 0;

    for (Edge edge : edges) {
        adj[edge.from].push_back(edge.to);
        adj_rev[edge.to].push_back(edge.from);
    }
    std::function<void(void)> mark_unvisited = [&]() {
        for (const auto& [label, v] : sccs) {
            if (mx < label) mx = label;
            visited[label] = false;
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
    std::function<void(long int)> dfs1 = [&](long int u) {
        st.push(u);
        // 0 - not called, 1 - called and not finished, 2 - finished
        while(!st.empty()) {
            long int v = st.top();
            if(visited[v] > 0) {
                if(visited[v] != 2) order.push_back(v);
                visited[v] = 2;
                st.pop();
                continue;
            }
            visited[v] = 1;
            for (long int w : adj[v]) {
                if (!visited[w]) {
                    st.push(w);
                }
            }
        }
    };
    mark_unvisited();
    for (const auto& [label,v]: sccs) {
        if (!visited[label]) {
            dfs1(label);
        }
    }
    mark_unvisited();
    for (long int i = order.size() - 1; i >= 0; i--) {
        if (!visited[order[i]]) {
            dfs2(order[i]);
            mx++;
        }
    }
    
}

void TreeNode::removeEdge(const Edge &edge)
{
    //erase using key
    corresponds_to.erase(edge);
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
            // removing the split node
            if (edge.to < 0)
                edge.to = -edge.to;
            corresponds_to[edge] = e;
            // corresponds_to.emplace_back(std::make_pair(edge, e));
        }
    }
}

void TreeNode::checkUnreachable(std::unordered_set<long int> &unreachable)
{
    // get edge list
    std::vector<Edge> edges;
    long int split_on = -1;
    for (auto edge : corresponds_to)
    {
        if (edge.second.to < 0)
        {
            split_on = -edge.second.to;
            edge.second.to = split_on;
        }
        edges.emplace_back(edge.second);
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
    for (auto &correspond_to : corresponds_to)
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
    std::unordered_map<Edge, Edge, EdgeHash> temp;
    for(auto &correspond_to: corresponds_to)
    {
        // std::pair<Edge, Edge> &correspond_to = corresponds_to.back();
        // corresponds_to.pop_back();
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
            parent.corresponds_to[correspond_to.first] = correspond_to.second;
        }
        else
        {
            temp[correspond_to.first] = correspond_to.second;
        }
    }
    corresponds_to.swap(temp);
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
    for (auto &edge : edges)
    {
        if (edge.to == node)
        {
            edge.to = -node;
        }
    }
}

// Divide edges by groups based on SCCs (O(mlog(n)))
void MaintainSCC::divideEdgesByScc(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs, std::unordered_map<long int, std::vector<Edge>> &sccEdges, std::vector<Edge> &inter_scc_edges, bool split = false)
{
    for (auto &edge : edges)
    {
        if (sccs[edge.from] == sccs[edge.to])
        {
            sccEdges[sccs[edge.from]].push_back(edge);
        } else if(split) {
            inter_scc_edges.push_back(edge);
        }
    }
}

/*** Functions to make SCC Tree ***/
void MaintainSCC::makeSccTreeInternals(std::vector<Edge> &edge, TreeNode &currentNode)
{

    // std::cout << "Current Node: " << currentNode.label << std::endl;

    scc_tree_labels[currentNode.label] = currentNode.label;

    //get split on vertex (we choose vertex with highest degree (total))
    long int split_on = -1;
    std::unordered_map<long int, int> degree;
    // O(m)
    for (auto &e : edge)
    {
        degree[e.from]++;
        degree[e.to]++;
    }
    int max_degree = -1;
    for (auto &d : degree)
    {
        if (d.second > max_degree)
        {
            max_degree = d.second;
            split_on = d.first;
        }
    }

    // O(m)
    splitGraphOnNode(edge, split_on);

    // O(m + n)
    // actual node -> temp label
    std::unordered_map<long int, long int> sccs;
    for(auto &[from, to] : edge) {
        sccs[from] = from;
        sccs[to] = to;
    }

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
            makeSccTreeInternals(sccEdges[node.first], child);
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
    // Create a TreeNode
    int label = 0;
    roots.push_back(label);
    TreeNode &root = scc_tree_nodes[label];
    root.label = label;
    root.parent = label;
    root.dept = 0;
    for (const auto &edge : edges)
    {
        if (sccs[edge.from] != sccs[edge.to])
        {
            Edge e;
            e.from = sccs[edge.from];
            e.to = sccs[edge.to];
            root.corresponds_to[edge] = e;
            // root.corresponds_to.emplace_back(std::make_pair(edge, e));
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
    if (type == 0)
    {
        for (auto &correspond_to : scc_tree_nodes[0].corresponds_to)
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
        for (const auto &edge : edges)
        {
            if (sccs[edge.from] != sccs[edge.to])
            {
                Edge e;
                e.from = sccs[edge.from];
                e.to = sccs[edge.to];
                scc_tree_nodes[0].corresponds_to[edge] = e;
                // scc_tree_nodes[0].corresponds_to.emplace_back(std::make_pair(edge, e));
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
        for(auto &edge_pair: parent.corresponds_to) {
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
    for (auto &correspond_to : temp.corresponds_to)
    {
        edges.emplace_back(correspond_to.first);
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
    TreeNode &root = scc_tree_nodes[0];
    for(auto &edge : edges) {
        Edge label_edge(sccs[edge.from], sccs[edge.to]);
        root.corresponds_to[edge] = label_edge;
        // root.corresponds_to.emplace_back(std::make_pair(edge, label_edge));
    }
    edges.clear();
    // filling the edges for computing sccs
    for(const auto &correspond_to: root.corresponds_to)
        edges.emplace_back(correspond_to.second);
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
    std::unordered_map<long int, std::vector<std::pair<const Edge, Edge>>> new_scc_root;
    std::unordered_map<Edge, Edge, EdgeHash> temp;
    for(auto &correspond_to: root.corresponds_to) {
        // std::pair<Edge, Edge> &correspond_to = root.corresponds_to.end();
        // root.corresponds_to.pop_back();
        actual_node[correspond_to.second.from] = correspond_to.first.from;
        actual_node[correspond_to.second.to] = correspond_to.first.to;
        if(sccs[correspond_to.second.from] == sccs[correspond_to.second.to]) {
            new_scc_root[sccs[correspond_to.second.from]].emplace_back(correspond_to);
        } else {
            temp.emplace_hint(temp.end(), correspond_to.first, correspond_to.second);
        }
    }
    root.corresponds_to.swap(temp);
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
    // add the edge to the lca
    lca.corresponds_to[edge] = e;
    // add the lca to the insert cache
    insert_cache.insert({lca.dept, lca.label});
}

void MaintainSCC::clearInsertCache(TreeNode &node)
{
    std::vector<Edge> edges;
    for (auto &correspond_to : node.corresponds_to)
        edges.emplace_back(correspond_to.second);

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

    std::unordered_map<Edge, Edge, EdgeHash> temp;
    for(auto &correspond_to: node.corresponds_to)
    {
        // std::pair<Edge, Edge> &correspond_to = node.corresponds_to.back();
        // node.corresponds_to.pop_back();
        int scc_label_from = sccs[correspond_to.second.from], scc_label_to = sccs[correspond_to.second.to];
        if (scc_label_from != scc_label_to)
            temp[correspond_to.first] = correspond_to.second;
        else {
            if(correspond_to.second.to == split_on[scc_label_to])
                correspond_to.second.to = -correspond_to.second.to;
            scc_tree_nodes[sccs[correspond_to.second.from]].corresponds_to[correspond_to.first] = correspond_to.second;
        }
    }
    node.corresponds_to.swap(temp);

    // changing connections in node
    for (auto &correspond_to : node.corresponds_to)
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
        else if (type == MessageType::SCC_TREE)
        {
            std::vector<long int> nodes;
            std::vector<Edge> edgeList;
            world.recv(0, 1, SCC_LABEL);
            world.recv(0, 2, nodes);
            world.recv(0, 3, edgeList);
            makeSccTree(edgeList, nodes);
        }
        else if (type == MessageType::EDGE_DELETE)
        {
            Edge edge;
            world.recv(0, 1, edge);
            deleteEdge(edge);
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
    world_size = world.size();

    if (world_rank == 0)
    {
        SCC_LABEL = n + 1;
        sccs.reserve(n+1);
        for (long int i = 1; i <= n; i++)
            sccs[i] = i;
        findScc(edges, sccs);

        // std::cout << "Done FIND SCC\n" << std::endl;

        std::unordered_map<long int, std::vector<Edge>> sccEdges;
        std::vector<Edge> inter_scc_edges;
        divideEdgesByScc(edges, sccs, sccEdges, inter_scc_edges, true);

        // std::cout << "Done DIVIDE EDGES\n" << std::endl;

        std::unordered_map<long int, std::vector<long int>> sccNodes;
        for (auto &scc : sccs)
        {
            sccNodes[scc.second].push_back(scc.first);
            which_rank[scc.first] = scc.second % (world.size() - 1) + 1;
        }

        // std::cout << "Done DIVIDE NODES\n" << std::endl;

        // std::cout << "Number of SCCs: " << sccNodes.size() << std::endl;
        // std::cout << "Set of Edge: " << sccEdges.size() << std::endl;


        for(auto &scc_node: sccNodes) {
            int to_rank = which_rank[scc_node.second[0]];
            world.send(to_rank, 0, MessageType::SCC_TREE);
            world.send(to_rank, 1, SCC_LABEL + to_rank);
            world.send(to_rank, 2, scc_node.second);
            world.send(to_rank, 3, sccEdges[scc_node.first]);
        }

        changeSccLabels(sccs, sccNodes);
        constructMasterNode(edges, sccs);

        // std::cout << "Done CONSTRUCT MASTER NODE\n" << std::endl;
    }
    else
        processMessage();
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
        for (int i = 0; i < world.size(); i++)
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

        for (int i = 1; i < world.size(); i++)
        {
            if (status[i] == STATUS::DONE_NEW)
            {
                std::vector<Edge> edges;
                std::unordered_map<long int, long int> node_list;
                // std::vector<std::pair<long int, int>> transfer_list;
                world.recv(i, 1, edges);
                world.recv(i, 2, node_list);
                for (auto edge : edges)
                {
                    dEdge(edge);
                }
                for (auto &node : node_list)
                {
                    sccs[node.first] = node.second;
                    int to_rank = node.second % (world.size() - 1) + 1;
                    which_rank[node.first] = to_rank;
                    // transfer_list.emplace_back(std::make_pair(node.first, to_rank));
                }
                updateMasterNode(0, edges, sccs);
                // world.send(i, 0, STATUS::TRANSFER);
                // world.send(i, 1, transfer_list);
                dTreeNode(scc_tree_nodes[0]);
            }
        }

        // std::cout << "Done DELETE EDGES\n" << std::endl;
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