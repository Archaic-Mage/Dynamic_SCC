#include "graph.hpp"
#include "debug.hpp"

int MaintainSCC::getLabel()
{
    int ret = SCC_LABEL;
    SCC_LABEL += world_size;
    return ret;
}

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
            // removing the split node
            if (edge.to < 0)
                edge.to = -edge.to;
            corresponds_to.emplace_back(std::make_pair(edge, e));
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

void TreeNode::getNewLabels(std::unordered_set<long int> &unreachable, std::unordered_map<long int, long int> &new_labels)
{
    for (auto &correspond_to : corresponds_to)
    {
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

void TreeNode::updateLabels(std::unordered_map<long int, long int> &new_labels)
{
    for (auto &correspond_to : corresponds_to)
    {
        Edge &label_edge = correspond_to.second;
        Edge &actual_edge = correspond_to.first;

        if ( new_labels.find ( actual_edge.from ) != new_labels.end() )
            label_edge.from = new_labels[actual_edge.from];
        if ( new_labels.find ( actual_edge.to ) != new_labels.end() )
            label_edge.to = new_labels[actual_edge.to];
    }
}

void TreeNode::exposeToParent(TreeNode &parent, std::unordered_set<long int> unreachable)
{
    std::vector<std::pair<Edge, Edge>> temp;
    while (!corresponds_to.empty())
    {
        std::pair<Edge, Edge> &correspond_to = corresponds_to.back();
        corresponds_to.pop_back();
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
            parent.corresponds_to.push_back(correspond_to);
        }
        else
        {
            temp.push_back(correspond_to);
        }
    }
    corresponds_to = temp;
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
void MaintainSCC::divideEdgesByScc(std::vector<Edge> &edges, std::unordered_map<long int, long int> &sccs, std::unordered_map<long int, std::vector<Edge>> &sccEdges)
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
void MaintainSCC::makeSccTreeInternals(std::vector<Edge> &edge, TreeNode &currentNode)
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
            root.corresponds_to.emplace_back(std::make_pair(edge, e));
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
            Edge &actual_edge = correspond_to.first;
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
                scc_tree_nodes[0].corresponds_to.emplace_back(std::make_pair(edge, e));
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
    {
        return 0;
    }
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
    TreeNode &lca = scc_tree_nodes.at(getLCANode(edge.from, edge.to));
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
    lca.corresponds_to.emplace_back(std::make_pair(edge, e));
    // add the lca to the insert cache
    insert_cache.insert({lca.dept, lca.label});
}

void MaintainSCC::clearInsertCache(TreeNode &node)
{
    std::vector<Edge> edges;
    for (auto &correspond_to : node.corresponds_to)
    {
        edges.emplace_back(correspond_to.second);
    }
    std::unordered_map<long int, long int> sccs;
    findScc(edges, sccs);

    std::unordered_map<long int, std::vector<long int>> sccNodes;
    for (auto &scc : sccs)
    {
        sccNodes[scc.second].push_back(scc.first);
    }

    bool to_end = true;
    for (auto &scc_nodes : sccNodes)
    {
        if (scc_nodes.second.size() > 1)
        {
            to_end = false;
        }
    }
    if (to_end)
        return;

    changeSccLabels(sccs, sccNodes);

    sccNodes.clear();
    for (auto &scc : sccs)
    {
        sccNodes[scc.second].push_back(scc.first);
    }

    std::vector<long int> new_tree_nodes;
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
            node.contains.insert(sccNodes[scc.second].begin(), sccNodes[scc.second].end());
            new_tree_nodes.push_back(child.label);
        }
    }

    std::vector<std::pair<Edge, Edge>> temp;
    while (!node.corresponds_to.empty())
    {
        std::pair<Edge, Edge> &correspond_to = node.corresponds_to.back();
        node.corresponds_to.pop_back();
        if (sccs[correspond_to.second.from] != sccs[correspond_to.second.to])
        {
            temp.emplace_back(correspond_to);
        }
        else
        {
            scc_tree_nodes[sccs[correspond_to.second.from]].corresponds_to.emplace_back(correspond_to);
        }
    }
    node.corresponds_to = temp;

    // changing connections in node
    for (auto &correspond_to : node.corresponds_to)
    {
        Edge &label_edge = correspond_to.second;
        Edge &actual_edge = correspond_to.first;
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
    dTreeNode(node);
    for (auto &new_tree_node : new_tree_nodes)
    {
        TreeNode &child = scc_tree_nodes[new_tree_node];
        dTreeNode(child);
        // insertCache.insert({child.dept, child.label});
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
            for (auto &root : roots)
            {
                dSccTree(root, scc_tree_nodes);
            }
        }
        else if (type == MessageType::EDGE_DELETE)
        {
            Edge edge;
            world.recv(0, 1, edge);
            deleteEdge(edge);
            std::cout << "Deleted Edge: " << edge.from << " " << edge.to << std::endl;
        }
        else if (type == MessageType::CLR_DEC_CACHE)
        {
            while (!delete_cache.empty())
            {
                dInfo(world, "Deleting Cache");
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
        }
        else if (type == MessageType::EDGE_INSERT)
        {
            Edge edge;
            world.recv(0, 1, edge);
            std::cout << "Inserted Edge: " << edge.from << " " << edge.to << std::endl;
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
        for (long int i = 1; i <= n; i++)
            sccs[i] = i;
        findScc(edges, sccs);

        std::unordered_map<long int, std::vector<Edge>> sccEdges;
        divideEdgesByScc(edges, sccs, sccEdges);

        std::unordered_map<long int, std::vector<long int>> sccNodes;
        for (auto &scc : sccs)
        {
            sccNodes[scc.second].push_back(scc.first);
            which_rank[scc.first] = scc.second % (world.size() - 1) + 1;
        }

        for (auto &scc_edges : sccEdges)
        {
            std::vector<long int> &nodes = sccNodes[scc_edges.first];
            int to_rank = which_rank[nodes[0]];
            world.send(to_rank, 0, MessageType::SCC_TREE);
            world.send(to_rank, 1, SCC_LABEL + to_rank);
            world.send(to_rank, 2, nodes);
            world.send(to_rank, 3, scc_edges.second);
        }

        changeSccLabels(sccs, sccNodes);
        constructMasterNode(edges, sccs);
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
    }
}

void MaintainSCC::insertEdges(std::vector<Edge> &increament)
{
    boost::mpi::communicator world;
    int world_rank = world.rank();
    if (world_rank == 0)
    {
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
                // insert edge in master node
            }
        }
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