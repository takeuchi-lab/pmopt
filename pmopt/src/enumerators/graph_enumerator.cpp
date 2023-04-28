#include "pmopt/enumerators/graph_enumerator.hpp"

#include <iostream>
#include "pmopt/utils/containers.hpp"
#include "pmopt/logger.hpp"

using namespace PMOpt;

static const unsigned int DUMMY = 0;
static const unsigned int LOG = 100;
static const std::string PREFIX = "  - ";


unsigned int GraphEnumerator::enumerate(std::vector<PDB> & children, const PDB & parent) noexcept
{
    _edge_enumerator.initialize(children);

    if (parent._length == 0)
    {
        // if parent is ROOT pattern then enumerating single vertex pattern
        for (unsigned int i_instance = 0; i_instance < _database.size(); ++i_instance)
            search_vertex(_edge_enumerator, i_instance);
    }
    else 
    {
        // if parent is not ROOT
        Logger logger(std::cout);

        // Build rightmost path
        // In gspan, edge is expanded from the vertices on the `rightmost path`.
        // Here, the ids of vertices and edges in pattern will be recoded.
        unsigned int backward_start_index;
        build_rightmost_path(_rightmost_path_pattern, backward_start_index, parent._pattern.get());

        logger(LOG, PREFIX, "RMP_P: ", Utils::to_string(_rightmost_path_pattern));

        // For each instance in which the parent occurs
        auto support = parent._occurrence.size();
        for (unsigned int i_occ = 0; i_occ < support; ++i_occ)
        {
            // For each position of instance in which the parent occurs
            auto i_instance = parent._occurrence[i_occ];
            const auto & graph = _database[i_instance];
            for (const auto & position : parent._positions[i_occ])
            {
                // initialize visited flag for edges and vertices
                initialize_visited(graph);

                // Trace the built rightmost path.
                // Here, the ids of rightmost edges and rightmost vertices in instance will be recorded.
                trace_rightmost_path(
                    _rightmost_path_instance,
                    _vertex_is_visited,
                    _edge_is_visited,
                    _rightmost_path_pattern,
                    position.get(),
                    parent._pattern.get()
                );


                logger(LOG, PREFIX, "RMP_S: ", Utils::to_string(_rightmost_path_instance));


                // Search edges
                search_backward_edge(_edge_enumerator, position, i_instance, backward_start_index);
                search_forward_edge(_edge_enumerator, position, i_instance);
            }
        }
    }


    // Postprocess for enumerated edges
    for (const auto & edge_and_child_index : _edge_enumerator._edge2index)
    {
        auto & child = children[edge_and_child_index.second];

        // If parent is ROOT (i.e. child is single vertex) then edge_id is also ROOT
        unsigned int edge_id;
        if (parent._length == 0)
            edge_id = 0;
        else
            edge_id = parent._pattern->_symbol._edge_pattern._edge_id + 1;

        // set length and pattern
        child._length = parent._length + 1;
        child._pattern = std::make_shared<const Pattern<Symbol>>(Pattern<Symbol>{
            Symbol{
                child._occurrence[0],
                child._positions[0][0]->_symbol,
                edge_and_child_index.first.edge(edge_id),
            },
            parent._pattern
        });

        // check dfs-code minimality
        // single vertex pattern obviously has minimum dfs-code
        if (parent._length != 0 && !check_dfscode(child))
            continue;
        
        // Passed dfs-code check
        child._is_valid = true;
    }

    return _edge_enumerator._n_children;
}




void GraphEnumerator::build_rightmost_path(
    std::vector<Edge> & rightmost_path,
    unsigned int & backward_start_index,
    const Pattern<Symbol> * pattern
)
const noexcept
{
    // reset rightmost path
    rightmost_path.clear();

    // set initial position
    auto start = pattern->_symbol._edge_pattern._start_vertex_id;
    auto end = pattern->_symbol._edge_pattern._end_vertex_id;
    auto rightmost_vertex = std::max(start, end);
    auto end_of_rightmost_edge = end;

    // determine from which index to start the search of backward-edge
    // whether rightmost edge is forward-edge or not
    // if true then searching backward-edge will be started from the end of rightmost path
    // else searching backward-edge will be started from the end of forward-edge in rightmost path
    auto rightmost_edge_is_forward = pattern->_symbol._edge_pattern.is_forward();

    // this will store index of backward-edge where the search should start
    backward_start_index = 0;

    // iterate for each prefix
    // Do not add ROOT edge because there is no ROOT edge in instance
    while (pattern->_parent)
    {
        // update starting vertex and ending one
        start = pattern->_symbol._edge_pattern._start_vertex_id;
        end = pattern->_symbol._edge_pattern._end_vertex_id;

        auto is_forward = pattern->_symbol._edge_pattern.is_forward();
        if (is_forward && rightmost_vertex == end)
        {
            // if the edge of the current pos is forward-edge
            // and it is connected to the ending vertex of the current rightmost edge
            // then, add the edge into the rightmost path
            rightmost_path.push_back(pattern->_symbol._edge_pattern);

            // update rightmost vertex
            rightmost_vertex = start;

            // set backward start index
            if (!rightmost_edge_is_forward && end_of_rightmost_edge == start)
            {
                backward_start_index = rightmost_path.size() + 1;
            }
        }

        // Next
        pattern = pattern->_parent.get();
    }

    // if the rightmost edge is forward
    // then search of backward-edge will be started from end of the edge of rightmost path
    if (rightmost_edge_is_forward)
    {
        backward_start_index = rightmost_path.size() - 1;
    }
}




void GraphEnumerator::trace_rightmost_path(
    std::vector<Edge> & rightmost_path_instance,
    std::vector<unsigned int> & vertex_is_visited,
    std::vector<unsigned int> & edge_is_visited,
    const std::vector<Edge> & rightmost_path_pattern,
    const Pattern<Edge> * position,
    const Pattern<Symbol> * pattern
)
const noexcept
{
    // Note
    // The rightmost path in pattern and right most path in instance is different
    // because vertex and edge id numbering is different between in pattern and instance
    // even if the pointing vertex is the same one.
    rightmost_path_instance.clear();

    // Trace rightmost path according to rightmost path in pattern.
    unsigned int index_rightmost_pattern = 0;
    while (position != nullptr)
    {
        // Check visited edges and vertices
        vertex_is_visited[position->_symbol._start_vertex_id] = true;
        vertex_is_visited[position->_symbol._end_vertex_id] = true;
        edge_is_visited[position->_symbol._edge_id] = true;

        // Check consistency between the edge in built rightmost path
        // and the edge found in tracing position in the pattern.
        if (pattern->_symbol._edge_pattern.is_forward()
            && rightmost_path_pattern[index_rightmost_pattern] == pattern->_symbol._edge_pattern)
        {
            rightmost_path_instance.push_back(position->_symbol);
            ++index_rightmost_pattern;
        }

        // finish tracing
        if (index_rightmost_pattern >= rightmost_path_pattern.size())
            break;

        // Next
        position = position->_parent.get();
        pattern = pattern->_parent.get();
    }

    assert(rightmost_path_pattern.size() == rightmost_path_instance.size());
}




bool GraphEnumerator::check_dfscode(const PDB & projected) noexcept
{
    // If the length is 0 or 1, the pattern has obviously minimum dfs-code.
    assert(projected._length > 1);

    const auto & pattern = * projected._pattern.get();
    const auto & instance = _database[pattern._symbol._instance_index];

    // If the length is 2, checking is done simply comparing vertex labels.
    if (projected._length == 2)
    {
        const auto & edge = pattern._symbol._edge_instance;
        return instance.vertex_label(edge._start_vertex_id)
            <= instance.vertex_label(edge._end_vertex_id);
    }
    
    // Create dfs-code of the pattern
    build_dfscode(_dfscode, pattern);

    // Create Graph-typed pattern of the pattern
    build_pattern_graph(_dfs_check_graph, pattern);

    // Initialize visited flags of edges and vertices
    initialize_visited(_dfs_check_graph);

    // For single-vertex edges
    // Initialize checker
    unsigned int edge_id = 0;
    _dfscode_checker.initialize(_dfscode[edge_id]);
    if (!search_vertex(_dfscode_checker, DUMMY))
        return false;

    // For non-single-vertex edges
    for (edge_id = 1; edge_id < _dfscode.size(); ++edge_id)
    {
        // Initialize checker
        _dfscode_checker.initialize(_dfscode[edge_id]);

        // TODO here need to set pattern
        // Build rightmost path
        unsigned int backward_start_index;
        build_rightmost_path(_rightmost_path_pattern, backward_start_index, nullptr);

        // For each positions
        for (auto position : _dfscode_checker._current_positions)
        {
            initialize_visited(_dfs_check_graph);

            // Trace rightmost path
            trace_rightmost_path(
                _rightmost_path_instance,
                _vertex_is_visited,
                _edge_is_visited,
                _rightmost_path_pattern,
                position.get(),
                nullptr // TODO
            );

            if (!search_backward_edge(_dfscode_checker, position, DUMMY, backward_start_index))
                return false;

            if (!search_forward_edge(_dfscode_checker, position, DUMMY))
                return false;
        }
    }

    return true;
}




void GraphEnumerator::build_dfscode(
    std::vector<DFSEdge> & dfscode,
    const Pattern<Symbol> & pattern
)
const noexcept
{
    const auto * p = & pattern;
    const Graph * instance;

    // # of edges except ROOT edge
    // The ROOT and single-vertex edge have id = 0
    // and the subsequent edges have consecutive ids,  like id = 1,2,3,...
    // So # of edges with the single-vertex one is maximun edge id + 1
    auto n_edges = p->_symbol._edge_pattern._edge_id + 1;
    assert(n_edges >= 2);

    dfscode.clear();
    dfscode.resize(n_edges);

    // Here add the edges iterativly (except ROOT edge and single-vertex edge)
    // edge_id = 0,1,...,n_edges-1
    // edge_id = 0 if the edge is a single-vertex edge
    unsigned int edge_id = n_edges - 1;
    const Edge * edge_pattern;
    const Edge * edge_instance;
    for (; edge_id > GraphDatabase::ROOT; --edge_id)
    {
        edge_pattern = & p->_symbol._edge_pattern;
        edge_instance = & p->_symbol._edge_instance;
        instance = & _database[p->_symbol._instance_index];

        // Create dfs-edge from the edge of pattern
        dfscode[edge_id] = DFSEdge{
            edge_pattern->_start_vertex_id,
            edge_pattern->_end_vertex_id,
            instance->vertex_label(edge_instance->_start_vertex_id),
            instance->edge_label(edge_instance->_edge_id),
            instance->vertex_label(edge_instance->_end_vertex_id)
        };

        // Next
        p = p->_parent.get();
    }

    // Here add the single-veretx edge
    edge_pattern = & p->_symbol._edge_pattern;
    edge_instance = & p->_symbol._edge_instance;
    instance = & _database[p->_symbol._instance_index];

    assert(edge_id == GraphDatabase::ROOT);
    assert(edge_pattern->_start_vertex_id == GraphDatabase::ROOT);
    assert(edge_pattern->_end_vertex_id == GraphDatabase::ROOT + 1);
    assert(edge_instance->_start_vertex_id == GraphDatabase::ROOT);

    // Single-vertex edge
    dfscode[0] = DFSEdge{
        GraphDatabase::ROOT,
        GraphDatabase::ROOT + 1,
        GraphDatabase::ROOT,
        GraphDatabase::ROOT,
        instance->vertex_label(edge_instance->_end_vertex_id)
    };
}




void GraphEnumerator::build_pattern_graph(
    Graph & graph,
    const Pattern<Symbol> & pattern
)
const noexcept
{
    const Graph * instance;
    const auto * p = & pattern;
    const auto * edge_pattern = & p->_symbol._edge_pattern;
    const auto * edge_instance = & p->_symbol._edge_instance;

    // # of vertices and edges
    // Here, the ROOT vertex, ROOT edge are excluded.
    auto n_vertices = std::max(edge_pattern->_start_vertex_id, edge_pattern->_end_vertex_id);
    auto n_edges = edge_pattern->_edge_id + 1;
    assert(n_vertices >= 1);
    assert(n_edges >= 1);

    // Initialize graph
    // ROOT and single-vertex are excluded
    graph.clear();
    graph.resize(n_vertices);
    graph._vertex_labels.resize(n_vertices);
    graph._edge_labels.resize(n_edges - 1); // minus one to exclude a single-vertex edge

    // Here add the edges iterativly (except a single-vertex edge)
    // edge_id = 0,1,...,n_edges-1
    // edge_id = 0 if the edge is a single-vertex edge
    unsigned int start_id_pattern;
    unsigned int end_id_pattern;
    auto edge_id = n_edges - 1;
    for (; edge_id > GraphDatabase::ROOT; --edge_id)
    {
        instance = & _database[p->_symbol._instance_index];
        edge_pattern = & p->_symbol._edge_pattern;
        edge_instance = & p->_symbol._edge_instance;
        start_id_pattern = edge_pattern->_start_vertex_id;
        end_id_pattern = edge_pattern->_end_vertex_id;

        // Add edges and its label
        assert(start_id_pattern != GraphDatabase::ROOT);
        assert(end_id_pattern != GraphDatabase::ROOT);

        // Add edges
        Edge edge{start_id_pattern, edge_id, end_id_pattern};
        if (start_id_pattern != end_id_pattern)
            graph[end_id_pattern].emplace_back(edge.reverse());
        graph[start_id_pattern].emplace_back(std::move(edge));

        // Add edge label
        graph.edge_label(edge_id) = instance->edge_label(edge_instance->_edge_id);

        // Add vertex label
        // Only add the ending vertex
        // because each vertex appears more than once at ecah of the starting and ending ones.
        graph.vertex_label(end_id_pattern) = instance->vertex_label(edge_instance->_end_vertex_id);

        // Next
        p = p->_parent.get();
    }

    instance = & _database[p->_symbol._instance_index];
    edge_pattern = & p->_symbol._edge_pattern;
    edge_instance = & p->_symbol._edge_instance;
    start_id_pattern = edge_pattern->_start_vertex_id;
    end_id_pattern = edge_pattern->_end_vertex_id;
    assert(start_id_pattern == GraphDatabase::ROOT);
    assert(end_id_pattern == GraphDatabase::ROOT + 1);

    // Add the first starting vertex (the ending of the single-vertex edge)
    graph.vertex_label(end_id_pattern) = instance->vertex_label(edge_instance->_end_vertex_id);
}


void GraphEnumerator::DFSCodeChecker::initialize(
    const DFSEdge & dfs_edge
)
noexcept
{
    _current_edge = dfs_edge;
    std::swap(_next_positions, _current_positions);
    _next_positions.clear();
}


bool GraphEnumerator::DFSCodeChecker::push(
    const DFSEdge & dfs_edge,
    const Edge & edge_instance,
    std::shared_ptr<const Pattern<Edge>> position,
    unsigned int // unused
)
noexcept
{
    // the case that there is smaller dfs edge than the current edge
    // this means that the current edge is not minimum dfs-code edge.
    // so it is not needed to check further dfs-code with respect to the current edge
    if (dfs_edge < _current_edge)
        return false;

    // the case that there is larger dfs edge than the current edge
    // this dfs edge is not achieve minimum dfs-code
    // so it is needed to search more to find the edge achieving minimum dfs-code
    if (_current_edge < dfs_edge)
        return true;

    
    // the case that there is an edge with same dfs-code order.
    // it is needed further dfs-code checking in this sub-graph
    _next_positions.push_back(
        std::make_shared<const Pattern<Edge>>(Pattern<Edge>{edge_instance,  position})
    );

    // continue searching
    return true;
}




void GraphEnumerator::EdgeEnumerator::initialize(std::vector<PDB> & children) noexcept
{
    _children = & children;
    _n_children = 0;
    _edge2index.clear();
}




bool GraphEnumerator::EdgeEnumerator::push(
    const DFSEdge & dfs_edge,
    const Edge & edge_instance,
    std::shared_ptr<const Pattern<Edge>> position,
    unsigned int instance_index
)
noexcept
{
    unsigned int index;
    PDB * child;

    // The case that new edge found
    if (_edge2index.find(dfs_edge) == _edge2index.end())
    {
        // register dfs edge
        _edge2index[dfs_edge] = _n_children;
        index = _edge2index[dfs_edge];

        // create new child
        ++_n_children;
        _children->resize(std::max((std::size_t) _n_children, _children->size()));
        child = & (* _children)[index];
        child->_is_valid = false;
        child->_occurrence.clear();
        child->_positions.clear();
    }

    index = _edge2index[dfs_edge];
    child = & (* _children)[index];

    // Add new instance index
    // The case that the edge is found in new instance
    if (child->_occurrence.empty() || child->_occurrence.back() != instance_index)
    {
        child->_occurrence.push_back(instance_index);
        child->_positions.push_back(Positions());
    }

    // Add new position
    child->_positions.back().push_back(
        std::make_shared<const Pattern<Edge>>(Pattern<Edge>{edge_instance,  position})
    );

    return true;
}




Graph GraphEnumerator::pat_to_struct(const Pattern<Symbol> & pattern) const noexcept
{
    Graph graph;
    build_pattern_graph(graph, pattern);
    return graph;
}




std::string GraphEnumerator::pat_to_str(const Pattern<Symbol> & pattern) const noexcept
{
    return _database.struct_to_str(pat_to_struct(pattern));
}




void GraphEnumerator::initialize_visited(const Graph & graph) noexcept
{
    _vertex_is_visited.resize(graph.n_vertices());
    _edge_is_visited.resize(graph.n_edges());
    std::fill(_vertex_is_visited.begin(), _vertex_is_visited.end(), false);
    std::fill(_edge_is_visited.begin(), _edge_is_visited.end(), false);
}