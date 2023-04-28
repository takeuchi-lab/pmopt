#include "pmopt/structures/graph.hpp"

using namespace PMOpt;

void Graph::build_boost() noexcept
{
    // convert the graph into boost compatible graph
    auto n_vertices = _vertex_labels.size();

    // add vertices
    std::vector<boostGraph::vertex_descriptor> boost_vertices(n_vertices);
    for (unsigned int i_vertex = 0; i_vertex < n_vertices; ++i_vertex)
    {
        // set vertex_labels[i] into _boost_graph
        // and store the result in boost_vertices[i]
        boost_vertices[i_vertex] = boost::add_vertex(_vertex_labels[i_vertex], _boost_graph);
    }

    // add edges
    std::vector<unsigned int> edges_added(_edge_labels.size(), false);
    for (const auto & edges : * this)
    {
        for (const auto & edge : edges)
        {
            if (edges_added[edge._edge_id])
                continue;

            // set _edge_labels[edge._edge_id] into _boost_graph
            auto v_1 = boost_vertices[edge._start_vertex_id];
            auto v_2 = boost_vertices[edge._end_vertex_id];
            boost::add_edge(v_1, v_2, _edge_labels[edge._edge_id], _boost_graph);

            edges_added[edge._edge_id] = true;
        }
    }
}