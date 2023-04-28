#include "pmopt/databases/graph_database.hpp"

#include <iostream>

using namespace PMOpt;


void GraphDatabase::push(std::string && buffer)
{
    if (_format == FORMAT_DFS)
        push_dfs(std::move(buffer));
}


void GraphDatabase::push_dfs(std::string && buffer)
{
    Graph graph;
    
    // relabel vertex id
    // this is needed to make vertices consecutively numbered 
    std::unordered_map<unsigned int, unsigned int> vertex_id_old2new;

    std::stringstream ss(std::move(buffer));
    while (ss >> buffer)
    {
        if (buffer == "v") // for vertex entry
        {
            unsigned int vertex_id;
            std::string vertex_symbol;

            ss >> vertex_id;
            ss >> vertex_symbol;

            // if vertex_id is already registered in vertex_id_old2new
            // this means that the vertiex is duplicated
            if (vertex_id_old2new.find(vertex_id) != vertex_id_old2new.end())
                throw std::runtime_error("Found duplicate vertex id: " + std::to_string(vertex_id));

            // Set new vertex_id
            // Plus one because of avoiding conflicting ROOT vertex id
            vertex_id_old2new[vertex_id] = vertex_id_old2new.size() + 1;

            // convert vid_old to vid_new
            vertex_id = vertex_id_old2new[vertex_id];

            // set vertex label
            if (_vertex_symbol2label.find(vertex_symbol) == _vertex_symbol2label.end())
            {
                // if the vertex symbol is newly found
                _vertex_symbol2label[vertex_symbol] = _vertex_symbols.size();
                _vertex_symbols.push_back(vertex_symbol);
            }

            auto vertex_label = _vertex_symbol2label[vertex_symbol];

            // push the vertex into the graph
            graph._vertex_labels.push_back(vertex_label); // minus one because id=0 is for ROOT
            graph.push_back(std::vector<Edge>());
        }
        else if (buffer == "e") // for edge
        {
            unsigned int start_vertex_id, end_vertex_id;
            std::string edge_symbol;

            // set vertex id and edge symbol
            ss >> start_vertex_id;
            ss >> end_vertex_id;
            ss >> edge_symbol;

            if (vertex_id_old2new.find(start_vertex_id) == vertex_id_old2new.end())
                throw std::runtime_error("Undefined vertex id: " + std::to_string(start_vertex_id));
            if (vertex_id_old2new.find(end_vertex_id) == vertex_id_old2new.end())
                throw std::runtime_error("Undefined vertex id: " + std::to_string(end_vertex_id));

            // vertex and edge ids
            start_vertex_id = vertex_id_old2new[start_vertex_id];
            end_vertex_id = vertex_id_old2new[end_vertex_id];

            // set edge label
            if (_edge_symbol2label.find(edge_symbol) == _edge_symbol2label.end())
            {
                // if the edge symbol is newly found
                _edge_symbol2label[edge_symbol] = _edge_symbols.size();
                _edge_symbols.push_back(edge_symbol);
            }

            auto edge_label = _edge_symbol2label[edge_symbol];
            unsigned int edge_id = graph._edge_labels.size() + 1; // id = 0 for ROOT


            // push the edge into the graph
            // minus one because id = 0 is researved for ROOT
            graph[start_vertex_id].push_back(Edge{start_vertex_id, edge_id, end_vertex_id});
            if (start_vertex_id != end_vertex_id)
                graph[end_vertex_id].push_back(Edge{end_vertex_id, edge_id, start_vertex_id});
            graph._edge_labels.push_back(edge_label);
        }
        else
            throw std::runtime_error("Undefined symbol identifier: " + buffer);
    }

    emplace_back(std::move(graph));
}

std::string GraphDatabase::struct_to_str(const Graph & graph) const noexcept
{
    std::stringstream ss;
    
    // for vertices
    for (unsigned int vertex_id = 1; vertex_id < graph.n_vertices(); ++vertex_id)
    {
        if (vertex_id > 1)
            ss << ' ';

        ss << "v " << vertex_id << ' ' << _vertex_symbols[graph.vertex_label(vertex_id)];
    }

    // for edges
    for (const auto & edges : graph)
    {
        for (const auto & edge : edges)
        {
            if (!edge.is_forward())
                continue;

            ss << " e " << edge._start_vertex_id
                << ' ' << edge._end_vertex_id
                << ' ' << _edge_symbols[graph.edge_label(edge._edge_id)];
        }
    }

    return ss.str();
}