#pragma once

#include <sstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include "pmopt/utils/containers.hpp"

namespace PMOpt
{

    struct Edge
    {
        unsigned int _start_vertex_id;
        unsigned int _edge_id;
        unsigned int _end_vertex_id;

        inline bool operator==(const Edge & edge) const noexcept
        {
            return _start_vertex_id == edge._start_vertex_id
                && _edge_id == edge._edge_id
                && _end_vertex_id == edge._end_vertex_id;
        }


        /**
         * @brief Check if the edge is `forward-edge` or not
         * @return bool 
         */
        inline bool is_forward() const noexcept
        {
            return _start_vertex_id <= _end_vertex_id;
        }


        /**
         * @brief construct reversed edge
         * @return Edge 
         */
        inline Edge reverse() const noexcept
        {
            return Edge{_end_vertex_id, _edge_id, _start_vertex_id};
        }
    };


    struct Graph : public std::vector<std::vector<Edge>>
    {
        using boostGraph = boost::adjacency_list<
            boost::vecS, // container type for edges for each vertex
            boost::vecS, // container type for vertices
            boost::undirectedS, // directed or undirecetd
            unsigned int, // label type of vertices
            unsigned int, // label type of edges
            boost::vecS // container type for edges for the graph
        >;

        std::vector<unsigned int> _vertex_labels;
        std::vector<unsigned int> _edge_labels;
        boostGraph _boost_graph;


        inline const std::vector<Edge> & operator[](std::size_t id) const noexcept
        {
            // const version
            // minus one because id=0 is researved for ROOT
            assert(id > 0);
            return std::vector<std::vector<Edge>>::operator[](id - 1);
        }


        inline std::vector<Edge> & operator[](std::size_t id) noexcept
        {
            // unconst version
            // minus one because id=0 is researved for ROOT
            assert(id > 0);
            return std::vector<std::vector<Edge>>::operator[](id - 1);
        }


        unsigned int & vertex_label(unsigned int id) noexcept
        {
            // unconst version
            // minus one because id==0 is researved for ROOT
            // Do not use graph._vertex_labels[id] !!!
            assert(id > 0);
            return _vertex_labels[id - 1];
        }


        unsigned int  vertex_label(unsigned int id) const noexcept
        {
            // const version
            // minus one because id==0 is researved for ROOT
            // Do not use graph._vertex_labels[id] !!!
            assert(id > 0);
            return _vertex_labels[id - 1];
        }


        unsigned int & edge_label(unsigned int id) noexcept
        {
            // unconst version
            // minus one because id==0 is researved for ROOT
            // Do not use graph._edge_labels[id] !!!
            assert(id > 0);
            return _edge_labels[id - 1];
        }


        unsigned int edge_label(unsigned int id) const noexcept
        {
            // const version
            // minus one because id==0 is researved for ROOT
            // Do not use graph._edge_labels[id] !!!
            assert(id > 0);
            return _edge_labels[id - 1];
        }


        unsigned int n_vertices() const noexcept
        {
            // Do not use graph._vertex_labels.size() !!!
            return _vertex_labels.size() + 1;
        }

        unsigned int n_edges() const noexcept
        {
            // Do not use graph._edge_labels.size() !!!
            return _edge_labels.size() + 1;
        }

        /**
         * @brief build boost library graph type
         */
        void build_boost() noexcept;


        /**
         * @brief Check if the graph contains subgraph 
         * @param subgraph
         * @return bool
         */
        bool contains(const Graph &) const noexcept;

    };


    inline std::ostream & operator<<(std::ostream & os, const Edge & edge) noexcept
    {
        return os 
            << '(' << edge._start_vertex_id
            << ',' << edge._edge_id
            << ',' << edge._end_vertex_id << ')';
    }


    inline std::ostream & operator<<(std::ostream & os, const Graph & graph) noexcept
    {
        os << "v: " << Utils::to_string(graph._vertex_labels) << '\n';
        os << "e: " << Utils::to_string(graph._edge_labels) << '\n';
        os << "adj:\n";
        for (const auto & edges : graph)
            os << ' ' << Utils::to_string(edges) << '\n';
        return os;
    }
}