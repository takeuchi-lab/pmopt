#pragma once

#include <iostream>
#include <stack>
#include "pmopt/databases/graph_database.hpp"
#include "pmopt/pattern.hpp"

namespace PMOpt
{
    /**
     * @brief Edge struct used in mining
     */
    struct DFSEdge
    {
        unsigned int _start_vertex_id;
        unsigned int _end_vertex_id;
        unsigned int _start_vertex_label;
        unsigned int _edge_label;
        unsigned int _end_vertex_label;

        inline bool operator==(const DFSEdge & edge) const noexcept;
        inline bool operator<(const DFSEdge & edge) const noexcept;

        /**
        * @brief Check if the edge is forward or not
        * @return bool 
        */
        inline bool is_forward() const noexcept
        {
            return _start_vertex_id <= _end_vertex_id;
        }


        /**
        * @brief get Edge instance
        * @return Edge 
        */
        inline Edge edge(unsigned int edge_id) const noexcept
        {
            return Edge{_start_vertex_id, edge_id, _end_vertex_id};
        }


        /**
        * @brief Convert string
        * @return std::string 
        */
        inline std::string str() const noexcept
        {
            std::stringstream ss;
            ss 
                << '(' << _start_vertex_id << ',' << _end_vertex_id
                << ',' << _start_vertex_label
                << ',' << _edge_label
                << ',' << _end_vertex_label << ')';
            return ss.str();
        }
    };
}


// This struct is required to use DFSEdge as the key of unordered_map
// To use DFSEdge as the key, it is also required DFSEdge::operator== to handle hash conflicts
template<>
struct std::hash<PMOpt::DFSEdge>
{
    inline std::size_t operator()(PMOpt::DFSEdge const & edge) const noexcept
    {
        return std::hash<std::string>()(std::string{
            (char) edge._start_vertex_id,
            (char) edge._end_vertex_id,
            (char) edge._start_vertex_label,
            (char) edge._end_vertex_label,
            (char) edge._edge_label
        });
    }
};


namespace PMOpt
{
    struct GraphEnumerator
    {
        using Database = GraphDatabase;

        // In graph mining,
        // the position of the pattern in the graphs is also represented as pattern type.
        // Patern type as edge position and pattern type as pattern itself differ in that:
        // - edge position records ids in the ORIGINAL graph which has the pattern
        // - pattern itself records ids in the SUB graph
        using Positions = std::vector<std::shared_ptr<const Pattern<Edge>>>;

        struct Symbol
        {
            unsigned int _instance_index;
            Edge _edge_instance;
            Edge _edge_pattern;
        };
        using PDB = ProjectedDatabase<Symbol, Positions>;


        static constexpr Symbol ROOT_SYMBOL = Symbol{0, Edge{0, 0, 0}, Edge{0, 0, 0}};


        /**
         * @brief Edge collector used in dfs-code minimality check
         */
        struct DFSCodeChecker
        {
            DFSEdge _current_edge;
            std::vector<std::shared_ptr<const Pattern<Edge>>> _next_positions;
            std::vector<std::shared_ptr<const Pattern<Edge>>> _current_positions;


            /**
             * @brief 
             * 
             * @param edge 
             */
            void initialize(const DFSEdge & edge) noexcept;
            

            /**
             * @brief 
             * @param dfs_edge
             * @param edge
             * @param position
             * @param i_instance
             * @return bool 
             */
            bool push(
                const DFSEdge &,
                const Edge &,
                std::shared_ptr<const Pattern<Edge>>,
                unsigned int
            ) noexcept;
        };


        /**
         * @brief Edge collector used in enumerate new edges
         */
        struct EdgeEnumerator
        {
            std::vector<PDB> * _children;
            const PDB * _parent;
            std::unordered_map<DFSEdge, unsigned int> _edge2index;
            unsigned int _n_children;

            /**
             * @brief 
             * @param[out] children 
             * @param parent
             */
            void initialize(std::vector<PDB> &, const PDB &) noexcept;

            /**
             * @brief 
             * @param dfs_edge
             * @param edge
             * @param position
             * @param i_instance
             * @return bool 
             */
            bool push(
                const DFSEdge &,
                const Edge &,
                std::shared_ptr<const Pattern<Edge>>,
                unsigned int
            ) noexcept;
        };


        const Database & _database;

        DFSCodeChecker _dfscode_checker;
        EdgeEnumerator _edge_enumerator;
        std::vector<Edge> _rightmost_path_pattern;
        std::vector<Edge> _rightmost_path_instance;

        Graph _dfs_check_graph;
        std::vector<DFSEdge> _dfscode;
        std::vector<const Pattern<Symbol> *> _dfscode_patterns;
        std::vector<unsigned int> _vertex_is_visited;
        std::vector<unsigned int> _edge_is_visited;

        /**
         * @brief Construct a new Graph Enumerator object
         * @param database 
         */
        GraphEnumerator(const Database & database) : _database(database) {}


        /**
         * @brief enumerate child patterns of the parent
         * @param[out] children
         * @param parent
         * @return unsigned int 
         */
        unsigned int enumerate(std::vector<PDB> &, const PDB &) noexcept;


        /**
         * @brief build rightmost path with respect to edges constructing the pattern
         * @param[out] rightmost_path
         * @param[out] backward_start_index
         * @param pattern
         */
        void build_rightmost_path(
            std::vector<Edge> &,
            unsigned int &,
            const Pattern<Symbol> *
        ) const noexcept;
        

        /**
         * @brief trace built rightmost path with respect to edges constructiing the instance
         */
        void trace_rightmost_path(
            std::vector<Edge> &,
            std::vector<unsigned int> &,
            std::vector<unsigned int> &,
            const std::vector<Edge> &,
            const Pattern<Edge> *,
            const Pattern<Symbol> *
        ) const noexcept;


        /**
         * @brief Check dfscode minimality
         * @param projected
         * @return bool 
         */
        bool check_dfscode(const PDB &) noexcept;


        /**
         * @brief Build dfscode of pattern
         * @param[out] dfscode
         * @param[out] dfscode_patterns
         * @param pattern
         * @param
         */
        void build_dfscode(
            std::vector<DFSEdge> &,
            std::vector<const Pattern<Symbol> *> &,
            const Pattern<Symbol> &
        ) const noexcept;


        /**
         * @brief 
         * @param[out] graph
         * @param pattern
         */
        void build_pattern_graph(Graph &, const Pattern<Symbol> &) const noexcept;


        /**
         * @brief Search backward-edge
         * @tparam EdgeCollector 
         * @param[out] collector
         * @param graph
         * @param i_instance
         * @return bool 
         */
        template < class EdgeCollector > 
        bool search_backward_edge(
            EdgeCollector &,
            std::shared_ptr<const Pattern<Edge>>,
            const Graph &,
            unsigned int,
            unsigned int
        ) const noexcept;


        /**
         * @brief Search forward-edge
         * @tparam EdgeCollector 
         * @param[out] collector,
         * @param graph
         * @param i_instance
         * @return bool 
         */
        template < class EdgeCollector > 
        bool search_forward_edge(
            EdgeCollector &,
            std::shared_ptr<const Pattern<Edge>>,
            const Graph &,
            unsigned int
        ) const noexcept;


        /**
         * @brief Search single vertex
         * @tparam EdgeCollector 
         * @param[out] collector
         * @param graph
         * @param i_instance
         * @return bool 
         */
        template < class EdgeCollector >
        bool search_vertex(EdgeCollector &, const Graph &, unsigned int) const noexcept;


        /**
         * @brief convert pattern into structure type
         * @return Structure 
         */
        Graph pat_to_struct(const Pattern<Symbol> &) const noexcept;


        /**
         * @brief convert pattern into string
         * @return std::string 
         */
        std::string pat_to_str(const Pattern<Symbol> &) const noexcept;


        /**
         * @brief initialize visited flags vector of vertices and edges
         * @param graph
         */
        void initialize_visited(const Graph &) noexcept;
    };


    inline std::ostream & operator<<(
        std::ostream & os,
        const GraphEnumerator::Symbol & symbol
    )
    noexcept
    {
        return os << symbol._edge_pattern;
    }
}



template < class EdgeCollector >
bool PMOpt::GraphEnumerator::search_backward_edge(
    EdgeCollector & collector,
    std::shared_ptr<const Pattern<Edge>> position,
    const Graph & graph,
    unsigned int i_instance,
    unsigned int backward_start_index
)
const noexcept
{
    // const auto & graph = _database[i_instance];

    // std::cout << _rightmost_path_instance.size() << "/" << backward_start_index << "\n";

    // get rightmost vertex
    auto rightmost_vertex_instance = _rightmost_path_instance[0]._end_vertex_id;
    auto rightmost_vertex_pattern = _rightmost_path_pattern[0]._end_vertex_id;

    // TODO check if self-loop edge is backward ????

    // Must be include index == 0
    // because there may be the other edge that:
    // its starting vertex is equal to the starting vertex of the rightmost edge
    // and its ending vertex is equal to the ending vertex of the rightmost edge.
    // It is not needed to check index such that: index > _backward_start_index
    // because the sub-graph that added these edges dose not achieve minimum dfs-code.
    for (auto i_edge_sign = (int) backward_start_index; i_edge_sign >= 0; --i_edge_sign)
    {
        // Convert to unsigned
        auto i_edge = (unsigned int) i_edge_sign;

        // Enumerate backward edges:
        // backward edge is the edges that its starting vertex is on the rightmost path
        // and its ending vertex is the rightmost vertex.
        auto end_vertex_id_instance = _rightmost_path_instance[i_edge]._end_vertex_id;

        for (const auto & edge_instance : graph[end_vertex_id_instance])
        {
            // the starting vertex of the edge must be equal to the rightmost vertex
            if (edge_instance._end_vertex_id != rightmost_vertex_instance)
                continue;

            // avoid adding duplicate edges
            if (_edge_is_visited[edge_instance._edge_id])
                continue;

            // std::cout << "echeck:" << edge_instance << "/" <<  << "\n";

            // New edge is that:
            // its starting vertex is the rightmost vertex, and
            // its ending vertex is one of the vertices on the rightmost path.
            // Note that label of starting and ending vertex is reversed
            // because the ending vertex of `edge_instance` is the rightmost vertex.
            DFSEdge dfs_edge{
                rightmost_vertex_pattern,
                _rightmost_path_pattern[i_edge]._end_vertex_id,
                graph.vertex_label(rightmost_vertex_instance),
                graph.edge_label(edge_instance._edge_id),
                graph.vertex_label(edge_instance._start_vertex_id)
                // graph.vertex_label(end_vertex_id_instance)
            };

            // std::cout << "NEW_BKEDGE:" << dfs_edge.str() << "\n";
            
            // if failed to push new edge, then return
            // this return value is used if EdgeCollector == DFSCodeCollector
            // in order to check dfscode is minimum or not
            if (!collector.push(dfs_edge, edge_instance, position, i_instance))
                return false;
        }
    }

    return true;
}




template < class EdgeCollector >
bool PMOpt::GraphEnumerator::search_forward_edge(
    EdgeCollector & collector,
    std::shared_ptr<const Pattern<Edge>> position,
    const Graph & graph,
    unsigned int i_instance
)
const noexcept
{
    // const auto & graph = _database[i_instance];

    auto rightmost_vertex_pattern = _rightmost_path_pattern[0]._end_vertex_id;

    for (unsigned int i_edge = 0; i_edge < _rightmost_path_pattern.size(); ++i_edge)
    {
        // The ending vertex of the current edge
        // will be the starting vertex of the newly added edge.
        auto start_vertex_id_instance = _rightmost_path_instance[i_edge]._end_vertex_id;
        auto start_vertex_id_pattern = _rightmost_path_pattern[i_edge]._end_vertex_id;

        for (const auto & edge_instance : graph[start_vertex_id_instance])
        {
            // avoid adding duplicate edges
            if (_edge_is_visited[edge_instance._edge_id])
                continue;

            // avoid adding duplicate vertices
            // because the edge whose ending vertex is already visisted is not forward-edge.
            if (_vertex_is_visited[edge_instance._end_vertex_id])
                continue;

            // New edge
            DFSEdge dfs_edge{
                start_vertex_id_pattern,
                rightmost_vertex_pattern + 1,
                graph.vertex_label(edge_instance._start_vertex_id),
                graph.edge_label(edge_instance._edge_id),
                graph.vertex_label(edge_instance._end_vertex_id)
            };

            // check if successfully pushed
            if (!collector.push(dfs_edge, edge_instance, position, i_instance))
                return false;
        }
    }

    return true;
}




template < class EdgeCollector >
bool PMOpt::GraphEnumerator::search_vertex(
    EdgeCollector & collector,
    const Graph & graph,
    unsigned int i_instance
)
const noexcept
{
    // const auto & graph = _database[i_instance];

    for (unsigned int vertex_id = 1; vertex_id < graph.n_vertices(); ++vertex_id)
    {
        // single vertex edge
        // starting vertex : ROOT vertex
        // edge : ROOT edge
        DFSEdge dfs_edge{
            GraphDatabase::ROOT,
            GraphDatabase::ROOT + 1,
            GraphDatabase::ROOT,
            GraphDatabase::ROOT,
            graph.vertex_label(vertex_id)
        };

        Edge edge{
            GraphDatabase::ROOT,
            GraphDatabase::ROOT,
            vertex_id
        };

        auto position = std::make_shared<const Pattern<Edge>>(Pattern<Edge>{
            Edge{GraphDatabase::ROOT, GraphDatabase::ROOT, GraphDatabase::ROOT},
            std::shared_ptr<const Pattern<Edge>>()
        });

        if (!collector.push(dfs_edge, edge, position, i_instance))
            return false;
    }
    return true;
}




// This operator is required to use DFSEdge as the key of unordered_map
inline bool PMOpt::DFSEdge::operator==(const DFSEdge & edge) const noexcept
{
    return _start_vertex_id == edge._start_vertex_id 
        && _end_vertex_id == edge._end_vertex_id 
        && _start_vertex_label == edge._start_vertex_label 
        && _edge_label == edge._edge_label 
        && _end_vertex_label == edge._end_vertex_label;
}




inline bool PMOpt::DFSEdge::operator<(const DFSEdge & edge) const noexcept
{
    if (!is_forward())
    {
        if (edge.is_forward()) return true;
        
        if (_end_vertex_id < edge._end_vertex_id) return true;
        
        if (_end_vertex_id != edge._end_vertex_id) return false;
        
        if (_edge_label < edge._edge_label) return true;

        return false;
    }
    else if (edge.is_forward())
    {
        if (_start_vertex_id > edge._start_vertex_id) return true;
        
        if (_start_vertex_id != edge._start_vertex_id) return false;

        if (_start_vertex_label < edge._start_vertex_label) return true;

        if (_start_vertex_label != edge._start_vertex_label) return false;

        if (_edge_label < edge._edge_label) return true;

        if (_edge_label != edge._edge_label) return false;

        if (_end_vertex_label < edge._end_vertex_label) return true;
        
        return false;
    }
    return false;
}
