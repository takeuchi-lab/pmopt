#pragma once

#include <vector>
#include "pmopt/structures/graph.hpp"

namespace PMOpt
{

    struct GraphDatabase : public std::vector<Graph>
    {
        const std::string FORMAT_DFS = "graph";

        static const unsigned int ROOT = 0;

        std::string _format;
        std::vector<std::string> _vertex_symbols;
        std::vector<std::string> _edge_symbols;
        std::unordered_map<std::string, unsigned int> _vertex_symbol2label;
        std::unordered_map<std::string, unsigned int> _edge_symbol2label;


        /**
         * @brief Construct a new Graph Database object
         * @param format name of column that stores strutuerd instance in input csv
         */
        GraphDatabase(const std::string & format) : _format(format)
        {
            // Both vertex and edge, ROOT symbol will be prepared.
            // They are used to represent ROOT pattern and single vertex pattern.
            // In the implementation, ROOT pattern which is empty sub-graph has an edge:
            //  - starting vertex symbol: ROOT
            //  - ending vertex symbol: ROOT
            //  - edge symbol: ROOT
            // and, single vertex pattern has an edge:
            //  - starting vertex symbol: ROOT
            //  - ending vertex symbol: vertex symbol represented by the pattern 
            //  - edge symbol: ROOT
            _vertex_symbols.push_back(std::string());
            _vertex_symbol2label[std::string()] = ROOT;
            _edge_symbols.push_back(std::string());
            _edge_symbol2label[std::string()] = ROOT;
        }


        /**
         * @brief push database to new instance that represented by string
         * @param string
         */
        void push(std::string &&);


        /**
         * @brief push database to new instance that represented by `dfs` format
         * @param string
         */
        void push_dfs(std::string &&);


        /**
         * @brief convert structure instance into string
         * @param graph
         * @return std::string 
         */
        std::string struct_to_str(const Graph &) const noexcept;
    };
}