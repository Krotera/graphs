// uugraph.hpp: an interface for an undirected, unweighted graph supporting generics
// Very useful page: https://en.wikipedia.org/wiki/Glossary_of_graph_theory_terms
#ifndef UUGRAPH_HPP
#define UUGRAPH_HPP

#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <utility> // std::pair

template <typename T>
class UUGraph {
private:
    std::map<T, std::set<const T *>> vertices; // Map keys are vertices, values are edges (adjacency sets)
    std::unordered_map<const T *, unsigned int> vert_degrees; // Mapping of vertex degrees (no reason for order due to pointer-based keys)
    unsigned int number_of_edges;

public:
    // Tier 1 functions
    UUGraph();
    unsigned int vertSize();
    unsigned int edgeSize();
    bool vertExists(const T&);
    std::pair<bool, bool> edgeExists(const T&, const T&);
    int vertDeg(const T&);
    bool addVert(const T&);
    bool addEdge(const T&, const T&);
    bool delEdge(const T&, const T&);
    bool delVert(const T&);
    // Tier 2 functions
    std::pair<std::unordered_map<const T *, unsigned int>, bool> BFS_DistMapFrom(const T&);
    unsigned int componentCount();
    bool isConnected();
    bool isCyclic();
    bool isEulerian();
};

#include "uugraph.tpp"

#endif
