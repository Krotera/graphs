// uugraph.tpp: undirected, unweighted graph implementation
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <stack>
#include <queue>
#include <string>
#include <limits>

#include "uugraph.hpp"

/*** Initializes an empty graph ***/
template <typename T>
UUGraph<T>::UUGraph() {
    number_of_edges = 0;
};

/*** Returns number of vertices in this graph ***/
template <typename T>
unsigned int UUGraph<T>::vertSize() {
    return vert_degrees.size(); // Using the unordered_map *may* be faster (likely not)
};

/*** Returns number of edges in this graph ***/
template <typename T>
unsigned int UUGraph<T>::edgeSize() {
    return number_of_edges;
};

/*** Returns whether a vertex exists in this graph ***/
template <typename T>
bool UUGraph<T>::vertExists(const T& target_vert) {
    return (vertices.find(target_vert) != vertices.end());
};

/***
Returns whether an undirected edge mapping between vert_A and vert_B exists in vertices,
checking both for vert_B in vert_A's adjacency set and vice versa,
    returning a bool pair encoded with whether the edge exist (first bool) and, if not,
    whether that was because one or both vertices don't exist in the graph (second bool = false)
    or because the vertices are valid but simply lack an edge to one another (second bool = true)
***/
template <typename T>
std::pair<bool, bool> UUGraph<T>::edgeExists(const T& vert_A, const T& vert_B) {
    // Encode existence of params in second bool
    if (!(vertExists(vert_A) && vertExists(vert_B))) {
        return std::make_pair(false, false);
    }
    else {
        std::pair<bool, bool> ret;
        typename std::map<T, std::set<const T *>>::iterator target_A, target_B;

        ret = std::make_pair(false, true);
        target_A = vertices.find(vert_A);
        target_B = vertices.find(vert_B);
        const T * ptrA = &target_A->first;
        const T * ptrB = &target_B->first;
        // Encode existence of edge in first bool
        if ((target_A->second.find(ptrB) != target_A->second.end()) &&
        (target_B->second.find(ptrA) != target_B->second.end())) {
            ret.first = true;
        }
        return ret;
    }
};

/***
Returns the degree of target_vert,
    returning an int < 0
        if target_vert doesn't exist in this graph
    or an int >= 0
        if target_vert exist and its degree is returned
***/
template <typename T>
int UUGraph<T>::vertDeg(const T& target_vert) {
    if (vertExists(target_vert)) {
        return vert_degrees.at(&vertices.find(target_vert)->first);
    }
    else {
        return -1;
    }
};

/***
Attempts to insert new_vert in vertices as well as inserting a pointer to
it in vertices inside of vert_degrees with degree 0,
    returning false
        if it failed because the vertex was already in vertices
    or true
        if it succeeded
***/
template <typename T>
bool UUGraph<T>::addVert(const T& new_vert) {
    if (!vertExists(new_vert)) {
        vertices.insert(std::pair<T, std::set<const T *>>(new_vert, std::set<const T *>()));
        if (vertices.find(new_vert) != vertices.end()) {
            vert_degrees.insert(std::pair<const T *, unsigned int>(&vertices.find(new_vert)->first, 0));
            if (vert_degrees.find(&vertices.find(new_vert)->first) != vert_degrees.end()) {
                return true;
            }
        }
    }
    return false;
};

/***
Attempts to add an edge between vert_A and vert_B, meaning that a pointer to vert_B is added to
vert_A's adjacency set and vice versa, and then increment number_of_edges and both vertices' degrees
in vert_degrees accordingly
    Returns false
        if it failed because one or both of vert_A and vert_B do not exist in vertices or
        one or both were found in the other's adjacency set, impying the edge already exists
    or true
        if it succeeded
Note: this graph implementation cannot support redundant edges.
 ***/
template <typename T>
bool UUGraph<T>::addEdge(const T& vert_A, const T& vert_B) {
    typename std::map<T, std::set<const T *>>::iterator target_A = vertices.find(vert_A);
    typename std::map<T, std::set<const T *>>::iterator target_B = vertices.find(vert_B);

    // Check if edge is valid via .second and whether it already exists via .first
    if ((edgeExists(vert_A, vert_B).second) &&
        !(edgeExists(vert_A, vert_B).first)) {
            // Insert vert_B in vert_A's adjancecy set
            vertices.at(vert_A).insert(&target_B->first);
            // Insert vert_A in vert_B's adjancecy set
            vertices.at(vert_B).insert(&target_A->first);
            if ((vertices.at(vert_A).find(&target_B->first) != vertices.at(vert_A).end()) &&
            (vertices.at(vert_B).find(&target_A->first) != vertices.at(vert_B).end())) {
                vert_degrees.at(&target_A->first) += 1;
                vert_degrees.at(&target_B->first) += 1;
                number_of_edges++;
                return true;
            }
    }
    return false;
};

/***
Deletes an edge between vert_A and vert_B by removing the pointer to vert_B from vert_A's
adjacency set and vice versa and decrementing number_of_edges and both vertices' degrees,
    returning false
        if the edge didn't exist in the graph
    or true
        if the vertices and their edge existed and the edge was deleted
***/
template <typename T>
bool UUGraph<T>::delEdge(const T& vert_A, const T& vert_B) {
    typename std::map<T, std::set<const T *>>::iterator target_A = vertices.find(vert_A);
    typename std::map<T, std::set<const T *>>::iterator target_B = vertices.find(vert_B);

    if (edgeExists(vert_A, vert_B).first) {
            target_A->second.erase(&target_B->first);
            target_B->second.erase(&target_A->first);
            vert_degrees.at(&target_A->first) -= 1;
            vert_degrees.at(&target_B->first) -= 1;
            number_of_edges--;
            return true;
    }
    else {
        return false;
    }
};

/***
Deletes target_vert from the graph by first deleting all of the edges between it
and its neighbors and then removing it from vertices and vert_degrees,
    returning false
        if target_vert didn't exist in the graph
    or true
        if target_vert existed and was deleted
 ***/
template <typename T>
bool UUGraph<T>::delVert(const T& target_vert) {
    if (vertExists(target_vert)) {
        // Iterate through target_vert's adjacency set (if non-empty), deleting target_vert from
        // the adj sets of all of its neighbors,
        if (!vertices.at(target_vert).empty()) {
            for (const T * adj_vert : vertices.at(target_vert)) {
                delEdge(target_vert, *adj_vert);
            }
        }
        // then delete target_vert from vertices and vert_degrees.
        vert_degrees.erase(&vertices.find(target_vert)->first);
        vertices.erase(target_vert);
        return true;
    }
    else {
        return false;
    }
};

/***
Returns the distances (numbers of edge hops) from target_vert to every other reachable vertex in this graph
as a pair object with an unordered_map and a bool.
    The map is empty and the bool is false
        if the graph is empty or target_vert does not exist in the graph
    or the map is populated and the bool is true
        otherwise
Note: a BFS traversal is done solely on the component containing target_vert, so vertices located
on other components will not have entries in the returned map.
***/
template <typename T>
std::pair<std::unordered_map<const T *, unsigned int>, bool> UUGraph<T>::BFS_DistMapFrom(const T& target_vert) {
    if (vertices.empty() || !vertExists(target_vert)) {
        return std::make_pair(std::unordered_map<const T *, unsigned int>(), false);
    }
    else {
        std::queue<const T *> Q;
        std::unordered_map<const T *, unsigned int> dist; // Distances to other vertices (the map to return)
        std::unordered_map<const T *, std::string> state; // Marks vertices as "unvisited", "visited", and "processed"
        //std::map<const T *, const T *> pred; // Traces traversal paths from the starting vertex to any other (unnecessary here)
        const T * target_ptr = &vertices.find(target_vert)->first; // Get ptr to target_vert in vertices

        // Initialize values
        for (const auto & vert_set_pair : vertices) {
            // vert_set_pair.first holds the vertex keys whose pointers we want to populate dist and state with.
            dist.insert(std::pair<const T *, unsigned int>(&vert_set_pair.first, std::numeric_limits<unsigned int>::max()));
            state.insert(std::pair<const T *, std::string>(&vert_set_pair.first, "unvisited"));
        }

        // Initialize starting vertex's values and add it to queue
        dist.at(target_ptr) = 0;
        state.at(target_ptr) = "visited";
        Q.push(target_ptr);
        // Traverse BFS
        while (!Q.empty()) {
            const T * curr_ptr = Q.front();
            const unsigned int curr_dist = dist.at(Q.front());

            Q.pop();
            for (const T * adj_vert : vertices.at(*curr_ptr)) {
                if (state.at(adj_vert) == "unvisited") {
                    dist.at(adj_vert) = curr_dist + 1;
                    state.at(adj_vert) = "visited";
                    Q.push(adj_vert);
                }
            }
        }
        // dist includes all the vertices in the graph, and the size of dist is compared with the total size of the graph
        // in isConnected() to determine whether the BFS traversed the entire graph and not one of several isolated components.
        // This means that we must scrub dist of all vertices this BFS didn't traverse and whose distances
        // remain std::numeric_limits<unsigned int>::max().
        for (auto vert_dist_pair : dist) {
            if (vert_dist_pair.second == std::numeric_limits<unsigned int>::max()) {
                dist.erase(vert_dist_pair.first);
            }
        }
        return std::make_pair(dist, true);
    }
};

/***
Returns the number of components in the graph (defined here as subgraphs of the graph
for which isConnected() would return true),
    from 0
        for an empty graph
    to beyond
        for a non-empty graph
***/
template <typename T>
unsigned int UUGraph<T>::componentCount() {
    unsigned int component_count = 0;

    if (!vertices.empty()) {
        // Similar to the DFS in isCyclic(), initiate a BFS on each connected component
        std::queue<const T *> Q;
        std::unordered_map<const T *, std::string> state; // Marks vertices as "unvisited" or "visited"

        // Initialize values
        for (const auto & vert_set_pair : vertices) {
            // vert_set_pair.first holds the vertex keys whose pointers we want to populate state and pred with.
            state.insert(std::pair<const T *, std::string>(&vert_set_pair.first, "unvisited"));
        }
        for (const auto & vert_str_pair : state) {
            if (vert_str_pair.second == "unvisited") {
                component_count++; // Entering here means we've stepped into a new component.
                state.at(vert_str_pair.first) = "visited";
                Q.push(vert_str_pair.first);
            }
            while (!Q.empty()) {
                const T * curr_ptr = Q.front();

                Q.pop();
                for (const T * adj_vert : vertices.at(*curr_ptr)) {
                    if (state.at(adj_vert) == "unvisited") {
                        state.at(adj_vert) = "visited";
                        Q.push(adj_vert);
                    }
                }
            }
        }
    }
    return component_count;
};

/***
Determines whether the graph is connected (i.e., whether every pair of any two vertices
forms the endpoints of a path with path defined as a series of edges and vertices with
no repetitions thereof),
    returning false
        if all the vertices in the graph don't belong to a single connected component
        (i.e., if there are islands of vertices unconnected by edges to the rest of the graph)
    or true
        otherwise
***/
template <typename T>
bool UUGraph<T>::isConnected() {
    // Initiate a BFS iteration over the graph from the first vertex in vertices,
    std::pair<std::unordered_map<const T *, unsigned int>, bool> traversal = BFS_DistMapFrom(vertices.cbegin()->first);

    if (traversal.second == false) {
        return true; // I guess a null graph is "connected".
    }
    else {
        // and return whether or not the map compiled from our traversal contains as many vertices as the graph.
        return traversal.first.size() == vertices.size();
    }
};

/***
Determines whether the graph contains cycles (a cycle is defined here as a path of edges
and vertices by which a vertex is reachable from itself but with no repetitions of edges or
vertices except, implicitly, for the starting and ending vertex),
    returning false
        if the graph contains no cycles
    or true
        otherwise
***/
template <typename T>
bool UUGraph<T>::isCyclic() {
    // Following a handy rule of thumb (ONLY FOR CONNECTED GRAPHS!), a cycle exists if |E| > |V| - 1
    if (isConnected()) {
        return number_of_edges > vertices.size() - 1;
    }
    else {
        // Do, for each connected component, a DFS that hunts for a cycle and stops prematurely upon finding one
        std::stack<const T *> S;
        std::unordered_map<const T *, std::string> state; // Marks vertices as "unvisited" and "visited"
        std::unordered_map<const T *, const T *> pred; // Tracks vertices' predecessors to know which visited neighbors to exclude from cycles

        // Initialize values
        for (const auto & vert_set_pair : vertices) {
            // vert_set_pair.first holds the vertex keys whose pointers we want to populate state and pred with.
            state.insert(std::pair<const T *, std::string>(&vert_set_pair.first, "unvisited"));
            pred.insert(std::pair<const T *, const T *>(&vert_set_pair.first, nullptr));
        }
        // Initiate DFS traversals of all components
        // For all unvisited vertices (functionally, the "first" vertex of a component), push them to the stack and DFS traverse their component
        for (const auto & vert_str_pair : state) {
            if (vert_str_pair.second == "unvisited") {
                pred.at(vert_str_pair.first) = vert_str_pair.first; // Set the first vertex's predecessor as itself
                state.at(vert_str_pair.first) = "visited";
                S.push(vert_str_pair.first);
            }
            while (!S.empty()) {
                const T * curr_ptr = S.top();

                S.pop();
                for (const T * adj_vert : vertices.at(*curr_ptr)) {
                    if (state.at(adj_vert) == "unvisited") {
                        pred.at(adj_vert) = curr_ptr;
                        state.at(adj_vert) = "visited";
                        S.push(adj_vert);
                    }
                    // If we have visited neighbors that are not this vertex's predecessor,
                    else if (pred.at(curr_ptr) != adj_vert) {
                        // we have a cycle!
                        return true;
                    }
                }
            }
        }
        return false;
    }
};

/***
Determines whether the graph has an Eulerian cycle (defined, for an undirected graph,
as a graph which is connected and in which every vertex has even degree),
    returning false
        if the graph does not have an Eulerian circuit
    or true
        otherwise
Note: this is distinct from a graph that contains an Eulerian trail (the heuristic
for which is being connected and having exactly zero or two vertices with odd degree).
***/
template <typename T>
bool UUGraph<T>::isEulerian() {
    if (isConnected()) {
        for (const auto & vert_deg_pair : vert_degrees) {
            if (vert_deg_pair.second % 2 != 0) {
                return false;
            }
        }
        return true;
    }
    else {
        return false;
    }
};
