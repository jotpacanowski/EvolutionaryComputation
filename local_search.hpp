#pragma once

#include <utility>
#include <vector>

#include "headers.hpp"

vector<int> steepestLocalSearch(vector<int> solution,
                                const vector<vector<int>> &distanceMatrix,
                                const vector<int> &costs, bool edges = false,
                                int *iterations = nullptr);

vector<int> greedyLocalSearch(vector<int> solution,
                              const vector<vector<int>> &distanceMatrix,
                              const vector<int> &costs, bool edges = false,
                              int *iterations = nullptr);

class SteepestLocalSearchWithCandidateMoves {
   public:
    // IDs of 10 nearest vertices for i-th vertex in terms of distance
    vector<vector<int>> nearest10_only_distance;
    // IDs of 10 nearest vertices for i-th vertex in terms of objective function
    vector<vector<int>> nearest10_objective;

   public:
    // constructor pre-computing the candidates for given instance
    SteepestLocalSearchWithCandidateMoves(const vector<vector<int>> &distanceMatrix,
                                          const vector<int> &costs,
                                          int candidate_count = 10);

    [[nodiscard]] vector<int> do_local_search(
        vector<int> solution, const vector<vector<int>> &distanceMatrix,
        const vector<int> &costs,
        // defaulting to edges-based type of intra-route move
        bool edges = true) const;
};
