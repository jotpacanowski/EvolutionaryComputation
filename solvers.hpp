#pragma once

#include "headers.hpp"

vector<int> randomTSP(const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
                      int seed);

// static int findNearestNeighbor(const vector<vector<int>>& distanceMatrix,
//                                const vector<int>& costs, int node,
//                                const vector<int>& solution,
//                                const vector<uint8_t>& is_in_sol);

vector<int> nearestNeighborTSP(const vector<vector<int>>& distanceMatrix,
                               const vector<int>& costs, int starting_node);

vector<int> nearestNeighborAnyTSP(const vector<vector<int>>& distanceMatrix,
                                  const vector<int>& costs, int starting_node);

vector<int> greedyCycleTSP(const vector<vector<int>>& distanceMatrix,
                           const vector<int>& costs, int starting_node);

vector<int> regret2TSP(const vector<vector<int>>& distanceMatrix,
                       const vector<int>& costs, int starting_node);