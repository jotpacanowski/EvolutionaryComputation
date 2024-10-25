#pragma once

#include "headers.hpp"

int intraSwapTwoNodesImpact(const vector<int> &solution, int id1, int id2,
                            const vector<vector<int>> &distanceMatrix,
                            const int solution_size);

int intraSwapTwoEdgesImpact(const vector<int> &solution, int id1, int id2,
                            const vector<vector<int>> &distanceMatrix,
                            const int solution_size);

int interSwapTwoNodesImpact(const vector<int> &solution, int idx, int external_node,
                            const vector<vector<int>> &distanceMatrix,
                            const vector<int> &costs, const int solution_size);

vector<int> localSearch(vector<int> solution, const vector<vector<int>> &distanceMatrix,
                        const vector<int> &costs, const int solution_size,
                        bool edges = false);

vector<int> localSearchGreedy(vector<int> solution,
                              const vector<vector<int>> &distanceMatrix,
                              const vector<int> &costs, const int solution_size,
                              bool edges = false);