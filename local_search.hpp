#pragma once

#include "headers.hpp"

vector<int> localSearch(vector<int> solution, const vector<vector<int>> &distanceMatrix,
                        const vector<int> &costs, const int solution_size,
                        bool edges = false);

vector<int> localSearchGreedy(vector<int> solution,
                              const vector<vector<int>> &distanceMatrix,
                              const vector<int> &costs, const int solution_size,
                              bool edges = false);
