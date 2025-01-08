#pragma once

#include <utility>
#include <vector>

#include "headers.hpp"

vector<int> evolutionarySolver(const vector<vector<int>>& distanceMatrix,
                               const vector<int>& costs, int seed);

vector<vector<int>> getCommonSubsequences(vector<int> v1, vector<int> v2,
                                          vector<int> common_nodes);