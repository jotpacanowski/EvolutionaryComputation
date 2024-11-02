#pragma once

#include "headers.hpp"

vector<int> steepestLocalSearch(vector<int> solution,
                                const vector<vector<int>> &distanceMatrix,
                                const vector<int> &costs, bool edges = false);

vector<int> greedyLocalSearch(vector<int> solution,
                              const vector<vector<int>> &distanceMatrix,
                              const vector<int> &costs, bool edges = false);
